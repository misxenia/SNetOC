%% Sparse Networks with Overlapping Communities (SNetOC) package: MMSB_usairport
%
% This Matlab script performs posterior inference on a network of airports 
% under the mixed membership stochastic blockmodel.
%
% For downloading the package and information on installation, visit the
% <https://github.com/OxCSML-BayesNP/SNetOC SNetOC webpage>.
% 
% References:
%
% * <http://jmlr.csail.mit.edu/papers/volume9/airoldi08a/airoldi08a.pdf>
% * A. Todeschini, X. Miscouridou and F. Caron (2017) <https://arxiv.org/abs/1602.02114 Exchangeable Random Measures for Sparse and Modular Graphs with Overlapping Communities>. arXiv:1602.02114.
%
% Authors: 
%
% * <http://adrien.tspace.fr/ A. Todeschini>, Inria
% * <http://csml.stats.ox.ac.uk/people/miscouridou/ X. Miscouridou>, University of Oxford
% * <http://www.stats.ox.ac.uk/~caron/ F. Caron>, University of Oxford
% 
% Tested on Matlab R2017a. Requires the Statistics toolbox.
%
% Last Modified: 2017-09-15
%%

%% General settings

clear
close all;

istest = true; % enable testing mode: quick run with small nb of iterations

root = '.';
if istest
    outpath = fullfile(root, 'results', 'MMSB_usairport', 'test');
else
    outpath = fullfile(root, 'results', 'MMSB_usairport', date);
end

if ~isdir(outpath)
    mkdir(outpath)
end

% Add path
addpath ./GGP/ ./CGGP/ ./utils/ ./MMSB/

set(0, 'DefaultAxesFontSize', 14)

% Set the seed
rng default

%% Load  USairports network 2010
% 

load ./data/usairport/usairport

titlenetwork = 'US airport network in 2010';
name = 'usairport';
labels = {'Airports','Airports'};

G = G | G'; % make undirected graph
G = logical(G-diag(diag(G))); % remove self edges (#164)

meta.degree = num2cell(full(sum(G,2)));

fn = fieldnames(meta);

% Remove nodes with no edge (#0)
ind = any(G);
G = G(ind, ind);
for i=1:length(fn)
    meta.(fn{i}) = meta.(fn{i})(ind,:);
end

% Plot adjacency matrix
figure
spy(G);
xlabel(labels{2})
ylabel(labels{1})


%% Posterior inference
%

% Model parameters
n = size(G, 1);
p = 4;
alpha = [];
W = [0.01, 0.01];
rho = [];

% MCMC parameters
if istest
    niter = 200;
    nsamples = 50;
else
    niter = 200000;
    nsamples = 500;    
end
nchains = 3;
nburn = floor(niter*2/4);
thin = ceil((niter-nburn)/nsamples);

verbose = true;

% Run MCMC
objprior =  graphmodel('MMSB', n, p, alpha, W, rho);
objmcmc = graphmcmc(objprior, niter, nburn, thin, nchains);
objmcmc = graphmcmcsamples(objmcmc, G, verbose);

% Get estimates of pis
[estimates, ~] = graphest(objmcmc);

%% Plots
%
prefix = sprintf('%s_%df_', name, p);
suffix = '';

% Assign max feature
[~, nodefeat] = max(estimates.pi, [],2);

% order features - features are not really interpretable here
[~, ind_features] = sort(sum(estimates.pi, 1), 'descend');
featnames = {'Feat. 1', 'Feat. 2', 'Feat. 3', 'Feat. 4'};
plot_sortedgraph(G, nodefeat, nodefeat , ind_features);

figure
for k=1:p
    subplot(2,2,k)
    plot(sum(G), estimates.pi(:,k), '.');
    xlabel('Degree node');
    ylabel(['\pi_' num2str(k)]);
    title(['Feature' num2str(k)]);
end


% Plot traces and histograms

% Block matrix
graph_title='W matrix';figure('name', graph_title)
for i=1:p
    for j=1:p
        subplot(p,p, p*(i-1)+j)
        for k=1:nchains
            plot(squeeze(objmcmc.samples(k).W(i,j,:)));ylim([0,1])
            hold on
        end
    end
end
legend('chain 1', 'chain 2','chain 3')


% dirichlet parameter alpha
graph_title = 'Dirichlet parameter alpha'; figure('name',graph_title)
for k=1:nchains
    plot(squeeze(objmcmc.samples(k).alpha));ylim([0,0.5])
    hold on
end
legend('chain 1', 'chain 2', 'chain 2')

% sparsity parameter rho
graph_title = 'sparsity parameter rho';figure('name',graph_title)
for k=1:nchains
    plot(squeeze(objmcmc.samples(k).rho));ylim([0,1])
    hold on
end
legend('chain 1', 'chain 2','chain 3')

plot_GOF_graphs=false;
if plot_GOF_graphs
    % Plot posterior predictive of degrees
    cond = false;
    verbose = true;
    [~, ~, degsamp, sd_rowmean, triad_dep] = plot_degreepostpredMMSB(G, objmcmc.samples, nsamples, cond, outpath, prefix, suffix, verbose);

    % empirical goodness of fit statistics (see P Hoff 200?)
    [TRUEsd_rowmean,TRUEsd_colmean,TRUEdyad_dep,TRUEtriad_dep] = GOFstats(G); 

    %  transitivity
    figure('name', 'cluster coefficient');
    myhist_fit(triad_dep,10, 'normal','-+g')
    hold on;
    y1=get(gca,'ylim');
    hold on
    plot([TRUEtriad_dep TRUEtriad_dep],y1,'--b')
    legend('MMSB ','empirical','location','northeast')
    box off;
    savefigs(gcf,  [prefix 'cluster coefficient' suffix], outpath);


    %sd deviation of degree
    figure('name', 'sd_deviation');
    myhist_fit(sd_rowmean,10, 'normal','-+g')
    hold on;
    y1=get(gca,'ylim');
    hold on
    plot([TRUEtriad_dep TRUEtriad_dep],y1,'--b')
    legend('MMSB ','empirical','location','northeast')
    box off;
    savefigs(gcf,  [prefix 'sd_degree' suffix], outpath);
end