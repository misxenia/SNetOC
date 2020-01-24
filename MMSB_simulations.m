%% Sparse Networks with Overlapping Communities (SNetOC) package: MMSB_simulations
%
% This Matlab script performs posterior inference on a network simulated from the mixed-membership stochastic blockmodel, 
% and performs posterior inference using the same model (well-specified
% case) and the model of Todeschini et al. (misspecified case).
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
%

clear
close all

istest = true; % enable testing mode: quick run with small nb of iterations

root = '.';
if istest
    outpath = fullfile(root, 'results', 'MMSB_simulations', 'test');
else
    outpath = fullfile(root, 'results', 'MMSB_simulations', date);
end

if ~isdir(outpath)
    mkdir(outpath)
end


% Add path
addpath ./GGP/ ./CGGP/ ./utils/ ./MMSB

set(0, 'DefaultAxesFontSize', 14)

% Set the seed
rng default

%% Sample graph from MMSB
%
n = 300;
p = 3;
alpha_true = .1;
W_true = [.7, 0.3, .05;.3, .7, .1; .05, .1, .7];
rho_true = 0;
objtrue = graphmodel('MMSB', n, p, alpha_true, W_true, rho_true);
[G, s_true, pi_true] = graphrnd(objtrue);
labels = {'Nodes', 'Nodes'};


% Plot adjacency matrix
figure
spy(G);
xlabel(labels{2})
ylabel(labels{1})

% Plot sorted nodes
[~, ind_com] = max(pi_true, [], 2);
plot_sortedgraph(G, ind_com, ind_com, 1:p, labels);


%% Posterior inference using the MMSB as prior (well-specified)
%

% Model parameters
n = size(G, 1);
alpha = [];
W = [0.01, 0.01];
rho = 0;

% MCMC parameters
if istest
    niter = 2000;
    nsamples = 50;
else
    niter = 200000;  
    nsamples = 500;
end
nburn = floor(niter*2/4);
thin = ceil((niter-nburn)/nsamples);
nchains = 3;
verbose = true;

% Run MCMC
objprior =  graphmodel('MMSB', n, p, alpha, W,rho);
objmcmc = graphmcmc(objprior, niter, nburn, thin, nchains);
objmcmc = graphmcmcsamples(objmcmc, G, verbose);

% Get estimates of pis
estimates = graphest(objmcmc);

%%

% Plots
%

prefix = '';
suffix = '';

% Assign max feature
[~, nodefeat] = max(estimates.pi, [],2);
plot_sortedgraph(G, nodefeat, nodefeat , 1:p, labels);

% Plot traces and histograms

%%

% Block matrix W
graph_title='W matrix';
figure('name', graph_title)
for i=1:p
    for j=1:p
        subplot(p,p, p*(i-1)+j)
        for k=1:nchains
            plot(squeeze(objmcmc.samples(k).W(i,j,:)));
            ylim([0,1])
            hold on
        end
    end
end
legend('chain 1', 'chain 2','chain 3')

%%

% Dirichlet parameter alpha
graph_title = 'Dirichlet parameter alpha'; figure('name',graph_title)
for k=1:nchains
    plot(squeeze(objmcmc.samples(k).alpha));
    hold on
end
legend('chain 1', 'chain 2','chain 3')
xlabel('Samples')
ylabel('\alpha')

%% Posterior inference using the CGGP model as prior (mis-specified)
%

objprior_CGGP =  graphmodel('CGGP', p);
objmcmc_CGGP = graphmcmc(objprior_CGGP, niter, nburn, thin, nchains);
objmcmc_CGGP = graphmcmcsamples(objmcmc_CGGP, G, verbose);

estimates_CGGP = graphest(objmcmc_CGGP);

%%

% Plots
%

prefix = '';
suffix = '';

% Assign max feature
[~, nodefeat_CGGP] = max(estimates_CGGP.w, [],2);
plot_sortedgraph(G, nodefeat_CGGP, nodefeat_CGGP , 1:p, labels);
