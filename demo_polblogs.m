%% Sparse Networks with Overlapping Communities (SNetOC) package: demo_polblogs
%
% This Matlab script performs posterior inference on a network of political
% blogs to find latent overlapping communities, using a Bayesian
% nonparametric approach.
%
% For downloading the package and information on installation, visit the
% <https://github.com/OxCSML-BayesNP/SNetOC SNetOC webpage>.
% 
% Reference: 
%
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
% Last Modified: 01/2020
%%

%% General settings
%

clear
close all
tstart = clock; % Starting time

istest = true; % enable testing mode: quick run with smaller nb of iterations

%%
% In test mode, a smaller number of iterations is run. Although the sampler
% clearly has not converged yet, the method already recovers well the ground truth
% communities. To reproduce the results of the paper, set this value to
% false. 

root = '.';
if istest
    outpath = fullfile(root, 'results', 'demo_polblogs', 'test');
else
    outpath = fullfile(root, 'results', 'demo_polblogs', date);
end

if ~isdir(outpath)
    mkdir(outpath);
end

% Add path
addpath ./GGP/ ./CGGP/ ./utils/

set(0, 'DefaultAxesFontSize', 14)

% Set the seed
rng default

%% Load network of political blogs
% Data can be downloaded from
% <https://www.cise.ufl.edu/research/sparse/mat/Newman/polblogs.mat here>.

load ./data/polblogs/polblogs.mat

titlenetwork = 'Political blogosphere Feb. 2005';
name = 'polblogs';
labels = {'Blogs', 'Blogs'};
groupfield = 'name'; % meta field displayed for group plot

% Transform the graph to obtain a simple graph
G = Problem.A | Problem.A'; % make undirected graph
G = logical(G-diag(diag(G))); % remove self edges (#3)

% Collect metadata
clear meta
meta.name = cellstr(Problem.aux.nodename);
meta.source = cellstr(Problem.aux.nodesource);
meta.isright = logical(Problem.aux.nodevalue);
meta.degree = num2cell(full(sum(G,2)));
meta.groups = zeros(size(meta.isright));
meta.groups(~meta.isright) = 1;
meta.groups(meta.isright) = 2;
color_groups = [0 0 .8; .8 0 0];
label_groups = {'Left', 'Right'};
fn = fieldnames(meta);

% Remove nodes with no edge (#266)
ind = any(G);
G = G(ind, ind);
for i=1:length(fn)
    meta.(fn{i}) = meta.(fn{i})(ind);
end

% Plot adjacency matrix (sorted)
figure('name', 'Adjacency matrix (sorted by ground truth political leaning)');
spy(G);
xlabel(labels{2})
ylabel(labels{1})

%%

% Shuffle nodes: irrelevant due to exchangeability, just to check we do not cheat!
ind = randperm(size(G,1));
G = G(ind, ind);
for i=1:length(fn)
    meta.(fn{i}) = meta.(fn{i})(ind);
end

% Plot adjacency matrix (unsorted)
figure('name', 'Adjacency matrix (unsorted)');
spy(G);
xlabel(labels{2})
ylabel(labels{1})

%%

% Plot degree distribution
figure('name', 'Empirical degree distribution');
hdeg = plot_degree(G);
set(hdeg, 'markersize', 10, 'marker', 'o','markeredgecolor', 'none', 'markerfacecolor', [1, .75, .75]);

%% Posterior Inference using Markov chain Monte Carlo and point estimation
% Users needs to start the parallel pool by using the command *parpool* to run multiple chains in parallel.

% Define the parameters of the prior
p = 2; % Number of commmunities
objprior =  graphmodel('CGGP', p); % CGGP graph model with p communities

% Define parameters of the MCMC sampler
nchains = 3;
if istest
    niterinit = 1000;
    niter = 10000;
    nsamples = 100;
    ndraws = 100; 
else
    niterinit = 10000;
    niter = 2e6;
    nsamples = 1000;
    ndraws = 500;
end   
nburn = floor(niter/2);
thin = ceil((niter-nburn)/nsamples); 
verbose = true;

% Create the graphMCMC object
objmcmc = graphmcmc(objprior, niter, 0, thin, nchains);

% Run initialisation
init = graphinit(objmcmc, G, niterinit);

%%

% Run MCMC sampler
objmcmc = graphmcmcsamples(objmcmc, G, verbose, init);

%%

% Print summary in text file
print_summary(['summary_' num2str(p) 'f.txt'], titlenetwork, G, niter, nburn, nchains, thin, p, outpath, tstart)

% Save workspace
save(fullfile(outpath, ['workspace_' num2str(p) 'f.mat']), '-v7.3')

%% 

% Log posterior approximation
[lp_nonlat, lp_lat, ll_nonlat, ll_lat] = logpost_approx(objmcmc, G);

%% 

% Compute identifiable parameters
for ch = 1:size(objmcmc.samples, 2)
    S = objmcmc.samples(1,ch);
    objmcmc.samples(1,ch).varsigma1 = -S.alpha .* S.tau.^S.sigma ./ S.sigma;
    a_sigma = S.Fparam.a.*S.sigma;
    objmcmc.samples(1,ch).varsigma2 = -a_sigma./S.Fparam.b2;
    objmcmc.samples(1,ch).varsigma3 = a_sigma.*(S.sigma-S.Fparam.a-1)./(S.Fparam.b2.^2);
end

    
%% discard burnin
objmcmc_noburn = objmcmc;
objmcmc_noburn.samples = discard(objmcmc.samples, floor(nburn/objmcmc.settings.thin));
objmcmc_noburn.settings.nburn = nburn;

%% Point estimation of the model parameters
[estimates, C_st] = graphest(objmcmc_noburn);

%% Plots 

prefix = sprintf('%s_%df_', name, p);
suffix = '';

%%

% Plot Log posterior approximation
iter = (1:size(lp_nonlat,1))*thin;
plot_logpost(lp_nonlat, iter, [], 'Log-posterior', outpath, prefix, '_nonlat');
plot_logpost(lp_lat, iter, [], 'Log-posterior', outpath, prefix, '_lat');

% Plot log-posterior autocorr
lp_nonlat_noburn = lp_nonlat(floor(nburn/niter*size(lp_nonlat, 1)):end, :);
lp_lat_noburn = lp_lat(floor(nburn/niter*size(lp_lat, 1)):end, :);
plot_autocorr_logpost(lp_nonlat_noburn, thin, 'Log-posterior', outpath, prefix, '_nonlat');
plot_autocorr_logpost(lp_lat_noburn, thin, 'Log-posterior', outpath, prefix, '_lat');

% Plot cost
if ~isempty(C_st)
    plot_cost(C_st, outpath, prefix, suffix);
end

%%

% Identify each feature to left/right wing using ground truth 
% (This step would normally require a human interpretation of the features)
[~, ind_features] = sort(median(estimates.w(meta.isright,:), 1)./median(estimates.w, 1));
featnames = {'Liberal', 'Conservative'};

% Print classification performance with ground truth
[~, nodefeat] = max(estimates.w, [],2); % Assign each node to the feature with highest weight
[confmat] = print_classif(fullfile(outpath, ['classif_' num2str(p) 'f.txt']), ...
    nodefeat, meta.groups, ind_features, label_groups);

%%

% Plot traces and histograms
variables = {'varsigma1', 'varsigma2', 'varsigma3', 'mean_w_rem'};
namesvar = {'$\varsigma_1$', '$\varsigma_2$', '$\varsigma_3$', '$\overline{w}_{\ast}$'};
plot_trace(objmcmc.samples, objmcmc.settings, variables, namesvar, [], outpath, prefix, suffix);
plot_hist(objmcmc_noburn.samples, variables, namesvar, [], ind_features, [], outpath, prefix, suffix);

%%

% Plot the graph by sorting the nodes by max feature to see block structure
plot_sortedgraph(G, nodefeat, nodefeat, ind_features, labels, outpath, prefix, suffix, {'png'});

%%

% Show the proportion in each features for a few nodes
if p==2
    color = color_groups;
elseif p==5
    color = [0 0 .5; .3 .3 1 ; 0.8 0.8 0.8 ; 1 .3 .3; .5 0 0];
end

if isfield(meta, 'groups')
    % Plots by groups right vs left
    plot_groups(estimates.w, meta.groups, meta.(groupfield), ind_features, label_groups, featnames, ...
        color_groups, outpath, prefix, suffix);
    
    y = sum(estimates.w,2);
    x = estimates.w(:,ind_features(2))./y;
    figure; hold on
    for i=1:numel(label_groups)
        ind = meta.groups==i;
        stem(x(ind), y(ind), 'o', 'markersize', 5, 'color', color(i,:))
    end
    xlabel('$w_{2}/\vert w \vert$', 'interpreter', 'latex', 'fontsize', 20)
    ylabel('$\vert w \vert$', 'interpreter', 'latex', 'fontsize', 20)
    text([0,.75], [-.5, -.5], featnames, 'fontsize', 16)
    legend(label_groups, 'fontsize', 16)
    legend boxoff
    axis tight
    xlim([0,1])
end

%%

% Plot normalized weights for a subset of blogs
prop_nodes = estimates.w(:,1)./sum(estimates.w, 2);
names = {'blogsforbush.com'
    'instapundit.com'
    'drudgereport.com'
    'tagorda.com'
    'danieldrezner.com/blog'
    'andrewsullivan.com'
    'iwantmycountryback.org'
    'democraticunderground.com'
    'wonkette.com'
    'washingtonmonthly.com'
    'dailykos.com'
    'atrios.blogspot.com'};
ind = zeros(size(names,1),1);
lab = cell(size(names,1), 1);
for i=1:size(names,1)
    ind(i) = find(strcmp(meta.name, names{i}) & strcmp(meta.name, names{i,1}));
    lab{i} = sprintf('%s (%s)', names{i}, label_groups{meta.groups(ind(i))}(1));
%     fprintf('%2d. (%.2f) #%d, %s\n', i, prop_nodes(ind(i)), meta.degree{ind(i)}, lab{i});
end
plot_nodesfeatures(estimates.w, ind, ind_features, lab, featnames, color, outpath, prefix, suffix);

%%

% Show blogs with highest weight in each feature
meta.wing = cell(meta.name);
meta.wing(meta.isright) = repmat({'R'}, sum(meta.isright), 1);
meta.wing(~meta.isright) = repmat({'L'}, sum(~meta.isright), 1);

fnames = {'degree', 'name', 'wing'};
formats = {'#%d,', '%s', '(%s).'}; % 
% Nodes with highest weight in each feature
fprintf('-----------------------------------\n')
fprintf('Nodes with highest weights in each feature\n')
fprintf('#edges, name of blog (Political leaning)\n')
fprintf('-----------------------------------\n')
print_features( [outpath 'features_' num2str(p) 'f.txt'], ...
    estimates.w, ind_features, [], meta, fnames, formats );
fprintf('-----------------------------------\n')

%%

% Correlation between the features
figure('name', 'Correlation between features');
imagesc(corr(estimates.w(:,ind_features)));
colormap('gray')
caxis([0,1])
colorbar
xlabel('Feature')
ylabel('Feature')
title('Correlation between features')

%%

% Plot posterior predictive of degree distribution
plot_degreepostpred(G, objmcmc_noburn, ndraws, 1e-6, outpath, prefix, suffix);
