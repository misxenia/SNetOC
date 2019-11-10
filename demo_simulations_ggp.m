%% Sparse Networks with Overlapping Communities (SNetOC) package: demo_simulations
%
% This Matlab script performs posterior inference on a (sparse) simulated graph with overlapping communities.
%
% For downloading the package and information on installation, visit the
% <https://github.com/misxenia/SNetOC SNetOC webpage>.
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
% Tested on Matlab R2016a. Requires the Statistics toolbox.
%
% Last Modified: 10/10/2017
%%

%% General settings
%

close all
clearvars

tstart = clock; % Starting time

istest = true; % enable testing mode: quick run with small nb of iterations

root = '.';
if istest
    outpath = fullfile(root, 'results', 'GGP_simulations', 'test');
else
    outpath = fullfile(root, 'results', 'GGP_simulations', date);
end

if ~isdir(outpath)
    mkdir(outpath)
end

% Add path
addpath ./GGP/ ./utils/ 

set(0, 'DefaultAxesFontSize', 14)

% Set the seed
rng default

%% Sample a GGP graph

titlenetwork = 'Simulated simple graph';
name = 'simugraph';
labels = {'Nodes', 'Nodes'};

% Set parameters of the GGP graph
alpha = 100; sigma = 0.1; tau = 1;
obj = graphmodel('GGP', alpha, sigma, tau);
% Sample a GGP graph
[G, w_true, w_rem_true] = graphrnd(obj, 1e-6);

nnodes = size(G, 1);
nedges = nnz(G);
fprintf('GGP graph with %d nodes and %d edges sampled\n', size(G, 1), nnz(triu(G)));

figure('name', 'adjacency matrix')
spy(G)
xlabel('Nodes')
ylabel('Nodes')
title('Adjacency Matrix')

figure('name', 'Empirical degree distribution')
hdeg = plot_degree(G);
set(hdeg, 'markersize', 10, 'marker', 'o','markeredgecolor', 'none', 'markerfacecolor', [1, .75, .75]);

%% Prior distribution

objprior =  graphmodel('GGP');

%% Posterior inference

% Parameters of the MCMC algorithm
if istest
    niterinit = 100;
    niter = 20000;
    nsamples = 100; % Nb of Monte Carlo samples to return
else
    niterinit = 10000;
    niter = 200000;
    nsamples = 500;
end   
nburn = floor(3*niter/4); nchains = 3;
thin = ceil((niter-nburn)/nsamples);
verbose = true;

% Create the graphMCMC object
objmcmc = graphmcmc(objprior, niter, nburn, thin, nchains); 
% Note: nburn is set to zero here in order to store samples in the transient regime of the MCMC

%% Run MCMC sampler
objmcmc = graphmcmcsamples(objmcmc, G, verbose);

%%

% Print summary in text file
print_summary(['summary_' num2str(p) 'f.txt'], titlenetwork, G, niter, nburn, nchains, thin, p, outpath, tstart)

% discard burnin to compute estimates
objmcmc_all = objmcmc;
objmcmc.samples = discard(objmcmc_all.samples, floor(nburn/objmcmc_all.settings.thin));
objmcmc.settings.nburn = nburn;

% Get estimates and cost
[estimates, C_st] = graphest(objmcmc);

% Save workspace
save(fullfile(outpath, ['workspace_' num2str(p) 'f.mat']))

%% Plots
%

prefix = sprintf('%s_%df_', name, p);
suffix = '';

% Plot cost
plot_cost(C_st, outpath, prefix, suffix);

% Assign max feature
[~, nodefeat] = max(estimates.w, [],2);

% order features 
ind_features = 1:p;

%%

% Plot traces and histograms
variables = {'logalpha2', 'sigma', 'Fparam.a', 'Fparam.b2', 'mean_w_rem'};
namesvar = {'$\log \tilde\alpha$', '$\sigma$', '$a$', '$\tilde b$', '$\overline{w}_{\ast}$'};
trueval = {log(alpha_true)+sigma_true*log(tau_true), sigma_true,Fdist_true.param.a, tau_true*Fdist_true.param.b, mean(w_rem_true)};
plot_trace(objmcmc_all.samples, objmcmc_all.settings, variables, namesvar, trueval, outpath, prefix, suffix);
plot_hist(objmcmc.samples, variables, namesvar, trueval, ind_features, [], outpath, prefix, suffix);

%%

% Plot the graph by sorting the nodes by max feature
[~,indk] = max(w_true,[],2);
[~,ind] = sort(indk, 'descend');
figure
spy(G(ind, ind))
xlabel('Nodes')
ylabel('Nodes')
title('Sorted Adjacency Matrix')

%%

% Plot posterior predictive of degrees
if istest
    ndraws = 100;
else
    ndraws = 5000;
end
plot_degreepostpred(G, objmcmc, ndraws, 1e-6, outpath, prefix, suffix);

%%

% Plot credible intervals for the weights
[~, ind] = sort(sum(G, 2)+sum(G, 1)', 'descend');
mean_w_true = mean(w_true,2);

% High degree nodes
figure('name','Credible intervals - high degree nodes'); hold on
for i=1:min(nnodes(1), 50)
    plot([i, i],quantile(objmcmc.samples(1).mean_w(ind(i), :, :),[.025,.975]), 'r','linewidth', 3);
    hold on
    plot(i, mean_w_true(ind(i)), 'xg', 'linewidth', 2)
end
xlim([0.1, min(nnodes(1), 50)+.5])
box off
ylabel('Mean sociability parameters', 'fontsize', 16)
xlabel('Index of node (sorted by dec. degree)', 'fontsize', 16)
legend('95% credible intervals', 'True value', 'location', 'northeast')
legend boxoff
[~, ind] = sort(sum(G, 2)+sum(G, 1)', 'descend');

% Low degree nodes
figure('name','Credible intervals - low degree nodes'); hold on
for i=max(1,nnodes(1)-50+1):nnodes(1)
    plot([i, i],...
        quantile(log(objmcmc.samples(1).mean_w(ind(i), :, :)),[.025,.975]), 'r', ...
        'linewidth', 3);
    hold on
    plot(i, log(mean_w_true(ind(i))), 'xg', 'linewidth', 2)
end
xlim([max(1,nnodes(1)-50+1)-.5, nnodes(1)+.5])
ylim([-12, -4])
box off
ylabel('Log mean sociability parameters', 'fontsize', 16)
xlabel('Index of node (sorted by dec. degree)', 'fontsize', 16)
legend('95% credible intervals', 'True value', 'location', 'southeast')
legend boxoff