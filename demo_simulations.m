%% Sparse Networks with Overlapping Communities (SNetOC) package: demo_simulations
%
% This Matlab script performs posterior inference on a (sparse) simulated graph with overlapping communities.
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

close all
clearvars

tstart = clock; % Starting time

istest = true; % enable testing mode: quick run with small nb of iterations

root = '.';
if istest
    outpath = fullfile(root, 'results', 'CGGP_simulations', 'test');
else
    outpath = fullfile(root, 'results', 'CGGP_simulations', date);
end

if ~isdir(outpath)
    mkdir(outpath)
end

% Add path
addpath ./GGP/ ./CGGP/ ./utils/ 

set(0, 'defaultAxesFontSize', 16)
set(0, 'defaultTextInterpreter', 'latex');
set(0, 'defaultAxesTickLabelInterpreter', 'latex');

% Set the seed
rng(22)

%% Sample a CGGP graph
%

titlenetwork = 'Simulated simple graph';
name = 'simugraph';
labels = {'Nodes', 'Nodes'};

% Set parameters of the CGGP graph
p = 2;
alpha_true = 200; sigma_true = 0.2; tau_true = 1; observe_all_true = false;
Fdist_true.name = 'gamma'; Fdist_true.param.b = 1/p; Fdist_true.param.a = 0.2;
gamma_true = zeros(p,1);
obj = graphmodel('CGGP', p, alpha_true, sigma_true, tau_true, Fdist_true, gamma_true, 'undirected', observe_all_true);
% Sample a CGGP graph
[G, w_true, w_rem_true] = graphrnd(obj, 1e-9);

nnodes = size(G, 1);
nedges = nnz(G);
fprintf('CGGP graph with %d nodes and %d edges sampled\n', size(G, 1), nnz(triu(G)));

% shuffle nodes
indperm = randperm(size(w_true,1));
w_true = w_true(indperm,:);
G = G(indperm, indperm);


figure('name', 'adjacency matrix')
spy(G)
xlabel('Nodes')
ylabel('Nodes')
title('Adjacency Matrix')

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

% Plot degree distribution
figure('name', 'Empirical degree distribution')
hdeg = plot_degree(G);
set(hdeg, 'markersize', 10, 'marker', 'o','markeredgecolor', 'none', 'markerfacecolor', [1, .75, .75]);

%% Prior distribution
% 

objprior =  graphmodel('CGGP', p);

%% Posterior inference
%

% Parameters of the MCMC algorithm
if istest
    niterinit = 100;
    niter = 20000;
    nsamples = 100; % Nb of Monte Carlo samples to return
    ndraws = 100;
else
    niterinit = 10000;
    niter = 200000;
    nsamples = 500;
    ndraws = 500;
end   
nburn = floor(3*niter/4); nchains = 3;
thin = ceil((niter-nburn)/nsamples);
verbose = true;

% Create the graphMCMC object
objmcmc = graphmcmc(objprior, niter, 0, thin, nchains); 
% Note: nburn is set to zero here in order to store samples in the transient regime of the MCMC

% Run initialisation
init = graphinit(objmcmc, G, niterinit);

%%

T = 1e-3;

% Run MCMC sampler
objmcmc = graphmcmcsamples(objmcmc, G, verbose, init, 'T', T);

%%

% True log-posterior
if sigma_true>0
    objprior.param.observe_all = false;
end
[lp_nonlat_true, ll_nonlat_true] = logpostcggp_approx_true(G, objprior, ...
    w_true, alpha_true, sigma_true, tau_true, Fdist_true, gamma_true);

% compute log-posterior
[lp_nonlat, lp_lat, ll_nonlat, ll_lat] = logpost_approx(objmcmc, G);

%%

% discard burnin to compute estimates
objmcmc_noburn = objmcmc;
objmcmc_noburn.samples = discard(objmcmc.samples, floor(nburn/objmcmc.settings.thin));
objmcmc_noburn.settings.nburn = nburn;

% Get estimates and cost
[estimates, C_st] = graphest(objmcmc_noburn);

%%

% Print summary in text file
print_summary(['summary_' num2str(p) 'f.txt'], titlenetwork, G, niter, nburn, nchains, thin, p, outpath, tstart)

% Save workspace
objmcmc.stats = [];
save(fullfile(outpath, ['workspace_' num2str(p) 'f.mat']), '-v7.3')

%% Plots
%

prefix = sprintf('%s_%df_', name, p);
suffix = '';

% Plot cost
plot_cost(C_st, outpath, prefix, suffix);

%%

% Plot log-posterior
iter = (1:size(lp_nonlat,1))*thin;
plot_logpost(lp_nonlat, iter, lp_nonlat_true, 'Log-posterior', outpath, prefix, '_nonlat');
plot_logpost(lp_lat, iter, [], 'Log-posterior', outpath, prefix, '_lat');

% Plot log-posterior autocorr
lp_nonlat_noburn = lp_nonlat(floor(nburn/niter*size(lp_nonlat, 1)):end, :);
lp_lat_noburn = lp_lat(floor(nburn/niter*size(lp_lat, 1)):end, :);
plot_autocorr_logpost(lp_nonlat_noburn, thin, 'Log-posterior', outpath, prefix, '_nonlat');
plot_autocorr_logpost(lp_lat_noburn, thin, 'Log-posterior', outpath, prefix, '_lat');


%%

% Plot traces and histograms

% order features 
ind_features = 1:p;

if sigma_true<0
    varsigma1_true = -alpha_true.*tau_true.^sigma_true./sigma_true;
    varsigma2_true = -sigma_true.*Fdist_true.param.a./(tau_true.*Fdist_true.param.b);
    varsigma3_true = sigma_true.*Fdist_true.param.a.*(sigma_true-Fdist_true.param.a-1)/(tau_true.*Fdist_true.param.b).^2;
    mean_w_rem_true = mean(w_rem_true);
    
    variables = {'varsigma1', 'varsigma2', 'varsigma3', 'mean_w_rem'};
    namesvar = {'$\varsigma_1$', '$\varsigma_2$', '$\varsigma_3$', '$\overline{w}_{\ast}$'};
    trueval = {varsigma1_true, varsigma2_true, varsigma3_true, mean_w_rem_true};
    plot_trace(objmcmc.samples, objmcmc.settings, variables, namesvar, trueval, outpath, prefix, suffix);
    plot_hist(objmcmc_noburn.samples, variables, namesvar, trueval, ind_features, [], outpath, prefix, suffix);
else
    logalpha2_true = log(alpha_true)+sigma_true*log(tau_true);
    b2_true = tau_true*Fdist_true.param.b;
    mean_w_rem_true = mean(w_rem_true);

    variables = {'logalpha2', 'sigma', 'Fparam.a', 'Fparam.b2', 'mean_w_rem'};
    namesvar = {'$\log \tilde\alpha$', '$\sigma$', '$a$', '$\tilde b$', '$\overline{w}_{\ast}$'};
    trueval = {logalpha2_true, sigma_true, Fdist_true.param.a, b2_true, mean_w_rem_true};
    plot_trace(objmcmc.samples, objmcmc.settings, variables, namesvar, trueval, outpath, prefix, suffix);
    plot_hist(objmcmc_noburn.samples, variables, namesvar, trueval, ind_features, [], outpath, prefix, suffix);
end

%%
% Plot credible intervals for the weights
plot_mean_w_ci(G, objmcmc.samples(1).mean_w, w_true, outpath, prefix, suffix)

%%

% Plot posterior predictive of degrees
plot_degreepostpred(G, objmcmc_noburn, ndraws, 1e-6, outpath, prefix, suffix);
