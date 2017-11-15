% Test script for the Matlab package SNetOC

clear
close all

fprintf('=======================================\n')
fprintf('    TEST Matlab Package SNetOC    \n')
fprintf('=======================================\n')

% Add path
addpath ./GGP/ ./CGGP/ ./utils/

% Sample CGGP graph
p = 2; % nb of communities
alpha = 20; sigma = .5; tau = 1; 
Fdist.name = 'gamma';
Fdist.param.a = .2;
Fdist.param.b = 1/p;
gamma = zeros(p,1);
obj = graphmodel('CGGP', p, alpha, sigma, tau, Fdist, gamma);
G = graphrnd(obj, 1e-6);
plot_degree(G, 'or', 2);

% MCMC sampler for CGGP graph
hyper_alpha = [0.01,0.01]; hyper_sigma = [0.01, 0.01]; hyper_tau = [0.01,0.01];
hyper_Fdist.name = 'gamma';
hyper_Fdist.param.a = [0.01, 0.01];
hyper_Fdist.param.b = [0.01, 0.01];
objprior =  graphmodel('CGGP', p, hyper_alpha, hyper_sigma, hyper_tau, hyper_Fdist, gamma);
niter = 100; nburn = 50; nadapt = niter/4; thin = 2; nchains = 1; verbose = true;
objmcmc = graphmcmc(objprior, niter, nburn, thin, nchains, nadapt);
objmcmc = graphmcmcsamples(objmcmc, G, verbose);

fprintf('=======================================\n')
fprintf('             TEST COMPLETED            \n')
fprintf('=======================================\n')