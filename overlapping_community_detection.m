function [degree_correction, community_affiliation, community_detection] = overlapping_community_detection(G, p, niter)

% OVERLAPPING_COMMUNITY_DETECTION finds latent (overlapping) communities in a network 
% -------------------------------------------------------------------------
% INPUTS
%   - G: sparse logical adjacency matrix of size n x n
%   - p: positive integer. number of features
% Optional input:
%       - niter: number of MCMC iterations
%
% OUTPUTS
%   - degree_correction: vector of size n x 1 of positive reals corresponding to the degree correction
%   - community_affiliation: matrix of size n x p; entry (i,k) in (0,1) corresponds to
%                         the level of affiliation of node i to community k
%   - community_detection: vector of integers of size n; index of the
%               community for which the node has highest level of affiliation
% -------------------------------------------------------------------------
% EXAMPLE
% n= 5;
% G = sparse(logical([ones(n)-eye(n), zeros(n); zeros(n), ones(n)-eye(n)]));
% p = 2; niter = 2000;
% [degree_correction, community_affiliation, community_detection] = overlapping_community_detection(G, p, niter);

% Copyright (c) A. Todeschini (Inria), X. Miscouridou (University of Oxford)
% and F. Caron (University of Oxford)% 
% adrien.todeschini@gmail.com
% xenia.miscouridou@spc.ox.ac.uk
% caron@stats.ox.ac.uk
% November 2017
%--------------------------------------------------------------------------


if nargin < 3
    niter=100000;
end

niterinit = 2000;
nsamples = 500;   
nburn = floor(niter/2);
thin = ceil((niter-nburn)/nsamples); 
verbose = true;
nchains = 1;
% CGGP graph model with p communities
objprior =  graphmodel('CGGP', p); 

% Create the graphMCMC object
objmcmc = graphmcmc(objprior, niter, nburn, thin, nchains);

% Run initialisation
init = graphinit(objmcmc, G, niterinit);

% Run MCMC sampler
objmcmc = graphmcmcsamples(objmcmc, G, verbose, init);

% Point estimation of the model parameters
[estimates, ~] = graphest(objmcmc);

degree_correction = sum(estimates.w,2);
community_affiliation = bsxfun(@rdivide, estimates.w, degree_correction);

[~, community_detection] = max(estimates.w, [],2); % Assign each node to the feature with highest weight
end
