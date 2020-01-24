function [lp_nonlat, lp_lat, ll_nonlat, ll_lat ] = logpostcggp_approx(G, sample, prior)

% logpostcggp_nonlat computes an approxaimtion to the log posterior p(w|G) (up to a constant)
% without using latent variables
% it is the sum of the log likelihood p(G|w) and the log prior p(w) (up to a constant)


% INPUTS:
%   - G: observed binary adjacency matrix
%   - sample: structure; sample from the MCMC posterior distribution
%   - prior: object of class graphmodel; 
% ---------------------------------------------
% OUTPUT:
%   - lp: logposterior
%   - ll: loglikelihood 


% Approximation with and without latent variables for loglikelihood p(G|w)
ll_nonlat = loglikcggp_nonlat(G, sample, prior);
ll_lat = loglikcggp_lat(G, sample, prior);

lp_lat = ll_lat;
lp_nonlat = ll_nonlat;
nsamples  = length(sample) ;

isparallel = false;
if getpoolsize()>0 % If parallel pool initiated
    isparallel = true; % Will run MCMC chains in parallel using parallel computing toolbox        
end

% Prior
if isparallel  
    fprintf('Sampling graphs in parallel \n ')       
    parfor (i = 1:nsamples, getpoolsize())
        fprintf('draw %d\n',i)
	    lprior = logpriorcggp_approx(sample(i), prior.param(i));
	    lp_nonlat = lp_nonlat + lprior; 
	    lp_lat = lp_lat +lprior;
	end
else

	for i=1:nsamples 
	    lprior = logpriorcggp_approx(sample(i), prior.param(i));
	    lp_nonlat = lp_nonlat + lprior; 
	    lp_lat = lp_lat +lprior;
	end

end


