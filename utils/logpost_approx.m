function [logpost_nonlat, logpost_lat, loglik_nonlat, loglik_lat] = logpost_approx(objmcmc, G)

% logpost_nonlat computes an approximation to the log posterior (up to a constant) 
% without using augmented/latent variables for all 
% samples of parameters from the MCMC posterior 
% given the observed adjacency matrix G
% ---------------------------------------------
% INPUTS:
%   - objmcmc: object of class graphmcmc
%   - G: observed binary adjacency matrix
% ---------------------------------------------
% OUTPUT:
%   - logpost: logposterior
%   - loglik: loglikelihood
% -------------------------------------------------------------------------

% Copyright (c) F. Caron (University of Oxford), A. Todeschini (Inria), and 
% X. Miscouridou (University of Oxford)
% caron@stats.ox.ac.uk
% adrien.todeschini@gmail.com
% xenia.miscouridou@spc.ox.ac.uk
% August 2018
%--------------------------------------------------------------------------


% Computes the logposterior (up to a constant)

poolsize = getpoolsize();

if poolsize
    nchains = size(objmcmc.samples, 2);
    nsamples = size(objmcmc.samples(1).w, 3);
    
    samples_all = arraystruct2structarray(combine(objmcmc.samples));
    
    logpost_nonlat = zeros(nsamples*nchains,1);
    logpost_lat = zeros(nsamples*nchains,1);
    loglik_nonlat = zeros(nsamples*nchains,1);
    loglik_lat = zeros(nsamples*nchains,1);

    objmcmc_prior = objmcmc.prior;
    
    parfor(i=1:nsamples*nchains, poolsize)
        [logpost_nonlat(i), logpost_lat(i), loglik_nonlat(i), loglik_lat(i)] = logpostcggp_approx(G, samples_all(:,:,i), objmcmc_prior);
    end
    
    logpost_nonlat = reshape(logpost_nonlat, [nsamples nchains]);
    logpost_lat = reshape(logpost_lat, [nsamples nchains]);
    
    loglik_nonlat = reshape(loglik_nonlat, [nsamples nchains]);
    loglik_lat = reshape(loglik_lat, [nsamples nchains]);
else    
    samples = arraystruct2structarray(objmcmc.samples);

    [~, nchains, nsamples] = size(samples);
    logpost_nonlat = zeros(nsamples,nchains);
    logpost_lat = zeros(nsamples,nchains);
    loglik_nonlat = zeros(nsamples,nchains);
    loglik_lat = zeros(nsamples,nchains);

    for ch=1:nchains
        for i=1:nsamples
            [logpost_nonlat(i,ch), logpost_lat(i,ch), loglik_nonlat(i,ch), loglik_lat(i,ch)] = logpostcggp_approx(G, samples(:,ch,i), objmcmc.prior);
        end
    end
end
