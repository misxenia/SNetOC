function ll = loglikcggp_lat(G, sample, prior)

% loglikcggp_nonlat computes the loglikelihood (ll) of the parameters
% obtained from the MCMC posterior distribution
% given the observed graph G
% ---------------------------------------------
% INPUTS:
%   - G: observed binary adjacency matrix
%   - sample: structure; parameter sample from the MCMC posterior distribution
%   - prior: object of graphmodel class; (only needed for the typegraph property )
% ---------------------------------------------
% OUTPUT:
%   - ll: loglikelihood
% -------------------------------------------------------------------------

% Copyright (c) F. Caron (University of Oxford), A. Todeschini (Inria), and 
% X. Miscouridou (University of Oxford)
% caron@stats.ox.ac.uk
% adrien.todeschini@gmail.com
% xenia.miscouridou@spc.ox.ac.uk
% September 2017
%--------------------------------------------------------------------------


if length(sample)==2
    w1 = sample(1).w;
    w2 = sample(2).w;
    rate = w1*w2';
    ok = G>0;
    %ll =  sum(log(1-exp(-rate(ok)))) + sum(-rate(~ok)); fix me
else
    w = sample.w;
    m = sample.m;
    K = size(w, 1);
    if strcmp(prior.typegraph, 'issimple')
        rate = 2*(w*w').*(1-eye(K));
        
    else
        rate = (w*w').*(2-eye(K));
        
    end
    ok = G>0;
    %ll =  sum(log(1-exp(-rate(ok)))) + sum(-rate(~ok));
    ll =  sum(sum(gammaln(sample.n+1))) + sum(sum(m .* log(w))) - sum(sum(w,1).^2);
end