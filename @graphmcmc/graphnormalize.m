function objmcmc = graphnormalize(objmcmc)

% graphnormalize normalizes feature importance of mcmc samples for solving 
% scale identifiability
% [objmcmc] = GRAPNORMALIZE(objmcmc)
%
% -------------------------------------------------------------------------
% INPUT
%   - objmcmc: an object of the class graphmcmc
%
% OUTPUT
%   - objmcmc: normalized object of the class graphmcmc
%
% See also GRAPHMCMC, GRAPHMCMC.GRAPHMCMC, GRAPHMCMCSAMPLES
% -------------------------------------------------------------------------

% Copyright (c) F. Caron (University of Oxford), A. Todeschini (Inria), and 
% X. Miscouridou (University of Oxford)
% caron@stats.ox.ac.uk
% adrien.todeschini@gmail.com
% xenia.miscouridou@spc.ox.ac.uk

nchains = size(objmcmc.samples, 2);
if size(objmcmc.samples, 1)==1 % unipartite graph
    for ch=1:numel(objmcmc.samples)
        sum_b = sum(objmcmc.samples(ch).Fparam.b , 1);

        objmcmc.samples(ch).logalpha = objmcmc.samples(ch).logalpha - objmcmc.samples(ch).sigma .* log(sum_b);
        objmcmc.samples(ch).alpha = exp(objmcmc.samples(ch).logalpha);
        objmcmc.samples(ch).tau = objmcmc.samples(ch).tau .* sum_b;
        objmcmc.samples(ch).Fparam.b = bsxfun(@rdivide, objmcmc.samples(ch).Fparam.b, sum_b);
    end
elseif size(objmcmc.samples, 1)==2 % bipartite graph
    for ch=1:nchains
        rate2 = bsxfun(@times, objmcmc.samples(2, ch).tau, objmcmc.samples(2, ch).Fparam.b);
        rate2 = permute(rate2, [2 1 3]);

        objmcmc.samples(1, ch).w = bsxfun(@rdivide, objmcmc.samples(1, ch).w, rate2);
        objmcmc.samples(1, ch).w_rem = bsxfun(@rdivide, objmcmc.samples(1, ch).w_rem, rate2);
        objmcmc.samples(2, ch).w = bsxfun(@times, objmcmc.samples(2, ch).w, rate2);
        objmcmc.samples(2, ch).w_rem = bsxfun(@times, objmcmc.samples(2, ch).w_rem, rate2);

        % sum b_k * sum b_k^prime
        prod_b = bsxfun(@times, objmcmc.samples(1, ch).Fparam.b, objmcmc.samples(2, ch).Fparam.b);
        sum_b = sum(prod_b ,1);

        objmcmc.samples(1, ch).logalpha = objmcmc.samples(1, ch).logalpha - objmcmc.samples(1, ch).sigma .* (log(objmcmc.samples(2, ch).tau) + log(sum_b));
        objmcmc.samples(1, ch).alpha = exp(objmcmc.samples(1, ch).logalpha);
        objmcmc.samples(1, ch).tau = objmcmc.samples(1, ch).tau .* objmcmc.samples(2, ch).tau .* sum_b;
        objmcmc.samples(1, ch).Fparam.b = bsxfun(@rdivide, prod_b, sum_b);

        objmcmc.samples(2, ch).logalpha = objmcmc.samples(2, ch).logalpha + objmcmc.samples(2, ch).sigma .* (log(objmcmc.samples(2, ch).tau));
        objmcmc.samples(2, ch).alpha = exp(objmcmc.samples(2, ch).logalpha);
        objmcmc.samples(2, ch).tau = ones(size(objmcmc.samples(2, ch).tau));
        objmcmc.samples(2, ch).Fparam.b = ones(size(objmcmc.samples(2, ch).Fparam.b));
    end
else
    error('not implemented yet');
end
