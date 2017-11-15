function [G, s, pi] = mmsbrnd(n, alpha, W, rho)

% MMSBRND Samples an undirected graph from a mixed membership stochastic blockmodel(mmsb)
%
% -------------------------------------------------------------------------
% INPUTS
%   - n:     positive integer: number of nodes
%   - alpha: positive scalar; Dirichlet parameter 
%   - W:     matrix of size [p,p]; block matrix of between blocks link probabilities 
%   - rho:   scalar in (0,1); sparsity parameter 
%
% OUTPUT
%   - G: symmetric binary matrix 
% -------------------------------------------------------------------------
% Example
% n = 50; alpha = 0.5; W = rand(3,3); rho = 0.5;
% G = mmsbrnd(n, alpha, W, rho);

% Reference:
% E. M Airoldi, D. Blei, S. E Fienberg, and E. Xing. 
% Mixed membership stochastic blockmodels. The Journal of Machine Learning Research, 2008
%
% Copyright (c) F. Caron (University of Oxford), A. Todeschini (Inria), and 
% X. Miscouridou (University of Oxford)
% caron@stats.ox.ac.uk
% adrien.todeschini@gmail.com
% xenia.miscouridou@spc.ox.ac.uk
% September 2017
%--------------------------------------------------------------------------

p = size(W, 1);
pi = gamrnd(alpha*ones(n, p), 1);
pi = bsxfun(@times, pi, 1./sum(pi, 2));

s = NaN*zeros(n);
G = zeros(n);
for i=1:n
    for j=i+1:n
        s(i,j) = find(rand<cumsum(pi(i, :)), 1);
        s(j, i) = find(rand<cumsum(pi(j, :)), 1);
        G(i, j) = (rand< ((1-rho)*W(s(i,j), s(j, i))));        
    end
end
G = G + G';


