function G = mmsbrnd2(s, W, rho)

% MMSBRND2 Samples an undirected graph from a mixed membership stochastic blockmodel
% conditionally on the group memberships of the nodes.
%
% -------------------------------------------------------------------------
% INPUTS
%   - s:   matrix of size n x n; group membership indicators
%   - W:     matrix of size [p,p]; block matrix of between blocks link probabilities 
%   - rho:   scalar in (0,1); sparsity parameter 

% OUTPUT
%   - G: symmetric binary matrix

% Reference:
% E. M Airoldi, D. Blei, S. E Fienberg, and E. Xing. 
% Mixed membership stochastic blockmodels. The Journal of Machine Learning Research, 2008

% Copyright (c) F. Caron (University of Oxford), A. Todeschini (Inria), and 
% X. Miscouridou (University of Oxford)
% caron@stats.ox.ac.uk
% adrien.todeschini@gmail.com
% xenia.miscouridou@spc.ox.ac.uk
% September 2017
%--------------------------------------------------------------------------

n = size(s, 1);
p = size(W, 1);

for i=1:n
    s(i,i) = 1;
end
s_tr = s';
% keyboard
prob = (1-rho) * W(sub2ind([p,p], s(:), s_tr(:)));
G = rand(n^2, 1)<prob;
G = triu(reshape(G, n, n), 1);
G = G + G';


