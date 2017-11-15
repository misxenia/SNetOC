function r = mn1rnd2(proba)

% MN1RND2 Samples a vector r of length n such that Pr(r(i)=k) = proba(i,k)
% -------------------------------------------------------------------------
% INPUT
%   - proba: matrix of size [n,p]; with probability vectors of size r

% OUTPUT
%   - r: a vector of size [1,n]
% -------------------------------------------------------------------------


out = mn1rnd(proba);
[r, ~] = find(out');

end