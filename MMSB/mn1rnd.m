function r = mn1rnd(p)

% MN1RND Samples a random logical vector r of size 1 from the multinomial distribution with n=1
% -------------------------------------------------------------------------
% INPUT
%   - p: probability

% OUTPUT
%   - r: logical indicator
% -------------------------------------------------------------------------


u = rand(size(p,1), 1); % uniform random var
r = indic(u, p); % logical indicator
end

function x = indic(u, p)
%INDIC Bin indicator vectors
cs = cumsum(p, 2);
under = bsxfun(@lt, u, cs);
cs(:,end) = [];
over = bsxfun(@ge, u, cs);
x = [under(:,1), under(:,2:end) & over];
end