function [samples, stats] = mmsbmcmc(G, objmodel, objmcmc, verbose)

% MCMC sampler for the mixed membership stochastic blockmodel (mmsb) 
% [samples, stats] = mmsbmcmc(G, objmodel, objmcmc, verbose)
%
% -------------------------------------------------------------------------
% INPUTS
%   - G: binary adjacency matrix
%   - objmodel: an object of class graphmodel with properties
%       - name: char with the full model name
%       - type: char with the abbreviated model name
%       - param: struct array of length 2, corresponding to each type of node
%           constaining model parameters with the following fields:
%           - n: positive integer. number of nodes
%           - p: positive integer. number of groups/blocks
%           - alpha: If empty, then alpha is estimated via an Improper Jeffrey prior. 
%             If scalar, the value of alpha is fixed. If vector of length 2,
%             parameters of the gamma prior.
%           - W: vector of length p, for the parameters of the beta prior
%             over the upper diagonal entries. (W symmetric) 
%           - rho: If empty, then rho is estimated via an Improper Jeffrey prior. If scalar,
%           - the initial rho for the MH update
%           - epsilon: If empty then not used, is scalar, then used as a
%             weight on the entries of W.
%       - typegraph: char ( default = 'undirected')
%   - objmcmc: an object of class graphmcmc with properties
%       - prior: graphmodel object
%       - settings: structure of mcmc parameters with the following fields:
%           - store_s: logical. If true, returns MCMC draws of s
%           - hyper.rw_std: standard deviation of the random walk
%           - hyper.nmh: number of MH iterations
%           - niter: number of MCMC iterations
%           - nburn: number of burn-in iterations
%           - thin: thinning of the MCMC output
%           - nchains: number of MCMC chains
% -------------------------------------------------------------------------
% OUTPUTS
%   - samples: MCMC samples 
%       - alpha: Dirichlet parameter
%       - W:     block matrix size pxp
%       - rho:   scalar in (0,1)
%       - pi:   matrix n x p with the p dimensional mixed membership probability vectors 
%         for each node
%       - s:  matrix n x n  with the vectors of emission indicator variable of the nodes'
%         group memberhsips  
%   - stats: 
%       -rateW:      acceptance rate in MH update for each entry in w
%       -rate_alpha: acceptance rate in MH update for the parameter alpha
%       -rate_rho:   (avg) acceptance rate in MH update for the parameter rho
% -------------------------------------------------------------------------

% Copyright (c) F. Caron (University of Oxford), A. Todeschini (Inria), and 
% X. Miscouridou (University of Oxford)
% caron@stats.ox.ac.uk
% adrien.todeschini@gmail.com
% xenia.miscouridou@spc.ox.ac.uk
% September 2017
%--------------------------------------------------------------------------

if nargin<4
    verbose=true;
end

samples = [];
stats = [];
mcmcparam = objmcmc.settings;
modelparam = objmodel.param;
niter = mcmcparam.niter;
nburn = mcmcparam.nburn;
thin = mcmcparam.thin;
nmh = mcmcparam.hyper.nmh;
rw_std = mcmcparam.hyper.rw_std;
n = size(G, 1);

ip = inputParser;
addRequired(ip, 'G', @(x) islogical(x) || isnumeric(x));
addRequired(ip, 'modelparam', @isstruct);
addRequired(ip, 'mcmcparam', @isstruct);
addOptional(ip, 'verbose', true, @islogical);
parse(ip, G, modelparam, mcmcparam, verbose);

if size(G, 2)~=n
    error('G must be a square matrix');
end
p = modelparam.p; % Number of latent communities
hyperW = modelparam.W;
if isempty(modelparam.alpha)
    % Estimate alpha - improper Jeffreys prior
    estimate_alpha = true;
    alpha = 1;
    hyperalpha=[0,0];
elseif(isnumeric(modelparam.alpha) && length(modelparam.alpha)==1)
    alpha = modelparam.alpha;
    estimate_alpha = false;
    rate_alpha = 0;    
elseif(isnumeric(modelparam.alpha) && length(modelparam.alpha)==2)
    % Estimate alpha - gamma prior
    estimate_alpha = true;
    alpha = 1;
    hyperalpha = modelparam.alpha;
else
    error('modelparam.alpha must be a scalar or empty')
end

%  Stuff for updating s
mask = triu(ones(n), 1);
[indup1, ~] = find(mask);
[inddown1, ~] = find(mask');
induplin = find(mask);
inddownlin = find(mask');

if isempty(modelparam.rho)
    % Estimate rho - Improper Jeffreys prior
    estimate_rho = true;
    rho = .5;
else
    estimate_rho = false;
    rho = modelparam.rho;
end

aW = hyperW(1);
bW = hyperW(2);


% Initialize W matrix
W = .1*ones(p)+.7*eye(p);% betarnd(hyperW(1), hyperW(2), p, p);
s = randi(p, n, n);
s = s + diag(NaN*ones(n,1));% Add NaN on the diagonal
% Compute summary statistics
[count, pp, nn] = get_counts(s, G, p);

% Initialize output
nsamples = floor((niter -nburn)/thin);
samples.W = zeros(p, p, nsamples);
samples.pi = zeros(n, p, nsamples);
samples.s = zeros(n, n, nsamples, 'int8');
samples.nn = zeros(p, p, nsamples);
samples.mm = zeros(p, p, nsamples);
samples.alpha = zeros(1,1, nsamples);
samples.rho = zeros(1,1, nsamples);
stats.rateW = zeros(p, p, nsamples);
stats.rate_alpha = zeros(1,1,nsamples);
stats.rate_rho = zeros(1,1,nsamples);

tic
% MCMC iterations
for iter=1:niter
    % Print current state
    printind = min(2000, round(niter/5)); % print output every 2000 or less
    if verbose && rem(iter, printind)==0
        fprintf('iter=%d/%d alpha=%.3f rho=%.3f ', iter, niter,alpha, rho);
        fprintf('W=');
        fprintf('%.3f ', W);fprintf('\n')
    end

    % Update W
    [W, rateW] = update_W(W, pp, nn, rho, aW, bW, nmh, rw_std);
    if ~isempty(modelparam.epsilon)
        W = diag(diag(W))  + modelparam.epsilon*ones(p) - modelparam.epsilon*eye(p);
    end

    % Update rho
    if estimate_rho
        [rho, rate_rho] = update_rho(rho, W, pp, nn, nmh, rw_std);
    else
    rate_rho =0;    
    end    

    % Update pi
    pi = gamrnd(alpha + count, 1);
    pi = bsxfun(@times, pi, 1./sum(pi, 2));

    % Update s
    s = update_s(s, G, rho, W, pi, indup1, inddown1, induplin, inddownlin);

    % Compute summary statistics
    [count, pp, nn] = get_counts(s, G, p);

    % Update alpha
    if estimate_alpha % We use improper Jeffreys prior here
        [alpha,rate_alpha] = update_alpha(alpha, pi, n, p, nmh, hyperalpha);
    end

       % Store output
    if (iter>nburn && rem((iter-nburn),thin)==0)
        ind = ((iter-nburn)/thin);
        samples.W(:,:, ind) = W;
        samples.pi(:,:, ind) = pi;
        samples.s(:, :, ind) = s;
        samples.alpha(1,1,ind) = alpha;
        samples.rho(1,1,ind) = rho;
        samples.pp(:, :, ind) = pp;
        samples.nn(:, :, ind) = nn;
        stats.rateW(:,:,ind) = rateW;
        stats.rate_alpha(1,1,ind) = rate_alpha;
        stats.rate_rho(1,1,ind) = rate_rho;
    end
    
end
end


%%%%%%%%%%%%%%%%%%%
% SUBFUNCTIONS
%%%%%%%%%%%%%%%%%%%

%%
function s = update_s(s, G, rho, W, pi, indup1, inddown1, induplin, inddownlin)

%n = size(G, 1);
p = size(W, 1);
s_tr = s'; 
% Update s(i,j), j>i
Gtemp = repmat(G(induplin), 1, p);
Wtemp = W(s_tr(induplin), :);
prob = pi(indup1, :) .* ( Gtemp .* (1-rho).* Wtemp...
    + (1-Gtemp) .* (1 - (1-rho).* Wtemp));
prob = bsxfun(@times, prob, 1./sum(prob, 2));
s(induplin) =  mn1rnd2(prob);

s_tr = s';

% Update s(j,i), j>i
Gtemp = repmat(G(inddownlin), 1, p);
Wtemp = W(s_tr(inddownlin), :);
prob = pi(inddown1, :) .* ( Gtemp .* (1-rho).* Wtemp...
    + (1-Gtemp) .* (1 - (1-rho).* Wtemp ));
prob = bsxfun(@times, prob, 1./sum(prob, 2));
s(inddownlin) =  mn1rnd2(prob);

end


function [W, rateW] = update_W(W, pp, nn, rho, aW, bW, nmh, rw_std)

p = size(W, 1);

if rho==0
    rateW = ones(p,p);
    W = betarnd(aW + pp, bW + nn);
    W = triu(W) + triu(W, 1)'; %keep only upper right and symmetrize
else
    for nitermh=1:nmh
        Wnew = logistic(logit(W) + rw_std*randn(p));
        llold = (aW + pp - 1) .*log(W) + nn.*log(1-(1-rho)*W) + (bW - 1)*log(1-W);
        llnew = (aW + pp - 1) .*log(Wnew) + nn.*log(1-(1-rho)*Wnew) + (bW - 1)*log(1-Wnew);
        qnew = -log(Wnew) - log(1-Wnew); % logitnormal proposal
        qold = -log(W) - log(1-W);
        logaccept_W = (llnew - llold - qnew + qold );
        ind = log(rand(p))< logaccept_W;
        W(ind) = Wnew(ind);
        W = triu(W) + triu(W, 1)'; %keep only upper right and symmetrize
        rateW = min(1, exp(logaccept_W));
    end
end
end

%%
function [rho,rate_rho] = update_rho(rho, W, pp, nn, nmh, rw_std)
for nitermh=1:nmh
    rhonew = logistic(logit(rho) + rw_std*randn);
    llold = pp .*log(1-rho) + nn.*log(1-(1-rho)*W);
    llnew = pp .*log(1-rhonew) + nn.*log(1-(1-rhonew)*W);
    qnew = -log(rhonew) - log(1-rhonew); % logitnormal proposal
    qold = -log(rho) - log(1-rho);
    logacc_rho = (llnew - llold - qnew + qold);
    if log(rand) < logacc_rho
        rho = rhonew;
    end
    rate_rho = min(1,mean(mean((exp(logacc_rho)))));
end

end

%%
function [alpha,rate_alpha] = update_alpha(alpha, pi, n, p, nmh, hyperalpha)
for nitermh=1:nmh
    alphanew = exp(log(alpha) + .02*randn);
    llold = (alpha-1)*sum(log(pi(:))) + n * gammaln(alpha*p) - n*p*gammaln(alpha);
    llnew = (alphanew-1)*sum(log(pi(:))) + n * gammaln(alphanew*p) - n*p*gammaln(alphanew);
    priornew = hyperalpha(1)*log(alphanew) - hyperalpha(2)*alphanew;
    priorold = hyperalpha(1)*log(alpha)    - hyperalpha(2)*alpha;
    logacc_alpha = (llnew-llold + priornew - priorold);
    if log(rand)< logacc_alpha
        alpha = alphanew;
    end
    rate_alpha = min(1,exp(logacc_alpha));

end
end

%%
function [count, pp, nn] = get_counts(s, G, p)
%     count = zeros(n, p);
    s_tr = s';
    ind = s.*(1-G)>0;
    temp = [s(ind), s_tr(ind)];
    nn = hist3(temp, {1:p, 1:p});
    nn = nn - diag(diag(nn))/2;

    ind2 = s.*G>0;
    temp2 = [s(ind2), s_tr(ind2)];
    pp = hist3(temp2, {1:p, 1:p});
    pp = pp - diag(diag(pp))/2;
    count = hist(s_tr, 1:p)';
end

function x = logit(p)
    x = log(p) - log(1-p);
end

function p  = logistic(x)
    p = 1./(1+exp(-x));
end
