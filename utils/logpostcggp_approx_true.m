function [lp_nonlat, ll_nonlat] = logpostcggp_approx_true(G, objprior, w_true, alpha_true, sigma_true, tau_true, Fdist_true, gamma_true)

if size(objprior.param, 1)>1
    error('Only valid for unipartite graphs')
end

p = objprior.param(1).p;

if nargin<8
    gamma_true=zeros(p,1);
end

% create sample with true values
sample_true.alpha = alpha_true;
sample_true.sigma = sigma_true;
sample_true.tau = tau_true;
sample_true.gamma = gamma_true;
sample_true.Fparam = Fdist_true.param;
if size(sample_true.Fparam.a, 1)==1
    sample_true.Fparam.a = sample_true.Fparam.a*ones(p,1);
end
if size(sample_true.Fparam.b, 1)==1
    sample_true.Fparam.b = sample_true.Fparam.b*ones(p,1);
end
sample_true.w = w_true;

% likelihood
ll_nonlat = loglikcggp_nonlat(G, sample_true, objprior);

% Prior
lprior = logpriorcggp_approx(sample_true, objprior.param);

lp_nonlat = ll_nonlat + lprior;
