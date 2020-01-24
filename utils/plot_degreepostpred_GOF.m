function [quantile_freq, quantile_freq2, degsamp, sd_rowmean, sd_colmean, dyad_dep, triad_dep] = plot_degreepostpred_GOF(G, objmcmc, ndraws, T, rep, prefix, suffix, verbose)

% Function that  
% 1) plots the degree distribution from the posterior predictive distribution
%    either conditioning on the hyperparameters only (cond = false)or
%    using the estimated parameter values (cond = true)
% 2) calculates the goodness of fit statistics:
%      sd_rowmean: standard deviation of rowmeans
%      triad_dep: triadic dependence (or cluster soefficient)
% -------------------------------------------------------------
% INPUTS
%   - G: undirected binary adjacency matrix
%   - samples: MCMC output 
%
% OPTIONAL INPUT
%   - ndraws: number of draws from the posterior predictive
%   - cond: indicates whether the graphs will be sampled from the posterior
%   of hyperparameters (cond = false) or parameters (cond = true)
%   - rep:  output directory
%   - prefix: output filename prefix
%   - suffix: output filename suffix
%   - verbose: logical to indicate whether to print progress info
    
% OUTPUT 
%   - quantile_freq:  quantiles of the distribution  for matrix G
%   - quantile_freq2: quantiles of the distribution for matrix G'
%   - degsamp: sampled degree for each node of the network (for cond = false)
%   - sd_rowmean: standard deviation of the degree row means
%   - triad_dep: triadic dependence or cluster coefficient statistic
% -------------------------------------------------------------------------

% Copyright (c) F. Caron (University of Oxford), A. Todeschini (Inria), and 
% X. Miscouridou (University of Oxford)
% caron@stats.ox.ac.uk
% adrien.todeschini@gmail.com
% xenia.miscouridou@spc.ox.ac.uk
% September 2017
%--------------------------------------------------------------------------

fprintf('Sampling %d samples from predictive posterior distribution\n', ndraws)
fprintf('----------------------------------------------------------\n')

if nargin<3
    ndraws = 500;
end
if nargin<4
    T = 1e-6;
end
if nargin<5
    rep = './results/';
end
if nargin<6
    prefix = '';
end
if nargin<7
    suffix = '';
end
if nargin<8
    verbose = true;
end

freq1 = zeros(ndraws, length(0:16));
freq2 = zeros(ndraws, length(0:16));
htemp = figure('Visible', 'off');

p = objmcmc.prior.param(1).p;
F1name = objmcmc.prior.param(1).Fdist.name;
F2name = '';
if numel(objmcmc.prior.param)>1
    F2name = objmcmc.prior.param(2).Fdist.name;
end
typegraph = objmcmc.prior.typegraph;
observe_all = {objmcmc.prior(:).param.observe_all};

samples = combine(objmcmc.samples);
nsamples = size(samples, 3);
ind = floor(linspace(1, nsamples, ndraws));

N = size(G,1);
degsamp = zeros(ndraws,(N+50));
sd_rowmean = zeros(ndraws,1);
sd_colmean = sd_rowmean;
dyad_dep = sd_colmean;
triad_dep = dyad_dep;

if verbose
    fprintf('%s\n', repmat('-',1,35))
    fprintf('Start degree predictive posterior estimation for CGGP graphs: %d draws\n', ndraws)
end
tstart = tic;
Fdist.name='gamma';

for i=1:ndraws
        alpha = samples.alpha(:,:,ind(i));
        sigma = samples.sigma(:,:,ind(i));
        tau = samples.tau(:,:,ind(i));
        Fparam = struct('a', samples.Fparam.a(:,:,ind(i)), 'b', samples.Fparam.b(:,:,ind(i)));
        Fdist = struct('name', F1name, 'param', Fparam);
        gamma = samples(:).gamma(:,:,ind(i))'; 
        postmodel = graphmodel('CGGP', p, alpha, sigma, tau, Fdist, gamma, typegraph, observe_all{1});
    if verbose
        progress(i, tstart, ndraws, 35);
    end
    Gsamp = graphrnd(postmodel, T);
    ind = any(Gsamp);  %%%remove nodes with no connection
    Gsamp = Gsamp(ind, ind);
    
    [~, ~, freq1(i, :)] = plot_degree(Gsamp);
    [~, ~, freq2(i, :)] = plot_degree(Gsamp');
    N = size(Gsamp,1);
    degsamp(i,1:N)=sum(Gsamp,2);
    sd_rowmean(i) = std(mean(Gsamp,1));
%    sd_colmean(i) = std(mean(Gsamp,2));
%   dyad_dep(i) = corr(reshape(Gsamp,N^2,1),reshape(Gsamp',N^2,1));
    E = Gsamp'- mean(mean(Gsamp,2)); D = 1*(~isnan(E)); E(isnan(E)) = 0;
    triad_dep(i) = sum(diag(E*E*E)) / (sum(diag(D*D*D)) * std(reshape(Gsamp',N^2,1)).^3);

end



if verbose
    fprintf('|\n')
    fprintf('End degree predictive posterior estimation for CGGP graphs\n')
    fprintf('Computation time: %.1f hours\n', toc(tstart)/3600);
    fprintf('%s\n', repmat('-',1,35))
end

[~, centerbins1, ~] = plot_degree(G);
[~, centerbins2, ~] = plot_degree(G');
close(htemp);

% keyboard

quantile_freq = plot_figure(G, freq1, centerbins1, rep, prefix, ['1' suffix]);
quantile_freq2 = plot_figure(G', freq2, centerbins2, rep, prefix, ['2' suffix]);
end

%%
function progress(i, tstart, niter, s)
if i==5
    est_time = toc(tstart) * niter/i/3600;
    fprintf('Estimated end of computation: %s (%.1f hours)\n', datestr(now + toc(tstart) * (niter-i)/i/3600/24), est_time);
    fprintf('|%s|\n|', repmat('-',1,(min(niter,s-2))))
    fprintf(repmat('*', 1, floor(i/ceil(niter/min(niter,s-2)))))
end
if i>5 && mod(i, ceil(niter/min(niter,s)))==0
    fprintf('*')
end
end
%%
function quantile_freq = plot_figure(G, freq, centerbins, rep, prefix, suffix)

quantile_freq = quantile(freq, [.025, .975]);
plot_variance = @(x,lower,upper,color) fill([x,x(end:-1:1)],[upper,lower(end:-1:1)],color, 'EdgeColor', color);

figure; hold on
plot(centerbins, quantile_freq, 'color', [.8, .8, 1], 'linewidth', 2.5);
ind = quantile_freq(1,:)>0;
ha = plot_variance(centerbins(ind), quantile_freq(1,ind),quantile_freq(2,ind), [.8, .8, 1] );
set(gca,'XScale','log')
set(gca,'YScale','log')

hb = plot_degree(G);
set(hb, 'markersize', 10, 'marker', 'o',...
    'markeredgecolor', 'none', 'markerfacecolor', [1, .75, .75])

legend([ha, hb], {'$95\%$ posterior predictive', 'Data'})
legend boxoff
xlim([.8, 1e3])
box off
set(gca,'XMinorTick','on','YMinorTick','on')

if ~isempty(rep)
    savefigs(gcf,  [prefix 'degreepostpred' suffix], rep);
end

end
