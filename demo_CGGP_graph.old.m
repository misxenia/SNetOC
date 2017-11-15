%% BNPcommunity_graph package: demo_CGGP_graph
%
% This Matlab script shows how to sample an undirected graph with overlapping community
% structure from a generalized gamma process and gamma scores 
% and how to perform posterior inference.
%
% For downloading the package and information on installation, visit the
% <https://github.com/misxenia/BNPcommunity_graph BNPcommunity_graph webpage>.
% 
% Reference: 
% A. Todeschini, X. Miscouridou and F. Caron (2017) <https://arxiv.org/abs/1602.02114 Exchangeable Random Measures for Sparse and Modular Graphs with Overlapping Communities>. arXiv:1602.02114.
%
%% Authors: 
% * <http://www.stats.ox.ac.uk/~caron/ F. Caron>, University of Oxford
% * <http://adrien.tspace.fr/ A. Todeschini>, Inria
% * <http://csml.stats.ox.ac.uk/people/miscouridou/ X. Miscouridou>, University of Oxford
% 
% Tested on Matlab R2016a.
%
% Last Modified: 10/10/2017
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% General settings

close all
clear all

istest = false; % enable testing mode : quick run with small nb of iterations
% root = '.';
% % root = '~/demo_cggp/';
% 
% if istest
%     outpath = fullfile(root, 'results/demo_CGGP_graph/test/');
% else
%     outpath = fullfile(root, 'results/demo_CGGP_graph/', date);
%     if ~initp1
%         outpath = [outpath '_noinit'];
%     end
% end
% 
% if ~isdir(outpath)
%     mkdir(outpath)
% end

% Add paths
addpath ./GGP/ ./CGGP/ ./utils/ 

% Set the fontsize
set(0, 'DefaultAxesFontSize', 14)

% Set the seed
rng default
% rng shuffle

%% Simulation of a GGP graph

% Sample graph
%title = 'Simulated simple CGGP graph';
name = 'simugraph';
labels = {'Nodes', 'Nodes'};

p = 2; % Number of communities
alpha = 100; sigma = 0.1; tau = 1; % Parameters tuning the degree correction
Fdist.name = 'gamma'; Fdist.param.b = 1/p; Fdist.param.a = .2; % Parameters tuning the community structure
obj = graphmodel('CGGP', p, alpha, sigma, tau, Fdist);
[G, w_true, w_rem_true] = graphrnd(obj, 1e-6);

% Plot the adjacency matrix
figure; spy(G)
xlabel('nodes', 'fontsize', 16);ylabel('nodes', 'fontsize', 16);
title('Simulated simple CGGP graph');

% Plot the empirical degree distribution
deg = sum(G);
N = size(G, 1);
figure; h = plot_degree(G, 'o');
set(h, 'markersize', 6, 'color',  [.8, .3, .3],  'markerfacecolor', [.8, .3, .3]);
xlim([.9, max(deg)]);title('Empirical Degree Distribution');


%% MCMC inference for the CGGP graph
%
nchains = 2; % Number of MCMC chains
niterinit = 1000; % Number of iterations with GGP model for init
niter = 20000; % Number of MCMC iterations
nburn = floor(niter/2); % Number of burn-in iterations
nsamples = 250; thin = ceil((niter-nburn)/nsamples);  % Number of samples and thinning
verbose =  true;

objprior =  graphmodel('CGGP', p);
objmcmc = graphmcmc(objprior, niter, nburn, thin, nchains);
init = graphinit(objmcmc, G, niterinit);
objmcmc = graphmcmcsamples(objmcmc, G, verbose, init);


%% Summary statistics of the posterior

% % Normalize outputs
%  objmcmc = graphnormalize(objmcmc); % Standardize parameters for identifiability

% Get estimates and cost
 [estimates, C_st] = graphest(objmcmc);

 for ch=1:nchains
     objmcmc.samples(ch).mean_w_rem = mean(objmcmc.samples(ch).w_rem, 2);
     objmcmc.samples(ch).mean_w = mean(objmcmc.samples(ch).w, 2);
 end


%% Plots
%
prefix = sprintf('%s_%df_', name, p);
suffix = '';
outpath = './results/demo/';

% Plot cost
plot_cost(C_st,outpath, prefix, suffix);

% Assign max feature
[~, nodefeat] = max(estimates.w, [],2);

[~, ind_features] = sort(sum(estimates.w), 'descend');
%ind_features = 1:p;

% Plot traces and histograms
variables = {'logalpha', 'sigma', 'tau', 'Fparam.a', 'Fparam.b', 'mean_w_rem'};
namesvar = {'$\log \alpha$', '$\sigma$', '$\tau$', '$a$', '$b$', '$\overline{w}_{\ast}$'};
truevalues = {log(obj.param.alpha), obj.param.sigma, obj.param.tau, obj.param.Fdist.param.a, obj.param.Fdist.param.b, mean(w_rem_true) };

plot_trace(objmcmc.samples, objmcmc.settings, variables, namesvar, truevalues,outpath, prefix, suffix);

plot_hist(objmcmc.samples, variables, namesvar, truevalues, ind_features, [],outpath, prefix, suffix);


% Plot the graph by sorting the nodes by max feature
% plot_sortedgraph(G, nodefeat, nodefeat, ind_features, labels, rep, prefix, suffix, {'png'});
[~,indk] = max(w_true,[],2);
[~,ind] = sort(indk, 'descend');
figure
spy(G(ind, ind));title('Feature Sorted Graph');
xlabel('Nodes')
ylabel('Nodes')

% % % Correlation between the features
% figure;
% imagesc(corr(estimates.w(:,ind_features)))
% colormap('gray');title('Correlation between Features - A');
% 
% corr(estimates.w)
% for k=p:-1:1
%     temp = 1-exp(-estimates.w(:,ind_features(k))*estimates.w(:,ind_features(k))');
%     out(:,k) = temp(:);
% end
% corr(out); 
% figure; imagesc(corr(out));
% colormap('gray');title('Correlation between Features - B');
% 
% figure; plot(out(:,1),out(:,end), '.');
% xlabel('Feature 1 Component of Link Probability');
% ylabel('Feature 2 Component of Link Probability');
% 
% figure;imagesc(out'*out); colormap('gray');title('Correlation between Feature Components in Link Probabilities');


% Plot posterior predictive of degrees
ndraws = 100;
plot_degreepostpred(G, objmcmc, ndraws, 1e-6,outpath, prefix, suffix);
title('Posterior Predictive Degree Distribution');

%% Plot credible intervals for the mean weights over features
[~, ind] = sort(sum(G, 2)+sum(G, 1)', 'descend');

mean_w_true = mean(w_true,2);
mean_w_samples = cat(3, objmcmc.samples(1).mean_w, objmcmc.samples(2).mean_w);

figure; hold on
for i=1:min(N(1), 50)
    plot([i, i],quantile(mean_w_samples(ind(i), :, :),[.025,.975]), 'r','linewidth', 3);
    hold on
    plot(i, mean_w_true(ind(i)), 'xg', 'linewidth', 2)
end
xlim([0.1, min(N(1), 50)+.5])
box off
ylabel('Mean sociability parameters', 'fontsize', 16)
xlabel('Index of node (sorted by dec. degree)', 'fontsize', 16)
legend('95% credible intervals', 'True value', 'location', 'northeast')
legend boxoff
[~, ind] = sort(sum(G, 2)+sum(G, 1)', 'descend');

figure; hold on
for i=max(1,N(1)-50+1):N(1)
    plot([i, i],...
        quantile(log(mean_w_samples(ind(i), :, :)),[.025,.975]), 'r', ...
        'linewidth', 3);
    hold on
    plot(i, log(mean_w_true(ind(i))), 'xg', 'linewidth', 2)
end
xlim([max(1,N(1)-50+1)-.5, N(1)+.5])
ylim([-12, -4])
box off
ylabel('Log mean sociability parameters', 'fontsize', 16)
xlabel('Index of node (sorted by dec. degree)', 'fontsize', 16)
legend('95% credible intervals', 'True value', 'location', 'southeast')
legend boxoff

