%% Sparse Networks with Overlapping Communities (SNetOC) package: demo_usairport
%
% This Matlab script performs posterior inference on a network of airports 
% to find latent overlapping communities, using a Bayesian
% nonparametric approach.
%
% For downloading the package and information on installation, visit the
% <https://github.com/OxCSML-BayesNP/SNetOC SNetOC webpage>.
% 
% Reference: 
%
% * A. Todeschini, X. Miscouridou and F. Caron (2017) <https://arxiv.org/abs/1602.02114 Exchangeable Random Measures for Sparse and Modular Graphs with Overlapping Communities>. arXiv:1602.02114.
%
% Authors: 
%
% * <http://adrien.tspace.fr/ A. Todeschini>, Inria
% * <http://csml.stats.ox.ac.uk/people/miscouridou/ X. Miscouridou>, University of Oxford
% * <http://www.stats.ox.ac.uk/~caron/ F. Caron>, University of Oxford
% 
% Tested on Matlab R2017a. Requires the Statistics toolbox.
%
% Last Modified: 01/2020
%%

%% General settings
%

clear
close all
tstart = clock; % Starting time

istest = true; % enable testing mode: quick run with smaller nb of iterations

%%
% In test mode, a smaller number of iterations is run. Although the sampler
% clearly has not converged yet, the method already recovers well
% interpretable latent communities. To reproduce the results of the paper, set this value to
% false. 

root = '.';
if istest
    outpath = fullfile(root, 'results', 'demo_usairport', 'test');
else
    outpath = fullfile(root, 'results', 'demo_usairport', date);
end

if ~isdir(outpath)
    mkdir(outpath);
end


% Add path
addpath ./GGP/ ./CGGP/ ./utils/

% Default fontsize
set(0, 'DefaultAxesFontSize', 14)

% Set the seed
rng default


%% Load Network of airports with a connection to a US airport
% 

load ./data/usairport/usairport

titlenetwork = 'US airport network in 2010';
name = 'usairport';
labels = {'Airports','Airports'};

G = G | G'; % make undirected graph
G = logical(G-diag(diag(G))); % remove self edges (#164)

meta.degree = num2cell(full(sum(G,2)));
fn = fieldnames(meta);

% Plot adjacency matrix
figure;
spy(G);
xlabel(labels{2})
ylabel(labels{1})


%% Posterior Inference using Markov chain Monte Carlo and point estimation
% Users needs to start the parallel pool by using the command *parpool* to run multiple chains in parallel.

% Define the parameters of the prior
p = 4; % Number of commmunities
objprior =  graphmodel('CGGP', p); % CGGP graph model with p communities

% Define parameters of the MCMC sampler
nchains = 3;
if istest
    niterinit = 1000;
    niter = 20000;
    nsamples = 100;
    ndraws = 100;
else
    niterinit = 10000;
    niter = 1e7;
    nsamples = 1000;
    ndraws = 500;
end   
nburn = floor(niter/2);
thin = ceil((niter-nburn)/nsamples); 
verbose = true;

% Create the graphMCMC object
objmcmc = graphmcmc(objprior, niter, 0, thin, nchains);

% Run initialisation
init = graphinit(objmcmc, G, niterinit);

%%

% Run MCMC sampler
objmcmc = graphmcmcsamples(objmcmc, G, verbose, init);

%%

% Print summary in text file
print_summary(['summary_' num2str(p) 'f.txt'], titlenetwork, G, niter, nburn, nchains, thin, p, outpath, tstart)

% Save workspace
save(fullfile(outpath, ['workspace_' num2str(p) 'f.mat']))

%% 

% Log posterior approximation
[lp_nonlat, lp_lat, ll_nonlat, ll_lat] = logpost_approx(objmcmc, G);

%% discard burnin
objmcmc_noburn = objmcmc;
objmcmc_noburn.samples = discard(objmcmc.samples, floor(nburn/objmcmc.settings.thin));
objmcmc_noburn.settings.nburn = nburn;

%%

% Point estimation of the model parameters
[estimates, C_st] = graphest(objmcmc_noburn);

%% Plots
%

prefix = sprintf('%s_%df_', name, p);
suffix = '';
%%

% Plot Log posterior approximation
iter = (1:size(lp_nonlat,1))*thin;
plot_logpost(lp_nonlat, iter, [], 'Log-posterior', outpath, prefix, '_nonlat');
plot_logpost(lp_lat, iter, [], 'Log-posterior', outpath, prefix, '_lat');

% Plot log-posterior autocorr
lp_nonlat_noburn = lp_nonlat(floor(nburn/niter*size(lp_nonlat, 1)):end, :);
lp_lat_noburn = lp_lat(floor(nburn/niter*size(lp_lat, 1)):end, :);
plot_autocorr_logpost(lp_nonlat_noburn, thin, 'Log-posterior', outpath, prefix, '_nonlat');
plot_autocorr_logpost(lp_lat_noburn, thin, 'Log-posterior', outpath, prefix, '_lat');

% Plot cost
if ~isempty(C_st)
    plot_cost(C_st, outpath, prefix, suffix);
end

% Assign max feature for community detection
[~, nodefeat] = max(estimates.w, [],2);

% Identify each feature/community as Hub/East/West/Alaska 
% (This step would normally require a human interpretation of the features)
code_airports = {'JFK', 'LAN', 'DEN', 'BET'};
for k=1:length(code_airports)
    [~, ind_features(k)] = max(estimates.w(strcmp(meta.code, code_airports{k}), :));
end
if length(unique(ind_features))~=4
    warning('Problem with the interpretation of features/communities');
    ind_features = 1:4;
    featnames = {'Feature 1', 'Feature 2', 'Feature 3', 'Feature 4'};
else
    featnames = {'Hub', 'East', 'West', 'Alaska'};
end

%%

% Plot traces and histograms
variables = {'logalpha2', 'sigma', 'Fparam.a', 'Fparam.b2', 'mean_w_rem'};
namesvar = {'$\log \tilde\alpha$', '$\sigma$', '$a$', '$\tilde b$', '$\overline{w}_{\ast}$'};
plot_trace(objmcmc.samples, objmcmc.settings, variables, namesvar, [], outpath, prefix, suffix);
plot_hist(objmcmc_noburn.samples, variables, namesvar, [], ind_features, [], outpath, prefix, suffix);

%%

% Plot the graph by sorting the nodes by max feature
plot_sortedgraph(G, nodefeat, nodefeat, ind_features, labels, outpath, prefix, suffix, {'png'});

if isfield(meta, 'groups')
    % Plots by groups right vs left
    plot_groups(estimates.w, meta.groups, meta.(groupfield), ind_features, label_groups, featnames, color_groups, outpath, prefix, suffix);
end


%%

% Show the proportion in each features for a few nodes
% https://www.mapcustomizer.com/
names = {
    'New York, NY'
%     'Washington, DC'
    'Miami, FL'
%     'Detroit, MI'
%     'Knoxville, TN'
%     'Atlanta, GA'
%     'Louisville, KY'
%     'Indianapolis, IN'
    'Raleigh/Durham, NC'
    'Nashville, TN'
%     'Chicago, IL'
%     'Fayetteville, NC'
    'Lansing, MI'
    'Louisville, KY'
%     'Memphis, TN'
%     'Cleveland, OH'
    'Minneapolis, MN'
%     'Charleston/Dunbar, WV'
%     'Baltimore, ML'
%     'Tallahassee, FL'
%     'Portland, ME'
%     'Flint, MI'
%     'Champaign/Urbana, IL'
%     'Oklahoma City, OK'
%     'Des Moines, IA'
%     'Houston, TX'
%     'Dallas, TX'
    'Denver, CO'
%     'Fort Wayne, IN'
%     'Tyler, TX'
%     'Salt Lake City, UT'
%     'Phoenix, AZ'
    'Los Angeles, CA'
    'Seattle, WA'
%     'San Francisco, CA'
%     'Fairbanks, AK'
    'Anchorage, AK'
    'Bethel, AK'
};
ind = zeros(size(names,1),1);
for i=1:size(names,1)
    I = find(strcmp(meta.city, names{i}));
    [~, imax] = max([meta.degree{I}]);
    ind(i) = I(imax);
end
[~, ind2] = sort(meta.lon(ind), 'descend');
ind = ind(ind2);
names = names(ind2);
color = hsv(p);
plot_nodesfeatures(estimates.w, ind, ind_features, names, featnames, color, outpath, prefix, suffix);

%%

% Show some of the nodes in each feature
fnames = {'degree', 'city'}; % meta fields displayed for features exploration
formats = {'#%d,', '%s.'}; % 
print_features( fullfile(outpath, ['features_' num2str(p) 'f.txt']), ...
    estimates.w, ind_features, featnames, meta, fnames, formats)
print_features( fullfile(outpath, ['featuresnorm_' num2str(p) 'f.txt']), ...
    bsxfun(@rdivide, estimates.w, sum(estimates.w,2)),...
    ind_features, featnames, meta, fnames, formats)

fnames = {'city'}; % meta fields displayed for features exploration
formats = {'%s.'}; % 
print_features( fullfile(outpath, ['features_' num2str(p) 'f_tex.txt']), ...
    estimates.w, ind_features, featnames, meta, fnames, formats)
print_features( fullfile(outpath, ['featuresnorm_' num2str(p) 'f_tex.txt']), ...
    bsxfun(@rdivide, estimates.w, sum(estimates.w,2)),...
    ind_features, featnames, meta, fnames, formats)


%%

% Plot posterior predictive of degrees
plot_degreepostpred(G, objmcmc_noburn, ndraws, 1e-6, outpath, prefix, suffix);
 