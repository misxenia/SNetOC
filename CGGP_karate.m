%% BNPcommunity_graph package: CGGP_karate
%
% This Matlab script performs posterior inference on a network of political
% blogs to find latent overlapping communities, using a Bayesian
% nonparametric approach.
%
% For downloading the package and information on installation, visit the
% <https://github.com/misxenia/BNPcommunity_graph BNPcommunity_graph webpage>.
% 
% *Reference*:
%
% * A. Todeschini, X. Miscouridou and F. Caron.
% <https://arxiv.org/abs/1602.02114 Exchangeable Random Measures for Sparse
% and Modular Graphs with Overlapping Communities>. arXiv:1602.02114, 2017 (First version: Jan. 2016).
%
% *Authors*: 
%
% * <http://adrien.tspace.fr/ A. Todeschini>, Inria
% * <http://csml.stats.ox.ac.uk/people/miscouridou/ X. Miscouridou>, University of Oxford
% * <http://www.stats.ox.ac.uk/~caron/ F. Caron>, University of Oxford
%
% 
% Tested on Matlab R2016a. Requires the statistics toolbox.
%
% Last Modified: 10/2017


%% General settings
%

clear
close all
tstart = clock; % Starting time
root='.';
outpath = fullfile(root, 'results/CGGP_karate', date);

if ~isdir(outpath)
    mkdir(outpath);
end

% Add path
addpath ./GGP/ ./CGGP/ ./utils/

set(0, 'DefaultAxesFontSize', 14)

% Set the seed
rng default

%% Load network of karate networks
% Data can be downloaded from
% <http://networkdata.ics.uci.edu/data/karate/ here>.

load ./data/karate/karate.mat

titlenetwork = 'Karate';
labels = {'athletes', 'athletes'};

% Remove nodes with no edge
ind = any(G);
G = G(ind, ind);
fn=fieldnames(meta);
for i=1:length(fn)
    meta.(fn{i}) = meta.(fn{i})(ind);
end


% Plot adjacency matrix (sorted)
figure('name', 'Adjacency matrix (sorted by ground truth political leaning)')
spy(G);
xlabel(labels{2})
ylabel(labels{1})

% Shuffle nodes: irrelevant due to exchangeability, just to check we do not cheat!
ind = randperm(size(G,1));
G = G(ind, ind);
for i=1:length(fn)
    meta.(fn{i}) = meta.(fn{i})(ind);
end

% Plot adjacency matrix (unsorted)
figure('name', 'Adjacency matrix (unsorted)')
spy(G);
xlabel(labels{2})
ylabel(labels{1})

p=3;
niter=10000;
[degree_corr, memberships] = overlapping_community_detection(G, p, niter);

