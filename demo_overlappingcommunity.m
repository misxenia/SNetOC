%% Sparse Networks with Overlapping Communities (SNOC) package: demo_overlappingcommunity
%
% This Matlab script finds overlapping communities in a network of
% political blogs, using the wrapper function
% overlapping_community_detection.m. For a full analysis of this dataset,
% see the Matlab demo demo_polblogs.m.
%
% For downloading the package and information on installation, visit the
% <https://github.com/OxCSML-BayesNP/SNOC SNOC webpage>.
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
% Tested on Matlab R2017a. Requires the statistics toolbox.
%
% Last Modified: 01/2020
%%

close all
clear all

% Add path
addpath ./GGP/ ./CGGP/ ./utils/

%% Load network of political blogs
% Data can be downloaded from
% <https://www.cise.ufl.edu/research/sparse/mat/Newman/polblogs.mat here>.

load ./data/polblogs/polblogs.mat

titlenetwork = 'Political blogosphere Feb. 2005';
name = 'polblogs';
labels = {'Blogs', 'Blogs'};
groupfield = 'name'; % meta field displayed for group plot

% Transform the graph to obtain a simple graph
G = Problem.A | Problem.A'; % make undirected graph
G = logical(G-diag(diag(G))); % remove self edges (#3)

% Collect metadata
meta.name = cellstr(Problem.aux.nodename);
meta.source = cellstr(Problem.aux.nodesource);
meta.isright = logical(Problem.aux.nodevalue);
meta.degree = num2cell(full(sum(G,2)));
meta.groups = zeros(size(meta.isright));
meta.groups(~meta.isright) = 1;
meta.groups(meta.isright) = 2;
color_groups = [0 0 .8; .8 0 0];
label_groups = {'Left', 'Right'};
fn = fieldnames(meta);

% Remove nodes with no edge (#266)
ind = any(G);
G = G(ind, ind);
for i=1:length(fn)
    meta.(fn{i}) = meta.(fn{i})(ind);
end

% Plot adjacency matrix (sorted)
figure('name', 'Adjacency matrix (sorted by ground truth political leaning)')
spy(G);
xlabel(labels{2})
ylabel(labels{1})

%%

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

%% Find overlapping communities
%

p = 2; % Number of communities
niter = 20000; % Number of iterations
[degree_correction, community_affiliation, community_detection] = overlapping_community_detection(G, p, niter);


%% Some plots
%

% Identify each feature as liberal or conservative using ground truth 
% (This step would normally require a human interpretation of the features)
[~, ind_features] = sort(median(community_affiliation(meta.isright,:), 1)./median(community_affiliation, 1));
featnames = {'Liberal', 'Conservative'}; % Name of the interpreted features

% Print classification performance with ground truth
[confmat] = print_classif(fullfile('./', ['classif_' num2str(p) 'f.txt']), ...
    community_detection, meta.groups, ind_features, label_groups);


%%

% Plot levels of affiliation to each community for a subset of blogs
color = [0 0 .8; .8 0 0];
names = {'blogsforbush.com'
    'instapundit.com'
    'drudgereport.com'
    'tagorda.com'
    'danieldrezner.com/blog'
    'andrewsullivan.com'
    'iwantmycountryback.org'
    'democraticunderground.com'
    'wonkette.com'
    'washingtonmonthly.com'
    'dailykos.com'
    'atrios.blogspot.com'};
ind = zeros(size(names,1),1);
lab = cell(size(names,1), 1);
for i=1:size(names,1)
    ind(i) = find(strcmp(meta.name, names{i}) & strcmp(meta.name, names{i,1}));
    lab{i} = sprintf('%s (%s)', names{i}, label_groups{meta.groups(ind(i))}(1));
end
plot_nodesfeatures(community_affiliation, ind, ind_features, lab, featnames, color);

%%

% Plot the graph by clustering the nodes by community with maximum level of affiliation to see block structure
plot_sortedgraph(G, community_detection, community_detection, ind_features, labels);
