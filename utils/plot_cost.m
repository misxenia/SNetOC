function h = plot_cost(C_st, rep, prefix, suffix, fontsize)

% plot_cost plots the cost of each MCMC sample 
%  the cost function is a permutation-invariant absolute loss on the weights
% see munkres.m
%
% INPUT
%   - C_st: matrix of cost values for MCMC output for every chain
%   - rep: output directory to save the figure
%   - prefix: output filename prefix
%   - suffix: output filename suffix
% OUTPUT
%
% - h: figure handle vector

if nargin<5
    fontsize=22;
end

h = figure;
plot(C_st)
xlabel('MCMC iterations', 'fontsize', fontsize, 'interpreter', 'latex');
ylabel('Cost', 'fontsize', fontsize, 'interpreter', 'latex');
nchains = size(C_st, 2);
leg = cell(nchains, 1);
for k=1:nchains
    leg{k} = ['Chain ' num2str(k)];
end
legend(leg, 'fontsize', fontsize, 'location', 'Best', 'Interpreter','latex');
legend boxoff
box off
axis tight
savefigs(gcf, [prefix 'cost' suffix], rep)
