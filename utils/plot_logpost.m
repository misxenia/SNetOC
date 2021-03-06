function h_trace = plot_logpost(logp, iter, trueval, namevar, rep, prefix, suffix, leg, fontsize)
% plot_logpost plots the given log posterior probability distributions 
%
% INPUT
%   - logp: log posterior probabilities
%   - rep:  output directory
%   - prefix: output filename prefix
%   - suffix: output filename suffix
%
% OUTPUT
% h = vector of handles to the plotted lines.

[nsamples, nchains] = size(logp);

if nargin<2
    iter = 1:nsamples;
end

colour = get(gca,'ColorOrder');
nbcolour = size(colour, 1);

if nargin<3
    trueval = {};
end
if nargin<4
    namevar = 'Log-posterior';
end
if nargin<5
    rep = [];
end
if nargin<6
    prefix = '';
end
if nargin<7
    suffix = '';
end
if nargin<8
    leg = cell(nchains, 1);
    for k=1:nchains
        leg{k} = ['Chain ' num2str(k)];
    end
end
if nargin<9
    fontsize=22;
end


h_trace = figure;
for k=1:nchains
    colourk = colour(mod(k-1, nbcolour)+1,:);
    h(k,:) = plot(iter, logp(:,k), 'color', colourk);
    hold on
end

if ~isempty(trueval)
    h_true = plot([iter(1); iter(end)], [trueval(:)'; trueval(:)'], 'g--', 'linewidth', 3);
%     legend([h(:,1); h_true(1)], [leg; 'True'], 'fontsize', fontsize, 'location', 'Best','Interpreter','latex');
end
legend(h(:,1), leg, 'fontsize', fontsize, 'location', 'Best', 'Interpreter','latex');

legend boxoff
xlabel('MCMC iterations', 'fontsize', fontsize,'Interpreter','latex');
ylabel(namevar, 'fontsize', fontsize, 'interpreter', 'latex');
box off
axis tight
xlim([iter(1), iter(end)])
clear h;
if ~isempty(rep)
    savefigs(gcf, [prefix 'trace_logpost' suffix], rep);
end
