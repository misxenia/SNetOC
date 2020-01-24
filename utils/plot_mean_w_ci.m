function plot_mean_w_ci(G, samples_mean_w, w_true, rep, prefix, suffix, fontsize)
% Plot credible intervals for the weights

if nargin<7
    fontsize=22;
end

nnodes = size(G, 1);

[~, ind] = sort(sum(G, 2)+sum(G, 1)', 'descend');
mean_w_true = mean(w_true,2);

% High degree nodes
figure('name','Credible intervals - high degree nodes'); hold on
first = 1:min(nnodes(1), 50);
for i=first
    q = quantile(samples_mean_w(ind(i), :, :),[.025,.975]);
    h1 = plot([i, i], q, 'r','linewidth', 3, 'Marker', '.', 'Markersize', 8);
end
h2 = plot(first, mean_w_true(ind(first)), 'xg', 'linewidth', 2);
xlim([0.1, min(nnodes(1), 50)+.5])
axis tight
box off
ylabel('Mean sociability parameters', 'fontsize', fontsize, 'Interpreter','latex')
xlabel('Index of node (sorted by dec. degree)', 'fontsize', fontsize, 'Interpreter','latex')
% legend([h1, h2], {'$95\%$ credible intervals', 'True value'}, 'location', 'northeast', 'Interpreter','latex')
legend boxoff
savefigs(gcf,  [prefix 'w_highdeg' suffix], rep);

% Low degree nodes
figure('name','Credible intervals - low degree nodes'); hold on
last = max(1,nnodes(1)-50+1):nnodes(1);
for i=last
    q = quantile(log(samples_mean_w(ind(i), :, :)),[.025,.975]);
    h1 = plot([i, i], q, 'r', 'linewidth', 3, 'Marker', '.', 'Markersize', 8);
    plot(i, log(mean_w_true(ind(i))), 'xg', 'linewidth', 2)
end
h2 = plot(last, log(mean_w_true(ind(last))), 'xg', 'linewidth', 2);
xlim([max(1,nnodes(1)-50+1)-.5, nnodes(1)+.5])
% ylim([-12, -4])
axis tight
box off
ylabel('Log mean sociability parameters', 'fontsize', fontsize, 'Interpreter','latex')
xlabel('Index of node (sorted by dec. degree)', 'fontsize', fontsize, 'Interpreter','latex')
% legend([h1, h2], {'$95\%$ credible intervals', 'True value'}, 'location', 'southeast', 'Interpreter','latex')
legend boxoff
savefigs(gcf,  [prefix 'w_lowdeg' suffix], rep);
