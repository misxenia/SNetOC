function h_trace = plot_autocorr_logpost(logp, thin, namevar, rep, prefix, suffix, leg, fontsize)
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
    thin = 1;
end

colour = get(gca,'ColorOrder');
nbcolour = size(colour, 1);

if nargin<3
    namevar = 'Log-posterior';
end
if nargin<4
    rep = [];
end
if nargin<5
    prefix = '';
end
if nargin<6
    suffix = '';
end
if nargin<7
    leg = cell(nchains, 1);
    for k=1:nchains
        leg{k} = ['Chain ' num2str(k)];
    end
end

if nargin<8
    fontsize=22;
end


h_trace = figure;
hold on
for k=1:nchains
    colourk = colour(mod(k-1, nbcolour)+1,:);
    
    [acf, lags, ~] = autocorr(logp(:,k), floor(nsamples/2));
    lags = lags * thin;
    
    h(k,:) = plot(lags, acf, 'color', colourk);
    hold on
end

legend(h(:,1), leg, 'fontsize', fontsize, 'location', 'Best', 'Interpreter','latex');

legend boxoff
xlabel('Lag', 'fontsize', fontsize,'Interpreter','latex');
ylabel([namevar, ' Autocorrelation'], 'fontsize', fontsize, 'interpreter', 'latex');
box off
axis tight
xlim([lags(1), lags(end)])
clear h;
if ~isempty(rep)
    savefigs(gcf, [prefix 'autocorr_logpost' suffix], rep);
end
