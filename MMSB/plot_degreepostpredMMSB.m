function [quantile_freq, quantile_freq2, degsamp, sd_rowmean, triad_dep] = plot_degreepostpredMMSB(G, samples, ndraws, cond, rep, prefix, suffix, verbose)

% plot_degreepostpredMMSB   
%
% 1) plots the degree distribution from samples drawn from the posterior predictive distribution
%    either conditioning on the hyperparameters only (cond = false)or
%    using the estimated parameter values (cond = true)
% 2) calculates the goodness of fit (GOF) statistics: (see GOFstats.m)
% -------------------------------------------------------------
% INPUTS
%   - G: symmetric binary adjacency matrix
%   - samples: structure; MCMC output 
%
% OPTIONAL INPUT
%   - ndraws: number of draws from the posterior predictive
%   - cond: logical; indicates whether the graphs will be sampled from the posterior
%     of hyperparameters (cond = false) or parameters (cond = true)
%   - rep:  output directory
%   - prefix: output filename prefix
%   - suffix: output filename suffix
%   - verbose: logical; indicates whether to print progress
%    
% OUTPUT 
%   - quantile_freq:  quantiles of the of the values in data G for the cumulative probability distribution (for matrix G)
%   - quantile_freq2: quantiles of the of the values in data G for the cumulative probability distribution (for matrix G')
%   - degsamp: vector; sampled degree for each node of the network (for cond = false)
%   - sd_rowmean: vector; standard deviation of the degree row means
%   - triad_dep: vector; triadic dependence or cluster coefficient statistic
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
      cond=false;
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

chain_num=1;
  n=size(G,1);
 
  freq = zeros(ndraws, length(0:16));
  freq2 = zeros(ndraws, length(0:16));
  htemp = figure('Visible', 'off');
   
  degsamp=zeros(ndraws,n);
    NS = zeros(ndraws,1);

  if verbose
      fprintf('%s\n', repmat('-',1,35))
      fprintf('Start degree predictive posterior estimation for mmsb graphs: %d draws\n', ndraws)
  end
  tstart = tic;


  for i=1:ndraws
    fprintf('iteration %d \n', i)

    W_est = samples(chain_num).W(:,:,i);
    rho_est = samples(chain_num).rho(i);
    alpha_est = samples(chain_num).alpha(i);

    if ~cond    %conditionally on hyperparameters only
      [G_samp, ~, ~] = mmsbrnd(n, alpha_est, W_est, rho_est);
    else        %conditionally on parameters
      s_est = samples{chain_num}.s(:,:,i);
      G_samp = mmsbrnd2(s_est, W_est, rho_est);
    end

    if verbose
   %   progress(i, tstart, ndraws, 35);
    end

    ind = any(G_samp);
    G_samp = G_samp(ind, ind);

    [~, ~, freq(i, :)] = plot_degree(G_samp);
    [~, ~, freq2(i, :)] = plot_degree(G_samp');


    nsamp = size(G_samp,1);
    NS(i)=nsamp;
    degsamp(i,1:nsamp)= sum(G_samp);
    % GOF statistics for matrix G_samp
    sd_rowmean(i)=std(mean(G_samp,2));
    E = G_samp - mean(mean(G_samp)); D=1*(~isnan(E)); E(isnan(E))=0;
    triad_dep(i)=sum(diag(E*E*E))/(sum(diag(D*D*D)) * std(reshape(G_samp,nsamp^2,1)).^3);

  end

  if verbose
      fprintf('\n')
      fprintf('End degree predictive posterior estimation for mmsb graphs\n')
      fprintf('Computation time: %.1f hours\n', toc(tstart)/3600);
      fprintf('%s\n', repmat('-',1,35))
  end

  [~, centerbins1, ~] = plot_degree(G);
  [~, centerbins2, ~] = plot_degree(G');
  close(htemp);

  % keyboard
  quantile_freq = plot_figure(G, freq, centerbins1, rep, prefix, ['1' suffix]);
  quantile_freq2 = plot_figure(G', freq2, centerbins2, rep, prefix, ['2' suffix]);


end


%% SUBFUNCTIONS %%
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

  legend([ha, hb], {'95% posterior predictive', 'Data'})
  legend boxoff
  xlim([.8, 1e3])
  box off
  set(gca,'XMinorTick','on','YMinorTick','on')

  if ~isempty(rep)
      savefigs(gcf,  [prefix 'degree_post_pred' suffix], rep);
  end

end
