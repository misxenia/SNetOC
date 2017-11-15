function [sd_rowmean,sd_colmean,dyad_dep,triad_dep] = GOFstats(Y)

% GOFstats calculates the goodness of fit statistics for a given adjacency
% matrix Y
% -----------------------------------------------------
% INPUT:
%   - Y: adjacency matrix 
% 
% OUTPUTS:
%   - sd_rowmean: vector; standard deviation of the degree row means
%   - sd_colmean: vector; standard deviation of the degree column means
%     (sd_colmean = sd_rowmean if Y is symmetric)
%   - dyad_dep: vector; dyadic dependence in Y data 
%      (dyad_dep =1 if Y symmetric)
%   - triad_dep: vector; triadic dependence or cluster coefficient is a statistic
%     that depends on the proportion of stars in the graph in every
%     occurence of a triplet connection
% 
% Reference:
% Peter D. Hoff,(2015). amen: Additive and Multiplicative Effects Models for Networks and
% Relational Data, R package
% version 1.1.
% -------------------------------------------------------------------------

% Copyright (c) F. Caron (University of Oxford), A. Todeschini (Inria), and 
% X. Miscouridou (University of Oxford)
% caron@stats.ox.ac.uk
% adrien.todeschini@gmail.com
% xenia.miscouridou@spc.ox.ac.uk
% September 2017
%--------------------------------------------------------------------------


    dims=size(Y);
    sd_rowmean=std(mean(Y,2));
    sd_colmean=std(mean(Y,1));
    dyad_dep = corr(reshape(Y,dims(1)^2,1),reshape(Y',dims(1)^2,1));
    E=Y-mean(mean(Y)); D=1*(~isnan(E)); E(isnan(E))=0;
    triad_dep=sum(diag(E*E*E))/(sum(diag(D*D*D)) * std(reshape(Y,dims(1)^2,1)).^3);
end
