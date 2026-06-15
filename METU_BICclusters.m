function [bic,best_bic] = METU_BICclusters(all_nClusters1,X,f1,f,matchedFeatures,plotBICorNot)
% This function defines BIC for a defined number of clusters.
%
% NOTE: all_nClusters1 has to be defined, e.g.:
% minNclusters1 = 1;
% maxNClusters1 = 100;
% gapNClusters1 = 1; 
% all_nClusters1 = (minNclusters1:gapNClusters1:maxNClusters1)';

bic = [];
for nClusters_idx = 1:length(all_nClusters1)
    nClusters = all_nClusters1(nClusters_idx);
    RRR = METU_defineClustersInA_visualiseInB(X,f1,f,matchedFeatures,nClusters,0);
    p = size(X(f).UMAPres.reduction,2)*nClusters; % parameters = 2 coords × k centroids
    current_bic = log(length(matchedFeatures(f1,f).Xt_connIdx))*p + length(matchedFeatures(f1,f).Xt_connIdx)*log(sum(RRR.Dist_set2_real)/length(matchedFeatures(f1,f).Xt_connIdx)); 
    disp('iteration  nClusters   BIC')
    fprintf('    %d\t      %d\t  %.4e\n\n', nClusters_idx, nClusters, current_bic)
    bic(nClusters_idx,1) = current_bic;
end
% Find best number of clusters (minimum value of BIC)
[sorted_bic, sorted_bic_idx] = sort(bic);
best_bic.nClusters = all_nClusters1(sorted_bic_idx(1));
best_bic.nClusters_BIC = sorted_bic(1);

% Plot
ResVisualisation = [];
if plotBICorNot > 0
    % Plot BIC per number of clusters
    figure, plot(all_nClusters1,bic,'k.')
    hold on, plot(all_nClusters1(sorted_bic_idx(1)), sorted_bic(1),'or')
    axis tight, grid on
    text(all_nClusters1(sorted_bic_idx(1)), sorted_bic(1),string(all_nClusters1(sorted_bic_idx(1))),'horizontalalignment','center','verticalalignment','bottom','fontsize',12,'color','r')
    title('Bayesian Information Criterion with increased number of clusters')
    xlabel('Number of clusters'); ylabel('BIC')
end
% if plotBICorNot > 1
%     % Plot UMAP with N clusters
%     nClusters = all_nClusters1(sorted_bic_idx(1))
%     ResVisualisation = METU_defineClustersInA_visualiseInB(X,f1,f,matchedFeatures,best_bic.nClusters,1)
% end