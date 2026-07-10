%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% This script needs the following packages:

% M2S 
% HELP_RPtools
% umap 
% METU_UMAPtools


% Load a structure with the datasets of interest and another with the 
% matched features

cd('C:\Users\rjdossan\OneDrive - Imperial College London\METU_UMAPtools\Data')
load AW2_AW3_100samples_plasmaXCMS_imputed

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Define the datasets in X to use

f1 = 1 % reference
f = 2 % target

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Create UMAP of dataset1 coloured by subclass

% Reference dataset
annotationForPlot_ref = X(f1).set.VarInfo.AbbreviatedAnnotation;
groupColorForPlot_ref = X(f1).set.VarInfo.SubClass;

figure; plot(X(f1).UMAPres.reduction(:,1), X(f1).UMAPres.reduction(:,2),'.k')
METU_plotIsodensityLines(gcf, X(f1).UMAPres.reduction); axis tight, grid on
METU_scatterPlot_AnnotColor_grpLabel(gcf,X(1).UMAPres.reduction,annotationForPlot_ref,groupColorForPlot_ref,'jet',1,45,8);
title('Reference dataset coloured by class')

% Target dataset
annotationsForPlot_target = repmat("",size(X(f).UMAPres.reduction,1),1);
annotationsForPlot_target(matchedFeatures(f1,f).Xt_connIdx) = annotationForPlot_ref(matchedFeatures(f1,f).Xr_connIdx);
groupColourForPlot_target = repmat("",size(X(f).UMAPres.reduction,1),1);
groupColourForPlot_target(matchedFeatures(f1,f).Xt_connIdx) = X(f1).set.VarInfo.SubClass(matchedFeatures(f1,f).Xr_connIdx);

figure; plot(X(f).UMAPres.reduction(:,1), X(f).UMAPres.reduction(:,2),'.k')
METU_plotIsodensityLines(gcf, X(f).UMAPres.reduction); axis tight, grid on
METU_scatterPlot_AnnotColor_grpLabel(gcf,X(f).UMAPres.reduction,annotationsForPlot_target,groupColourForPlot_target,'jet',1,45,8);
title('Target dataset coloured by class')


%% Calculate BIC

%%%%%%%%%%%%%%%%
% Define the number of clusters to calculate BIC:
minNclusters1 = 1;
maxNClusters1 = 100;
gapNClusters1 = 1; 

all_nClusters1 = (minNclusters1:gapNClusters1:maxNClusters1)';

%% BIC

[bic, best_bic] = METU_BICclusters(all_nClusters1,X,f1,f,matchedFeatures,1);

% Plot number of clusters, with best number
figure, plot(all_nClusters1,bic,'k.')
% Find best number of clusters
[sorted_bic, sorted_bic_idx] = sort(bic);
hold on, plot(all_nClusters1(sorted_bic_idx(1)), sorted_bic(1),'or')
axis tight, grid on; xlabel('Number of clusters'); ylabel('BIC')
text(all_nClusters1(sorted_bic_idx(1)), sorted_bic(1),string(all_nClusters1(sorted_bic_idx(1))),'horizontalalignment','center','verticalalignment','bottom','fontsize',12,'color','r')
title('Bayesian Information Criterion with increased number of clusters')

% Calculate inter-dataset clusters with N clusters
nClusters = all_nClusters1(sorted_bic_idx(1))
Res_interDatasetClusters = METU_defineClustersInA_visualiseInB(X,f1,f,matchedFeatures,nClusters,1);

disp('Percentages of good features with distances to cluster median within threshold')
disp(struct2table(Res_interDatasetClusters.percentGoodFeaturesWithinThresh))

% PLOTS

% Plot UMAP of target dataset with distantce to median of cluster
figure; plot(X(f).UMAPres.reduction(:,1), X(f).UMAPres.reduction(:,2),'.k')
METU_plotIsodensityLines(gcf, X(f).UMAPres.reduction); axis tight, grid on
% hold on, plot(X(2).UMAPres.plottingData(:,1), X(2).UMAPres.plottingData(:,2),'.k')
scatter(X(f).UMAPres.reduction(matchedFeatures(f1,f).Xt_connIdx,1), X(f).UMAPres.reduction(matchedFeatures(f1,f).Xt_connIdx,2),40,Res_interDatasetClusters.Dist_set2_real,'filled')
colormap jet; colorbar; title('UMAP of target dataset coloured by distance to median of cluster')

% Calculate and plot kernel density functions
[f1real,xi1real] = ksdensity(Res_interDatasetClusters.Dist_set1_real);
[f2real,xi2real] = ksdensity(Res_interDatasetClusters.Dist_set2_real);
[f2random,xi2random] = ksdensity(Res_interDatasetClusters.Dist_set2_random);

% Plot kernel density functions of distances
figure, 
h_real1 = area(xi1real, f1real);
h_real1.FaceColor = 'r';
h_real1.FaceAlpha = 0.8; 
hold on

h_real2 = area(xi2real, f2real);
h_real2.FaceColor = 'b';
h_real2.FaceAlpha = 0.5; 

h_random2 = area(xi2random, f2random);
h_random2.FaceColor = 'c';
h_random2.FaceAlpha = 0.7; 
legend('ref real', 'target real', 'target random','Location','northeast')

axis tight, grid on
xlabel('Euclidean distance to centroid'); ylabel('Normal kernel density')
title('Kernel density for each distance')

M2S_plotaxes('-k',[Res_interDatasetClusters.percentGoodFeaturesWithinThresh.threshold,NaN]);


