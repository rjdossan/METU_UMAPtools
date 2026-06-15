%% MetU_defineClustersInA_visualiseInB
% previously called UMAP_defineClustersInA_visualiseInB_function
% Res_interDatasetClusters = MetU_defineClustersInA_visualiseInB(X,f1,matchedFeatures,nClusters)
% Show the clusters in set 2 as defined in set 1. It requires that UMAP is already done.
% 
% INPUT
% X - struct of type X(f) containing:
% datasetName; set; ; UMAPres.reduction; 
% X(f).set - contains Data; VarInfo; VarInfo.MZRT_str; SampleInfo; 

% f1 - index of the reference dataset in X (e.g., 1)
% f - index of the target dataset in X (e.g., 2)
% matchedFeatures - structure of type  matchedFeatures(f1,f).Xx_connIdx
% containing the matched indices of reference and target for each position of X as columns: Xr_connIdx ; Xt_connIdx
% nClusters - number of dataset1 clusters to plot (e.g., 20)
%
% plotOrNot - 0 No 1 Yes
% 
% Rui Pinto 2025


function Res_interDatasetClusters = MetU_defineClustersInA_visualiseInB(X,f1,f,matchedFeatures,nClusters,plotOrNot)

if nargin<6
    plotOrNot == 1;
end
% Initialise the output
Res_interDatasetClusters = [];


%% THIS MUST BE CALCULATED AND PART OF THE INPUT 
%{
%% Variable containing the matched features in both sets
% This case it is a full match between all datasets, sorted in same order.
% For each set there is a Xr_connIdx and a Xt_connIdx
matchedFeatures = struct;
for f=1:length(X)
    if f==f1
        disp("do nothing")
    else
    matchedFeatures(f1,f).Xr_connIdx = (1:size(X(f1).set.Data,2))';
    matchedFeatures(f1,f).Xt_connIdx = (1:size(X(f).set.Data,2))';
    end
end
%}


%% Calculate the clusters on **dataset 1**

opt.nClusters = nClusters;

% Define set1
opt.dataset1Name = X(f1).datasetName;
opt.reduction1 = X(f1).UMAPres.reduction;
opt.dataset1_MZRTstr = X(f1).set.VarInfo.MZRT_str;
%opt.h_figUMAP_1 = X(f1).UMAPres.h_figUMAP;

% for f=:size(X)
    % Define set2
    opt.dataset2Name = X(f).datasetName;
    opt.reduction2 = X(f).UMAPres.reduction;
    opt.dataset2_MZRTstr = X(f).set.VarInfo.MZRT_str;
    %opt.h_figUMAP_2 = X(f).UMAPres.h_figUMAP;
    
    % Find the MZRTstr in the other dataset
    % dataset2_MZRTstrInDataset1 = opt.dataset1_MZRTstr(matchedFeatures(f1,f).Xr_connIdx); 
    % nDataset2MZstr_existingInDataset1 = length(dataset2_MZRTstrInDataset1); 

    % Find the vars in dataset2 with equivalent in dataset1
    % matchedFeatures(f1,f).Xt_connIdx = matchedFeatures(f1,f).Xt_connIdx;
    % matchedFeatures(f1,f).Xr_connIdx = matchedFeatures(f1,f).Xr_connIdx;
    %matchedFeatures(f1,f).Xr_connIdx = matchedFeatures(f1,f).Xr_connIdx;

    % Define z-scored data using only common features
    
    mean1 = nanmean(opt.reduction1(matchedFeatures(f1,f).Xr_connIdx,:));
    std1 = nanstd(opt.reduction1(matchedFeatures(f1,f).Xr_connIdx,:));
    reduction1_zscored = opt.reduction1./repmat(std1,size(opt.reduction1,1),1);
    reduction1_zscored = reduction1_zscored-repmat(mean1,size(opt.reduction1,1),1);
    
    mean2 = nanmean(opt.reduction2(matchedFeatures(f1,f).Xt_connIdx,:));
    std2 = nanstd(opt.reduction2(matchedFeatures(f1,f).Xt_connIdx,:));
    reduction2_zscored = opt.reduction2./repmat(std2,size(opt.reduction2,1),1);
    reduction2_zscored = reduction2_zscored-repmat(mean2,size(opt.reduction2,1),1);

    %% Define groups in dataset1
    clusterRes = kmeans(zscore(opt.reduction1),opt.nClusters);
    clusterRes_inDataset2 = NaN(size(X(f).set.VarInfo,1),1);
    clusterRes_inDataset2(matchedFeatures(f1,f).Xt_connIdx) = clusterRes(matchedFeatures(f1,f).Xr_connIdx);
    %groupPerNclusters1(:,nClusters_idx) =clusterRes;
    % calculate average euclidean distance for features also in dataset2
    dataset1_current_distancesToGroupMedian = calculateDistancesToGroupAverageX(reduction1_zscored, string(clusterRes));
    dataset1_current_distancesToGroupMedian = dataset1_current_distancesToGroupMedian.^2;
    dataset1_current_distancesToGroupMedian_common12 = dataset1_current_distancesToGroupMedian(matchedFeatures(f1,f).Xr_connIdx);
    mean_sum_squaredDist_groupPerNclusters1 = (1/length(matchedFeatures(f1,f).Xr_connIdx))* ((nansum(dataset1_current_distancesToGroupMedian(matchedFeatures(f1,f).Xr_connIdx)))');
    
    % These plots work but were copied to the end
    %{
    plotUMAPorNot = 0
    if plotUMAPorNot
        % Plot the clusters of dataset1
        figUMAP1 = METU_plotDatasetByReferenceColour(X(f1).UMAPres.reduction, string(clusterRes), X(f1).UMAPres.h_figUMAP, string(clusterRes))
        title(X(f1).datasetName,'interpreter','none')
        METU_addClassMedianLabel(gcf,X(f1).UMAPres.reduction,string(clusterRes),'k',15)
        
        % Plot the clusters of dataset1 in the current dataset2
        groupColorForPlot_2 = repmat("",size(X(f).UMAPres.reduction,1),1);
        groupColorForPlot_2(matchedFeatures(f1,f).Xt_connIdx) = string(clusterRes(matchedFeatures(f1,f).Xr_connIdx));
        figUMAP2 = METU_plotDatasetByReferenceColour(X(f).UMAPres.reduction, groupColorForPlot_2, X(f).UMAPres.h_figUMAP, string(clusterRes))
        title(X(f).datasetName,'interpreter','none')
        METU_addClassMedianLabel(gcf,X(f).UMAPres.reduction,groupColorForPlot_2,'k',15)
    end
    % Plot the distances
    % figure, plot(dataset1_current_distancesToGroupMedian_common12,'.')
    % title('Dist to centre of N clusters dataset1')
    %}
    %% Plot a subplot in dataset2 per group of dataset1 ***BEST***
    % N_grps_inDataset1 = opt.nClusters
    %dataset1grps_i = clusterRes;
    uniqueG1 = unique(str2double(clusterRes));
    % all_D1idx = matchedFeatures(f1,f).Xr_connIdx(matchedFeatures(f1,f).Xt_connIdx);
    % all_D2idx = matchedFeatures(f1,f).Xt_connIdx;
    all_D1idx = matchedFeatures(f1,f).Xr_connIdx;
    all_D2idx = matchedFeatures(f1,f).Xt_connIdx;
    % Define if real dataset1 groups or randomised ones
    for iterationNr = 1:2
        if iterationNr == 1
            useReal1_or_Random0 = 1;
        else
            useReal1_or_Random0 = 0;
        end
        if useReal1_or_Random0 == 1 % Real
            % Use the real vector with groups in dataset1
            dataset1grps{f,iterationNr} = clusterRes;
        else % Random
            % Use random groups from the groups in dataset1
            rng(3)
            randidx = randperm(numel(clusterRes));
            dataset1grps_random = clusterRes(randidx);
            dataset1grps{f,iterationNr} = dataset1grps_random;
        end
        all_D1clusters = dataset1grps{f,iterationNr}(all_D1idx);
        %{
        % Plot a defined number of dataset1 groups in dataset2
        
        if iterationNr == 1
            clusterSeq = (string(greedyClusterSequence(opt.reduction1(:,1), opt.reduction1(:,2), dataset1grps{f,iterationNr})))';
        end
        
        [r,c] = subplotDim(length(clusterSeq));
        figure('Name',opt.dataset2Name)
        for g_idx = 1:length(clusterSeq)
            subplot(r,c,g_idx), plot(opt.reduction2(:,1), opt.reduction2(:,2),'.k','markersize',3); axis tight, grid on
            local_currentGroup_idx = find(string(all_D1clusters) == clusterSeq(g_idx));
            hold on, plot(opt.reduction2(all_D2idx(local_currentGroup_idx),1), opt.reduction2(all_D2idx(local_currentGroup_idx),2),'or','markersize',3);
            title(strcat("Cluster ",clusterSeq(g_idx)) );
        end  
        %}
        
        % Get the distances to average of cluster for each feature, per cohort
        Dist_set1{f,iterationNr} = [];
        Dist_set2{f,iterationNr} = [];
        all_local_currentGroup_idx{f,iterationNr} = [];
        all_grp1string{f,iterationNr} = [];
        for g_idx = 1:nClusters
            % local_currentGroup_idx = find(string(all_D1clusters)== clusterSeq(g_idx));
            local_currentGroup_idx = find(string(all_D1clusters)== string(g_idx));
            all_local_currentGroup_idx{f,iterationNr} = [all_local_currentGroup_idx{f,iterationNr}; local_currentGroup_idx];
            all_grp1string{f,iterationNr} = [all_grp1string{f,iterationNr}; repmat(string(g_idx),length(local_currentGroup_idx),1)];
        
            % Compute Euclidean distance of each point to that centroid USING MEDIAN:
            Dist_set1{f,iterationNr} = [Dist_set1{f,iterationNr}; sqrt( (reduction1_zscored(all_D1idx(local_currentGroup_idx),1) - median(reduction1_zscored(all_D1idx(local_currentGroup_idx),1))).^2 +...
                (reduction1_zscored(all_D1idx(local_currentGroup_idx),2) - median(reduction1_zscored(all_D1idx(local_currentGroup_idx),2))).^2 )];
            Dist_set2{f,iterationNr} = [Dist_set2{f,iterationNr}; sqrt( (reduction2_zscored(all_D2idx(local_currentGroup_idx),1) - median(reduction2_zscored(all_D2idx(local_currentGroup_idx),1))).^2 +...
                (reduction2_zscored(all_D2idx(local_currentGroup_idx),2) - median(reduction2_zscored(all_D2idx(local_currentGroup_idx),2))).^2 )];
        end
        % MAYBE THE FIGURES BELOW NEED TO BE HERE
    end

% end

percentGoodFeaturesWithinThresh = struct;
if plotOrNot
    %% Define a threshold for the distances, based on reference distances
    
    % Fit a gamma distribution ot the reference dataset distances
    % threshold_all = [];
    d = Dist_set1{f,1};
    pd = fitdist(d,'Gamma');
    threshold = icdf(pd,0.99);
    %threshold_all(f,1) = threshold;
    
    % Find number of matched features within cluster thresholds (good features)
    percentGoodFeaturesWithinThresh.ref = sum((Dist_set1{f,1} < threshold) / length(Dist_set1{f,1}));
    percentGoodFeaturesWithinThresh.target = sum((Dist_set2{f,1} < threshold) / length(Dist_set2{f,1}));
    percentGoodFeaturesWithinThresh.targetRandom = sum((Dist_set2{f,2} < threshold) / length(Dist_set2{f,2}));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PLOTS

if plotOrNot
%% UMAP plots with clusters

% create UMAP plots for both datasets
X(f1).UMAPres.h_figUMAP = figure; plot(X(f1).UMAPres.reduction(:,1), X(f1).UMAPres.reduction(:,2),'.k')
METU_plotIsodensityLines(gcf, X(f1).UMAPres.reduction); axis tight, grid on

X(f).UMAPres.h_figUMAP = figure; plot(X(f).UMAPres.reduction(:,1), X(f).UMAPres.reduction(:,2),'.k')
METU_plotIsodensityLines(gcf, X(f).UMAPres.reduction); axis tight, grid on


% Plot the clusters of dataset1
[figUMAP1, figUMAP1_FigLegend] = METU_plotDatasetByReferenceColour(X(f1).UMAPres.reduction, string(dataset1grps{f,1}), X(f1).UMAPres.h_figUMAP, string(dataset1grps{f,1}));
title(X(f1).datasetName,'interpreter','none')
% METU_addClassMedianLabel(h_f,reduction,classLabels,textColor,textSize,medianOrMean,visibility_onOff)
  METU_addClassMedianLabel(gcf,X(f1).UMAPres.reduction,string(dataset1grps{f,1}),'k',15,'median','on'); 
close(figUMAP1_FigLegend);
set(gcf,'Name',strcat("UMAP on reference dataset colored by its clusters"));
xlabel('UMAP dimension 1'); ylabel('UMAP dimension 2'); 

% Plot the clusters of dataset1 in the current dataset2
groupColorForPlot_2 = repmat("",size(X(f).UMAPres.reduction,1),1);
groupColorForPlot_2(matchedFeatures(f1,f).Xt_connIdx) = string(dataset1grps{f,1}(matchedFeatures(f1,f).Xr_connIdx));
[figUMAP2, figUMAP2_FigLegend] = METU_plotDatasetByReferenceColour(X(f).UMAPres.reduction, groupColorForPlot_2, X(f).UMAPres.h_figUMAP, string(dataset1grps{f,1}));
set(gcf,'Name',strcat("UMAP on target dataset colored by clusters in the reference dataset")); title(X(f).datasetName,'interpreter','none')
METU_addClassMedianLabel(gcf,X(f).UMAPres.reduction,groupColorForPlot_2,'k',15,'median','on'); 
close(figUMAP2_FigLegend);
xlabel('UMAP dimension 1'); ylabel('UMAP dimension 2'); 

% Plot the RANDOMISED clusters of dataset1 in the current dataset2
groupColorForPlot_2RAND = repmat("",size(X(f).UMAPres.reduction,1),1);
groupColorForPlot_2RAND(matchedFeatures(f1,f).Xt_connIdx) = string(dataset1grps{f,2}(matchedFeatures(f1,f).Xr_connIdx));
[figUMAP2rand,  figUMAP2rand_FigLegend] = METU_plotDatasetByReferenceColour(X(f).UMAPres.reduction, groupColorForPlot_2RAND, X(f).UMAPres.h_figUMAP, string(dataset1grps{f,2}));
set(gcf,'Name',strcat("UMAP on dataset ", X(f).datasetName," colored by RANDOMISED clusters in the reference dataset")); title(strcat(X(f).datasetName,"_RANDOMISED"),'interpreter','none')
METU_addClassMedianLabel(gcf,X(f).UMAPres.reduction,groupColorForPlot_2RAND,'k',15,'median','on'); 
close(figUMAP2rand_FigLegend);
xlabel('UMAP dimension 1'); ylabel('UMAP dimension 2'); 

clusterSeq = (string(METU_greedyClusterSequence(opt.reduction1(:,1), opt.reduction1(:,2), dataset1grps{f,1})))';

% Subplot a defined number of dataset1 groups in dataset1
for iterationNr = 1:2
    all_D1clusters = dataset1grps{f,iterationNr}(all_D1idx);
    if iterationNr == 1
        figure('Name',strcat("UMAP on reference dataset ", X(f1).datasetName," with subplots colored by clusters in the reference dataset"));
    else
        figure('Name',strcat("UMAP on reference dataset ", X(f1).datasetName," with subplots colored by RANDOMISED clusters in the reference dataset"));
    end
    [r,c] = M2S_subplotDim(length(clusterSeq));
    for g_idx = 1:length(clusterSeq)
        subplot(r,c,g_idx), plot(X(f1).UMAPres.reduction(:,1), X(f1).UMAPres.reduction(:,2),'.k','markersize',3); axis tight, grid on
        local_currentGroup_idx = find(string(all_D1clusters) == clusterSeq(g_idx));
        hold on, plot(X(f1).UMAPres.reduction(all_D1idx(local_currentGroup_idx),1), X(f1).UMAPres.reduction(all_D1idx(local_currentGroup_idx),2),'or','markersize',3);
        title(strcat("RefCluster ",clusterSeq(g_idx)) );
    end  
end

% Subplot a defined number of dataset1 groups in dataset2
for iterationNr = 1:2
    all_D1clusters = dataset1grps{f,iterationNr}(all_D1idx);
    if iterationNr == 1
        figure('Name',strcat("UMAP on target dataset ", X(f).datasetName," with subplots colored by clusters in the reference dataset"));
    else
        figure('Name',strcat("UMAP on target dataset ", X(f).datasetName," with subplots colored by RANDOMISED clusters in the reference dataset"));
    end
    [r,c] = M2S_subplotDim(length(clusterSeq));
    for g_idx = 1:length(clusterSeq)
        subplot(r,c,g_idx), plot(X(f).UMAPres.reduction(:,1), X(f).UMAPres.reduction(:,2),'.k','markersize',3); axis tight, grid on
        local_currentGroup_idx = find(string(all_D1clusters) == clusterSeq(g_idx));
        hold on, plot(X(f).UMAPres.reduction(all_D2idx(local_currentGroup_idx),1), X(f).UMAPres.reduction(all_D2idx(local_currentGroup_idx),2),'or','markersize',3);
        title(strcat("TargetCluster ",clusterSeq(g_idx)) );
    end  
end
%xlabel('UMAP dimension 1'); ylabel('UMAP dimension 2'); 

%% PLOTS of distances

%f = 2
% Plot the distances to dataset1 centroids in dataset1 *** BEST ***
% scatterPlot_AnnotColor(figure,[(1:length(all_grp1string))',Dist_set2{iterationNr}],F.Xtomatch.VarInfo.AbbreviatedAnnotation(matchedFeatures(f1,f).Xt_connIdx),all_grp1string,'lines');
distFig(f).datasetName1 = X(f1).datasetName;
ap_fig = figure;
[distFig(f).distDataset1, distFig(f).distDataset1_legend] = METU_plotDatasetByReferenceColour([(1:length(Dist_set1{f,1}))',Dist_set1{f,1}], all_grp1string{f,1}, ap_fig, string(dataset1grps{f,1}));
title(strcat("Euclidean distance to dataset1 centroids in dataset ",distFig(f1).datasetName1),'interpreter','none'); ylim1 = ylim; close(distFig(f).distDataset1_legend); close(ap_fig);
set(gcf, 'Name',strcat("Euclidean distance to reference data centroids in reference dataset ",distFig(f1).datasetName1));
xlabel('Common variable number'); ylabel('Euclidean distance to centroid'); 

% Plot the distances to dataset1 centroids in dataset2 *** BEST ***
ap_fig = figure;
distFig(f).datasetName2 = X(f).datasetName;
[distFig(f).distDataset2, distFig(f).distDataset2_legend]  = METU_plotDatasetByReferenceColour([(1:length(Dist_set2{f,1}))',Dist_set2{f,1}], all_grp1string{f,1}, ap_fig, string(dataset1grps{f,1}));
title(strcat("Euclidean distance to dataset1 centroids in dataset ",distFig(f).datasetName2),'interpreter','none'); 
ylim2=ylim; close(distFig(f).distDataset2_legend); close(ap_fig);
set(gcf, 'Name',strcat("Euclidean distance to reference data centroids in target dataset ",distFig(f).datasetName2));
xlabel('Common variable number'); ylabel('Euclidean distance to centroid'); 

% Plot the distances to dataset1 centroids in dataset1 after randomisation *** BEST ***
ap_fig = figure;
distFig(f).datasetName1RAND = strcat(string(X(f1).datasetName), "_RANDOMISED");
[distFig(f).distDataset1RAND, distFig(f).distDataset1RAND_legend]  = METU_plotDatasetByReferenceColour([(1:length(Dist_set1{f,2}))',Dist_set1{f,2}], all_grp1string{f,2}, ap_fig, string(dataset1grps{f,1}));
title(strcat("Euclidean distance to RANDOMISED dataset1 centroids in dataset ",distFig(f).datasetName1RAND),'interpreter','none'); 
ylim1RAND=ylim; close(distFig(f).distDataset1RAND_legend); close(ap_fig);
set(gcf, 'Name',strcat("Euclidean distance to RANDOMISED dataset1 centroids in dataset ",distFig(f).datasetName1RAND));
xlabel('Common variable number'); ylabel('Euclidean distance to centroid'); 

% Plot the distances to dataset1 centroids in dataset2 after randomisation *** BEST ***
ap_fig = figure;
distFig(f).datasetName2RAND = strcat(string(X(f).datasetName), "_RANDOMISED");
[distFig(f).distDataset2RAND, distFig(f).distDataset2RAND_legend] = METU_plotDatasetByReferenceColour([(1:length(Dist_set2{f,2}))',Dist_set2{f,2}], all_grp1string{f,2}, ap_fig, string(dataset1grps{f,1}))
title(strcat("Euclidean distance to dataset1 centroids in dataset ",distFig(f).datasetName2RAND),'interpreter','none'); 
ylim2RAND=ylim; close(distFig(f).distDataset2RAND_legend); close(ap_fig);
set(gcf, 'Name',strcat("Euclidean distance to dataset1 centroids in dataset ",distFig(f).datasetName2RAND));
xlabel('Common variable number'); ylabel('Euclidean distance to centroid'); 

% Homogenise the ylim
ylim_minMax = [min([ylim1(1); ylim2(1); ylim1RAND(1); ylim2RAND(1)]), max([ylim1(2); ylim2(2); ylim1RAND(2); ylim2RAND(2)])];
figure(distFig(f).distDataset1); ylim(ylim_minMax);
figure(distFig(f).distDataset2); ylim(ylim_minMax);
figure(distFig(f).distDataset1RAND); ylim(ylim_minMax);
figure(distFig(f).distDataset2RAND); ylim(ylim_minMax);


drawnow
end

%% Get the output
Res_interDatasetClusters.clustersInReference = clusterRes;
Res_interDatasetClusters.clustersInTarget = clusterRes_inDataset2;

Res_interDatasetClusters.Dist_set1_real = Dist_set1{f,1};
Res_interDatasetClusters.Dist_set1_random = Dist_set1{f,2};
Res_interDatasetClusters.Dist_set2_real = Dist_set2{f,1};
Res_interDatasetClusters.Dist_set2_random = Dist_set2{f,2};

Res_interDatasetClusters.percentGoodFeaturesWithinThresh = struct;
if plotOrNot
    Res_interDatasetClusters.percentGoodFeaturesWithinThresh = percentGoodFeaturesWithinThresh;
    Res_interDatasetClusters.percentGoodFeaturesWithinThresh.threshold = threshold;
end


if plotOrNot
    %% Kernel density plot

    % Calculate kernel density functions
    [f1real,xi1real] = ksdensity(Res_interDatasetClusters.Dist_set1_real);
    [f2real,xi2real] = ksdensity(Res_interDatasetClusters.Dist_set2_real);
    [f2random,xi2random] = ksdensity(Res_interDatasetClusters.Dist_set2_random);
    
    % Plot kernel density functions as areas BEST!!!
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

    M2S_plotaxes('-k',[Res_interDatasetClusters.percentGoodFeaturesWithinThresh.threshold,NaN]);

    legend('ref real', 'target real', 'target random','Location','northeast');
    
    axis tight, grid on
    xlabel('Euclidean distance to centroid'); ylabel('Normal kernel density');
    title('Kernel density for each distance');


    
    drawnow
end
end


%% HELP FUNCTIONS

function distances = calculateDistancesToGroupAverageX(reduction, groupDefinition)
    % Get unique groups
    uniqueGroups = unique(groupDefinition);
    numPoints = size(reduction, 1);
    distances = zeros(numPoints, 1);
    
    % Loop through each group to calculate the median and distances
    for i = 1:length(uniqueGroups)
        group = uniqueGroups{i};
        groupIdx = strcmp(groupDefinition, group);
        groupData = reduction(groupIdx, :);
        
        % Calculate the median of the group
        groupMean = mean(groupData, 1);
        
        % Calculate the Euclidean distance to the group median
        distances(groupIdx) = sqrt(sum((groupData - groupMean).^2, 2));
    end
    distances(groupDefinition == "") = NaN;
end


