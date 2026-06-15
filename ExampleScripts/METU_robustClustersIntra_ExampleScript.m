% I AM WORKING ON THIS!!!!!!!!!!!!!!!!!!!!!

% error('This function was incorporated into the package METU as METU_scatterPlot_AnnotColor_grpLabel, for the UMAP article')


%% FUNCTIONS FOR THE MetU package
% Visualise list
deps = matlab.codetools.requiredFilesAndProducts('METU_interDatasetMatching_ExampleScript.m');
for i = 1:numel(deps)
    disp(deps{i})
end

% Find all files containing a function 
functionStringToFind = 'DMA_plotPercentOfDiagonal('; % DEFINE 
for i = 1:numel(deps)
    if string(deps{i}(end-3:end)) ~= ".mat"
        txt = fileread(deps{i});
        if contains(txt,functionStringToFind)
            fprintf('%s\n',deps{i});
        end
    end
end



%% This script needs the following packages:

% M2S 
% HELP
% umap 
% CanlabCore-master
% METU

%% Load AW1 old

cd('C:\Users\rjdossan\OneDrive - Imperial College London\Imperial projects\UMAPmetabolomics\Data');
load 220219_AirwaveDataImputed.mat;
Xtomatch.VarInfo = HELP_tableCellToString(Xtomatch.VarInfo);
X = Xtomatch;
X.UMAPres = UMAPres;

idxToKeep= (1:500)';
idxToKeep = ceil(linspace(1,size(X.Data,1)-1,500)')

X.Data = X.Data(idxToKeep,:);
X.SampleInfo = X.SampleInfo(idxToKeep,:);

X.Data = zscore(X.Data);

% save AW1_plasmaXCMS_500samples_imputed X
%UMAPres.reduction = plottingData;
%UMAPres.h_figUMAP = h_figUMAP;


%% Load AW1 new
% No missing values; It already contains UMAP coordinates
cd('C:\Users\rjdossan\OneDrive - Imperial College London\METU_UMAPtools\Data')
load('AW1_plasmaXCMS_500samples_imputed.mat')
X.VarInfo = HELP_tableCellToString(X.VarInfo);% make all columns string or numerical

%% START

% Plot UMAP with a diagonal
figure, plot(X.UMAPres.reduction(:,1), X.UMAPres.reduction(:,2),'.k')
METU_plotIsodensityLines(gcf, X.UMAPres.reduction); axis tight, grid on
METU_plotDiagonalORcircle(gcf); % close(h_ap)
METU_plotDiagonalORcircle(gcf, X.UMAPres.reduction(1500,:),0.075);

% Plot UMAP coloured by RefMet classes
annotationsForPlot = X.VarInfo.AbbreviatedAnnotation;
groupColorForPlot = X.VarInfo.RefMet_Sub_class;

figure, plot(X.UMAPres.reduction(:,1), X.UMAPres.reduction(:,2),'.k')
METU_plotIsodensityLines(gcf, X.UMAPres.reduction); axis tight, grid on
[h_figFinal,h_FigLegend] = METU_scatterPlot_AnnotColor_grpLabel(gcf,X.UMAPres.reduction,annotationsForPlot,groupColorForPlot,'jet',1,45,8)


%% Function 1: find all neighbours in multiple iterations

% Define settings for function 1 (and 2)
n_iterations = 10;
nSamples = 250;
louvain_gamma = 10;
onlyAdductsIsotopes=0;
UMAPthreshDist_betweenNeighbours=0.075;
nNeighUMAP = round(0.01*size(X.Data,2));
plotType = 1; % only one plot per iteration

[sortedEdges_all_str,TF_edges,all_neigh_edges, G_neigh, stats] = METU_findUMAPneighbours_part1_withCommunities...
    (X,onlyAdductsIsotopes,UMAPthreshDist_betweenNeighbours,nNeighUMAP, n_iterations, nSamples, louvain_gamma,plotType);


%% Display the results of function 1

disp('RESULTS')
disp('Number of edges with a specific frequency:')
disp(HELP_tabulateFormatted(TF_edges.Count))

% for each edge frequency, count the number of their nodes BEST
nNodes_atEdgeFrequency = []; nNodes_atEdgeFrequency_cum=[];
for c_i=1:n_iterations
    c=n_iterations-(c_i-1)% start by the last
    tempEdgesValue = TF_edges.Value(TF_edges.Count==c);
    tempNodes = str2double(split(tempEdgesValue,'_'));
    if ~isempty(tempNodes)
        uniqueNodesIdx = unique([tempNodes(:,1);tempNodes(:,2)]);
        nNodes_atEdgeFrequency(c,1) = length(uniqueNodesIdx);
        if c_i==1
            nodes_atEdgeFrequency_cum = uniqueNodesIdx;
            nNodes_atEdgeFrequency_cum(c,1) = nNodes_atEdgeFrequency(c,1);
        else        
            nodes_atEdgeFrequency_cum = unique([nodes_atEdgeFrequency_cum;uniqueNodesIdx]);
            nNodes_atEdgeFrequency_cum(c,1) = length(nodes_atEdgeFrequency_cum);
        end
    else
        nNodes_atEdgeFrequency(c,1) = 0;
        %nodes_atEdgeFrequency_cum = 0;
        nNodes_atEdgeFrequency_cum(c,1) = 0;
    end
end

% Frequency of edges and number of nodes for that edge frequency - not cumulative
M2S_figureH(0.4,0.7), subplot(2,1,1), [histogramInfo] = histogram(TF_edges.Count,'BinMethod','integers')
axis tight, grid on, title('Number of edges with x frequency')
text((1:length(histogramInfo.Values))', histogramInfo.Values', string(histogramInfo.Values'),...
    'fontsize',8,'VerticalAlignment','bottom','HorizontalAlignment','center')
xlabel('number of iterations'); ylabel('frequency of x edges')
subplot(2,1,2), bar(nNodes_atEdgeFrequency)
text((1:length(nNodes_atEdgeFrequency))', nNodes_atEdgeFrequency, string(nNodes_atEdgeFrequency'),...
    'fontsize',8,'VerticalAlignment','bottom','HorizontalAlignment','center')
axis tight, grid on, title('Number of nodes for edges with x frequency')
xlabel('number of iterations'); ylabel('number of features')


% Number of iterations an edge exists, e.g. edge 1_6, or edge 1_4
figure, bar(TF_edges.Count)
xlabel('number of edges'); ylabel('number of iterations in which the edge is reproducible')

%% Function 2: decide which features are reproducible neighbours

nNeighThreshold = 10 % Number of times features have to be neighbours
plotType = 0

% For Airwave1 using a previous UMAP

% nNeighUMAP_or_UMAP = round(0.01*size(X.Data,2)); % Use number of neighbours to define a UMAP, otherwise if a UMAP exists, use it
nNeighUMAP_or_UMAP = UMAPres;% Use UMAP if it already exists, otherwise use number of neighbours to define a UMAP

% NOTE: this will show all clusters, even small ones
[R] = METU_findUMAPneighbours_part2...
    (sortedEdges_all_str, TF_edges,nNeighThreshold,X, nNeighUMAP_or_UMAP, louvain_gamma, onlyAdductsIsotopes, UMAPthreshDist_betweenNeighbours, plotType) 

% Put ALL the edges on the plot
% METU_putEdgesOnPlot(gcf,X.UMAPres.reduction,R.FINAL_edges_idx(:,1),R.FINAL_edges_idx(:,2),'k-');

% Get the edges within threshold (number of times they are found in iterations)
R_TFedges = TF_edges(R.FINAL_edges_idx,:);

%% Plot results of function 2 - using CCs

% NOTE: the edges have now been filtered according to frequency, but their
% nodes can stil have created clusters that are "too small". 

% Define the minimum number of features per cluster
minNumberOfFeaturesPerCC = 10 % DEFINE THIS!!

% Plot the clusters within this threshold
figure
T_FINAL_CC = HELP_tabulateFormatted(R.FINAL_selectedNodes_CC)
bar(T_FINAL_CC.Count)
xlabel('Clusters sorted by size'); ylabel('Number of features')
axis tight, grid on
title('number of members of each FINAL CC')
M2S_plotaxes('r-',[NaN,minNumberOfFeaturesPerCC])


%% Get the results with this defined minimum number of features per CC - BEST ***

T_FINAL_CC_selByCC = T_FINAL_CC(T_FINAL_CC.Count >= minNumberOfFeaturesPerCC,:);
selByCC_LocalIdx_all = M2S_find_idxInReference(string(R.FINAL_selectedNodes_CC), T_FINAL_CC_selByCC.Value);
selByCC_LocalIdx = find(~isnan(selByCC_LocalIdx_all));
FINAL_selectedNodes_CC_selByCC = R.FINAL_selectedNodes_CC(selByCC_LocalIdx);
FINAL_selectedNodes_idx_inVarInfo_selByCC = R.FINAL_selectedNodes_idx_inVarInfo(selByCC_LocalIdx);

% Create a final UMAP, coloured by CC
groupColorForPlot_FinalCC = repmat("",size(annotationsForPlot));
groupColorForPlot_FinalCC(FINAL_selectedNodes_idx_inVarInfo_selByCC) = FINAL_selectedNodes_CC_selByCC;
FinalUMAPfig = figure;
plot(X.UMAPres.reduction(:,1), X.UMAPres.reduction(:,2),'.k')
METU_plotIsodensityLines(gcf, X.UMAPres.reduction); axis tight, grid on
% Plot the robust clusters
METU_scatterPlot_AnnotColor_grpLabel(FinalUMAPfig,X.UMAPres.reduction,annotationsForPlot,groupColorForPlot_FinalCC,'jet',1,45,8);
title('FINAL coloured by cluster')

% Create a table with final results

classForTable = X.VarInfo.RefMet_Sub_class(FINAL_selectedNodes_idx_inVarInfo_selByCC)
T_FINAL_selByCC = table(annotationsForPlot(FINAL_selectedNodes_idx_inVarInfo_selByCC) , classForTable, FINAL_selectedNodes_CC_selByCC,FINAL_selectedNodes_idx_inVarInfo_selByCC,...
    'VariableNames',{'Annotations','Class','clusterNr','idx_inVarInfo'})
T_FINAL_selByCC = sortrows(T_FINAL_selByCC, {'clusterNr'})

nClustersFINAL = length(unique(T_FINAL_selByCC.clusterNr))

VarInfo_tableToSave = [T_FINAL_selByCC, X.VarInfo(T_FINAL_selByCC.idx_inVarInfo,:)]


