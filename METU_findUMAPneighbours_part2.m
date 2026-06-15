%function [sortedEdges_sel,final_neigh_idx] = findUMAPneighbours_v2part2(sortedEdges_all_str, TF_edges,nNeighThreshold, X, nNeighUMAP, louvain_gamma, onlyAdductsIsotopes, UMAPthreshDist_betweenNeighbours) 

function [R] = METU_findUMAPneighbours_part2(sortedEdges_all_str, TF_edges,nNeighThreshold, X, nNeighUMAP_or_UMAP, louvain_gamma, onlyAdductsIsotopes, UMAPthreshDist_betweenNeighbours,plotType) 



%% DEFINITIONS
% INPUT
% sortedEdges_all_str: all the edges found in all iterations, as string
% with nodes e.g. 3_6. Comes from previous function (findUMAPneighbours_v2part1_withCommunities.m).
% TF_edges is a table counting the number of times each edge exists in the
% iterations calculated.
% nNeighThreshold = 5 is the minimum Count an edge needs to have on
% TF_edges to be accepted as reproducible.
% X contains X.VarInfo, X.SampleInfo, X.Data
% nNeighUMAP = 0.01*size(X.VarInfo,1): number of UMAP neighbours e.g. 1% of total of objects to cluster
% OUTPUT
% TF_edges_sel is a table with the selected edges (indices of VarInfo)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

%% FOR ONLY ADDUCTS AND ISOTOPES
if onlyAdductsIsotopes == 1
    RTthresh_neighbours_iterations =  0.25/60;
    corrThresh = 0.75;
else
%% FOR EVERYTHING
    RTthresh_neighbours_iterations =  1000;
    corrThresh = -1;
end

% Only calculate correlations if needed
if corrThresh>-1
    calculate_corr = true;
else
    calculate_corr = false;
end
% In case correlations are to be calculated: ****************
calculate_corr = false;
if calculate_corr
    allCorr = corr(X.Data','type','Pearson'); 
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% WORK ON THE FULL DATASET
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create a UMAP of the entire dataset
% Find neighbours and Louvain communities, delete edges when both their 
% nodes are not in the same Louvain community.
% Intersect those edges with the edges calculated in the previous function
% (above)

%% Part 1: Create a UMAP plot using the entire dataset
if isnumeric(nNeighUMAP_or_UMAP)% if a number of neighbours is supplied
    metricType = 'cosine';
    UMAPres = METU_UMAPregular(zscore(X.Data)', nNeighUMAP_or_UMAP,metricType);
else % if a UMAPres is supplied
    UMAPres = nNeighUMAP_or_UMAP;
end



% Normalise and collect UMAP distances in percentage of diagonal
reduction_temp = UMAPres.reduction - repmat(min(UMAPres.reduction),size(UMAPres.reduction,1),1);
reduction_between01 = reduction_temp ./ repmat(max(reduction_temp),size(reduction_temp,1),1);
temp_UMAPcoordPerc = reduction_between01;
clear reduction_between01;

%% Part 2: find neighbours in the entire dataset, plot UMAP with edges

neighbours_idx_all = []; %cell(n_iterations,1);
% all_neigh_edges_twoColumns_all = cell(n_iterations,1);
all_neigh_edges = []; % cell(n_iterations,1);
% all_neigh_edges_Louvain= cell(n_iterations,1);

%temp_UMAPcoordPerc_iterations = UMAPcoordPerc_iterations{iterationNr,1};
neighbours_idx = cell(size(temp_UMAPcoordPerc,1),1);
% all_neigh_edges_twoColumns = [];
temp_all_neigh_edges = [];
for varNr = 1:size(temp_UMAPcoordPerc,1)
    %[varNr, size(X.Data,2)]
    tempDiff1 = temp_UMAPcoordPerc(:,1) - temp_UMAPcoordPerc(varNr,1);
    tempDiff2 = temp_UMAPcoordPerc(:,2) - temp_UMAPcoordPerc(varNr,2);
    tempDist = sqrt(tempDiff1.^2 + tempDiff2.^2);
    tempDist(varNr) = 10000; % to not match itself
    if calculate_corr
        current_neighbours_idx = find(tempDist < UMAPthreshDist_betweenNeighbours & allCorr(:,varNr)>= corrThresh); % This could be calculated using correlations
    else
        current_neighbours_idx = find(tempDist < UMAPthreshDist_betweenNeighbours);
    end
    local_acceptedIdx = find(abs(X.VarInfo.rtmed(current_neighbours_idx) - X.VarInfo.rtmed(varNr)) < RTthresh_neighbours_iterations);
    neighbours_idx{varNr} = current_neighbours_idx(local_acceptedIdx);
    %all_neigh_edges = [all_neigh_edges;[repmat(varNr,length(current_neighbours_idx(local_acceptedIdx)),1),current_neighbours_idx(local_acceptedIdx)]];
    temp_all_neigh_edges = [temp_all_neigh_edges;[repmat(varNr,length(current_neighbours_idx(local_acceptedIdx)),1),current_neighbours_idx(local_acceptedIdx)]];
end
all_neigh_edges = temp_all_neigh_edges(temp_all_neigh_edges(:,1) < temp_all_neigh_edges(:,2),:);
% all_neigh_edges_twoColumns = [all_neigh_edges_twoColumns; all_neigh_edges{iterationNr,1}];

% Plot all the edges between neighbours in the full dataset
if plotType > 1
    UMAPfig_toEditFullDataset = figure; 
    plot(UMAPres.reduction(:,1), UMAPres.reduction(:,2),'.k')
    METU_plotIsodensityLines(gcf, UMAPres.reduction); axis tight, grid on
    % UMAPfig_toEditFullDataset = copyfig(UMAPres.h_figUMAP);
    %figure(UMAPfig_toEditFullDataset); 
    hold on
    METU_putEdgesOnPlot(UMAPfig_toEditFullDataset,UMAPres.reduction,all_neigh_edges(:,1),all_neigh_edges(:,2),'r-');
end

neighbours_idx_all = neighbours_idx;  


% Plot all the edges distances THIS WORKS!
% distancesBetweenNodesInEdges=[];for a=1:length(all_neigh_edges_twoColumns_all{iterationNr,1} (:,1)); distancesBetweenNodesInEdges(a,1) = pdist([temp_UMAPcoordPerc_iterations(all_neigh_edges_twoColumns_all{iterationNr,1} (a,1));temp_UMAPcoordPerc_iterations(all_neigh_edges_twoColumns(a,2))],'euclidean');end
% figure, histogram(distancesBetweenNodesInEdges,201), axis tight, grid on
% title(strcat("edges distances in iteration number ",string(iterationNr)))

%% Part 3: Create a directed graph for the FULL DATASET

G_neighFullDataset = digraph(X.VarInfo.MZRT_str(all_neigh_edges(:,1)), X.VarInfo.MZRT_str(all_neigh_edges(:,2)));
G_neighFullDataset.Nodes.CC = (conncomp(G_neighFullDataset,'Type','weak'))';
G_neighFullDataset.Nodes.idx_inVarInfo = M2S_find_idxInReference(G_neighFullDataset.Nodes.Name, X.VarInfo.MZRT_str);
%figure, plot(G_neighFullDataset,'Layout','force')

%% Plot a UMAP with edges, coloured by connected component using the FULL DATASET

if plotType > 1
    % Plot UMAP with color by CC
    % annotationsForPlot_FullDataset = X.VarInfo.MZRT_str;
    annotationsForPlot_FullDataset = X.VarInfo.CompoundName;
    groupColorForPlot_FullDataset_CC = repmat("",size(X.VarInfo,1),1);
    groupColorForPlot_FullDataset_CC(G_neighFullDataset.Nodes.idx_inVarInfo) = (string(G_neighFullDataset.Nodes.CC))';
    
    hhh0 = METU_scatterPlot_AnnotColor_grpLabel(UMAPfig_toEditFullDataset,UMAPres.reduction,annotationsForPlot_FullDataset,groupColorForPlot_FullDataset_CC,'jet',1,45,8);
    title('UMAP on full dataset with selected edges coloured by connected component');
end

%% Part 4: FIND Louvain COMMUNITIES in the full dataset (all edges)

[M_FullDataset,Q_FullDataset]=community_louvain(adjacency(G_neighFullDataset),louvain_gamma);
G_neighFullDataset.Nodes.ModuleNrs_Louvain = M_FullDataset;
fprintf('\nLouvain communities in the full dataset\n');
disp(HELP_tabulateFormatted(G_neighFullDataset.Nodes.ModuleNrs_Louvain));


% Colour UMAP according to Louvain communities
if plotType > 1
    % annotationsForPlot_FullDataset = X.VarInfo.MZRT_str;
    annotationsForPlot_FullDataset = X.VarInfo.CompoundName;
    groupColorForPlot_FullDataset_Louvain = repmat("",size(X.VarInfo,1),1);
    % Choose the type of cluster/community
    groupColorForPlot_FullDataset_Louvain(G_neighFullDataset.Nodes.idx_inVarInfo) = string(G_neighFullDataset.Nodes.ModuleNrs_Louvain);
    hhh = METU_scatterPlot_AnnotColor_grpLabel(UMAPfig_toEditFullDataset,UMAPres.reduction,annotationsForPlot_FullDataset,groupColorForPlot_FullDataset_Louvain,'jet',1,45,8);
    title('All edges coloured by Louvain community in full dataset');
end


%% Part 5: Keep only edges in which the nodes are in the same Louvain community FULL DATASET

all_neigh_edges_sel = all_neigh_edges;
node1Idx_inGraphNodes = M2S_find_idxInReference(all_neigh_edges(:,1), G_neighFullDataset.Nodes.idx_inVarInfo) ;
node2Idx_inGraphNodes = M2S_find_idxInReference(all_neigh_edges(:,2), G_neighFullDataset.Nodes.idx_inVarInfo) ;
node1LouvainCommunity = G_neighFullDataset.Nodes.ModuleNrs_Louvain(node1Idx_inGraphNodes);
node2LouvainCommunity = G_neighFullDataset.Nodes.ModuleNrs_Louvain(node2Idx_inGraphNodes);

% delete the edges which nodes are not in the same Louvain community
all_neigh_edges_sel(node1LouvainCommunity ~= node2LouvainCommunity,:) = [];




%% Create a plot without the bad edges, and with the Louvain communities

if plotType > 1
% UMAPfig_toEditFullDataset2 = copyfig(UMAPres.h_figUMAP);
UMAPfig_toEditFullDataset2 = figure;
plot(UMAPres.reduction(:,1), UMAPres.reduction(:,2),'.k')
METU_plotIsodensityLines(gcf, UMAPres.reduction); axis tight, grid on
axis tight; hold on

    METU_putEdgesOnPlot(UMAPfig_toEditFullDataset2,UMAPres.reduction,all_neigh_edges_sel(:,1),all_neigh_edges_sel(:,2),'r-');

% annotationsForPlot_FullDataset = X.VarInfo.MZRT_str;
annotationsForPlot_FullDataset = X.VarInfo.CompoundName;
groupColorForPlot_FullDataset_Louvain = repmat("",size(X.VarInfo,1),1);
% Choose the type of cluster/community
groupColorForPlot_FullDataset_Louvain(G_neighFullDataset.Nodes.idx_inVarInfo) = string(G_neighFullDataset.Nodes.ModuleNrs_Louvain);
METU_scatterPlot_AnnotColor_grpLabel(UMAPfig_toEditFullDataset2,UMAPres.reduction,annotationsForPlot_FullDataset,groupColorForPlot_FullDataset_Louvain,'jet',1,45,8);
title('Only edges with nodes in same Louvain commmunity in full dataset');
end
% close(UMAPfig_toEditFullDataset2)%% THIS IS CLOSED ONLY TO DELETE MEMORY SPACE
% 
% close(UMAPfig_toEditFullDataset)%% THIS IS CLOSED ONLY TO DELETE MEMORY SPACE


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Part 6: SELECT EDGES ACCORDING TO THE MIN NUMBER OF REPEATED EDGES IN ITERATIONS
TF_edges_sel = TF_edges(TF_edges.Count >= nNeighThreshold,:);
TF_edges = TF_edges_sel;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Intersect the selected edges in the entire dataset calculated now
% % with the table TF_edges calculated previously with the first function

all_neigh_edges_sel_str_fullDataset = strcat(string(all_neigh_edges_sel(:,1)),"_",string(all_neigh_edges_sel(:,2)));

[all_neigh_edges_sel_str_fullDataset_commonFINAL ,all_neigh_edges_sel_str_fullDataset_common_localIdx,TF_edges_common_localIdx] = intersect(all_neigh_edges_sel_str_fullDataset, TF_edges.Value,'stable');
%figure, plot(TF_edges.Count(TF_edges_common_localIdx),'.k') %% THIS IS CLOSED ONLY TO DELETE MEMORY SPACE

[~,notCommonEdges_idx] = setdiff(TF_edges.Value,TF_edges.Value(TF_edges_common_localIdx),'stable');

T_commonEdges = TF_edges(TF_edges_common_localIdx,:);
T_notcommonEdges = TF_edges(notCommonEdges_idx,:);


% HELP_tabulateFormatted(T_commonEdges.Count);
% HELP_tabulateFormatted(T_notcommonEdges.Count);

%% Part 7: Create a graph with edges common with table FINAL 

FINAL_edges_idx = str2double(split(T_commonEdges.Value, "_"));
G_FINAL = digraph(X.VarInfo.MZRT_str(FINAL_edges_idx(:,1)), X.VarInfo.MZRT_str(FINAL_edges_idx(:,2)));
G_FINAL.Nodes.CC = (conncomp(G_FINAL,'Type','weak'))';
G_FINAL.Nodes.idx_inVarInfo = M2S_find_idxInReference(G_FINAL.Nodes.Name, X.VarInfo.MZRT_str);

G_FINAL.Edges.idx_inVarInfo1 = M2S_find_idxInReference(G_FINAL.Edges.EndNodes(:,1), X.VarInfo.MZRT_str);
G_FINAL.Edges.idx_inVarInfo2 = M2S_find_idxInReference(G_FINAL.Edges.EndNodes(:,2), X.VarInfo.MZRT_str);

[M_FINAL_Dataset,Q_FINAL_Dataset]=community_louvain(adjacency(G_FINAL),louvain_gamma);
G_FINAL.Nodes.ModuleNrs_Louvain = M_FINAL_Dataset;
fprintf('\nLouvain communities in the FINAL dataset\n');
disp(HELP_tabulateFormatted(G_FINAL.Nodes.ModuleNrs_Louvain));
length(unique(G_FINAL.Nodes.ModuleNrs_Louvain))

% Colour the graph by connected components
if plotType >= 0
    UMAPfig_toEditFullDataset4 = figure; 
    plot(UMAPres.reduction(:,1), UMAPres.reduction(:,2),'.k')
    METU_plotIsodensityLines(gcf, UMAPres.reduction); axis tight, grid on
    %UMAPfig_toEditFullDataset4 = copyfig(UMAPres.h_figUMAP);
    %UMAPfig_toEditFullDataset4 = copyfig(gcf);
    axis tight; hold on;
    if plotType > 0
        METU_putEdgesOnPlot(UMAPfig_toEditFullDataset4,UMAPres.reduction,FINAL_edges_idx(:,1),FINAL_edges_idx(:,2),'k-');
    end
    annotationsForPlot_FullDataset = X.VarInfo.MZRT_str;
    %annotationsForPlot_FullDataset = X.VarInfo.CompoundName;
    groupColorForPlot_FullDataset_CCs_withEdgesFINAL = repmat("",size(X.VarInfo,1),1);
    groupColorForPlot_FullDataset_CCs_withEdgesFINAL(G_FINAL.Nodes.idx_inVarInfo) = string(G_FINAL.Nodes.CC);
    METU_scatterPlot_AnnotColor_grpLabel(UMAPfig_toEditFullDataset4,UMAPres.reduction,annotationsForPlot_FullDataset,groupColorForPlot_FullDataset_CCs_withEdgesFINAL,'jet',1,45,8);
    title('Only edges with nodes in same Connected Component  common with iteration table edges');
    close(UMAPfig_toEditFullDataset4);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Part 8: Create UMAP with edges common with table FINAL coloured by Louvain communities
% Colour by Louvain communities, but only nodes with edges are coloured
if plotType > 1
    % UMAPfig_toEditFullDataset3 = copyfig(UMAPres.h_figUMAP);
    UMAPfig_toEditFullDataset3 = figure; 
    plot(UMAPres.reduction(:,1), UMAPres.reduction(:,2),'.k')
    METU_plotIsodensityLines(gcf, UMAPres.reduction); axis tight, grid on
    axis tight; hold on
    % METU_putEdgesOnPlot(UMAPfig_toEditFullDataset3,UMAPres.reduction,all_neigh_edges_sel(all_neigh_edges_sel_str_fullDataset_common_localIdx,1),all_neigh_edges_sel(all_neigh_edges_sel_str_fullDataset_common_localIdx,2),'k-');
    % if plotType > 0
    METU_putEdgesOnPlot(UMAPfig_toEditFullDataset3,UMAPres.reduction,G_FINAL.Edges.idx_inVarInfo1,G_FINAL.Edges.idx_inVarInfo2,'k-');
    % end
    % annotationsForPlot_FullDataset = X.VarInfo.MZRT_str;
    annotationsForPlot_FullDataset = X.VarInfo.CompoundName;
    groupColorForPlot_FullDataset_Louvain_withEdges = repmat("",size(X.VarInfo,1),1);
    groupColorForPlot_FullDataset_Louvain_withEdges(G_FINAL.Nodes.idx_inVarInfo) = string(G_FINAL.Nodes.ModuleNrs_Louvain);
    %[~,final_selectedNodes_idx_inG_neighFullDatasetNodes,~] = intersect(G_neighFullDataset.Nodes.idx_inVarInfo, unique([all_neigh_edges_sel(all_neigh_edges_sel_str_fullDataset_common_localIdx,1);all_neigh_edges_sel(all_neigh_edges_sel_str_fullDataset_common_localIdx,2)]));
    % groupColorForPlot_FullDataset_Louvain_withEdges(R.FINAL_selectedNodes_idx_inVarInfo) = R.FINAL_selectedNodes_LouvainCommunity;
    METU_scatterPlot_AnnotColor_grpLabel(UMAPfig_toEditFullDataset3,UMAPres.reduction,annotationsForPlot_FullDataset,groupColorForPlot_FullDataset_Louvain_withEdges,'jet',1,45,8);
    title('Only edges with nodes in same Louvain commmunity in full dataset common with iteration table edges');
    close(UMAPfig_toEditFullDataset3);
end

%% Part 9: Collect the results
% R.FINAL_selectedNodes_idx_inVarInfo = G_neighFullDataset.Nodes.idx_inVarInfo(final_selectedNodes_idx_inG_neighFullDatasetNodes);% RESULT
% R.FINAL_selectedNodes_LouvainCommunity = string(G_neighFullDataset.Nodes.ModuleNrs_Louvain(final_selectedNodes_idx_inG_neighFullDatasetNodes));% RESULT
% [ap1,ap2,ap3] = intersect(X.VarInfo.MZRT_str(R.FINAL_selectedNodes_idx_inVarInfo),G_FINAL.Nodes.Name, 'stable');
% R.FINAL_selectedNodes_CC = G_FINAL.Nodes.CC(ap2);
% R.FINAL_edges_idx = all_neigh_edges_sel(all_neigh_edges_sel_str_fullDataset_common_localIdx,:);
% R.UMAPres = UMAPres;
R.FINAL_selectedNodes_idx_inVarInfo = G_FINAL.Nodes.idx_inVarInfo;% RESULT
R.FINAL_selectedNodes_LouvainCommunity = string(G_FINAL.Nodes.ModuleNrs_Louvain);% RESULT
%[ap1,ap2,ap3] = intersect(X.VarInfo.MZRT_str(R.FINAL_selectedNodes_idx_inVarInfo),G_FINAL.Nodes.Name, 'stable');
R.FINAL_selectedNodes_CC = G_FINAL.Nodes.CC;
R.FINAL_edges_idx = [G_FINAL.Edges.idx_inVarInfo1, G_FINAL.Edges.idx_inVarInfo2];
R.UMAPres = UMAPres;
%{

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% WORK ON THE SELECTED EDGES FROM THE PREVIOUS FUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Create a directed graph for the selected nodes
G_neigh = digraph(X.VarInfo.MZRT_str(edgeNode1_sel), X.VarInfo.MZRT_str(edgeNode2_sel))
G_neigh.Nodes.CC = (conncomp(G_neigh,'Type','weak'))';
G_neigh.Nodes.idx_inVarInfo = M2S_find_idxInReference(G_neigh.Nodes.Name, X.VarInfo.MZRT_str)
%figure, plot(G_neigh,'Layout','force')

%% Create a UMAP plot using the entire dataset
% metricType = 'cosine'
% UMAPres = METU_UMAPregular(zscore(X.Data_imputed)', nNeighUMAP,metricType)
% UMAPfig_toEditFullDataset = copyfig(gcf);

METU_putEdgesOnPlot(UMAPfig_toEditFullDataset,UMAPres.reduction,edgeNode1_sel,edgeNode2_sel,'r-')
% tempUMAPdist = UAU_dist_inPercent(UMAPres.reduction,(1:size(UMAPres.reduction,1))',(1:size(UMAPres.reduction,1))');

% Plot UMAP
annotationsForPlot = X.VarInfo.MZRT_str;
groupColorForPlot = repmat("",size(X.VarInfo,1),1)
groupColorForPlot(G_neigh.Nodes.idx_inVarInfo) = (string(G_neigh.Nodes.CC))';

METU_scatterPlot_AnnotColor_grpLabel(UMAPfig_toEditFullDataset,UMAPres.reduction,annotationsForPlot,groupColorForPlot,'jet',1,45,8)
title('UMAP on full dataset with selected edges and Louvain communities')
close(UMAPfig_toEditFullDataset)

%% FIND Louvain COMMUNITIES 

[M,Q]=community_louvain(adjacency(G_neigh),louvain_gamma)
G_neigh.Nodes.ModuleNrs_Louvain = M;

figure, histogram(G_neigh.Nodes.ModuleNrs_Louvain,'BinMethod','integers')
axis tight, grid on, title('Number of members of each Louvain community in full dataset')
HELP_tabulateFormatted(G_neigh.Nodes.ModuleNrs_Louvain)


figure, histogram(G_neigh.Nodes.CC,'BinMethod','integers')
axis tight, grid on, title('Number of members of each conncomp in full dataset')

%% Decompose all edges in indices

parts = split(sortedEdges_all_str, "_");             % split into 2 columns
edgeNode1 = str2double(parts(:,1));       % first column -> numbers
edgeNode2 = str2double(parts(:,2)); 


%{
sortedEdges_sel = NaN(size(TF_edges_sel,1),2);
for e=1:size(TF_edges_sel,1)
    temp = double(strsplit(TF_edges.Value(e), '_'));
    sortedEdges_sel(e,:) = [temp(1), temp(2)];
end
final_neigh_idx = unique([sortedEdges_sel(:,1); sortedEdges_sel(:,2)]);
%}

%% Previous version

idx_edgesGood = M2S_find_idxInReference(sortedEdges_all_str,TF_edges_sel.Value);
%{
sortedEdges_sel = sortedEdges_all(~isnan(idx_edgesGood),:);
sortedEdges_sel = unique(sortedEdges_sel,'rows'); % RESULT: These are the indexes of neighbours!!!
% Get the idx of features that are neighbours
final_neigh_idx = unique([sortedEdges_sel(:,1);sortedEdges_sel(:,2)]);
%}
sortedEdges_str_sel = unique(sortedEdges_all_str(~isnan(idx_edgesGood)));
% sortedEdges_sel = arrayfun(@strsplit, sortedEdges_str_sel, repmat('_',length(sortedEdges_str_sel),1),'UniformOutput', false)

sortedEdges_sel = NaN(length(sortedEdges_str_sel),2);
for e=1:length(sortedEdges_str_sel)
    temp = double(strsplit(sortedEdges_str_sel(e), '_'));
    sortedEdges_sel(e,:) = [temp(1), temp(2)];
end

final_neigh_idx = unique([sortedEdges_sel(:,1); sortedEdges_sel(:,2)]);
%}