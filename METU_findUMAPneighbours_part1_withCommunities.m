function [sortedEdges_all_str,TF_edges,all_neigh_edges,G_neigh,stats] = METU_findUMAPneighbours_part1_withCommunities...
    (X,onlyAdductsIsotopes,UMAPthreshDist_betweenNeighbours,nNeighUMAP, n_iterations, nSamples, louvain_gamma,plotType)

stats=[];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% THIS IS DEFINED BECAUSE THE FUNCTION REPEATS "FOR" ITERATIONS IN CASE THERE ARE 
%% OUTLIERS, AND COULD CONTINUE DOING SO INDEFINITELY
maxNrForIterations = 20;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Define if it is to use whole dataset or only adducts and isotopes
%% FOR ONLY ADDUCTS AND ISOTOPES
if onlyAdductsIsotopes == 1
    RTthresh_neighbours_iterations =  0.25/60;
    corrThresh = 0.75;
else
%% FOR EVERYTHING there is no restriction of RT or correlation
    RTthresh_neighbours_iterations =  1000;
    corrThresh = -1;
end

% Only calculate correlations if needed
if corrThresh>-1
    calculate_corr = true;
else
    calculate_corr = false;
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% Define the samples to use in each iteration
% withResampling = 1;
% % create an index of random samples
% randSamples_idx = randperm(size(X.Data,1))';
% if withResampling == 0
%     %% Check if the number of samples and iterations will allow the calculations
%     % If not, give ERROR
% 
%     if n_iterations * nSamples> size(X.Data,1)
%         fprintf('\n**************************************************************\n')
%         disp('***ERROR*** Nothing will be calculated.')
%         disp('This number of samples and iterations cannot be used.')
%         fprintf('\nWith this number of samples per iteration, you can only do %d iterations.\nOR\n',floor(size(X.Data,2)/nSamples))
%         fprintf('With this number of iterations, you can only use %d samples in each iteration.\n',floor(size(X.Data,2)/n_iterations))
%         fprintf('**************************************************************\n')
%     error('CustomError:TooManySamplesIterations', ...
%               'There are too many samples or iterations, choose other values');
%     end
% 
%     idxLeft = 1;
%     idxRight = nSamples;
%     % Get the idx of each sample in each of the iterations
%     for iterationNr=1:n_iterations
%         samplesIdx_iteration{iterationNr,1} = randSamples_idx(idxLeft:idxRight);
%         idxLeft = idxRight+1;
%         idxRight = idxRight+nSamples;
%     end
% else
%     for iterationNr=1:maxNrForIterations %%*************************************************************************************
%         % create a random index vector
%         randSamples_idx2 = randperm(size(X.Data,1))';
%         samplesIdx_iteration{iterationNr,1} = randSamples_idx(randSamples_idx2(1:nSamples));
%     end
% end
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialisations
neighbours_idx_all = cell(n_iterations,1);
all_neigh_edges_twoColumns_all = cell(n_iterations,1);
all_neigh_edges = cell(n_iterations,1);

sortedEdges_all = [];
sortedEdges_all_str = repmat("",size(sortedEdges_all,1),1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part 1: For each iteration, create a UMAP
disp('PART 1: UMAP plots and get coordinates')
% rng(3)
iterationNr = 0
real_iterationNr=0;
numberOfForIterations = 0;
jumpIteration = 0;
for iterationNr_i = 1:maxNrForIterations
    iterationNr = iterationNr+1;
    real_iterationNr = real_iterationNr+1;
    numberOfForIterations = numberOfForIterations+1;
    
    % Check if the iteration needed to be repeated, use same iterationNr
    if jumpIteration == 1 & numberOfForIterations <= maxNrForIterations;
        iterationNr = iterationNr-1;
        jumpIteration = 0;
        disp('Repeated the iteration because UMAP had big outliers')
        real_iterationNr = real_iterationNr-1;
    end
    % Check if the desired number of iterations has been reached. In case
    % yes, jump over all the next iterations
    if iterationNr > n_iterations
        real_iterationNr = real_iterationNr-1;
        continue
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Define the samples to use in each iteration
    withResampling = 1; % it is always with resampling, otherwise there may not be samples enough

    % create an index of random samples
    randSamples_idx = randperm(size(X.Data,1))';
    if withResampling == 0
        %% Check if the number of samples and iterations will allow the calculations
        % If not, give ERROR
        
        if n_iterations * nSamples> size(X.Data,1)
            fprintf('\n**************************************************************\n')
            disp('***ERROR*** Nothing will be calculated.')
            disp('This number of samples and iterations cannot be used.')
            fprintf('\nWith this number of samples per iteration, you can only do %d iterations.\nOR\n',floor(size(X.Data,2)/nSamples))
            fprintf('With this number of iterations, you can only use %d samples in each iteration.\n',floor(size(X.Data,2)/n_iterations))
            fprintf('**************************************************************\n')
        error('CustomError:TooManySamplesIterations', ...
                  'There are too many samples or iterations, choose other values');
        end
    
        idxLeft = 1;
        idxRight = nSamples;
        % Get the idx of each sample in each of the iterations
        %for iterationNr=1:n_iterations
            samplesIdx_iteration{iterationNr,1} = randSamples_idx(idxLeft:idxRight);
            idxLeft = idxRight+1;
            idxRight = idxRight+nSamples;
        %end
    else % without Resampling
        %for iterationNr=1:maxNrForIterations %%*************************************************************************************
            % create a random index vector
            randSamples_idx2 = randperm(size(X.Data,1))';
            samplesIdx_iteration{iterationNr,1} = randSamples_idx(randSamples_idx2(1:nSamples));
        %end
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    
    randSamples_idx = samplesIdx_iteration{iterationNr,1};
    selData_Zscored = X.Data(randSamples_idx,:);
    selData_Zscored = zscore(selData_Zscored);
    
    % Variables that could not be zscored receive random values
    idxColNaN = find(sum(isnan(selData_Zscored))>0);
    selData_Zscored(:,idxColNaN) = randn(size(selData_Zscored,1),length(idxColNaN));
    if ~isempty(idxColNaN);
        disp('The following variables indices have NaN:')
        disp(idxColNaN)
    end
    dataForUMAP_iterations = selData_Zscored';   
    

    % Perform UMAP
    UMAPres.opt.U_neighbors=nNeighUMAP;
    UMAPres.opt.U_components = 2;
    UMAPres.opt.U_metric = 'cosine';
    UMAPres.opt.U_epochs = 1000;
    UMAPres.opt.min_dist = 0.3;% default is 0.3
    [UMAPres.reduction,UMAPres.Xumap] = run_umap(dataForUMAP_iterations,...
        'n_neighbors',UMAPres.opt.U_neighbors,...
        'n_components',UMAPres.opt.U_components,...
        'metric',UMAPres.opt.U_metric,...
        'min_dist',UMAPres.opt.min_dist,...
        'n_epochs',UMAPres.opt.U_epochs);%,...
        %'verbose','none');
    % close
    

% Check if there are outliers, and in case yes, skip this iteration but run
% again the same iteration index.
    nMAD=6;
    [outliersIdx] = METU_findUMAPoutliers(UMAPres.reduction,0,nMAD);
    if ~isempty(outliersIdx) 
        jumpIteration = 1;
        close
        continue % skips this iteration of FOR cycle
    end

    UMAPfig(iterationNr).h_UMAPfig = gcf; % *** FIGURE ***
    plottingData_iterations_all{iterationNr,1} = UMAPres.reduction;

    %% Duplicate the UMAP figure to write over it
    
    if plotType > 1
        
        UMAPfig_toEdit(iterationNr).figure = figure, plot(UMAPres.reduction(:,1), UMAPres.reduction(:,2),'.k');
        METU_plotIsodensityLines(gcf, UMAPres.reduction); axis tight, grid on

        
        %UMAPfig_toEdit(iterationNr).figure = copyfig(UMAPfig(iterationNr).h_UMAPfig); % *** FIGURE ***
        axis tight; 
        title(strcat("iteration number ", string(iterationNr))); 
        
        % Create the diagonal line in the UMAP figure
        % h_ap = figure, plot(plottingData_iterations(:,1), plottingData_iterations(:,2),'.k'), axis tight, grid on;
        [diagLength, h] = METU_plotDiagonalORcircle(UMAPfig_toEdit(iterationNr).figure); % close(h_ap)
    end
    
    % Normalise and collect UMAP distances in percentage of diagonal
    plottingData_temp = plottingData_iterations_all{iterationNr,1} - repmat(min(plottingData_iterations_all{iterationNr,1}),size(plottingData_iterations_all{iterationNr,1},1),1);
    plottingData_between01 = plottingData_temp ./ repmat(max(plottingData_temp),size(plottingData_temp,1),1);
    UMAPcoordPerc_iterations{iterationNr,1} = plottingData_between01;
% end
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % In case correlations are to be calculated:
    if calculate_corr
        %for iterationNr = 1:n_iterations
            allCorr{iterationNr,1} = corr(UMAPcoordPerc_iterations{iterationNr,1}','type','Pearson'); 
        %end
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Part 2: Find neighbours at a certain max distance from each other


%for iterationNr = 1:n_iterations
    temp_UMAPcoordPerc_iterations = UMAPcoordPerc_iterations{iterationNr,1};
    neighbours_idx = cell(size(temp_UMAPcoordPerc_iterations,1),1);
    %all_neigh_edges_twoColumns = [];
    temp_all_neigh_edges = [];
    for varNr = 1:size(temp_UMAPcoordPerc_iterations,1)
        %[iterationNr,n_iterations, varNr, size(X.Data,2)]
        tempDiff1 = temp_UMAPcoordPerc_iterations(:,1) - temp_UMAPcoordPerc_iterations(varNr,1);
        tempDiff2 = temp_UMAPcoordPerc_iterations(:,2) - temp_UMAPcoordPerc_iterations(varNr,2);
        tempDist = sqrt(tempDiff1.^2 + tempDiff2.^2);
        tempDist(varNr) = 1000; % to not match itself
        if calculate_corr
            current_neighbours_idx = find(tempDist < UMAPthreshDist_betweenNeighbours & allCorr{iterationNr,1}(:,varNr)>= corrThresh); % This could be calculated using correlations
        else
            current_neighbours_idx = find(tempDist < UMAPthreshDist_betweenNeighbours);
        end
        local_acceptedIdx = find(abs(X.VarInfo.rtmed(current_neighbours_idx) - X.VarInfo.rtmed(varNr)) < RTthresh_neighbours_iterations);
        neighbours_idx{varNr} = current_neighbours_idx(local_acceptedIdx);
        %all_neigh_edges{iterationNr,1} = [all_neigh_edges{iterationNr,1};[repmat(varNr,length(current_neighbours_idx(local_acceptedIdx)),1),current_neighbours_idx(local_acceptedIdx)]];
        temp_all_neigh_edges = [temp_all_neigh_edges;[repmat(varNr,length(current_neighbours_idx(local_acceptedIdx)),1),current_neighbours_idx(local_acceptedIdx)]];
    end
    all_neigh_edges{iterationNr,1} = temp_all_neigh_edges(temp_all_neigh_edges(:,1) < temp_all_neigh_edges(:,2),:);
    %all_neigh_edges_twoColumns = [all_neigh_edges_twoColumns; all_neigh_edges{iterationNr,1}];

    if plotType > 1
        % Plot the edges between neighbours
        figure(UMAPfig_toEdit(iterationNr).figure); hold on % *** FIGURE ***
        METU_putEdgesOnPlot(UMAPfig_toEdit(iterationNr).figure,plottingData_iterations_all{iterationNr,1},all_neigh_edges{iterationNr,1}(:,1),all_neigh_edges{iterationNr,1}(:,2),'r-');
    end

    neighbours_idx_all{iterationNr,1} = neighbours_idx;  
    %all_neigh_edges_twoColumns_all{iterationNr,1} = all_neigh_edges_twoColumns;

    % Plot all the edges distances THIS WORKS!
    % distancesBetweenNodesInEdges=[];for a=1:length(all_neigh_edges_twoColumns_all{iterationNr,1} (:,1)); distancesBetweenNodesInEdges(a,1) = pdist([temp_UMAPcoordPerc_iterations(all_neigh_edges_twoColumns_all{iterationNr,1} (a,1));temp_UMAPcoordPerc_iterations(all_neigh_edges_twoColumns(a,2))],'euclidean');end
    % figure, histogram(distancesBetweenNodesInEdges,201), axis tight, grid on
    % title(strcat("edges distances in iteration number ",string(iterationNr)))
%end % XXXXXXXXXXXXXXXXXXXX

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
%% Part 3: Define graphs with neighbours and detect Louvain communities within the graphs

%for iterationNr = 1:n_iterations % XXXXXXXXXXXXXXXXXXXX
    % all_neigh_edges_unique{iterationNr} = all_neigh_edges{iterationNr};
    % all_neigh_edges_unique{iterationNr}  = (sort(all_neigh_edges_unique{iterationNr}'))';
    % all_neigh_edges_unique{iterationNr}  = unique(all_neigh_edges_unique{iterationNr} ,'rows');
    G_neigh(iterationNr).graph = digraph(X.VarInfo.MZRT_str(all_neigh_edges{iterationNr}(:,1)), X.VarInfo.MZRT_str(all_neigh_edges{iterationNr}(:,2)));
    G_neigh(iterationNr).graph.Nodes.CC = (conncomp(G_neigh(iterationNr).graph,'Type','weak'))';
    G_neigh(iterationNr).graph.Nodes.idx_inVarInfo = M2S_find_idxInReference(G_neigh(iterationNr).graph.Nodes.Name, X.VarInfo.MZRT_str);
    
    
    [M,Q]=community_louvain(adjacency(G_neigh(iterationNr).graph),louvain_gamma);
    G_neigh(iterationNr).graph.Nodes.ModuleNrs_Louvain = M;
    
    % Show the size of the communities
    %figure, histogram(G_neigh(iterationNr).graph.Nodes.ModuleNrs_Louvain,'BinMethod','integers')
    fprintf('\nLouvain communities of iteration %d\n',num2str(iterationNr));
    disp(HELP_tabulateFormatted(G_neigh(iterationNr).graph.Nodes.ModuleNrs_Louvain));

    if plotType > 1
        annotationsForPlot = X.VarInfo.MZRT_str; % *** FIGURE ***
        groupColorForPlot = repmat("",size(X.VarInfo,1),1);
        % Choose the type of cluster/community
        groupColorForPlot(G_neigh(iterationNr).graph.Nodes.idx_inVarInfo) = string(G_neigh(iterationNr).graph.Nodes.ModuleNrs_Louvain);
        F(iterationNr).hUMAPfig_withBadEdges_Louvain = METU_scatterPlot_AnnotColor_grpLabel(UMAPfig_toEdit(iterationNr).figure,plottingData_iterations_all{iterationNr,1},annotationsForPlot,groupColorForPlot,'jet',1,45,8);
        title(strcat("iteration number ", string(iterationNr)))
    end
    %close(UMAPfig_toEdit(iterationNr).figure)%% THIS IS CLOSED ONLY TO DELETE MEMORY SPACE
    %close(F(iterationNr).hUMAPfig_withBadEdges_Louvain)%% THIS IS CLOSED ONLY TO DELETE MEMORY SPACE
% end % XXXXXXXXXXXXXXXXXXXX

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Part 4: Keep only edges in which the nodes are in the same Louvain community

    % all_neigh_edges_sel{iterationNr} = all_neigh_edges{iterationNr};
%for iterationNr=1:length(all_neigh_edges) % XXXXXXXXXXXXXXXXXXXX
    node1Idx_inGraphNodes = M2S_find_idxInReference(all_neigh_edges{iterationNr}(:,1), G_neigh(iterationNr).graph.Nodes.idx_inVarInfo) ;
    node2Idx_inGraphNodes = M2S_find_idxInReference(all_neigh_edges{iterationNr}(:,2), G_neigh(iterationNr).graph.Nodes.idx_inVarInfo) ;
    node1LouvainCommunity = G_neigh(iterationNr).graph.Nodes.ModuleNrs_Louvain(node1Idx_inGraphNodes);
    node2LouvainCommunity = G_neigh(iterationNr).graph.Nodes.ModuleNrs_Louvain(node2Idx_inGraphNodes);

    % delete the edges which nodes are not in the same Louvain community
    all_neigh_edges{iterationNr}(node1LouvainCommunity ~= node2LouvainCommunity,:) = [];
%end % XXXXXXXXXXXXXXXXXXXX

%% Part 5: Create a plot without the bad edges, but with the Louvain communities

    if plotType > 0
                                                                    % *** FIGURE ***
    %for iterationNr = 1:n_iterations % XXXXXXXXXXXXXXXXXXXX
        annotationsForPlot = X.VarInfo.MZRT_str;
        groupColorForPlot = repmat("",size(X.VarInfo,1),1);
        % Choose the type of cluster/community
        groupColorForPlot(G_neigh(iterationNr).graph.Nodes.idx_inVarInfo) = string(G_neigh(iterationNr).graph.Nodes.ModuleNrs_Louvain);
        if plotType > 1
            tempfigureSelectedEdges = figure; plot(UMAPres.reduction(:,1), UMAPres.reduction(:,2),'.k')
            METU_plotIsodensityLines(gcf, UMAPres.reduction); axis tight, grid on
            %tempfigureSelectedEdges = copyfig(UMAPfig(iterationNr).h_UMAPfig);
            METU_putEdgesOnPlot(tempfigureSelectedEdges,plottingData_iterations_all{iterationNr,1},all_neigh_edges{iterationNr}(:,1),all_neigh_edges{iterationNr}(:,2),'k-');
            F(iterationNr).hUMAPfig_FINAL_Louvain = METU_scatterPlot_AnnotColor_grpLabel(tempfigureSelectedEdges,plottingData_iterations_all{iterationNr,1},annotationsForPlot,groupColorForPlot,'jet',1,45,8);
    
        else
            F(iterationNr).hUMAPfig_FINAL_Louvain = METU_scatterPlot_AnnotColor_grpLabel(UMAPfig(iterationNr).h_UMAPfig,plottingData_iterations_all{iterationNr,1},annotationsForPlot,groupColorForPlot,'jet',1,45,8);
        end
        axis tight; hold on
        title(strcat("Only edges with nodes in same Louvain commmunity - iteration number ", string(iterationNr)));
        % close(tempfigureSelectedEdges);
    end
    %end

    close(UMAPfig(iterationNr).h_UMAPfig) %% THIS IS CLOSED ONLY TO DELETE MEMORY SPACE

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Part 6: Prepare to create a string with the edges
%{
sortedEdges_all = [];
for iterationNr=1:n_iterations
    sortedEdges_all = [sortedEdges_all;all_neigh_edges_sel{iterationNr}];  
end
sortedEdges_all_str = repmat("",size(sortedEdges_all,1),1);
for idx = 1:size(sortedEdges_all,1)
    sortedEdges_all_str(idx,1) = strcat(string(sortedEdges_all(idx,1)),"_",string(sortedEdges_all(idx,2)));
end
%}
%% Prepare to create a string with the edges
    sortedEdges_all = [sortedEdges_all;all_neigh_edges{iterationNr}];  
end

sortedEdges_all_str = strcat(string(sortedEdges_all(:,1)),"_",string(sortedEdges_all(:,2)));
TF_edges = HELP_tabulateFormatted(sortedEdges_all_str);

%% Collect stats
stats.numberOfForIterations = real_iterationNr;
%stats.nCalculatedIterations = iterationNr;


