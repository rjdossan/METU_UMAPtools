function seq = greedyClusterSequence(x, y, clusterIDs)


% GREEDYCLUSTERSEQUENCE  Return an ordering of clusters by starting at
% the cluster with the highest mean-y, then repeatedly jumping to the
% nearest unvisited cluster (Euclidean distance of centroids).
%
%  seq = greedyClusterSequence(x, y, clusterIDs)
%
%  Inputs:
%    x, y         – vectors of point coordinates (N×1)
%    clusterIDs   – vector of cluster labels (N×1), numeric or categorical
%
%  Output:
%    seq          – 1×K vector of cluster labels in the visitation order

    % 1) compute each cluster’s centroid
    [G, uniqueClusters] = findgroups(clusterIDs);
    centroids = [ splitapply(@mean, x, G), ...
                  splitapply(@mean, y, G) ];  % K×2
    % for g=1:max(uniqueClusters)
    %     centroids = mean()

    K = numel(uniqueClusters);

    % 2) find the “start” cluster = one with highest y‐centroid
    [~, startIdx] = max(centroids(:,2));
    seq = uniqueClusters(startIdx);      % initialize output
    visited = false(K,1);
    visited(startIdx) = true;
    currentPoint = centroids(startIdx,:);

    % 3) greedy loop: pick nearest unvisited centroid each time
    for i = 2:K
        unv = find(~visited);           % indices of unvisited clusters
        diffs = centroids(unv,:) - currentPoint;      % M×2
        dists = sqrt(sum(diffs.^2,2));                % M×1
        [~, loc] = min(dists);         % index in “unv”
        nextIdx = unv(loc);            % global index of that cluster

        % update
        seq(end+1) = uniqueClusters(nextIdx);
        visited(nextIdx) = true;
        currentPoint = centroids(nextIdx,:);
    end
end
