%% Find outliers in UMAP
% Outliers are the ones far from the median (a number of MADs away). 

function [outliersIdx] = METU_findUMAPoutliers(plottingData,plotOrNot,nMAD)

if nargin == 1
    plotOrNot = 0;
    nMAD = 6;
elseif nargin == 2
    nMAD = 6;
end

% normUMAPdist = UAU_dist_inPercent(plottingData,(1:size(plottingData,1))',(1:size(plottingData,1))');
% The median value for all points is calculated for each dimension
% normMedianVal = median(normUMAPdist.data_std01);

plottingData_temp = plottingData - repmat(min(plottingData),size(plottingData,1),1);
plottingData_between01 = plottingData_temp ./ repmat(max(plottingData_temp),size(plottingData_temp,1),1);
normMedianVal = median(plottingData_between01);

normDistToMedian = plottingData_between01 - repmat(normMedianVal,size(plottingData_temp,1),1);
normEuclideanDist = sqrt(normDistToMedian(:,1).^2 + normDistToMedian(:,2).^2);

% Outliers
outliersIdx = find(normEuclideanDist >  median(normEuclideanDist) + nMAD * mad(normEuclideanDist));

if plotOrNot == 1
    figure('Position',[414         378        1377         420])
    
    % subplot(1,2,1), plot(normUMAPdist.data_std01(:,1),normUMAPdist.data_std01(:,2),'.k');
    subplot(1,2,1), plot(plottingData_between01(:,1),plottingData_between01(:,2),'.k');
    hold on, axis tight, grid on
    plot(normMedianVal(1),normMedianVal(2),'+r','markersize',20,'linewidth',4)
    plot(plottingData_between01(outliersIdx,1),plottingData_between01(outliersIdx,2),'or');
    title('UMAP points and outliers')
    
    % Draw circle with threshold
    th = 0:pi/50:2*pi;
    xunit = nMAD * mad(normEuclideanDist) * cos(th) + normMedianVal(1);
    yunit = nMAD * mad(normEuclideanDist) * sin(th) + normMedianVal(2);
    plot(xunit, yunit,'-r');
    
    subplot(1,2,2), plot((normEuclideanDist),'.k'), axis tight, grid
    hold on, plot(outliersIdx,normEuclideanDist(outliersIdx),'or')
    plot(xlim',[nMAD * mad(normEuclideanDist);nMAD * mad(normEuclideanDist)],'-r')
    title('Normalized euclidean distance to median')
    
    
end
