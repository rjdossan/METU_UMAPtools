% addClassMedianLabel(gcf ,plottingData,classLabels,'r')

function [classLabelsUnique, classMedianVals, tH] = METU_addClassMedianLabel(h_f,plottingData,classLabels,textColor,textSize,medianOrMean,visibility_onOff)
if nargin<4
    textColor = 'k'
    textSize = 12
elseif nargin<5
    textSize = 12
elseif nargin<6
    medianOrMean = 'median';
elseif nargin<7
    visibility_onOff = 'on';
end

if iscell(classLabels)
    classLabels=string(classLabels);
end
%{
classLabelsUnique = unique(classLabels);
classLabelsUnique(classLabelsUnique=="") = [];
classMedianVals = NaN(length(classLabelsUnique),2);
plottingData_zscored = zscore(plottingData);
for u=1:length(classLabelsUnique)
    idx=find(classLabels == classLabelsUnique(u));
    % classMedianVals(u,:) = median(plottingData(idx,:));
    currentGroupData = plottingData_zscored(idx,:);
    currentGroupData_distToZero = vecnorm(currentGroupData, 2, 2);
    if strcmp(medianOrMean,'median')
        medianDistValue = median(currentGroupData_distToZero);
    else
        medianDistValue = mean(currentGroupData_distToZero);
    end
    % Find the index of the median value
    [~, medianIndex] = min(abs(currentGroupData_distToZero - medianDistValue));
    classMedianVals(u,:) = plottingData(idx(medianIndex),:);
end
hold on
%text(classMedianVals(:,1),classMedianVals(:,2),classLabelsUnique,'fontsize',textSize,'color',textColor)
text(classMedianVals(:,1),classMedianVals(:,2),classLabelsUnique,'fontsize',textSize,'color',textColor,'verticalalignment','top')
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classLabelsUnique = unique(classLabels);
classLabelsUnique(classLabelsUnique=="") = [];
classMedianVals = NaN(length(classLabelsUnique),2);
plottingData_zscored = zscore(plottingData);
for u=1:length(classLabelsUnique)
    idx=find(classLabels == classLabelsUnique(u));
    % classMedianVals(u,:) = median(plottingData(idx,:));
    currentGroupData = plottingData_zscored(idx,:);  
    if strcmp(medianOrMean,'median')
        currentGroupData_distToZero = vecnorm(currentGroupData, 2, 2);
        medianDistValue = median(currentGroupData_distToZero);
        [~, medianIndex] = min(abs(currentGroupData_distToZero - medianDistValue));
    else
        medianDistValue = mean(currentGroupData);
        %centroid = mean(currentGroupData, 1);   % 1×2 vector: [meanX, meanY]
        diffs = currentGroupData - medianDistValue;    % N×2 matrix of [x−meanX, y−meanY]
        distances = sqrt(sum(diffs.^2, 2));   
        [minDist, medianIndex] = min(distances);
    end
    % Find the index of the median value
    classMedianVals(u,:) = plottingData(idx(medianIndex),:);
end
hold on
%text(classMedianVals(:,1),classMedianVals(:,2),classLabelsUnique,'fontsize',textSize,'color',textColor)
if strcmp(medianOrMean,'median')
    tH = text(classMedianVals(:,1),classMedianVals(:,2),classLabelsUnique,'fontsize',textSize,'color',textColor,'verticalalignment','top','Visible',visibility_onOff);
else
    tH = text(classMedianVals(:,1),classMedianVals(:,2),classLabelsUnique,'fontsize',textSize,'color',textColor,'horizontalalignment','center','verticalalignment','middle','Visible',visibility_onOff);
end

