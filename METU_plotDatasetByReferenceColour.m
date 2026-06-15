%% Plot a scatterplot with colors according to a labels vector common to all desired plots




function [h_fig1, h_FigLegend] = METU_plotDatasetByReferenceColour(data1, groupColorForPlot, h_figure1, uniqueColors)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
% EXAMPLE WITH N DATASETS (run N times)
data1 = plottingData_FINGER;
groupColorForPlot = FINGER_Lipidomics_VarInfo_updated.predicted_RefMet_Main_class;
h_figure1 = h_figUMAP_FINGER;
uniqueColors = ALLdatasetsLabels; % unique labels across all datasets, including
the labels in groupColorForPlot
%}
%%%%%%%%%%%%%%%%%%%%%%%%%

%% DEFINITIONS
randomizeColours = 0;
transparencyValue = 1;
LegendMarkerSize = 20;
%LegendFontSize = 10;
LegendFontSize = 8;
markerSize = 50;
markerSizeBlackDots = 5;

%%%%%%%%%%%%%%%%%%%%%%%%%

% uniqueColors = unique([groupColorForPlot;uniqueColors]);
uniqueColors = unique(uniqueColors);
uniqueColors(uniqueColors=="")= [];
uniqueColors_noEmpty = uniqueColors;



isIntegerString = false(size(uniqueColors_noEmpty));

% Loop through each string in the vector
for i = 1:length(uniqueColors_noEmpty)
    % Check if the string can be converted to a number and if it's an integer
    num = str2double(uniqueColors_noEmpty{i});
    if ~isnan(num) && mod(num, 1) == 0
        isIntegerString(i) = true;
    end
end
if sum(isIntegerString) == length(isIntegerString)
    sortedUniqueColors = sort(str2double(uniqueColors_noEmpty),'descend');
    sortedUniqueColors = string(sortedUniqueColors);
else
    sortedUniqueColors = sort(uniqueColors_noEmpty,'descend');
end





rng(10);
colors = jet(length(sortedUniqueColors));
if randomizeColours == 1
    colors = colors(randperm(size(colors,1)),:);
end

possibleShapes = ('osdvp^h')';
nRepmatShapes = length(sortedUniqueColors)/length(possibleShapes);
uniqueShapes = repmat(possibleShapes,ceil(nRepmatShapes),1);
uniqueShapes = uniqueShapes(1:length(sortedUniqueColors));
uniqueShapes = repmat(possibleShapes,1000,1);

% Plot 

h_fig1 = HELP_duplicateFigure(h_figure1); hold on;
tempidx1 = find(groupColorForPlot == "");
hold on, plot(data1(tempidx1, 1), data1(tempidx1, 2), '.k','markersize',markerSizeBlackDots);

% h_fig2 = duplicateFigure(h_figure2); hold on;
% tempidx2 = find(uniqueColors == "");
% hold on, plot(data2(tempidx2, 1), data2(tempidx2, 2), '.k');

for i = 1:length(sortedUniqueColors)
    
    figure(h_fig1);
    tempidx1 = find(groupColorForPlot == sortedUniqueColors(i));
    scatter(data1(tempidx1, 1), data1(tempidx1, 2), markerSize*ones(size(tempidx1)), colors(i, :), uniqueShapes(i),'filled','MarkerEdgeColor','k','MarkerEdgeAlpha',0.5,'MarkerFaceAlpha',transparencyValue);

    % figure(h_fig2);
    % tempidx2 = find(uniqueColors == uniqueColors(i));
    % scatter(data2(tempidx2, 1), data2(tempidx2, 2), [], colors(i, :), uniqueShapes(i),'filled','MarkerEdgeColor','k','MarkerEdgeAlpha',0.5,'MarkerFaceAlpha',transparencyValue);
end


%% get the legend 


legendTable = table(sortedUniqueColors,colors,'VariableNames',{'label','color'});
%axes('Position',[.01 .7 .1 .23]), box on
h_FigLegend = figure;
if ~isempty(legendTable)
    %axes('Position',[0.01 0.1 0.16 0.83]); box on;
    plot(0.25,1,'.k'), hold on
    text(0.25,1,"   unknown",'fontsize',LegendFontSize,'verticalalignment','middle','interpreter','none');
    % scatter(1,1,50,colors(1,:),uniqueShapes(1),'filled'); hold on;
    if length(sortedUniqueColors)>1
        for c=1:length(sortedUniqueColors)
            scatter(0.25,c+1,LegendMarkerSize,colors(c,:),uniqueShapes(c),'filled','MarkerEdgeColor','k','MarkerEdgeAlpha',0.5,'MarkerFaceAlpha',transparencyValue);
        end
    end
    text(0.25*ones(length(sortedUniqueColors),1),(2:length(sortedUniqueColors)+1)',...
        strcat("   ",sortedUniqueColors),'fontsize',LegendFontSize,'verticalalignment','middle','interpreter','none');
    set(gca,'xticklabels',[]);set(gca,'yticklabels',[]);
    set(gca,'xtick',[]);set(gca,'ytick',[]);
    xlim([0,4]);ylim([0,length(sortedUniqueColors)+2]);
else
    h_FigLegend = 0;
end

% give the handle to the actual figure, not the legend
figure(h_fig1);



