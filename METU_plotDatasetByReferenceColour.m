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
% markerSizeBlackDots = 1;

% This defines if legend only contains strings for the current dataset
onlyExistingAnnotationsForLabel = 1

%%%%%%%%%%%%%%%%%%%%%%%%%

%% Define uniqueColors_noEmpty from uniqueColors

% uniqueColors = unique([groupColorForPlot;uniqueColors]);
uniqueColors = unique(uniqueColors);
uniqueColors(uniqueColors=="")= [];
uniqueColors_noEmpty = uniqueColors;

% % Define if the legend only contains groups existing in the current dataset
% % or all the groups in the "uniqueColors", in general common to many datasets
% if onlyExistingAnnotationsForLabel == 1
%     [uniqueColors_noEmpty,uniqueColors_noEmpty_initialIdx] = intersect(uniqueColors_noEmpty, groupColorForPlot,'stable');
% else
%     uniqueColors_noEmpty_initialIdx = (i:length(uniqueColors_noEmpty))';
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%



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


if onlyExistingAnnotationsForLabel == 1 & h_FigLegend ~=0
    groupColorForPlot_unique = unique(groupColorForPlot);
    groupColorForPlot_unique(groupColorForPlot_unique=="")=[];
    %[uniqueColors_noEmpty] = intersect(uniqueColors_noEmpty, groupColorForPlot_unique,'stable');
    [h_FigLegend_Final, h_FigLegend_axFinal] = editLegendEntriesLocal(h_FigLegend, uniqueColors_noEmpty,groupColorForPlot_unique)
end


% give the handle to the actual figure, not the legend
figure(h_fig1);

end

function [figFinal, axFinal] = editLegendEntriesLocal(figP, S, Sfinal)
%% NOTE: There is also a function editLegendEntries
%
% figP    : handle to figure containing P (or use gcf)
% S       : string array (Nx1) labels in original plot order
% Sfinal  : string array subset of S
%
% Returns: new figure + axes with the subset plotted using original marker/color.

disp('You may want to use addClassMedianLabel function')

S = string(S(:));
Sfinal = string(Sfinal(:));

% Map Sfinal -> indices in S (keeps Sfinal order)
notFoundS = setdiff(Sfinal, S);
[tf, idx] = ismember(Sfinal, S);
idx = idx(tf);
Skeep = Sfinal(tf);

% Get axes from original figure
axP = findall(figP, 'Type','axes');
axP = axP(1);  % if multiple axes, adjust as needed

% Get plotted objects (lines/scatters) from original axes
% We assume one plotted object per entry in S, in the same order they were created.
% objs = flipud(findall(axP, 'Type','line'));      % common for plot() points
% scat = flipud(findall(axP, 'Type','scatter'));   % common for scatter()
scat = (findall(axP, 'Type','scatter'));   % common for scatter()

if ~isempty(scat)
    % If scatter objects exist, prefer them
    objs = scat;
end

% Safety check
if numel(objs) < numel(S)
    warning('Found fewer plot objects (%d) than labels in S (%d). Order may not match.', numel(objs), numel(S));
end

% Extract styles for kept items
mk = strings(numel(idx),1);
col = zeros(numel(idx),3);

for i = 1:numel(idx)
    h = objs(idx(i));
    %h = objs((i));
    % Marker
    if isprop(h,'Marker') && ~isempty(h.Marker)
        mk(i) = string(h.Marker);
    else
        mk(i) = "o";
    end

    % Color: scatter uses MarkerEdgeColor/MarkerFaceColor; line uses Color
    col(i,:) = h.CData;
    %{
    if isprop(h,'Color')
        col(i,:) = h.Color;
    elseif isprop(h,'MarkerEdgeColor') && isnumeric(h.MarkerEdgeColor)
        col(i,:) = h.MarkerEdgeColor;
    else
        col(i,:) = [0 0 0];
    end
    %}
end

% Create final plot
figFinal = figure;
axFinal = axes(figFinal); hold(axFinal,'on');

x = ones(numel(idx),1);
y = (1:numel(idx))';
text(1,1,"   unknown",'FontSize',8)
plot(axFinal, 1, 1,'.k')
for i = 1:numel(idx)
    plot(axFinal, x(i), i+1, ...
        'LineStyle','none', ...
        'Marker', char(mk(end-i+1)), ...
        'MarkerSize', 6, ...
        'Color','k',...
        'MarkerFaceColor',col(end-i+1,:));
        %'faceColor', col(i,:));
end

text(x,flipud(y+1),strcat("   ",Skeep),'FontSize',8)
% Add labels on y-axis
%axFinal.YTick = y;
%axFinal.YTickLabel = cellstr(Skeep);
%axFinal.XLim = [0.5 1.5];
axFinal.Box = 'on';
xlim([0.9,2]); axFinal.YLim = [0 numel(idx)+2];
set(axFinal,'xtick',[]);set(axFinal,'ytick',[]);
figureTightened_andHighQualityTiffExport_RP_B(gcf)
end




