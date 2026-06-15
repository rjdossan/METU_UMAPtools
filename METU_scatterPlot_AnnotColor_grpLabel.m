%% function to color two dimensional maps
% [h_figFinal,h_FigLegend] = scatterPlot_AnnotColor(h_fig,plottingData,annotationsForPlot,groupColorForPlot,colormapName,transparencyValue,markerSize_value,legend_font_size)
function [h_figFinal,h_FigLegend] = MetU_scatterPlot_AnnotColor_grpLabel(h_fig,plottingData,annotationsForPlot,groupColorForPlot,colormapName,transparencyValue,markerSize_value,legend_font_size)

% INPUTS:
% h_fig,plottingData,annotationsForPlot,groupColorForPlot,colormapName,transparencyValue,markerSize_value,legend_font_size

%% PLOT 3: plot with class colors

% colormaps that work well: jet, hsv, colorcube, prism
colormap('default');
if nargin <=4
    colormapPlot = colormap(hsv); colormapPlot(end-20:end,:)=[];
    transparencyValue = 1;
    markerSize_value = 50;
    markSize = 50*ones(size(plottingData,1),1);
    legend_font_size = 11;
elseif nargin == 5
    colormapPlot = colormap(colormapName); colormapPlot(end-20:end,:)=[];
    transparencyValue=1;
    markerSize_value = 50;
    markSize = 50*ones(size(plottingData,1),1);
    legend_font_size = 11;
elseif nargin == 6
    colormapPlot = colormap(colormapName); %colormapPlot(end-20:end,:)=[];
    markerSize_value = 50;
    markSize = 50*ones(size(plottingData,1),1);
    legend_font_size = 11;
elseif nargin == 7
    colormapPlot = colormap(colormapName); %colormapPlot(end-20:end,:)=[];
    markSize = markerSize_value*ones(size(plottingData,1),1);
    legend_font_size = 11;
elseif nargin == 8
    colormapPlot = colormap(colormapName); %colormapPlot(end-20:end,:)=[];
    markSize = markerSize_value*ones(size(plottingData,1),1);
    %legend_font_size = 11;
end

%% Run the following to plot

% get the unique classes
if string(class(groupColorForPlot)) == 'double'
    uniqueClass= unique(groupColorForPlot(~isnan(groupColorForPlot)));
else 
    groupColorForPlot(ismissing(groupColorForPlot)) = "";
    uniqueClass= unique(groupColorForPlot,'stable');
end

% get the colours for each class
annotationsForNet.color_idx = M2S_find_idxInReference(groupColorForPlot,uniqueClass);
if length(unique(annotationsForNet.color_idx))==1
    annotationsForNet.color_idx256 = ones(length(annotationsForNet.color_idx),1);
else
    annotationsForNet.color_idx256 = round((((size(colormapPlot,1)-1)/(length(uniqueClass)-1))) * (annotationsForNet.color_idx-1))+1;
end
annotationsForNet.color = colormapPlot(annotationsForNet.color_idx256,:);


% Define which points will have small dots (e.g. without annotation, or
% without groupClass). It should be the ones without groupClass

sizeByClassOrAnnotation = 'groupColor'; % {'annot','groupColor'}
if strcmp(sizeByClassOrAnnotation, 'annot')    
    smallDotsIdx = find(annotationsForPlot == "");    
elseif strcmp(sizeByClassOrAnnotation, 'groupColor')
    if string(class(groupColorForPlot)) == 'double'
        smallDotsIdx = find(groupColorForPlot == 0);
    else
        smallDotsIdx = find(groupColorForPlot == "");
    end
end   
if ~isempty(smallDotsIdx)
markSize(smallDotsIdx) = 5;% unknown annotation are small dots
annotationsForNet.color(smallDotsIdx,:) = repmat([0,0,0],length(smallDotsIdx),1);% black dots for undefined class   
end


%% Main figure with objects of different shapes

possibleShapes = ('osdvp^h')';
uniqueClass_noEmpty = uniqueClass;
uniqueClass_noEmpty(uniqueClass_noEmpty=="") =[];
nRepmatShapes = length(uniqueClass_noEmpty)/length(possibleShapes);
uniqueShapes = repmat(possibleShapes,ceil(nRepmatShapes),1);
uniqueShapes = uniqueShapes(1:length(uniqueClass_noEmpty));
h_figFinal = HELP_duplicateFigure(h_fig); hold on;
if sum(groupColorForPlot=="")>0
    scatter(plottingData(groupColorForPlot=="",1),plottingData(groupColorForPlot=="",2),markSize(groupColorForPlot==""),annotationsForNet.color(groupColorForPlot=="",:),'filled');
end
for c = 1:length(uniqueClass_noEmpty)
    tempIdx = groupColorForPlot == uniqueClass_noEmpty(c);
    scatter(plottingData(tempIdx,1),plottingData(tempIdx,2),markSize(tempIdx),annotationsForNet.color(tempIdx,:),uniqueShapes(c),'filled','MarkerEdgeColor','k','MarkerEdgeAlpha',0.5,'MarkerFaceAlpha',transparencyValue); 
end
grid on;
set(gca,'YAxisLocation','right');
h_axisScatter=gca;


%% Create a legend 

groupColorForPlot_noSmallDots_idx = setdiff((1:length(groupColorForPlot))',smallDotsIdx);
[uniqueLegendLabel,uniqueLegendLabel_idx] = unique(groupColorForPlot(groupColorForPlot_noSmallDots_idx),'stable');
if string(class(uniqueLegendLabel)) == "double"
    [uniqueLegendLabel,sortedIdx] = sort(uniqueLegendLabel);
    uniqueLegendLabel_idx = uniqueLegendLabel_idx(sortedIdx);
    uniqueLegendLabel = string(uniqueLegendLabel);
end
legendTable = table(uniqueLegendLabel,annotationsForNet.color(groupColorForPlot_noSmallDots_idx(uniqueLegendLabel_idx),:),'VariableNames',{'label','color'});

if ~isempty(legendTable)

    axes('Position',[0.01 0.11 0.16 0.817]); box on;
    h_FigLegend = scatter(1,1,markerSize_value,legendTable.color(1,:),uniqueShapes(1),'filled','MarkerEdgeColor','k'); hold on;
    if length(uniqueClass_noEmpty)>1
        for c=1:length(uniqueClass_noEmpty)
            scatter(1,c,markerSize_value,legendTable.color(c,:),uniqueShapes(c),'filled','MarkerEdgeColor','k');
        end
    end
    % text(ones(length(uniqueLegendLabel_idx),1),(1:length(uniqueLegendLabel_idx))',...
    %     uniqueLegendLabel,'fontsize',8,'verticalalignment','middle','interpreter','none');
    text(ones(length(uniqueLegendLabel_idx),1),(1:length(uniqueLegendLabel_idx))',...
        strcat("  ",uniqueLegendLabel),'fontsize',legend_font_size,'verticalalignment','middle','interpreter','none');
    set(gca,'xticklabels',[]);set(gca,'yticklabels',[]);
    set(gca,'xtick',[]);set(gca,'ytick',[]);
    xlim([0.9,1.6]);ylim([0,length(uniqueLegendLabel_idx)+1]);
else
    h_FigLegend = 0;
end

% give back the handle to the main figure
figure(h_figFinal);

%% add labels with visibility off

axes(h_axisScatter); % get the axes of the main figure
hT=text(plottingData(:,1),plottingData(:,2),annotationsForPlot,'FontSize',legend_font_size,'interpreter','none','VerticalAlignment','top','visible','off');
% Create button to turn on/off  annotations
uicontrol('Parent',h_figFinal,'Style','pushbutton','String','Annotation on/off',...
        'Units','normalized','Position',[0.75 0.01 0.2 0.05],'Visible','on','Callback',{@functionForFig_addText,hT});


%% add labels for median of each group with visibility off 

axes(h_axisScatter); % get the axes of the main figure
[classLabelsUnique, classMedianVals, tH] = METU_addClassMedianLabel(gcf,plottingData,groupColorForPlot,[0.45, 0.45, 0.45],round(1.3*legend_font_size),'median','off');
% Create button to turn on/off  annotations
uicontrol('Parent',h_figFinal,'Style','pushbutton','String','Group labels on/off',...
        'Units','normalized','Position',[0.05 0.01 0.2 0.05],'Visible','on','Callback',{@functionForFig_addText,tH});

end

% create the ButtonH figure
function functionForFig_addText(source, event, handlesText);
        length(handlesText);
        currentOnOff = (get(handlesText,'Visible'));
        if strcmp(currentOnOff{1},'off')
            set(handlesText,'visible','on');
            disp('Annotations ON');
        else
            set(handlesText,'visible','off');
            disp('Annotations OFF');
        end        
end