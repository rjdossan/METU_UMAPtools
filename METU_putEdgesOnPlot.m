%% Put edges on a graph

function H_fig = METU_putEdgesOnPlot(H_fig,plottingData,refIdx,targetIdx,lineTypeColour)
figure(H_fig); hold on
for metabEdge = 1:length(refIdx)
    
    plot([plottingData(refIdx(metabEdge),1);plottingData(targetIdx(metabEdge),1)],...
        [plottingData(refIdx(metabEdge),2);plottingData(targetIdx(metabEdge),2)],lineTypeColour);
% plot([plottingData(refIdx(metabEdge),1);plottingData(refIdx(metabEdge),2)],...
%         [plottingData(targetIdx(metabEdge),1);plottingData(targetIdx(metabEdge),2)],lineTypeColour);
    
end