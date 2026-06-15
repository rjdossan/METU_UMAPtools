% UMAPregular
% Calculates simple UMAP with nNeighbours and defined metric.
% UMAPres = UMAPregular(UMAPdata, nNeighbours,metricType)
% UMAPres = UMAPregular(UMAPdata, nNeighbours) metric is cosine
% UMAPres = UMAPregular(UMAPdata) nNeighbours is 1% nFeatures, metric is cosine

function UMAPres = METU_UMAPregular(UMAPdata, nNeighbours,metricType)

if nargin<2
    nNeighbours = round(0.01*size(UMAPdata,1));
    metricType = 'cosine';
elseif nargin < 3
    metricType = 'cosine';
end
rng(3) 
UMAPres.opt.U_neighbors=nNeighbours;
UMAPres.opt.U_components = 2;
UMAPres.opt.U_metric = metricType;
%UMAPres.opt.U_metric = 'cosine'
%UMAPres.opt.U_metric = 'euclidean'
% UMAPres.opt.U_metric = 'correlation'
UMAPres.opt.U_epochs = 1000;
UMAPres.opt.min_dist = 0.3;% default is 0.3
[UMAPres.reduction,UMAPres.Xumap] = run_umap(UMAPdata,...
    'n_neighbors',UMAPres.opt.U_neighbors,...
    'n_components',UMAPres.opt.U_components,...
    'metric',UMAPres.opt.U_metric,...
    'min_dist',UMAPres.opt.min_dist,...
    'n_epochs',UMAPres.opt.U_epochs);
UMAPres.h_figUMAP = gcf;
%UMAPres.plottingData = UMAPres.reduction;

disp('To plot with labels and colours')
disp('annotationsForPlotProt = Xtomatch.VarInfo.AbbreviatedAnnotation;')
disp('groupColorForPlotProt = Xtomatch.VarInfo.RefMet_Sub_class;')
disp("[h_figUMAPannotatedProt,h_figUMAPannotated_legendProt] = METU_scatterPlot_AnnotColor_grpLabel(h_figUMAP_Prot,plottingDataProt,annotationsForPlotProt,groupColorForPlotProt,'jet');")
