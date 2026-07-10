%% This script needs the following packages:
% HELP_RPtools
% umap

%% Load AW1 data with 500 samples and 4089 metabolomic features
% NOTE: It has been imputed, there are already no missing values; 

cd('C:\Users\rjdossan\OneDrive - Imperial College London\METU_UMAPtools\Data')
load('AW1_plasmaXCMS_500samples_imputed.mat')
X.VarInfo = HELP_tableCellToString(X.VarInfo);% make all columns string or numerical

%% Log scale and z-score the data by columns (metabolites) - OPTIONAL
% X.Data = log(X.Data+1);
% X.Data = zscore(X.Data);

%% Create UMAP
% NOTE: there was already a UMAP representation calculated previously, in
% X.UMAPres, but we will calculate UMAP again, for this example.

X.UMAPres = METU_UMAPregular(X.Data',round(0.01*size(X.VarInfo,1)),'cosine');
close

%% Plot UMAP with diagonal and circle

figure, plot(X.UMAPres.reduction(:,1), X.UMAPres.reduction(:,2),'.k')
xlabel('UMAP dimension 1'); xlabel('UMAP dimension 2');
METU_plotIsodensityLines(gcf, X.UMAPres.reduction); axis tight, grid on
METU_plotDiagonalORcircle(gcf); % plot diagonal
METU_plotDiagonalORcircle(gcf, X.UMAPres.reduction(1500,:),0.075); % plot circle at point number 1500

%% Plot UMAP coloured by RefMet classes

% Define annotations and classes
annotationsForPlot = X.VarInfo.AbbreviatedAnnotation;
groupColorForPlot = X.VarInfo.RefMet_Sub_class;

% Plot
figure, plot(X.UMAPres.reduction(:,1), X.UMAPres.reduction(:,2),'.k')
METU_plotIsodensityLines(gcf, X.UMAPres.reduction); axis tight, grid on
[h_figFinal,h_FigLegend] = METU_scatterPlot_AnnotColor_grpLabel(gcf,X.UMAPres.reduction,annotationsForPlot,groupColorForPlot,'jet',1,45,8)
xlabel('UMAP dimension 1'); xlabel('UMAP dimension 2');