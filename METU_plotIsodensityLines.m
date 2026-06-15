function h_contour = MetU_plotIsodensityLines(hFigure, plottingData, bandwidthFactor)
% This function creates isodensity lines aroung a set of coordinates.
% Inputs
% hFigure: e.g., gcf
% plottingData: two column matrix with coordinates
% bandwidthFactor regulates tightness of lines (default = 0.3)

if nargin == 2
    bandwidthFactor = 0.3;         % <--- smaller = tighter contours
end

% x,y = UMAP coordinates (column vectors)
x = plottingData(:,1);
y = plottingData(:,2);

figure(hFigure); hold on
axis tight
expandAxesPadding(gca, 0.05)

%--- 1) Plot your coloured points  ---------------------------------------
% scatter(x, y, 15, 'b', 'filled');   % replace by your colouring
% colormap parula

%--- 2) Finer grid -------------------------------------------------------
nx = 400;          % was 200
ny = 400;


xlin = linspace(min(xlim), max(xlim), nx);
ylin = linspace(min(ylim), max(ylim), ny);
[Xg, Yg] = meshgrid(xlin, ylin);

%--- 3) 2-D kernel density with *tighter* bandwidth ----------------------
XY  = [x y];
pts = [Xg(:) Yg(:)];

% Get MATLAB’s default bandwidth, then shrink it
[~,~,bw] = ksdensity(XY);      % bw is a 1x2 vector

bw_tight = bw * bandwidthFactor;

f = ksdensity(XY, pts, 'Bandwidth', bw_tight);
F = reshape(f, size(Xg));

%--- 4) Contour levels concentrated on higher density --------------------
% ignore very low-density tail
low  = prctile(F(:), 20);
high = prctile(F(:), 99);
levels = linspace(low, high, 12);

h_contour = contour(Xg, Yg, F, levels, ...
        'LineColor', [0.5 0.5 0.5], ...
        'LineWidth', 0.3);

%--- 5) Cosmetics ---------------------------------------------------------
%axis equal
box on
%set(gca,'Layer','top')

end


function expandAxesPadding(ax, pct)
% expandAxesPadding(ax, pct)
% Adds pct (e.g., 0.05 for 5%) extra space to both sides of X and Y axes.
% ax  - axes handle (use gca if omitted)
% pct - fraction padding (default 0.05)

if nargin < 1 || isempty(ax), ax = gca; end
if nargin < 2 || isempty(pct), pct = 0.05; end

% --- X limits ---
xl = xlim(ax);
if isdatetime(xl)
    span = seconds(xl(2) - xl(1));
    if span == 0, span = 1; end
    pad  = pct * span;
    xlim(ax, [xl(1) - seconds(pad), xl(2) + seconds(pad)]);
else
    span = xl(2) - xl(1);
    if span == 0, span = max(abs(xl))*0.1 + 1; end
    pad  = pct * span;
    xlim(ax, [xl(1) - pad, xl(2) + pad]);
end

% --- Y limits ---
yl = ylim(ax);
if isdatetime(yl)
    span = seconds(yl(2) - yl(1));
    if span == 0, span = 1; end
    pad  = pct * span;
    ylim(ax, [yl(1) - seconds(pad), yl(2) + seconds(pad)]);
else
    span = yl(2) - yl(1);
    if span == 0, span = max(abs(yl))*0.1 + 1; end
    pad  = pct * span;
    ylim(ax, [yl(1) - pad, yl(2) + pad]);
end
end
