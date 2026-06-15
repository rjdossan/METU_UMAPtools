%% Function to plot two things:
% - A diagonal line with 10%s of diagonal marked in it 
% OR
% - A circle with a defined radius on a selected point around it
%
% Input example
%
% For diagonal line only
% [diagLength, h] = M2S_plotPercentDiagonal(gcf)
% For circle only
% [diagLength, h] = M2S_plotPercentDiagonal(gcf, xy, diagPercent)
% in which:
% xy = [2.4, 3.8];
% diagPercent = 0.1;

function [diagLength, h] = METU_plotDiagonalORcircle(handleFig, xy, diagPercent);
disp('In function METU_plotDiagonalORcircle the inputs xy diagPercent changed order')


if nargin  == 2
    diagPercent = 0.1;
end

figure(handleFig);
axis tight

% Get the coordinates from the plot
ax = get(gcf,'children');
[xData, yData] = HELP_getDotsFromAxes(ax);
xData = xData';
yData = yData';

xyLim = [min(xData) max(xData) min(yData) max(yData)];

diagLength = sqrt((xyLim(2)-xyLim(1))^2+(xyLim(4)-xyLim(3))^2);


if nargin == 1 
    xPlot = 0.1*(0:1:10)*(xyLim(2)-xyLim(1))+xyLim(1);
    yPlot = 0.1*(0:1:10)*(xyLim(4)-xyLim(3))+xyLim(3);
    hold on, plot(xPlot,yPlot,'k.-');
    text(xPlot(2:10), yPlot(2:10), strcat(string((10:10:90)'),"%"),'fontsize',7);
    h = [];
else
    x = xy(1);
    y = xy(2);
    h = XplotCircle(x, y, diagPercent * diagLength,'k');
end


function h = XplotCircle(x,y,r,lineColour);
hold on
th = 0:pi/50:2*pi;
xunit = r * cos(th) + x;
yunit = r * sin(th) + y;
h.circle = plot(xunit, yunit,lineColour);
h.centre = plot(x, y,[lineColour,'o']);
hold off