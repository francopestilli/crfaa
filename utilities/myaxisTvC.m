% myaxisTvC.m
%
%      usage: myaxisTvC()
%         by: franco
%       date: 02/03/2009
%    purpose: put up publication quality axis
%             plot data on log log
%             to use this function specify 
%             exactly the x- and y-ticks like: 
%             set(gca, ...
%                'XScale','log', ...
%                'YScale','log', ...
%                'XTick', [.001 .01 .1 1], ...
%                'YTick', [.001 .01 .1 1], ...
%                'XLim',[.001 1], ...
%                'YLim',[.001 1]);
%
% developed to plot TvC functions.
function [xphyaxis yphyaxis] = myaxisTvC(putzero,multiplier)

if ~(strcmp(get(gca,'xscale'),'log') && strcmp(get(gca,'yscale'),'log'))
 % axis must be log-log for this function
 keyboard
end

if ieNotDefined('putzero')
 putzero = 1;
end

if ieNotDefined('multiplier')
 multiplier = 100;
end

xphyaxis = [];
yphyaxis = [];

% get the edgecolor, this should be
% either 'k' or white = [0.99 0.99 0.99]
% if it is set to 'w' it will be reversed
% by matlab to black when you export the figure
global edgecolor;
if (isempty(edgecolor))
 edgecolor = 'k';
end
if (isequal(edgecolor,[0.99 0.99 0.99]))
 set(gcf,'Color','k');
end

% this should turn off clipping, but doesn't seem to work
hold on, axis manual;
set(gcf,'Clipping','off');
set(gca,'Clipping','off');

axis off;

%%%%%%%%%%%%%%%%%%%%%%%%%
%% get axis info first %%
%%%%%%%%%%%%%%%%%%%%%%%%%
tickaspect = get(gca,'TickLength');
tickaspect = tickaspect(2)/tickaspect(1);

xTick = get(gca,'XTick');
xLim = get(gca,'XLim');
xaxissize = max(xTick)-min(xTick);

yTick = get(gca,'YTick');
yLim = get(gca,'YLim');
yaxissize = max(yTick)-min(yTick);

%%%%%%%%%%%%%%%%%%%%%%
%% deal with X-axis %%
%%%%%%%%%%%%%%%%%%%%%%
% get tick values and reset min max and axissize accordingly
xticks = suggesttick(min(xTick),max(xTick),xTick);
xlabelOffset = yLim(1)/10;

% make the X-axis a variable and draw them
xphyaxis.x = xLim(1);
xphyaxis.y = yLim(1);
xphyaxis.len = xaxissize;
xphyaxis.orient = 0;

% tick lenght depends on y-axis limits
xphyaxis.ticklen = yLim(1)/10;
xphyaxis.offset = yLim(1)/10;
xphyaxis.ticks = xticks;

% set labels
for i = 1:length(xticks.major)
 if putzero && i == 1 % plot 0 at first place
  xphyaxis.ticks.labels{i} = '0';
 elseif i == length(xticks.major) % do not plot 1
  xphyaxis.ticks.labels{i} = '';
 else
  xphyaxis.ticks.labels{i} = num2str(multiplier*xticks.major(i));
 end
end


%%%%%%%%%%%%%%%%%%%%%%
%% deal with Y-axis %%
%%%%%%%%%%%%%%%%%%%%%%
% get tick values and reset min max and axissize accordingly
yticks = suggesttick(min(yTick),max(yTick),yTick);
yticks.minor = [.004375 .00875 .0175 .035 .07 .14 .28 .56];
yticks.minortick = [];
yticks.major = [.004375 .00875 .0175 .035 .07 .14 .28 .56];
yticks.majortick = [];

ylabelOffset = xLim(1)/10;

% make the axis variable and draw them
yphyaxis.x = xLim(1);
yphyaxis.y = yLim(1);
yphyaxis.len = yaxissize;
yphyaxis.orient = 90;

% tick lenght depends on x-axis limits
yphyaxis.ticklen = xLim(1)/10;
yphyaxis.offset = xLim(1)/10;
yphyaxis.ticks = yticks;

% set labels
for i = 1:length(yticks.major)
  yphyaxis.ticks.labels{i} = num2str(multiplier*yticks.major(i));
end

% draw the axis now computed:

drawaxis(xphyaxis,edgecolor);
drawaxis(yphyaxis,edgecolor);

% set axis limits:
xaxis(xLim(1)*.2,xLim(2));
yaxis(yLim(1)*.5,yLim(2)+.06);


% turn back on the labels
set(get(gca,'XLabel'),'Visible','on');
set(get(gca,'YLabel'),'Visible','on');
set(get(gca,'XLabel'),'FontSize',14);
set(get(gca,'XLabel'),'FontName','Helvetica');
set(get(gca,'YLabel'),'FontSize',14);
set(get(gca,'YLabel'),'FontName','Helvetica');
set(get(gca,'Title'),'FontName','Helcetica');
set(get(gca,'Title'),'FontSize',14);
set(get(gca,'XLabel'),'Position',[xlabelOffset .5 1]);
set(get(gca,'YLabel'),'Position',[ylabelOffset .1 1]);
drawnow





