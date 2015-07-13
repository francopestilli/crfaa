function plotBestHrf(varargin)
% function plotBestHrf(varargin)
%
% this function is created for crfaa, to use data at riken.
% to make a pretty plot for figure 2.
%
% it plots the hrf resulting from a GLM and a Deconvolution,
% analysis as performaed with fitTimecourse.m
%
% the deconvolution is plotted as data, with ste as errobars
% the GLM as the model fit.
%
% it can plot any number fo conditions among the existing ones.
%
% franco pestilli 2009.07.04
%
% best hrf result from these calls:
% close all, plotBestHrf('observer','fp','adaptation',0,'conditions',[25 23 24 28 29 30 31 32],'timeseries','v1')
% close all, plotBestHrf('observer','fp','adaptation',0,'conditions',[25 26 27 28 29 30 31 32],'timeseries','v1')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
conditions=[];                   %
observer=[];                     %
adaptation=[];                   %
timeseries=[];                   %
day=[];                          %
saveFig=[];                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set default variables:
getArgs(varargin,{...
 'conditions',[25 26 27 28 29 30 31 32], ...
 'observer','fp', ...
 'adaptation', 0, ...
 'timeseries','v1', ...
 'day','05-Jul-2009',...% old data set: '15-Feb-2009', ...
 'saveFig',1});

% if one of the original roi timeseries is requested switch to the
% otirginal date automatically
if strcmp(timeseries(1),'v')
 day = '29-Dec-2008'; 
end

% het the filename given, observer, adaptation, day etc... 
filename = getfilename(observer,adaptation,timeseries,day);

% loading a file, just a set one for the moment:
d1 = load(filename,'d1');
d1=d1.d1;

% strip out what we want to plot:
time = d1.time;

% results of deconvolution:
deconvResponses = d1.deconv.ehdr;
deconvErrors    = d1.deconv.ehdrste;

% glm results
glmResponses = d1.ehdr;

% make the plot:
plotHRFlocal(time,deconvResponses,deconvErrors, glmResponses, conditions, saveFig)
keyboard



%%%%%%%%%%%
% plotHrf %
%%%%%%%%%%%
function plotHRFlocal(time,deconvResponses,deconvErrors,glmResponses, conditions, savePlots)

h = smartfig('plotBestHrf','reuse');

ncontrasts = length(conditions);

% generates nice colors:
c = makeColor(ncontrasts,0);

% set formatting for the plots:
plotInfo = plotFormat;


% plot the deconvolution result
for i = 1:ncontrasts
 myerrorbar(time,deconvResponses(conditions(i),:), ...
  'yError',deconvErrors(conditions(i),:), ...
  'Symbol','o', ...
  'Color', c{i}, ...
  'MarkerFaceColor',c{i}, ...
  'MarkerEdgeColor','k', ...
  'MarkerSize',plotInfo.MarkerSize);
 hold on
end

% plot the glm result
for i = 1:ncontrasts
 
 % resample the glm response with a higher resolution
 resapledTime = linspace(min(time),max(time),100);
 r = interp1(time,glmResponses(conditions(i),:),resapledTime,'linear');
 
 myerrorbar(resapledTime,r, ...
  'Symbol','-', ...
  'LineWidth',plotInfo.LineSize,...
  'Color',c{i});
end



xlabel('Contrast (%)','FontSize',plotInfo.Fsize);
ylabel(sprintf('Response\n(% signal change)','%'),'FontSize',plotInfo.Fsize);

xlabel('Time (s)','FontSize',plotInfo.Fsize);
ylabel(sprintf('fMRI response\n(%s signal change)','%'),'FontSize',plotInfo.Fsize);


% figure formatting:
set(gca,...
 'FontName','Helvetica', ...
 'FontSize',plotInfo.Fsize, ...
 'PlotBoxAspectRatio',plotInfo.PlotBoxAspectRatio, ...
 'XLim', plotInfo.Xlim, ...
 'xScale','lin', ...
 'XTick', plotInfo.XTicks, ...
 'XTickLabel', plotInfo.XTicks ,...
 'XColor',plotInfo.XYColor, ...
 'YLim', plotInfo.YLim,...
 'yScale','lin', ...
 'YTick', plotInfo.YTicks, ...
 'YTickLabel', plotInfo.YTicksLabel, ...
 'YColor',plotInfo.XYColor, ...
 'LineWidth',plotInfo.LineWidth, ...
 'TickLength',plotInfo.TickLength, ...
 'Box','off');

myaxisHRF;

drawnow

% save figure
if savePlots
 set(h,'PaperPosition',[.25 .25 8 10.5]);
 set(h,'PaperOrientation','Portrait');
 
 figName    = 'bestHRF';
 dataFolder = '/Volumes/data/riken/crfaa/fmridata/crfaa/';
 figDir     = 'figureR1';

 disp(sprintf('Saving figure %s.eps',figName));
 savefig(h,'figName',figName,'defaultDataFolder',dataFolder,'figDir',figDir,'verbose',1);
end


  
%%%%%%%%%%%%%%
% plotFormat %
%%%%%%%%%%%%%%
function plotInfo = plotFormat

plotInfo.XYColor = 'k';
plotInfo.Fsize = 20;
plotInfo.LineWidth = 1.5;
plotInfo.PlotBoxAspectRatio = [1.5 1 1];
plotInfo.YLim = [-.25  1.25];
plotInfo.Xlim = [0.4 20.4];
plotInfo.YTicks = [-.25 0 0.25 0.5 0.75 1 1.25];
plotInfo.YTicksLabel = [-.25 0 0.25 0.5 0.75 1 1.25];
plotInfo.XTicks = [0 5 10 15 20];
plotInfo.XTicksLabel = [0 5 10 15 20 25];
plotInfo.MarkerSize = 8;
plotInfo.TickLength = [0.025 .01];
plotInfo.LineSize = 1;



%%%%%%%%%%%%%%
% makeColors %
%%%%%%%%%%%%%%
function c = makeColor(numC,dispColors)
% generate some nice colors:

for i = 1:numC
 c{i} = getSmoothColor(i,numC,'hsv');
 
 if dispColors
  plot(i,1,'ko','MarkerFaceColor',c{i},'MarkerSize',26);hold on
  text(i,1,sprintf('%i',i),'HorizontalAlignment','Center','Color',[0 0 0]);
 end
end


%%%%%%%%%%%%%%%%%
%  getFileName  %
%%%%%%%%%%%%%%%%%
function filename = getFileName(observer,visualarea,adaptation,day)
% set right r2 string
if strcmp(visualarea,'v1')
 r2 = '0_7';
else
 r2 = '0_5';
end

filename = deblank(sprintf('%s_A%s_%s_%s_%sroi.mat',upper(observer),num2str(adaptation),day,lower(visualarea),r2));



%%%%%%%%%%%%%%%%
%  getsession  %
%%%%%%%%%%%%%%%%
function session = getsession(observer,adaptation)

% check which computer are we using:
[notok c] = system('hostname');
if ~(~notok && strcmp('riken.jp',deblank(c(end-8:end))))
 % check if there is a drive connected which one is it ...
 if isdir('/Volumes/homosacer_data')
  hd = '/Volumes/homosacer_data';
 elseif isdir('/Volumes/homosacer')
  hd = '/Volumes/homosacer';
 elseif isdir('/Volumes/homosacer_backup')
  hd = '/Volumes/homosacer_backup';
 else
  disp('No data hard drive found.')
  keyboard
 end
elseif (~notok && strcmp('riken.jp',deblank(c(end-8:end))))
 hd = '/Volumes/data/riken/crfaa';
end

switch lower(observer)
 case {'fm'}
  switch adaptation
   case {0}
    mrsession = 'fm20080209';
   case 28
    mrsession = 'fm20080406';
   otherwise
    disp(sprintf('This condition is NOT set (%s, adp-%i).',observer,adaptation))
    mrsession = [];
  end
  
 case {'fp'}
  switch adaptation
   case {0}
    mrsession = 'fp20071019';
   case {28}
    mrsession = 'fp20080402';
   case {100}
    mrsession = 'fp20080415';
   otherwise
    disp(sprintf('This condition is NOT set (%s, adp-%i).',observer,adaptation))
    keyboard
  end
  
 case {'jg'}
  switch adaptation
   case {0}
    mrsession = 'jg20070919';
   case {28}
    mrsession = 'jg20080402';
   case {100}
    mrsession = 'jg20080414';
   otherwise
    disp(sprintf('This condition is NOT set (%s, adp-%i).',observer,adaptation))
    keyboard
  end
  
 otherwise
  disp(sprintf('Observer (%s) NOT found.',observer))
  keyboard
end

if ~isempty(mrsession)
 session = deblank(sprintf('%s/fmridata/crfaa/%s_A%i_testStream/%s',hd,upper(observer),adaptation,mrsession));
else
 session = '';
end
disp(sprintf('Loading: %s',session))



%%%%%%%%%%%%%%%%
% switchTSfile %
%%%%%%%%%%%%%%%%
function timeseries = switchTSfile(timeseries,day)

switch timeseries
 case 'v1'
  ts = 'v1_0_7roi';
 case 'v2'
  ts = 'v2_0_5roi';
 case 'v3'
  ts = 'v3_0_5roi';
 case 'v4'
  ts = 'v4_0_5roi';
 case 'combined' % restircted to a stricter r2
  ts = 'v1_TestsSubtractTSnoiseGroupsCombV1_3';
 case '025' % restircted to a stricter r2
  ts = 'v1_subTSnoiseSTD_combinedR2025';
 case '05' % restircted to a stricter r2
  ts = 'v1_subTSnoiseSTD_combinedR205';
 case '075' % restircted to a stricter r2
  ts = 'v1_subTSnoiseSTD_combinedR2075';
 case '085' % restircted to a stricter r2
  ts = 'v1_subTSnoiseSTD_combinedR2085';
 case '088' % restircted to a stricter r2
  ts = 'v1_subTSnoiseSTD_combinedR2088';

 otherwise
  keyboard
end

% add the day the file was created:
timeseries = sprintf('%s_%s',day,ts);



%%%%%%%%%%%%%%%
% getfilename %
%%%%%%%%%%%%%%%
function filename = getfilename(observer,adaptation,timeseries,day)

% get the session:
session = getsession(observer,adaptation);

% get the right filename:
ts = switchTSfile(timeseries,day);


% construct the filename
filename = sprintf('%s/%s',session,ts);


