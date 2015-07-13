function plotBestHrf_R1(varargin)
% function plotBestHrf_R1(varargin)
%
% this function is created for crfaa, to use data at columbia for 
% NN Rebuttal.
%
% it plots the hrf resulting from a GLM and a Deconvolution,
% analysis as performaed with fitTimecourse.m
%
% the deconvolution is plotted as data, with ste as errobars
% the GLM as the model fit.
%
% it can plot any number fo cueCondition among the existing ones.
%
% franco pestilli 2010.06.03
%
% best hrf result from these calls:

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cueCondition=[];                 %
observer=[];                     %
adaptation=[];                   %
visualArea=[];                   %
day=[];                          %
saveFig=[];                      %
subtractionType=[];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set default variables:
getArgs(varargin,{...
 'cueCondition','dt', ...
 'observer','fp',    ...
 'adaptation', 0,    ...
 'visualArea','v1',  ...
 'day','2010-30-May',...
 'saveFig',1,        ...
 'r2','070',         ...
 'subtractionType', 'ORIG'});

% get the filename given, observer, adaptation, day etc... 
defaultDataFolder = '/data2/crfaa/crfaa/data_used_files_R1/';
[filename] = makeFileName(observer,adaptation,visualArea,defaultDataFolder, r2, subtractionType);
 
% loading a file, just a set one for the moment:
disp(sprintf('[%s] Loading %s...',mfilename,filename))

switch subtractionType
 case {'ORIG'} % load d1 and show the fitTimecourse result
  d1 = load(filename, 'd1');
  d1=d1.d1;
  
  % strip out what we want to plot:
  time = d1.time;
  
  % results of deconvolution:
  deconvResponses = d1.deconv.ehdr;
  deconvErrors    = d1.deconv.ehdrste;
  
  % glm results
  glmResponses = d1.ehdr;
  
 case {'BYTRIAL'} % load d
  d = load(filename, 'd');
  d=d.d.BeforeSubtraction_d;
  
  % strip out what we want to plot:
  time = d.time;
  
  % results of deconvolution:
  deconvResponses = d.ehdr;
  deconvErrors    = d.ehdrste;
  
  % glm results
  glmResponses = [];
 otherwise
  keyboard
end



% get the cue condition
[condition colorIndex] = getcondition(cueCondition);

% make the plot:
plotHRFlocal(time,deconvResponses,deconvErrors, glmResponses, condition, saveFig, ...
             observer,defaultDataFolder,r2,visualArea,subtractionType,cueCondition, colorIndex)
return



%%%%%%%%%%%
% plotHrf %
%%%%%%%%%%%
function plotHRFlocal(time,deconvResponses,deconvErrors,glmResponses, condition, savePlots, ...
                      observer,defaultDataFolder,r2,visualArea,subtractionType, cueCondition, colorIndex)

h = smartfig('plotBestHrf_R1','reuse');

ncontrasts = length(condition);

% generates nice colors:
% c = makeColor(ncontrasts,0);

% set formatting for the plots:
plotInfo = plotFormat;

colorScaleFactor = linspace(.38,1,8);


% plot the deconvolution result
for i = 1:ncontrasts
 myerrorbar(time,deconvResponses(condition(i),:), ...
  'yError',deconvErrors(condition(i),:), ...
  'Symbol','o', ...
  'Color', plotInfo.thisColor{1}{colorIndex}*colorScaleFactor(i), ...
  'MarkerFaceColor',plotInfo.thisColor{1}{colorIndex}*colorScaleFactor(i), ...
  'MarkerEdgeColor',plotInfo.MarkerEdgeColor, ...
  'MarkerSize',plotInfo.MarkerSize);
 hold on
end

% plot the glm result
for i = 1:ncontrasts
 
 % resample the glm response with a higher resolution
 resapledTime = linspace(min(time),max(time),100);
 r = interp1(time,glmResponses(condition(i),:),resapledTime,'linear');
 
 myerrorbar(resapledTime,r, ...
  'Symbol','-', ...
  'LineWidth',plotInfo.LineSize,...
  'Color',plotInfo.thisColor{1}{colorIndex}*colorScaleFactor(i));
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
 
 figName    = sprintf('%s_%s_%s_%s_%s',observer, visualArea, r2, subtractionType,cueCondition);
 figDir     = 'fig_best_hrf';

 disp(sprintf('Saving figure %s.eps',figName));
 savefig(h,'figName',figName,'defaultDataFolder',defaultDataFolder,'figDir',figDir,'verbose',1);
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

% generate some nice colors:
numC =60;
for i = 1:numC
 c{i} = getSmoothColor(i,numC,'hsv');
 %  plot(i,1,'ko','MarkerFaceColor',c{i},'MarkerSize',26);hold on
 %  text(i,1,sprintf('%i',i),'HorizontalAlignment','Center','Color',[0 0 0]);
end

% choose the ones i like:
reds = {c{2} c{5} c{7} c{11}};
greens = {c{23} c{17} c{15} c{13}};
blues = {c{40} c{37} c{34} c{30}};
purples = {c{46} c{49} c{54} c{50}};

for i = 1:4 % visual area
 plotInfo.thisColor{i} = {reds{i} blues{i} purples{i} greens{i}};
 plotInfo.MarkerFaceColor{i} = plotInfo.thisColor{i};
end

plotInfo.MarkerEdgeColor = [.4 .4 .4];



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


%%%%%%%%%%%%%%%%
% makeFileName %
%%%%%%%%%%%%%%%%
function [filename vars2load] = makeFileName(observer,adaptation,visualArea,defaultDataFolder, r2, subtractionType)
% NB whereas the main analysis files are saved into the special folder
% <data_used_files> so that changes to that folder do not screw up the
% files...
% the files for the other analyses are instad in each folder/session. so
% for construtcting the filename we use getsession.

if strcmpi(observer,'avrg')
 keyboard
 % not used
 
else
 % choose the type of analysis files to load
 vars2load = {'d1' 'events'};
 this_date = '2010-30-May';
 if (~strcmp(visualArea,'v1') && strcmp(r2,'050')) || (strcmp(visualArea,'v1') && strcmp(r2,'070'))
  filename = fullfile(defaultDataFolder,sprintf('%s_A%s_%s_%s_era_R1_subtr%s_%s.mat',upper(observer),num2str(adaptation),this_date,visualArea, ...
   subtractionType,visualArea));
  
 else
  filename = fullfile(defaultDataFolder,sprintf('%s_A%s_%s_%s_era_R1_subtr%s_%sr%s.mat',upper(observer),num2str(adaptation),this_date,visualArea, ...
   subtractionType,visualArea,r2));
 end
end


%%%%%%%%%%%%%%%%
% getcondition %
%%%%%%%%%%%%%%%%
function [condition colorIndex] = getcondition(cueCondition)
switch cueCondition
 case {'dnt', 'distributed-non-target'}
  condition = [25 26 27 28 29 30 31 32];
  colorIndex = 3;
 case {'dt', 'distributed-target'}
  condition = [17 18 19 20 21 22 23 24];
  colorIndex = 2;
 case {'fnt' 'focal-non-target'}
  condition = [ 9 10 11 12 13 14 15 16];
  colorIndex = 4;
 case {'ft' 'focal-target'}
  condition = [1 2 3 4 5 6 7 8];
  colorIndex = 1;
 otherwise
  keyboard
end