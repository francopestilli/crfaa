function plotCRF(varargin)
% function plotCRF(varargin)
%
% this function is created for crfaa, to use data at riken.
% to plot a crf from any analysis file in any roi or observer.
%
% it loads a file e.g.,
% 23-Jul-2009_v1_firstSecondIntervalAnalysis_STD_v1.mat
% runAllAnalyses.m and meanERA_ROI_remove0_script3.m
%
%
% it can plot any number fo conditions among the existing ones.
% if you need to change or add an alanysis set up day and filename in the
% subfunction: 'getfilename'
% 
% franco pestilli 2009.07.26
%
% best hrf result from these calls:
% close all, plotCRF('observer','fp','adaptation',0,'conditions',[25 23 24 28 29 30 31 32],'timeseries','v1')
% close all, plotCRF('observer','fm','adaptation',0,'conditions',[25 26 27 28 29 30 31 32],'timeseries','v1ci')
% close all, plotCRF('observer','fm','adaptation',0,'conditions',[25 26 27 28 29 30 31 32],'timeseries','v112')

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
 'conditions',{[1 2 3 4 5 6 7 8] ...
               [9 10 11 12 13 14 15 16]  ...
               [17 18 19 20 21 22 23 24] ...
               [25 26 27 28 29 30 31 32] ...
               [33 34 35 36 37 38 39 40] ...
               [41 42 43 44 45 46 47 48]}, ...
 'observer','jg', ...
 'adaptation', 0, ...
 'timeseries','v412', ...
 'saveFig',0});


% get the filename given, observer, adaptation, day etc... 
filename = getfilename(observer,adaptation,timeseries);

% loading a file, just a set one for the moment:
disp(sprintf('[plotCRF] loading file: <%s>',filename))
d1 = load(filename,'d1');
d1 = d1.d1;

% strip out what we want to plot:
d.contrasts = [0.875 1.75 3.5 7 14 28 56 84]./100;

% generates nice colors:
color = makeColor(length(conditions),0);

% extract the responses, fit them with a naka-rushton and plot them:
for c = 1:length(conditions)
 if max(conditions{c}) <= length(d1.amplitude)
  % glm results
  d.glmResponses = d1.amplitude(conditions{c});
  d.glmErrors = d1.amplitudeSTE(conditions{c});
  
  % fit a nakarushton to the data:
  crfFittype = 'naka';
  testType   = '1111';
  crf = fitCRF(d,crfFittype, testType);
  
  % make the plot:
  h = plotCRFlocal(d.contrasts,d.glmResponses,d.glmErrors, crf.fit.fitx, crf.fit.fity, color{c});
  
  % make a legend
  L{c} = sprintf('Condition: %i',c);
 end
end

% % make a legend NOT WORKING
% % 1. get the axis
% ax = get(h,'CurrentAxes');
% % 2. get children and their properties
% ch = get(ax,'Child');
% sy_properties = get(ch);
% c = 0;
% % 3. save the indexes of the data:
% for i = 1:length(sy_properties)
%  if strcmp(sy_properties(i).Marker,'o')
%   c = c + 1;
%   legendIndex{c} = ch(i);
%  end
% end
% 
% legend(ax,'String',L);

% save figure
if saveFig
 set(h,'PaperPosition',[.25 .25 8 10.5]);
 set(h,'PaperOrientation','Portrait');
 
 figName    = 'bestCRF';
 dataFolder = '/Volumes/data/riken/crfaa/fmridata/crfaa/';
 figDir     = 'fig_plotCRF';
 
 disp(sprintf('Saving figure %s.eps',figName));
 savefig(h,'figName',figName,'defaultDataFolder',dataFolder,'figDir',figDir,'verbose',1);
end



%%%%%%%%%%%%%%%%
% plotCRFlocal %
%%%%%%%%%%%%%%%%
function h = plotCRFlocal(contrast,deconvResponses,deconvErrors, fitx, glmResponses, c)

h = smartfig('plotCRF','reuse');

% set formatting for the plots:
plotInfo = plotFormat;


% plot the deconvolution result
myerrorbar(contrast,deconvResponses, ...
 'yError',deconvErrors, ...
 'Symbol','o', ...
 'Color', c, ...
 'MarkerFaceColor',c, ...
 'MarkerEdgeColor','k', ...
 'MarkerSize',plotInfo.MarkerSize);
hold on

% plot the glm result
[x index] = find(fitx>.00875);

myerrorbar(fitx(index),glmResponses(index), ...
 'Symbol','-', ...
 'LineWidth',plotInfo.LineSize,...
 'Color',c);

xlabel('Contrast (%)','FontSize',plotInfo.Fsize);
ylabel(sprintf('Response\n(% signal change)','%'),'FontSize',plotInfo.Fsize);

xlabel('Time (s)','FontSize',plotInfo.Fsize);
ylabel(sprintf('fMRI response\n(%s signal change)','%'),'FontSize',plotInfo.Fsize);

axis('square');

% figure formatting:
set(gca,...
 'FontName','Helvetica', ...
 'FontSize',plotInfo.Fsize, ...
 'PlotBoxAspectRatio',plotInfo.PlotBoxAspectRatio, ...
 'XLim', plotInfo.Xlim, ...
 'xScale','log', ...
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

drawnow


%%%%%%%%%%%%%%
%   fitCRF   %
%%%%%%%%%%%%%%
function crf = fitCRF(d,crfFittype, testType)
disp(sprintf('[plotCRF:fitCRF] start...'));

switch crfFittype
 case {'naka' 'NAKA' 'nk'}
  % fit a naka-rushton sigmodal function:
  crf.fit = fitsigmoidtest(d.contrasts,d.glmResponses, testType);
  
 case {'exp' 'EXP' 'e'}
  % this is the fuction i need to fit an exponential
  % if neededget this function from plotBestCRF.m
  crf.fit = fitexptest(d.contrasts,d.glmResponses, testType);

 otherwise
  disp(sprintf('[plotCRF:fitCRF] %s fit type not defined',crfFittype))
  keyboard
end

% set the fittype
% crf.fit.type = crfFittype;
disp(sprintf('[plotCRF:fitCRF] done'));


%%%%%%%%%%%%%%%%%%%%
%  fitsigmoidtest  %
%%%%%%%%%%%%%%%%%%%%
function bestfit = fitsigmoidtest(contrasts,responses, testType)

% check arguments
if (nargin ~= 3)
 help fitsigmoidtest
 bestfit = nan;
 return
end

% make sure we have a column vector
if (size(contrasts,1) == 1)
 contrasts = contrasts';
end
if (size(responses,1) == 1)
 responses = responses';
end

% set initial params
numConditions = size(responses,2);

[initparams minfit maxfit] = switchInitialParamsNK(testType,numConditions);


% set optimization parameters
optimizationParams = optimset( ...
 'LevenbergMarquardt','on', ...
 'MaxIter',inf, ...
 'TolFun',10^-10, ...
 'MaxFunEvals', 1000);

global numIters;
numIters = 0;

% fit function using lsqnonlin in LevenbergMarquardt mode.
[fitparams resnorm residual exitflag output lambda jacobian] = ...
 lsqnonlin(@sigmoiderr,initparams,minfit,maxfit,optimizationParams, ...
 contrasts,responses,testType,numConditions);

% Taken from Numerical Recipies,
% the leastsq function seems to return the transposed gradient
% instead of the jacobian...
hessian = jacobian'*jacobian;
reducedChiSquared = (residual'*residual)/(numel(responses)-length(initparams));
covar = (reducedChiSquared * inv(hessian));
M = max(diag(covar)); % get the maximum value on the
covar = covar/M; % normalized covariance matrix

% set params to the size of the conditions (e.g., 3*4*4):
[Rmax c50 n offset] = setNakaParams(fitparams,testType,numConditions);

[bestfit.err bestfit.yfit] = sigmoiderr(fitparams,contrasts,responses,testType,numConditions);

% compute R2 values and save them in the fit:
% (1) reshape the errors by condition and contrast
bestfit.err = reshape(bestfit.err,size(responses));

% (2) run through conditions and compute r2
% for each one of them
thisError = bestfit.err;
thisResponses = responses;
bestfit.r2 = 1-var(thisError)/var(thisResponses);

% if we have the best fit then keep it.
bestfit.params = fitparams;
bestfit.eachparam.Rmax = Rmax;
bestfit.eachparam.c50 = c50;
bestfit.eachparam.n = n;
bestfit.eachparam.offset = offset;
bestfit.covar = covar;
bestfit.output = output;
bestfit.testType = testType;
bestfit.numConditions = numConditions;

% compute the function returned by the best fit:
bestfit.fitx = .0001:.0001:1;

% compute th current model estimate:
[bestfit.fity Rmax c50 n offset] = testnaka(bestfit.params,bestfit.fitx,testType,numConditions);


%%%%%%%%%%%%
% testnaka %
%%%%%%%%%%%%
function [fitfun Rmax c50 n offset] = testnaka(fitparams,contrasts,testType,numConditions)

% set params of the naka rushton depending on the type of test being
% performed:
[Rmax c50 n offset] = setNakaParams(fitparams,testType,numConditions);


% calculate function
thisRmax = Rmax;
thisc50 = c50;
thisn = n;
thisoffset = offset;
thiscontrasts = contrasts;
% now compute responses for this CRF
fitfun = thisRmax*(thiscontrasts.^thisn)./(thiscontrasts.^thisn+thisc50.^thisn)+thisoffset;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function for fitting naka-rushton to crf %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [err fitfun] = sigmoiderr(fitparams,contrasts,responses,testType,numConditions)

% compute th current model estimate:
[fitfun Rmax c50 n offset] = testnaka(fitparams,contrasts,testType,numConditions);

% calculate error:
err = responses-fitfun;

% display the current fit (FIX: This doesn't yet work for displaying simultaneous fit)
global numIters 
numIters = numIters+1;
if 0
 % display the fit as we go along
 if ~isnan(numIters)
  cla
  % extract the current values:
  thisC = 100*contrasts;
  thisR = responses;
  thisFitFun = fitfun;
  this_c50 = 100*c50;
  this_Rmax = Rmax;
  this_n = n;
  this_offset = offset;
  [sortcontrasts sortindex] = sort(thisC);
  
  plot(thisC,thisR,'ko');hold on
  plot(sortcontrasts,thisFitFun(sortindex),'r-');
  vline(this_c50);
  title(sprintf('Rmax=%0.2f,\n c50=%0.2f,\n n=%0.2f,\n offset=%0.2f\n numIters=%i', ...
   this_Rmax,this_c50,this_n,this_offset,numIters));
  axis([1 100 0 1.5])
  drawnow
 end
end


%%%%%%%%%%%%%%%%%
% setNakaParams %
%%%%%%%%%%%%%%%%%
function [Rmax c50 n offset] = setNakaParams(fitparams,testType,numConditions)
% init values
Rmax   = zeros(numConditions);
c50    = zeros(numConditions);
n      = zeros(numConditions);
offset = zeros(numConditions);

switch testType
 case {'1111'} % fit all params
  % make fitparams back into a nAttention x nAdaptation x nParams matrix;
  fitparams = reshape(fitparams,[numConditions 4]);
  % and grab the values
  Rmax   = fitparams(1);
  c50    = fitparams(2);
  n      = fitparams(3);
  offset = fitparams(4);

 otherwise
  keyboard
end


%%%%%%%%%%%%%%%%%%%%%%%%%
% switchInitialParamsNK %
%%%%%%%%%%%%%%%%%%%%%%%%%
function [initparams minfit maxfit] = switchInitialParamsNK(testType,numConditions)
% Rmax = fitparams(1);
% c50 = fitparams(2);
% n = fitparams(3);
% offset = fitparams(4);

switch testType
 case {'1111'} % fit all params
  % init the params to 0
  initparams = zeros([numConditions 4]);
  minfit = zeros([numConditions 4]);
  maxfit = zeros([numConditions 4]);
  % Rmax
  initparams(1) = 18.0;
  minfit(1) = .5;
  maxfit(1) = inf;
  % c50
  initparams(2) = 5000;
  minfit(2) = 0;
  maxfit(2) = inf;
  % n
  initparams(3) = .5;
  minfit(3) = 0;
  maxfit(3) = inf;
  % offset
  initparams(4) = .15;
  minfit(4) = 0;
  maxfit(4) = inf;
  
 otherwise
  keyboard
end

  
%%%%%%%%%%%%%%
% plotFormat %
%%%%%%%%%%%%%%
function plotInfo = plotFormat

plotInfo.XYColor = 'k';
plotInfo.Fsize = 20;
plotInfo.LineWidth = .5;
plotInfo.PlotBoxAspectRatio = [1.5 1 1];
plotInfo.YLim = [0 1.25];
plotInfo.YTicks = [0 0.25 .5 .75 1 1.25];
plotInfo.YTicksLabel = [0 0.25 .5 .75 1 1.25];
plotInfo.Xlim = [.00875  1];
plotInfo.XTicks = [.00875 .0175 .035 .07 .14 .28 .56 .84 1];
plotInfo.XTicksLabel = [.00875 .0175 .035 .07 .14 .28 .56 .84 1]*100;
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
function timeseries = switchTSfile(timeseries)

% set up day of analysis
switch timeseries
 case {'v1' 'v2' 'v3' 'v4'}
  day = '29-Dec-2008';
  
 case {'combined' '025' '05' '075' '085' '088'}
  day = '05-Jul-2009';
  
 case {'v1ci' 'v2ci' 'v3ci' 'v4ci'}
  day = '28-Jul-2009';
  
 case {'v112' 'v212' 'v312' 'v412'}
  day = '28-Jul-2009';
  
end

% set up analysis folename:
switch timeseries
 % original rois used for the main analysis of crfs, and noise
 case 'v1'
  ts = 'v1_0_7roi';
 case 'v2'
  ts = 'v2_0_5roi';
 case 'v3'
  ts = 'v3_0_5roi';
 case 'v4'
  ts = 'v4_0_5roi';
  
  % combined roi (for STD analysis) files
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
 
  % correct and incorrect analysis files
 case 'v1ci'
  ts = 'v1_correctIncorrectAnalysis_STD_v1';
 case 'v2ci'
  ts = 'v2_correctIncorrectAnalysis_STD_v2';
 case 'v3ci'
  ts = 'v3_correctIncorrectAnalysis_STD_v3';
 case 'v4ci'
  ts = 'v4_correctIncorrectAnalysis_STD_v4';
  
  % 1st and 2nd interval analysis files
 case 'v112'
  ts = 'v1_firstSecondIntervalAnalysis_STD_v1';
 case 'v212'
  ts = 'v2_firstSecondIntervalAnalysis_STD_v2';
 case 'v312'
  ts = 'v3_firstSecondIntervalAnalysis_STD_v3';
 case 'v412'
  ts = 'v4_firstSecondIntervalAnalysis_STD_v4';
 
 otherwise
  keyboard
end

% add the day the file was created:
timeseries = sprintf('%s_%s',day,ts);


%%%%%%%%%%%%%%%
% getfilename %
%%%%%%%%%%%%%%%
function filename = getfilename(observer,adaptation,timeseries)

% get the session:
session = getsession(observer,adaptation);

% get the right filename:
ts = switchTSfile(timeseries);


% construct the filename
filename = sprintf('%s/%s',session,ts);


