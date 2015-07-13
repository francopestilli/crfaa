function adp_makenoisestatisticsfigs2_r2_compareNoise(varargin)
%
% function adp_makenoisestatisticsfigs2_r2_compareNoise({plotType},{observer},{visualarea},{adaptation})
%
% this function calls adp_dprimefit3_r2_compareNoise.m
% it does three things:
% (1) computes noise and CRF for distributed
% (2) computes a CRF for attended using the distributed noise estimates in (1)
% (3) computes noise and CRF fot attended
% (4) plots graphs for each one of the above points
%
% the whole thing is done for each observer, visual areas and adaptation condition.
% (A) for the Attention paper. only adapt-0 is used.
% (B) for the Adaptation paper... to be decided
%
% call type e.g.,
% adp_makenoisestatisticsfigs2_r2_compareNoise('observer',{'jg' 'fp' 'fm'},'visualarea',{'v1'},'adaptation',[0])
%
% data files are loaded from:
% './crfaa/fmridata/crfaa/data_used_files/';
%
% results are saved here:
%   here: datadir =  './crfaa/fmridata/crfaa/';
%   as:  filename = [date,'CompareNoise.mat'];
%
% figures are saved here:
%   TvC = './crfaa/fmridata/crfaa/fig_tvc_noiseTest';
%   CRF = './crfaa/fmridata/crfaa/fig_crf_noiseTest';
%
% if savePlot the function will save figures in the current dir.
%
% '_r' means riken, it has been modified in the last stages of the paper
% to use a naka-rushton and its derivative to do the fits
%
% franco pestilli 2010.01.10

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
observer  =[];
visualarea=[];
dprime    =[];     % dprime to use for fits
savePlots =[];     % 1 saves the plot in the current dierectory
dispFit   =[];     % if 1 the plot is displayed
recomputeNoise=[]; % if 1 recomputes noise for the passed conditions otherwise
% it loads up the values in the file in the current dir allnoise.mat
saveFig   =[];     % if 1 saves a figure
dprimefittype=[];  % this is used to test the naka-rushton fit OR the polinomial fit
adaptation=[];     % which adaptation condition to run this on
plotLog   =[];     % choose whether to plot crfs on log scale or on linear scale
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set default variables:
getArgs(varargin,{...
 'observer',{'avrg', 'jg', 'fp', 'fm' }, ... 'avrg', 'jg', 'fp', 'fm'
 'visualarea',{'v1' 'v2' 'v3' 'v4'},    ...
 'dprime',1,                            ...
 'savePlots',1,                         ...
 'dispFit',0,                           ...
 'recomputeNoise',1,                    ...
 'saveFig',1,                           ...
 'filename',[date,'allNoise.mat'],      ...
 'adaptation', 0,                     ...
 'doIndividualPlots',1,                 ... % THIS IS NOT USED ANYMORE LEAVE AT 0
 'dprimefittype','nkn2',                ...
 'plotLog',1});

% find adaptation index:
if adaptation == 0
 adp = 1;
elseif adaptation == 28
 adp = 2;
elseif adaptation == 100
 adp = 3;
else
 keyboard;
end


fitType = {'tvc2crf'}; % fit type used to be either tvc2crf or crf2tvc, for the papers we only use crf2tvc

if recomputeNoise
 % set some parameters:
 % date of the file to be used:
 defautDataFolder = adp_defaultDataFolder;
 
 doBootstrap = 1;
 numBootstraps = 2;
 noisetype = 'additive';
 whichCRF = {'crf_distributed_nontarget' 'crf_attended'};
 data = [];
 noiseTestType = {'D2D' 'D2A' 'A2A'}; % (1) estimate noise for D and plot the results for D,
 % (2) use the noise estimated on D to plot A
 % (3) estimate noise for A and plot the results for A
 
 for iObserver = 1:length(observer) % do it for each observer independently
  disp(sprintf('[adp_makenoisestatisticsfigs2_r2_compareNoise] OBS: %s',observer{iObserver}))
  if ~(strcmpi(observer{iObserver},'fm') && adaptation == 100)
   disp(sprintf('[adp_runNestedHypothesisTest_r] WARNING! Subject fm does not have the 100%% adaptation, exiting.\n[adp_runNestedHypothesisTest_r] To incur in no problems please run fit on observer fm as last in a series of observers.'));
   
   for iArea = 1:length(visualarea)
    disp(sprintf('[adp_makenoisestatisticsfigs2_r2_compareNoise] AREA: %s',visualarea{iArea}))
    for iFit = 1:length(noiseTestType) % compute noise from one condition (D or A) and test it to the other
     % initialize the noise struct to nan:
     % noise{iFit}.k = nan.*ones(length(observer),length(visualarea));
     % noise{iFit}.responseOffset = nan.*ones(length(observer),length(visualarea));
     disp(sprintf('[adp_makenoisestatisticsfigs2_r2_compareNoise] FIT: %s',noiseTestType{iFit}))
     
     % load the proper file.
     [d meanC events thisfullfile doMean] = loadData(observer{iObserver},visualarea{iArea},defautDataFolder,adaptation);
     
     
     data{iFit, iObserver, iArea}{1}.d1     = d;
     data{iFit, iObserver, iArea}{1}.events = events;
     data{iFit, iObserver, iArea}{1}.meanC  = meanC;
     clear d meanC events;
     
     data{iFit, iObserver, iArea}{1}.thisfullfile = thisfullfile;
     data{iFit, iObserver, iArea}{1}.mfilename = mfilename;
     disp(sprintf('[adp_makenoisestatisticsfigs2_r2_compareNoise] Loaded %s',thisfullfile))
     
     % set the minOffset and maxk for fitting:
     minOffset = 0;
     maxk = inf;
     
     switch noiseTestType{iFit}
      case {'D2D'} % estimate noise from distributed and plot the crf
       disp(sprintf('[adp_makenoisestatisticsfigs2_r2_compareNoise] FITTING <Distributed Noise>'))
       
       % set up general structures and info
       fig = setUpPlot(observer{iObserver},whichCRF{1},adaptation,savePlots, dispFit);
       
       % launch the analysis code
       [n exptData] = adp_dprimefit3_r2_compareNoise(data{iFit, iObserver, iArea}, ...
        'doNoise',      fitType{1},    ...
        'noise',        [],            ... % noise is not passed so it will be recomputed
        'dprimefittype',dprimefittype, ...
        'doBootstrap',  doBootstrap,   ...
        'numBootstraps',numBootstraps, ...
        'fudgeIndex',1,           ...
        'whichCRF',  whichCRF{1},      ... % this is set to 1 here to use distributed
        'noisetype', noisetype,        ...
        'dataTitle', [],               ...
        'VisualArea',visualarea{iArea},...
        'figureInfo',fig,              ...
        'dispFit',   dispFit,          ...
        'dprime',    dprime,           ...
        'conditionIndex',[1 1],        ... % NB the second index here was following adaptation, it is now hard coded to 1 (it might screw up the colors)
        'minOffset', minOffset,         ...
        'maxk',      maxk,              ...
        'doMean',    doMean,                           ...
        'adaptation',adp);
       
       
       % save noise info
       noise{iFit,iObserver, iArea}.k              = n.k;
       noise{iFit,iObserver, iArea}.responseOffset = n.responseOffset;
       noise{iFit,iObserver, iArea}.noisetype      = n.noisetype;
       noise{iFit,iObserver, iArea}.bestFit        = n.bestFit;
       
       % store exptdata:
       results{iObserver, iArea, 1} = exptData;
       clear exptData;
       
      case  {'D2A'} % noise was from distributed and it is now used to plot the crf
       disp(sprintf('[adp_makenoisestatisticsfigs2_r2_compareNoise] FITTING <Attended crf> using <Distributed Noise>'))
       % set up general structures and info
       fig = setUpPlot(observer{iObserver},whichCRF{2},adaptation,savePlots, dispFit);
       
       [n exptData] = adp_dprimefit3_r2_compareNoise(data{iFit, iObserver, iArea}, ...
        'doNoise',        fitType{1},                 ...
        'noise',          noise{1, iObserver, iArea}, ... % this line is the one that passes the noise (n) in...
        'dprimefittype',  dprimefittype,              ...
        'doBootstrap',    doBootstrap,                ...
        'numBootstraps',  numBootstraps,              ...
        'fudgeIndex',1,                          ...
        'whichCRF',  whichCRF{2},                     ...  % this is set to 2 here to use attended
        'noisetype', noisetype,                       ...
        'dataTitle', [],                              ...
        'VisualArea',visualarea{iArea},               ...
        'figureInfo',fig,                             ...
        'dispFit',   dispFit,                         ...
        'dprime',    dprime,                          ...
        'conditionIndex',[2 1],                       ... % NB the second index here was following adaptation, it is now hard coded to 1 (it might screw up the colors)
        'minOffset',minOffset,                        ...
        'maxk',     maxk,                             ...
        'doMean',   doMean,                           ...
        'adaptation',adp);
       
       % save noise info
       noise{iFit, iObserver, iArea}.k              = n.k;
       noise{iFit, iObserver, iArea}.responseOffset = n.responseOffset;
       noise{iFit,iObserver, iArea}.noisetype       = n.noisetype;
       noise{iFit, iObserver, iArea}.bestFit        = n.bestFit;
       
       
       % store exptdata:
       results{iObserver, iArea, 2} = exptData;
       clear exptData;
       
       
      case  {'A2A'} % estimate noise from attended and plot the crf
       disp(sprintf('[adp_makenoisestatisticsfigs2_r2_compareNoise] FITTING <Attended Noise>'))
       
       % set up general structures and info
       fig = setUpPlot(observer{iObserver},whichCRF{1},adaptation,savePlots, dispFit);
       clear n
       [n exptData] = adp_dprimefit3_r2_compareNoise(data{iFit, iObserver, iArea}, ...
        'doNoise',       fitType{1},       ...
        'noise',         [],               ... % here noise is not passed, so it will be recomputed
        'dprimefittype', dprimefittype,    ...
        'crfFittype',    'Naka',           ...
        'doBootstrap',   doBootstrap,      ...
        'numBootstraps', numBootstraps,    ...
        'fudgeIndex',    1,                ...
        'whichCRF',      whichCRF{2},      ... % this is set to 2 here to use attended
        'noisetype',     noisetype,        ...
        'dataTitle',     [],               ...
        'VisualArea',    visualarea{iArea},...
        'figureInfo',    fig,              ...
        'dispFit',       dispFit,          ...
        'dprime',        dprime,           ...
        'conditionIndex',[2 1],            ... % NB the second index here was following adaptation, it is now hard coded to 1 (it might screw up the colors)
        'minOffset',     minOffset,        ...
        'maxk',          maxk,             ...
        'doMean',        doMean,                           ...
        'adaptation',adp);
       
       % save noise info
       noise{iFit,iObserver, iArea}.k = n.k;
       noise{iFit,iObserver, iArea}.responseOffset = n.responseOffset;
       noise{iFit,iObserver, iArea}.noisetype = n.noisetype;
       noise{iFit,iObserver, iArea}.bestFit = n.bestFit;
       
       % store exptdata:
       results{iObserver, iArea, 3} = exptData;
       clear exptData;
       
      otherwise
       keyboard
     end
     
     % save results:
     % it saves the file at every iteretion so that if it crashes it will
     % still have saved everything up to here.
     noise{iFit,iObserver, iArea}.mfilename = mfilename;
     
     datadir = adp_defaultDataFolder;
     
     filename = [date,'CompareNoiseFit_',dprimefittype, '_A',num2str(adaptation),'.mat'];
     save(sprintf('%s',fullfile(datadir,filename)),'noise', 'results','data');
    end
   end
  else
   % do nothing
  end
 end
 
else
 datadir = adp_defaultDataFolder;
 filename = [date,'CompareNoiseFit_',dprimefittype,'.mat'];
 filename = '03-Sep-2009CompareNoiseFit_nkn2.mat';
 s = load(sprintf('%s',fullfile(datadir,filename)),'noise','results');
 noise = s.noise;
 results = s.results;
 clear s;
end

if ~isempty(results)
 % plot TvCs and CRFs
 makeIndAttentionPlot(results,noise,observer,fitType,saveFig,plotLog,adaptation);
 
 % make scatter plot
 noiseScatterPlot(noise,saveFig,adaptation)
end
% end main call



%%%%%%%%%%%%%%%%%%%%%%%%
% makeIndAttentionPlot %
%%%%%%%%%%%%%%%%%%%%%%%%
function makeIndAttentionPlot(results,noise,observer,fitType,saveFig,plotLog,adaptation)
% makes individual plots for attention effect, tvc and crf

if ~strcmpi(fitType,'tvc2crf')
 keyboard
end

% set up basic plot info:
plotInfo = setUpPlotInfo;

% plot only one TvC disregarding visual areas:
areas = {'v1' 'v2' 'v3' 'v4' };
tests = {'D2D' 'D2A' 'A2A'};
testColorIndex = [2 1 1];

% plot TvCs
for iObserver = 1:size(results,1)
 for iArea = 1% DO thi only once the rest of the TvC are identical... :size(results,2)
  for itest = 1:3
   % open up a new figure for each observer:
   figurename = deblank(sprintf('noisetest_TvC_testType%s_OBS%s_%s_A%i',tests{itest},observer{iObserver},areas{iArea},adaptation));
   h = smartfig(figurename,'reuse');
   set(h,'color','w')
   Contrast = results{iObserver, iArea, itest}.pedestalsTvC;
   if Contrast(1) == 0
    Contrast(1) = results{iObserver, iArea, itest}.pedestalsTvC(2)/2;
   end
   
   TvCsmoothX1 = results{iObserver, iArea, itest}.behavior.tvcfit.fitParams.fit.x;
   TvCsmoothY1 = results{iObserver, iArea, itest}.behavior.tvcfit.fitParams.fit.y;
   
   % plot the smooth function only starting at the first pedestal:
   TvCsmoothX = TvCsmoothX1(find(TvCsmoothX1 > Contrast(1)));
   TvCsmoothY = TvCsmoothY1(find(TvCsmoothX1 > Contrast(1)));
   TvCsmoothX = TvCsmoothX1(find(TvCsmoothX1 <= Contrast(end)));
   TvCsmoothY = TvCsmoothY1(find(TvCsmoothX1 <= Contrast(end)));
   
   TvC = results{iObserver, iArea, itest}.behavior.tvc.thisTvC;
   % this is to distinguish between idnvidual observers data
   % they have a 8 x 2 ste, and average who has a 8x1 ste
   if size(results{iObserver, iArea, itest}.behavior.tvc.thisTvCste,2) == 2
    keyboard
    % this section should not be used now.
    TvC_low = .3*ones(size(results{iObserver, iArea, itest}.behavior.tvc.thisTvCste(:,1)));
    TvC_high = .5*ones(size(results{iObserver, iArea, itest}.behavior.tvc.thisTvCste(:,2)));
    % now plot the data:
    myerrorbar(Contrast,TvC, ...
     'yLow',TvC.*TvC_low', 'yHigh',TvC.*TvC_high', ...
     'Symbol',plotInfo.thisSymbol{itest}(1),...
     'MarkerFaceColor',plotInfo.MarkerFaceColor{1}{testColorIndex(itest)}, ...
     'MarkerEdgeColor','k',...
     'MarkerSize',plotInfo.MarkerSize, ...
     'Color',plotInfo.thisColor{1}{testColorIndex(itest)});
    
   else % THIS IS TO TAKE CARE OF THE AVERAGE
    TvC_low = results{iObserver, iArea, itest}.behavior.tvc.thisTvCste;
    TvC_high = results{iObserver, iArea, itest}.behavior.tvc.thisTvCste;
    
    % now plot the data:
    myerrorbar(Contrast',TvC, ...
     'yLow',TvC_low, 'yHigh',TvC_high, ...
     'Symbol',plotInfo.thisSymbol{itest}(1),...
     'MarkerFaceColor',plotInfo.MarkerFaceColor{1}{testColorIndex(itest)}, ...
     'MarkerEdgeColor','k',...
     'MarkerSize',plotInfo.MarkerSize, ...
     'Color',plotInfo.thisColor{1}{testColorIndex(itest)});
    
   end
   
   
   % now plot the fit:
   myerrorbar(TvCsmoothX,TvCsmoothY, ...
    'Symbol','-',...
    'yLow',0, 'yHigh',0, ...
    'Color',plotInfo.thisColor{1}{testColorIndex(itest)}, ...
    'LineWidth',plotInfo.LineWidth);
  
  % axis formatting:
  set(gca,...
    'tickDir','out', ...
    'tickLen',[0.05 0.01], ...
    'XScale','log', ...
    'XTick', [.00875 .0175 .035 .07 .14 .28 .56 .84 1], ...
    'XLim',[.00875 1], ...
    'XTickLabel', [0 .0175 .035 .07 .14 .28 .56 .84 1], ...
    'YScale','log', ...
    'YLim',[0.0021875  1], ...
    'YTick',[0.0021875 .004375 .00875 .0175 .035 .07 .14 .28],...
    'YTickLabel', [0.0021875 .004375 .00875 .0175 .035 .07 .14 .28]);
  
  axis('square')
   drawnow
   
   if saveFig
    set(h,'PaperPosition',[.25 .25 8 10.5]);
    set(h,'PaperOrientation','Portrait');
    
    AxLabel = 'LogLog';
    figurename = sprintf('%s_%s',figurename,AxLabel);
    disp(sprintf('Saving figure %s',figurename));
    defaultDataFolder = adp_defaultDataFolder;
    figDir = sprintf('fig_tvc_noiseTest_adp%i',adaptation);
    figDir = savefig(h,'figName',figurename,'defaultDataFolder', defaultDataFolder,'figDir',figDir);
   end
  end
 end
end
drawnow


% crf fits
for itest = 1:3
 for iObserver = 1:size(results,1)
  for iArea = 1:size(results,2)
   % open up a new figure for each observer:
   figurename = deblank(sprintf('noisetest_CRF_testType%s_OBS%s_%s_A%i',tests{itest},observer{iObserver},areas{iArea},adaptation));
   h = smartfig(figurename,'reuse'); set(h,'Name',figurename,'color','w');
   set(h,'color','w')

   % first plot the fitted functions:
   Contrast = results{iObserver, iArea, itest}.pedestalsCRF;
   
   CRFsmoothX1 = results{iObserver, iArea, itest}.noise.contrast;
   CRFsmoothY1 = results{iObserver, iArea, itest}.noise.response;
   
   % plot the smooth function only starting at the first pedestal:
   CRFsmoothX = CRFsmoothX1;
   CRFsmoothY = CRFsmoothY1;
   
   % now plot the fit:
   myerrorbar(CRFsmoothX,CRFsmoothY, ...
    'Symbol','-',...
    'yLow',0, 'yHigh',0, ...
    'Color',plotInfo.thisColor{1}{testColorIndex(itest)}, ...
    'LineWidth',plotInfo.LineWidth);
   
   % now make a legend with the noie info:
   n = noise{itest,iObserver, iArea}.k;
   offset = noise{itest,iObserver, iArea}.responseOffset;
   L{itest} = sprintf('Noise: %1.3f\n Offset: %1.3f',n,offset);
   disp(sprintf('Noise: %1.3f\n Offset: %1.3f',n,offset))
   legend({L{itest}},'Location','NorthEastOutside');
   
  end
 end
 
 % now plot the data
 for iObserver = 1:size(results,1)
  for iArea = 1:size(results,2)
   % open up a new figure for each observer:
   figurename = deblank(sprintf('noisetest_CRF_testType%s_OBS%s_%s_A%i',tests{itest},observer{iObserver},areas{iArea},adaptation));
   h = smartfig(figurename,'reuse'); set(h,'Name',figurename,'color','w');
   set(h,'color','w')

   crf_data = results{iObserver, iArea, itest}.use_crf;
   crf_ste = results{iObserver, iArea, itest}.use_crf_ste;
   crf_contrast = results{iObserver, iArea, itest}.pedestalsCRF;
   
   
   myerrorbar(crf_contrast,crf_data, ...
    'yLow',crf_ste, 'yHigh',crf_ste, ...
    'Symbol',plotInfo.thisSymbol{itest}(1),...
    'MarkerFaceColor',plotInfo.MarkerFaceColor{1}{testColorIndex(itest)}, ...
    'MarkerEdgeColor','k',...
    'MarkerSize',plotInfo.MarkerSize, ...
    'Color',plotInfo.thisColor{1}{testColorIndex(itest)});
   
  % axis formatting:
  if plotLog
    AxLabel = 'log';
    set(gca,...
      'tickDir','out', ...
      'tickLen',[0.05 0.01], ...
      'XScale','log', ...
      'YScale','lin', ...
      'XTick', [0.0021875 .004375 .00875 .0175 .035 .07 .14 .28 .56 .84 1], ...
      'XTickLabel', [0 .004375 .00875 .0175 .035 .07 .14 .28 .56 .84 1], ...
      'YTick', [0 .25 .5 .75 1 1.25 1.5] ,...
      'YTickLabel', [0 .25 .5 .75 1 1.25 1.5] ,...
      'XLim',[0.0021875 1], ...
      'YLim',[0 1.5]);
    axis('square')
    
  else
    AxLabel = 'lin';
    set(gca,...
      'tickDir','out', ...
      'tickLen',[0.05 0.01], ...
      'XScale','lin', ...
      'YScale','lin', ...
      'XTick', [.004375 .00875 .0175 .035 .07 .14 .28 .56 .84 1], ...
      'XTickLabel', [0 .00875 .0175 .035 .07 .14 .28 .56 .84 1], ...
      'YTick', [0 .25 .5 .75 1 1.25 1.5] ,...
      'YTickLabel', [0 .25 .5 .75 1 1.25 1.5] ,...
      'XLim',[0 1], ...
      'YLim',[0 1.5]);
    axis('square')
    
  end
  
   if saveFig
    set(h,'PaperPosition',[.25 .25 8 10.5]);
    set(h,'PaperOrientation','Portrait');
    figurename = sprintf('%s_%s',figurename,AxLabel);
    disp(sprintf('Saving figure %s',[defaultDataFolder,figurename]));
    defaultDataFolder = adp_defaultDataFolder;
    
    figDir = sprintf('fig_crf_noiseTest_adp%i',adaptation);
    figDir = savefig(h,'figName',figurename,'defaultDataFolder', defaultDataFolder,'figDir',figDir);
   end
  end
 end
end


%%%%%%%%%%%%%%%%%%%%
% noiseScatterPlot %
%%%%%%%%%%%%%%%%%%%%
function noiseScatterPlot(noise,saveFig,adaptation)

% set up basic plot info:
plotInfo = setUpPlotInfo;
plotInfo.thisSymbol = {'o' 's'  '^' 'd'};
plotInfo.XTicks = [0 .25 .5 .75 1];
plotInfo.YTicks = [0 .25 .5 .75 1];
plotInfo.title = {'v1' 'v2' 'v3' 'v4'};
plotInfo.MarkerSize = [14 10 10 10 10];
plotInfo.YLim = [.001 1];
plotInfo.Xlim = [.001 1];

% plot sigma
% attention figure
figurename = 'attention_NoiseScatterPlot_sigma';
h = smartfig(figurename,'reuse');
set(h,'Name',figurename,'color','w');

loglog([.045/8,.045/8;.18,.18],[.045/8,.045/8;.18,.18],'k-','Color',plotInfo.XYColor)
hold on;
for iObserver = size(noise,2):-1:1
 for iVis = 1:size(noise,3) % visual areas
  loglog(noise{1,iObserver,iVis}.k,noise{3,iObserver,iVis}.k, ...
   plotInfo.thisSymbol{iObserver}, ...
   'LineWidth',plotInfo.LineWidth, ...
   'Color',plotInfo.thisColor{1}{iVis}, ...
   'MarkerFaceColor',plotInfo.MarkerFaceColor{1}{iVis}, ...
   'MarkerEdgeColor','w', ...
   'MarkerSize', plotInfo.MarkerSize(iObserver));
 end
end

title(sprintf('Noise: Sigma'))
xlabel(sprintf('Distributed\n(std of the internal response)'),'FontSize',plotInfo.Fsize);
ylabel(sprintf('Attended\n(std of the internal response)'),'FontSize',plotInfo.Fsize);

set(gca,...
 'FontName','Helvetica','FontSize',plotInfo.Fsize, ...
 'PlotBoxAspectRatio',plotInfo.PlotBoxAspectRatio, ...
 'XLim', [.045/8 .18],'YLim', [.045/8 .18],... 
 'tickDir','out', ...
 'tickLen',[0.05 0.01], ...
 'LineWidth',plotInfo.LineWidth,'TickLength',plotInfo.TickLength, ...
 'XTick', [ .045/8  .045/4 .045/2  .045 .09 .18], ...
 'XTickLabel', [ .045/8  .045/4 .045/2  .045 .09 .18] ,...
 'YTick', [ .045/8  .045/4 .045/2  .045 .09 .18], ...
 'YTickLabel', [ .045/8  .045/4 .045/2  .045 .09 .18] ,...
 'XColor',plotInfo.XYColor,'YColor',plotInfo.XYColor,'Box','off');

if saveFig
 set(h,'PaperPosition',[.25 .25 8 10.5]);
 set(h,'PaperOrientation','Portrait');
 
 disp(sprintf('Saving figure %s',figurename));
 defaultDataDir = adp_defaultDataFolder;
 figDir = sprintf('fig_scatter_noiseTest_adp%s',num2str(adaptation));
 figDir = savefig(h,'defaultDataFolder',defaultDataDir,'figDir',figDir,'figName',figurename);
end


% plot offset
figurename = 'attention_NoiseScatterPlot_offset';
h = smartfig(figurename,'reuse');
set(h,'Name',figurename, 'color','w');

plot([0,0;1,1],[0,0;1,1],'k-','Color',plotInfo.XYColor)
hold on;
for iObserver = size(noise,2):-1:1
 for iVis = 1:size(noise,3) % visual areas
  plot(noise{1,iObserver,iVis}.responseOffset,noise{3,iObserver,iVis}.responseOffset, ...
   plotInfo.thisSymbol{iObserver}, ...
   'LineWidth',plotInfo.LineWidth, ...
   'Color',plotInfo.thisColor{1}{iVis}, ...
   'MarkerFaceColor',plotInfo.MarkerFaceColor{1}{iVis}, ...
   'MarkerEdgeColor','w', ...
   'MarkerSize', plotInfo.MarkerSize(iObserver));
 end
end

title(sprintf('Noise: offset'))
xlabel(sprintf('Distributed\n(fMRI response, percent change image intensity)'),'FontSize',plotInfo.Fsize);
ylabel(sprintf('Attended\n(fMRI response, percent change image intensity)'),'FontSize',plotInfo.Fsize);
set(gca,...
 'FontName','Helvetica','FontSize',plotInfo.Fsize, ...
 'PlotBoxAspectRatio',plotInfo.PlotBoxAspectRatio, ...
 'XLim', [0 1],'YLim', [0 1],...
 'tickDir','out', ...
 'tickLen',[0.05 0.01], ...
 'LineWidth',plotInfo.LineWidth,'TickLength',plotInfo.TickLength, ...
 'XTick', [0 .5 1], 'XTickLabel', [0 .5 1] ,...
 'YTick', [0 .5 1], 'YTickLabel', [0 .5 1] ,...
 'XColor',plotInfo.XYColor,'YColor',plotInfo.XYColor,'Box','off');

if saveFig
 set(h,'PaperPosition',[.25 .25 8 10.5]);
 set(h,'PaperOrientation','Portrait');
 
 disp(sprintf('Saving figure %s',figurename));
 defaultDataDir = adp_defaultDataFolder;
 figDir = sprintf('fig_scatter_noiseTest_adp%i',adaptation);
 figDir = savefig(h,'defaultDataFolder',defaultDataDir,'figDir',figDir,'figName',figurename);
end


%%%%%%%%%%%%%%%%%
% setUpPlotInfo %
%%%%%%%%%%%%%%%%%
function plotInfo = setUpPlotInfo
% figure formatting info:
plotInfo.MarkerEdgeColor{1} = [0 0 0];
plotInfo.MarkerEdgeColor{2} = [1 1 1];

plotInfo.thisSymbol = {'o-' 'o-' 'o-' 'o-'};
plotInfo.plotPeds = [0 1.75 3.5 7 14 28 56 84];

% generate some nice colors:
numC =60;
for i = 1:numC
 c{i} = getSmoothColor(i,numC,'hsv');
 %  plot(i,1,'ko','MarkerFaceColor',c{i},'MarkerSize',26);hold on
 %  text(i,1,sprintf('%i',i),'HorizontalAlignment','Center','Color',[0 0 0]);
end

%  special colors for scatter plot
plotInfo.customColor{1} = {c{1} c{5} c{8}};% DISTRIBUTED
plotInfo.customColor{2} = {c{42} c{37} c{34}};% ATTENDED
plotInfo.customColor{3} = {c{20} c{26} c{16}};% ADAPT-100


% choose the ones i like:
reds = {c{2} c{5} c{7}};
greens = {c{23} c{17} c{15}};
% greens = {[13 130 73]./255 [40 255 73]./255 [141 255 73]./255};
blues = {c{40} c{37} c{34}};
purples = {c{46} c{49} c{54}};

for i = 1:3 % adaptation
 plotInfo.thisColor{i} = {reds{i} blues{i} purples{i} greens{i}};
 plotInfo.MarkerFaceColor{i} = plotInfo.thisColor{i};
end

% plot basic set up:
plotInfo.XYColor = [0 0 0];
plotInfo.Fsize = 6;
plotInfo.MarkerSize = 10;
plotInfo.LineWidth = 1;
plotInfo.TickLength = [0.025 .01];
plotInfo.PlotBoxAspectRatio = [1 1 1];
plotInfo.YLim = [0 1.3];
plotInfo.Xlim = [.001,1];
plotInfo.plotPeds = [plotInfo.Xlim(1) 0.0175 0.035 0.07 0.14 0.28 0.56 0.84];


%%%%%%%%%%%%%
% setUpPlot %
%%%%%%%%%%%%%
function figs = setUpPlot(observer,whichCRF,adaptation,savePlots, dispFit)
% this function makes a figure name
% and a figure handle depending on the type of test
% observer etc.

figs.name = sprintf('%s_A%i_%s',observer,adaptation,whichCRF);
if dispFit
 figs.handle = smartfig(figs.name,'reuse');
 set(figs.handle,'Name',figs.name);
end
figs.savePlots = savePlots;


%%%%%%%%%%%%
% loadData %
%%%%%%%%%%%%
function [d meanC events filename doMean] = loadData(observer,visualAreas,defautDataFolder,adaptation)
d = []; meanC = []; filename = []; doMean = []; events = [];

% find adaptation index:
if adaptation == 0
 adp = 1;
elseif adaptation == 28
 adp = 2;
elseif adaptation == 100
 adp = 3;
else
 keyboard;
end


% load each file:
if strcmpi(observer,'avrg')
  
 % (0) load data file
 [filename vars2load] = makeFileName(observer,[],visualAreas,defautDataFolder);
 [data] = load(sprintf('%s',filename),vars2load{1}, vars2load{2});
 d = data.d{adp};
 % NB for some reason v1 data of the average has a different threshold at the first
 % pedestal contrast, i am using the v2 mean contrast instead
%  if strcmpi(visualAreas,'v1')
%   [filename vars2load] = makeFileName(observer,[],'v2',defautDataFolder);
%   [data] = load(sprintf('%s',filename),vars2load{1}, vars2load{2});
%   meanC = data.meanC{adp};
%  else
  meanC = data.meanC{adp};
 %end
 clear data;
 events = [];
 doMean = 1;
else % individual observers
 if ~(strcmpi('fm',observer) && adaptation == 100)
  [filename vars2load] = makeFileName(observer,adaptation,visualAreas,defautDataFolder);
  data = load(sprintf('%s',filename),vars2load{1}, vars2load{2});
  
  % now make the average contrast
  meanC = computedMeanC(data.events,.5);
  meanC.quantileAcontrast = 100*meanC.quantileAcontrast';
  meanC.quantileDcontrast = 100*meanC.quantileDcontrast';
  meanC.pedC = 100*meanC.pedC;
  
  d.A_amplitude   = data.d1.amplitude(1:8);
  d.U_amplitude   = data.d1.amplitude(9:16);
  d.Dt_amplitude  = data.d1.amplitude(17:24);
  d.Dnt_amplitude = data.d1.amplitude(25:32);
  
  d.A_ste         = data.d1.amplitudeSTE(1:8);
  d.U_ste         = data.d1.amplitudeSTE(9:16);
  d.Dt_ste        = data.d1.amplitudeSTE(17:24);
  d.Dnt_ste       = data.d1.amplitudeSTE(25:32);
  
  events = data.events;
 end
 doMean = 0;
end

if isempty(d)
 return
end

%%%%%%%%%%%%%%%%
% makeFileName %
%%%%%%%%%%%%%%%%
function [filename vars2load] = makeFileName(observer,adaptation,visualArea,defautDataFolder)

if strcmpi(observer,'avrg')
 filename = sprintf('%sADAPTATION-EFFECT-%s-19-Jan-2009.mat',defautDataFolder,visualArea);
 vars2load = {'d' 'meanC'};
else
 vars2load = 'events';
 if strcmpi(visualArea,'v1')
  filename = sprintf('%s%s_A%s_29-Dec-2008_%s_0_7roi.mat',defautDataFolder,upper(observer),num2str(adaptation),visualArea);
 else
  filename = sprintf('%s%s_A%s_29-Dec-2008_%s_0_5roi.mat',defautDataFolder,upper(observer),num2str(adaptation),visualArea);
 end
 vars2load = {'d1' 'events'};
end


