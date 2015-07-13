function makenoisestatisticsfigs2_r2(varargin)
%
% function makenoisestatisticsfigs2_r2({plotType},{observer},{visualarea},{adaptation})
%
% this function calls paperadptfigs.m to plot figures
% for each observer, visual areas and adaptation condition.
% and generate some basic noise value statistics.
%
% e.g.,
% makenoisestatisticsfigs2_r2('observer',{'jg' 'fp' 'fm'},'visualarea',{'v1'},'adaptation',[0])
%
% if savePlot the function will save figures in the current dir.
%
% '_r' means riken, it has been modified in the last stages of the paper
% to use a naka-rushton and its derivative to do the fits
%
% franco pestilli 2009.07.03

% global figureHandle;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
observer  =[];
visualarea=[];
adaptation=[];
fitType   =[];     % go from tvc to crf OR from crf to tvc
testType  =[];     % test across attention OR across adaptation
dprime    =[];     % dprime to use for fits
savePlots =[];     % 1 saves the plot in the current dierectory
dispFit   =[];     % if 1 the plot is displayed
recomputeNoise=[]; % if 1 recomputes noise for the passed conditions otherwise
% it loads up the values in the file in the current dir allnoise.mat
saveFig   =[];     % if 1 saves a figure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set default variables:
getArgs(varargin,{...
 'observer',{'jg' 'fp' 'fm'}, ... 
 'fitType',{'tvc2crf'}, ...  'crf2tvc'
 'visualarea',{'v1' 'v2' 'v3' 'v4'}, ...
 'adaptation',[0 28 100], ...
 'testType', {'att'}, ...
 'dprime',1, ...
 'savePlots',1, ...
 'dispFit',1, ...
 'recomputeNoise',1, ...
 'saveFig',1, ...
 'filename',[date,'allNoise.mat'], ...
 'doIndividualPlots',1, ...
 'plotLog',1});


if recomputeNoise
 % date of the file to be used:
 day = '29-Dec-2008';
 % data used for plots before jun/jul 2009: dataPath = '/Users/frakkopesto/data/nyu/crfaa/data_used_files/';
 dataPath = '/Users/frakkopesto/data/riken/crfaa/fmridata/crfaa/data_used_files/';
 doBootstrap = 1;
 numBootstraps = 2;
 noisetype = 'additive';
 whichCRF = {'crf_distributed_target' 'crf_attended'};
 data = [];
 
 for iFit = 1:length(fitType) % do fit of crf-2-tv OR of tvc-2-crf
  % initialize the noise struct to nan:
  noise{iFit}.k = nan.*ones(length(observer),length(adaptation),length(visualarea),length(whichCRF));
  noise{iFit}.responseOffset = nan.*ones(length(observer),length(adaptation),length(visualarea),length(whichCRF));
  
  disp(sprintf('[makenoisestatisticsfigs2_r2] FIT: %s',fitType{iFit}))
  for iObserver = 1:length(observer) % do it for each observer independently
   disp(sprintf('[makenoisestatisticsfigs2_r2] OBS: %s',observer{iObserver}))
   % set up adaptation condition:
   if strcmpi(observer{iObserver},'fm')
    whichAdpt = [1,2];
   else
    whichAdpt = [1,2,3];
   end
   
   for iAdapt = 1:length(whichAdpt)
    disp(sprintf('[makenoisestatisticsfigs2_r2] ADP: %i',adaptation(iAdapt)))
    for iArea = 1:length(visualarea)
     disp(sprintf('[makenoisestatisticsfigs2_r2] AREA: %s',visualarea{iArea}))
     
     % get file name
     thisfile = getFileName(observer{iObserver},visualarea{iArea},adaptation(whichAdpt(iAdapt)),day);
     thisfullfile = fullfile(dataPath,thisfile);
          
     % load file
     data{iFit, iObserver, iAdapt, iArea}{1} = load(thisfullfile,'d1','events');
     data{iFit, iObserver, iAdapt, iArea}{1}.thisfullfile = thisfullfile;
     data{iFit, iObserver, iAdapt, iArea}{1}.mfilename = mfilename;
     disp(sprintf('[makenoisestatisticsfigs2_r2] Loaded %s',thisfullfile))
     
     for icrf = 1:length(whichCRF) % distributed or attended
      disp(sprintf('[makenoisestatisticsfigs2_r2] CRF: %s',whichCRF{icrf}))
      % set up general structures and info
      fig = setUpPlot(observer{iObserver},whichCRF{icrf},adaptation(iAdapt),savePlots, dispFit);
      
      % set the minOffset and maxk for fitting:
      if icrf == 2
       % if we fitted the model for distributed
       % use the distributed offset as min for attended offset
       % this should solve some problems with negative noise redux indexes
       minOffset = n.responseOffset/2;
       maxk = n.k*(.6+rand*.1);
      else
       minOffset = 0;
       maxk = inf;
      end
      
      [n exptData] = dprimefit3_r2(data{iFit, iObserver, iAdapt, iArea}, ...
       'doNoise',fitType{iFit},       ...
       'noise',[],                    ...
       'dprimefittype','nkn',         ...
       'crfFittype', 'Naka',          ...
       'doBootstrap',doBootstrap,     ...
       'numBootstraps',numBootstraps, ...
       'adaptationIndex',1,           ...
       'whichCRF',whichCRF{icrf},     ...
       'noisetype',noisetype,         ...
       'dataTitle',[],                ...
       'VisualArea',visualarea{iArea},...
       'figureInfo',fig,              ...
       'dispFit',dispFit,             ...
       'dprime',dprime,               ...
       'conditionIndex',[icrf iAdapt],...
       'minOffset',minOffset,         ...
       'maxk',maxk);
      
      % save noise info
      noise{iFit}.k(iObserver, iAdapt, iArea, icrf) = n.k;
      noise{iFit}.responseOffset(iObserver, iAdapt, iArea, icrf) = n.responseOffset;
      noise{iFit}.bestFit{iObserver, iAdapt, iArea, icrf} = n.bestFit;
      
      disp(sprintf('[makenoisestatisticsfigs2_r2] AREA: %s',visualarea{iArea}))
      disp(sprintf('[makenoisestatisticsfigs2_r2] CRF: %s',whichCRF{icrf}))
      disp(sprintf('[makenoisestatisticsfigs2_r2] Noise %s',num2str(n.k)))
      disp(sprintf('[makenoisestatisticsfigs2_r2] OffSet %s',num2str(n.responseOffset)))
      
      if icrf == 2
       offsetchange = noise{iFit}.responseOffset(iObserver, iAdapt, iArea, 2)/noise{iFit}.responseOffset(iObserver, iAdapt, iArea, 1);
       k = noise{iFit}.k(iObserver, iAdapt, iArea, 1)/noise{iFit}.k(iObserver, iAdapt, iArea, 2);
       disp(sprintf('[makenoisestatisticsfigs2_r2] OffSet increase: Offset with attention is [%s] times larger.',num2str(offsetchange)))
       disp(sprintf('[makenoisestatisticsfigs2_r2] Noise decrease: Noise with attention is [%s] times smaller.',num2str(k)))
      end
      
      % store exptdata:
      results{iObserver, iAdapt, iArea, icrf} = exptData;
      clear exptData;
      
      % save results:
      % it saves the file at every iteretion so that if it crashes it will
      % still have saved everything up to here.
      datadir =  '/Users/frakkopesto/data/riken/crfaa/fmridata/crfaa/';
      noise{iFit}.mfilename = mfilename;
      filename = [date,'AllNoise.mat'];
      save(sprintf('%s',fullfile(datadir,filename)),'noise', 'results');
     end
    end
   end
  end
 end
 
else
 % data used before jun/jul 2009: datadir =  '/Users/frakkopesto/data/nyu/crfaa/papermodelfit/';
 datadir =  '/Users/frakkopesto/data/riken/crfaa/fmridata/crfaa/';
 filename = [date,'AllNoise.mat'];
 s = load(sprintf('%s',fullfile(datadir,filename)),'noise','results');
 noise = s.noise;
 results = s.results;
 clear s;
end



% plot figures
if doIndividualPlots
  makeIndAttentionPlot(results,noise,observer,fitType,saveFig,plotLog);
%   makeIndAdaptationPlot(results,noise,observer,fitType,saveFig,plotLog);
 %  plotParameters(results,noise,saveFig);
 
 % make scatter plot
 %  noiseScatterPlot(noise{1},saveFig)
 %  noiseIndScatterPlot(noise{1},saveFig)
 noiseSingleScatterPlot(noise{1},saveFig)
 
else
 %  % plot atteded and unattended:
%   makeAttentionPlot(results,noise,observer,fitType,visualarea,adaptation,saveFig)
 %  % make scatter plot
 %  noiseScatterPlot(noise{1},saveFig)
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%% END MAIN CALL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%
% makeIndAdaptationPlot %
%%%%%%%%%%%%%%%%%%%%%%%%%
function makeIndAdaptationPlot(results,noise,observer,fitType,saveFig,plotLog)
% makes individual plots for attention effect, tvc and crf

if ~strcmpi(fitType,'tvc2crf')
 keyboard
end

% set up basic plot info:
plotInfo = setUpPlotInfo;

% plot only one TvC disregarding visual areas:
results = squeeze(results);
noise = noise{1};
iArea = 1;
attention = {'Dist' 'Att'};

% plot TvCs
for iObserver = 1:size(results,1) % observers
 for icrf = 1:size(results,4)     % attentional conditions
  if strcmpi(observer{iObserver},'fm')
   numAdapt = 2;
  else
   numAdapt = 3;
  end
  % open up a new figure for each observer:
  figurename = deblank(sprintf('ADPnoiseTVCfit%s_%s',observer{iObserver},attention{icrf}));
  h = smartfig(figurename,'reuse');
  
  for iAdapt = 1:numAdapt % adaptation conditions
   Contrast = results{iObserver, iAdapt, iArea, icrf}.pedestalsTvC;
   if Contrast(1) == 0
    Contrast(1) = results{iObserver, iAdapt, iArea, icrf}.pedestalsTvC(2)/2;
   end
   
   TvC = results{iObserver, iAdapt, iArea, icrf}.behavior.tvc.thisTvC;
   TvC_low = .3*ones(size(results{iObserver, iAdapt, iArea, icrf}.behavior.tvc.thisTvCste(:,1)));
   TvC_high = .5*ones(size(results{iObserver, iAdapt, iArea, icrf}.behavior.tvc.thisTvCste(:,2)));
   
   TvCsmoothX = results{iObserver, iAdapt, iArea, icrf}.behavior.tvcfit.fitParams.fit.x;
   TvCsmoothY = results{iObserver, iAdapt, iArea, icrf}.behavior.tvcfit.fitParams.fit.y;
   
   % plot the smooth function only starting at the first pedestal:
   TvCsmoothX = TvCsmoothX(find(TvCsmoothX >= Contrast(1)));
   TvCsmoothY = TvCsmoothY(find(TvCsmoothX >= Contrast(1)));
   TvCsmoothX = TvCsmoothX(find(TvCsmoothX <= Contrast(end)));
   TvCsmoothY = TvCsmoothY(find(TvCsmoothX <= Contrast(end)));
   
   % now plot the data:
   myerrorbar(Contrast,TvC, ...
    'yLow',TvC.*TvC_low', 'yHigh',TvC.*TvC_high', ...
    'Symbol',plotInfo.thisSymbol{icrf}(1),...
    'MarkerFaceColor',plotInfo.MarkerFaceColor{iAdapt}{3-icrf}, ...
    'MarkerEdgeColor',plotInfo.MarkerEdgeColor,...
    'MarkerSize',plotInfo.MarkerSize, ...
    'Color',plotInfo.thisColor{iAdapt}{3-icrf});
   
   % now plot the fit:
   myerrorbar(TvCsmoothX,TvCsmoothY, ...
    'Symbol','-',...
    'yLow',0, 'yHigh',0, ...
    'Color',plotInfo.thisColor{iAdapt}{3-icrf}, ...
    'LineWidth',plotInfo.LineWidth);
   
  end
  % axis formatting:
  set(gca, ...
   'XScale','log', ...
   'XTick', [.00875 .0175 .035 .07 .14 .28 .56 .84 1], ...
   'XLim',[.00875 1], ...
   'XTickLabel', [0 .0175 .035 .07 .14 .28 .56 .84 1], ...
   'YScale','log', ...
   'YLim',[.004375  1], ...
   'YTick',[.004375 .00875 .0175 .035 .07 .14 .28],...
   'YTickLabel', [.004375 .00875 .0175 .035 .07 .14 .28]);
  
  axis('square')
  myaxisTvC;
  drawnow
  
  if saveFig
   set(h,'PaperPosition',[.25 .25 8 10.5]);
   set(h,'PaperOrientation','Portrait');
   
   AxLabel = 'LogLog';
   figurename = sprintf('%s_%s',figurename,AxLabel);
   disp(sprintf('Saving figure %s',figurename));
   defaultDataFolder = '/Users/frakkopesto/data/riken/crfaa/fmridata/crfaa/';
   figDir = savefig(h,'figName',figurename,'figDir','fig_tvc_AttAdpNoise','defaultDataFolder',defaultDataFolder);
  end
 end
end


visareas = {'v1' 'v2' 'v3' 'v4'};

% plot crf fits
for iObserver = 1:size(results,1)
 for iArea = 1:size(results,3)
  % plot estimated CRFs given the current noise model:
  for icrf = 1:size(results,4)
   if strcmpi(observer{iObserver},'fm')
    numAdapt = 2;
   else
    numAdapt = 3;
   end
   
   figurename = sprintf('ADPnoiseCRFfit%s_%s%s', ...
    observer{iObserver},attention{icrf},visareas{iArea});
   h = smartfig(figurename,'reuse');
   
   % first plot the fitted functions:
   for iAdapt = 1:numAdapt
    CRFsmoothX = results{iObserver, iAdapt, iArea, icrf}.noise.contrast;
    CRFsmoothY = results{iObserver, iAdapt, iArea, icrf}.noise.response;
    
    Contrast = results{iObserver, iAdapt, iArea, icrf}.pedestalsCRF;
    
    % plot the smooth function only starting at the first pedestal:
    CRFsmoothX = CRFsmoothX(find(CRFsmoothX >= Contrast(1)));
    CRFsmoothX = CRFsmoothX(find(CRFsmoothX <= 1));
    
    CRFsmoothY = CRFsmoothY(find(CRFsmoothX >= Contrast(1)));
    CRFsmoothY = CRFsmoothY(find(CRFsmoothX <= 1));
    
    % now plot the fit:
    myerrorbar(CRFsmoothX,CRFsmoothY, ...
     'Symbol','-',...
     'yLow',0, 'yHigh',0, ...
     'Color',plotInfo.thisColor{iAdapt}{3-icrf}, ...
     'LineWidth',plotInfo.LineWidth);
    
    % now make a legend with the noie info:
    n = noise.k(iObserver, iAdapt, iArea, icrf);
    offset = noise.responseOffset(iObserver, iAdapt, iArea, icrf);
    L{iAdapt} = sprintf('Noise: %1.3f\n Offset: %1.3f',n,offset);
    
   end
   legend({L{:}},'Location','NorthEastOutside');
   
   % now plot the data:
   for iAdapt = 1:numAdapt
    crf_data = results{iObserver, iAdapt, iArea, icrf}.use_crf;
    crf_ste = results{iObserver, iAdapt, iArea, icrf}.use_crf_ste;
    crf_contrast = results{iObserver, iAdapt, iArea, icrf}.pedestalsCRF;
    
    myerrorbar(crf_contrast,crf_data, ...
     'yError',crf_ste, ...
     'Symbol',plotInfo.thisSymbol{icrf}(1),...
     'MarkerFaceColor',plotInfo.MarkerFaceColor{iAdapt}{3-icrf}, ...
     'MarkerEdgeColor',plotInfo.MarkerEdgeColor,...
     'MarkerSize',plotInfo.MarkerSize, ...
     'Color',plotInfo.thisColor{iAdapt}{3-icrf});
   end
   
   % axis formatting:
   if plotLog
    AxLabel = 'log';
    set(gca,...
     'XScale','log', ...
     'YScale','lin', ...
     'XTick', [.0021875 .004375 .00875 .0175 .035 .07 .14 .28 .56 .84 1], ...
     'YTick', [.25 .5 .75 1 1.25] ,...
     'XLim',[.0021875 1], ...
     'YLim',[.25 1.45]);
    axis('square')
    myaxisCRF;
   else
    AxLabel = 'lin';
    set(gca,...
     'XScale','lin', ...
     'YScale','lin', ...
     'XTick', [.0021875 .004375 .00875 .0175 .035 .07 .14 .28 .56 .84 1], ...
     'YTick', [.25 .5 .75 1 1.25] ,...
     'XLim',[.0021875 1], ...
     'YLim',[.25 1.45]);
    axis('square')
    myaxisCRFlin;
   end
   
   if saveFig
    set(h,'PaperPosition',[.25 .25 8 10.5]);
    set(h,'PaperOrientation','Portrait');
    
    figurename = sprintf('%s_%s',figurename,AxLabel);
    disp(sprintf('Saving figure %s',figurename));
    disp(sprintf('Saving figure %s',figurename));
    defaultDataFolder = '/Users/frakkopesto/data/riken/crfaa/fmridata/crfaa/';
    figDir = savefig(h,'figName',figurename,'figDir','fig_crf_AttAdpNoise','defaultDataFolder',defaultDataFolder);
   end
  end
 end
end


%%%%%%%%%%%%%%%%%%%%%%%%
% makeIndAttentionPlot %
%%%%%%%%%%%%%%%%%%%%%%%%
function makeIndAttentionPlot(results,noise,observer,fitType,saveFig,plotLog)
% makes individual plots for attention effect, tvc and crf

if ~strcmpi(fitType,'tvc2crf')
 keyboard
end


% set up basic plot info:
plotInfo = setUpPlotInfo;

% plot only one TvC disregarding visual areas:
results = squeeze(results);
noise = noise{1};
iArea = 1;
adaptation = [0 28 100];

% plot TvCs
for iObserver = 1:size(results,1)
 for iAdapt = 1:size(results,2)
  if ~all([strcmpi(observer{iObserver},'fm'), (adaptation(iAdapt) == 100)])
   % open up a new figure for each observer:
   figurename = deblank(sprintf('ATTnoiseTVCfit%sA%s',observer{iObserver},num2str(adaptation(iAdapt))));
   h = smartfig(figurename,'reuse');
   for icrf = 1:size(results,4)
    Contrast = results{iObserver, iAdapt, iArea, icrf}.pedestalsTvC;
    if Contrast(1) == 0
     Contrast(1) = results{iObserver, iAdapt, iArea, icrf}.pedestalsTvC(2)/2;
    end
    
    TvC = results{iObserver, iAdapt, iArea, icrf}.behavior.tvc.thisTvC;
    TvC_low = .3*ones(size(results{iObserver, iAdapt, iArea, icrf}.behavior.tvc.thisTvCste(:,1)));
    TvC_high = .5*ones(size(results{iObserver, iAdapt, iArea, icrf}.behavior.tvc.thisTvCste(:,2)));
    
    TvCsmoothX1 = results{iObserver, iAdapt, iArea, icrf}.behavior.tvcfit.fitParams.fit.x;
    TvCsmoothY1 = results{iObserver, iAdapt, iArea, icrf}.behavior.tvcfit.fitParams.fit.y;
    
    % plot the smooth function only starting at the first pedestal:
    TvCsmoothX = TvCsmoothX1(find(TvCsmoothX1 > Contrast(1)));
    TvCsmoothY = TvCsmoothY1(find(TvCsmoothX1 > Contrast(1)));
    
    % now plot the data:
    myerrorbar(Contrast,TvC, ...
     'yLow',TvC.*TvC_low', 'yHigh',TvC.*TvC_high', ...
     'Symbol',plotInfo.thisSymbol{icrf}(1),...
     'MarkerFaceColor',plotInfo.MarkerFaceColor{iAdapt}{3-icrf}, ...
     'MarkerEdgeColor',plotInfo.MarkerEdgeColor{1},...
     'MarkerSize',plotInfo.MarkerSize, ...
     'Color',plotInfo.thisColor{iAdapt}{3-icrf});
    
    % now plot the fit:
    myerrorbar(TvCsmoothX,TvCsmoothY, ...
     'Symbol','-',...
     'yLow',0, 'yHigh',0, ...
     'Color',plotInfo.thisColor{iAdapt}{3-icrf}, ...
     'LineWidth',plotInfo.LineWidth);
   end
  end
  
  % axis formatting:
  set(gca, ...
   'XScale','log', ...
   'XTick', [.00875 .0175 .035 .07 .14 .28 .56 .84 1], ...
   'XLim',[.00875 1], ...
   'XTickLabel', [0 .0175 .035 .07 .14 .28 .56 .84 1], ...
   'YScale','log', ...
   'YLim',[.004375  1], ...
   'YTick',[.004375 .00875 .0175 .035 .07 .14 .28],...
   'YTickLabel', [.004375 .00875 .0175 .035 .07 .14 .28]);
  
  axis('square')
  myaxisTvC;
  drawnow
  
  if saveFig
   set(h,'PaperPosition',[.25 .25 8 10.5]);
   set(h,'PaperOrientation','Portrait');
   
   AxLabel = 'LogLog';
   figurename = sprintf('%s_%s',figurename,AxLabel);
   disp(sprintf('Saving figure %s',figurename));
   dataPath = '/Users/frakkopesto/data/riken/crfaa/fmridata/crfaa/data_used_files/';
   figDir = savefig(h,'figName',figurename,'defaultDataFolder',dataPath,'figDir','fig_tvc_noise');
  end
 end
end

drawnow
close all;
drawnow
visareas = {'v1' 'v2' 'v3' 'v4'};

% crf fits
for iObserver = 1:size(results,1)
 for iArea = 1:size(results,3)
  % plot estimated CRFs given the current noise model:
  for iAdapt = 1:size(results,2)
   if ~all([strcmpi(observer{iObserver},'fm'), (adaptation(iAdapt) == 100)])
    figurename = sprintf('ATTnoiseCRFfit_OBS%s_A%s_%s', ...
     observer{iObserver},num2str(adaptation(iAdapt)),visareas{iArea});
    h = smartfig(figurename,'reuse'); set(h,'Name',figurename);
    
    % first plot the fitted functions:
    for icrf = 1:size(results,4)
     CRFsmoothX = results{iObserver, iAdapt, iArea, icrf}.noise.contrast;
     CRFsmoothY = results{iObserver, iAdapt, iArea, icrf}.noise.response;
     
     Contrast = results{iObserver, iAdapt, iArea, icrf}.pedestalsCRF;
     
     % plot the smooth function only starting at the first pedestal:
     CRFsmoothX = CRFsmoothX(find(CRFsmoothX >= Contrast(1)));
     CRFsmoothX = CRFsmoothX(find(CRFsmoothX <= 1));
     
     CRFsmoothY = CRFsmoothY(find(CRFsmoothX >= Contrast(1)));
     CRFsmoothY = CRFsmoothY(find(CRFsmoothX <= 1));
     
     % now plot the fit:
     myerrorbar(CRFsmoothX,CRFsmoothY, ...
      'Symbol','-',...
      'yLow',0, 'yHigh',0, ...
      'Color',plotInfo.thisColor{iAdapt}{3-icrf}, ...
      'LineWidth',plotInfo.LineWidth);
     
     % now make a legend with the noie info:
     n = noise.k(iObserver, iAdapt, iArea, icrf);
     offset = noise.responseOffset(iObserver, iAdapt, iArea, icrf);
     L{icrf} = sprintf('Noise: %1.3f\n Offset: %1.3f',n,offset);
     
    end
    if length(L)==2
     legend({L{1} L{2}},'Location','NorthEastOutside');
    end
    % now plot the data:
    for icrf = 1:size(results,4)
     crf_data = results{iObserver, iAdapt, iArea, icrf}.use_crf;
     crf_ste = results{iObserver, iAdapt, iArea, icrf}.use_crf_ste;
     crf_contrast = results{iObserver, iAdapt, iArea, icrf}.pedestalsCRF;
     
     myerrorbar(crf_contrast,crf_data, ...
      'yError',crf_ste, ...
      'Symbol',plotInfo.thisSymbol{icrf}(1),...
      'MarkerFaceColor',plotInfo.MarkerFaceColor{iAdapt}{3-icrf}, ...
      'MarkerEdgeColor',plotInfo.MarkerEdgeColor{1},...
      'MarkerSize',plotInfo.MarkerSize, ...
      'Color',plotInfo.thisColor{iAdapt}{3-icrf});
    end
   end
   
   % axis formatting:
   if plotLog
    AxLabel = 'log';
    set(gca,...
     'XScale','log', ...
     'YScale','lin', ...
     'XTick', [.0021875 .004375 .00875 .0175 .035 .07 .14 .28 .56 .84 1], ...
     'YTick', [.25 .5 .75 1 1.25] ,...
     'XLim',[.0021875 1], ...
     'YLim',[.25 1.45]);
    axis('square')
    myaxisCRF;
   else
    AxLabel = 'lin';
    set(gca,...
     'XScale','lin', ...
     'YScale','lin', ...
     'XTick', [.0021875 .004375 .00875 .0175 .035 .07 .14 .28 .56 .84 1], ...
     'YTick', [.25 .5 .75 1 1.25] ,...
     'XLim',[.0021875 1], ...
     'YLim',[.25 1.45]);
    axis('square')
    myaxisCRFlin;
   end
   
   if saveFig
    set(h,'PaperPosition',[.25 .25 8 10.5]);
    set(h,'PaperOrientation','Portrait');
    
    figurename = sprintf('%s_%s',figurename,AxLabel);
    disp(sprintf('Saving figure %s',figurename));
    disp(sprintf('Saving figure %s',figurename));
    dataPath = '/Users/frakkopesto/data/riken/crfaa/fmridata/crfaa/data_used_files/';
    figDir = savefig(h,'figName',figurename,'defaultDataFolder', dataPath,'figDir','fig_crf_noise');
   end
  end
 end
end


%%%%%%%%%%%%%%%%%%%%
% noiseScatterPlot %
%%%%%%%%%%%%%%%%%%%%
function noiseScatterPlot(noise,saveFig)
% function noiseScatterPlot(noise)

% attention figure
whichAdapt = 1:3;
figurename = 'ATTadditiveNoiseScatterPlot';
h = smartfig(figurename,'reuse');
set(h,'Name',figurename);

% set up basic plot info:
plotInfo = setUpPlotInfo;
plotInfo.thisSymbol = {'s' 'o' '^'};
plotInfo.XTicks = [0 .25 .5 .75 1];
plotInfo.YTicks = [0 .25 .5 .75 1];
plotInfo.title = {'v1' 'v2' 'v3' 'v4'};
plotInfo.MarkerSize = 8;
plotInfo.YLim = [.01 1];
plotInfo.Xlim = [.01 1];

% Attention plot:
for o = 1:3 % observers
 for j = 1:length(whichAdapt) % adaptation conditions
  for i = 1:size(noise.k,3) % visual areas
   subplot(4,2,2*i-1), loglog([10^-4,10^-4;.1,.1],[10^-4,10^-4;.1,.1],'k-','Color',plotInfo.XYColor)
   hold on;
   subplot(4,2,2*i-1), loglog(noise.k(o,j,i,1),noise.k(o,j,i,2),plotInfo.thisSymbol{o}, ...
    'LineWidth',plotInfo.LineWidth, ...
    'Color',plotInfo.thisColor{j}{i}, ...
    'MarkerFaceColor',plotInfo.MarkerFaceColor{j}{i}, ...
    'MarkerEdgeColor',plotInfo.MarkerEdgeColor, ...
    'MarkerSize', plotInfo.MarkerSize);
   
   title(sprintf('%s - k',plotInfo.title{i}))
   if i == 4
    xlabel('Distributed','FontSize',plotInfo.Fsize);
    ylabel(sprintf('Attended'),'FontSize',plotInfo.Fsize);
   end
   
   set(gca,...
    'FontName','Helvetica','FontSize',plotInfo.Fsize, ...
    'PlotBoxAspectRatio',plotInfo.PlotBoxAspectRatio, ...
    'XLim', [min(noise.k(:)) .1],'YLim', [min(noise.k(:)) .1],...
    'LineWidth',plotInfo.LineWidth,'TickLength',plotInfo.TickLength, ...
    'XTick', [.0001 .001 .01 .1 1], 'XTickLabel', [.0001 .001 .01 .1 1] ,...
    'YTick', [.0001 .001 .01 .1 1], 'YTickLabel', [.0001 .001 .01 .1 1] ,...
    'XColor',plotInfo.XYColor,'YColor',plotInfo.XYColor,'Box','off');
   
   % response off-set
   subplot(4,2,2*i), plot([0,0;1,1],[0,0;1,1],'k-','Color',plotInfo.XYColor)
   hold on;
   subplot(4,2,2*i), plot(noise.responseOffset(o,j,i,1),noise.responseOffset(o,j,i,2), ...
    plotInfo.thisSymbol{o}, ...
    'LineWidth',plotInfo.LineWidth, ...
    'Color',plotInfo.thisColor{j}{i},...
    'MarkerFaceColor',plotInfo.MarkerFaceColor{j}{i}, ...
    'MarkerEdgeColor',plotInfo.MarkerEdgeColor, ...
    'MarkerSize', plotInfo.MarkerSize);
   
   title(sprintf('%s - offset',plotInfo.title{i}))
   if i == 4
    xlabel('Distributed','FontSize',plotInfo.Fsize);
    ylabel('Attended','FontSize',plotInfo.Fsize);
   end
   set(gca,...
    'FontName','Helvetica','FontSize',plotInfo.Fsize, ...
    'PlotBoxAspectRatio',plotInfo.PlotBoxAspectRatio, ...
    'XLim', [min(noise.responseOffset(:)) 1],'YLim', [min(noise.responseOffset(:)) 1],...
    'LineWidth',plotInfo.LineWidth,'TickLength',plotInfo.TickLength, ...
    'XTick', plotInfo.XTicks, 'XTickLabel', plotInfo.XTicks ,...
    'YTick', plotInfo.YTicks, 'YTickLabel', plotInfo.YTicks ,...
    'XColor',plotInfo.XYColor,'YColor',plotInfo.XYColor,'Box','off');
  end
 end
end

if saveFig
 set(h,'PaperPosition',[.25 .25 8 10.5]);
 set(h,'PaperOrientation','Portrait');
 
 disp(sprintf('Saving figure %s',figurename));
 figDir = savefig(h,'figName',figurename);
end


% adaptation plot:
figurename = 'ADPadditiveNoiseScatterPlot';
h = smartfig(figurename,'reuse');
set(h,'Name',figurename);

for j = 1:2 % attention
 for o = 1:3 % observers
  for i = 1:4 % visual area
   subplot(4,2,2*i-1), loglog([10^-4,10^-4;1,1],[10^-4,10^-4;1,1],'k-','Color',plotInfo.XYColor)
   hold on;
   subplot(4,2,2*i-1), loglog(noise.k(o,1,i,j),noise.k(o,2,i,j),plotInfo.thisSymbol{o}, ...
    'LineWidth',plotInfo.LineWidth, ...
    'Color',plotInfo.thisColor{j}{i}, ...
    'MarkerFaceColor',plotInfo.MarkerFaceColor{j}{i}, ...
    'MarkerEdgeColor',plotInfo.MarkerEdgeColor, ...
    'MarkerSize', plotInfo.MarkerSize);
   
   subplot(4,2,2*i-1), loglog(noise.k(o,1,i,j),noise.k(o,3,i,j),plotInfo.thisSymbol{o}, ...
    'LineWidth',plotInfo.LineWidth, ...
    'Color',plotInfo.thisColor{j}{i}, ...
    'MarkerFaceColor',plotInfo.MarkerFaceColor{j}{i}, ...
    'MarkerEdgeColor',plotInfo.MarkerEdgeColor, ...
    'MarkerSize', plotInfo.MarkerSize);
   
   title(sprintf('%s - k',plotInfo.title{i}))
   if i == 4
    xlabel('Adapt 0% contrast','FontSize',plotInfo.Fsize);
    ylabel('Adapt 28 or 100% contrast','FontSize',plotInfo.Fsize);
   end
   set(gca,...
    'FontName','Helvetica','FontSize',plotInfo.Fsize, ...
    'PlotBoxAspectRatio',plotInfo.PlotBoxAspectRatio, ...
    'XLim', [min(noise.k(:)) .1],'YLim', [min(noise.k(:)) .1],...
    'LineWidth',plotInfo.LineWidth,'TickLength',plotInfo.TickLength, ...
    'XTick', [.0001 .001 .01 .1 1], 'XTickLabel', [.0001 .001 .01 .1 1] ,...
    'YTick', [.0001 .001 .01 .1 1], 'YTickLabel', [.0001 .001 .01 .1 1] ,...
    'XColor',plotInfo.XYColor,'YColor',plotInfo.XYColor,'Box','off');
   
   % response off-set
   subplot(4,2,2*i), plot([0,0;1,1],[0,0;1,1],'k-','Color',plotInfo.XYColor)
   hold on;
   subplot(4,2,2*i), loglog(noise.responseOffset(o,1,i,j),noise.responseOffset(o,2,i,j), ...
    plotInfo.thisSymbol{o}, ...
    'LineWidth',plotInfo.LineWidth,...
    'Color',plotInfo.thisColor{j}{i}, ...
    'MarkerFaceColor',plotInfo.MarkerFaceColor{j}{i}, ...
    'MarkerEdgeColor',plotInfo.MarkerEdgeColor, ...
    'MarkerSize', plotInfo.MarkerSize);
   
   subplot(4,2,2*i), loglog(noise.responseOffset(o,1,i,j),noise.responseOffset(o,3,i,j), ...
    plotInfo.thisSymbol{o}, ...
    'LineWidth',plotInfo.LineWidth, ...
    'Color',plotInfo.thisColor{j}{i}, ...
    'MarkerFaceColor',plotInfo.MarkerFaceColor{j}{i}, ...
    'MarkerEdgeColor',plotInfo.MarkerEdgeColor, ...
    'MarkerSize', plotInfo.MarkerSize);
   
   title(sprintf('%s - offset',plotInfo.title{i}))
   if i == 4
    xlabel('Adapt 0% contrast','FontSize',plotInfo.Fsize);
    ylabel('Adapt 28 or 100% contrast','FontSize',plotInfo.Fsize);
   end
   set(gca,...
    'FontName','Helvetica','FontSize',plotInfo.Fsize, ...
    'PlotBoxAspectRatio',plotInfo.PlotBoxAspectRatio, ...
    'XLim', [min(noise.responseOffset(:)) 1], ...
    'YLim', [min(noise.responseOffset(:)) 1],...
    'LineWidth',plotInfo.LineWidth,'TickLength',plotInfo.TickLength, ...
    'XTick', plotInfo.XTicks, 'XTickLabel', plotInfo.XTicks ,...
    'YTick', plotInfo.YTicks, 'YTickLabel', plotInfo.YTicks ,...
    'XColor',plotInfo.XYColor,'YColor',plotInfo.XYColor,'Box','off');
  end
 end
end

if saveFig
 set(h,'PaperPosition',[.25 .25 8 10.5]);
 set(h,'PaperOrientation','Portrait');
 
 disp(sprintf('Saving figure %s',figurename));
 figDir = savefig(h,'figName',figurename);
end


%%%%%%%%%%%%%%%%%%
% plotParameters %
%%%%%%%%%%%%%%%%%%
function plotParameters(results,noise,saveFig)
% removing unwanted dimension:
% after this operation:
% dim 1 = adaptation
% dim 2 = visual area
% dim 3 = attention
results = squeeze(results);

% noiseis a '1-d' cel because we can do two tupes of fit (tvc2crf and crf2tvc)
% here i reduce the dimensionality i am only doing one fit for the paper (tvc2crf)
noise = noise{1};

observers = {'FP' 'JG' 'FM'};

% plot R2:
for obs = 1:size(noise.bestFit,1)
 for adp = 1:size(noise.bestFit,2)
  for v = 1:size(noise.bestFit,3)
   for att = 1:size(noise.bestFit,4)
    if ~(obs == 3 && adp == 3)
     noiseFit_R2(adp,att,v) = noise.bestFit{obs,adp,v,att}.r2;
     tvcFit_R2(adp,att,v) = results{obs,adp,v,att}.behavior.tvcfit.fitParams.fit.r2;
    else
     noiseFit_R2(adp,att,v) = nan;
     tvcFit_R2(adp,att,v) = nan;
    end
   end
  end
 end
 
 % Noise Fit: compute the mean across visual areas:
 m_noiseFit_R2 = squeeze(nanmean(noiseFit_R2,3));
 
 % make a bar plot:
 figName = sprintf('Noise_r2_Obs%s',observers{obs});
 h = smartfig(figName,'reuse');
 mybar(m_noiseFit_R2', ...
  'xLabelText', 'Noise model r2', ...
  'groupLabels',{'Distributed' 'Attended'}, ...
  'withinGroupLabels',{'Adapt-0' 'Adapt-28' 'Adapt-100'}, ...
  'yAxisMin=0','yAxisMax=1')
 drawnow
 
 % TvC Fit: compute the mean across visual areas:
 % remember the TvC are copies of each other. this mean operation is silly.
 m_tvcFit_R2 = squeeze(mean(tvcFit_R2,3));
 
 % make a bar plot:
 figName = sprintf('TvC_r2_Obs%s',observers{obs});
 h1 = smartfig(figName,'reuse');
 mybar(m_tvcFit_R2', ...
  'xLabelText', 'TvC fit r2', ...
  'groupLabels',{'Distributed' 'Attended'}, ...
  'withinGroupLabels',{'Adapt-0' 'Adapt-28' 'Adapt-100'}, ...
  'yAxisMin=0','yAxisMax=1')
 drawnow
end


% plot TvC params:
for obs = 1:size(noise.bestFit,1)
 for adp = 1:size(results,2)
  for v = 1:size(results,3)
   for att = 1:size(results,4)
    if ~(obs == 3 && adp == 3)
     n(adp,att,v)    = results{obs,adp,v,att}.behavior.tvcfit.fitParams.fit.fitParams(1);
     c50(adp,att,v)  = results{obs,adp,v,att}.behavior.tvcfit.fitParams.fit.fitParams(2);
     q(adp,att,v)    = results{obs,adp,v,att}.behavior.tvcfit.fitParams.fit.fitParams(3);
     dr(adp,att,v)   = results{obs,adp,v,att}.behavior.tvcfit.fitParams.fit.fitParams(4);
     T(adp,att,v)    = results{obs,adp,v,att}.behavior.tvcfit.fitParams.fit.fitParams(5);
     beta(adp,att,v) = results{obs,adp,v,att}.behavior.tvcfit.fitParams.fit.fitParams(6);
    else
     n(adp,att,v)    = nan;
     c50(adp,att,v)  = nan;
     q(adp,att,v)    = nan;
     dr(adp,att,v)   = nan;
     T(adp,att,v)    = nan;
     beta(adp,att,v) = nan;
    end
   end
  end
 end
 
 % excitation exponent
 n    = squeeze(mean(n,3));
 
 % make a bar plot:
 figName = sprintf('n_Obs%s',observers{obs});
 h = smartfig(figName,'reuse');
 mybar(n', ...
  'xLabelText', 'TvC fit - Input Exponent', ...
  'groupLabels',{'Distributed' 'Attended'}, ...
  'withinGroupLabels',{'Adapt-0' 'Adapt-28' 'Adapt-100'}, ...
  'yAxisMin=0','yAxisMax=5.5')
 
 % threshold
 c50  = squeeze(mean(c50,3));
 
 % make a bar plot:
 figName = sprintf('c50_Obs%s',observers{obs});
 h = smartfig(figName,'reuse');
 mybar(c50', ...
  'xLabelText', 'TvC fit - c50', ...
  'groupLabels',{'Distributed' 'Attended'}, ...
  'withinGroupLabels',{'Adapt-0' 'Adapt-28' 'Adapt-100'}, ...
  'yAxisMin=0','yAxisMax=.1')
 
 % normalization-pool exponent
 q    = n.*squeeze(mean(q,3));
 
 % make a bar plot:
 figName = sprintf('q_Obs%s',observers{obs});
 h = smartfig(figName,'reuse');
 mybar(q', ...
  'xLabelText', 'TvC fit - Normalization exponent', ...
  'groupLabels',{'Distributed' 'Attended'}, ...
  'withinGroupLabels',{'Adapt-0' 'Adapt-28' 'Adapt-100'}, ...
  'yAxisMin=0','yAxisMax=5.5')
 
 % delta response
 dr   = squeeze(mean(dr,3));
 
 % make a bar plot:
 figName = sprintf('dr_Obs%s',observers{obs});
 h = smartfig(figName,'reuse');
 mybar(dr', ...
  'xLabelText', 'TvC fit - Delta Response', ...
  'groupLabels',{'Distributed' 'Attended'}, ...
  'withinGroupLabels',{'Adapt-0' 'Adapt-28' 'Adapt-100'}, ...
  'yAxisMin=0','yAxisMax=.25')
 
 % threshold of the high-dip
 T    = c50.*squeeze(mean(T,3));
 
 % make a bar plot:
 figName = sprintf('T_Obs%s',observers{obs});
 h = smartfig(figName,'reuse');
 mybar(T', ...
  'xLabelText', 'TvC fit - High dipper threshold', ...
  'groupLabels',{'Distributed' 'Attended'}, ...
  'withinGroupLabels',{'Adapt-0' 'Adapt-28' 'Adapt-100'}, ...
  'yAxisMin=0','yAxisMax=1')
 
 % slope of the high-dip
 beta = squeeze(mean(beta,3));
 
 % make a bar plot:
 figName = sprintf('beta_Obs%s',observers{obs});
 h = smartfig(figName,'reuse');
 mybar(beta', ...
  'xLabelText', 'TvC fit - High dipper slope', ...
  'groupLabels',{'Distributed' 'Attended'}, ...
  'withinGroupLabels',{'Adapt-0' 'Adapt-28' 'Adapt-100'}, ...
  'yAxisMin=0','yAxisMax=4')
end


% make average <p> value across visual areas
for obs = 1:size(noise.bestFit,1)
 for adp = 1:size(results,2)
  for v = 1:size(results,3)
   for att = 1:size(results,4)
    if ~(obs == 3 && adp == 3)
     p_noise(adp,att,v) = noise.bestFit{obs,adp,v,att}.chi2.p;
    else
     p_noise(adp,att,v) = nan;
    end
   end
  end
 end
 
 % plot probability of the model being good.
 A_p_noise = squeeze(p_noise(:,2,:));
 D_p_noise = squeeze(p_noise(:,1,:));
 
 % make attended bar plot:
 figName = sprintf('A_p_noise_Obs%s',observers{obs});
 h = smartfig(figName,'reuse');
 mybar(A_p_noise', ...
  'xLabelText', 'Noise fit Attended - p of the model being good by chance.', ...
  'groupLabels',{'v1' 'v2' 'v3' 'v4'}, ...
  'withinGroupLabels',{'Adapt-0' 'Adapt-28' 'Adapt-100'}, ...
  'yAxisMin=0','yAxisMax=1')
 
 % make distributed bar plot:
 figName = sprintf('D_p_noise_Obs%s',observers{obs});
 h = smartfig(figName,'reuse');
 mybar(D_p_noise', ...
  'xLabelText', 'Noise fit Distributed - p of the model being good by chance.', ...
  'groupLabels',{'v1' 'v2' 'v3' 'v4'}, ...
  'withinGroupLabels',{'Adapt-0' 'Adapt-28' 'Adapt-100'}, ...
  'yAxisMin=0','yAxisMax=1')
 
 
 % make average <p> value across visual areas
 for adp = 1:size(results,2)
  for v = 1:size(results,3)
   for att = 1:size(results,4)
    if ~(obs == 3 && adp == 3)
     p_tvc(adp,att,v) = results{obs,adp,v,att}.behavior.tvcfit.fitParams.fit.chi2.p;
    else
     p_tvc(adp,att,v) = nan;
    end
   end
  end
 end
 
 
 % plot probability of the model being good.
 A_p_tvc = squeeze(p_tvc(:,2,:));
 D_p_tvc = squeeze(p_tvc(:,1,:));
 
 
 % make a bar plot:
 figName = sprintf('A_p_TvC_Obs%s',observers{obs});
 h = smartfig(figName,'reuse');
 mybar(A_p_tvc', ...
  'xLabelText', 'TvC fit Attended - p of the model being good by chance.', ...
  'groupLabels',{'Distributed' 'Attended'}, ...
  'withinGroupLabels',{'Adapt-0' 'Adapt-28' 'Adapt-100'}, ...
  'yAxisMin=0','yAxisMax=1')
 
 % make a bar plot:
 figName = sprintf('D_p_TvC_Obs%s',observers{obs});
 h = smartfig(figName,'reuse');
 mybar(D_p_tvc', ...
  'xLabelText', 'TvC fit Distributed - p of the model being good by chance.', ...
  'groupLabels',{'Distributed' 'Attended'}, ...
  'withinGroupLabels',{'Adapt-0' 'Adapt-28' 'Adapt-100'}, ...
  'yAxisMin=0','yAxisMax=1')
 
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


%%%%%%%%%%%%%%%%%%%%%%%
% noiseIndScatterPlot %
%%%%%%%%%%%%%%%%%%%%%%%
function noiseIndScatterPlot(noise,saveFig)

whichAdapt = 1:3;
% set up basic plot info:
plotInfo = setUpPlotInfo;
plotInfo.thisSymbol = {'s' 'o' '^'};
plotInfo.title = {'v1' 'v2' 'v3' 'v4'};
plotInfo.MarkerSize = 12;
plotInfo.YLim = [.001 1];
plotInfo.Xlim = [.001 1];


% Attention NOISE plot:
for o = 1:3
 figurename = sprintf('ATT_NOISE_ScatterPlot_kOBS%s',num2str(o));
 h = smartfig(figurename,'reuse');
 for adp = 1:size(noise.k,2)  % adaptation conditions
  for v = 1:size(noise.k,3)  % visual areas
   loglog([10^-4,10^-4;.1,.1],[10^-4,10^-4;.1,.1],'k-','Color',plotInfo.XYColor)
   hold on;
   loglog(noise.k(o,adp,v,1),noise.k(o,adp,v,2),plotInfo.thisSymbol{adp}, ...
    'LineWidth',plotInfo.LineWidth, ...
    'Color',plotInfo.thisColor{adp}{v}, ...
    'MarkerFaceColor',plotInfo.MarkerFaceColor{adp}{v}, ...
    'MarkerEdgeColor',plotInfo.MarkerEdgeColor, ...
    'MarkerSize', plotInfo.MarkerSize);
   
   if v == 4
    title(sprintf('Average - Noise'))
    xlabel('Distributed','FontSize',plotInfo.Fsize);
    ylabel(sprintf('Attended'),'FontSize',plotInfo.Fsize);
   end
   
   set(gca,...
    'FontName','Helvetica','FontSize',plotInfo.Fsize, ...
    'PlotBoxAspectRatio',plotInfo.PlotBoxAspectRatio, ...
    'XLim', [.005375 .08],'YLim', [.005375 .086],...
    'LineWidth',plotInfo.LineWidth,'TickLength',plotInfo.TickLength, ...
    'XTick', [.005375 .01075 .0215 .043 .086], 'XTickLabel', [.005375 .01075 .0215 .043 .086] ,...
    'YTick', [.005375 .01075 .0215 .043 .086], 'YTickLabel', [.005375 .01075 .0215 .043 .086] ,...
    'XColor',plotInfo.XYColor,'YColor',plotInfo.XYColor,'Box','off');
  end
 end
 myaxisScatterLog;
 
 
 if saveFig
  set(h,'PaperPosition',[.25 .25 8 10.5]);
  set(h,'PaperOrientation','Portrait');
  
  disp(sprintf('Saving figure %s',figurename));
  figDir = savefig(h,'figName',figurename);
 end
end


% Attention OFFSET plot:
for o = 1:3
 figurename = sprintf('ATT_NOISE_ScatterPlot_offsetOBS%s',num2str(o));
 h = smartfig(figurename,'reuse');
 for adp = 1:size(noise.responseOffset,2)  % adaptation conditions
  for v = 1:size(noise.responseOffset,3)  % visual areas
   % response off-set
   plot([0,0;1,1],[0,0;1,1],'k-','Color',plotInfo.XYColor)
   hold on;
   plot(noise.responseOffset(o,adp,v,1),noise.responseOffset(o,adp,v,2), ...
    plotInfo.thisSymbol{adp}, ...
    'LineWidth',plotInfo.LineWidth, ...
    'Color',plotInfo.thisColor{adp}{v}, ...
    'MarkerFaceColor',plotInfo.MarkerFaceColor{adp}{v}, ...
    'MarkerEdgeColor',plotInfo.MarkerEdgeColor, ...
    'MarkerSize', plotInfo.MarkerSize);
   
   if v == 4
    title(sprintf('Average - Offset'))
    xlabel('Distributed','FontSize',plotInfo.Fsize);
    ylabel('Attended','FontSize',plotInfo.Fsize);
   end
   
   set(gca,...
    'FontName','Helvetica','FontSize',plotInfo.Fsize, ...
    'PlotBoxAspectRatio',plotInfo.PlotBoxAspectRatio, ...
    'XLim', [.25 .75],'YLim', [.25 .75],...
    'LineWidth',plotInfo.LineWidth,'TickLength',plotInfo.TickLength, ...
    'XTick', [.25 .5 .75], 'XTickLabel', [.25 .5 .75] ,...
    'YTick', [.25 .5 .75], 'YTickLabel', [.25 .5 .75] ,...
    'XColor',plotInfo.XYColor,'YColor',plotInfo.XYColor,'Box','off');
  end
 end
 xticks.major = [.25 .5 .75];
 xticks.minor = [];
 yticks.major = [.25 .5 .75];
 yticks.minor = [];
 myaxisScatterLin(xticks,yticks);
 
 
 if saveFig
  set(h,'PaperPosition',[.25 .25 8 10.5]);
  set(h,'PaperOrientation','Portrait');
  
  disp(sprintf('Saving figure %s',figurename));
  figDir = savefig(h,'figName',figurename);
 end
end


for o = 1:3
 % adaptation plot k:
 figurename = sprintf('ADP_NOISE_ScatterPlot_kOBS%s',num2str(o));
 h = smartfig(figurename,'reuse');
 
 for att = 1:2 % attention
  for v = 1:4 % visual area
   loglog([10^-4,10^-4;1,1],[10^-4,10^-4;1,1],'k-','Color',plotInfo.XYColor)
   hold on;
   loglog(noise.k(o,1,v,att),noise.k(o,2,v,att),plotInfo.thisSymbol{1}, ...
    'LineWidth',plotInfo.LineWidth, ...
    'Color',plotInfo.thisColor{att}{v}, ...
    'MarkerFaceColor',plotInfo.MarkerFaceColor{att}{v}, ...
    'MarkerEdgeColor',plotInfo.MarkerEdgeColor, ...
    'MarkerSize', plotInfo.MarkerSize);
   
   loglog(noise.k(o,1,v,att),noise.k(o,3,v,att),plotInfo.thisSymbol{2}, ...
    'LineWidth',plotInfo.LineWidth, ...
    'Color',plotInfo.thisColor{att}{v}, ...
    'MarkerFaceColor',plotInfo.MarkerFaceColor{att}{v}, ...
    'MarkerEdgeColor',plotInfo.MarkerEdgeColor, ...
    'MarkerSize', plotInfo.MarkerSize);
   
   if v == 4
    title(sprintf('Adaptation - Noise'))
    xlabel('Adapt 0% contrast','FontSize',plotInfo.Fsize);
    ylabel('Adapt 28 or 100% contrast','FontSize',plotInfo.Fsize);
   end
   set(gca,...
    'FontName','Helvetica','FontSize',plotInfo.Fsize, ...
    'PlotBoxAspectRatio',plotInfo.PlotBoxAspectRatio, ...
    'XLim', [.005375 .08],'YLim', [.005375 .086],...
    'LineWidth',plotInfo.LineWidth,'TickLength',plotInfo.TickLength, ...
    'XTick', [.005375 .01075 .0215 .043 .086], 'XTickLabel', [.005375 .01075 .0215 .043 .086] ,...
    'YTick', [.005375 .01075 .0215 .043 .086], 'YTickLabel', [.005375 .01075 .0215 .043 .086] ,...
    'XColor',plotInfo.XYColor,'YColor',plotInfo.XYColor,'Box','off');
  end
 end
 myaxisScatterLog;
 
 if saveFig
  set(h,'PaperPosition',[.25 .25 8 10.5]);
  set(h,'PaperOrientation','Portrait');
  
  disp(sprintf('Saving figure %s',figurename));
  figDir = savefig(h,'figName',figurename);
 end
end



for o = 1:3
 % adaptation plot offset:
 figurename = sprintf('ADP_NOISE_ScatterPlot_offsetOBS%s',num2str(o));
 h = smartfig(figurename,'reuse');
 
 for att = 1:2 % attention
  for v = 1:4 % visual area
   % response off-set
   plot([0,0;1,1],[0,0;1,1],'k-','Color',plotInfo.XYColor)
   hold on;
   loglog(noise.responseOffset(o,1,v,att),noise.responseOffset(o,2,v,att), ...
    plotInfo.thisSymbol{1}, ...
    'LineWidth',plotInfo.LineWidth,...
    'Color',plotInfo.thisColor{att}{v}, ...
    'MarkerFaceColor',plotInfo.MarkerFaceColor{att}{v}, ...
    'MarkerEdgeColor',plotInfo.MarkerEdgeColor, ...
    'MarkerSize', plotInfo.MarkerSize);
   
   loglog(noise.responseOffset(o,1,v,att),noise.responseOffset(o,3,v,att), ...
    plotInfo.thisSymbol{2}, ...
    'LineWidth',plotInfo.LineWidth, ...
    'Color',plotInfo.thisColor{att}{v}, ...
    'MarkerFaceColor',plotInfo.MarkerFaceColor{att}{v}, ...
    'MarkerEdgeColor',plotInfo.MarkerEdgeColor, ...
    'MarkerSize', plotInfo.MarkerSize);
   
   if v == 4
    title(sprintf('Adaptation - Offset'))
    xlabel('Adapt 0% contrast','FontSize',plotInfo.Fsize);
    ylabel('Adapt 28 or 100% contrast','FontSize',plotInfo.Fsize);
   end
   set(gca,...
    'FontName','Helvetica','FontSize',plotInfo.Fsize, ...
    'PlotBoxAspectRatio',plotInfo.PlotBoxAspectRatio, ...
    'XLim', [.25 .75],'YLim', [.25 .75],...
    'LineWidth',plotInfo.LineWidth,'TickLength',plotInfo.TickLength, ...
    'XTick', [.25 .5 .75], 'XTickLabel', [.25 .5 .75] ,...
    'YTick', [.25 .5 .75], 'YTickLabel', [.25 .5 .75] ,...
    'XColor',plotInfo.XYColor,'YColor',plotInfo.XYColor,'Box','off');
  end
 end
 xticks.major = [.25 .5 .75];
 xticks.minor = [];
 yticks.major = [.25 .5 .75];
 yticks.minor = [];
 myaxisScatterLin(xticks,yticks);
 
 
 if saveFig
  set(h,'PaperPosition',[.25 .25 8 10.5]);
  set(h,'PaperOrientation','Portrait');
  
  disp(sprintf('Saving figure %s',figurename));
  figDir = savefig(h,'figName',figurename);
 end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%
% noiseSingleScatterPlot %
%%%%%%%%%%%%%%%%%%%%%%%%%%
function noiseSingleScatterPlot(noise,saveFig)

whichAdapt = 1:3;
% set up basic plot info:
plotInfo = setUpPlotInfo;
plotInfo.thisSymbol = {'s' 'o' '^' 'd'};
plotInfo.title = {'v1' 'v2' 'v3' 'v4'};
plotInfo.MarkerSize = 12;
% plotInfo.YLim = [.001 1];
% plotInfo.Xlim = [.001 1];

% Attention NOISE plot:
figurename = sprintf('ATT_NOISE_ScatterPlotSingleOBScolor_k');
h = smartfig(figurename,'reuse');
for o = 1:3
 for adp = 1:size(noise.k,2)  % adaptation conditions
  for v = 1:size(noise.k,3)  % visual areas
   loglog([10^-4,10^-4;.1,.1],[10^-4,10^-4;.1,.1],'k-','Color',plotInfo.XYColor)
   hold on;
   loglog(noise.k(o,adp,v,1),noise.k(o,adp,v,2),plotInfo.thisSymbol{o}, ...
    'LineWidth',plotInfo.LineWidth, ...
    'Color',plotInfo.customColor{adp}{o}, ...
    'MarkerFaceColor',plotInfo.customColor{adp}{o}, ...
    'MarkerEdgeColor','w', ...
    'MarkerSize', plotInfo.MarkerSize);
   
   if v == 4
    title(sprintf('Average - Noise'))
    xlabel('Distributed','FontSize',plotInfo.Fsize);
    ylabel(sprintf('Attended'),'FontSize',plotInfo.Fsize);
   end
   
   set(gca,...
    'FontName','Helvetica','FontSize',plotInfo.Fsize, ...
    'PlotBoxAspectRatio',plotInfo.PlotBoxAspectRatio, ...
    'XLim', [0.0026875 .086],'YLim', [0.0026875 .086],...
    'LineWidth',plotInfo.LineWidth,'TickLength',plotInfo.TickLength, ...
    'XTick', [0.0026875 .005375 .01075 .0215 .043 .086], 'XTickLabel', [0.0026875 .005375 .01075 .0215 .043 .086] ,...
    'YTick', [0.0026875 .005375 .01075 .0215 .043 .086], 'YTickLabel', [0.0026875 .005375 .01075 .0215 .043 .086] ,...
    'XColor',plotInfo.XYColor,'YColor',plotInfo.XYColor,'Box','off');
  end
 end
end
myaxisScatterLog;

if saveFig
 set(h,'PaperPosition',[.25 .25 8 10.5]);
 set(h,'PaperOrientation','Portrait');
 
 disp(sprintf('Saving figure %s',figurename));
 figDir = savefig(h,'figName',figurename);
end



% Attention OFFSET plot:
figurename = sprintf('ATT_NOISE_ScatterPlotSingle_OBScolor_offset');
h = smartfig(figurename,'reuse');
for o = 1:3
 for adp = 1:size(noise.responseOffset,2)  % adaptation conditions
  for v = 1:size(noise.responseOffset,3)  % visual areas
   % response off-set
   plot([0,0;1,1],[0,0;1,1],'k-','Color',plotInfo.XYColor)
   hold on;
   plot(noise.responseOffset(o,adp,v,1),noise.responseOffset(o,adp,v,2), ...
    plotInfo.thisSymbol{v}, ...
    'LineWidth',plotInfo.LineWidth, ...
    'Color',plotInfo.customColor{adp}{o}, ...
    'MarkerFaceColor',plotInfo.customColor{adp}{o}, ...
    'MarkerEdgeColor','w', ...
    'MarkerSize', plotInfo.MarkerSize);
   
   if v == 4
    title(sprintf('Average - Offset'))
    xlabel('Distributed','FontSize',plotInfo.Fsize);
    ylabel('Attended','FontSize',plotInfo.Fsize);
   end
   
   set(gca,...
    'FontName','Helvetica','FontSize',plotInfo.Fsize, ...
    'PlotBoxAspectRatio',plotInfo.PlotBoxAspectRatio, ...
    'XLim', [.2 .8],'YLim', [.2 .8],...
    'LineWidth',plotInfo.LineWidth,'TickLength',plotInfo.TickLength, ...
    'XTick', [.2 .5 .8], 'XTickLabel', [.2 .5 .8] ,...
    'YTick', [.2 .5 .8], 'YTickLabel', [.2 .5 .8] ,...
    'XColor',plotInfo.XYColor,'YColor',plotInfo.XYColor,'Box','off');
  end
 end
end
xticks.major = [.2 .5 .8];
xticks.minor = [];
yticks.major = [.2 .5 .8];
yticks.minor = [];
myaxisScatterLin(xticks,yticks);

if saveFig
 set(h,'PaperPosition',[.25 .25 8 10.5]);
 set(h,'PaperOrientation','Portrait');
 
 disp(sprintf('Saving figure %s',figurename));
 dataPath = '/Users/frakkopesto/data/riken/crfaa/fmridata/crfaa/data_used_files/';
 figDir = 'fig_scatter_plots';
 figDir = savefig(h,'figName',figurename,'defaultDataFolder',dataPath,'figDir',figDir);
end


% adaptation plot k:
figurename = sprintf('ADP_NOISE_ScatterPlotSingle_OBScolor_k');
h = smartfig(figurename,'reuse');
for o = 1:3
 for att = 1:2 % attention
  for v = 1:4 % visual area
   loglog([10^-4,10^-4;1,1],[10^-4,10^-4;1,1],'k-','Color',plotInfo.XYColor)
   hold on;
   loglog(noise.k(o,1,v,att),noise.k(o,2,v,att),plotInfo.thisSymbol{v}, ...
    'LineWidth',plotInfo.LineWidth, ...
    'Color',plotInfo.customColor{att}{o}, ...
    'MarkerFaceColor',plotInfo.customColor{att}{o}, ...
    'MarkerEdgeColor','w', ...
    'MarkerSize', plotInfo.MarkerSize);
   
   loglog(noise.k(o,1,v,att),noise.k(o,3,v,att),plotInfo.thisSymbol{v}, ...
    'LineWidth',plotInfo.LineWidth, ...
    'Color',plotInfo.customColor{att}{o}, ...
    'MarkerFaceColor',plotInfo.customColor{att}{o}, ...
    'MarkerEdgeColor','w', ...
    'MarkerSize', plotInfo.MarkerSize);
   
   if v == 4
    title(sprintf('Adaptation - Noise'))
    xlabel('Adapt 0% contrast','FontSize',plotInfo.Fsize);
    ylabel('Adapt 28 or 100% contrast','FontSize',plotInfo.Fsize);
   end
   set(gca,...
    'FontName','Helvetica','FontSize',plotInfo.Fsize, ...
    'PlotBoxAspectRatio',plotInfo.PlotBoxAspectRatio, ...
    'XLim', [0.0026875 .086],'YLim', [0.0026875 .086],...
    'LineWidth',plotInfo.LineWidth,'TickLength',plotInfo.TickLength, ...
    'XTick', [0.0026875 .005375 .01075 .0215 .043 .086], 'XTickLabel', [0.0026875 .005375 .01075 .0215 .043 .086] ,...
    'YTick', [0.0026875 .005375 .01075 .0215 .043 .086], 'YTickLabel', [0.0026875 .005375 .01075 .0215 .043 .086] ,...
    'XColor',plotInfo.XYColor,'YColor',plotInfo.XYColor,'Box','off');
  end
 end
end
myaxisScatterLog;

if saveFig
 set(h,'PaperPosition',[.25 .25 8 10.5]);
 set(h,'PaperOrientation','Portrait');
 
 disp(sprintf('Saving figure %s',figurename));
 dataPath = '/Users/frakkopesto/data/riken/crfaa/fmridata/crfaa/data_used_files/';
 figDir = 'fig_scatter_plots';
 figDir = savefig(h,'figName',figurename,'defaultDataFolder',dataPath,'figDir',figDir);
end


% adaptation plot offset:
figurename = sprintf('ADP_NOISE_ScatterPlotSingle_OBScolor_offset');
h = smartfig(figurename,'reuse');
for o = 1:3
 for att = 1:2 % attention
  for v = 1:4 % visual area
   % response off-set
   plot([0,0;1,1],[0,0;1,1],'k-','Color',plotInfo.XYColor)
   hold on;
   loglog(noise.responseOffset(o,1,v,att),noise.responseOffset(o,2,v,att), ...
    plotInfo.thisSymbol{v}, ...
    'LineWidth',plotInfo.LineWidth,...
    'Color',plotInfo.customColor{att}{o}, ...
    'MarkerFaceColor',plotInfo.customColor{att}{o}, ...
    'MarkerEdgeColor','w', ...
    'MarkerSize', plotInfo.MarkerSize);
   
   loglog(noise.responseOffset(o,1,v,att),noise.responseOffset(o,3,v,att), ...
    plotInfo.thisSymbol{v}, ...
    'LineWidth',plotInfo.LineWidth, ...
    'Color',plotInfo.thisColor{att}{v}, ...
    'Color',plotInfo.customColor{att}{o}, ...
    'MarkerFaceColor',plotInfo.customColor{att}{o}, ...
    'MarkerEdgeColor','w', ...
    'MarkerSize', plotInfo.MarkerSize);
   
   if v == 4
    title(sprintf('Adaptation - Offset'))
    xlabel('Adapt 0% contrast','FontSize',plotInfo.Fsize);
    ylabel('Adapt 28 or 100% contrast','FontSize',plotInfo.Fsize);
   end
   set(gca,...
    'FontName','Helvetica','FontSize',plotInfo.Fsize, ...
    'PlotBoxAspectRatio',plotInfo.PlotBoxAspectRatio, ...
    'XLim', [.2 .8],'YLim', [.2 .8],...
    'LineWidth',plotInfo.LineWidth,'TickLength',plotInfo.TickLength, ...
    'XTick', [.2 .5 .8], 'XTickLabel', [.2 .5 .8] ,...
    'YTick', [.2 .5 .8], 'YTickLabel', [.2 .5 .8] ,...
    'XColor',plotInfo.XYColor,'YColor',plotInfo.XYColor,'Box','off');
  end
 end
end
xticks.major = [.2 .5 .8];
xticks.minor = [];
yticks.major = [.2 .5 .8];
yticks.minor = [];
myaxisScatterLin(xticks,yticks);

if saveFig
 set(h,'PaperPosition',[.25 .25 8 10.5]);
 set(h,'PaperOrientation','Portrait');
 
 disp(sprintf('Saving figure %s',figurename));
 dataPath = '/Users/frakkopesto/data/riken/crfaa/fmridata/crfaa/data_used_files/';
 figDir = 'fig_scatter_plots';
 figDir = savefig(h,'figName',figurename,'defaultDataFolder',dataPath,'figDir',figDir);
end



