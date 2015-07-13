function makenoisestatisticsfigs2_rAVRG(varargin)
%
% function makenoisestatisticsfigs2_rAVRG({plotType},{observer},{visualarea},{adaptation})
%
% this function calls paperadptfigs.m to plot figures
% for each observer, visual areas and adaptation condition.
% and generate some basic noise value statistics.
%
% e.g.,
% makenoisestatisticsfigs2_rAVRG('observer',{'jg' 'fp' 'fm'},'visualarea',{'v1'},'adaptation',[0])
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
plotLog   =[];      % if '0' makes linear axis crf plots
filename  =[];
doIndividualPlots=[];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set default variables:
getArgs(varargin,{...
 'observer',{'average'}, ...
 'fitType',{'tvc2crf'}, ...  'crf2tvc'
 'visualarea',{'v1' 'v2' 'v3' 'v4'}, ...
 'adaptation',[0 28 100], ...
 'testType', {'att'}, ...
 'dprime',1, ...
 'savePlots',0, ...
 'dispFit',0, ...
 'recomputeNoise',0, ...
 'saveFig',1, ...
 'filename',[date,'allNoise.mat'], ...
 'doIndividualPlots',1, ...
 'plotLog',1});


if recomputeNoise
 % set default data folder:
 [u computer] = system('hostname');
 
 % if i am running on the laptop:
 if strcmp(strtrim(computer),'ticonzero.local')
 % the i have the external HD homosacer_data mounted
 dataPath = '/Volumes/homosacer_data/riken/crfaa/fmridata/crfaa/data_used_files/';
 
 else % assuming this is run on homoacer at riken
 dataPath = '/Users/frakkopesto/data/riken/crfaa/fmridata/crfaa/data_used_files/';
 end
 
 noisetype = 'additive';
 whichCRF = {'crf_distributed_target' 'crf_attended'};
 data = [];
 
 for iFit = 1:length(fitType) % do fit of crf-2-tv OR of tvc-2-crf
  % initialize the noise struct to nan:
  noise{iFit}.k = nan.*ones(length(adaptation),length(visualarea),length(whichCRF));
  noise{iFit}.responseOffset = nan.*ones(length(adaptation),length(visualarea),length(whichCRF));
  
  disp(sprintf('[makenoisestatisticsfigs2_rAVRG] FIT: %s',fitType{iFit}))
  for iArea = 1:length(visualarea)
   disp(sprintf('[makenoisestatisticsfigs2_rAVRG] AREA: %s',visualarea{iArea}))
   
   whichAdpt = [1,2,3];
   for iAdapt = 1:length(whichAdpt)
    disp(sprintf('[makenoisestatisticsfigs2_rAVRG] ADP: %i',adaptation(iAdapt)))
    
    % get file name
    [thisfile vars2load] = getFileName(visualarea{iArea});
    thisfullfile = fullfile(dataPath,thisfile);
    
    % load file
    data{iFit, iArea} = load(thisfullfile,vars2load{1},vars2load{2});
    data{iFit, iArea}.thisfullfile = thisfullfile;
    data{iFit, iArea}.mfilename = mfilename;
    disp(sprintf('[makenoisestatisticsfigs2_rAVRG] Loaded file: %s.',thisfullfile))
    
    for icrf = 1:length(whichCRF) % distributed or attended
     disp(sprintf('[makenoisestatisticsfigs2_rAVRG] CRF: %s',whichCRF{icrf}))
     % set up general structures and info
     %      fig = setUpPlot(whichCRF{icrf},adaptation(iAdapt),savePlots, dispFit);
     fig = 1;
     % set the minOffset and maxk for fitting:
     if icrf == 2
      % if we fitted the model for distributed
      % use the distributed offset as min for attended offset
      % this should solve some problems with negative noise redux indexes
      minOffset = n.responseOffset;
      maxk = n.k*(.6+rand*.1);
     else
      minOffset = 0;
      maxk = inf;
     end
     
     [n exptData] = dprimefit3_rAVRG(data{iFit, iArea}, ...
      'doNoise',fitType{iFit},       ...
      'noise',[],                    ...
      'dprimefittype','nkn',         ...
      'crfFittype', 'Naka',          ...
      'adaptationIndex',iAdapt,      ...
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
     noise{iFit}.k(iAdapt, iArea, icrf) = n.k;
     noise{iFit}.responseOffset(iAdapt, iArea, icrf) = n.responseOffset;
     noise{iFit}.bestFit{iAdapt, iArea, icrf} = n.bestFit;
     
     disp(sprintf('[makenoisestatisticsfigs2_rAVRG] AREA: %s',visualarea{iArea}))
     disp(sprintf('[makenoisestatisticsfigs2_rAVRG] CRF: %s',whichCRF{icrf}))
     disp(sprintf('[makenoisestatisticsfigs2_rAVRG] Noise %s',num2str(n.k)))
     disp(sprintf('[makenoisestatisticsfigs2_rAVRG] OffSet %s',num2str(n.responseOffset)))
     
     if icrf == 2
      offsetchange = noise{iFit}.responseOffset(iAdapt, iArea, 2)/noise{iFit}.responseOffset(iAdapt, iArea, 1);
      k = noise{iFit}.k(iAdapt, iArea, 1)/noise{iFit}.k(iAdapt, iArea, 2);
      disp(sprintf('[makenoisestatisticsfigs2_rAVRG] OffSet increase: Offset with attention is [%s] times larger.',num2str(offsetchange)))
      disp(sprintf('[makenoisestatisticsfigs2_rAVRG] Noise decrease: Noise with attention is [%s] times smaller.',num2str(k)))
     end
     
     % store exptdata:
     results{iFit, iAdapt, iArea, icrf} = exptData;
     clear exptData;
     
    end
   end
   % save results:
   noise{iFit}.mfilename = mfilename;
   filename = '12_Jul_2009_AVERAGE_Noise.mat';
   save(sprintf('%s',filename),'noise', 'results');
  end
 end
 
 
else
 % load noise file:
  % set default data folder:
 [u computer] = system('hostname');
 
 % if i am running on the laptop:
 if strcmp(strtrim(computer),'ticonzero.local')
  % the i have the external HD homosacer_data mounted
  datadir = '/Volumes/homosacer_data/riken/crfaa/fmridata/crfaa/data_used_files/';

 else % otherwise assuming we are on homosacer at riken
  datadir =  '/Volumes/data/riken/crfaa/fmridata/crfaa/';
 end

 filename = '12_Jul_2009_AVERAGE_Noise.mat'; % '04-Apr-2009AllNoise.mat';
 s = load(sprintf('%s',fullfile(datadir,filename)),'noise','results');
 noise = s.noise;
 results = s.results;
 clear s;
end


% plot figures
if doIndividualPlots
 % ATTENTION PLOT: plot atteded and distributed crf and tvc:
  makeIndAttentionPlot(results,noise,fitType,saveFig,   plotLog)
 
 % ADAPTATION PLOT: plot atteded and distributed crf and tvc:
 %  makeIndAdaptationPlot(results,noise,fitType,saveFig,  plotLog)
 
 % scatter plot of adaptation and attention for both offset and k:
 %  noiseIndScatterPlot(noise{1},saveFig)
 noiseScatterPlot(noise{1},saveFig)
 %  plotParameters(results,noise,saveFig)
else
 keyboard
end

% END MAIN CALL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%
% makeIndAdaptationPlot %
%%%%%%%%%%%%%%%%%%%%%%%%%
function makeIndAdaptationPlot(results,noise,fitType,saveFig,plotLog)
% makes individual plots for attention effect, tvc and crf
if ~strcmpi(fitType,'tvc2crf')
 keyboard
end

% set up basic plot info:
plotInfo = setUpPlotInfo;

% plot only one TvC disregarding visual areas:
iFit = 1;
iArea = 1;
attention = {'Dist' 'Att'};

% plot TvCs
for icrf = 1:size(results,4)     % attentional conditions
 figurename = deblank(sprintf('ADPnoiseAVRG_TVCfit_%s',attention{icrf}));
 h = smartfig(figurename,'reuse');
 
 for iAdapt = 1:3 % adaptation conditions
  Contrast = results{iFit, iAdapt, iArea, icrf}.pedestalsTvC;
  
  TvC = squeeze(results{iFit, iAdapt, iArea, icrf}.behavior.tvc.meanTvC(iAdapt,icrf, :));
  TvC_low = squeeze(results{iFit, iAdapt, iArea, icrf}.behavior.tvc.steTvC(iAdapt,icrf, :));
  TvC_high = squeeze(results{iFit, iAdapt, iArea, icrf}.behavior.tvc.steTvC(iAdapt,icrf, :));
  
  TvCsmoothX = results{iFit, iAdapt, iArea, icrf}.behavior.tvcfit.fitParams.fit.x;
  TvCsmoothY = results{iFit, iAdapt, iArea, icrf}.behavior.tvcfit.fitParams.fit.y;
  
  % plot the smooth function only starting at the first pedestal:
  TvCsmoothX = TvCsmoothX(find(TvCsmoothX >= Contrast(1)));
  TvCsmoothY = TvCsmoothY(find(TvCsmoothX >= Contrast(1)));
  TvCsmoothX = TvCsmoothX(find(TvCsmoothX <= Contrast(end)));
  TvCsmoothY = TvCsmoothY(find(TvCsmoothX <= Contrast(end)));
  
  % now plot the data:
  myerrorbar(Contrast,TvC, ...
   'yLow',TvC_low, 'yHigh',TvC_high, ...
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
 
 if saveFig
  set(h,'PaperPosition',[.25 .25 8 10.5]);
  set(h,'PaperOrientation','Portrait');
  
  AxLabel = 'LogLog';
  figurename = sprintf('%s_%s',figurename,AxLabel);
  disp(sprintf('Saving figure %s',figurename));
  figDir = savefig(h,'figName',figurename,'figDir','fig_tvc_noise');
 end
end
drawnow

% plot crf fits
visareas = {'v1' 'v2' 'v3' 'v4'};

for iArea = 1:size(results,3)
 % plot estimated CRFs given the current noise model:
 for icrf = 1:size(results,4)
  figurename = sprintf('ADPnoiseAVRG_CRFfit%s_%s', attention{icrf},visareas{iArea});
  h = smartfig(figurename,'reuse');
  
  % first plot the fitted functions:
  for iAdapt = 1:3
   CRFsmoothX = results{iFit, iAdapt, iArea, icrf}.noise.contrast;
   CRFsmoothY = results{iFit, iAdapt, iArea, icrf}.noise.response;
   
   Contrast = results{iFit, iAdapt, iArea, icrf}.pedestalsCRF;
   
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
   n = noise{iFit}.k(iAdapt, iArea, icrf);
   offset = noise{iFit}.responseOffset(iAdapt, iArea, icrf);
   L{iAdapt} = sprintf('Noise: %1.3f\n Offset: %1.3f',n,offset);
   
  end
  legend({L{:}},'Location','NorthEastOutside');
  
  % now plot the data:
  for iAdapt = 1:3
   crf_data = results{iFit, iAdapt, iArea, icrf}.use_crf;
   crf_ste = results{iFit, iAdapt, iArea, icrf}.use_crf_ste;
   crf_contrast = results{iFit, iAdapt, iArea, icrf}.pedestalsCRF;
      
   myerrorbar(crf_contrast,crf_data, ...
    'yError',crf_ste, ...
    'Symbol',plotInfo.thisSymbol{icrf}(1),...
    'MarkerFaceColor',plotInfo.MarkerFaceColor{iAdapt}{3-icrf}, ...
    'MarkerEdgeColor',plotInfo.MarkerEdgeColor,...
    'MarkerSize',plotInfo.MarkerSize, ...
    'Color',plotInfo.thisColor{iAdapt}{3-icrf});
  end
  axis('square')
  
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
   myaxisCRFlin;
  end
  %   axis('square')
  
  if saveFig
   set(h,'PaperPosition',[.25 .25 8 10.5]);
   set(h,'PaperOrientation','Portrait');
   
   figurename = sprintf('%s_%s',figurename,AxLabel);
   disp(sprintf('Saving figure %s',figurename));
   disp(sprintf('Saving figure %s',figurename));
   figDir = savefig(h,'figName',figurename,'figDir','fig_crf_noise');
  end
 end
end


%%%%%%%%%%%%%%%%%%%%%%%%%
% makeIndAdaptationPlot %
%%%%%%%%%%%%%%%%%%%%%%%%%
function makeIndAttentionPlot(results,noise,fitType,saveFig,plotLog)
% makes individual plots for attention effect, tvc and crf

if ~strcmpi(fitType,'tvc2crf')
 keyboard
end

% set up basic plot info:
plotInfo = setUpPlotInfo;

% plot only one TvC disregarding visual areas:
iFit = 1;
iArea = 1;
adaptation = {'0' '28' '100'};

% plot TvC %
for iAdapt = 1:size(results,2)
 % open up a new figure for each observer:
 figurename = deblank(sprintf('ATTnoiseAVRG_TVCfit_A%s',num2str(adaptation{iAdapt})));
 h = smartfig(figurename,'reuse');
 
 for icrf = 1:size(results,4)
  Contrast = results{iFit, iAdapt, iArea, icrf}.pedestalsTvC;
  if Contrast(1) == 0
   Contrast(1) = results{iFit, iAdapt, iArea, icrf}.pedestalsTvC(2)/2;
  end
  
  TvC = squeeze(results{iFit, iAdapt, iArea, icrf}.behavior.tvc.meanTvC(iAdapt,icrf, :));
  TvC_low = squeeze(results{iFit, iAdapt, iArea, icrf}.behavior.tvc.steTvC(iAdapt,icrf, :));
  TvC_high = squeeze(results{iFit, iAdapt, iArea, icrf}.behavior.tvc.steTvC(iAdapt,icrf, :));
  
  TvCsmoothX1 = results{iFit, iAdapt, iArea, icrf}.behavior.tvcfit.fitParams.fit.x;
  TvCsmoothY1 = results{iFit, iAdapt, iArea, icrf}.behavior.tvcfit.fitParams.fit.y;
  
  % plot the smooth function only starting at the first pedestal:
  TvCsmoothX = TvCsmoothX1(find(TvCsmoothX1 >= Contrast(1)));
  TvCsmoothY = TvCsmoothY1(find(TvCsmoothX1 >= Contrast(1)));
  TvCsmoothX = TvCsmoothX(find(TvCsmoothX1 <= Contrast(end)));
  TvCsmoothY = TvCsmoothY(find(TvCsmoothX1 <= Contrast(end)));
  
  % now plot the data:
  myerrorbar(Contrast,TvC, ...
   'yLow',TvC_low, 'yHigh',TvC_high, ...
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
 
 if saveFig
  set(h,'PaperPosition',[.25 .25 8 10.5]);
  set(h,'PaperOrientation','Portrait');
  
  AxLabel = 'LogLog';
  figurename = sprintf('%s_%s',figurename,AxLabel);
  disp(sprintf('Saving figure %s',figurename));
  figDir = savefig(h,'figName',figurename,'figDir','fig_tvc_noise');
 end
end



% crf fits
visareas = {'v1' 'v2' 'v3' 'v4'};
for iArea = 1:size(results,3)
 % plot estimated CRFs given the current noise model:
 for iAdapt = 1:size(results,2)
  
  figurename = sprintf('ATTnoiseCRFfitAVRG_A%s%s', num2str(adaptation{iAdapt}),visareas{iArea});
  h = smartfig(figurename,'reuse');
  
  % first plot the fitted functions:
  for icrf = 1:size(results,4)
   CRFsmoothX = results{iFit, iAdapt, iArea, icrf}.noise.contrast;
   CRFsmoothY = results{iFit, iAdapt, iArea, icrf}.noise.response;
   
   Contrast = results{iFit, iAdapt, iArea, icrf}.pedestalsCRF;
   
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
   hold on
   
   % now make a legend with the noie info:
   n = noise{iFit}.k(iAdapt, iArea, icrf);
   offset = noise{iFit}.responseOffset(iAdapt, iArea, icrf);
   
   L{icrf} = sprintf('Noise: %1.3f\n Offset: %1.3f',n,offset);
  end
  
  if icrf == 2
   legend({L{1} L{2}},'Location','NorthEastOutside');
  end
  
  % now plot the data:
  for icrf = 1:size(results,4)
   crf_data = results{iFit, iAdapt, iArea, icrf}.use_crf;
   crf_ste = results{iFit, iAdapt, iArea, icrf}.use_crf_ste;
   crf_contrast = results{iFit, iAdapt, iArea, icrf}.pedestalsCRF;

   myerrorbar(crf_contrast,crf_data, ...
    'yError',crf_ste, ...
    'Symbol',plotInfo.thisSymbol{icrf}(1),...
    'MarkerFaceColor',plotInfo.MarkerFaceColor{iAdapt}{3-icrf}, ...
    'MarkerEdgeColor',plotInfo.MarkerEdgeColor,...
    'MarkerSize',plotInfo.MarkerSize, ...
    'Color',plotInfo.thisColor{iAdapt}{3-icrf});
  end
  axis('square')
  
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
   myaxisCRFlin;
  end
  %   axis('square')
  
  if saveFig
   set(h,'PaperPosition',[.25 .25 8 10.5]);
   set(h,'PaperOrientation','Portrait');
   
   figurename = sprintf('%s_%s',figurename,AxLabel);
   disp(sprintf('Saving figure %s',figurename));
   figDir = savefig(h,'figName',figurename,'figDir','fig_crf_noise');
  end
 end
end


%%%%%%%%%%%%%%%%%%%%%
% makeAttentionPlot %
%%%%%%%%%%%%%%%%%%%%%
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
figurename = 'ATT_NOISE_ScatterPlot_k';
h = smartfig(figurename,'reuse');
for adp = 1:size(noise.k,1)  % adaptation conditions
 for v = 1:size(noise.k,2)  % visual areas
  loglog([10^-4,10^-4;.1,.1],[10^-4,10^-4;.1,.1],'k-','Color',plotInfo.XYColor)
  hold on;
  loglog(noise.k(adp,v,1),noise.k(adp,v,2),plotInfo.thisSymbol{adp}, ...
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


% Attention OFFSET plot:
figurename = 'ATT_NOISE_ScatterPlot_offset';
h = smartfig(figurename,'reuse');
for adp = 1:size(noise.responseOffset,1)  % adaptation conditions
 for v = 1:size(noise.responseOffset,2)  % visual areas
  % response off-set
  plot([0,0;1,1],[0,0;1,1],'k-','Color',plotInfo.XYColor)
  hold on;
  plot(noise.responseOffset(adp,v,1),noise.responseOffset(adp,v,2), ...
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


% adaptation plot k:
figurename = 'ADP_NOISE_ScatterPlot_k';
h = smartfig(figurename,'reuse');

for att = 1:2 % attention
 for v = 1:4 % visual area
  loglog([10^-4,10^-4;1,1],[10^-4,10^-4;1,1],'k-','Color',plotInfo.XYColor)
  hold on;
  loglog(noise.k(1,v,att),noise.k(2,v,att),plotInfo.thisSymbol{1}, ...
   'LineWidth',plotInfo.LineWidth, ...
   'Color',plotInfo.thisColor{att}{v}, ...
   'MarkerFaceColor',plotInfo.MarkerFaceColor{att}{v}, ...
   'MarkerEdgeColor',plotInfo.MarkerEdgeColor, ...
   'MarkerSize', plotInfo.MarkerSize);
  
  loglog(noise.k(1,v,att),noise.k(3,v,att),plotInfo.thisSymbol{2}, ...
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


% adaptation plot offset:
figurename = 'ADP_NOISE_ScatterPlot_offset';
h = smartfig(figurename,'reuse');

for att = 1:2 % attention
 for v = 1:4 % visual area
  % response off-set
  plot([0,0;1,1],[0,0;1,1],'k-','Color',plotInfo.XYColor)
  hold on;
  loglog(noise.responseOffset(1,v,att),noise.responseOffset(2,v,att), ...
   plotInfo.thisSymbol{1}, ...
   'LineWidth',plotInfo.LineWidth,...
   'Color',plotInfo.thisColor{att}{v}, ...
   'MarkerFaceColor',plotInfo.MarkerFaceColor{att}{v}, ...
   'MarkerEdgeColor',plotInfo.MarkerEdgeColor, ...
   'MarkerSize', plotInfo.MarkerSize);
  
  loglog(noise.responseOffset(1,v,att),noise.responseOffset(3,v,att), ...
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


%%%%%%%%%%%%%%%%%%%%
% getFigureHandles %
%%%%%%%%%%%%%%%%%%%%
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

% plot R2:
for adp = 1:size(results,1)
 for v = 1:size(results,2)
  for att = 1:size(results,3)
   noiseFit_R2(adp,att,v) = noise.bestFit{adp,v,att}.r2;
   tvcFit_R2(adp,att,v) = results{adp,v,att}.behavior.tvcfit.fitParams.fit.r2;
  end
 end
end

% Noise Fit: compute the mean across visual areas:
m_noiseFit_R2 = squeeze(mean(noiseFit_R2,3));

% make a bar plot:
h = smartfig('Noise_r2','reuse');
mybar(m_noiseFit_R2', ...
 'xLabelText', 'Noise model r2', ...
 'groupLabels',{'Distributed' 'Attended'}, ...
 'withinGroupLabels',{'Adapt-0' 'Adapt-28' 'Adapt-100'}, ...
 'yAxisMin=0','yAxisMax=1')


% TvC Fit: compute the mean across visual areas:
% remember the TvC are copies of each other. this mean operation is silly.
m_tvcFit_R2 = squeeze(mean(tvcFit_R2,3));

% make a bar plot:
h1 = smartfig('TvC_r2','reuse');
mybar(m_tvcFit_R2', ...
 'xLabelText', 'TvC fit r2', ...
 'groupLabels',{'Distributed' 'Attended'}, ...
 'withinGroupLabels',{'Adapt-0' 'Adapt-28' 'Adapt-100'}, ...
 'yAxisMin=0','yAxisMax=1')


% plot TvC params:
for adp = 1:size(results,1)
 for v = 1:size(results,2)
  for att = 1:size(results,3)
   n(adp,att,v) = results{adp,v,att}.behavior.tvcfit.fitParams.fit.fitParams(1);
   c50(adp,att,v) = results{adp,v,att}.behavior.tvcfit.fitParams.fit.fitParams(2);
   q(adp,att,v) = results{adp,v,att}.behavior.tvcfit.fitParams.fit.fitParams(3);
   dr(adp,att,v) = results{adp,v,att}.behavior.tvcfit.fitParams.fit.fitParams(4);
   T(adp,att,v) = results{adp,v,att}.behavior.tvcfit.fitParams.fit.fitParams(5);
   beta(adp,att,v) = results{adp,v,att}.behavior.tvcfit.fitParams.fit.fitParams(6);
  end
 end
end

% excitation exponent
n    = squeeze(mean(n,3));
% make a bar plot:
h1 = smartfig('n','reuse');
mybar(n', ...
 'xLabelText', 'TvC fit - Input Exponent', ...
 'groupLabels',{'Distributed' 'Attended'}, ...
 'withinGroupLabels',{'Adapt-0' 'Adapt-28' 'Adapt-100'}, ...
 'yAxisMin=0','yAxisMax=5.5')

% threshold
c50  = squeeze(mean(c50,3));
% make a bar plot:
h1 = smartfig('c50','reuse');
mybar(c50', ...
 'xLabelText', 'TvC fit - c50', ...
 'groupLabels',{'Distributed' 'Attended'}, ...
 'withinGroupLabels',{'Adapt-0' 'Adapt-28' 'Adapt-100'}, ...
 'yAxisMin=0','yAxisMax=.1')

% normalization-pool exponent
q    = n.*squeeze(mean(q,3));
% make a bar plot:
h1 = smartfig('q','reuse');
mybar(q', ...
 'xLabelText', 'TvC fit - Normalization exponent', ...
 'groupLabels',{'Distributed' 'Attended'}, ...
 'withinGroupLabels',{'Adapt-0' 'Adapt-28' 'Adapt-100'}, ...
 'yAxisMin=0','yAxisMax=5.5')

% delta response
dr   = squeeze(mean(dr,3));
% make a bar plot:
h1 = smartfig('dr','reuse');
mybar(dr', ...
 'xLabelText', 'TvC fit - Delta Response', ...
 'groupLabels',{'Distributed' 'Attended'}, ...
 'withinGroupLabels',{'Adapt-0' 'Adapt-28' 'Adapt-100'}, ...
 'yAxisMin=0','yAxisMax=.25')

% threshold of the high-dip
T    = c50.*squeeze(mean(T,3));
% make a bar plot:
h1 = smartfig('T','reuse');
mybar(T', ...
 'xLabelText', 'TvC fit - High dipper threshold', ...
 'groupLabels',{'Distributed' 'Attended'}, ...
 'withinGroupLabels',{'Adapt-0' 'Adapt-28' 'Adapt-100'}, ...
 'yAxisMin=0','yAxisMax=1')

% slope of the high-dip
beta = squeeze(mean(beta,3));
% make a bar plot:
h1 = smartfig('beta','reuse');
mybar(beta', ...
 'xLabelText', 'TvC fit - High dipper slope', ...
 'groupLabels',{'Distributed' 'Attended'}, ...
 'withinGroupLabels',{'Adapt-0' 'Adapt-28' 'Adapt-100'}, ...
 'yAxisMin=0','yAxisMax=4')


% make average <p> value across visual areas
for adp = 1:size(results,1)
 for v = 1:size(results,2)
  for att = 1:size(results,3)
   p_noise(adp,att,v) = noise.bestFit{adp,v,att}.chi2.p;
  end
 end
end


% plot probability of the model being good.
A_p_noise = squeeze(p_noise(:,2,:));
D_p_noise = squeeze(p_noise(:,1,:));

% make attended bar plot:
h1 = smartfig('A_p_noise','reuse');
mybar(A_p_noise', ...
 'xLabelText', 'Noise fit Attended - p of the model being good by chance.', ...
 'groupLabels',{'v1' 'v2' 'v3' 'v4'}, ...
 'withinGroupLabels',{'Adapt-0' 'Adapt-28' 'Adapt-100'}, ...
 'yAxisMin=0','yAxisMax=1')

% make distributed bar plot:
h1 = smartfig('D_p_noise','reuse');
mybar(D_p_noise', ...
 'xLabelText', 'Noise fit Distributed - p of the model being good by chance.', ...
 'groupLabels',{'v1' 'v2' 'v3' 'v4'}, ...
 'withinGroupLabels',{'Adapt-0' 'Adapt-28' 'Adapt-100'}, ...
 'yAxisMin=0','yAxisMax=1')


% make average <p> value across visual areas
for adp = 1:size(results,1)
 for v = 1:size(results,2)
  for att = 1:size(results,3)
   p_tvc(adp,att,v) = results{adp,v,att}.behavior.tvcfit.fitParams.fit.chi2.p;
  end
 end
end


% plot probability of the model being good.
A_p_tvc = squeeze(p_tvc(:,2,:));
D_p_tvc = squeeze(p_tvc(:,1,:));


% make a bar plot:
h1 = smartfig('A_p_tvc','reuse');
mybar(A_p_tvc', ...
 'xLabelText', 'TvC fit Attended - p of the model being good by chance.', ...
 'groupLabels',{'Distributed' 'Attended'}, ...
 'withinGroupLabels',{'Adapt-0' 'Adapt-28' 'Adapt-100'}, ...
 'yAxisMin=0','yAxisMax=1')

% make a bar plot:
h1 = smartfig('D_p_tvc','reuse');
mybar(D_p_tvc', ...
 'xLabelText', 'TvC fit Distributed - p of the model being good by chance.', ...
 'groupLabels',{'Distributed' 'Attended'}, ...
 'withinGroupLabels',{'Adapt-0' 'Adapt-28' 'Adapt-100'}, ...
 'yAxisMin=0','yAxisMax=1')


%%%%%%%%%%%%%
% setUpPlot %
%%%%%%%%%%%%%
function plotInfo = setUpPlotInfo
% figure formatting info:
plotInfo.MarkerEdgeColor = [0 0 0];
plotInfo.thisSymbol = {'o-' 'o-' 'o-' 'o-'};
plotInfo.plotPeds = [.00875 1.75 3.5 7 14 28 56 84];

% generate some nice colors:
numC =60;
for i = 1:numC
 c{i} = getSmoothColor(i,numC,'hsv');
 %  plot(i,1,'ko','MarkerFaceColor',c{i},'MarkerSize',26);hold on
 %  text(i,1,sprintf('%i',i),'HorizontalAlignment','Center','Color',[0 0 0]);
end

% choose the ones i like:
reds = {c{2} c{5} c{7}};
greens = {c{23} c{17} c{15}};
blues = {c{40} c{37} c{34}};
purples = {c{46} c{49} c{54}};

for i = 1:3 % adaptation
 plotInfo.thisColor{i} = {reds{i} blues{i} purples{i} greens{i}};
 plotInfo.MarkerFaceColor{i} = plotInfo.thisColor{i};
end

% plot basic set up:
plotInfo.XYColor = [0 0 0];
plotInfo.Fsize = 6;
plotInfo.MarkerSize = 6;
plotInfo.LineWidth = 1;
plotInfo.TickLength = [0.025 .01];
plotInfo.PlotBoxAspectRatio = [1 1 1];
plotInfo.YLim = [0 1.3];
plotInfo.Xlim = [.001,1];
plotInfo.plotPeds = [plotInfo.Xlim(1) 0.0175 0.035 0.07 0.14 0.28 0.56 0.84];


%%%%%%%%%%%%%%%
% getFileName %
%%%%%%%%%%%%%%%
function [filename  vars2load] = getFileName(visualArea)

filename = sprintf('ADAPTATION-EFFECT-%s-19-Jan-2009.mat',visualArea);
vars2load = {'d' 'meanC'};


%%%%%%%%%%%%%%%%%%%%
% noiseScatterPlot %
%%%%%%%%%%%%%%%%%%%%
function noiseScatterPlot(noise,saveFig)
% function noiseScatterPlot(noise)
keyboard

% fix the indexes now it is not plotting right.

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
for j = 1:length(whichAdapt) % adaptation conditions
 for i = 1:size(noise.k,3) % visual areas
  subplot(4,2,2*i-1), loglog([10^-4,10^-4;.1,.1],[10^-4,10^-4;.1,.1],'k-','Color',plotInfo.XYColor)
  hold on;
  subplot(4,2,2*i-1), loglog(noise.k(j,i,1),noise.k(j,i,2),plotInfo.thisSymbol{1}, ...
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
  subplot(4,2,2*i), plot(noise.responseOffset(j,i,1),noise.responseOffset(j,i,2), ...
   plotInfo.thisSymbol{1}, ...
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
 for i = 1:4 % visual area
  subplot(4,2,2*i-1), loglog([10^-4,10^-4;1,1],[10^-4,10^-4;1,1],'k-','Color',plotInfo.XYColor)
  hold on;
  subplot(4,2,2*i-1), loglog(noise.k(1,i,j),noise.k(2,i,j),plotInfo.thisSymbol{1}, ...
   'LineWidth',plotInfo.LineWidth, ...
   'Color',plotInfo.thisColor{j}{i}, ...
   'MarkerFaceColor',plotInfo.MarkerFaceColor{j}{i}, ...
   'MarkerEdgeColor',plotInfo.MarkerEdgeColor, ...
   'MarkerSize', plotInfo.MarkerSize);
  
  subplot(4,2,2*i-1), loglog(noise.k(1,i,j),noise.k(3,i,j),plotInfo.thisSymbol{1}, ...
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
  subplot(4,2,2*i), loglog(noise.responseOffset(1,i,j),noise.responseOffset(2,i,j), ...
   plotInfo.thisSymbol{1}, ...
   'LineWidth',plotInfo.LineWidth,...
   'Color',plotInfo.thisColor{j}{i}, ...
   'MarkerFaceColor',plotInfo.MarkerFaceColor{j}{i}, ...
   'MarkerEdgeColor',plotInfo.MarkerEdgeColor, ...
   'MarkerSize', plotInfo.MarkerSize);
  
  subplot(4,2,2*i), loglog(noise.responseOffset(1,i,j),noise.responseOffset(3,i,j), ...
   plotInfo.thisSymbol{1}, ...
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


if saveFig
 set(h,'PaperPosition',[.25 .25 8 10.5]);
 set(h,'PaperOrientation','Portrait');
 
 disp(sprintf('Saving figure %s',figurename));
 figDir = savefig(h,'figName',figurename);
end


