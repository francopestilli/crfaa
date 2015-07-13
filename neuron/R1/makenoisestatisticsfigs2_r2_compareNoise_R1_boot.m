function makenoisestatisticsfigs2_r2_compareNoise_R1_boot(varargin)
%
% function makenoisestatisticsfigs2_r2_compareNoise_R1_boot({plotType},{observer},{visualarea},{adaptation})
%
% R1 is the version modified at columbia for the rebuttal.
%
% this function calls dprimefit3_r2_compareNoise_R1.m
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
% makenoisestatisticsfigs2_r2_compareNoise_R1_boot('observer',{'jg' 'fp' 'fm'},'visualarea',{'v1'},'adaptation',[0])
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
% franco pestilli 2009.07.25

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
nBoots    =[];     % number of bootstraps
r2        =[];     % the r2 level of the roi, this is used to load the correct data file and to save the new results
subtractionType  = []; % ORIG or BYTRIAL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set default variables:
getArgs(varargin,{...
 'observer',{ 'jg', 'fm', 'fp'}         ... % 'jg', 'fm', 'fp'no average
 'visualarea',{'v1' 'v2', 'v3' 'v4'},   ...
 'dprime',1,                            ...
 'savePlots',0,                         ...
 'dispFit',0,                           ...
 'recomputeNoise',1,                    ...
 'saveFig',1,                           ...
 'filename',[date,'allNoise_boot_fm.mat'], ...
 'nBoots',100,                            ...
 'adaptation', 0,                       ...
 'doIndividualPlots',0,                 ... % THIS IS NOT USED ANYMORE LEAVE AT 0
 'dprimefittype','nkn2',                  ... % 'nkn2'
 'plotLog',1,                           ...
 'r2','070', ...
 'subtractionType', 'BOOT'});


fitType = {'tvc2crf'}; % fit type used to be either tvc2crf or crf2tvc, for the papers we only use crf2tvc

if recomputeNoise
 % load files from this dir:
 defaultDataFolder = '/data2/crfaa/crfaa/data_used_files_R1/';
 
 doBootstrap = 1;
 numBootstraps = 2;
 noisetype = 'additive';
 whichCRF = {'crf_distributed_target' 'crf_attended'};
 data = [];
 noiseTestType = {'D2D' 'D2A' 'A2A'}; % (1) estimate noise for D and plot the results for D,
 % (2) use the noise estimated on D to plot A
 % (3) estimate noise for A and plot the results for A
 
 
 for iObserver = 1:length(observer) % do it for each observer independently
  disp(sprintf('[makenoisestatisticsfigs2_r2_compareNoise_R1_boot] OBS: %s',observer{iObserver}))
  for iArea = 1:length(visualarea)
   disp(sprintf('[makenoisestatisticsfigs2_r2_compareNoise_R1_boot] AREA: %s',visualarea{iArea}))
   for iFit = 1:length(noiseTestType) % compute noise from one condition (D or A) and test it to the other
    % initialize the noise struct to nan:
    % noise{iFit}.k = nan.*ones(length(observer),length(visualarea));
    % noise{iFit}.responseOffset = nan.*ones(length(observer),length(visualarea));
    disp(sprintf('[makenoisestatisticsfigs2_r2_compareNoise_R1_boot] FIT: %s',noiseTestType{iFit}))
    for iBoot = 1:nBoots
     disp(sprintf('[makenoisestatisticsfigs2_r2_compareNoise_R1_boot] BOOT: %i/%i',iBoot,nBoots))
     
     % load the proper file.
     [d meanC events thisfullfile doMean] = loadData(observer{iObserver},visualarea{iArea}, ...
      defaultDataFolder,adaptation,r2,subtractionType,iBoot);
     
     data{iFit, iObserver, iArea}{1}.d1 = d;
     data{iFit, iObserver, iArea}{1}.events = events;
     data{iFit, iObserver, iArea}{1}.meanC = meanC;
     clear d meanC events;
     
     data{iFit, iObserver, iArea}{1}.thisfullfile = thisfullfile;
     data{iFit, iObserver, iArea}{1}.mfilename = mfilename;
     disp(sprintf('[makenoisestatisticsfigs2_r2_compareNoise_R1_boot] Loaded %s',thisfullfile))
     
     
     % set the minOffset and maxk for fitting:
     minOffset = 0;
     maxk = inf;
     
     switch noiseTestType{iFit}
      case {'D2D'} % estimate noise from distributed and plot the crf
       disp(sprintf('[makenoisestatisticsfigs2_r2_compareNoise_R1_boot] FITTING <Distributed Noise>'))
       
       % set up general structures and info
       fig = setUpPlot(observer{iObserver},whichCRF{1},adaptation,savePlots, dispFit);
       [n exptData] = dprimefit3_r2_compareNoise_R1(data{iFit, iObserver, iArea}, ...
        'doNoise',fitType{1},          ...
        'noise',[],                    ... % noise is not passed so it will be recomputed
        'dprimefittype',dprimefittype, ...
        'doBootstrap',doBootstrap,     ...
        'numBootstraps',numBootstraps, ...
        'adaptationIndex',1,           ...
        'whichCRF',whichCRF{1},        ... % this is set to 1 here to use distributed
        'noisetype',noisetype,         ...
        'dataTitle',[],                ...
        'VisualArea',visualarea{iArea},...
        'figureInfo',fig,              ...
        'dispFit',dispFit,             ...
        'dprime',dprime,               ...
        'conditionIndex',[1 1],        ... % NB the second index here was following adaptation, it is now hard coded to 1 (it might screw up the colors)
        'minOffset',minOffset,         ...
        'maxk',maxk,                   ...
        'doMean', doMean);
       
       
       % save noise info
       
       sensoryNoise(iBoot,iFit,iObserver, iArea)   = n.k;
       responseOffset(iBoot,iFit,iObserver, iArea) = n.responseOffset;
       
       noise{iBoot,iFit,iObserver, iArea}.k              = n.k;
       noise{iBoot,iFit,iObserver, iArea}.responseOffset = n.responseOffset;
       noise{iBoot,iFit,iObserver, iArea}.noisetype      = n.noisetype;
       noise{iBoot,iFit,iObserver, iArea}.bestFit        = n.bestFit;
       
       % store exptdata:
       results{iBoot,iObserver, iArea, 1} = exptData;
       clear exptData;
       
      case  {'D2A'} % noise was from distributed and it is now used to plot the crf
       disp(sprintf('[makenoisestatisticsfigs2_r2_compareNoise_R1_boot] FITTING <Attended crf> using <Distributed Noise>'))
       % set up general structures and info
%        fig = setUpPlot(observer{iObserver},whichCRF{2},adaptation,savePlots, dispFit);
%        
%        [n exptData] = dprimefit3_r2_compareNoise_R1(data{iFit, iObserver, iArea}, ...
%         'doNoise',fitType{1},               ...
%         'noise',noise{1, iObserver, iArea}, ... % this line is the one that passes the noise (n) in...
%         'dprimefittype',dprimefittype,      ...
%         'doBootstrap',doBootstrap,          ...
%         'numBootstraps',numBootstraps,      ...
%         'adaptationIndex',1,                ...
%         'whichCRF',whichCRF{2},             ...  % this is set to 2 here to use attended
%         'noisetype',noisetype,              ...
%         'dataTitle',[],                     ...
%         'VisualArea',visualarea{iArea},     ...
%         'figureInfo',fig,                   ...
%         'dispFit',dispFit,                  ...
%         'dprime',dprime,                    ...
%         'conditionIndex',[2 1],             ... % NB the second index here was following adaptation, it is now hard coded to 1 (it might screw up the colors)
%         'minOffset',minOffset,              ...
%         'maxk',maxk,                        ...
%         'doMean', doMean);
       
       % save noise info
       sensoryNoise(iBoot,iFit,iObserver, iArea)   = nan;%n.k;
       responseOffset(iBoot,iFit,iObserver, iArea) = nan;%n.responseOffset;
       
       noise{iBoot,iFit,iObserver, iArea}.k              = nan;%n.k;
       noise{iBoot,iFit,iObserver, iArea}.responseOffset = nan;%n.responseOffset;
       noise{iBoot,iFit,iObserver, iArea}.noisetype      = nan;%n.noisetype;
       noise{iBoot,iFit,iObserver, iArea}.bestFit        = nan;%n.bestFit;
       
       
       % store exptdata:
       results{iBoot, iObserver, iArea, 2} = nan;%exptData;
       clear exptData;
       
       
      case  {'A2A'} % estimate noise from attended and plot the crf
       disp(sprintf('[makenoisestatisticsfigs2_r2_compareNoise_R1_boot] FITTING <Attended Noise>'))
       
       % set up general structures and info
       fig = setUpPlot(observer{iObserver},whichCRF{1},adaptation,savePlots, dispFit);
       clear n
       [n exptData] = dprimefit3_r2_compareNoise_R1(data{iFit, iObserver, iArea}, ...
        'doNoise',fitType{1},          ...
        'noise',[],                    ... % here noise is not passed, so it will be recomputed
        'dprimefittype',dprimefittype, ...
        'crfFittype', 'Naka',          ...
        'doBootstrap',doBootstrap,     ...
        'numBootstraps',numBootstraps, ...
        'adaptationIndex',1,           ...
        'whichCRF',whichCRF{2},        ... % this is set to 2 here to use attended
        'noisetype',noisetype,         ...
        'dataTitle',[],                ...
        'VisualArea',visualarea{iArea},...
        'figureInfo',fig,              ...
        'dispFit',dispFit,             ...
        'dprime',dprime,               ...
        'conditionIndex',[2 1],        ... % NB the second index here was following adaptation, it is now hard coded to 1 (it might screw up the colors)
        'minOffset',minOffset,         ...
        'maxk',maxk,                   ...
        'doMean', doMean);
       
       % save noise info
       sensoryNoise(iBoot,iFit,iObserver, iArea)   = n.k;
       responseOffset(iBoot,iFit,iObserver, iArea) = n.responseOffset;
       
       noise{iBoot,iFit,iObserver, iArea}.k              = n.k;
       noise{iBoot,iFit,iObserver, iArea}.responseOffset = n.responseOffset;
       noise{iBoot,iFit,iObserver, iArea}.noisetype      = n.noisetype;
       noise{iBoot,iFit,iObserver, iArea}.bestFit        = n.bestFit;
       
       % store exptdata:
       results{iBoot, iObserver, iArea, 3} = exptData;
       clear exptData;
       
      otherwise
       keyboard
     end
     
     noise{iFit,iObserver, iArea}.mfilename = mfilename;
     datadir = '/data2/crfaa/crfaa/data_used_files_R1/';
     

    end
   end
   filename = ['2010-50-Aug','CompareNoiseFit_Dnt_bootNum',num2str(nBoots),'_',dprimefittype,'_',r2,'_',subtractionType,'.mat'];
   disp(sprintf('[makenoisestatisticsfigs2_r2_compareNoise_R1_boot] SAVING interim Results: %s',filename))
   save(sprintf('%s',fullfile(datadir,filename)),'noise', 'results','data','sensoryNoise', 'responseOffset');
  end
 end
 
else
 datadir = '/data2/crfaa/crfaa/data_used_files_R1/';
 
 filename = '2010-10-JulyCompareNoiseFit_Dnt_bootNumd_nkn2_070_BOOT.mat';
 filename = '2010-20-JulyCompareNoiseFit_Dnt_bootNum500_nkn2_070_BOOT.mat';
 disp(sprintf('[makenoisestatisticsfigs2_r2_compareNoise_R1_boot] LOADING  Results: %s',filename))

 s = load(sprintf('%s',fullfile(datadir,filename)),'noise','results','sensoryNoise', 'responseOffset');
 noise = s.noise;
 sensoryNoise = s.sensoryNoise;
 responseOffset = s.responseOffset;
 if isfield(s,'results')
  results = s.results;
 end
 clear s;
end

% plot the ratio plots
noiseRatiPlot(sensoryNoise,responseOffset,saveFig,r2,subtractionType,nBoots)

return
keyboard
% plot histograms 
noiseHist(sensoryNoise,responseOffset,saveFig,r2,subtractionType,nBoots)

% make scatter plot
noiseScatterPlot(sensoryNoise,responseOffset,saveFig,r2,subtractionType,nBoots)

% end main call



%%%%%%%%%%%%%%%%%%
% noiseRatioPlot %
%%%%%%%%%%%%%%%%%%
function noiseRatiPlot(sensoryNoise,responseOffset,saveFig,r2,subtractionType,nBoots)

% set up basic plot info:
plotInfo = setUpPlotInfo;
plotInfo.thisSymbol = {'s'  '^' 'd'};
plotInfo.XTicks = [0 .25 .5 .75 1];
plotInfo.YTicks = [0 .25 .5 .75 1];
plotInfo.title = {'v1' 'v2' 'v3' 'v4'};
plotInfo.MarkerSize = [6 6 6 6];
plotInfo.YLim = [.001 1];
plotInfo.Xlim = [.001 1];

% plot sigma
% attention figure
figurename = sprintf('attention_NoiseRatioBarPlot_sigma_r%s_%s%i',r2,subtractionType,nBoots);
h = smartfig(figurename,'reuse');
set(h,'Name',figurename);

for iObserver = size(sensoryNoise,3):-1:1
 for iVis = 1:size(sensoryNoise,4) % visual areas
  noise_f = squeeze(sensoryNoise(:,3,iObserver,iVis));
  noise_d = squeeze(sensoryNoise(:,1,iObserver,iVis));
  ratio_noise = noise_d./noise_f;
  mratio_noise(iObserver,iVis) = mean(ratio_noise);
  sdratio_noise(iObserver,iVis) = std(ratio_noise);
 end
end

mybar(mratio_noise, 'yError', sdratio_noise, ...
 'hline', 1, ...
 'yAxisMax',10, ... 
 'yAxisMin',0, ...
 'groupLabels', {'fp' 'fm', 'jg'}, ...
 'withinGroupLabels', {'V1' 'V2', 'V3', 'hV4'}, ...
 'xLabelText','Visual Areas and Observers', ...
 'yLabelText','sigma_d/sigma_f');

if saveFig
 set(h,'PaperPosition',[.25 .25 8 10.5]);
 set(h,'PaperOrientation','Portrait');
 
 disp(sprintf('Saving figure %s',figurename));
 defaultDataDir = '/data2/crfaa/crfaa/';
 figDir = 'fig_scatter_noiseTest';
 figDir = savefig(h,'defaultDataFolder',defaultDataDir,'figDir',figDir,'figName',figurename);
end


% plot offset
figurename = sprintf('attention_NoiseBarPlot_offset_r%s_%s%i',r2,subtractionType,nBoots);
h = smartfig(figurename,'reuse');
set(h,'Name',figurename);

for iObserver = size(responseOffset,3):-1:1
 for iVis = 1:size(responseOffset,4) % visual areas
  b_f = squeeze(responseOffset(:,3,iObserver,iVis));
  b_d = squeeze(responseOffset(:,1,iObserver,iVis));
  ratio_b = b_f./b_d;
  mratio_b(iObserver,iVis) = mean(ratio_b);
  sdratio_b(iObserver,iVis) = std(ratio_b);
 end
end

mybar(mratio_b, 'yError', sdratio_b, ...
 'hline', 1, ...
 'yAxisMax',1.5, ... 
 'yAxisMin',0, ...
 'groupLabels', {'fp' 'fm', 'jg'}, ...
 'withinGroupLabels', {'V1' 'V2', 'V3', 'hV4'}, ...
 'xLabelText','Visual Areas and Observers', ...
 'yLabelText','b_f/b_d');


if saveFig
 set(h,'PaperPosition',[.25 .25 8 10.5]);
 set(h,'PaperOrientation','Portrait');
 
 disp(sprintf('Saving figure %s',figurename));
 defaultDataDir = '/data2/crfaa/crfaa/';
 figDir = 'fig_scatter_noiseTest';
 figDir = savefig(h,'defaultDataFolder',defaultDataDir,'figDir',figDir,'figName',figurename);
end




%%%%%%%%%%%%%
% noiseHist %
%%%%%%%%%%%%%
function noiseHist(sensoryNoise,responseOffset,saveFig,r2,subtractionType,nBoots)

% set up basic plot info:
plotInfo = setUpPlotInfo;
plotInfo.thisSymbol = {'s'  '^' 'd'};
plotInfo.XTicks = [0 .25 .5 .75 1];
plotInfo.YTicks = [0 .25 .5 .75 1];
plotInfo.title = {'v1' 'v2' 'v3' 'v4'};
plotInfo.MarkerSize = [6 6 6 6];
plotInfo.YLim = [.001 1];
plotInfo.Xlim = [.001 1];

% plot sigma
% attention figure
figurename = sprintf('attention_sensoryNoiseHist_sigma_r%s_%s%i',r2,subtractionType,nBoots);
h = smartfig(figurename,'reuse');
set(h,'Name',figurename);

c = 1;
for iObserver = size(sensoryNoise,3):-1:1
 for iVis = 1:size(sensoryNoise,4) % visual areas
  noise_f = squeeze(sensoryNoise(:,3,iObserver,iVis));
  noise_d = squeeze(sensoryNoise(:,1,iObserver,iVis));
  % restricting the distributions to p=0.05
  q_f = 10.^quantile(log10(noise_f),[.05 .95]);
  noise_f_i = and(noise_f > q_f(1), noise_f < q_f(2));
  noise_f = noise_f(noise_f_i);
    
  q_d = 10.^quantile(log10(noise_d),[.05 .95]);
  noise_d_i = and(noise_d > q_d(1), noise_d < q_d(2));
  noise_d = noise_d(noise_d_i);
  
  subplot(size(sensoryNoise,3),size(sensoryNoise,4),c)
  myhist(noise_d,ceil(length(noise_d)/10),'b',0);
  hold on
  myhist(noise_f,ceil(length(noise_f)/10),'r',0);
  axis([.004 .16 0 (length(noise_f)/10)*4])
  set(gca,'xScale','lin')
  c = c + 1;
  title(sprintf('Sigma - Area %s Observer %i',plotInfo.title{iVis},iObserver))
 end
end


xlabel('Sigma','FontSize',plotInfo.Fsize);
ylabel('number of occurrencies','FontSize',plotInfo.Fsize);


if saveFig
 set(h,'PaperPosition',[.25 .25 8 10.5]);
 set(h,'PaperOrientation','Portrait');
 
 disp(sprintf('Saving figure %s',figurename));
 defaultDataDir = '/data2/crfaa/crfaa/';
 figDir = 'fig_scatter_noiseTest';
 figDir = savefig(h,'defaultDataFolder',defaultDataDir,'figDir',figDir,'figName',figurename);
end


% plot offset
figurename = sprintf('attention_sensoryNoiseHist_offset_r%s_%s%i',r2,subtractionType,nBoots);
h = smartfig(figurename,'reuse');
set(h,'Name',figurename);

c = 1;
for iObserver = size(sensoryNoise,3):-1:1
 for iVis = 1:size(sensoryNoise,4) % visual areas
  off_f = squeeze(responseOffset(:,3,iObserver,iVis));
  off_d = squeeze(responseOffset(:,1,iObserver,iVis));
  
  % restricting the distributions to p=0.05
  q_f = quantile(off_f,[.05 .95]);
  off_f_i = and(off_f > q_f(1), off_f < q_f(2));
  off_f = off_f(off_f_i);
    
  q_d = quantile(off_d,[.05 .95]);
  off_d_i = and(off_d > q_d(1), off_d < q_d(2));
  off_d = off_d(off_d_i);
  
  subplot(size(responseOffset,3),size(responseOffset,4),c)
  myhist(off_d,ceil(length(off_f)/10),'b',0);
  hold on
  myhist(off_f,ceil(length(off_f)/10),'r',0);
  axis([0 1 0 (length(off_f)/10)*4])
  title(sprintf('Offset - Area %s Observer %i',plotInfo.title{iVis},iObserver))
  c = c + 1;
 end
end
 
  xlabel('Response offset','FontSize',plotInfo.Fsize);
  ylabel('Numbe rof occurrences','FontSize',plotInfo.Fsize);

if saveFig
 set(h,'PaperPosition',[.25 .25 8 10.5]);
 set(h,'PaperOrientation','Portrait');
 
 disp(sprintf('Saving figure %s',figurename));
 defaultDataDir = '/data2/crfaa/crfaa/';
 figDir = 'fig_scatter_noiseTest';
 figDir = savefig(h,'defaultDataFolder',defaultDataDir,'figDir',figDir,'figName',figurename);
end


%%%%%%%%%%%%%%%%%%%%
% noiseScatterPlot %
%%%%%%%%%%%%%%%%%%%%
function noiseScatterPlot(sensoryNoise,responseOffset,saveFig,r2,subtractionType,nBoots)

% set up basic plot info:
plotInfo = setUpPlotInfo;
plotInfo.thisSymbol = {'s'  '^' 'd'};
plotInfo.XTicks = [0 .25 .5 .75 1];
plotInfo.YTicks = [0 .25 .5 .75 1];
plotInfo.title = {'v1' 'v2' 'v3' 'v4'};
plotInfo.MarkerSize = [6 6 6 6];
plotInfo.YLim = [.001 1];
plotInfo.Xlim = [.001 1];

% plot sigma
% attention figure
figurename = sprintf('attention_NoiseScatterPlot_sigma_r%s_%s%i',r2,subtractionType,nBoots);
h = smartfig(figurename,'reuse');
set(h,'Name',figurename);

loglog([.007,.007;.08,.08],[.007,.007;.08,.08],'k-','Color',plotInfo.XYColor)
hold on;
for iObserver = size(sensoryNoise,3):-1:1
 for iVis = 1:size(sensoryNoise,4) % visual areas
  noise_f = squeeze(sensoryNoise(:,3,iObserver,iVis));
  noise_d = squeeze(sensoryNoise(:,1,iObserver,iVis));
  mf_noise = mean(noise_f);
  md_noise = mean(noise_d);
  sdf_noise = std(noise_f);
  sdd_noise = std(noise_d);
  myerrorbar(md_noise,mf_noise, 'yError', sdf_noise, 'xError', sdd_noise, ...
   'Symbol',plotInfo.thisSymbol{iObserver}, ...
   'LineWidth',plotInfo.LineWidth, ...
   'Color',plotInfo.thisColor{1}{iVis}, ...
   'MarkerFaceColor',plotInfo.MarkerFaceColor{1}{iVis}, ...
   'MarkerEdgeColor','w', ...
   'MarkerSize', plotInfo.MarkerSize(iObserver));
 end
end

title(sprintf('Noise: Sigma'))
xlabel('Distributed','FontSize',plotInfo.Fsize);
ylabel(sprintf('Attended'),'FontSize',plotInfo.Fsize);

set(gca,...
 'FontName','Helvetica','FontSize',plotInfo.Fsize, ...
 'PlotBoxAspectRatio',plotInfo.PlotBoxAspectRatio, ...
 'XLim', [.007 .18],'YLim', [.007 .18],...
 'LineWidth',plotInfo.LineWidth,'TickLength',plotInfo.TickLength, ...
 'XTick', [.007  .012 .023  .045 .09 .18], 'XTickLabel', [.007  .012 .023 .045 .09 .18] ,...
 'YTick', [.007  .012 .023  .045 .09 .18], 'YTickLabel', [.007  .012 .023  .045 .09 .18] ,...
 'XColor',plotInfo.XYColor,'YColor',plotInfo.XYColor,'Box','off');

myaxisScatterLog;

if saveFig
 set(h,'PaperPosition',[.25 .25 8 10.5]);
 set(h,'PaperOrientation','Portrait');
 
 disp(sprintf('Saving figure %s',figurename));
 defaultDataDir = '/data2/crfaa/crfaa/';
 figDir = 'fig_scatter_noiseTest';
 figDir = savefig(h,'defaultDataFolder',defaultDataDir,'figDir',figDir,'figName',figurename);
end


% plot offset
figurename = sprintf('attention_NoiseScatterPlot_offset_r%s_%s%i',r2,subtractionType,nBoots);
h = smartfig(figurename,'reuse');
set(h,'Name',figurename);

plot([.1,.1;.8,.8],[.1,.1;.8,.8],'k-','Color',plotInfo.XYColor)
hold on;
for iObserver = size(responseOffset,3):-1:1
 for iVis = 1:size(responseOffset,4) % visual areas
  rOffset_f = squeeze(responseOffset(:,3,iObserver,iVis));
  rOffset_d = squeeze(responseOffset(:,1,iObserver,iVis));
  mf_off = mean(rOffset_f);
  md_off = mean(rOffset_d);
  sdf_off = std(rOffset_f);
  sdd_off = std(rOffset_d);
 myerrorbar(md_off,mf_off, 'yError', sdf_off, 'xError', sdd_off, ...
  'Symbol',plotInfo.thisSymbol{iObserver}, ...
  'LineWidth',plotInfo.LineWidth, ...
  'Color',plotInfo.thisColor{1}{iVis}, ...
  'MarkerFaceColor',plotInfo.MarkerFaceColor{1}{iVis}, ...
  'MarkerEdgeColor','w', ...
  'MarkerSize', plotInfo.MarkerSize(iObserver));
 end
end
 
 title(sprintf('Noise: offset'))
  xlabel('Distributed','FontSize',plotInfo.Fsize);
  ylabel('Attended','FontSize',plotInfo.Fsize);
set(gca,...
 'FontName','Helvetica','FontSize',plotInfo.Fsize, ...
 'PlotBoxAspectRatio',plotInfo.PlotBoxAspectRatio, ...
 'XLim', [0 1],'YLim', [0 1],...
 'LineWidth',plotInfo.LineWidth,  'TickLength',plotInfo.TickLength, ...
 'XTick', [0 .1 .2 .3 .4 .5 .6 .7 .8], 'XTickLabel', [0 .1 .2 .3 .4 .5 .6 .7 .8] ,...
 'YTick', [0 .1 .2 .3 .4 .5 .6 .7 .8], 'YTickLabel', [0 .1 .2 .3 .4 .5 .6 .7 .8] ,...
 'XColor',plotInfo.XYColor,'YColor',plotInfo.XYColor,'Box','off');

% this is the origianl call, this function was lost at a certain point
% i am not making it now unless we need it: myaxisScatterLin;
myaxisScatterLin;

if saveFig
 set(h,'PaperPosition',[.25 .25 8 10.5]);
 set(h,'PaperOrientation','Portrait');
 
 disp(sprintf('Saving figure %s',figurename));
 defaultDataDir = '/data2/crfaa/crfaa/';
 figDir = 'fig_scatter_noiseTest';
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
function [d meanC events filename doMean] = loadData(observer,visualAreas,defaultDataFolder,adaptation,r2,subtractionType,boot)

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
 keyboard
 clear data;
 events = [];
 doMean = 1;

else % individual observers

 if ~(strcmpi('fm',observer) && adaptation == 100)
  [filename vars2load] = makeFileName(observer,adaptation,visualAreas,defaultDataFolder,r2,subtractionType);
  data = load(sprintf('%s',filename),vars2load{1}, vars2load{2});
  
  % now make the average contrast
  meanC = computedMeanC(data.events,.5);
  meanC.quantileAcontrast = 100*meanC.quantileAcontrast';
  meanC.quantileDcontrast = 100*meanC.quantileDcontrast';
  meanC.pedC = 100*meanC.pedC;
  
  d.A_amplitude   = data.d2.amplitude(1:8,boot)';
  d.U_amplitude   = data.d2.amplitude(9:16,boot)';
  d.Dt_amplitude  = data.d2.amplitude(17:24,boot)';
  d.Dnt_amplitude = data.d2.amplitude(25:32,boot)';
  
  d.A_ste         = data.d2.amplitudeSTE(1:8,boot)';
  d.U_ste         = data.d2.amplitudeSTE(9:16,boot)';
  d.Dt_ste        = data.d2.amplitudeSTE(17:24,boot)';
  d.Dnt_ste       = data.d2.amplitudeSTE(25:32,boot)';
  
  events = data.events;
 end
 doMean = 0;
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
 filename =  fullfile(defaultDataFolder,sprintf('ADAPTATION-EFFECT-%s-%s-%s.mat', visualArea,r2,subtractionType));
 vars2load = {'d' 'meanC'};
 
else
 % choose the type of analysis files to load
 vars2load = {'d2' 'events'};
 this_date = '2010-20-July';
 filename = fullfile(defaultDataFolder,sprintf('%s_A%s_%s_%s_era_R1_subtr%s_%s.mat',upper(observer), ...
                     num2str(adaptation),this_date,visualArea, subtractionType,visualArea));
  
end


%%%%%%%%%%%%%%%%%%%%%%%%
% makeIndAttentionPlot %
%%%%%%%%%%%%%%%%%%%%%%%%
function makeIndAttentionPlot(results,noise,observer,fitType,saveFig,plotLog,r2,subtractionType)
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
   for itest = 1:2
   % open up a new figure for each observer:
   figurename = deblank(sprintf('noisetest_TvC_testType%s_OBS%s_%s_A0',tests{itest},observer{iObserver},areas{iArea}));
   h = smartfig(figurename,'reuse');
   
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
    TvC_low  = .3*ones(size(results{iObserver, iArea, itest}.behavior.tvc.thisTvCste(:,1)));
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
   set(gca, ...
    'XScale','log', ...
    'XTick', [.00875 .0175 .035 .07 .14 .28 .56 .84 1], ...
    'XLim',[.00875 1], ...
    'XTickLabel', [0 .0175 .035 .07 .14 .28 .56 .84 1], ...
    'YScale','log', ...
    'YLim',[0.0021875  1], ...
    'YTick',[0.0021875 .004375 .00875 .0175 .035 .07 .14 .28],...
    'YTickLabel', [0.0021875 .004375 .00875 .0175 .035 .07 .14 .28]);
   
   axis('square')
   myaxisTvC;
   drawnow
   
   if saveFig
    set(h,'PaperPosition',[.25 .25 8 10.5]);
    set(h,'PaperOrientation','Portrait');
    
    AxLabel = 'LogLog';
    figurename = sprintf('%s_%s_r%s_%s',figurename,AxLabel,r2,subtractionType);
    disp(sprintf('Saving figure %s',figurename));
    defaultDataFolder = '/data2/crfaa/crfaa/';

    
    figDir = 'fig_tvc_noiseTest';
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
   figurename = deblank(sprintf('noisetest_CRF_testType%s_OBS%s_%s_A0',tests{itest},observer{iObserver},areas{iArea}));
   h = smartfig(figurename,'reuse'); set(h,'Name',figurename);
   
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
   figurename = deblank(sprintf('noisetest_CRF_testType%s_OBS%s_%s_A0',tests{itest},observer{iObserver},areas{iArea}));
   h = smartfig(figurename,'reuse'); set(h,'Name',figurename);

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
     'XScale','log', ...
     'YScale','lin', ...
     'XTick', [0.0021875 .004375 .00875 .0175 .035 .07 .14 .28 .56 .84 1], ...
     'XTickLabel', [0 .004375 .00875 .0175 .035 .07 .14 .28 .56 .84 1], ...
     'YTick', [.25 .5 .75 1 1.25] ,...
     'YTickLabel', [.25 .5 .75 1 1.25] ,...
     'XLim',[0.0021875 1], ...
     'YLim',[.25 1.45]);
    axis('square')
    myaxisCRF;
   else
    AxLabel = 'lin';
    set(gca,...
     'XScale','lin', ...
     'YScale','lin', ...
     'XTick', [.004375 .00875 .0175 .035 .07 .14 .28 .56 .84 1], ...
     'XTickLabel', [0 .00875 .0175 .035 .07 .14 .28 .56 .84 1], ...
     'YTick', [.25 .5 .75 1 1.25] ,...
     'YTickLabel', [.25 .5 .75 1 1.25] ,...
     'XLim',[0 1], ...
     'YLim',[.25 1.45]);
    axis('square')
    myaxisCRFlin;
   end

   if saveFig
    set(h,'PaperPosition',[.25 .25 8 10.5]);
    set(h,'PaperOrientation','Portrait');
    figurename = sprintf('%s_%s_r%s_%s',figurename,AxLabel,r2,subtractionType);
    disp(sprintf('Saving figure %s',[defaultDataFolder,figurename]));
    defaultDataFolder = '/data2/crfaa/crfaa/';


    figDir = 'fig_crf_noiseTest';
    figDir = savefig(h,'figName',figurename,'defaultDataFolder', defaultDataFolder,'figDir',figDir);
   end
  end
 end
end

