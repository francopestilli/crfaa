function r = measureWeberFraction(varargin)
% function r = measureWeberFraction(varargin)
%
% fits a power law to the log-linearportion of the TvC in the attention
% paper.
%
% data files are:
%  for the average          filename = '/Volumes/data/riken/crfaa/fmridata/crfaa/data_used_files/12-Jul-2009averageTvC.mat';
%  for individual observers filename = '/Volumes/data/riken/crfaa/fmridata/crfaa/31-Jul-2009CompareNoiseFit_nkn2.mat';
%
% and figures are saved here:
% figDir = 'fig_weber_fraction';
% defaultDataFolder = '/Users/frakkopesto/data/riken/crfaa/fmridata/crfaa/';
%
% franco pestilli 2009.07.31

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
observer  =[];
visualarea=[];
adaptation=[];
fitType   =[];     % go from tvc to crf OR from crf to tvc
savePlots =[];     % 1 saves the plot in the current dierectory
dispFit   =[];     % if 1 the plot is displayed
numBoots  =[];     % number of bootstrap samples
saveFig   =[];     % if 1 saves a figure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set default variables:
getArgs(varargin,{...
 'observer',{'fm' 'fp' 'jg' 'average' }, ...
 'adaptation',[0 28 100], ...
 'dispFit',0, ...
 'saveFig',1, ...
 'numBoots', 10000, ...
 'filename',[date,'allNoise.mat']});

% default folder to load the data:
defautDataFolder = '/Users/frakkopesto/data/riken/crfaa/fmridata/crfaa/data_used_files/';


for iObserver = 1:length(observer)
 disp(sprintf('[measureWeberFraction] ADP: %i',observer{iObserver}))
 for iAdapt = 1:length(adaptation)
  disp(sprintf('[measureWeberFraction] ADP: %i',adaptation(iAdapt)))
  
  % load the proper file.
  [t] = loadData(observer{iObserver},'v2',defautDataFolder,adaptation(iAdapt));
  
  if ~isempty(t)
  % get the threshold (data)
  tvc = getDeltaContrast(t);

  % now fit a power law to the actual data:
  for a = 1:size(tvc.thresholds,1)
   w(iObserver,iAdapt,a,1:2) = fitPowerFun(tvc.pedestals,tvc.thresholds(a,:));
  end
  
  % now generate numBoot TvC functions
  % resapling the data from a distribution centered at the observed
  % threshold with the standard deviation given by the variability in th
  % data (tvc.thresholds_ste):
  numCueConditions = size(tvc.thresholds,1);
  numContrasts = size(tvc.thresholds,2);
  bt_threshold = nan.*ones(numCueConditions,numContrasts,numBoots);
  for a = 1:numCueConditions
   disp(sprintf('[measureWeberFraction] Bootstrap cue condition %i',a))

   % (1) compute the interval around the mean threshold
   %     given the ste of the data:
   for c = 1:numContrasts
    low = log10(tvc.thresholds(a,c)) - tvc.thresholds_ste(a,c);
    high = log10(tvc.thresholds(a,c)) + tvc.thresholds_ste(a,c);
    bt_threshold(a,c,:) = 10.^(low + (high-low).*rand(1,numBoots));
   end
   % fit a power function here and get the w parameter distribution
   % directly OR figure out the indexes of this matrices
   for bt = 1:numBoots
%       disp(sprintf('[measureWeberFraction] Fit power function to bootstrap #%i',bt))
      thisTvC = squeeze(bt_threshold(a,:,bt));
    bt_w(iObserver,iAdapt,a,bt,1:2) = fitPowerFun(tvc.pedestals,thisTvC);
   end
   % SAVE the result
   this_w.mean(iObserver,iAdapt,a,1) = squeeze(mean(bt_w(iObserver,iAdapt,a,:,1),4)); 
   this_w.mean(iObserver,iAdapt,a,2) = squeeze(mean(bt_w(iObserver,iAdapt,a,:,2),4));
   this_w.std(iObserver,iAdapt,a,1)  = squeeze(std(bt_w(iObserver,iAdapt,a,:,1), [], 4));
   this_w.std(iObserver,iAdapt,a,2)  = squeeze(std(bt_w(iObserver,iAdapt,a,:,2), [], 4));

  end

  smartfig(sprintf('ThisPowerFitA%i',iAdapt),'reuse');
  plotPowerFunFit(tvc,squeeze(this_w.mean(iObserver,iAdapt,:,:)),squeeze(this_w.std(iObserver,iAdapt,:,:)))
  end
 end
end

w(w==0)=nan;

% compute the mean and ste of the power exponent:
r.w = w;
thismean = squeeze(nanmean(w,1));
thisste = squeeze(nanstd(w,1))./sqrt(length(observer));

r.offset_mean = squeeze(thismean(:,:,1));
r.offset_ste = squeeze(thisste(:,:,1));
r.w_mean = squeeze(thismean(:,:,2));
r.w_ste = squeeze(thisste(:,:,2));

% plot the results
plotBarExponent(r,saveFig)




%%%%%%%%%%%%%%%%%%%
% plotBarExponent %
%%%%%%%%%%%%%%%%%%%
function plotBarExponent(r,saveFig)

figurename = 'powerExponent';
h = smartfig(figurename,'reuse');
mybar(r.w_mean, ...
 'yError',r.w_ste, ...
 'groupLabels',{'0','28' '100'}, ...
 'withinGroupLabels',{'Distributed cue, target' 'Focal cue, target'}, ...
 'yAxisMin',.5, ...
 'yAxisMax',1, ...
 'xLabelText','Adapter contrast (%)','dispValues',1);

if saveFig
 set(h,'PaperPosition',[.25 .25 8 10.5]);
 set(h,'PaperOrientation','Portrait');
 
 disp(sprintf('Saving figure %s',figurename));
 figDir = 'fig_weber_fraction';
 defaultDataFolder = '/Users/frakkopesto/data/riken/crfaa/fmridata/crfaa/';
 figDir = savefig(h,'figName',figurename,'defaultDataFolder', defaultDataFolder,'figDir',figDir);
end

% plot the vertical offset on thrContrast 
figurename = 'offestExponent';
h = smartfig(figurename,'reuse');
mybar(r.offset_mean, ...
 'yError',r.offset_ste, ...
 'groupLabels',{'0','28' '100'}, ...
 'withinGroupLabels',{'Distributed cue, target' 'Focal cue, target'}, ...
 'yAxisMin',0, ...
 'yAxisMax',1, ...
 'xLabelText','Adapter contrast (%)','dispValues',1);

if saveFig
 set(h,'PaperPosition',[.25 .25 8 10.5]);
 set(h,'PaperOrientation','Portrait');
 
 disp(sprintf('Saving figure %s',figurename));
 figDir = 'fig_weber_fraction';
 defaultDataFolder = '/Users/frakkopesto/data/riken/crfaa/fmridata/crfaa/';
 figDir = savefig(h,'figName',figurename,'defaultDataFolder', defaultDataFolder,'figDir',figDir);
end


%%%%%%%%%%%%%%%%%%%
% plotPowerFunFit %
%%%%%%%%%%%%%%%%%%%
function plotPowerFunFit(tvc,w,wstd)
plotColor = { [.02 .1 .9] [.9 .1 .02]};
col = {'c' 'm'};
for a = 1:2
 % get the powerfunction:;
 low_t = powerLaw(w(a,:)-wstd(a,:),tvc.pedestals);
 t = powerLaw(w(a,:),tvc.pedestals);
 high_t = powerLaw(w(a,:)+wstd(a,:),tvc.pedestals);

 % plot the fit:
 myerrorbar(tvc.pedestals,t, ...
   'yLow',low_t, ...
   'yHigh',high_t, ...
   'Symbol','-',...
   'Color',col{a});
 
 % plot the data:
 myerrorbar(tvc.pedestals,tvc.thresholds(a,:), ...
   'yLow',0, ...
   'yHigh',0, ...
   'Symbol','o',...
   'MarkerFaceColor',plotColor{a}, ...
   'MarkerEdgeColor','k',...
   'MarkerSize',10, ...
   'Color',plotColor{a});
end

% axis formatting:
set(gca, ...
 'XScale','log', ...
 'XTick', [.0175 .035 .07 .14 .28], ...
 'XLim',[.0175 1], ...
 'XTickLabel', [.0175 .035 .07 .14 .28], ...
 'YScale','log', ...
 'YLim',[.002175   .56], ...
 'YTick',[.002175 .004375 .00875  .0175  .035 .07 .14 .28 .56],...
 'YTickLabel', [.002175 .004375 .00875 .0175 .035 .07 .14 .28 .56]);
axis('square')
drawnow

%%%%%%%%%%%%%%%
% fitPowerFun %
%%%%%%%%%%%%%%%
function w = fitPowerFun(pedestals,thresholds)
% the power function has the form:
% threshold = pedestals^omega
% taken from gorea and sagi, NN 2001, p. 1147

% optimization parameters
maxIter = 1000;
MaxFunEvals = 1000;
optimParams = optimset('LevenbergMarquardt','on','MaxIter',maxIter,'MaxFunEvals',MaxFunEvals,'Display','off','Diagnostics','off');

% init the noise params
initParams = [.1 .5];
minParams = [0 0];
maxParams = [1 1];

[w resnorm residual exitflag output lambda jacobian] = lsqnonlin(@errPowerFun, ...
 initParams,minParams,maxParams,optimParams,pedestals,thresholds);


%%%%%%%%%%%%%%%
% errPowerFun %
%%%%%%%%%%%%%%%
function err = errPowerFun(params,pedestals,thresholds)

% estimated threshold given the current exponent:
t_hat = powerLaw(params,pedestals);

err = log10(thresholds) - log10(t_hat);


%%%%%%%%%%%%
% powerLaw %
%%%%%%%%%%%%
function t = powerLaw(w,c)
% the power function has the form:
% threshold = pedestals^omega
% taken from gorea and sagi, NN 2001, p. 1147

t = w(1) .* c.^w(2);


%%%%%%%%%%%%%%%
% getContrast %
%%%%%%%%%%%%%%%
function tvc = getDeltaContrast(t)

if ~isempty(t)
 % drop 0 pedestal and last threshold
 tvc.thresholds    = t.threshold(:,2:end-2);
 tvc.thresholds_ste = t.threshold_ste(:,2:end-2);
 tvc.pedestals   = [.0175 .035 .07 .14 .28];
 
else
 keyboard
end


%%%%%%%%%%%%
% loadData %
%%%%%%%%%%%%
function [tvc filename] = loadData(observer,visualAreas,defautDataFolder,adaptation)

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
switch observer
 case {'avrg' 'average'}
  % (0) load data file
  [filename vars2load] = makeFileName(observer,[],visualAreas,defautDataFolder);
  [data] = load(sprintf('%s',filename),vars2load{1},vars2load{2});
  events = [];
  
  tvc.threshold = squeeze(data.meanTvC(adp,:,:));
  tvc.threshold_ste = squeeze(data.meanTvC(adp,:,:));
  
 case {'fm' 'FM' 'fp' 'FP' 'jg' 'JG'}
 % this is from compareNoise the fucntio who made the files being loaded 
 % 'observer',{'avrg','jg' 'fm' 'fp'}, ...
  if strcmpi(observer,'fm')
   o = 3;
  elseif  strcmpi(observer,'fp')
   o = 4;
  elseif strcmpi(observer,'jg')
   o = 2;
  end
  
  if ~(strcmpi('fm',observer) && adaptation == 100)
   [filename vars2load] = makeFileName(observer,adaptation,visualAreas,defautDataFolder);
   data = load(sprintf('%s',filename),vars2load{1});
   tvc.threshold(1,1:8) = data.results{o,1,1}.behavior.tvc.thisTvC;
   tvc.threshold(2,1:8) = data.results{o,1,2}.behavior.tvc.thisTvC;
   tvc.threshold_ste(1,1:8) = data.results{o,1,2}.behavior.tvc.thisTvCste;
   tvc.threshold_ste(2,1:8) = data.results{o,1,2}.behavior.tvc.thisTvCste;
  
  else
   tvc = [];
  end
  
 otherwise
  keyboard
end


%%%%%%%%%%%%%%%%
% makeFileName %
%%%%%%%%%%%%%%%%
function [filename vars2load] = makeFileName(observer,adaptation,visualArea,defautDataFolder)

switch observer
 case {'avrg' 'average'}
  filename = '/Volumes/data/riken/crfaa/fmridata/crfaa/data_used_files/12-Jul-2009averageTvC.mat';
  vars2load = {'meanTvC' 'steTvC'};
  
 case {'fm' 'FM' 'fp' 'FP' 'jg' 'JG'}
  % load a file result of any dprimefit3.m so that the
  % thresholds are already computed
  
  filename = '/Volumes/data/riken/crfaa/fmridata/crfaa/31-Jul-2009CompareNoiseFit_nkn2.mat';
  vars2load = {'results'};
  
 otherwise
  keyboard
end

