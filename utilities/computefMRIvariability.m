function err = computefMRIvariability(varargin)
% function err = computefMRIvariability(varargin)
%
% this function computes the variability of the variability of the fMRI
% response. It takes the STD analysis for each observer, bootstraps
% the data, and computes the average.
%
% These bootstrapped averages are used to estimate the variability of the
% STD of the fMRI signal (figure 6 of the attention paper).
%
% the function also computes the probability that the observerved difference in 
% fMRI STD between the Focal cue, target and Distributed cue, targed
% conditions is not due to chance. This is the test reported in the
% attention paper.
%
% franco pestilli 2007.07.29
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% visualAreas = []; % selects which visual area to use (to load the file)                                                                     %
% analysisType= []; % choose which analysis files to load, standard is the main results, 1st/2nd and corr/incorrect are the other options     %
% saveData          = []; % if 1, saves the results                                                                                           %
% observer          = []; % selects the observers to load the data (average also)                                                             %
% adaptation        = []; % which adaptation contiion for the current observer, default for the attention paper is '0'                        %
% numBootStraps     = []; % number of bootstrap samples take                                                                                  %
% restrictToPositive = [];% restricts the amplitude of the responses only to the positive ones                                                %
% normalizeByTrialNum = []; % normalize the estimated standardeviation by the number fo trials in the condition (basically compute the ste)   %
% saveFig           = []; % saves the figure just plotted (check functionality it has been inherited)                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
visualAreas = []; % selects which visual area to use (to load the file)                                                                     %
analysisType= []; % choose which analysis files to load, standard is the main results, 1st/2nd and corr/incorrect are the other options     %
saveData          = []; % if 1, saves the results                                                                                           %
observer          = []; % selects the observers to load the data (average also)                                                             %
adaptation        = []; % which adaptation contiion for the current observer, default for the attention paper is '0'                        %
numBootStraps     = []; % number of bootstrap samples take                                                                                  %
restrictToPositive  = []; % restricts the amplitude of the responses only to the positive ones                                                %
normalizeByTrialNum = []; % normalize the estimated standardeviation by the number fo trials in the condition (basically compute the ste)   %
saveFig             = []; % saves the figure just plotted (check functionality it has been inherited)                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


getArgs(varargin,{...
 'visualAreas',{'v1','v2','v3','v4'}, ... 
 'analysisType', 'standard', ...
 'saveData',1, ...
 'observer',{'fp' 'fm' 'jg'}, ...
 'adaptation',0, ...
 'numBootStraps', 100000, ...
 'restrictToPositive', 0, ...
 'normalizeByTrialNum', 0, ...
 'saveFig',1});

global defautDataFolder;

% set figure filder (to save) and data folder (to save):
figDir = ['figs_nestedAttTests_',date];
if isdir('/Users/frakkopesto/data/riken/crfaa/fmridata/crfaa')
  defautDataFolder = '/Users/frakkopesto/data/riken/crfaa/fmridata/crfaa/data_used_files/';
else
  defaultDataFolder = '/Volumes/frakkopesto_HomeDir/data/riken/crfaa/fmridata/crfaa/data_used_files/';
end
 
%%%%%%%%%%%%%
% load data %
%%%%%%%%%%%%%
for o = 1:length(observer)      % do the whole thing for each observer
 for v = 1:length(visualAreas) % v1-v4 
  d{o,v} = loadData(observer{o},visualAreas{v},adaptation,analysisType);
 end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% bootstrap the STD estimates %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numConditions = length(d{1,1}.groupSTD.amplitude);
bootSTD = nan.*ones(length(observer),length(visualAreas),numConditions,numBootStraps);
for o = 1:length(observer)      % do the whole thing for each observer
 for v = 1:length(visualAreas)  % v1-v4
  % do a bootstrap for each condition:
  for icond = 1:numConditions
   allAmplitudes = d{o,v}.groupSTD.amplitude{icond};% here make a test for the estimates of STD on only positive trials...
   numAmplitudes = length(allAmplitudes);
   disp(sprintf('[computefMRIvariability] now boot strapping observer %s visual area %s condition: %i, num of Amplitudes: %s', ...
                 observer{o},visualAreas{v},icond,num2str(numAmplitudes)));

   % save the STDs estimated from the original data:
   originalSTD(o,v,icond) = d{o,v}.groupSTD.amplitudeSTD(icond);
   
   % randsample with replacement the amplitudes:
   for bt = 1:numBootStraps
    withReplacement = 1; % if '1' rand samples are taken with replacement
    index = randsample(numAmplitudes,numAmplitudes,withReplacement);
    bootSTD(o,v,icond,bt) = std(allAmplitudes(index));
   end
  end
 end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Compute a permutation test %
%   This test is used in the   %
%   figure for statistical     %
%   significance               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (1) compute the empirical difference between 
%     <focal cue, target> and <distributed cue, target>
avrg_originalSTD    = squeeze(mean(originalSTD,1));
empiricalDifference = avrg_originalSTD(:,3) - avrg_originalSTD(:,1);

% (2) compute the difference between the std for <focal cue, target> and <distributed cue, target> for each bootstrap sample.
avrg_stdbootSTD_A   = squeeze(mean(squeeze(bootSTD(:,:,1,:)),1));
avrg_stdbootSTD_D   = squeeze(mean(squeeze(bootSTD(:,:,3,:)),1));
bootstrapDifference = avrg_stdbootSTD_A - avrg_stdbootSTD_D;


% (4) now find the probability (the quantile) that corresponds to the 
%     <empiricalDifference> given the distribution <bootstrappedDifference> 
for v = 1:length(visualAreas)  % v1-v4
 p(v) = 1 - fitBootStrapP(empiricalDifference(v),bootstrapDifference(v,:));
end

disp(sprintf('[computefMRIvariability]\nThe probability that the observed difference in STD between\n<Focal cue, Target> and <Distributed cue, Target> is:\n[%s]',num2str(p)));


%%%%%%%%%%%%%%%%%%%%
% make a bar graph %
%%%%%%%%%%%%%%%%%%%%
saveFig = 1;
bootSTD = squeeze(bootSTD(:,1:3,:,:));
originalSTD = originalSTD(:,1:3,:);
fmriNoiseBarGraph(bootSTD,originalSTD,saveFig)


% plot the results of the bootstrap in individual figures %
% plotEachBootStrap(bootSTD,originalSTD,observer,visualAreas,numConditions,numBootStraps,normalizeByTrialNum);


% plot the results of the bootstrap in individual figures %
% plotAvrgBootStrap(bootSTD,originalSTD,observer,visualAreas,numConditions,numBootStraps,normalizeByTrialNum);


% make a scatter plot like the noise estimates of the model %
% noiseScatterPlot(bootSTD,originalSTD,saveFig)
 

keyboard

%%%%%%%%%%%%%%%%%
% END main call %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%
% fitBootStrapP %
%%%%%%%%%%%%%%%%%
function p = fitBootStrapP(observation,bootStrapDistribution)
% find the probability value that corresponds to the observed difference between
% distributed- and attended-target fMRI STDs 

% optimization parameters
maxIter = inf;
MaxFunEvals = inf;
optimParams = optimset('LevenbergMarquardt','on','MaxIter',maxIter,'MaxFunEvals',MaxFunEvals,'Display','off','Diagnostics','off');

% init the noise params
initParams = 1;
minParams = 0;
maxParams = 100;

[probability resnorm residual exitflag output lambda jacobian] = lsqnonlin(@getErrProbability,   ...
                                                                 initParams,minParams,maxParams, ...
                                                                 optimParams,observation,bootStrapDistribution);
p = probability/100;
                                                  

%%%%%%%%%%%%%%%%%%%%%
% getErrProbability %
%%%%%%%%%%%%%%%%%%%%%
function err = getErrProbability(p, observation, bootStrapDistribution)
% computes the error difference between the observed data and the
% percentile given the empirical distribution at the current probability
% value (p)

currentPercentile = prctile(bootStrapDistribution,p);
err = observation - currentPercentile;
disp(sprintf('[getErrProbability] current error <%s> - current percentile <%s> - observer percentile <%s>.', ...
              num2str(err),num2str(currentPercentile),num2str(observation)));


%%%%%%%%%%%%%%%%%%%%%
% fmriNoiseBarGraph %
%%%%%%%%%%%%%%%%%%%%%
function fmriNoiseBarGraph(bootSTD,originalSTD,saveFig)
% function noiseScatterPlot(noise)
% this function is able to plot bar-graphs across visual areas bu tit is
% currently not used.
%

if size(bootSTD,2) == 1
 bootSTD = squeeze(bootSTD);
 originalSTD = squeeze(originalSTD);
 doSingelROI = 1;
else
 doSingelROI = 0;
end

figurename = sprintf('fmriNoiseBar');
label = {'Focal cue, target','Focal cue, non-target','Distributed cue, target','Distributed cue, non-target'};
ylabel = sprintf('STD of fMRI response\n (%s signal change)','%');

% standard deviation averaged across observers and visual areas
meanSTD = squeeze(mean(bootSTD,1));% average across observers

% standard deviation averaged across observers and visual areas
originalSTD = squeeze(mean(originalSTD,1)); % use the original data
% originalSTD = squeeze(mean(meanSTD,3)); % use the bootstrap average

if doSingelROI
 visualArea = {'Single'};
 % 95% confidence intervals of the STD
 CI_meanSTD = zeros(1,size(meanSTD,1));
 for icond = 1:size(meanSTD,1)
  this_meanSTD = squeeze(meanSTD(icond,:));
  CI_meanSTD(icond) = quantile(this_meanSTD,.95);
  CI_meanSTD(icond) = CI_meanSTD(icond) - originalSTD(icond);
 end
 
else
 visualArea = {'v1' 'v2' 'v3' 'v4'};
 % 95% confidence intervals of the STD
 CI_meanSTD = zeros(size(meanSTD,1),size(meanSTD,2));
 for v = 1:size(meanSTD,1)
  for icond = 1:size(meanSTD,2)
   this_meanSTD = squeeze(meanSTD(v,icond,:));
   CI_meanSTD(v,icond) = quantile(this_meanSTD,.95);
   CI_meanSTD(v,icond) = CI_meanSTD(v,icond) - originalSTD(v,icond);
  end
 end
end


for v = 1:size(CI_meanSTD,1)
 h = smartfig([figurename,visualArea{v}],'reuse');
 mybar(originalSTD(v,:), ...
  'yError',CI_meanSTD(v,:), ...
  'withinGroupLabels',label, ...
  'yAxisMin',0.5, ...
  'yAxisMax',2, ...
  'yLabelText',ylabel,'dispValues',1);
 
 if saveFig
  set(h,'PaperPosition',[.25 .25 8 10.5]);
  set(h,'PaperOrientation','Portrait');
  
  disp(sprintf('Saving figure %s',[figurename,visualArea{v}]));
  figDir = 'fig_attention_fMRInoise';
  defaultDataFolder = '/Users/frakkopesto/data/riken/crfaa/fmridata/crfaa/';
  figDir = savefig(h,'figName',[figurename,visualArea{v}],'defaultDataFolder', defaultDataFolder,'figDir',figDir);
 end
end

if ~doSingelROI
 h = smartfig(figurename,'reuse');
 mybar(mean(originalSTD,1), ...
  'yError',mean(CI_meanSTD,1), ...
  'withinGroupLabels',label, ...
  'yAxisMin',0.5, ...
  'yAxisMax',2, ...
  'yLabelText',ylabel,'dispValues',1)
end

if saveFig
 set(h,'PaperPosition',[.25 .25 8 10.5]);
 set(h,'PaperOrientation','Portrait');
 
 disp(sprintf('Saving figure %s',figurename));
 figDir = 'fig_attention_fMRInoise';
 defaultDataFolder = '/Users/frakkopesto/data/riken/crfaa/fmridata/crfaa/';
 figDir = savefig(h,'figName',figurename,'defaultDataFolder', defaultDataFolder,'figDir',figDir);
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


%%%%%%%%%%%%%%%%%%%%%
% plotAvrgBootStrap %
%%%%%%%%%%%%%%%%%%%%%
function plotAvrgBootStrap(bootSTD, originalSTD,observer,visualAreas,numConditions,numBootStraps,normalizeByTrialNum)
conditions = {'Focal cue, target','Focal cue, non-target','Distributed cue, target','Distributed cue, non-target'};

% average across visual areas:
% h = smartfig(sprintf('AvrgBootstraphistVisArea'),'reuse');
% VavrgBootSTD = squeeze(mean(bootSTD,2));

% average across observers
h = smartfig(sprintf('AvrgBootstraphistObservers'),'reuse');
OavrgBootSTD = squeeze(mean(bootSTD,1));
OavrgOriginalSTD = squeeze(mean(originalSTD,1));

if (length(size(OavrgBootSTD)) == 2)
 for icond = 1:numConditions
  subplot(4,1,icond),
  thisBT = squeeze(OavrgBootSTD(icond,:));
  hist(thisBT'), hold on,
  
  % set the axis limits proptperly
  axis([.4 1.5 0 size(OavrgBootSTD,2)/3])
  
  vline(OavrgOriginalSTD(icond),'g--');
  vline(mean(thisBT(:)),'r--');
  legend(conditions{icond});
  title(sprintf('%s',visualAreas{1}));
 end
 
else % multiple rois
 c = 0;
 for icond = 1:numConditions
  for v = 1:length(visualAreas) % v1-v4
   c = c + 1;
   subplot(size(OavrgBootSTD,2),size(OavrgBootSTD,1),c),
   thisBT = squeeze(OavrgBootSTD(v,icond,:));
   hist(thisBT), hold on,
   
   % set the axis limits proptperly
   if v == 4
    axis([.5 2.7 0 size(OavrgBootSTD,3)/3])
   else
    axis([.4 1.5 0 size(OavrgBootSTD,3)/3])
   end
   
   vline(OavrgOriginalSTD(v,icond),'g--');
   vline(mean(thisBT(:)),'r--');
   legend(conditions{icond});
   title(sprintf('%s',visualAreas{v}));
  end
 end
end


%%%%%%%%%%%%%%%%%%%%
% noiseScatterPlot %
%%%%%%%%%%%%%%%%%%%%
function noiseScatterPlot(bootSTD,originalSTD,saveFig)

% set up basic plot info:
plotInfo = setUpPlotInfo;
plotInfo.thisSymbol = {'o' 's'  '^' 'd'};
plotInfo.XTicks = [0 .25 .5 .75 1];
plotInfo.YTicks = [0 .25 .5 .75 1];
plotInfo.title = {'v1' 'v2' 'v3' 'v4'};
plotInfo.MarkerSize = [10 10 10 10 10];
plotInfo.YLim = [.001 1];
plotInfo.Xlim = [.001 1];

% plot sigma
% attention figure
figurename = 'attention_fMRINoiseScatterPlot';
h = smartfig(figurename,'reuse');
set(h,'Name',figurename);


% take the mean across observers and visual areas:
meanOriginalSTD = squeeze(mean(originalSTD,1));

% % 95% confidence intervals of the STD
meanSTD = squeeze(mean(bootSTD,1));
CI_meanSTD = zeros(size(meanSTD,1),size(meanSTD,2));
for v = 1:size(meanSTD,1)
 for icond = 1:size(meanSTD,2)
  this_meanSTD = squeeze(meanSTD(v,icond,:));
  CI_meanSTD(v,icond,1:2) = quantile(this_meanSTD,[.05 .95]);
  for ci = 1:2
   if ci == 1
    CI_meanSTD(v,icond,ci) = meanOriginalSTD(v,icond) - CI_meanSTD(v,icond,ci);
   else
    CI_meanSTD(v,icond,ci) = CI_meanSTD(v,icond,ci) - meanOriginalSTD(v,icond);
   end
   end
 end
end

% plot mean across observers
plot([.5,.5;2.5,2.5],[.5,.5;2.5,2.5],'k-','Color',plotInfo.XYColor)
hold on;
for iVis = 1:size(originalSTD,2) % visual areas
 myerrorbar(meanOriginalSTD(iVis,3),meanOriginalSTD(iVis,1), ...
  'yError', CI_meanSTD(v,3,:), ...
  'xError', CI_meanSTD(v,1,:), ...
  'Symbol',plotInfo.thisSymbol{1}, ...
  'LineWidth',plotInfo.LineWidth, ...
  'Color',plotInfo.thisColor{1}{iVis}, ...
  'MarkerFaceColor',plotInfo.MarkerFaceColor{1}{iVis}, ...
  'MarkerEdgeColor','w', ...
  'MarkerSize', plotInfo.MarkerSize(1));
end

% plot each observer
loglog([.5,.5;4,4],[.5,.5;4,4],'k-','Color',plotInfo.XYColor)
hold on;
for iObserver = size(originalSTD,1):-1:1
 for iVis = 1:size(originalSTD,2) % visual areas
  myerrobar(originalSTD(iObserver,iVis,3),originalSTD(iObserver,iVis,1), ...
   plotInfo.thisSymbol{iObserver}, ...
   'LineWidth',plotInfo.LineWidth, ...
   'Color',plotInfo.thisColor{1}{iVis}, ...
   'MarkerFaceColor',plotInfo.MarkerFaceColor{1}{iVis}, ...
   'MarkerEdgeColor','w', ...
   'MarkerSize', plotInfo.MarkerSize(iObserver));
 end
end

title(sprintf('fMRI noise'))
% xlabel('Distributed','FontSize',plotInfo.Fsize);
% ylabel(sprintf('Attended'),'FontSize',plotInfo.Fsize);

set(gca,...
 'FontName','Helvetica','FontSize',plotInfo.Fsize, ...
 'PlotBoxAspectRatio',plotInfo.PlotBoxAspectRatio, ...
 'XLim', [.5 4],'YLim', [.5 4],...
 'LineWidth',plotInfo.LineWidth,'TickLength',plotInfo.TickLength, ...
 'XTick', [.5 1 2 4], 'XTickLabel', [.5 1 2 4] ,...
 'YTick', [.5 1 2 4], 'YTickLabel', [.5 1 2 4] ,...
 'XColor',plotInfo.XYColor,'YColor',plotInfo.XYColor,'Box','off');

myaxisScatterLog;

keyboard
if saveFig
 set(h,'PaperPosition',[.25 .25 8 10.5]);
 set(h,'PaperOrientation','Portrait');
 
 disp(sprintf('Saving figure %s',figurename));
 defaultDataDir = '/Users/frakkopesto/data/riken/crfaa/fmridata/crfaa/';
 figDir = 'fig_scatter_noiseTest';
 figDir = savefig(h,'defaultDataFolder',defaultDataDir,'figDir',figDir,'figName',figurename);
end



%%%%%%%%%%%%%%%%%%%%%
% plotEachBootStrap %
%%%%%%%%%%%%%%%%%%%%%
function plotEachBootStrap(bootSTD,originalSTD,observer,visualAreas,numConditions,numBootStraps,normalizeByTrialNum)
conditions = {'Focal cue, target','Focal cue, non-target','Distributed cue, target','Distributed cue, non-target'};

for o = 1:length(observer)      % do the whole thing for each observer
 for v = 1:length(visualAreas) % v1-v4
  for icond = 1:numConditions   
   h = smartfig(sprintf('BootstraphistOBS%sVisArea%s',observer{o},visualAreas{v}),'reuse');
   subplot(4,1,icond),
   thisBT = squeeze(bootSTD(o,v,icond,:));
   hist(thisBT);
   legend(conditions{icond});
   
   % set the axis limits properly
   if strcmp(visualAreas{v},'v4')
    if strcmp(observer{o},'fp')
     axis([0 4 0 numBootStraps/3]);
    else
     axis([0 2.5 0 numBootStraps/3]);
    end
   else
    axis([0 1.6 0 numBootStraps/3]);
   end
   
   vline(originalSTD(o,v,icond),'g-');hold on
   vline(mean(thisBT(:)),'r--');
  end
 end
end


%%%%%%%%%%%%
% loadData %
%%%%%%%%%%%%
function d = loadData(observer,visualAreas,adaptation,analysisType)

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
 [filename vars2load] = makeFileName(observer,[],visualAreas, analysisType);
 [data] = load(sprintf('%s',filename),vars2load{1}, vars2load{2});
 d = data.d1;
 clear data;
 
else % individual observers
 if ~(strcmpi('fm',observer) && adaptation == 100)
  [filename vars2load] = makeFileName(observer,adaptation,visualAreas,analysisType);
  data = load(sprintf('%s',filename),vars2load);
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % now get the STD estimates from the analysis loaded %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  switch analysisType
   case {'standard' '4conditions' '4c' '4' 's' 'CorrectIncorrect' 'corrincorr' 'CorrIncorr' 'ci' '1st2ndInterval' '1st2nd' '12'}
    d.groupSTD = data.d1.groupSTD;
        
   otherwise
    keyboard
  end
 end
end


%%%%%%%%%%%%%%%%
% makeFileName %
%%%%%%%%%%%%%%%%
function [filename vars2load] = makeFileName(observer,adaptation,visualArea,analysisType)
% NB whereas the main analysis files are saved into the special folder
% <data_used_files> so that changes to that folder do not screw up the
% files...
% the files for the other analyses are instad in each folder/session. so
% for construtcting the filename we use getsession.

% get the right session:
session = getsession(observer,adaptation);

if strcmpi(observer,'avrg')
keyboard % this is not operative here.

else
 % choose the type of analysis files to load
 switch analysisType
  case {'standard' '4conditions' '4c' '4' 's'} % this is the original analysis for the main results
   switch visualArea
    case {'v1' 'v2' 'v3' 'v4'}
     day = '15-Feb-2009';
     filename = deblank(sprintf('%s/%sNoiseAnalysis%s.mat',session,day,lower(visualArea)));

    case {'05' '06' '07' '075'}
     day = '05-Jul-2009';
     filename = deblank(sprintf('%s/%s_v1_subTSnoiseSTD_combinedR2%s.mat',session,day,lower(visualArea)));

    otherwise
     keyboard
   end
   
  case {'CorrectIncorrect' 'corrincorr' 'CorrIncorr' 'ci'} % this is the correct incorrect analysis
   filename = sprintf('%s/28-Jul-2009_%s_correctIncorrectAnalysis_STD_%s.mat',session,visualArea,visualArea);

  case {'1st2ndInterval' '1st2nd' '12'} % This is the first and second interval analysis
    filename = sprintf('%s/28-Jul-2009_%s_firstSecondIntervalAnalysis_STD_%s.mat',session,visualArea,visualArea);
  otherwise
   keyboard
 end
 vars2load = 'd1';
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


if ~isdir(hd)
  hd = '/Volumes/frakkopesto_HomeDir/data/riken/crfaa';
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

