function sigma = crfaaMaxPool(varargin)
% sigma = crfaaMaxPool.m
%
%        $Id:$ 
%      usage: crfaaMaxPool()
%         by: justin gardner, updated by franco pestilli
%       date: 2009/10/28
%    purpose: decision model using softmax/exponential sum pooling.
%
% The following re the poolingRules that are implemented:
% (poolingRule = 1) max, or Winner-takes-all, choose the stimulus that gives max response
% (poolingRule = 2) average, responses across locations
% (poolingRule = 3) softmax, weight differently each location, it goes from a perfect
%                   average to a perfect WTA
% (poolingRule = 4) targetLoc, always attend to the correct location (sensitivity model)
% (poolingRule = 5) targetLoc/average, FIX FIX FIX not sure 
% (poolingRule = 6) Exponent sum equivalent to softmax but with a biologically plausible
%                   implementation (sum of exponent, rectification, due to cell-response 
%                   normalization)
% (poolingRule = 7-9) Same as 6 - i.e. uses exponentSum - but with observer attending to 
%                   2, 3 or 4 targets respectively
% (poolingRule = 10) Same as 9, (exponent-sum with observer distributing attention across
%                   all 4 targets, but all distractors have a lower contrast)
% (poolingRule = 11) Exponent-sum non-strategy rule. The one used in the paper.
%
% The following are the type of crf that can be used to compute the stimulus response:
% (crftype = 0) raw for attended and unattended)
% (crftype = 1) tvc2crf fit for attended, raw for unattended)
% (crftype = 2) tvc2crf fit for attended, offset adjusted tvc2crf fit for unattended);
% (crftype = 3) tvc2crf average fit for attended/distributed, 
%               offset adjusted to fit independently attended and unattended curves)
% (crftype = 4) Use independently fit functions for attended and unattended
% (crftype = 5) tvc2crf fit for attended - unattended is offset adjusted to fit
%               tvc2crf fit for distributed -- offset adjusted to fit distributed target/non-target.
%               Should be used in paper
%
% Other options:
% (makeDistPlots = 0) shows plot of distributions of responses by quadrant
% (doSigma       = 0) shows scatter-plot of reduction in variability of the distribution 
%                     of resposes
% (doSearchTau   = 1) searches for optimal tau for each observer and
%                     visual area, makes histograms, NB it takes few days to run
% (doTvCFits     = 0) fits the TvC with the three models (average, WTA or
%                     optimal), generates graphs of each fit, NB it takes few minutes to run
% (recalc        = 0) Recomputes from scratch, rather than loading data
% (observerNum   = 1) Choses which observer data to use, {'Average' 'fp' 'jg' 'fm'} 
%
% (visuaArea     = 1) Select which visual are crf to use {'V1' 'V2' 'V3' 'V4'}
% (dispFig       = 0) Allows to display figure while computing
% (snrImprovementWithAttention = 0) Sets how much sensory SNR improvement there
%                     is with attention. For the Mitchell et al. paper this is 0.39 , for
%                     Cohen it should be 0.54. Defaults to 0, so that there is none.
% (crossValidate = 1) Fit with odd data points, test on even data points
%
% < end of help >

% display help with no arguments
if nargin == 0,help crfaaMaxPool;return,end

%%%%%%%%%%%%%%%%%%%%% set default variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
crftype       = [];
poolingRule   = [];
crffit        = [];
decisionModel = [];
getWeights    = []; 
makeDistPlots = [];
doSigma       = [];
doSearchTau   = [];
doTvCFits     = [];
doStdFit      = [];
recalc        = [];
observerNum   = [];
visualArea    = [];
dispFig       = [];
dataDir       = [];
bootstrapNum  = [];
snrImprovementWithAttention = [];
global crossValidate;

getArgs(varargin,{ ...
         'decisionModel','decisionModelData_Dnt.mat',      ...
         'crffit', '03-Aug-2009nestsedTests_standard.mat', ...
         'dataDir',[],...
         'poolingRule', 11,   ...
         'crftype', 5,       ...
         'getWeights', 0,    ...
         'makeDistPlots', 0, ...
         'doSigma', 0,       ... % this is VERY slow 
         'doSearchTau', 0,   ... % this is VERY slow
         'doTvCFits', 0,     ... % this is slow 
         'doStdFit', 0,      ...
         'recalc', 0,        ... % this is VERY slow
         'observerNum', 1,   ...
         'visualArea', 1,    ...
         'snrImprovementWithAttention',0, ...
         'crossValidate',1,...
         'bootstrapNum',[],...
         'dispFig', 0});
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% check on which computer this code is being called and
% set corrects paths to data folders and files:
if isempty(dataDir)
  dataFilenames = setDataDefault(decisionModel,crffit);
else
  dataFilenames = setDataDir(dataDir,bootstrapNum);
end

d = getData(dataFilenames,observerNum,visualArea,dispFig);

% recompute the weights
if (getWeights)
 tau = 0.23;
 [weightsAttended weightsDistributed] = getWeightsForTau(d,poolingRule,tau);
end

% make a plot of each distribution of responses by quadrants
if (makeDistPlots)
% poolingRule = 3; tau = 5;sigma = 0.0205;crftype = 3;
  tau = 13.38;sigma = 0.04;
  pedNum = 2;
  makeDistributionPlot(d,poolingRule,tau,sigma,crftype,pedNum,snrImprovementWithAttention)
  keyboard
end

% recompute sigma
if (doSigma)
  sigma = doSigmaPlot(dataFilenames,poolingRule,crftype,recalc,snrImprovementWithAttention,dispFig);
end

% search for the optimal tau
if (doSearchTau)
 doSearchTauPlot(dataFilenames,poolingRule,crftype,recalc,snrImprovementWithAttention,dispFig);
end

% display a plot of the fit to the TvC predicted by the model
if (doTvCFits)
 doTvCFitsPlot(dataFilenames,poolingRule,crftype,recalc,dispFig,observerNum,visualArea,snrImprovementWithAttention);
end

if (doStdFit)
  doStdFitPlot(dataFilenames);
end

%%%%%%%%%%%%%%%%%%%%%%
%    doStdFitPlot    %
%%%%%%%%%%%%%%%%%%%%%%
function doStdFitPlot(dataFilenames)

% load file
filename = sprintf('%s.mat',fullfile(dataFilenames.path,dataFilenames.tbyt));
if ~isfile(filename)
  disp(sprintf('(crfaaMaxPool:doStdFitPlot) Could not find file %s',filename));
  return
end
load(filename);

% get some sizes
nVisualAreas = size(d,2);
nObservers = size(d,1);
nConditions = 4;

clear a;
for v = 1:nVisualAreas
  % init amplitudes
  for i = 1:nConditions,a{v,nConditions} = [];end
  % collect amplitudes across observers
  for o = 1:nObservers
    for i = 1:nConditions
      a{v,i} = [a{v,i} d{o,v}.groupSTD.amplitude{i}];
    end
  end
end

focalCueTargetOffset = [0.29 0.39 0.52 0.51];
smartfig('crfaaMaxPoolVar','reuse');
clf;
for v = 1:nVisualAreas
  attended = a{v,2} + focalCueTargetOffset(v);
  unattended = a{v,1};
  distributedNonTarget = a{v,3};
  distributedTarget = a{v,4};
  
  n = 5000000;pAttended = 0.25;
  mix = [randsample(attended,round(n*pAttended),true) randsample(unattended,round(n*(1-pAttended)),true)];

  stds(v,:) = [std(attended),std(unattended),std(distributedNonTarget),std(distributedTarget),std(mix)];
  disp(sprintf('A: %0.3f U: %0.3f DNT: %0.3f DT: %0.3f MIX: %0.3f',std(attended),std(unattended),std(distributedNonTarget),std(distributedTarget),std(mix)));
end
    
mybar(stds,'groupLabels',{'V1','V2','V3','V4'},'withinGroupLabels',{'focal-cue target','focal-cue non-target','distributed-cue non-target','distributed-cue target','mix'},'yAxisMin=0',sprintf('yAxisMax=%0.2f',2.5));
keyboard



% Subsidiary functions that actually do the work.
%%%%%%%%%%%%%%%%%%%%%
%   doTvCFitsPlot   %
%%%%%%%%%%%%%%%%%%%%%
function doTvCFitsPlot(dataFilenames,poolingRule,crftype,recalc,dispFig,observerNum,visualArea,snrImprovementWithAttention)

% dispFig = 0;
% observerNum = 1;
% visualArea = 1;
d = getData(dataFilenames,observerNum,visualArea,dispFig);

if 0
  % both targetLoc
  thisPoolingRule = 4;
  [sigmaDistributed sigmaAttended] = CRF2TvCfitSigma(d,thisPoolingRule,nan,crftype,[],dataFilenames);
  sigma = sigmaAttended; %0.1
  dfit = CRF2TvCfitDeltaC(d,thisPoolingRule,nan,sigma,crftype,dataFilenames);
  dispModelFitToTvC(d,dfit); % new call: dispModelFitToTvC(d,dout,dataFolder,visArea)

  % both targetLoc -- with snrImprovement as found in Mitchell et al.
  thisPoolingRule = 4;
  [sigmaDistributed sigmaAttended] = CRF2TvCfitSigma(d,thisPoolingRule,nan,crftype,[],dataFilenames);
  sigma(2) = sigmaAttended; %0.1
  sigma(1) = sigmaAttended*(1+snrImprovementWithAttention);
  dfit = CRF2TvCfitDeltaC(d,thisPoolingRule,nan,sigma,crftype,dataFilenames);
  dispModelFitToTvC(d,dfit); % new call; dispModelFitToTvC(d,dout,dataFolder,visArea)

  % targetLoc versus averaging w/snrImprovement
  thisPoolingRule = 5;
  [sigmaDistributed sigmaAttended] = CRF2TvCfitSigma(d,thisPoolingRule,nan,crftype,[],dataFilenames);
  sigma(2) = sigmaAttended; %0.1
  sigma(1) = sigmaAttended*(1+snrImprovementWithAttention);
  dfit = CRF2TvCfitDeltaC(d,thisPoolingRule,nan,sigma,crftype,dataFilenames);
  dispModelFitToTvC(d,dfit); % new call: dispModelFitToTvC(d,dout,dataFolder,visArea)
end

% optimal rule
[tau sigma sigmad dfit] = fitTau(d,poolingRule,crftype,dataFilenames,recalc,snrImprovementWithAttention);

dfit = CRF2TvCfitDeltaC(d,poolingRule,tau,sigma,crftype,dataFilenames,recalc);
%dispModelFitToTvC(d,dfit); % new call: 
dispModelFitToTvC(d,dfit,dataFilenames,1)
if 0
  makeDistributionPlot(d,poolingRule,tau,sigma,crftype,5)
end
plotCRF(d,crftype)

%%%%%%%%%%%%%%%%%
%    plotCRF    %
%%%%%%%%%%%%%%%%%
function plotCRF(d,crftype)

smartfig('plotCRF');
CRF = getCRFFromD(d,crftype);
plot(d.unattendedContrast,d.unattendedResponse,'go');hold on
plot(d.attendedContrast,d.attendedResponse,'ro');
plot(d.attendedContrast,d.distributed.CRFtarget,'bo');
plot(d.attendedContrast,d.distributed.CRFnontarget,'mo');
plot(CRF.unattendedContrast,CRF.unattendedResponse,'g-');
plot(CRF.attendedContrast,CRF.attendedResponse,'r-');
plot(CRF.distributedContrast,CRF.distributedResponse,'b-');
plot(CRF.distributedContrast,CRF.distributedNontargetResponse,'m-');
xlabel('contrast');
ylabel('response');
legend('focal non-target','focal target','distributed target','distributed non-taret');

%%%%%%%%%%%%%%%%%%%%%%%
%   doSearchTauPlot   %
%%%%%%%%%%%%%%%%%%%%%%%
function doSearchTauPlot(dataFilenames,poolingRule,crftype,recalc,snrImprovementWithAttention,dispFig)

dispFig = 0;
filename = fullfile(dataFilenames.saveDataDir,fixBadChars(sprintf('optimalTauSigma_%s_crftype%i_%s',poolingRuleNum2Str(poolingRule),crftype,mynum2str(snrImprovementWithAttention,'doFixBadChars=1'))));
if isfile(sprintf('%s.mat',filename))
  disp(sprintf('(doSearchTauPlot) Loading %s',filename));
  eval(sprintf('load %s',filename));
else
  tau = 0;
end

% FIX FIX FIX start at 1 here
observerNumStart = 1;
observerNumEnd = 4;
% for bootstrap, we only do the observers individually
if ~isempty(dataFilenames.bootstrapNum)
  observerNumEnd = 3;
end
  
visualAreaStart = 1;

% now for each observer /visual area calculate the optimal tau
% that matches the sigmas
for observerNum = observerNumStart:observerNumEnd;
  for visualArea = visualAreaStart:4;
    d = getData(dataFilenames,observerNum,visualArea,dispFig);
    disp(sprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'));
    disp(sprintf('Fit tau for observer %i visualArea %s',observerNum,d.visualArea));
    disp(sprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'));
    [tau(observerNum,visualArea) sigma(observerNum,visualArea) sigmad(observerNum,visualArea)] = fitTau(d,poolingRule,crftype,dataFilenames,recalc,snrImprovementWithAttention);
    eval(sprintf('save %s tau sigma sigmad',filename));
  end
end

% get a d structure for calculate what weights look like in plot below
observerNum = 1;visualArea = 1;
d = getData(dataFilenames,observerNum,visualArea,dispFig);

% make plot
tau2 = tau(observerNumStart:observerNumEnd,:);
if dispFig
  smartfig('crfaaMaxPool:tau');
  if poolingRule == 6
    bins = 1:0.5:10;
  elseif poolingRule == 11
    bins = floor(min(tau2(:))):(ceil(max(tau2(:)))-floor(min(tau2(:))))/5:ceil(max(tau2(:)));
  else
    bins = 0.05:0.1:1;
  end
  myhist(tau2(:),bins);
  plotmean(tau2(:));%,'k',1);
  ylabel('Num observers/visual areas');
  if poolingRule == 6
    xlabel('Exponent');
    xaxis(0,7);
    yaxis(0,7);
    axis square
    [weightsAttended] = getWeightsForTau(d,poolingRule,1);
    plotWeights(weightsAttended,[0.5 1.5],[5.75 6.75],0);
    [weightsAttended] = getWeightsForTau(d,poolingRule,3.58);
    plotWeights(weightsAttended,[3 4],[5.75 6.75],0);
    [weightsAttended] = getWeightsForTau(d,poolingRule,2.25);
    plotWeights(weightsAttended,[1.75 2.75],[5.75 6.75],0);
    [weightsAttended] = getWeightsForTau(d,poolingRule,4.75);
    plotWeights(weightsAttended,[4.25 5.25],[5.75 6.75],0);
    [weightsAttended] = getWeightsForTau(d,poolingRule,6);
    plotWeights(weightsAttended,[5.5 6.5],[5.75 6.75],0);
  else
    xlabel('k');
    xaxis(min(0,min(bins)),max(bins));
    yaxis(0,7.5);
    [weightsAttended] = getWeightsForTau(d,poolingRule,0.1);
    plotWeights(weightsAttended,[0.05 0.15],[6 7],0);
    [weightsAttended] = getWeightsForTau(d,poolingRule,0.351);
    plotWeights(weightsAttended,[0.3 0.4],[6 7],0);
    [weightsAttended] = getWeightsForTau(d,poolingRule,0.6);
    plotWeights(weightsAttended,[0.55 0.65],[6 7],0);
    [weightsAttended] = getWeightsForTau(d,poolingRule,1);
    plotWeights(weightsAttended,[0.95 1.05],[6 7],0);
  end
end
  

%%%%%%%%%%%%%%%%%%%
%   doSigmaPlot   %
%%%%%%%%%%%%%%%%%%%
function sigma = doSigmaPlot(dataFilenames,poolingRule,crftype,recalc,snrImprovementWithAttention,dispFig)

if ieNotDefined('dispFig'),dispFig = 0;end
% show only the average of all three observers
averageObserverOnly = 0;

% use a separate k for each observet
kByObserver = 1;

% for the bar graphs, average over these observers
averageOver = 2:4;

% for averageObserver, we don't have observers to average over
if averageObserverOnly
  averageOver = 1;
end

disp(sprintf('(crfaaMaxPool) Calculating sigmas'));
filename = fullfile(dataFilenames.saveDataDir,fixBadChars(sprintf('optimalTauSigma_%s_crftype%i_%s',poolingRuleNum2Str(poolingRule),crftype,mynum2str(snrImprovementWithAttention,'doFixBadChars=1'))));
if ~isfile(sprintf('%s.mat',filename))
  disp(sprintf('(crfaaMaxPool:doSigmaPlot) Precompute file %s not found. Run doSearchTau with these parameters first',filename));
  keyboard
else
  eval(sprintf('load %s',filename));
end

originalTau = tau;

if averageObserverOnly
  tau = mean(tau(1,:));
  observers = 1;
else
  if kByObserver
    tau = mean(tau,2);
  else
    tau = mean(tau(:));
  end
  if isempty(dataFilenames.bootstrapNum )
    observers = 1:4;
  else
    observers = 1:3;
    if length(averageOver) == 3
      averageOver = 1:3;
    end
  end
end

% get the fit parameters when tau is individually fit
for iObserver = observers
  for iArea = 1:4
    [thisSigmaCrossVal thisDfitCrossVal] = getCrossValSigmaAcrossObserversAndAreas(dataFilenames,poolingRule,originalTau(iObserver,iArea),crftype,recalc,iObserver,iArea,snrImprovementWithAttention);
    [thisSigmaCrossValSensitivity thisDfitCrossValSensitivity] = getCrossValSigmaAcrossObserversAndAreas(dataFilenames,4,originalTau(iObserver,iArea),crftype,recalc,iObserver,iArea,snrImprovementWithAttention);
    % grab the right one
    sigmaCrossValFitTau(iObserver,iArea) = thisSigmaCrossVal(iObserver,iArea);
    dfitCrossValFitTau(iObserver,iArea) = thisDfitCrossVal(iObserver,iArea);
    sigmaCrossValSensitivityFitTau(iObserver,iArea) = thisSigmaCrossValSensitivity(iObserver,iArea);
    dfitCrossValSensitivityFitTau(iObserver,iArea) = thisDfitCrossValSensitivity(iObserver,iArea);
  end
end

% display statistics for each
for iObserver = 1:size(dfitCrossValFitTau,1)
  for iArea = 1:size(dfitCrossValFitTau,2)
    % get the fits
    selectionFit = dfitCrossValFitTau(iObserver,iArea);
    sensitivityFit = dfitCrossValSensitivityFitTau(iObserver,iArea);
    % grab some info
    k(iObserver,iArea) = selectionFit.tau;
    r2selection(iObserver,iArea) = selectionFit.crossValidate.r2;
    r2sensitivity(iObserver,iArea) = sensitivityFit.crossValidate.r2;
    AICdiff(iObserver,iArea) = selectionFit.overallAIC-sensitivityFit.overallAIC;
    % display info
    disp(sprintf('Observer: %i Area: %i Selection (Cross-val r2: %f k: %f) Sensitivity (Cross-val r2: %f) AIC diff: %f',iObserver,iArea,r2selection(iObserver,iArea),k(iObserver,iArea),r2sensitivity(iObserver,iArea),AICdiff(iObserver,iArea)));
    
  end
end

if dispFig
  for iArea = 1:4
    d = getData(dataFilenames,1,iArea,0);
    disp(sprintf('Plotting TvC FIT (V%i)',iArea));
    
    dispModelFitToTvC(d,dfitCrossValSensitivityFitTau(1,iArea),dataFilenames,iArea);
    dispModelFitToTvC(d,dfitCrossValFitTau(1,iArea),dataFilenames,iArea);
  end
end

% display average stats
for iArea = 1:4
  disp(sprintf('V%i: selection(r2=%f) sensitivity(r2=%f) AIC difference: %f',iArea,mean(r2selection(observers,iArea)),mean(r2sensitivity(observers,iArea)),mean(AICdiff(observers,iArea))));
end
disp(sprintf('Mean k for average: %f', mean(originalTau(1,:))));
k = originalTau(observers,:);
disp(sprintf('Mean k for subjects: %f', mean(k(:))));

% compute the sigma
if kByObserver
  % use a different tau for each observer
  for iObserver = observers
    [thisDistributedSigma thisAttendedSigma] = getSigmaAcrossObserversAndAreas(dataFilenames,poolingRule,tau(iObserver),crftype,recalc,iObserver);
    [thisSigmaCrossVal thisDfitCrossVal] = getCrossValSigmaAcrossObserversAndAreas(dataFilenames,poolingRule,tau(iObserver),crftype,recalc,iObserver);
    [thisSigmaCrossValSensitivity thisDfitCrossValSensitivity] = getCrossValSigmaAcrossObserversAndAreas(dataFilenames,4,tau(iObserver),crftype,recalc,iObserver);
    % grab the right one
    sigmaDistributedSoftmax(iObserver,:) = thisDistributedSigma(iObserver,:);
    sigmaAttendedSoftmax(iObserver,:) = thisAttendedSigma(iObserver,:);
    sigmaCrossVal(iObserver,:) = thisSigmaCrossVal(iObserver,:);
    dfitCrossVal(iObserver,:) = thisDfitCrossVal(iObserver,:);
    sigmaCrossValSensitivity(iObserver,:) = thisSigmaCrossValSensitivity(iObserver,:);
    dfitCrossValSensitivity(iObserver,:) = thisDfitCrossValSensitivity(iObserver,:);
  end
else
  % use same tau for each observer
  [sigmaDistributedSoftmax sigmaAttendedSoftmax] = getSigmaAcrossObserversAndAreas(dataFilenames,poolingRule,tau,crftype,recalc,observers);
  [sigmaCrossVal dfitCrossVal] = getCrossValSigmaAcrossObserversAndAreas(dataFilenames,poolingRule,tau,crftype,recalc,observers);
end

poolingRule = 4;tau = nan;
%recalc = 1;
[sigmaDistributedTargetLoc sigmaAttendedTargetLoc] = getSigmaAcrossObserversAndAreas(dataFilenames,poolingRule,tau,crftype,recalc,observers);
%poolingRule = 5;tau = nan;
%[sigmaDistributedTargetLocAverage sigmaAttendedTargetLocAverage] = getSigmaAcrossObserversAndAreas(dataFilenames,poolingRule,tau,crftype,recalc);

% save the sigma obtained for each observer and visual areas
% in a variable to output.
% -- sigmaSoftmax is the selection model sigma, 
%    one for focal one for distributed
% -- sigmaTargetLoc is the sensitivity model sigma, 
%    one for focal on for distributed
% by Franco doing the botostrap tests.
sigma.sensitivity = nan(3, size(sigmaDistributedTargetLoc,1), size(sigmaDistributedTargetLoc,2));
sigma.selection = sigma.sensitivity;
sigma.sensitivity(1,:,:) = sigmaDistributedTargetLoc;
sigma.sensitivity(3,:,:) = sigmaAttendedTargetLoc;

sigma.selection(1,:,:) = sigmaDistributedSoftmax;
sigma.selection(3,:,:) = sigmaAttendedSoftmax;

% make scatter plot
if dispFig
  figName = 'scatter';
  h = smartfig(figName,'reuse');
  loglog(sigmaDistributedTargetLoc(1,:),sigmaAttendedTargetLoc(1,:),'ko','MarkerFaceColor','k','MarkerSize',10);
  hold on
  loglog(sigmaDistributedSoftmax(1,:),sigmaAttendedSoftmax(1,:),'ro','MarkerFaceColor','r','MarkerSize',10);
  if ~averageObserverOnly
    loglog(sigmaDistributedTargetLoc(observers(1),:),sigmaAttendedTargetLoc(observers(1),:),'kx','MarkerFaceColor','k');
    loglog(sigmaDistributedTargetLoc(observers(2),:),sigmaAttendedTargetLoc(observers(2),:),'kd','MarkerFaceColor','k');
    loglog(sigmaDistributedTargetLoc(observers(3),:),sigmaAttendedTargetLoc(observers(3),:),'kp','MarkerFaceColor','k');
  end
  %loglog(sigmaDistributedTargetLocAverage,sigmaAttendedTargetLocAverage,'go','MarkerFaceColor','g');
  if ~averageObserverOnly
    loglog(sigmaDistributedSoftmax(observers(1),:),sigmaAttendedSoftmax(observers(1),:),'rx','MarkerFaceColor','r');
    loglog(sigmaDistributedSoftmax(observers(2),:),sigmaAttendedSoftmax(observers(2),:),'rd','MarkerFaceColor','r');
    loglog(sigmaDistributedSoftmax(observers(3),:),sigmaAttendedSoftmax(observers(3),:),'rp','MarkerFaceColor','r');
  end
  if 0
    clf
    plot(sigmaDistributedTargetLoc(1,:),sigmaAttendedTargetLoc(1,:),'ko','MarkerFaceColor','k','MarkerSize',10);
    hold on
    plot(sigmaDistributedTargetLoc(2,:),sigmaAttendedTargetLoc(2,:),'kx','MarkerFaceColor','k');
    plot(sigmaDistributedTargetLoc(3,:),sigmaAttendedTargetLoc(3,:),'kd','MarkerFaceColor','k');
    plot(sigmaDistributedTargetLoc(4,:),sigmaAttendedTargetLoc(4,:),'kp','MarkerFaceColor','k');

    plot(sigmaDistributedSoftmax(1,:),sigmaAttendedSoftmax(1,:),'rx','MarkerFaceColor','r');
    plot(sigmaDistributedSoftmax(2,:),sigmaAttendedSoftmax(2,:),'ro','MarkerFaceColor','r');
    plot(sigmaDistributedSoftmax(3,:),sigmaAttendedSoftmax(3,:),'rd','MarkerFaceColor','r');
    plot(sigmaDistributedSoftmax(4,:),sigmaAttendedSoftmax(4,:),'rp','MarkerFaceColor','r');
  end

  xaxis(0.01,0.15);
  yaxis(0.01,0.15);
  dline
  axis('square')

  xlabel('sigma distributed (% signal change)');
  ylabel('sigma attended (% signal change)');
  legend('Sensitivity model','Selection model');

  savefigure = 1;
  if savefigure
    savefig(h,'verbose',1,'figName',figName,'figDir','fig_sigmaratio','defaultDataFolder',dataFilenames.saveDataDir);
  end

  % make the bar graph.
  figName = 'ratio_plot';
  h = smartfig(figName,'reuse');

  % get the sensitivity and selection sigma ratios we need
  sensitivity = sigmaDistributedTargetLoc(averageOver,:)./sigmaAttendedTargetLoc(averageOver,:);
  selection = sigmaDistributedSoftmax(averageOver,:)./sigmaAttendedSoftmax(averageOver,:);

  % some things we need for doing small n correction
  n = length(averageOver);
  smallNcorrection = sqrt(2/(n-1))*gamma(n/2)/gamma((n-1)/2);

  % now calculate mean and standard error
  ratioMean = [ mean(sensitivity,1) ; mean(selection,1)];
  %ratioSTE = [std(selection,1) ; std(sensitivity,1)]/sqrt(length(averageOver));
  % calculate ratio STE using small n corection
  ratioSTE = [std(sensitivity)/smallNcorrection ;  std(selection)/smallNcorrection]/sqrt(n);

  % display the bar graph
  mybar(ratioMean,'yError',ratioSTE,'yMin',ratioMean-ratioSTE,'groupLabels',{'Sensitivity', 'Selection'}, ...
	'withinGroupLabels',{'V1','V2','V3','V4'}, ...
	'yAxisMin=0','yAxisMax=7', ...
	'yLabelText=Ratio of sigmaDistributed to sigmaAttended', ...
	'xLabelText=Model','dispValues=1','hline=1','hLineStyle=--');
            
  savefigure = 1;
  if savefigure
    savefig(h,'verbose',1,'figName',figName,'figDir','fig_sigmaratio','defaultDataFolder',dataFilenames.saveDataDir);
  end

  if 0
    poolingRule = 3;
    [tau sigma] = fitTau(d,poolingRule,crftype,dataFilenames,0,snrImprovementWithAttention);
    [sigmaDistributed sigmaAttended] = getSigmaAcrossObserversAndAreas(dataFilenames,poolingRule,tau,crftype,recalc);
    poolingRule = 2;tau = nan;
    [sigmaDistributed sigmaAttended] = getSigmaAcrossObserversAndAreas(dataFilenames,poolingRule,tau,crftype,recalc);
    poolingRule = 3;tau = 1000;
    [sigmaDistributed sigmaAttended] = getSigmaAcrossObserversAndAreas(dataFilenames,poolingRule,tau,crftype,recalc);
  end
end

%%%%%%%%%%%%%%
%   fitTau   %
%%%%%%%%%%%%%%
function [tau sigmaa sigmad dfit] = fitTau(d,poolingRule,crftype,dataFilenames,recalc,snrImprovementWithAttention)

dispFit = 0;
sigma=[];
% get filename for saving.
filename = fullfile(dataFilenames.saveDataDir,fixBadChars(sprintf('bestTauFor%s_%i_o%i_%s_%s',poolingRuleNum2Str(poolingRule,nan),crftype,d.observerNum,d.visualArea,mynum2str(snrImprovementWithAttention,'doFixBadChars=1'))));

% check to see whether we should just load precomputed
if isfile(sprintf('%s.mat',filename)) && ~recalc
  disp(sprintf('(crfaaMaxPool) LOADING precomputed tau for %s from %s',poolingRuleNum2Str(poolingRule,nan),filename));
  load(filename);
  % fix for when sigmas were calculated incorrectly
  if sigma(1) > sigma(2)*(1+snrImprovementWithAttention)
    newsigma = sigma;
    newsigma(1) = newsigma(2)*(1+snrImprovementWithAttention);
    disp(sprintf('(crfaaMaxPool:fitTau) FIXING sigma from %s to %s',mynum2str(sigma),mynum2str(newsigma)));
    sigma = newsigma;
    dfit = CRF2TvCfitDeltaC(d,poolingRule,tau,sigma,crftype,dataFilenames,1);
    save(filename);
  end
  % check whether we need to do fit
  if ~exist('dfit','var')
    dfit = CRF2TvCfitDeltaC(d,poolingRule,tau,sigma,crftype,dataFilenames,1);
    save(filename);
  end
  if ~isfield(dfit,'modelName'),
    dfit.modelName = sprintf('%s crftype=%i',poolingRuleNum2Str(poolingRule,tau),crftype);
  end
  if ieNotDefined('sigmaa')
    sigmaa = sigma(2);
  end
  return
else
  disp(sprintf('(crfaaMaxPool:fitTau) Searching for tau using %s',poolingRuleNum2Str(poolingRule,nan)));
end

% look for sigmas bi bisection across tau
% start with extreme values of tau
if poolingRule == 3
  tau = [0.1 5];
else
  tau = [5 1];
end
tau = [20 1];

% compute sigma and p for extreme values
[sigmaDistributed(1) sigmaAttended(1)] = CRF2TvCfitSigma(d,poolingRule,tau(1),crftype,0,dataFilenames,recalc);
% display fit if called for
if dispFit
%  dfit = CRF2TvCfitDeltaC(d,poolingRule,tau(1),[sigmaDistributed(1) sigmaAttended(1)],crftype,dataFilenames,recalc);
  dfit = CRF2TvCfitDeltaC(d,poolingRule,tau(1),[sigmaAttended(1) sigmaAttended(1)],crftype,dataFilenames,recalc);
%  dfit = CRF2TvCfitDeltaC(d,poolingRule,tau(1),[10 sigmaAttended(1)],crftype,dataFilenames,recalc);
  disp(sprintf('**** TAU: %f sigmaAttended: %f sigmaDistributed: %f',tau(1),sigmaAttended(1),sigmaDistributed(1)));
  dispModelFitToTvC(d,dfit); % new call: dispModelFitToTvC(d,dout,dataFolder,visArea)
end

% compute sigma and p for extreme values
[sigmaDistributed(2) sigmaAttended(2)] = CRF2TvCfitSigma(d,poolingRule,tau(2),crftype,[],dataFilenames,recalc);
% display fit if called for
if dispFit
  dfit = CRF2TvCfitDeltaC(d,poolingRule,tau(2),sigmaAttended(2),crftype,dataFilenames,recalc);
  disp(sprintf('**** TAU: %f sigmaAttended: %f sigmaDistributed: %f',tau(2),sigmaAttended(2),sigmaDistributed(2)));
  dispModelFitToTvC(d,dfit); % new call: dispModelFitToTvC(d,dout,dataFolder,visArea)
end

midSigmaDistributed = -inf;midSigmaAttended = inf;
epsilon = 0.0001;
nBisections = 0;maxBisections = 25;
% now do bisection search
while (abs(midSigmaDistributed-(1+snrImprovementWithAttention)*midSigmaAttended) > epsilon) && (nBisections < maxBisections)
  midTau = mean(tau);
  [midSigmaDistributed midSigmaAttended] = CRF2TvCfitSigma(d,poolingRule,midTau,crftype);
  if midSigmaDistributed-((1+snrImprovementWithAttention)*midSigmaAttended) > 0
    tau(2) = midTau;
    sigmaDistributed(2) = midSigmaDistributed;
    sigmaAttended(2) = midSigmaAttended;
  else
    tau(1) = midTau;
    sigmaDistributed(1) = midSigmaDistributed;
    sigmaAttended(1) = midSigmaAttended;
  end
  nBisections = nBisections+1;

  if dispFit
    % display fit as we go along
    dfit = CRF2TvCfitDeltaC(d,poolingRule,midTau,midSigmaAttended,crftype,dataFilenames,recalc);
    disp(sprintf('**** TAU: %f sigmaAttended: %f sigmaDistributed: %f',midTau,midSigmaAttended,midSigmaDistributed));
    dispModelFitToTvC(d,dfit); % new call: dispModelFitToTvC(d,dout,dataFolder,visArea)
  end
end

% output arguments
tau = midTau;
sigmaa = midSigmaAttended;
sigmad = midSigmaDistributed;

sigma = [sigmaa*(1+snrImprovementWithAttention) sigmaa];

% do the fit
dfit = CRF2TvCfitDeltaC(d,poolingRule,tau,sigma,crftype,dataFilenames,recalc);

% save file
eval(sprintf('save %s tau sigma sigmad sigmaa dfit',filename));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   getSigmaAcrossObserversAndAreas   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [sigmaDistributed sigmaAttended] = getSigmaAcrossObserversAndAreas(dataFilenames,poolingRule,tau,crftype,recalc,observers)

if ieNotDefined('observers'),observers = 1:4;,end
filename = fullfile(dataFilenames.saveDataDir,fixBadChars(sprintf('sigma%scrf%i',poolingRuleNum2Str(poolingRule,tau),crftype)));

if ~isfile(sprintf('%s.mat',filename)) || recalc
  for observerNum = observers
    for visualArea = 1:4
      disp(sprintf('(crfaaMaxPool:getSigmaAcrossObserversAndAreas) Observer: %i Visual Area: %i',observerNum,visualArea));
      d = getData(dataFilenames,observerNum,visualArea,0);
      [sigmaDistributed(observerNum,visualArea) sigmaAttended(observerNum,visualArea)] = CRF2TvCfitSigma(d,poolingRule,tau,crftype,[],dataFilenames);
      eval(sprintf('save %s sigmaDistributed sigmaAttended;',filename));
    end
  end
  eval(sprintf('save %s sigmaDistributed sigmaAttended;',filename));
else
  disp(sprintf('(crfaaMaxPool:getSigmaAcrossObserversAndAreas) Loading %s',filename));
  load(filename);
  % recompute if there is a problem with the data
  if (size(sigmaDistributed,1) < max(observers))
    [sigmaDistributed sigmaAttended] = getSigmaAcrossObserversAndAreas(dataFilenames,poolingRule,tau,crftype,1,observers);
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   getCrossValSigmaAcrossObserversAndAreas   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [sigmaCombined dfit] = getCrossValSigmaAcrossObserversAndAreas(dataFilenames,poolingRule,tau,crftype,recalc,observers,areas,snrImprovementWithAttention)

if ieNotDefined('observers'),observers = 1:4;,end
if ieNotDefined('areas'),areas = 1:4;,end
if ieNotDefined('snrImprovementWithAttention'),snrImprovementWithAttention=0;end

filename = fullfile(dataFilenames.saveDataDir,fixBadChars(sprintf('sigmaXval%scrf%io%sa%ssnr%s',poolingRuleNum2Str(poolingRule,tau),crftype,mynum2str(observers),mynum2str(areas),mynum2str(snrImprovementWithAttention))));

if ~isfile(sprintf('%s.mat',filename)) || recalc
  for observerNum = observers
    for visualArea = areas
      disp(sprintf('(crfaaMaxPool:getCrossValSigmaAcrossObserversAndAreas) Observer: %i Visual Area: %i',observerNum,visualArea));
      d = getData(dataFilenames,observerNum,visualArea,0);
      [sigmaDistributed sigmaAttended] = CRF2TvCfitSigma(d,poolingRule,tau,crftype,1,dataFilenames,recalc,snrImprovementWithAttention);
      % just return attended sigma (for snrImprovementWithAttention > 0), the two will be different, but
      % by the (1+snrImprovementWithAttention)
      sigmaCombined(observerNum,visualArea) = sigmaAttended;
      [sigmaDistributed sigmaAttended] = CRF2TvCfitSigma(d,poolingRule,tau,crftype,0,dataFilenames,recalc,snrImprovementWithAttention);
      fitSigma = [sigmaAttended*(1+snrImprovementWithAttention) sigmaAttended];
      dfit(observerNum,visualArea) = CRF2TvCfitDeltaC(d,poolingRule,tau,fitSigma,crftype,dataFilenames,recalc);
      eval(sprintf('save %s sigmaCombined dfit;',filename));
    end
  end
  eval(sprintf('save %s sigmaCombined dfit;',filename));
else
  disp(sprintf('(crfaaMaxPool:getCrossValSigmaAcrossObserversAndAreas) Loading %s',filename));
  load(filename);
  % recompute if there is a problem with the data
  if (size(sigmaCombined,1) < max(observers)) || (size(sigmaCombined,2) < max(areas))
    [sigmaCombined dfit] = getCrossValSigmaAcrossObserversAndAreas(dataFilenames,poolingRule,tau,crftype,1,observers,areas,snrImprovementWithAttention);
  end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   computeModelPCorrect   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [areaUnderROC dout] = computeModelPcorrect(responseMeanInt1,responseMeanInt2,sigma,poolingRule,tau,dispFig)

% number of samples in distributions -- the more the better the estimate
nSamples = 1000000;
% FIX FIX FIX
nSamples = 100000;
nSamples = 10000;

% number of criterion points to compute the ROC curve from -- the more the
% finer the ROC curve will be sampled
nCriterion = 10000; % make sure this is an even number

% number of different stimuli and number of different attended locations. i.e. there
% is one attended location for attended trials, for distributed trials there are 4 different
% attended locs as coded above
nStimuli = size(responseMeanInt1,1);
nAttendedLocations = size(responseMeanInt1,2);
% just a short cut 
nResponsesTotal = nStimuli*nAttendedLocations;

% make the response distributions for the first and second intervals, one row for each stimulus
responseDistInt1 = diag(responseMeanInt1(:))*ones(nResponsesTotal,nSamples)+sigma*randn(nResponsesTotal,nSamples);
responseDistInt2 = diag(responseMeanInt2(:))*ones(nResponsesTotal,nSamples)+sigma*randn(nResponsesTotal,nSamples);


% subtract off min for exponentSum (otherwise, exponent gives imaginary values)
if any(poolingRule == [6 7 8 9 10 11])
  minResponse = min([responseDistInt1(:) ; responseDistInt2(:)]);
  responseDistInt1 = responseDistInt1-minResponse;
  responseDistInt2 = responseDistInt2-minResponse;
end

% reshape the response distributions into a nStimuli x nAttendedLocations x nSamples
responseDistInt1 = reshape(responseDistInt1,nStimuli,nAttendedLocations,nSamples);
responseDistInt2 = reshape(responseDistInt2,nStimuli,nAttendedLocations,nSamples);

if nargout == 2
  dout.responseDistInt1 = responseDistInt1;
  dout.responseDistInt2 = responseDistInt2;
end

% now pool the responses across all 4 stimuli according to which pooling rule was selected
if poolingRule == 1
  % max rule
  responseDistInt1 = squeeze(max(responseDistInt1,[],1));
  responseDistInt2 = squeeze(max(responseDistInt2,[],1));
elseif poolingRule == 2
  % mean-rule
  responseDistInt1 = squeeze(mean(responseDistInt1,1));
  responseDistInt2 = squeeze(mean(responseDistInt2,1));
elseif poolingRule == 3
  % soft-max rule
  responseDistInt1 = softmax(responseDistInt1,tau).*responseDistInt1;
  responseDistInt1 = squeeze(sum(responseDistInt1,1));
  responseDistInt2 = softmax(responseDistInt2,tau).*responseDistInt2;
  responseDistInt2 = squeeze(sum(responseDistInt2,1));
elseif poolingRule == 4
  % taget-loc
  responseDistInt1 = squeeze(responseDistInt1(1,:,:));
  responseDistInt2 = squeeze(responseDistInt2(1,:,:));
elseif poolingRule == 5
  % taget-loc vs averaging
  responseDistInt1 = squeeze(mean(responseDistInt1,1));
  responseDistInt2 = squeeze(mean(responseDistInt2,1));
elseif any(poolingRule == [6 7 8 9 10 11])
  % exponent sum rule
  responseDistInt1 = exponentSum(responseDistInt1,tau);
  responseDistInt2 = exponentSum(responseDistInt2,tau);
else
  disp(sprintf('(crfaaMaxPool:computeModelPcorrect) Unknown poolingRule %i',poolingRule));
end

% check to see if we need to transpose
if size(responseDistInt1,2) ~= nSamples
  responseDistInt1 = responseDistInt1';
  responseDistInt2 = responseDistInt2';
end

if nargout == 2
  dout.pooledResponseDistInt1 = responseDistInt1;
  dout.pooledResponseDistInt2 = responseDistInt2;
end

% now compute min/max and step size for stepping through criterion
minResponse = min([responseDistInt1(:)' responseDistInt2(:)']);
maxResponse = max([responseDistInt1(:)' responseDistInt2(:)']);
responseStep = (maxResponse-minResponse)/(nCriterion-1);
criterions = minResponse:responseStep:maxResponse;

% Now get hits and false alarams for ROC function.
% i.e. compute the proportion of hist and false alarms as a function
% of the criterion (i.e. each element in hits/false alarms is for
% a criterion as specified in the criterions array). 
hits = max(0,1-cumsum(histc(responseDistInt1,criterions,2)/nSamples,2));
falseAlarms = max(0,1-cumsum(histc(responseDistInt2,criterions,2)/nSamples,2));

% return hits and false alarms
dout.hits = hits;
dout.falseAlarms = falseAlarms;

% compute area under ROC
for i = 1:nAttendedLocations
  areaUnderROC(i) = abs(trapz(falseAlarms(i,:),hits(i,:)));
end

% performance can be comupted as a difference between trials instead
% of area under the ROC
dout.areaUnderROC = sum((responseDistInt1-responseDistInt2)>0,2)/nSamples;

% return unaveraged area under ROC
%dout.areaUnderROC = areaUnderROC;

if dispFig
  % make some polts
  smartfig('crfaaMaxPool','reuse');
  clf;
  nPlots = nAttendedLocations+1;
  for i = 1:nAttendedLocations
    subplot(1,nPlots,i);
    myhist(responseDistInt1(i,:),50,'r');
    plotmean(responseDistInt1(i,:),'r',1);
    xlabel(sprintf('Response strength\n(Pooled %%signal change)'));
    ylabel('N');

    myhist(responseDistInt2(i,:),50);
    plotmean(responseDistInt2(i,:),'k');
%    if nAttendedLocations > 1
      title(sprintf('Signal and Noise (n=%i)\n[%0.3f vs %0.3f] %0.3f %0.3f %0.3f',nSamples,responseMeanInt1(1,i),responseMeanInt2(1,i),responseMeanInt1(2,i),responseMeanInt1(3,i),responseMeanInt1(4,i)));
%    else
%      title(sprintf('Signal and Noise (n=%i)\n[%0.3f vs %0.3f] %0.3f %0.3f %0.3f',nSamples,responseMeanInt1(1),responseMeanInt2(1),responseMeanInt1(2),responseMeanInt1(3),responseMeanInt1(4)));
%    end
    ylabel('N');
    xaxis(minResponse,maxResponse);
  end
  
  subplot(1,nPlots,nPlots);
  plot(falseAlarms',hits');
  hold on
  axis square
  dline;
  xlabel('False Alarm rate (proportion)');
  ylabel('Hit rate (proportion)');
  if nAttendedLocations > 1
    title(sprintf('Pooling rule: %s sigma: %0.5f\nArea under ROC: %s (%s)',poolingRuleNum2Str(poolingRule,tau),sigma,mynum2str(areaUnderROC,'sigfigs=3'),mynum2str(mean(areaUnderROC),'sigfigs=3')));
  else
    title(sprintf('Pooling rule: %s sigma: %0.5f\nArea under ROC: %s',poolingRuleNum2Str(poolingRule,tau),sigma,mynum2str(areaUnderROC,'sigfigs=3')));
  end
  drawnow
end

% compute the average areaUnderROC (this is only relevant for the distributed condition in which 
% we have 4 attended conditions -- one for attending to each stimulus). Then the performance
% is equal to the average performance across all 4 effective attention conditions.
areaUnderROC = mean(areaUnderROC);


%%%%%%%%%%%%%%%%%
%   softmax   %%
%%%%%%%%%%%%%%%%%
function p = softmax(q,tau)

% default tau setting
if nargin == 1,tau = 1;end

% softmax operates column-wise, that is, it computes the softmax operator 
% across each column -- converting each element into a weight such that
% the weights across each column sum to 1. With tau set small (say 0.1), this operates
% like max weighting (i.e. the entry in the row with the largest value approach 1 and
% all others approach 0). With tau set large (say 10), this acts like average weighting
% where each element in a row approaches the same value.
% make sure to pass in a *column* array e.g. [3 4 1 5]'
expq = exp(q/tau);
p = expq./ repmat(sum(expq,1),size(q,1),1);


%%%%%%%%%%%%%%%%%%
%   exponentSum  %
%%%%%%%%%%%%%%%%%%
function [r w] = exponentSum(r,exponent)
% this is the fit suggested by david
% it is equivalent to the softmax operator
% but it has a biologically plausible implementation
% Normalization.

% default exponent setting (does averaging)
if nargin == 1,exponent = 1;end

% compute what the effective weights are
if nargout >= 2
  w = ((r.^exponent/sum(r.^exponent,1))./abs((r./sum(r,1))))/4;
end

% compute response as sum of responses (column-wise) to the exponent.
% normalization is to divide by 4 and undo the exponent
r = (squeeze(sum(r.^exponent,1))/size(r,1)).^(1/exponent);


%%%%%%%%%%%%%%%%%%%%%%%
% interpolateResponse %
%%%%%%%%%%%%%%%%%%%%%%%
function response = interpolateResponse(modelContrast,modelResponse,contrast)
% get responses from a contrast response fit
% this uses linear interpolation

response = interp1(modelContrast,modelResponse,contrast);

% use interpolation for points that are outside linear interpolation bounds (i.e.
% contrasts below the lowest tested.
badResponse = find(isnan(response));

% spline interpolation gave some funny non-monotonic results, so commenting out.
%  response(badResponse) = interp1(modelContrast,modelResponse,contrast(badResponse),'spline');

% instead, linearlly interpolate using the line that fits the closest
% two points on the CRF
if ~isempty(badResponse)
  for badI = 1:length(badResponse)
    % find the closest contrast to the missing one
    [unusedTempValue closestContrastIndex] = sort(abs(modelContrast-contrast(badResponse(badI))));
    % fit with the closest two points
    fitIndex = closestContrastIndex(1:2);
    % now fit a line to the values closest to the unknown contrast
    response(badResponse(badI)) = linearlyExtrapolate(modelContrast(fitIndex)',modelResponse(fitIndex)',contrast(badResponse(badI)));
  end  
end


%%%%%%%%%%%%%%%%%%%%%%%
% linearlyExtrapolate %
%%%%%%%%%%%%%%%%%%%%%%%
function yi = linearlyExtrapolate(x,y,xi)
%   fit a line to the x,y points and then extrapolate the xi points   %%

% fit the line
A = ones(size(x,1),2);
A(:,1) = x;
b = pinv(A)*y;

A = ones(size(xi,1),2);
A(:,1) = xi;
yi = A*b;


%%%%%%%%%%%%%%%%%%%%%%
%   getResponses   %%
%%%%%%%%%%%%%%%%%%%%%%
function [response1 response2] = getResponses(type,d,pedestalNum,altdata,poolingRule,deltaC)
% altdata is passed in to make this function work differently. If altdata is absent,
% then we do "normal" processing, which means that we get the deltaC from the TvC
% and we get the unattended and attended responses from the d structure.

alttype = 'normal';
if (nargin >= 4) && ~isempty(altdata)
  % if altdata is a tvc, then we do tvc processing, which means we get the CRF that
  % is inferred from the tvc function
  if isfield(altdata,'istvc')
    alttype = 'tvc';
  % if altdata is crf then we get the CRF and the deltaC from the altdata structure
  elseif isfield(altdata,'iscrf')
    alttype = 'crf';
  else
    disp(sprintf('(crfaaMaxPool:getResponse) Unknown altdata type'));
    keyboard
  end
end

% get delta C
if ~strcmp(alttype,'crf')
  deltaC = d.(type).TvC(pedestalNum);
else
  if nargin < 6
    deltaC = altdata.deltaC;
  end
end

% get pedestalContrast
pedestalContrast = d.pedestalContrasts(pedestalNum);


% get distractor contrsats
distractorContrasts = d.distractorContrast{pedestalNum};

% for pooling rule 10, we just test what happens with distactors at lower contrast
if poolingRule == 10
  distractorContrasts(:) = min(distractorContrasts);
  disp(sprintf('(crfaaMaxPool:getResponses) Setting distractors to lower contrast %s',mynum2str(distractorContrasts)));
end

% now we can get all 5 contrasts we need responses for
contrasts = [pedestalContrast pedestalContrast+deltaC distractorContrasts];

if strcmp(alttype,'normal')
  % get the attended and unattended responses to those contrasts
  unattendedResponses = interpolateResponse(d.unattendedContrast,d.unattendedResponse,contrasts);
  attendedResponses = interpolateResponse(d.attendedContrast,d.attendedResponse,contrasts);
  % this is only used for poolingRule 4 which is the old type
  distributedResponses = interpolateResponse(d.distributed.pedestalsCRF,d.distributed.CRFtarget,contrasts);
elseif strcmp(alttype,'tvc')
  % get the values for the unattended and attended responses from the current tvc.
  unattendedResponses = interpolateResponseFromTvC(d.pedestalContrasts,tvc.unattended.deltaR,tvc.unattended.deltaC,tvc.unattended.responseOffset,contrasts);
  attendedResponses = interpolateResponseFromTvC(d.pedestalContrasts,tvc.attended.deltaR,tvc.attended.deltaC,tvc.attended.responseOffset,contrasts);
  if poolingRule == 4
    disp(sprintf('(crfaaMaxPool:getResponse) PoolingRule 4 not implemented yet for tvc2crf fit'));
    keyboard
  end
elseif strcmp(alttype,'crf')
  unattendedResponses = interpolateResponse(altdata.unattendedContrast,altdata.unattendedResponse,contrasts);
  attendedResponses = interpolateResponse(altdata.attendedContrast,altdata.attendedResponse,contrasts);
  distributedResponses = interpolateResponse(altdata.distributedContrast,altdata.distributedResponse,contrasts);
  distributedNontargetResponses = interpolateResponse(altdata.distributedContrast,altdata.distributedNontargetResponse,contrasts);
end

if strcmp(type,'distributed')
  % poolingRule 4 is targetLoc, so we use the distributed responses
  if any(poolingRule == [4 5])
    response1 = [distributedResponses(2) distributedResponses(3:5)]';
    response2 = [distributedResponses(1) distributedResponses(3:5)]';
  % pooling Rule selects concurrently from two locaations
  elseif poolingRule == 7
    % now construct response for unattended conditions. Same as for attended conditions
    % but each column now has an attention to a different stimulus
    response1(:,1) = [attendedResponses(2) attendedResponses(3) unattendedResponses(4) unattendedResponses(5)]';
    response2(:,1) = [attendedResponses(1) attendedResponses(3) unattendedResponses(4) unattendedResponses(5)]';

    response1(:,2) = [attendedResponses(2) unattendedResponses(3) attendedResponses(4) unattendedResponses(5)]';
    response2(:,2) = [attendedResponses(1) unattendedResponses(3) attendedResponses(4) unattendedResponses(5)]';

    response1(:,3) = [attendedResponses(2) unattendedResponses(3) unattendedResponses(4) attendedResponses(5)]';
    response2(:,3) = [attendedResponses(1) unattendedResponses(3) unattendedResponses(4) attendedResponses(5)]';

    response1(:,4) = [unattendedResponses(2) attendedResponses(3) attendedResponses(4) unattendedResponses(5)]';
    response2(:,4) = [unattendedResponses(1) attendedResponses(3) attendedResponses(4) unattendedResponses(5)]';

    response1(:,5) = [unattendedResponses(2) attendedResponses(3) unattendedResponses(4) attendedResponses(5)]';
    response2(:,5) = [unattendedResponses(1) attendedResponses(3) unattendedResponses(4) attendedResponses(5)]';

    response1(:,6) = [unattendedResponses(2) unattendedResponses(3) attendedResponses(4) attendedResponses(5)]';
    response2(:,6) = [unattendedResponses(1) unattendedResponses(3) attendedResponses(4) attendedResponses(5)]';
  elseif poolingRule == 8
    % now construct response for unattended conditions. Same as for attended conditions
    % but each column now has an attention to a different stimulus
    response1(:,1) = [attendedResponses(2) attendedResponses(3) attendedResponses(4) unattendedResponses(5)]';
    response2(:,1) = [attendedResponses(1) attendedResponses(3) attendedResponses(4) unattendedResponses(5)]';

    response1(:,2) = [attendedResponses(2) unattendedResponses(3) attendedResponses(4) attendedResponses(5)]';
    response2(:,2) = [attendedResponses(1) unattendedResponses(3) attendedResponses(4) attendedResponses(5)]';

    response1(:,3) = [attendedResponses(2) attendedResponses(3) unattendedResponses(4) attendedResponses(5)]';
    response2(:,3) = [attendedResponses(1) attendedResponses(3) unattendedResponses(4) attendedResponses(5)]';

    response1(:,4) = [unattendedResponses(2) attendedResponses(3) attendedResponses(4) attendedResponses(5)]';
    response2(:,4) = [unattendedResponses(1) attendedResponses(3) attendedResponses(4) attendedResponses(5)]';
  elseif any(poolingRule == [9 10])
    % now construct response for unattended conditions. Same as for attended conditions
    % but each column now has an attention to a different stimulus
    response1(:,1) = [distributedResponses(2) distributedResponses(3) distributedResponses(4) distributedResponses(5)]';
    response2(:,1) = [distributedResponses(1) distributedResponses(3) distributedResponses(4) distributedResponses(5)]';
  elseif any(poolingRule == [11])
    % now construct response for unattended conditions, based on distributed cue target and non-target responses 
    response1(:,1) = [distributedNontargetResponses(2) distributedNontargetResponses(3) distributedNontargetResponses(4) distributedNontargetResponses(5)]';
    response2(:,1) = [distributedNontargetResponses(1) distributedNontargetResponses(3) distributedNontargetResponses(4) distributedNontargetResponses(5)]';
  else
    % now construct response for unattended conditions. Same as for attended conditions
    % but each column now has an attention to a different stimulus
    response1(:,1) = [attendedResponses(2) unattendedResponses(3:5)]';
    response2(:,1) = [attendedResponses(1) unattendedResponses(3:5)]';
    
    response1(:,2) = [unattendedResponses(2) attendedResponses(3) unattendedResponses([4 5])]';
    response2(:,2) = [unattendedResponses(1) attendedResponses(3) unattendedResponses([4 5])]';
    
    response1(:,3) = [unattendedResponses(2) unattendedResponses(3) attendedResponses(4) unattendedResponses(5)]';
    response2(:,3) = [unattendedResponses(1) unattendedResponses(3) attendedResponses(4) unattendedResponses(5)]';
    
    response1(:,4) = [unattendedResponses(2) unattendedResponses([3 4]) attendedResponses(5)]';
    response2(:,4) = [unattendedResponses(1) unattendedResponses([3 4]) attendedResponses(5)]';
  end
elseif strcmp(type,'attended')
  if poolingRule ~= 5
    % now construct response model for attended condition. This is the attended response
    % to pedestal+deltaC in interval 1 + all the distractor responses
    % and for interval 2 is the pedestal response + all distractor responses
    response1 = [attendedResponses(2) unattendedResponses(3:5)]';
    response2 = [attendedResponses(1) unattendedResponses(3:5)]';
  else
    response1 = attendedResponses(2);
    response2 = attendedResponses(1);
  end
else
  disp(sprintf('(crfaaMaxPool:getResponse) Unknown type %s',type));
  keyboard
end


%%%%%%%%%%%%%%%
%   getData   %
%%%%%%%%%%%%%%%
function dout = getData(dataFilenames,observerNum,visualArea,dispFig)
% observerNum is 1 for average and 2-4 for others
% visualArea is 1-4 for which visual area

% display figure with data or not
if ieNotDefined('dispFig')
  dispFig = 1;
end

  
% load data file
decisionModelFilename = fullfile(dataFilenames.path,dataFilenames.decisionModel);
if ~isfile(decisionModelFilename)
  disp(sprintf('(crfaaMaxPool:getData) Could not find file %s',decisionModelFilename));
  keyboard
end
d = load(decisionModelFilename);

% set info about which data set we are using
dout.observerNum = observerNum;
visualAreas = {'V1','V2','V3','V4'};
dout.visualArea = visualAreas{visualArea};

% first get some results from the distributed condition
testType = 1;
% check bootstrap number
if isempty(dataFilenames.bootstrapNum)
  if length(size(d.results)) == 4
    disp(sprintf('(crfaaMaxPool:getData) Data file %s is a bootstrap file with %i boostraps, but no bootstapNum was specified',decisionModelFilename,size(d.results,1)));
    keyboard
  end
  results = d.results{observerNum,visualArea,testType};
else
  results = d.results{dataFilenames.bootstrapNum,observerNum,visualArea,testType};
end
dout.distributed.TvC = results.behavior.tvc.thisTvC;
dout.distributed.TvCste = results.behavior.tvc.thisTvCste;
dout.distributed.pedestalsTvC = results.pedestalsTvC;
dout.distributed.pedestalsCRF = results.pedestalsCRF;
dout.distributed.CRFtarget = results.crf_distributed_target;
dout.distributed.CRFnontarget = results.crf_distributed_nontarget;
dout.distributed.sigma = results.noise.k;
dout.distributed.responseOffset = results.noise.responseOffset;
dout.distributed.modelContrast = results.noise.contrast;
dout.distributed.modelResponse = results.noise.response;
if isfield(results.behavior.tvcfit.fitParams,'fit')
  dout.distributed.TvCFitY = results.behavior.tvcfit.fitParams.fit.y;
  dout.distributed.TvCFitX = results.behavior.tvcfit.fitParams.fit.x;
else
  dout.distributed.TvCFitY = [];
  dout.distributed.TvCFitX = [];
end
% check array orientation
if size(dout.distributed.TvC,1)  == 1
  dout.distributed.TvC = dout.distributed.TvC';
  dout.distributed.TvCste = dout.distributed.TvCste';
end

% now get some results from attended condition
testType = 3;
if isempty(dataFilenames.bootstrapNum)
  results = d.results{observerNum,visualArea,testType};
else
  results = d.results{dataFilenames.bootstrapNum,observerNum,visualArea,testType};
end
dout.attended.TvC = results.behavior.tvc.thisTvC;
dout.attended.TvCste = results.behavior.tvc.thisTvCste;
dout.attended.pedestalsTvC = results.pedestalsTvC;
dout.attended.pedestalsCRF = results.pedestalsCRF;
dout.attended.CRF = results.crf_attended;
dout.attended.sigma = results.noise.k;
dout.attended.responseOffset = results.noise.responseOffset;
dout.attended.modelContrast = results.noise.contrast;
dout.attended.modelResponse = results.noise.response;
% no fit to draw for polynomial fits
if isfield(results.behavior.tvcfit.fitParams,'fit')
  dout.attended.TvCFitY = results.behavior.tvcfit.fitParams.fit.y;
  dout.attended.TvCFitX = results.behavior.tvcfit.fitParams.fit.x;
else
  dout.attended.TvCFitY = [];
  dout.attended.TvCFitX = [];
end
% check array orientation
if size(dout.attended.TvC,1)  == 1
  dout.attended.TvC = dout.attended.TvC';
  dout.attended.TvCste = dout.attended.TvCste';
end

% unattended information
dout.unattended.pedestalsCRF = results.pedestals;
dout.unattended.CRF = results.crf_unattended;

% get the distractor contrasts
dout.distractorContrast = results.distractorContrast;

% get the pedestals -- remmeber that first pedestal should be 0 (and not the
% small number that it is set to to make the plotting work).
dout.pedestalContrasts = dout.attended.pedestalsTvC;
dout.pedestalContrasts(1) = 0;
dout.nPedestals = length(dout.pedestalContrasts);

% initialize unattended CRF
dout.unattendedContrast = dout.unattended.pedestalsCRF;
dout.unattendedResponse = dout.unattended.CRF;

% initialize attended CRF
dout.attendedContrast = dout.attended.pedestalsCRF;
dout.attendedResponse = dout.attended.CRF;

% load crffit data
crffitFilename = fullfile(dataFilenames.path,dataFilenames.crffit);
if isempty(dataFilenames.crffit)
  disp(sprintf('(crfaaMaxPool:getData) No crffit data - only used when using a crfypte that does not use model fits ',crffitFilename));
  dout.hasFits = 0;
elseif ~isfile(crffitFilename)
  disp(sprintf('(crfaaMaxPool:getData) Could not find file %s',crffitFilename));
  dout.hasFits = 0;
else
  dcrffit = load(crffitFilename);

  % grab the fits of the full model -- the first number in the cell array is irrelevant (coding
  % issue -- all thre cell arrays should contain the same info). If you want to get the
  % fit with just the offset shifting, you can look at dcrffit.results{3} -- this is the fit 
  % in the paper).
  crffit = dcrffit.fullResults{1}{observerNum,visualArea}.crf.fit;

  % now get each one of the crffits
  dout.attended.fitContrast = crffit.fitx(1,:);
  dout.attended.fitResponse = crffit.fity(1,:);

  dout.distributed.fitContrast = crffit.fitx(2,:);
  dout.distributed.fitResponse = crffit.fity(2,:);
  
  dout.unattended.fitContrast = crffit.fitx(4,:);
  dout.unattended.fitResponse = crffit.fity(4,:);
  dout.hasFits = 1;
end


% display figure
if dispFig
  smartfig('crfaaMaxPoolData','reuse');
  clf;
  if dout.hasFits,numPlots = 3;else, numPlots = 3; end
  subplot(1,numPlots,1);
  loglog(dout.distributed.pedestalsTvC,dout.distributed.TvC,'k.');
  hold on;
  loglog(dout.distributed.TvCFitX,dout.distributed.TvCFitY,'k-');
  loglog(dout.attended.pedestalsTvC,dout.attended.TvC,'r.');
  loglog(dout.attended.TvCFitX,dout.attended.TvCFitY,'r-');
  xlabel('Contrast');
  ylabel('deltaC');
  title(sprintf('Contrast discrimination functions (observer=%i)',dout.observerNum));

  subplot(1,numPlots,2);
  semilogx(dout.distributed.pedestalsCRF,dout.distributed.CRFtarget,'k.--');
  hold on
  semilogx(dout.distributed.pedestalsCRF,dout.distributed.CRFnontarget,'k.');
  semilogx(dout.distributed.modelContrast,dout.distributed.modelResponse,'k-');
  semilogx(dout.attended.pedestalsCRF,dout.attended.CRF,'r.');
  semilogx(dout.attended.modelContrast,dout.attended.modelResponse,'r-');
  semilogx(dout.unattended.pedestalsCRF,dout.unattended.CRF,'g.-');
  xlabel('Contrast');
  ylabel('Response');
  title(sprintf('Contrast response functions (%s)',dout.visualArea));

  if dout.hasFits 
    subplot(1,numPlots,3);
    semilogx(dout.distributed.pedestalsCRF,dout.distributed.CRFtarget,'k.');
    hold on
    semilogx(dout.distributed.fitContrast,dout.distributed.fitResponse,'k-');
    semilogx(dout.attended.pedestalsCRF,dout.attended.CRF,'r.');
    semilogx(dout.attended.fitContrast,dout.attended.fitResponse,'r-');
    semilogx(dout.unattended.pedestalsCRF,dout.unattended.CRF,'g.');
    semilogx(dout.unattended.fitContrast,dout.unattended.fitResponse,'g-');
    xlabel('Contrast');
    ylabel('Response');
    title(sprintf('Contrast response functions\nwith full model fit(%s)',dout.visualArea));
  end
end


%%%%%%%%%%%%%%%%%%%%
% bisectionTvC2CRF %
%%%%%%%%%%%%%%%%%%%%
function [unattendedResponses attendedResponses] = bisectionTvC2CRF(targetPCorrect,d,sigma,poolingRule,tau)
%   bisectionTvC2CRF: This function goes from TvC's to CRF's. Its a bit tricky
%   since you need to know the full CRF to compute the model, I abandoned it
%   after switching to going in the opposite direction, so it likely needs
%   to be debugged and toyed with if it is to be resurrected. It uses
%   the interpolateResponseFromTvC function below.

% we will search for deltaR's between the min and max
minDeltaR = 0.001;
maxDeltaR = 0.5;

% this function takes model parameters sigma, poolingRule, tau, responseOffset 
% and computes the implied CRF
deltaRattended(1:d.nPedestals,1) = minDeltaR;
deltaRattended(1:d.nPedestals,2) = maxDeltaR;
deltaRdistributed(1:d.nPedestals,1) = minDeltaR;
deltaRdistributed(1:d.nPedestals,2) = maxDeltaR;

% set how close we have to be the targetPCorrect to terminate
epsilon = 0.01;

% get parameters of the TvC
tvc.istvc = 1;
tvc.attended.deltaC = d.attended.TvC;
tvc.attended.responseOffset = d.attendedResponse(1);
tvc.unattended.deltaC = d.distributed.TvC;
tvc.unattended.responseOffset = d.unattendedResponse(1);
tvc.attended.deltaR = mean(deltaRattended,2);
tvc.unattended.deltaR = mean(deltaRdistributed,2);

% initialize the probability 
disppercent(-inf,'Initializing pCorrect for bisection search');
for pedestalNum = 1:d.nPedestals
  % set up deltaR's so that all the other pedestals are set to the mid point of their current
  % interval, but this pedestal is set to the lower end of the interval
  tvc.attended.deltaR = mean(deltaRattended,2);
  tvc.unattended.deltaR = mean(deltaRdistributed,2);
  tvc.unattended.deltaR(pedestalNum) = deltaRdistributed(pedestalNum,1);
  [responseDistributedInt1 responseDistributedInt2] = getResponses('distributed',d,pedestalNum,tvc,poolingRule);
  % compute model
  plower(1,pedestalNum) = computeModelPcorrect(responseDistributedInt1,responseDistributedInt2,sigma,poolingRule,tau,0);
  % same for attended
  tvc.attended.deltaR = mean(deltaRattended,2);
  tvc.attended.deltaR(pedestalNum) = deltaRattended(pedestalNum,1);
  tvc.unattended.deltaR = mean(deltaRdistributed,2);
  [responseAttendedInt1 responseAttendedInt2] = getResponses('attended',d,pedestalNum,tvc,poolingRule);
  plower(2,pedestalNum) = computeModelPcorrect(responseAttendedInt1,responseAttendedInt2,sigma,poolingRule,tau,0);
  
  % now compute the same thing for the upper bound (same as above, except that we save
  % teh result in pupper/ not plower and use the deltaR's second column and not its 1st column
  tvc.attended.deltaR = mean(deltaRattended,2);
  tvc.unattended.deltaR = mean(deltaRdistributed,2);
  tvc.unattended.deltaR(pedestalNum) = deltaRdistributed(pedestalNum,2);
  [responseDistributedInt1 responseDistributedInt2] = getResponses('distributed',d,pedestalNum,tvc,poolingRule);
  pupper(1,pedestalNum) = computeModelPcorrect(responseDistributedInt1,responseDistributedInt2,sigma,poolingRule,tau,0);
  tvc.attended.deltaR = mean(deltaRattended,2);
  tvc.attended.deltaR(pedestalNum) = deltaRattended(pedestalNum,2);
  tvc.unattended.deltaR = mean(deltaRdistributed,2);
  [responseAttendedInt1 responseAttendedInt2] = getResponses('attended',d,pedestalNum,tvc,poolingRule);
  pupper(2,pedestalNum) = computeModelPcorrect(responseAttendedInt1,responseAttendedInt2,sigma,poolingRule,tau,0);
  disppercent(pedestalNum/d.nPedestals);
end
disppercent(inf);

% make sure our target is between the p bounds computed here
if any(plower(:) > targetPCorrect) || any(pupper(:) < targetPCorrect)
  % show what the pCorrect were set to
  for pedestalNum = 1:d.nPedestals
    disp(sprintf('Pedestal %i distibuted: dR=[%0.8f %0.8f] p=[%0.4f %0.4f] attended: dR=[%0.8f %0.8f] p=[%0.4f %0.4f]',pedestalNum,deltaRdistributed(pedestalNum,1),deltaRdistributed(pedestalNum,2),plower(1,pedestalNum),pupper(1,pedestalNum),deltaRattended(pedestalNum,1),deltaRattended(pedestalNum,2),plower(2,pedestalNum),pupper(2,pedestalNum)));
  end
  % display error
  disp(sprintf('(crfaa:bisectionTvC2CRF) p-value out of range -- need to expand deltaR min max?'));
  % return nan
  unattendedResponses = nan(1,d.nPedestals)
  attendedResponses = nan(1,d.nPedestals);
  return
end

nBisections = 0;
% give up bisecting after we reach some very small difference in deltaR
maxBisections = round(log2((maxDeltaR-minDeltaR)/0.000001));
% now do bisection search until we are with epsilon of the targetPCorrect
while any(min([abs(plower(:)-targetPCorrect) abs(pupper(:)-targetPCorrect)],[],2)>epsilon) && (nBisections < maxBisections)
  % update current best estimate of deltaR's
  tvc.attended.deltaR = mean(deltaRattended,2);
  tvc.unattended.deltaR = mean(deltaRdistributed,2);
  disppercent(-inf,'Finding mid delta model p correct');
  for pedestalNum = 1:d.nPedestals
    % get the distributed responses for this pedestal
    [responseDistributedInt1 responseDistributedInt2] = getResponses('distributed',d,pedestalNum,tvc,poolingRule);
    % compute model
    p(1,pedestalNum) = computeModelPcorrect(responseDistributedInt1,responseDistributedInt2,sigma,poolingRule,tau,0);
    % get the attended responses for this pedestal
    [responseAttendedInt1 responseAttendedInt2] = getResponses('attended',d,pedestalNum,tvc,poolingRule);
    % compute model
    p(2,pedestalNum) = computeModelPcorrect(responseAttendedInt1,responseAttendedInt2,sigma,poolingRule,tau,0);
    disppercent(pedestalNum/d.nPedestals);
  end
  disppercent(inf);
  nBisections = nBisections+1;
  % now update the intervals
  for pedestalNum = 1:d.nPedestals
    % update distributed p's
    if p(1,pedestalNum) > targetPCorrect
      % set upper interval
      pupper(1,pedestalNum) = p(1,pedestalNum);
      deltaRdistributed(pedestalNum,2) = mean(deltaRdistributed(pedestalNum,1:2));
    else
      % set lower interval
      plower(1,pedestalNum) = p(1,pedestalNum);
      deltaRdistributed(pedestalNum,1) = mean(deltaRdistributed(pedestalNum,1:2));
    end
    % update attended p's
    if p(2,pedestalNum) > targetPCorrect
      % set upper interval
      pupper(2,pedestalNum) = p(2,pedestalNum);
      deltaRattended(pedestalNum,2) = mean(deltaRattended(pedestalNum,1:2));
    else
      % set lower interval
      plower(2,pedestalNum) = p(2,pedestalNum);
      deltaRattended(pedestalNum,1) = mean(deltaRattended(pedestalNum,1:2));
    end
    disp(sprintf('%i: Pedestal %i distibuted: dR=[%0.8f %0.8f] p=[%0.4f %0.4f] attended: dR=[%0.8f %0.8f] p=[%0.4f %0.4f]',nBisections,pedestalNum,deltaRdistributed(pedestalNum,1),deltaRdistributed(pedestalNum,2),plower(1,pedestalNum),pupper(1,pedestalNum),deltaRattended(pedestalNum,1),deltaRattended(pedestalNum,2),plower(2,pedestalNum),pupper(2,pedestalNum)));
  end
end

% get the best response for each pedestal
[minP minI] = min([abs(plower(:)-targetPCorrect) abs(pupper(:)-targetPCorrect)],[],2);
tvc.unattended.deltaR = deltaRdistributed(sub2ind([d.nPedestals 2],1:d.nPedestals,minI(1:d.nPedestals)'))';
tvc.attended.deltaR = deltaRattended(sub2ind([d.nPedestals 2],1:d.nPedestals,minI(d.nPedestals+1:d.nPedestals*2)'))';

% get unattended and attended Responses
unattendedResponses = interpolateResponseFromTvC(d.pedestalContrasts,tvc.unattended.deltaR,tvc.unattended.deltaC,tvc.unattended.responseOffset,d.pedestalContrasts);
attendedResponses = interpolateResponseFromTvC(d.pedestalContrasts,tvc.attended.deltaR,tvc.attended.deltaC,tvc.attended.responseOffset,d.pedestalContrasts);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% interpolateResponseFromTvC %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function response = interpolateResponseFromTvC(pedestalContrasts,deltaR,deltaC,responseOffset,contrast)
% get responses starting from the TvC data, i.e.
% we know the deltaR's that match certain deltaC's
% construct a pieceiwise linear function between
% the pedestal contrasts that matches this.
% Then use interpolateResponse to get the responses
% for the passed in contrast

nPedestals = length(pedestalContrasts);
% first find the piecewise linear function of contrast that matches
% the given deltaR's and deltaC's. We wish to find the difference
% in response between each pedestal contrast (deltaPedestalResponse) that
% does this (and not between each deltaC which is variable). 
% We thus are constructin a square A matrix that does the following
% deltaR = A * deltaPedestalResponse

% First, get the difference in contrast between each pedestal. For the last pedestal
% we (arbitrarily) use the difference between the pedestal contrast and the
% last deltaC. These will be the points upon which our new CRF will be defined.
pedestalContrastDifferences = diff([pedestalContrasts (pedestalContrasts(end)+deltaC(end))]);
% now init the A to zeros
A = zeros(nPedestals,nPedestals);
% multiply each one of these by the weighting factors for each deltaC
for i = 1:nPedestals-1
  % The first part of the delta response is the slope up to the
  % next pedestal contrast (if the deltaC is smaller then the difference
  % in contrast to the next pedestal, then this is all there is to calculating
  % what the deltaPedestalResponse should be.
  A(i,i) = min(deltaC(i)/pedestalContrastDifferences(i),1);
  % but if the deltaC is larger than the difference in contrast to the next
  % contrast we have to add in the slope of the next contrast difference
  A(i,i+1) = max(((deltaC(i)-pedestalContrastDifferences(i))/pedestalContrastDifferences(i+1)),0);
end
% for the last pedestal, we just have the slope define by the TvC
A(nPedestals,nPedestals) = 1;

% now use the A to compute the deltaPedestalResponse
deltaPedestalResponse = inv(A)*deltaR;

% now construct the CRF. We know the difference in responses between
% different pedestals as calculated above so we just take the cumsum
% and add that to the passed in responseOffset
modelContrast = [pedestalContrasts pedestalContrasts(end)+deltaC(end)];
modelResponse = [0 cumsum(deltaPedestalResponse')]+responseOffset;

% now go look up the contrasts using interpolateResponse
response = interpolateResponse(modelContrast,modelResponse,contrast);

% display
if 0
  smartfig('crfaaMaxPoolInterpolateResponseFromTvC','reuse');
  clf;
  plot(modelContrast,modelResponse,'k.-');
  hold on
  vline(deltaC+pedestalContrasts');
  hline(deltaR+modelResponse(1:nPedestals));
  hline(responseOffset,'g:');
  xlabel('Contrast');ylabel('Response');
  hline(response,'r:');
  vline(contrast,'r:');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%
%   poolingRuleNum2Str   %
%%%%%%%%%%%%%%%%%%%%%%%%%%
function poolingRule = poolingRuleNum2Str(poolingRule,tau)

if (nargin == 2) && ~isnan(tau)
  poolingRuleNames = {'max','average',sprintf('softmax %0.2f',tau),'targetLoc','targetLoc vs average',sprintf('exponent sum %0.2f',tau),sprintf('exponent sum %0.2f attend 2',tau),sprintf('exponent sum %0.2f attend 3',tau),sprintf('exponent sum %0.2f attend 4',tau),sprintf('exponent sum %0.2f attend 4 lo c',tau),sprintf('exponent sum no strategy %0.2f',tau)};
else
  poolingRuleNames = {'max','average','softmax','targetLoc','targetLoc vs average','exponent sum','exponent sum attend 2','exponent sum attend 3','exponent sum attend 4','exponent sum attend lo c','exponent sum no-strategy'};
end
poolingRule = poolingRuleNames{poolingRule};
		   

%%%%%%%%%%%%%%%%%
%   fitOffset   %
%%%%%%%%%%%%%%%%%
function offset = fitOffset(fit,data)

% least squares fit of offset is just the mean difference between the two.
% (pretty obvious I guess, but it took a bit  of calculus to convice myself
% this was correct. i.e. define f(k) = sum(fit-data-offset).^2 where offset
% is a single value and fit/data are arrays. Then compute df/dk and set
% to zero and you will see this comes out to be true)
offset = mean(fit-data);

%%%%%%%%%%%%%%%%%%%%%%%%%
%   CRF2TvCfitSigma   %%
%%%%%%%%%%%%%%%%%%%%%%%%%
function [sigmaDistributed sigmaAttended pDistributed pAttended] = CRF2TvCfitSigma(d,poolingRule,tau,crftype,fitSingleSigma,dataFilenames,recompute,snrImprovementWithAttention)

if ieNotDefined('snrImprovementWithAttention'),snrImprovementWithAttention=0;end
% don't preload if recompute is set
if ieNotDefined('recompute'),recompute = 0;end
 
% decided whether to fit a single sigma to both attended/distributed functions
if ieNotDefined('fitSingleSigma'), fitSingleSigma = 0;end
if fitSingleSigma
  disp(sprintf('(crfaaMaxPool:CRF2TvCfitSigma) Fitting a single sigma to attended/distributed'));
end

% name of preload file
filename = '';
if ~ieNotDefined('dataFilenames')
  filename = fullfile(dataFilenames.saveDataDir,fixBadChars(sprintf('fitSigma_%i_%s_%s_crftype%isnr%f_%i',d.observerNum,d.visualArea,poolingRuleNum2Str(poolingRule,tau),crftype,mynum2str(snrImprovementWithAttention),fitSingleSigma),[],{'.','p'}));
end

% look for preload file
if ~recompute
  if isfile(sprintf('%s.mat',filename))
    disp(sprintf('(crfaaMaxPool:CRF2TvCfitSigma) Loading %s',filename));
    load(filename);
    return
  end
end

% get which CRF we are using for this
CRF = getCRFFromD(d,crftype);

% get the pedestals which we are testing
pedestals = d.pedestalContrasts;
nPedestals = length(pedestals);

% some constats, to make things easier to read
DIST = 1;
ATT = 2;

% starting parameters
sigma = [0.5 0.00001];
sigma = [5 0.0000001];
%sigma = [0.03 0.01];
%sigma = [0.05 0.01];
epsilon = 0.001;
targetPCorrect = 0.7605;
dispFig = 0;

% give up bisecting after we reach some very small difference in sigma
maxBisections = round(log2(diff(sigma)/0.000001));

% cycle over pedestals, intializing min max of sigma
disppercent(-inf,sprintf('(crfaaMaxPool:CRF2TvCfitSigma) Computing model %s for extreme sigma values: %s',poolingRuleNum2Str(poolingRule,tau),mynum2str(sigma,'sigfigs=5')));
CRF.deltaC = d.distributed.TvC(4);

% test for stability and speed
if 0
  pedNum = 4;
  [responseDistributedInt1 responseDistributedInt2] = getResponses('distributed',d,pedNum,CRF,poolingRule,d.distributed.TvC(pedNum));
  disppercent(-inf,'testing')
  for i = 1:100
    x(i) = computeModelPcorrect(responseDistributedInt1,responseDistributedInt2,sigma(1),poolingRule,tau,dispFig);
    disppercent(i/100);
  end
  disppercent(inf);
  std(x)
end

% pedestals to compute model over
global crossValidate;
% for crossValidation, we fit parameters on only odd pedestals and test on the even pedestals
if crossValidate
  nPedTo = nPedestals/2;
  allPedNum = 1:2:nPedestals;
else
  nPedTo = nPedestals;
  allPedNum = 1:nPedestals;
end

for pedNumTo = 1:nPedTo
  pedNumFrom = allPedNum(pedNumTo);
  % set the deltaC for distributed and go grab the responses for the distributed trials
  [responseDistributedInt1 responseDistributedInt2] = getResponses('distributed',d,pedNumFrom,CRF,poolingRule,d.distributed.TvC(pedNumFrom));
  % set the deltaC for attended and go grab the responses for the attended trials
  [responseAttendedInt1 responseAttendedInt2] = getResponses('attended',d,pedNumFrom,CRF,poolingRule,d.attended.TvC(pedNumFrom));
  % and compute the model (use bisection search to find the sigma that fits. Start
  % by testing the extreme values of the interval. P has dimension 1:nPedestals, 1:2 (2 attention conditions, attended and distributed), 1:2 (2 sigma values)
  pDIST1(pedNumTo) = computeModelPcorrect(responseDistributedInt1,responseDistributedInt2,sigma(1),poolingRule,tau,dispFig);
  pATT1(pedNumTo) = computeModelPcorrect(responseAttendedInt1,responseAttendedInt2,sigma(1),poolingRule,tau,dispFig);
  pDIST2(pedNumTo) = computeModelPcorrect(responseDistributedInt1,responseDistributedInt2,sigma(2),poolingRule,tau,dispFig);
  pATT2(pedNumTo) = computeModelPcorrect(responseAttendedInt1,responseAttendedInt2,sigma(2),poolingRule,tau,dispFig);
%  disppercent(find(pedNum==pedestals)/length(pedestals));
end

% reconstruct p-array (for doesn't allow index and static array adressing at the same time.
p(:,ATT,1) = pATT1;
p(:,DIST,1) = pDIST1;
p(:,ATT,2) = pATT2;
p(:,DIST,2) = pDIST2;
disppercent(inf);


disp(sprintf('++ %s pedestals: (%s) ++',poolingRuleNum2Str(poolingRule,tau),mynum2str(pedestals,'sigfigs=0')));
% now maintain separate sigmas for distributed and attended conditions
sigma = repmat(sigma,[2 1]);
nBisections = 0;
pmiddist = inf;pmidatt = inf;

while (median(abs([pmiddist(:) ; pmidatt(:)]-targetPCorrect))>epsilon) && (nBisections < maxBisections)
  disppercent(-inf,sprintf('Recomputing model for sigma distributed: %0.4f sigma attended: %0.4f',mean(sigma(DIST,:)),mean(sigma(ATT,:))));
  % get mid sigma values
  midSigmaDist = mean(sigma(DIST,:));
  midSigmaAtt = mean(sigma(ATT,:));
  % can use parfor here || parallel
  parfor pedNumTo = 1:nPedTo
    pedNumFrom = allPedNum(pedNumTo);
    % set the deltaC for distributed and go grab the responses for the distributed trials
    [responseDistributedInt1 responseDistributedInt2] = getResponses('distributed',d,pedNumFrom,CRF,poolingRule,d.distributed.TvC(pedNumFrom));
    % calculate the p from the model given the mid sigma point
    pmiddist(pedNumTo) = computeModelPcorrect(responseDistributedInt1,responseDistributedInt2,midSigmaDist,poolingRule,tau,dispFig);
    % set the deltaC for attended and go grab the responses for the attended trials
    [responseAttendedInt1 responseAttendedInt2] = getResponses('attended',d,pedNumFrom,CRF,poolingRule,d.attended.TvC(pedNumFrom));
    % calculate the p from the model given the mid sigma point
    pmidatt(pedNumTo) = computeModelPcorrect(responseAttendedInt1,responseAttendedInt2,midSigmaAtt,poolingRule,tau,dispFig);
  end
  disppercent(inf);
  if fitSingleSigma
    % decide which interval the calculated p is on and set up 
    % the interval accordingly. Fit one sigma across both distributed/attended
%    if sqrt(sum((([pmiddist(:) ; pmidatt(:)]-targetPCorrect).^2).*sign([pmiddist(:) ; pmidatt(:)]-targetPCorrect)))>0
    if sqrt(sum((([pmidatt(:)]-targetPCorrect).^2).*sign([pmidatt(:)]-targetPCorrect)))>0
      p(1:nPedTo,DIST,2) = pmiddist(1:nPedTo);
      p(1:nPedTo,ATT,2) = pmidatt(1:nPedTo);
      sigma(DIST,2) = midSigmaAtt*(1+snrImprovementWithAttention);
      sigma(ATT,2) = midSigmaAtt;
    else
      p(1:nPedTo,DIST,1) = pmiddist(1:nPedTo);
      p(1:nPedTo,ATT,1) = pmidatt(1:nPedTo);
      sigma(DIST,1) = midSigmaAtt*(1+snrImprovementWithAttention);
      sigma(ATT,1) = midSigmaAtt;
    end
  else
    % decide which interval the calculated p is on and set up 
    % the interval accordingly. First for the distributed
    % trials. Select the interval according to which way most of the intervals
    % want to be set.
    if sqrt(sum(((pmiddist(:)-targetPCorrect).^2).*sign(pmiddist(:)-targetPCorrect))) > 0
      p(1:nPedTo,DIST,2) = pmiddist(1:nPedTo);
      sigma(DIST,2) = midSigmaDist;
    else
      p(1:nPedTo,DIST,1) = pmiddist(1:nPedTo);
      sigma(DIST,1) = midSigmaDist;
    end
    % now update the interval for the attended condition
    if sqrt(sum(((pmidatt(1:nPedTo)-targetPCorrect).^2).*sign(pmidatt(1:nPedTo)-targetPCorrect))) > 0
      p(1:nPedTo,ATT,2) = pmidatt(1:nPedTo);
      sigma(ATT,2) = midSigmaAtt;
    else
      p(1:nPedTo,ATT,1) = pmidatt(1:nPedTo);
      sigma(ATT,1) = midSigmaAtt;
    end
  end
  % update number of bisections done
  nBisections = nBisections + 1;
  % and display what we have so far
  disp(sprintf('%i: Distributed sigma=[%s %s] pmin=[%s] pmax=[%s]',nBisections,mynum2str(sigma(DIST,1),'sigfigs=4'),mynum2str(sigma(DIST,2),'sigfigs=4'),mynum2str(p(1:nPedTo,DIST,1)),mynum2str(p(1:nPedTo,DIST,2))));
  disp(sprintf('%i: Attended sigma=[%s %s] pmin=[%s] pmax=[%s]',nBisections,mynum2str(sigma(ATT,1),'sigfigs=4'),mynum2str(sigma(ATT,2),'sigfigs=4'),mynum2str(p(1:nPedTo,ATT,1)),mynum2str(p(1:nPedTo,ATT,2))));
end

% set the sigma to the last midpoint calculated
sigmaDistributed = midSigmaDist;
sigmaAttended = midSigmaAtt;
pDistributed = pmiddist(:);
pAttended = pmidatt(:);

if ~ieNotDefined('dataFilenames')
  eval(sprintf('save %s sigmaDistributed sigmaAttended pDistributed pAttended',filename));
end

%%%%%%%%%%%%%%%%%%%
%   getCRFFromD   %
%%%%%%%%%%%%%%%%%%%
function CRF = getCRFFromD(d,crftype)

% this sets up the CRF variable for extracting CRFs. Note that the way this is all set up
% you provided two arrays each for the attended and unattended condition. The responses
% and for what contrast those responses are for. The getResponse subroutine then pulls
% out any particular response by linear interpolation. So, if you want to use a functional
% form for a CRF, that is OK. Just evaluate that functional form very finely (say c = 0:0.001:1,
% r = functtionalForm(c,parameters)) and then stick c into the attendedContrast (or unattendedContrast) and r into
% attendedResponse (or unattendedResponse).

% is also sets the distributedContrast / distributedRespose this is only needed for the
% targetLoc model (i.e. simple old model).

% set this field so getResponse knows what to do with this structure
CRF.iscrf = 1;

% note, distributedResponse/distributedContrast is only used for the targetLoc poolingRule.

% TYPE 0: Uses the raw attended and the raw unattended
if crftype == 0
  disp(sprintf('(crfaaMaxPool:getCRFFromD) Using raw for attended and unattended'));
  CRF.attendedContrast = d.attended.pedestalsCRF;
  CRF.attendedResponse = d.attended.CRF;
  CRF.unattendedContrast = d.unattended.pedestalsCRF;
  CRF.unattendedResponse = d.unattended.CRF;
  CRF.distributedContrast = d.distributed.pedestalsCRF;
  CRF.distributedResponse = d.distributed.CRFtarget;
  CRF.distributedNontargetResponse = d.distributed.CRFnontarget;
  % TYPE 1: Uses the model fit for the attended and the raw unattended
elseif crftype == 1
  disp(sprintf('(crfaaMaxPool:getCRFFromD) Using model fit for attended and raw for unattended'));
  CRF.attendedContrast = d.attended.modelContrast;
  CRF.attendedResponse = d.attended.modelResponse;
  CRF.unattendedContrast = d.unattended.pedestalsCRF;
  CRF.unattendedResponse = d.unattended.CRF;
  CRF.distributedContrast = d.distributed.modelContrast;
  CRF.distributedResponse = d.distributed.modelResponse;
% TYPE 2: Uses the model fit for the attended and the the model fit for the
% distributed offset to fit the unattended CRF
elseif crftype == 2
  disp(sprintf('(crfaaMaxPool:getCRFFromD) Using model fit for attended and unattended (i.e. model for distributed is offset to fit unattended data)'));
  CRF.attendedContrast = d.attended.modelContrast;
  CRF.attendedResponse = d.attended.modelResponse;
  % fit distributed model to unattended responses
  % start by evaluating the model at the pedestal contrast
  unattendedContrast = d.unattendedContrast;
  unattendedResponse = interpolateResponse(d.distributed.modelContrast,d.distributed.modelResponse,unattendedContrast);
  % now compute offset
  offset = fitOffset(unattendedResponse,d.unattendedResponse);
  % now we set the CRF to be the modelResponse for distributed offset by the above
  CRF.unattendedContrast = unattendedContrast;
  CRF.unattendedResponse = unattendedResponse - offset;
  if 0
    plot(d.unattendedContrast,d.unattendedResponse,'ko');hold on
    plot(CRF.unattendedContrast,CRF.unattendedResponse,'k-');
    plot(CRF.unattendedContrast,CRF.unattendedResponse+offset,'g-');
    xlabel('contrast');
    ylabel('response');
    legend('unattended CRF','Fit','Fit w/out offset');
  end
  % use the model for the distributed
  CRF.distributedContrast = d.distributed.modelContrast;
  CRF.distributedResponse = d.distributed.modelResponse;
% use average of attended and distributed model fits to fit the attended and unattended curves
elseif crftype == 3
  disp(sprintf('(crfaaMaxPool:getCRFFromD) Using average attended/distributed model fit for attended and unattended'));
  % recmpute model for known contrasts
  modelContrast = 0:0.01:1;
  attendedModelResponse = interpolateResponse(d.attended.modelContrast,d.attended.modelResponse,modelContrast);
  distributedModelResponse = interpolateResponse(d.distributed.modelContrast,d.distributed.modelResponse,modelContrast);
  % subtract off means
  attendedModelResponse = attendedModelResponse-mean(attendedModelResponse);
  distributedModelResponse = distributedModelResponse-mean(distributedModelResponse);
  % and calculate mean.
  meanModelResponse = (attendedModelResponse + distributedModelResponse)/2;
  % get actual attended response at the modelContrast
  meanModelAttendedResponse = interpolateResponse(modelContrast,meanModelResponse,d.attendedContrast);
  % now compute offset
  offset = fitOffset(meanModelAttendedResponse,d.attendedResponse);
  CRF.attendedContrast = modelContrast;
  CRF.attendedResponse = meanModelResponse-offset;
  % fit mean model to unattended responses
  % start by evaluating the model at the pedestal contrast
  unattendedContrast = d.unattendedContrast;
  meanModelUnattendedResponse = interpolateResponse(modelContrast,meanModelResponse,unattendedContrast);
  % now compute offset
  offset = fitOffset(meanModelUnattendedResponse,d.unattendedResponse);
  % now we set the CRF to be the modelResponse, offset by what we calculate here.
  CRF.unattendedContrast = modelContrast;
  CRF.unattendedResponse = meanModelResponse - offset;
  if 0
    plot(d.unattendedContrast,d.unattendedResponse,'ko');hold on
    plot(d.attendedContrast,d.attendedResponse,'ro');
    plot(CRF.unattendedContrast,CRF.unattendedResponse,'k-');
    plot(CRF.unattendedContrast,CRF.unattendedResponse+offset,'g-');
    plot(CRF.attendedContrast,CRF.attendedResponse,'k-');
    plot(CRF.attendedContrast,CRF.attendedResponse+offset,'g-');
    xlabel('contrast');
    ylabel('response');
    legend('unattended CRF','attended CRF','Fit','Fit w/out offset');
  end
  % fit to distribued target responses
  distributedContrast = d.attendedContrast; % contrast are same as focal-cue case
  meanModelDistributedTargetResponse = interpolateResponse(modelContrast,meanModelResponse,distributedContrast);
  % now compute offset
  offset = fitOffset(meanModelDistributedTargetResponse,d.distributed.CRFtarget);
  % now we set the CRF to be the modelResponse, offset by what we calculate here.
  CRF.distributedContrast = modelContrast;
  CRF.distributedResponse = meanModelResponse - offset;

  % fit to distribued non-target responses
  distributedContrast = d.unattendedContrast; % contrast are same as focal-cue case
  meanModelDistributedTargetResponse = interpolateResponse(modelContrast,meanModelResponse,distributedContrast);
  % now compute offset
  offset = fitOffset(meanModelDistributedTargetResponse,d.distributed.CRFnontarget);
  % now we set the CRF to be the modelResponse, offset by what we calculate here.
  CRF.distributedNontargetResponse = meanModelResponse - offset;
  % use the model for the distributed
  %CRF.distributedContrast = d.distributed.modelContrast;
  %CRF.distributedResponse = d.distributed.modelResponse;
% get crf from fits
elseif crftype == 4
  % unattended CRF  
  CRF.unattendedContrast = d.unattended.fitContrast;
  CRF.unattendedResponse = d.unattended.fitResponse;
  % attended
  CRF.attendedContrast = d.attended.fitContrast;
  CRF.attendedResponse = d.attended.fitResponse;
  % distributed
  CRF.distributedContrast = d.distributed.fitContrast;
  CRF.distributedResponse = d.distributed.fitResponse;
  % distributed
  CRF.distributedContrast = d.distributed.fitContrast;
  CRF.distributedResponse = d.distributed.fitResponse;
% TYPE 5: Uses the model fit for the attended and the the model fit for the
% distributed offset to fit the unattended CRF
elseif crftype == 5
  % recmpute model for known contrasts
  modelContrast = 0:0.01:1;
  attendedModelResponse = interpolateResponse(d.attended.modelContrast,d.attended.modelResponse,modelContrast);
  distributedModelResponse = interpolateResponse(d.distributed.modelContrast,d.distributed.modelResponse,modelContrast);
  CRF.attendedContrast = modelContrast;
  CRF.attendedResponse = attendedModelResponse;
  % subtract off means
  attendedModelResponse = attendedModelResponse-mean(attendedModelResponse);
  distributedModelResponse = distributedModelResponse-mean(distributedModelResponse);
  % and calculate mean.
  meanModelResponse = (attendedModelResponse + distributedModelResponse)/2;
  % fit attendedModelResponse to unattended
  unattendedContrast = d.unattendedContrast;
  modelUnattendedResponse = interpolateResponse(modelContrast,attendedModelResponse,unattendedContrast);
  % now compute offset
  offset = fitOffset(modelUnattendedResponse,d.unattendedResponse);
  % now we set the CRF to be the modelResponse, offset by what we calculate here.
  CRF.unattendedContrast = modelContrast;
  CRF.unattendedResponse = attendedModelResponse - offset;

  % fit to distribued target responses
  distributedContrast = d.attendedContrast; % contrast are same as focal-cue case
  modelDistributedTargetResponse = interpolateResponse(modelContrast,distributedModelResponse,distributedContrast);
  % now compute offset
  offset = fitOffset(modelDistributedTargetResponse,d.distributed.CRFtarget);
  % now we set the CRF to be the modelResponse, offset by what we calculate here.
  CRF.distributedContrast = modelContrast;
  CRF.distributedResponse = distributedModelResponse - offset;

  % fit to distribued non-target responses
  distributedContrast = d.unattendedContrast; % contrast are same as focal-cue case
  modelDistributedTargetResponse = interpolateResponse(modelContrast,distributedModelResponse,distributedContrast);
  % now compute offset
  offset = fitOffset(modelDistributedTargetResponse,d.distributed.CRFnontarget);
  % now we set the CRF to be the modelResponse, offset by what we calculate here.
  CRF.distributedNontargetResponse = distributedModelResponse - offset;
  
end


%%%%%%%%%%%%%%%%%%%%%%%%
%   CRF2TvCfitDeltaC   %
%%%%%%%%%%%%%%%%%%%%%%%%
function dout = CRF2TvCfitDeltaC(d,poolingRule,tau,sigma,crftype,dataFilenames,recalc)

filename = fullfile(dataFilenames.saveDataDir,fixBadChars(sprintf('TvCFitfor%s_%i_o%i_%s_%s',poolingRuleNum2Str(poolingRule,tau),crftype,d.observerNum,d.visualArea,mynum2str(sigma,'sigfigs=5','doFixBadChars=1'))));

if ieNotDefined('recalc'),recalc = 0;end

% preload
if isfile(sprintf('%s.mat',filename)) && ~recalc
  disp(sprintf('(crfaaMaxPool:CRF2TvCfitDeltaC) Loading %s',filename));
  load(filename);
  return
end

% get which CRF we are using for this
CRF = getCRFFromD(d,crftype);

% get the pedestals which we are testing
pedestals = d.pedestalContrasts;
nPedestals = length(pedestals);

% some constats, to make things easier to read
DIST = 1;
ATT  = 2;

% starting parameters. DeltaC's get started with min/max set around the
% actual psychophysical measurements
deltaCdistributed = [d.distributed.TvC'/10 ; d.distributed.TvC'*10]';
deltaCattended = [d.attended.TvC'/10 ; d.attended.TvC'*10]';
epsilon = 0.001;
targetPCorrect = 0.7605;
dispFig = 0;

% give up bisecting after we reach some very small difference in deltaC
maxBisections = round(log2(max(diff(deltaCdistributed'))/0.001));

%deltaCdistributed = [d.distributed.TvC'/2 ; d.distributed.TvC'*2]';
%deltaCattended = [d.attended.TvC'/2 ; d.attended.TvC'*2]';
%maxBisections = 5;
disp(sprintf('(crfaaMaxPool:CRF2TvCfitDeltaC) Max bisections: %i',maxBisections));


% if we only have one sigma value then use it for both curves
if length(sigma) == 1
  sigma(2) = sigma;
end

% pedestals to compute model over
pedestals = 1:nPedestals;
nPedTo = nPedestals;
allPedNum = 1:nPedestals;


% cycle over pedestals, intializing min max of deltaC
disppercent(-inf,sprintf('(crfaaMaxPool:CRF2TvCfitDeltaC) Calculating model %s on extreme values',poolingRuleNum2Str(poolingRule,tau)));
% can use parfor here (though it's not working well for me). || parfor
parfor pedNumTo = 1:nPedTo
  pedNumFrom = allPedNum(pedNumTo);
  % set the deltaC for distributed using the minimum deltaC
  % and go grab the responses for the distributed trials
  [responseDistributedInt1 responseDistributedInt2] = getResponses('distributed',d,pedNumFrom,CRF,poolingRule,deltaCdistributed(pedNumFrom,1));
  % set the deltaC for attended and go grab the responses for the attended trials
  [responseAttendedInt1 responseAttendedInt2] = getResponses('attended',d,pedNumFrom,CRF,poolingRule,deltaCattended(pedNumFrom,1));
  % and compute the model (use bisection search to find the sigma that fits. Start
  % by testing the extreme values of the interval. P has dimension 1:nPedestals, 1:2 (2 attention conditions, attended and distributed), 1:2 (2 sigma values)
  pDIST1(pedNumTo) = computeModelPcorrect(responseDistributedInt1,responseDistributedInt2,sigma(1),poolingRule,tau,dispFig);
  pATT1(pedNumTo) = computeModelPcorrect(responseAttendedInt1,responseAttendedInt2,sigma(2),poolingRule,tau,dispFig);

  % set the deltaC for distributed using the maximum deltaC
  % and go grab the responses for the distributed trials
  [responseDistributedInt1 responseDistributedInt2] = getResponses('distributed',d,pedNumFrom,CRF,poolingRule,deltaCdistributed(pedNumFrom,2));
  % set the deltaC for attended and go grab the responses for the attended trials
  [responseAttendedInt1 responseAttendedInt2] = getResponses('attended',d,pedNumFrom,CRF,poolingRule,deltaCattended(pedNumFrom,2));
  pDIST2(pedNumTo) = computeModelPcorrect(responseDistributedInt1,responseDistributedInt2,sigma(2),poolingRule,tau,dispFig);
  pATT2(pedNumTo) = computeModelPcorrect(responseAttendedInt1,responseAttendedInt2,sigma(2),poolingRule,tau,dispFig);

%  disppercent(find(pedNum==pedestals)/length(pedestals));
end
p(:,ATT,1) = pATT1;
p(:,DIST,1) = pDIST1;
p(:,ATT,2) = pATT2;
p(:,DIST,2) = pDIST2;
disppercent(inf);

disp(sprintf('(crfaaMaxPool:CRF2TvCfitDeltaC) ++ %s pedestals: (%s) sigma=%s ++',poolingRuleNum2Str(poolingRule,tau),mynum2str(pedestals,'sigfigs=0'),mynum2str(sigma)));
nBisections = 0;
pmid = inf;

while (median(abs(pmid(:)-targetPCorrect))>epsilon) && (nBisections < maxBisections)
  % get the mid deltaC
  midDeltaCattended = mean(deltaCattended,2);
  midDeltaCdistributed = mean(deltaCdistributed,2);
  % and recompute the model for that mid point
  disppercent(-inf,sprintf('(crfaaMaxPool:CRF2TvCfitDeltaC) Recomputing model for deltaC distributed: %s deltaC attended: %s',mynum2str(midDeltaCdistributed,'sigfigs=4'),mynum2str(midDeltaCattended,'sigfigs=4')));
															  % can use parfor here (though it's not working well for me). | parallel
  parfor pedNumTo = 1:nPedTo
    pedNumFrom = allPedNum(pedNumTo);
    % set the deltaC for distributed and go grab the responses for the distributed trials
    [responseDistributedInt1 responseDistributedInt2] = getResponses('distributed',d,pedNumFrom,CRF,poolingRule,midDeltaCdistributed(pedNumFrom));
    % calculate the p from the model given the mid sigma point
    pmiddist(pedNumTo) = computeModelPcorrect(responseDistributedInt1,responseDistributedInt2,sigma(1),poolingRule,tau,dispFig);
    % set the deltaC for attended and go grab the responses for the attended trials
    [responseAttendedInt1 responseAttendedInt2] = getResponses('attended',d,pedNumFrom,CRF,poolingRule,midDeltaCattended(pedNumFrom));
    % calculate the p from the model given the mid deltaC point
    pmidatt(pedNumTo) = computeModelPcorrect(responseAttendedInt1,responseAttendedInt2,sigma(2),poolingRule,tau,dispFig);
%    disppercent(find(pedNum==pedestals)/length(pedestals));
  end
  disppercent(inf);
  % decide which interval the calculated p is on and set up 
  % the interval accordingly. First for the distributed
  % trials. Select the interval according to which way most of the intervals
  % want to be set.
  for pedNum = 1:nPedTo
    if pmiddist(pedNum)>targetPCorrect
      p(pedNum,DIST,2) = pmiddist(pedNum);
      deltaCdistributed(pedNum,2) = midDeltaCdistributed(pedNum);
    else
      p(pedNum,DIST,1) = pmiddist(pedNum);
      deltaCdistributed(pedNum,1) = midDeltaCdistributed(pedNum);
    end
  end
  % now update the interval for the attended condition
  for pedNum = 1:nPedTo
    if pmidatt(pedNum)>targetPCorrect
      p(pedNum,ATT,2) = pmidatt(pedNum);
      deltaCattended(pedNum,2) = midDeltaCattended(pedNum);
    else
      p(pedNum,ATT,1) = pmidatt(pedNum);
      deltaCattended(pedNum,1) = midDeltaCattended(pedNum);
    end
  end
  % update number of bisections done
  nBisections = nBisections + 1;
  % and display what we have so far
  disp(sprintf('%i: Distributed deltaC=[%s] p=[%s]',nBisections,mynum2str(deltaCdistributed,'sigfigs=4'),mynum2str(pmiddist(1:nPedTo))));
  disp(sprintf('%i: Attended deltaC=[%s] p=[%s]',nBisections,mynum2str(deltaCattended,'sigfigs=4'),mynum2str(pmidatt(1:nPedTo))));

end

% since we met criterion on the last bisection, return the best fitting deltaC
dout.deltaCattended = midDeltaCattended;
dout.deltaCdistributed = midDeltaCdistributed;

% compute % variance explained
dout.distributedResidual = d.distributed.TvC-dout.deltaCdistributed;
dout.distributedResidualSS = sum(dout.distributedResidual.^2);
dout.distributedDataSS = sum((d.distributed.TvC-mean(d.distributed.TvC)).^2);
dout.distributedr2 = 1-dout.distributedResidualSS/dout.distributedDataSS;
dout.attendedResidual = d.attended.TvC-dout.deltaCattended;
dout.attendedResidualSS = sum(dout.attendedResidual.^2);
dout.attendedDataSS = sum((d.attended.TvC-mean(d.attended.TvC)).^2);
dout.attendedr2 = 1-dout.attendedResidualSS/dout.attendedDataSS;
dout.overallResidualSS = dout.distributedResidualSS+dout.attendedResidualSS;
dout.overalDataSS = dout.distributedDataSS+dout.attendedDataSS;
dout.overallr2 = 1- dout.overallResidualSS/dout.overalDataSS;

% new calculation of variance explaines. Based on all of the TvC Data
dout.total.data = [d.attended.TvC(:)' d.distributed.TvC(:)'];
dout.total.SSdata = sum((dout.total.data-mean(dout.total.data)).^2);
dout.total.model = [dout.deltaCattended(:)' dout.deltaCdistributed(:)'];
dout.total.SSmodel = sum((dout.total.model-mean(dout.total.model)).^2);
dout.total.residual = dout.total.data-dout.total.model;
dout.total.SSres = sum((dout.total.residual-mean(dout.total.residual)).^2);
dout.total.r2 = 1 - dout.total.SSres/dout.total.SSdata;

global crossValidate;
% compute cross validated r2
if crossValidate
  allPedNum = 2:2:nPedestals;
  dout.crossValidate.data = [d.attended.TvC(allPedNum)' d.distributed.TvC(allPedNum)'];
  dout.crossValidate.SSdata = sum((dout.crossValidate.data-mean(dout.crossValidate.data)).^2);
  dout.crossValidate.model = [dout.deltaCattended(allPedNum)' dout.deltaCdistributed(allPedNum)'];
  dout.crossValidate.SSmodel = sum((dout.crossValidate.model-mean(dout.crossValidate.model)).^2);
  dout.crossValidate.residual = dout.crossValidate.data-dout.crossValidate.model;
  dout.crossValidate.SSres = sum((dout.crossValidate.residual-mean(dout.crossValidate.residual)).^2);
  dout.crossValidate.r2 = 1 - dout.crossValidate.SSres/dout.crossValidate.SSdata;
end

n = nPedTo * 2;
if poolingRule == 3
%  k = 2; (need to account for sigma)
  k = 3;
else
%  k = 1;
  k = 2;
end
dout.overallBIC = n *log (dout.overallResidualSS/n) + k*log(n);
dout.overallAIC = 2*k + n * log (dout.overallResidualSS/n) + 2*k*(k+1)/(n-k-1);
%dout.overallAIC = log((1-dout.overallr2)/n)+2*k;

% compute goodness-of-fit
chi2 = sum(((dout.distributedResidual).^2)./(d.distributed.TvCste.^2));
dof = length(dout.distributedResidual)-2;
chi2pdf(chi2,dof);
chi2 = sum(((dout.attendedResidual).^2)./(d.attended.TvCste.^2));
chi2pdf(chi2,dof);
chi2 = sum([dout.attendedResidual'.^2 dout.distributedResidual'.^2]./[d.attended.TvCste'.^2 d.distributed.TvCste'.^2]);
dof = n-k;
chi2pdf(chi2,dof);

% set model name for printing
dout.modelName = sprintf('%s crftype=%i',poolingRuleNum2Str(poolingRule,tau),crftype);

% and parameters
dout.sigma = sigma;
dout.tau = tau;
dout.poolingRule = poolingRule;
dout.crftype = crftype;

dispFig = 0;
if dispFig
  dispModelFitToTvC(d,dout); % new call: dispModelFitToTvC(d,dout,dataFolder,visArea)
end

% save file
eval(sprintf('save %s dout',filename));

%%%%%%%%%%%%%%%%%%%%%%%%%
%   dispModelFitToTvC   %
%%%%%%%%%%%%%%%%%%%%%%%%%
function dispModelFitToTvC(d,dout,dataFolder,visArea)

% set display contrasts so they can plot on the log axis
displayContrasts = d.pedestalContrasts;
displayContrasts(1) = 0.01;

% open figure
figname = sprintf('crfaaMaxPoolTvCFit_V%i_poolingRule%i',visArea,dout.poolingRule);
h = smartfig(figname,'reuse');
clf;
% to do here:
% (1) add errorbars for TvCs
% (2) make a smooth curve (maybe by interpolation)
% (3) change colors, format, add title/labels

% display distributed curve and fit
loglog(displayContrasts,d.distributed.TvC,'bo','MarkerFaceColor','b');
hold on
loglog(displayContrasts,d.attended.TvC,'ro','MarkerFaceColor','r');
loglog(displayContrasts,dout.deltaCdistributed,'b-');
loglog(displayContrasts,dout.deltaCattended,'r-');
xlabel('Contrast (%)');
ylabel('Threshold contrast (%)');
axis('square')

if isfield(dout,'crossValidate')
  title(sprintf('%s sigma: %s crftype: %i\ncrossValidate: r^2=%0.2f total: r^2=%0.2f AIC:%3f BIC:%0.3f',dout.modelName,mynum2str(dout.sigma,'sigfigs=4'),dout.crftype,dout.crossValidate.r2*100,dout.total.r2*100,dout.overallAIC,dout.overallBIC));
elseif isfield(dout,'total')
  title(sprintf('%s sigma: %s crftype: %i\ntotal: r^2=%0.2f AIC:%3f BIC:%0.3f',dout.modelName,mynum2str(dout.sigma,'sigfigs=4'),dout.crftype,dout.total.r2*100,dout.overallAIC,dout.overallBIC));
else
  title(sprintf('%s sigma: %s crftype: %i\noverall: r^2=%0.2f AIC:%3f BIC:%0.3f',dout.modelName,mynum2str(dout.sigma,'sigfigs=4'),dout.crftype,dout.overallr2*100,dout.overallAIC,dout.overallBIC));
end
drawnow

savefigure = 1;
if savefigure
 savefig(h,'verbose',1,'figName',figname,'figDir','fig_tvcfit','defaultDataFolder',dataFolder.saveDataDir);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   makeDistributionPlot   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function makeDistributionPlot(d,poolingRule,tau,sigma,crftype,pedNum,snrImprovementWithAttention)

if ieNotDefined('pedNum'),pedNum = 3;end
if ieNotDefined('snrImprovementWithAttention'),snrImprovementWithAttention=0;end

% get which CRF we are using for this
CRF = getCRFFromD(d,crftype);

% use same sigma for attended and distributed if only one passed in.
if length(sigma) == 1
  disp(sprintf('(crfaaMaxPool:makeDistributionPlot) Setting sigma(1) be %f times simga(2)',1+snrImprovementWithAttention));
  sigma(2) = sigma(1);
  sigma(1) = (1+snrImprovementWithAttention)*sigma(1);
end

doDistributed = 1;doAttended = 0;

if doDistributed
  % get distributions for distributed
  CRF.deltaC = d.distributed.TvC(pedNum);
  [responseDistributedInt1 responseDistributedInt2] = getResponses('distributed',d,pedNum,CRF,poolingRule);
  [p ddist] = computeModelPcorrect(responseDistributedInt1,responseDistributedInt2,sigma(1),poolingRule,tau,0);
  % plot figure
  for i = 1:size(ddist.responseDistInt1,2)
    smartfig(sprintf('crfaaMaxPoolAttDist%i_K',i),'reuse');clf;
    displayResponseAndPoolingFigure(ddist.responseDistInt1(:,i,:),ddist.responseDistInt2(:,i,:),ddist.pooledResponseDistInt1(i,:),ddist.pooledResponseDistInt2(i,:),ddist.hits(i,:),ddist.falseAlarms(i,:),ddist.areaUnderROC(i));
    currentTitle = get(get(gca,'Title'),'String');
    title(sprintf('%s\nDistributed %i: %s\nsigma: %0.4f crftype: %i',deblank(strtrim(currentTitle')),i,poolingRuleNum2Str(poolingRule,tau),sigma(1),crftype));
  end
end

if doAttended
  % get distributions for attended
  CRF.deltaC = d.attended.TvC(pedNum);
  [responseAttendedInt1 responseAttendedInt2] = getResponses('attended',d,pedNum,CRF,poolingRule);
  [p datt] = computeModelPcorrect(responseAttendedInt1,responseAttendedInt2,sigma(2),poolingRule,tau,0);
  smartfig('crfaaMaxPoolAttDist','reuse');clf;
  displayResponseAndPoolingFigure(datt.responseDistInt1,datt.responseDistInt2,datt.pooledResponseDistInt1,datt.pooledResponseDistInt2,datt.hits,datt.falseAlarms,datt.areaUnderROC);
  currentTitle = get(get(gca,'Title'),'String');
  title(sprintf('%s\nAttended %s\nsigma: %0.4f crftype: %i',deblank(strtrim(currentTitle')),poolingRuleNum2Str(poolingRule,tau),sigma(2),crftype));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   displayresponseAndPoolingFigure   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function displayResponseAndPoolingFigure(responseDistInt1,responseDistInt2,pooledResponseDistInt1,pooledResponseDistInt2,hits,falseAlarms,areaUnderROC)

% get min and max of axis
minx = 0;  %floor(100*min([responseDistInt1(:)' responseDistInt2(:)' pooledResponseDistInt1(:)' pooledResponseDistInt2(:)']))/100;
maxx = .75;%ceil(100*max([responseDistInt1(:)' responseDistInt2(:)' pooledResponseDistInt1(:)' pooledResponseDistInt2(:)']))/100;

% plotnums
plotnums = {[1 2 9 10],[3 4 11 12],[17 18 25 26],[19 20 27 28]};
stimulusLocation = [4 1 3 2]; % this is to match the clock-wise cycling of locations

% set up the general infor for the plot
% this is consistent across figures in the attention paper
plotInfo = setUpPlotInfo;
Int1useThisColor = plotInfo.MarkerFaceColor{3}{1};
Int2useThisColor = plotInfo.MarkerFaceColor{3}{2};
rocUseThisColor  = plotInfo.MarkerFaceColor{2}{1};

% display distributions for each stimulus quandrant
for stimNum = 1:4
  subplot(4,8,plotnums{stimNum});
  nbars = 20;
  myhist(squeeze(responseDistInt1(stimNum,:,:)),nbars,Int1useThisColor,1);
  hold on
  myhist(squeeze(responseDistInt2(stimNum,:,:)),nbars,Int2useThisColor,1);
  title(sprintf('Stimulus location %i',stimulusLocation(stimNum)));
  set(gca,...
     'FontName','Helvetica','FontSize',plotInfo.Fsize, ...
     'FontAngle','italic', ...
     'PlotBoxAspectRatio',plotInfo.PlotBoxAspectRatio, ...
     'LineWidth',plotInfo.LineWidth,'TickLength',plotInfo.TickLength, ...
     'YTick', [0 .25],        'YTickLabel', [0 0.25] ,...
     'XTick', [minx:0.25:maxx], 'XTickLabel', [minx:.15:maxx] ,...
     'XLim',[minx maxx], 'YLim', [0 .25], ...
     'XColor',plotInfo.XYColor,'YColor',plotInfo.XYColor,'Box','off');
  myaxishist;
end

% display distributions after pooling
subplot(4,8,[6 7 14 15]);
myhist(pooledResponseDistInt1,nbars,Int1useThisColor,1);
hold on
myhist(pooledResponseDistInt2,nbars,Int2useThisColor,1);
xaxis(minx,maxx);
yaxis(0,.25);
set(gca,...
 'FontName','Helvetica','FontSize',plotInfo.Fsize, ...
 'FontAngle','Italic', ...
 'PlotBoxAspectRatio',plotInfo.PlotBoxAspectRatio, ...
 'LineWidth',plotInfo.LineWidth,'TickLength',plotInfo.TickLength, ...
 'YTick', [0 .25], 'YTickLabel', [0 0.25] ,...
 'XTick', [minx:0.25:maxx], 'XTickLabel', [minx:.15:maxx] ,...
 'XColor',plotInfo.XYColor,'YColor',plotInfo.XYColor,'Box','off');
myaxishist;
  
subplot(4,8,[22 30 31]);
plot(falseAlarms,hits,'Color',rocUseThisColor, 'LineWidth',plotInfo.LineWidth);
axis square
ylabel('False alarms','FontSize',plotInfo.Fsize)
xlabel('Hits','FontSize',plotInfo.Fsize);
title(sprintf('Percent Correct\n(Area under ROC = %0.3f)',areaUnderROC), ...
      'FontSize',plotInfo.Fsize);
set(gca,...
 'FontName','Helvetica','FontSize',plotInfo.Fsize, ...
 'PlotBoxAspectRatio',plotInfo.PlotBoxAspectRatio, ...
 'LineWidth',plotInfo.LineWidth,'TickLength',plotInfo.TickLength, ...
 'XTick', [0 .5 1], 'XTickLabel', [0 0.5 1] ,...
 'YTick', [0 .5 1], 'YTickLabel', [0 0.5 1] ,...
 'XColor',plotInfo.XYColor,'YColor',plotInfo.XYColor,'Box','off');
myaxisroc;
drawnow


%%%%%%%%%%%%%%%%%%%%%%%%
%   getWeightsForTau   %
%%%%%%%%%%%%%%%%%%%%%%%%
function [weightsAttended weightsDistributed] = getWeightsForTau(d,poolingRule,tau)

% get which pooling rule to use
if poolingRule == 3
  type = 'softmax';
elseif any(poolingRule == [6:11])
  type = 'exponent';
else
  disp(sprintf('(crfaaMaxPool:getWeightsForTau) Unknown pooling rule %i',poolingRule));
  keyboard
end
  
% Set pooling rule to targetLoc so that we just get regular attended
% and distributed responses
poolingRule = 4;

% deltaC and altdata aren't used for getting the weights. (i.e. we just
% get them for the pedestal contrasts, using the interpolated responses
% from the CRF).
deltaC = nan;
altdata = [];

% now go through all pedestals.
for pedestalNum = 1:d.nPedestals
  % get the responses to target and distractors for attended
  response = getResponses('attended',d,pedestalNum,altdata,poolingRule,deltaC);
  % and compute their weights with softmax
  if isequal(type,'softmax');
    weightsAttended(:,pedestalNum) = softmax(response,tau);
  elseif isequal(type,'exponent')
    [r weightsAttended(:,pedestalNum)] = exponentSum(response,tau);
  end

  % do the same for the distributed condition.
  response = getResponses('distributed',d,pedestalNum,altdata,poolingRule,deltaC);
  if isequal(type,'softmax');
    weightsDistributed(:,pedestalNum) = softmax(response,tau);
  elseif isequal(type,'exponent')
    [r weightsDistributed(:,pedestalNum)] = exponentSum(response,tau);
  end
end

% get average weights
weightsAttended = mean(weightsAttended,2);
weightsDistributed = mean(weightsDistributed,2);


%%%%%%%%%%%%%%%%%%%
%   plotWeights   %
%%%%%%%%%%%%%%%%%%%
function plotWeights(w,x,y,dispText,titleText)

border = 0.1;
fix = border/2;
width = diff(x);
height = diff(y);
borderWidth = width*border;
borderHeight = height*border;
innerWidth = (width-borderWidth*4)/2;
innerHeight = (height-borderHeight*4)/2;
x1 = x(1)+borderWidth;
x2 = x(1)+borderWidth*3+innerWidth;
y1 = y(1)+borderHeight;
y2 = y(1)+borderHeight*3+innerHeight;
rectangle('Position',[x1 y2 innerWidth innerHeight]);
hold on
rectangle('Position',[x1 y2 innerWidth innerHeight*w(1)],'Facecolor','k');
rectangle('Position',[x2 y2 innerWidth innerHeight]);
rectangle('Position',[x2 y2 innerWidth innerHeight*w(2)],'FaceColor','k');
rectangle('Position',[x2 y1 innerWidth innerHeight]);
rectangle('Position',[x2 y1 innerWidth innerHeight*w(3)],'FaceColor','k');
rectangle('Position',[x1 y1 innerWidth innerHeight]);
if dispText
  text(x1+innerWidth/2,y2+innerHeight+borderHeight/5,sprintf('%0.2f',w(1)),'HorizontalAlignment','Center','Color','k','VerticalAlignment','bottom');
  text(x2+innerWidth/2,y2+innerHeight+borderHeight/5,sprintf('%0.2f',w(2)),'HorizontalAlignment','Center','Color','k','VerticalAlignment','bottom');
  text(x2+innerWidth/2,y1+innerHeight+borderHeight/5,sprintf('%0.2f',w(3)),'HorizontalAlignment','Center','Color','k','VerticalAlignment','bottom');
  text(x1+innerWidth/2,y1+innerHeight+borderHeight/5,sprintf('%0.2f',w(4)),'HorizontalAlignment','Center','Color','k','VerticalAlignment','bottom');
end
if nargin >= 5
  text(x(1)+width/2,y(2)+borderHeight/2,titleText,'HorizontalAlignment','Center','VerticalAlignment','bottom','Color','k');
end
rectangle('Position',[x1 y1 innerWidth innerHeight*w(4)],'FaceColor','k');
rectangle('Position',[x(1) y(1) width height]);
plot([x(1)+width/2-width*fix/2 x(1)+width/2+width*fix/2],[y(1)+height/2 y(1)+height/2],'k-');
plot([x(1)+width/2 x(1)+width/2],[y(1)+height/2-height*fix/2 y(1)+height/2+height*fix/2],'k-');


%%%%%%%%%%%%%%%%%%%%%%%%%
%   calcAkaikeWeights   %
%%%%%%%%%%%%%%%%%%%%%%%%%
function [AIC AICweights] = calcAkaikeWeights(AIC)

% get delta weights
AIC = AIC-min(AIC);

%compute weights
AICweights = exp(-0.5*AIC)/sum(exp(-0.5*AIC));


%%%%%%%%%%%%%%%%%%
% setDataDefault %
%%%%%%%%%%%%%%%%%%
function dataFilenames = setDataDefault(decisionModel,crffit)
% this function sets the default data path 
% and filenames for the data files, to load 

if isdir('~/data/crfaa/crfaaMaxPool')
  dataFilenames.path = '~/data/crfaa/crfaaMaxPool/data';
  dataFilenames.saveDataDir = '~/data/crfaa/crfaaMaxPool';
elseif isdir('~/data/riken/crfaa/crfaaMaxPool')
  dataFilenames.path = '~/data/riken/crfaa/crfaaMaxPool/data';
  dataFilenames.saveDataDir = '~/data/riken/crfaa/crfaaMaxPool';
elseif isdir('~/data/riken/crfaa/fmridata/crfaa/crfaaMaxPool')
  dataFilenames.path = '~/data/riken/crfaa/fmridata/crfaa/crfaaMaxPool/data';
  dataFilenames.saveDataDir = '~/data/riken/crfaa/fmridata/crfaa/crfaaMaxPool';
else
  disp(sprintf('(crfaaMaxPool:setDataDefault) Could not find data directory'));
  keyboard
end

dataFilenames.computer = gethostname;

% files names to load
dataFilenames.decisionModel = decisionModel;
dataFilenames.crffit        = crffit;
dataFilenames.tbyt = 'tbyt';

%%%%%%%%%%%%%%
% setDataDir %
%%%%%%%%%%%%%%
function dataFilenames = setDataDir(dataDir,bootstrapNum)
% this function sets the data path  and looks for the decisionmodel data in the
% dataDir/data

% check for existence of data
if ~isdir(dataDir)
  disp(sprintf('(crfaaMaxPool) Could not find dataDir: %s',dataDir));
  keyboard
end
if ~isdir(fullfile(dataDir,'data'))
  disp(sprintf('(crfaaMaxPool) Could not find data folder in dataDir: %s/data',dataDir));
  keyboard
end

% set directories
dataFilenames.path = fullfile(dataDir,'data');
if isempty(bootstrapNum)
  dataFilenames.saveDataDir = dataDir;
else
  dataFilenames.saveDataDir = fullfile(dataDir,sprintf('bootstrap%i',bootstrapNum));
  % make the directory if necessary
  if ~isdir(dataFilenames.saveDataDir)
    mkdir(dataFilenames.saveDataDir);
  end
end
disp(sprintf('(crfaaMaxPool) SaveDataDir set to: %s',dataFilenames.saveDataDir));

dataFilenames.computer = gethostname;

% files names to load. Look for a mat file in dataDir/data. This should be the decision model
d = dir(fullfile(dataFilenames.path,'*.mat'));
if length(d) ~= 1
  disp(sprintf('(crfaaMaxPool:setDataDir) Should be only a single file mat file in %s - which should be decision model data',dataFilenames.path));
  keyboard
else
  disp(sprintf('(crfaaMaxPool:setDataDir) Using %s as datafile for model',d(1).name));
end
dataFilenames.decisionModel = d(1).name;

% for now, we don't load crffit for alternate data directories.
dataFilenames.crffit        = '';

% what is this for? jg
dataFilenames.tbyt = 'tbyt';
dataFilenames.bootstrapNum = bootstrapNum;

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
plotInfo.Fsize = 10;
plotInfo.MarkerSize = 10;
plotInfo.LineWidth = 1;
plotInfo.TickLength = [0.025 .01];
plotInfo.PlotBoxAspectRatio = [1 1 1];
plotInfo.YLim = [0 1.3];
plotInfo.Xlim = [.001,1];
plotInfo.plotPeds = [plotInfo.Xlim(1) 0.0175 0.035 0.07 0.14 0.28 0.56 0.84];
