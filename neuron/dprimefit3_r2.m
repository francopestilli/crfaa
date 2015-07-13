% dprimefit3_r2.m
%
%        $Id: dprimefit3_r2.m
%      usage: dprimefit3_r2(events,d1)
%         by: franco pestilli
%       date: 01/23/09
%    purpose: compute the contrast
%             threshold at a certain d' value
%
%  This function computes estimates of noise params
%  across adaptation conditions. it works only for
%  distributed target and attended.
%  the rest of the conditions do not have behavior.
%
% this file works with some hand-work.
% dispFit needs to be turned on/off manually
% and only one type of figure can be done at the time.
% with some more work on passing around the right figure legends
% these issues should resolved.
%
% version 3_r uses the inv/deriv of a Naka-Rushton to do the fits
% _r stand for riken. this is done in the last stages of the paper.
%
% it is now called by: makenoisestatisticsfigs2_r.m
%
% franco pestilli 2009.07.03

function [noise exptData] = dprimefit3_r2(data,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% variables that define fit types.
crfFittype =[];    % Possible: {'linear' {'naka' 'Naka' 'NakaRushton' 'nk'}};
psychoFittype=[];  % Possible: {'mllfit' 'quest' 'average'}
dprimefittype=[];  % possible {'nakarushton' 'm'} {'poly' 'polynomial'} {'powerlaw' 'exp' 'exponent' 'power' 'e'} {'skewedgaussian' 'sg'}
adaptationIndex=[];% events and d1 are assumed to be structures for each adaptation condition (1=adapt-0, 2=adapt-28, 3=adapt-100)
doNoise=[];        % 'dual','crf2tvc','tvc2crf'
simulateData=[];   % 1 or 0
doBootstrap=[];    % 1 or 0
numBootstraps=[];  % number of bootstraps to do
dprime=[];         % dprime to compute thresholds
use0contrast=[];   % whether to use the "detection" contrast
whichCRF=[];       % 'crf_distributed_nontarget' 'crf_unattended' 'crf_attended' 'crf_distributed_target';
noisetype=[];      % 'square','additive','both','multiplicative'
noise=[];          % pass in already computed noise parameters
dispFit=[];        % 1 or 0, to display the fit (set to 2 for subsidary fits like TvC etc)
dataTitle=[];      % arbitrary string to use as data title (should be subject data filename)
VisualArea=[];     % which visual area are we working on
figureInfo=[];     % pass in a figure info for different plots to be made on the same figure
conditionIndex=[]; % this is an index for plot info the first one is
minOffset=[];      % this value passed in sets the minimum offset value, when fitting distributed frst and
% attended after i use the distributed offset as min
% value for attended. i hope this solves the problems
% with negative noise reduction index values in fp-v3
maxk=[];           % same thing as min offset
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set arguments
getArgs(varargin,{ ...
 'crfFittype=exp', ...
 'psychoFittype=mllfit', ...
 'dprimefittype=nk', ...
 'doNoise=1', ...
 'simulateData=0', ...
 'doBootstrap=0', ...
 'numBootstraps=200', ...
 'dprime=1', ...
 'use0contrast=1', ...
 'whichCRF=crf_distributed_nontarget', ...
 'noisetype=additive', ...
 'noise=[]', ...
 'dispFit=1', ...
 'dataTitle=[]', ...
 'adaptationIndex=1', ...
 'VisualArea','v1', ...
 'figureInfo',[], ...
 'conditionIndex',[1 1], ...
 'minOffset',0, ...
 'maxk',0},'verbose=0');

% check arguments
if nargin < 1
 help dprimefit3_r2
 return
end

% do a large bootstrap only for the last visual area:
if strcmpi(VisualArea,'v1')
 numBootstraps=numBootstraps;
else
 numBootstraps=2;
end

% setup exptData structure
exptData = setDataStruct(data{adaptationIndex}.events,data{adaptationIndex}.d1,whichCRF, ...
 psychoFittype,simulateData,doBootstrap,numBootstraps,dprime,use0contrast, ...
 VisualArea,adaptationIndex,conditionIndex,minOffset,maxk);


exptData.title = dataTitle;
clear data;

% 1. estimate threshold contrast from data
% by using a maxLikelihood fit:
exptData = computeXcontrast(exptData);

% 2. fit interpolating function to crf data
% exptData = fitCRF(exptData,crfFittype);
exptData.crf = []; 
% NB this field must be set always 
% because some of the calls of the noise model (e.g. multiplicative) require it

% 3. compute noise estimate given the crf and the TvC
switch doNoise
 case {'dual'}
  % fit noise: compute noise value given a noise model and the the behavioral data:
  if isempty(noise)
   noise = fitNoiseDual(exptData,noisetype);
  end
  % compute the tvc2crf fit here
  
 case {'tvc2crf'}
  % fit a tvc
  exptData = fitTVC(exptData,dprimefittype,dispFit>1);
  % fit the CRF from the TvC allowing noise parameters and offset to be adjusted
  if isempty(noise)
   noise = fitCRFfromTvC(exptData,noisetype);
  end
  % display the fit
  [contrast response] = makeModelCRFfromTvC(exptData,noise,noise.responseOffset,dispFit,figureInfo);
  
  % save the response estimated by this noise model:
  exptData.noise.response = response;
  exptData.noise.contrast = contrast;
  
 case {'crf2tvc'}
  % thsi part has not been used for so long that probably doesn't work
  keyboard
  
  % fit the TvC from the CRF allowing noise parameters to be adjusted
  if isempty(noise)
   noise = fitTvCfromCRF(exptData,noisetype);
  end
  % display the fit
  threshold = makeModelTvCfromCRF(exptData.pedestals,noise,exptData.crf,dprime, ...
   exptData,dispFit,figureInfo);
  
 otherwise
  keyboard
end

% %%%%%%%%%%%%%%%%%%%%% end main call %%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%
% makeModelCRFfromTvC %
%%%%%%%%%%%%%%%%%%%%%%%
function [contrast response] = makeModelCRFfromTvC(exptData,noise,responseOffset,dispFit,figureInfo)

if nargin < 4,dispFit = 1;end
% disp(sprintf('[dprimefit3_r2:makeModelCRFfromTvC] start'));

% grab the dprime
dprime = exptData.behavior.dprime;

% we always start the curve at a response of 0
contrast = exptData.pedestalsCRF(1);

response = responseOffset;
deltaContrast = 1;
while (contrast(end) < 1) %
 % get the threshold contrast for this pedestal from the fit
 % of the behaviorally measured TvC function
 deltaContrast = getThresholdContrast(contrast(end),exptData.behavior);
 
%  disp(sprintf('[dprimefit3_r2:makeModelCRFfromTvC] current DeltaC %s',num2str(deltaContrast)));

 % get the sigma
 sigma = sqrt(crfsigma(contrast(end),noise,exptData.crf)^2+crfsigma(contrast(end)+deltaContrast,noise,exptData.crf)^2);
 
 % now compute the response at the pedestal+delta contrast
 contrast(end+1) = contrast(end)+deltaContrast;
 response(end+1) = dprime*sigma+response(end);
 
end
% disp(sprintf('[dprimefit3_r2:makeModelCRFfromTvC] done'));


%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   makeModelTvCfromCRF   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%
function threshContrast = makeModelTvCfromCRF(contrasts,noise,crffit,dprime,exptData,dispFit,figureInfo)

% optimization parameters
maxIter = inf;
MaxFunEvals = inf;
optimParams = optimset('LevenbergMarquardt','on','MaxIter',maxIter,'MaxFunEvals',MaxFunEvals,'Display','off','Diagnostics','off');

% use lsqnonlin to find function minimum
minParams = zeros(1,length(contrasts));
initParams = 0.1*ones(1,length(contrasts));
maxParams = ones(1,length(contrasts));
global numCalls;numCalls = 0;
[threshContrast resnorm residual exitflag output lambda jacobian] = lsqnonlin(@findThreshold,initParams, ...
 minParams,maxParams,optimParams,contrasts,noise,crffit,dprime);

if dispFit
 %  minx = 0.01;
 contrastNum =  length(exptData.pedestalsTvC);
 if contrastNum == 8
  contrastStart = 2;
 else
  contrastStart = 1;
 end
 
 figureHandle =  figureInfo.handle;
 figureName = figureInfo.name;
 figure(figureHandle);
 %%%%%%%%%%% TvC
 figNum = exptData.plotInfo.subplot(4);
 rows = exptData.plotInfo.subplot(1);
 cols = exptData.plotInfo.subplot(2);
 
 subplot(rows,cols,figNum),plot(contrasts(contrastStart:end),threshContrast(contrastStart:end), ...
  'k-','LineWidth',exptData.plotInfo.LineWidth,'Color',exptData.plotInfo.thisColor);
 hold on;
 subplot(rows,cols,figNum),plot(.01,.001,'k*',1,1,'k*')
 for cNum = contrastStart:contrastNum
  yLow = exptData.behavior.tvcCI(cNum,1);
  yHigh = exptData.behavior.tvcCI(cNum,2);
  subplot(rows,cols,figNum)
  myerrorbar(exptData.pedestalsTvC(cNum),exptData.behavior.tvc(cNum), ...
   'Symbol',exptData.plotInfo.thisSymbol(1),'MarkerFaceColor',exptData.plotInfo.MarkerFaceColor, ...
   'MarkerEdgeColor',exptData.plotInfo.MarkerEdgeColor,'MarkerSize',exptData.plotInfo.MarkerSize, ...
   'Color',exptData.plotInfo.MarkerFaceColor);
 end
 set(gca,...
  'FontName','Helvetica','FontSize',exptData.plotInfo.Fsize, ...
  'PlotBoxAspectRatio',exptData.plotInfo.PlotBoxAspectRatio, ...
  'XLim', exptData.plotInfo.Xlim,'YLim', exptData.plotInfo.YLim,...
  'LineWidth',exptData.plotInfo.LineWidth,'TickLength',exptData.plotInfo.TickLength, ...
  'yScale','log','xScale','log','XTick', exptData.pedestals, 'XTickLabel', 100*exptData.pedestals ,...
  'YTick', exptData.plotInfo.YTicksTvC, 'YTickLabel', exptData.plotInfo.YTicksLabelTvC, ...
  'XColor',exptData.plotInfo.XYColor,'YColor',exptData.plotInfo.XYColor,'Box','off');
 
 if exptData.plotInfo.plotLabel
  xlabel('Contrast (%)','FontSize',exptData.plotInfo.Fsize);
  ylabel(sprintf('Threshold\n(% contrast)','%'),'FontSize',exptData.plotInfo.Fsize);
 end
 title(sprintf('Inferred TvC\n%s',getNoiseParamsStr(noise)),'FontSize',exptData.plotInfo.Fsize);
 % myaxis;drawnow
 
 %%%%%%%%%%% CRF
 % some info
 xindex = find(crffit.fit.fitx>exptData.plotInfo.Xlim(1));
 yindex = find(crffit.fit.fity>10);
 if ~isempty(yindex)
  yindex = yindex(1);
 else
  yindex = length(crffit.fit.fity);
 end
 contrastNum =  length(exptData.pedestalsCRF);
 
 figNum = exptData.plotInfo.subplot(3);
 rows = exptData.plotInfo.subplot(1);
 cols = exptData.plotInfo.subplot(2);
 
 subplot(rows,cols,figNum),plot(crffit.fit.fitx(xindex:yindex),crffit.fit.fity(xindex:yindex),'k-','LineWidth',exptData.plotInfo.LineWidth,'Color',exptData.plotInfo.thisColor);hold on
 subplot(rows,cols,figNum),plot(.01,0,'k*',1,1,'k*')
 for cNum = contrastStart:contrastNum
  subplot(rows,cols,figNum)
  myerrorbar(exptData.pedestalsCRF(cNum),exptData.use_crf(cNum), 'yError',exptData.use_crf_ste(cNum), ...
   'Symbol',exptData.plotInfo.thisSymbol(1),'MarkerFaceColor',exptData.plotInfo.MarkerFaceColor, ...
   'MarkerEdgeColor',exptData.plotInfo.MarkerEdgeColor,'MarkerSize',exptData.plotInfo.MarkerSize, ...
   'Color',exptData.plotInfo.MarkerFaceColor);
 end
 set(gca,...
  'FontName','Helvetica','FontSize',exptData.plotInfo.Fsize, ...
  'PlotBoxAspectRatio',exptData.plotInfo.PlotBoxAspectRatio, ...
  'XLim', exptData.plotInfo.Xlim,'YLim', exptData.plotInfo.YLim,...
  'LineWidth',exptData.plotInfo.LineWidth,'TickLength',exptData.plotInfo.TickLength, ...
  'xScale','log', 'XTick', exptData.pedestals, 'XTickLabel', 100*exptData.pedestals ,...
  'YTick', exptData.plotInfo.YTicksCRF, 'YTickLabel', exptData.plotInfo.YTicksLabelCRF ,...
  'XColor',exptData.plotInfo.XYColor,'YColor',exptData.plotInfo.XYColor,'Box','off');
 
 set(figureHandle,'PaperPosition',[.25 .25 8 10.5]);
 set(figureHandle,'PaperOrientation','Portrait');
 
 title(sprintf('CRF: %s %s\n%s',exptData.used_crf,exptData.plotInfo.VisualArea,exptData.title),'Interpreter','none','FontSize',exptData.plotInfo.Fsize);
 if exptData.plotInfo.plotLabel
  xlabel('Contrast (%)','FontSize',exptData.plotInfo.Fsize);
  ylabel(sprintf('Response\n(% signal change)','%'),'FontSize',exptData.plotInfo.Fsize);
 end
 drawnow
 if figureInfo.savePlots
  disp(sprintf('Saving figure %s',figureName));
  eval(sprintf('print(%s,''-depsc2'',''-tiff'',''-r300'', ''%s'')', num2str(figureHandle),figureName));
 end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%
%    computeXcontrast    %
%%%%%%%%%%%%%%%%%%%%%%%%%%
function exptData = computeXcontrast(exptData)
% function exptData = computeXcontrast(exptData)
%
% this function uses computeTvC to compute
% the contrast to be used for fitting and plotting
% te distributed-target and attended condition
%

% compute the contrast
exptData = computeTvC(exptData,exptData.behavior.psychoFittype);

% save the contrasts obtained
% adding 1/2 of it to the pedestal contrast:
% used to fit and plot the crf:
exptData.pedestalsCRF = exptData.pedestals+exptData.behavior.tvc.thisTvC/2;


%%%%%%%%%%%%%%%%%
% setDataStruct %
%%%%%%%%%%%%%%%%%
function exptData = setDataStruct(events,d1,whichCRF,psychoFittype,simulateData, doBootstrap,numBoostraps,dprime,use0contrast, VisualArea,a,conditionIndex,minOffset,maxk)

if strcmp(whichCRF,'crf_attended')
 indexAttention = 1;
elseif strcmp(whichCRF,'crf_distributed_target')
 indexAttention = 2;
else
 keyboard
end


% generate some nice colors:
numC =60;
for i = 1:numC
 c{i} = getSmoothColor(i,numC,'hsv');
 %  plot(i,1,'ko','MarkerFaceColor',c{i},'MarkerSize',26);hold on
 %  text(i,1,sprintf('%i',i),'HorizontalAlignment','Center','Color',[0 0 0]);
end

% choose the ones i like:
reds    = {c{2}  c{5}  c{7} };
greens  = {c{23} c{17} c{15}};
blues   = {c{40} c{37} c{34}};
purples = {c{46} c{49} c{54}};

for i = 1:3 % adaptation
 plotInfo.thisColors{i} = {reds{i} blues{i} purples{i} greens{i}};
end
exptData.plotInfo.thisSymbol = 'o-';

exptData.plotInfo.thisColor = plotInfo.thisColors{conditionIndex(1)}{conditionIndex(2)};
exptData.plotInfo.MarkerFaceColor = plotInfo.thisColors{conditionIndex(1)}{conditionIndex(2)};

exptData.plotInfo.XYColor    = [.4 .4 .4];
exptData.plotInfo.Fsize      = 8;
exptData.plotInfo.MarkerSize = 10;
plotInfo.MarkerEdgeColor     = [.4 .4 .4];
exptData.plotInfo.LineWidth  = 1;
exptData.plotInfo.TickLength = [0.025 .01];
exptData.plotInfo.PlotBoxAspectRatio = [1 1 1];
exptData.plotInfo.YLim       = [0 1.3];
exptData.plotInfo.Xlim       = [.01,1];
exptData.plotInfo.YTicksCRF  = [0 .25 .5 .75 1 1.25];
exptData.plotInfo.YTicksLabelCRF = [0 .25 .5 .75 1 1.25];
exptData.plotInfo.YTicksTvC  = [.01 .1 1 10 100];
exptData.plotInfo.YTicksLabelTvC = 100*[.01 .1 1 10 100];

% set up subplot info and title info:
exptData.plotInfo.VisualArea = VisualArea;
switch VisualArea
 case {'v1' 'V1'}
  exptData.plotInfo.plotLabel = 0;
  exptData.plotInfo.subplot = [4,2,1,2];
 case {'v2' 'V2'}
  exptData.plotInfo.plotLabel = 0;
  exptData.plotInfo.subplot = [4,2,3,4];
 case {'v3' 'V3'}
  exptData.plotInfo.plotLabel = 0;
  exptData.plotInfo.subplot = [4,2,5,6];
 case {'v4' 'V4'}
  exptData.plotInfo.plotLabel = 1;
  exptData.plotInfo.subplot = [4,2,7,8];
end

% behavior (tvc) parameters:
exptData.behavior.psychoFittype = psychoFittype;
exptData.behavior.simulateData = simulateData;
exptData.behavior.doBootstrap = doBootstrap;
exptData.behavior.numBoostraps = numBoostraps;
exptData.behavior.dprime = dprime;% set the d-prime value for the threshold wanted:

% get information about trial conditions and observer responses
% events.allEvents attention condition = 1 (attended) or = 2 (distributed)
useTrials = find(events.allEvents(1,:) == indexAttention);
exptData.n = length(useTrials);
exptData.pedContrast = events.allEvents(9,useTrials);
exptData.deltaContrast = events.allEvents(8,useTrials);
exptData.correctIncorrect = events.allEvents(4,useTrials);
exptData.use0contrast = use0contrast;

% get pedestal contrasts
exptData.pedestals = sort(unique(exptData.pedContrast));
exptData.pedestals(1) = exptData.plotInfo.Xlim(1);

exptData.pedestalsTvC = [.00875 .0175 .035 .07 .14 .28 .56 .84]; % used to compute and plot the tvc

exptData.fitcontrast = logspace(log10(.0001),log10(1),100);

% get contrast response function
exptData.crf_distributed_target = d1.amplitude(17:24);
exptData.crf_attended = d1.amplitude(1:8);

% and standard error
exptData.crf_distributed_target_ste = d1.amplitudeSTE(17:24);
exptData.crf_attended_ste = d1.amplitudeSTE(1:8);

% used crf and ste:
exptData.use_crf = exptData.(whichCRF);
exptData.use_crf_ste = exptData.([whichCRF,'_ste']);
exptData.used_crf = whichCRF;

% set the minimum value for noise offset during fit.
% for attended this value is set to be the value for distrbuted
exptData.minOffset = minOffset;
exptData.maxk = maxk;


%%%%%%%%%%%%%%%%%%%%%%%
%    nakarushtonTvC   %
%%%%%%%%%%%%%%%%%%%%%%%
function deltac = nakarushtonTvC(pedestals,params)
rmax  = 1;
n     = params(1);
c50   = params(2);
T     = params(3)*c50;
beta  = params(4);
d     = params(5); % differentiation factor

c     = pedestals;

% disp(sprintf('Current Params: Rmax:%s - n:%s - c50: %s\nd: %s - T: %s - Beta: %s.',...
%  num2str(rmax), ...
%  num2str(n),    ...
%  num2str(c50),  ...
%  num2str(d),    ...
%  num2str(T),    ...
%  num2str(beta)));

% call the function for the second
% (high contrast dip):
y = dip2(c,T,beta);

% breaking down the function into
% components for simplicity (kinda)
v1 = c.^n;
v2 = c50^n;
v3 = v1+v2;
v4 = v1./v3;
v5 = d/rmax;
v6 = (1./(v5+v4))-1;
v7 = v2./v6;
v8 = v7.^(1/n);
v9 = v8 - c;

% deltac = M.*(v9.*y);
deltac = (v9.*y);


%%%%%%%%%%%%%%
%   fitTVC   %
%%%%%%%%%%%%%%
function exptData = fitTVC(exptData,dprimefittype,dispFit)
% this function fits different type of TvC to the data:

% disp(sprintf('[dprimefit3_r2:fitTVC] start, (%s)',dprimefittype));
% global figLabel;

% fit the crf with the appropriate function
switch dprimefittype
 case {'nakarushton' 'nk'}
  
  % fit a naka-rushton sigmodal function:
  %            slope  thr                           scaleC50  dipSlope   d
  initParams = [2     min(exptData.pedestalsTvC)-min(exptData.pedestalsTvC)/10     4         .5        .00006];
  
  %            slope  thr      scaleC50   dipSlope   d
  minParams = [.1     0.001    2.5          .1         0];
  maxParams = [ 3.5     0.25   10           3.5       .1];
  
  
  pedestals = [exptData.pedestalsTvC(2)/2 exptData.pedestalsTvC(2:end)];
  thresholds = exptData.behavior.tvc;
  
  bestFit = fitnakarushtonTvC(pedestals,thresholds,initParams,minParams,maxParams);
  exptData.behavior.tvcfit.fitParams.fit = bestFit;
  
 case {'nk2' 'nkn' 'nakarushtonNumerical' 'nknum' 'NakaNum' 'nakanum'}
  
  % fit a naka-rushton sigmodal function:
  %             <n>    <thr>    <p>     <dr>    <scaleC50> <dipSlope>
  initParams = [1.6    0.0120    0.9   0.001     50       1];
  minParams  = [.2    0.0001    0.01  0.000001   1       .1];
  maxParams  = [6      0.75      1       1       100      4];
  
  pedestals = exptData.pedestalsTvC;
  thresholds = exptData.behavior.tvc.thisTvC;

  thresholds_ste = exptData.behavior.tvc.thisTvCste;
  
  bestFit = fitnakarushtonTvCNumerical(pedestals,thresholds,thresholds_ste,initParams,minParams,maxParams);
  exptData.behavior.tvcfit.fitParams.fit = bestFit;
  
  
 case {'poly' 'polynomial'}
  % fit a polynomial:
  exptData.behavior.tvcfit.polyorder = 2;
  fitParams = fitpolynomial(exptData.pedestalsTvC,exptData.behavior.tvc,exptData.behavior.tvcfit.polyorder);
  exptData.behavior.tvcfit.fitParams = fitParams;
  
 case {'powerlaw' 'exp' 'exponent' 'power' 'e'}
  % fit a powerlaw (exponential function (basically a straight line on loglog axis)):
  fitParams = fitexponent(exptData.pedestalsTvC,exptData.behavior.tvc,dispFit);
  exptData.behavior.tvcfit.fitParams = fitParams;
  
 case {'skewedgaussian' 'sg' 'skew' 'gaussian' 'g'}
  % fit a a skewed gaussian on log axis:
  fitParams = skewfit(exptData.pedestalsTvC(2:end),exptData.behavior.tvc(2:end),dispFit);
  exptData.behavior.tvcfit.fitParams = fitParams;
  
 otherwise
  %disp(sprintf('[dprimefit3_r2:fitTVC] Unknown TvC psychoFittype: %s',dprimefittype));
  keyboard
end

% set the fittype
exptData.behavior.tvcfit.type = dprimefittype;

% make sure we can compute the thresholds properly
dispFit = 1;
if dispFit
 fitContrast = 0.001:0.001:1;
 %disp(sprintf('[dprimefit3_r2:fitTVC] showingthe fit: %s',dprimefittype));
 for i = 1:length(fitContrast)
  threshold(i) = getThresholdContrast(fitContrast(i),exptData.behavior);
 end
 
 smartfig('dprimefit_tvcfit','reuse');
 cla
 loglog(fitContrast,threshold,'r-');
 hold on
 loglog(exptData.pedestalsTvC,exptData.behavior.tvc.thisTvC,'bo')
 % myaxis;
 %axis('square')
 title(sprintf('Fit TvC using:\n%s',dprimefittype));
 ylabel('Contrast threshold')
 xlabel('Pedestal contrast')
 legend(sprintf('TvC fitted as: %s\ndprime:%0.2f',dprimefittype,exptData.behavior.dprime),4)
 drawnow
end
disp(sprintf('[dprimefit3_r2:fitTVC] done'));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   getThresholdContrast   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function threshold = getThresholdContrast(pedestal,behavior)
% check that we have done the fit
if ~isfield(behavior,'tvcfit')
 %disp(sprintf('[dprimefit3_r2] tvcfit has not been done'));
 keyboard
end

% return the threshold contrast for the given pedestal. Depends
% on which fit type we did.
switch (behavior.tvcfit.type)
 case {'nakarushton' 'nk'}
  params = behavior.tvcfit.fitParams.fit.fitParams;
  threshold = nakarushtonTvC(pedestal,params);
  
 case {'nkn' 'nakarushtonNumerical' 'NakaNum' 'nknum' 'nakanum'}
  params = behavior.tvcfit.fitParams.fit.fitParams;
  threshold = nakarushtonTvCNumerical(pedestal,params);

 case {'poly' 'polynomial'}
  threshold = 0;
  for p = 0:behavior.tvcfit.polyorder
   threshold = threshold+behavior.tvcfit.fitParams(p+1)*pedestal^p;
  end
  
 case {'powerlaw' 'exp' 'exponent' 'power' 'e'}
  x = pedestal;
  a = behavior.tvcfit.fitParams.a;
  tau = behavior.tvcfit.fitParams.tau;
  offset = behavior.tvcfit.fitParams.offset;
  
  % use the model in fitexponent and the fit parameters to get the
  % model thresholds:
  threshold = a*exp(-x/tau)+offset;
  
 case {'skewedgaussian' 'sg' 'skew' 'gaussian' 'g'}
  threshold = skewFun(behavior.tvcfit.fitParams.params,pedestal);
  
 otherwise
  %disp(sprintf('[dprimefit3_r2:getThresholdContrast) Unknown TvC fittype: %s',behavior.tvcfit.type));
  keyboard
end


%%%%%%%%%%%%%%%%%%%%%%%
%    fitpolynomial    %
%%%%%%%%%%%%%%%%%%%%%%%
function fitParams = fitpolynomial(x,y,polyorder)

% get the glm
x = x';y = y';
glm = getPolynomialGLM(x,polyorder);

% fit the polynomial
fitParams = pinv(glm)*y;

% display the fit
dispFit = 1;
if dispFit
 smartfig('fitpolynomial','reuse');
 plotx = 0.001:0.001:1;
 plotglm = getPolynomialGLM(plotx,polyorder);
 plot([.001 x(2:end)'],y,'ko');
 hold on
 plot(plotx,plotglm*fitParams,'r-');
 axis([0 1 0 max(plotglm*fitParams)])
end


%%%%%%%%%%%%%%%%%%%%%%%%
%   getPolynomialGLM   %
%%%%%%%%%%%%%%%%%%%%%%%%
function glm = getPolynomialGLM(x,polyorder)

% create the glm
for i = 0:polyorder
 glm(:,i+1) = x.^i;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    fitnakarushtonTvC    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%
function bestFit = fitnakarushtonTvCNumerical(c,t,t_ste,initParams,minParams,maxParams)

% optimization parameters
maxIter = inf;
MaxFunEvals = inf;
optimParams = optimset('LevenbergMarquardt','on', ...
 'MaxIter',maxIter,         ...
 'MaxFunEvals',MaxFunEvals, ...
 'Display','off',           ...
 'Diagnostics','off',       ...
 'TolFun',10^-15);

% use lsqnonlin to find function minimum
[fitParams resnorm residual exitflag output lambda jacobian] = lsqnonlin(@resNakarushtonNumerical,initParams,minParams,maxParams,optimParams,c,t);

% compute R2 values and save them in the fit:
% (1) reshape the errors by condition and contrast
bestFit.err = 10.^resNakarushtonNumerical(fitParams,c,t);


% (2) run through conditions and compute r2
% for each one of them
thisError = bestFit.err;
thisdeltac = t;
bestFit.r2 = 1-var(thisError)/var(log10(thisdeltac));

% doing  chi-squared test of goodness of fit.
bestFit.chi2.chi2value = sum((thisError./(sum(t_ste')./2)).^2);
bestFit.chi2.M = 6;
bestFit.chi2.N = 8;
bestFit.chi2.nparams = bestFit.chi2.N - bestFit.chi2.M;
bestFit.chi2.dataVariability = t_ste; 
bestFit.chi2.modelFitError = thisError;
bestFit.chi2.p = 1 - gammainc(0.5*bestFit.chi2.chi2value,0.5*bestFit.chi2.nparams); 


% if we have the best fit then keep it.
bestFit.fitParams   = fitParams;

% add here a smoothfit function.
bestFit.x = c(1):.0001:.9;

bestFit.y = nakarushtonTvCNumerical(bestFit.x,bestFit.fitParams);

bestFit.resnorm = resnorm;
bestFit.residual = residual;
bestFit.exitflag = exitflag;
bestFit.output = output;
bestFit.lambda = lambda;
bestFit.jacobian = jacobian;


%%%%%%%%%%%%%%%%%%%%
% contrastResponse %
%%%%%%%%%%%%%%%%%%%%
function r = contrastResponse(c,params)
rmax = 1;
n   = params(1);
c50 = params(2);
p   = params(3);

% disp(sprintf('[contrastResponse] n: %s, c50: %s, p: %s.',num2str(n),num2str(c50),num2str(p)))
r = rmax.*((c.^n)./((c.^(n*p)) + (c50.^(n*p))));


%%%%%%%%%%%%%%%%%%%%%%%%%%%
% nakarushtonTvCNumerical %
%%%%%%%%%%%%%%%%%%%%%%%%%%%
function deltac = nakarushtonTvCNumerical(pedestals,otherparams)
% rmax   = 1; % fixed to '1';
% n      = params(1);
% c50    = params(2);
% p      = params(3);
% dr     = params(4); % differentiation factor

NKparams = otherparams(1:4);
T    = otherparams(5)*otherparams(2);
beta = otherparams(6);


% disp(sprintf('[nakarushtonTvCNumerical] n: %s, c50: %s, p: %s, d: %s, T: %s, beta: %s.', ...
%  num2str(NKparams(1)),num2str(NKparams(2)),num2str(NKparams(3)),num2str(NKparams(4)),num2str(T),num2str(beta)));

for i = 1:length(pedestals)
 c = pedestals(i);
%  disp(sprintf('[nakarushtonTvCNumerical] pedestals [%s], c: %s.',num2str(pedestals),num2str(c)))

 deltac(i) = fzero(@deltaContrast,[0 10],[],c,NKparams);
%  disp(sprintf('[nakarushtonTvCNumerical] deltac: %s, c: %s.',num2str(deltac(i)),num2str(c)))
end

% call the function for the second
% (high contrast dip):
y = dip2(pedestals,T,beta);

deltac = deltac.*y;
% disp(sprintf('[nakarushtonTvCNumerical] deltac: %s.',num2str(deltac)))


%%%%%%%%%%%%%%%%%
% deltaContrast %
%%%%%%%%%%%%%%%%%
function m = deltaContrast(dC,c,params)
% rmax   = 1; % fixed to '1';
% n      = params(1);
% c50    = params(2);
% p      = params(3);
dr       = params(4); % differentiation factor


rMIN = contrastResponse(c,params(1:3));
rMAX = contrastResponse(c+dC,params(1:3));

m = rMAX - rMIN - dr;
% disp(sprintf('[deltaContrast] m: %s, d: %s.',num2str(m),num2str(dr)))
% disp(sprintf('[deltaContrast] rMAX: %s, rMIN: %s.',num2str(rMIN),num2str(rMAX)))


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    resNakarushtonNumerical    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function residual = resNakarushtonNumerical(params,pedestal,data)

model = nakarushtonTvCNumerical(pedestal,params);
residual = log10(data)-log10(model);


dispFit = 1;
if dispFit
 smartfig('dprimefit_nakaTvCNumerical','reuse');
 cla;
 loglog(pedestal,data,'bs');hold on
 loglog(pedestal,model,'r.-');
 axis('square')
 
 rmax  = '1';%num2str(params(1)); % it is now fixed to 1
 n     = num2str(params(1));
 c50   = num2str(params(2));
 p     = num2str(params(3));
 d     = num2str(params(4));
 
 txt = sprintf('rmax: %s, n: %s,\nc50: %s, p: %s\nd: %s.',rmax,n,c50,p,d);
 title(txt)
end
% $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ END NUMERICAL $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ %%


%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    fitnakarushtonTvC    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%
function bestFit = fitnakarushtonTvC(pedestals,deltac,initParams,minParams,maxParams)

% optimization parameters
maxIter = 10^3;
MaxFunEvals = 10^3;
optimParams = optimset('LevenbergMarquardt','on', ...
 'MaxIter',maxIter,         ...
 'MaxFunEvals',MaxFunEvals, ...
 'Display','off',           ...
 'Diagnostics','off',       ...
 'TolFun',10^-8);

% use lsqnonlin to find function minimum
[fitParams resnorm residual exitflag output lambda jacobian] = lsqnonlin(@resNakarushton,initParams,minParams,maxParams,optimParams,pedestals,deltac);


% compute R2 values and save them in the fit:
% (1) reshape the errors by condition and contrast
bestFit.err = resNakarushton(fitParams,pedestals,deltac);

% (2) run through conditions and compute r2
% for each one of them
thisError = bestFit.err;
thisdeltac = deltac;
bestFit.r2 = 1-var(thisError)/var(thisdeltac);

% if we have the best fit then keep it.
bestFit.fitParams   = fitParams;

% add here a smoothfit function.
bestFit.x = 0.0001:.0001:1;
bestFit.y = nakarushtonTvC(bestFit.x,bestFit.fitParams);

bestFit.resnorm = resnorm;
bestFit.residual = residual;
bestFit.exitflag = exitflag;
bestFit.output = output;
bestFit.lambda = lambda;
bestFit.jacobian = jacobian;


%%%%%%%%%%%%%%%%%%%%%%%%
%    resNakarushton    %
%%%%%%%%%%%%%%%%%%%%%%%%
function residual = resNakarushton(params,pedestal,data)

model = nakarushtonTvC(pedestal,params);

residual = log10(data)-log10(model);

dispFit = 1;
if dispFit
 smartfig('dprimefit_nakaTvC','reuse');
 cla;
 loglog([.001 pedestal(2:end)],data,'bs');hold on
 loglog([.001 pedestal(2:end)],model,'r.-');
 axis('square')
 axis([.001 1 min(min([data(:) model(:)])) 1]);
 
 rmax  = '1'; % it is now fixed to 1
 n     = num2str(params(1));
 c50   = num2str(params(2));
 T     = num2str(params(3));
 beta  = num2str(params(4));
 d     = num2str(params(5));
 
%  txt = sprintf('rmax: %s, n: %s, c50: %s,\nd: %s, M: %s\nT: %s, beta: %s.',rmax,n,c50,d,M,T,beta);
 txt = sprintf('rmax: %s, n: %s, c50: %s,\nd: %s\nT: %s, beta: %s.',rmax,n,c50,d,T,beta);

title(txt)
end


%%%%%%%%
% dip2 %
%%%%%%%%
function y = dip2(c,T,beta)
% this one fixes the decrease in t at high contrast
% parameters of the copressive function at high contrast
% T center (threshold) of the function
% beta slope of the function
y = .015+exp(-(c./T).^beta-.025);


%%%%%%%%%%%%%%%%%%%%%%%%
%    findThreshold     %
%%%%%%%%%%%%%%%%%%%%%%%%
function residual = findThreshold(params,contrasts,noise,crffit,dprime)

for i = 1:length(contrasts)
 % get the pedestal contrast passed in:
 thisContrast = contrasts(i);
 
 % compute the response difference and standard deviation of response
 % responseDifference = mu and sigma is standard deviation
 responseDifference(i) = crf(thisContrast+params(i),crffit) - crf(thisContrast,crffit);
 sigma(i) = sqrt(crfsigma(thisContrast+params(i),noise,crffit)^2+crfsigma(thisContrast,noise,crffit)^2);
 
 % now compute dprime given signal and noise for each pedestal:
 modelDprime(i) = responseDifference(i)/sigma(i);
end

% compute the sum of squared difference between model and wanted d-prime:
residual = modelDprime-dprime;

dispFit = 0;
if dispFit
 global numCalls;
 smartfig('dprimefit_findThreshold','reuse');
 clf;
 subplot(1,2,1);
 loglog(contrasts,p);
 subplot(1,2,2);
 semilogx(contrasts,p);
end


%%%%%%%%%%%
%   crf   %
%%%%%%%%%%%
function response = crf(contrast,crffit)

% get the contrast response
switch crffit.fittype
 case {'linear'}
  %  linear fit
  response = crffit.fit.m*contrast+crffit.fit.b;
  
 case {'naka' 'Naka' 'NakaRushton' 'nr'}
  % recompose the naka-rushton sigmodal function:
  Rmax = crffit.fit.params(1);
  c50 = crffit.fit.params(2);
  n   = crffit.fit.params(3);
  offset   = crffit.fit.params(4);
  response = Rmax*(contrast.^n)./(contrast.^n+c50.^n)+offset;
  
 case {'exp' 'EXP' 'e'}
  % recompose the exponential function:
  [response c50 n offset] = testexp(fitparams,contrast,'111',[1 1]);
  
 case {'interp' }
  [c index] = sort(crffit.contrast);
  r = crffit.response(index);
  response = interp1q(c,r,contrast);
  
 otherwise
  %disp(sprintf('[dprimefit3_r2:crf] Unknown crf fittype: %s',crffit.fittype));
  keyboard
end


%%%%%%%%%%%%%%%%
%   crfsigma   %
%%%%%%%%%%%%%%%%
function sigma = crfsigma(contrast,noise,crffit)

switch (noise.noisetype)
 case {'additive'} % noise is fixed, independent of contrast
  sigma = noise.k;
  
 case {'multiplicative'} % noise scales with response
  sigma = noise.k*crf(contrast,crffit);
  
 case {'both'} % noise is aditive and scales with response
  sigma = noise.k+noise.k1*crf(contrast,crffit);
  
 case {'square'} % noise is additive and scales with response
  sigma = noise.k+noise.k1*crf(contrast,crffit)+noise.k2*(crf(contrast,crffit))^2;
  
 otherwise
  %disp(sprintf('[dprimefit3_r2:crfsigma] Unknown model type %s',p.noisetype));
end


%%%%%%%%%%%%%%
%   fitCRF   %
%%%%%%%%%%%%%%
function exptData = fitCRF(exptData,crfFittype)
disp(sprintf('[dprimefit3_r2:fitCRF] start... FitType: %s',crfFittype));

% fit the crf with the appropriate function
switch crfFittype
 case {'linear'}
  %  linear fit
  exptData.crf.fit = myregress(exptData.pedestalsCRF,exptData.use_crf);
  
 case {'naka' 'Naka' 'NakaRushton' 'NR'}
  % fit a naka-rushton sigmodal function:
  testType = '1111';
  exptData.crf.fit = fitsigmoidtest(exptData.pedestalsCRF,exptData.use_crf, testType);
  
 case {'exp' 'EXP' 'e'}
  % fit a exponential function:
  testType = '111';
  exptData.crf.fit = fitexptest(exptData.pedestalsCRF,exptData.use_crf, testType);
  
 case {'interp'}
  exptData.crf.contrast = exptData.pedestalsCRF;
  exptData.crf.response = exptData.use_crf;
 otherwise
  %disp(sprintf('[dprimefit3_r2:fitCRF] Unknown crf fittype: %s',crfFittype));
  keyboard
end

% set the fittype
exptData.crf.fittype = crfFittype;

% create a smooth function to display
if ~isfield(exptData.crf,'fit') || ~isfield(exptData.crf.fit,'fitx')
 exptData.crf.fit.fitx = 0:0.01:1;
 for i = 1:length(exptData.crf.fit.fitx)
  exptData.crf.fit.fity(i) = crf(exptData.crf.fit.fitx(i),exptData.crf);
 end
end

disp(sprintf('[dprimefit3_r2:fitCRF] done'));


%%%%%%%%%%%%%%%%%%
%   computeTvC   %
%%%%%%%%%%%%%%%%%%
function exptData = computeTvC(exptData,psychoFittype)
% this function computes a tvc
%disp(sprintf('[dprimefit3_r2:computeTvC] start...'));

% fit the crf with the appropriate function
switch psychoFittype
 case {'average'}
  %  compute mean contrast across trials
  for pedc = 1:length(exptData.pedestalsTvC)
   whichPedc = exptData.pedestalsTvC(pedc);
   exptData.behavior.tvc(pedc) = median(exptData.deltaContrast(find(exptData.pedContrast==whichPedc)));
   exptData.behavior.tvcSD(pedc) = std(exptData.deltaContrast(find(exptData.pedContrast==whichPedc)));
  end
  exptData.behavior.tvcCI = [];
 case {'quest' 'Quest'}
  % file name to load (make this smart):
  filename = '/Users/frakkopesto/Desktop/071030_stim09.mat';
  [results stimulus] = analyzeQuest(filename);
  exptData.behavior.tvc = 10.^mean(squeeze(results.meanThreshold(:,2,:)),2)';
  exptData.behavior.tvcSD = median(squeeze(results.thresholdSD(:,2,:)),2)';
  exptData.behavior.tvcCI = [];
  
 case {'mllfit'}
  % optimization parameters
  maxIter = inf;
  MaxFunEvals = inf;
  optimParams = optimset('MaxIter',maxIter,'MaxFunEvals',MaxFunEvals,'Display','off','Diagnostics','off');
  
  doBootstrap = exptData.behavior.doBootstrap;
  numBoostraps = exptData.behavior.numBoostraps;
  
  for pedc = 1:length(exptData.pedestalsTvC)
   % get mean delta contrast as threshold guess
   if exptData.pedestalsTvC(pedc) < .0175
    whichPedc = 0;
   else
    whichPedc = exptData.pedestalsTvC(pedc);
   end
   trials = find(exptData.pedContrast==whichPedc);
   deltaContrast = exptData.deltaContrast(trials);
   correctIncorrect = exptData.correctIncorrect(trials);
   numTrials = length(trials);
   
   if doBootstrap
    bootstrpWBLfitp = zeros(numBoostraps,3); % initialize array for bootstrap parameters
    % bootstraping some errorbars:
    for bootStrapN = 1:numBoostraps
     exitflag = 0;
     disp(sprintf('[dprimefit3_r2] Pedestal <%s> - percent bootstrap done <%s>',num2str(whichPedc),num2str(100*bootStrapN/numBoostraps)))
     while exitflag == 0 % repeat the bootstrap if the fit is not succesful
      % randsaple with replacement:
      index = randsample(1:numTrials,numTrials,true);
      bootDeltaContrast = deltaContrast(index);
      bootCorrectIncorrect = correctIncorrect(index);
      
      % init parameters for threshold fit:
      thresholdInit = median(bootDeltaContrast); % mll weight the number of occurrences the median is a good start
      initParams = [thresholdInit 3.5 .05];
      
      [bootstrpWBLfitp(bootStrapN,:) maxMLL exitflag fitOutput] = fminsearch(@mllweibull,initParams,optimParams, ...
       bootDeltaContrast,bootCorrectIncorrect,exptData.pedestalsTvC,exptData.fitcontrast,whichPedc);
      
      % get threshold at requested d-prime
      btThreshold(bootStrapN) = y2x(exptData.behavior.dprime,exptData.fitcontrast,'wbl2dprime',bootstrpWBLfitp(bootStrapN,:));
     end
    end
    
    % use fminsearch to find function for the actual data:
    % init parameters for threshold fit:
    thresholdInit = median(deltaContrast); % mll weight the number of occurrences the median is a good start
    initParams = [thresholdInit 3.5 .05];
    
    deltaContrast = exptData.deltaContrast(trials);
    correctIncorrect = exptData.correctIncorrect(trials);
    [wblfitp maxMLL exitflag fitOutput] = fminsearch(@mllweibull,initParams,optimParams,deltaContrast,correctIncorrect,exptData.pedestalsTvC,exptData.fitcontrast,whichPedc);
    
    % get threshold at requested d-prime
    th = y2x(exptData.behavior.dprime,exptData.fitcontrast,'wbl2dprime',wblfitp);
    
    % fix the threshold to 0 when is returned as 'nan'
    if isnan(th)
     th = 0;
    end
    
    % save info into structure
    exptData.behavior.tvc.thisTvC(pedc) = th;
    exptData.behavior.tvc.thisTvCste(pedc,1:2) = quantile(btThreshold,[.33 .66]);
    exptData.behavior.tvc.thistvcSD = [];
    exptData.behavior.tvc.BootMedian(pedc) =  quantile(btThreshold,.5);
    exptData.behavior.tvc.BootStrapWBLparams{pedc} = bootstrpWBLfitp;
    exptData.behavior.tvc.fitresults{pedc}.wblparams = wblfitp;
    exptData.behavior.tvc.fitresults{pedc}.maxMLL = maxMLL;
    exptData.behavior.tvc.fitresults{pedc}.exitflag = exitflag;
    exptData.behavior.tvc.fitresults{pedc}.fitOutput = fitOutput;
   else
    %disp(sprintf('[dprimefit3_r2:computeTvC] no bootstrap, pedestal %0.4f',whichPedc))
    % use fminsearch to find function minimum
    % init parameters for threshold fit:
    thresholdInit = median(deltaContrast); % mll weight the number of occurrences the median is a good start
    if isfield(exptData,'syntheticThresholdContrast')
     thresholdInit = exptData.syntheticThresholdContrast(pedc);
     % disp(sprintf('[dprimefit3_r2] Using synthetic threshold of %0.4f as initial threshold',exptData.syntheticThresholdContrast(pedc)));
    else
     % disp(sprintf('[dprimefit3_r2] Using %0.4f as initial threshold',thresholdInit));
    end
    initParams = [thresholdInit 3.5 .05];
    
    deltaContrast = exptData.deltaContrast(trials);
    correctIncorrect = exptData.correctIncorrect(trials);
    [wblfitp maxMLL exitflag fitOutput] = fminsearch(@mllweibull,initParams,optimParams, ...
     deltaContrast,correctIncorrect,exptData.pedestalsTvC,exptData.fitcontrast,whichPedc);
    
    % get threshold at requested d-prime
    th = y2x(exptData.behavior.dprime,exptData.fitcontrast,'wbl2dprime',wblfitp);
    
    % fix the threshold to 0 when is returned as 'nan'
    if isnan(th)
     th = 0;
    end
    
    % save info into structure
    exptData.behavior.tvc(pedc) = th;
    exptData.behavior.tvcSD(pedc) = 0;
    exptData.behavior.tvcCI = [];
    exptData.behavior.fitresults{pedc}.wblparams = wblfitp;
    exptData.behavior.fitresults{pedc}.maxMLL = maxMLL;
    exptData.behavior.fitresults{pedc}.exitflag = exitflag;
    exptData.behavior.fitresults{pedc}.fitOutput = fitOutput;
   end
  end
  exptData.behavior.fitparams.initial = initParams;
  exptData.behavior.fitparams.optimParams = optimParams;
  
 otherwise
  % disp(sprintf('[dprimefitmll:computeTvC] Unknown crf tvctype: %s',psychoFittype));
end
exptData.behavior.fittype = psychoFittype;
% disp(sprintf('[dprimefit3_r2:computeTvC] done'));


%%%%%%%%%%%%%%%%%%%%
%     mllweibull   %
%%%%%%%%%%%%%%%%%%%%
function loglike = mllweibull(p,deltaContrast,correctIncorrect,pedestals,fitcontrast,pedC)
% function loglike = mllweibull(p,exptData,pedC)
%
% this function computes a maximum likelihood estimate of
% the psychometric function at each pedestal contrast
% given all delta contrasts presented and all responses
%
% for each trial compute the probability
% of being correct given the paramaters
% and the data

for iTrial = 1:length(deltaContrast)
 % get the pedestal contrast
 pcorrect(iTrial) = weibullFun(deltaContrast(iTrial),p);
end
pincorrect = 1-pcorrect;

% compute the probablity with which the weibull
% generates the observers correct and incorrect
% responses for each trial
pObserverCorrectTrials = pcorrect(find(correctIncorrect==1));
pObserverIncorrectTrials = pincorrect(find(correctIncorrect==0));

% now compute the loglikelihood of the model
% generating the observers responses
loglike = sum(log([pObserverCorrectTrials pObserverIncorrectTrials]));

% we are using fminsearch, so need to minimize
if (p(3) < 0) || (p(3) > .075) || (p(2) < 1 || p(2) > 6.5) || (p(1) < 0 || p(1) > .28575)%p(1) > .575)
 loglike = inf;
else
 loglike = -loglike;
 % display the current fit
 dispFit = 0;
 
 if dispFit
  if pedC == 0
   pedNum = 1;
  else
   pedNum = find(pedestals == pedC);
  end
  behavior.tvc(pedNum) = p(1);
  behavior.tvcSD(pedNum) = 0;
  behavior.fitresults{pedNum}.wblparams = p;
  dispThisWBL(fitcontrast,behavior,deltaContrast,correctIncorrect,pedestals,pedNum);
 end
end


%%%%%%%%%%%%%%%%%%%%%%%%
%       weibullFun     %
%%%%%%%%%%%%%%%%%%%%%%%%
function pcorrect = weibullFun(contrasts,params)
% Weibull function
% pcorrect = weibullFun(contrasts,params)

% parameters
alpha = params(1); % threshold
beta  = params(2); % slope
lambda= params(3); % 1-asymptote

pcorrect = .5+(.5-lambda)*(1-exp(-(contrasts/alpha).^beta));


%%%%%%%%%%%%%%%%%%%
%   dispThisWBL   %
%%%%%%%%%%%%%%%%%%%
function dispThisWBL(fitcontrast,behavior,deltaContrast,correctIncorrect,pedestals,i)

numPedestals = length(behavior.tvc);
mycolor = repmat(linspace(0,.75,numPedestals),3,1);
smartfig(sprintf('wbl_tvc%i',i),'reuse');
cla;

% display fit
wblpc(i,:) = weibullFun(fitcontrast,behavior.fitresults{i}.wblparams);
if any(~isreal(wblpc(i,:)))
 keyboard
else % plot
 semilogx(fitcontrast,wblpc(i,:),'k-', 'Color',mycolor(:,i)');hold on
 semilogx([behavior.tvc(i);behavior.tvc(i)],[0.5;1],'k-', 'Color',mycolor(:,i)');
 text(behavior.tvc(i),1-behavior.fitresults{i}.wblparams(3)+.006,sprintf('p%i',i));
 hold on
end

% not displaying the data anymore.
thisparams = behavior.fitresults{i}.wblparams;
title(sprintf('Pedestal %0.2f Threshold: %0.4f Slope: %0.4f asymptote: %0.4f numTrials: %i',pedestals(i),thisparams(1),thisparams(2),1-thisparams(3),length(correctIncorrect)));
xlabel('Contrast');
ylabel('Probability correct');
% legend('TvC fit','distributed non target data',4);
%axis('square'), axis([fitcontrast(1) 1 .5 1]);
drawnow


%%%%%%%%%%
%  y2x   %
%%%%%%%%%%
function x = y2x(y,xvals,fun,p)
% find a x value (e.g., a threshold)
% given a function (e.g., a psychometric function)
% a set of values for x (xvals)
% and the value on the 'y' axis
%

eval(sprintf('yvals = %s(%s,%s);',fun,'xvals','p'));

% only use unique values otherwise interp1 complains:
[yvals i] = unique(yvals);
x = interp1(yvals,xvals(i),y,'pchip');


%%%%%%%%%%%%%%%%%%%
%    wbl2dprime   %
%%%%%%%%%%%%%%%%%%%
function dp = wbl2dprime(x,params)
% transform pcorrect retunred by
% a weibul function into dprime
% by dp = sqrt(2)*Z(PC);

dp = sqrt(2)*icdf('normal',weibullFun(x,params),0,1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   generateSyntheticData   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function exptData = generateSyntheticData(exptData)

% first get implied threshold contrast for
% the following dprime
testDprime = exptData.behavior.dprime;
testDeltaContrasts = exptData.fitcontrast;%0:0.001:1;
for pedNum = 1:length(exptData.pedestals)
 % first compute dprime as a function of delta contrast for
 % each pedstal
 for deltaContrastNum = 1:length(testDeltaContrasts)
  % get the contrast + delta contrast the trial was shown at
  pedestal = exptData.pedestals(pedNum);
  deltaContrast = testDeltaContrasts(deltaContrastNum);
  % get the response difference and sigma
  responseDiff = crf(pedestal+deltaContrast,exptData.crf)-crf(pedestal,exptData.crf);
  sigma = sqrt(crfsigma(pedestal+deltaContrast,exptData.noise,exptData.crffit)^2+crfsigma(pedestal,exptData.noise,exptData.crffit)^2);
  dprime(deltaContrastNum) = responseDiff/sigma;
 end
 % now use linear interpolation to get the delta contrast for
 % a given testDprime
 [dprime uniqueIndexes] = unique(dprime);
 exptData.syntheticThresholdContrast(pedNum) = interp1q(dprime,testDeltaContrasts(uniqueIndexes),testDprime);
 % now compute some contrasts limits that are between 1/2 and 2 of testDprime
 testMinContrast(pedNum) = interp1(dprime,testDeltaContrasts(uniqueIndexes),testDprime/2);
 testMaxContrast(pedNum) = interp1(dprime,testDeltaContrasts(uniqueIndexes),testDprime*2);
end

exptData.syntheticThresholdContrast(isnan(exptData.syntheticThresholdContrast)) = 1;
testMinContrast(isnan(testMinContrast)) = 1;
testMaxContrast(isnan(testMaxContrast)) = 1;

%disp(sprintf('[dprimefit3_r2:generateSyntheticData] Synthetic thresholds: [%s] - Noise: %1.4f',num2str(exptData.syntheticThresholdContrast),exptData.noise.k));

% now compute synthetic data, being careful to test delta contrasts
% near threshold
for trialNum = 1:exptData.n
 % get the contrast + delta contrast the trial was shown at
 pedestal = exptData.pedContrast(trialNum);
 
 % choose a delta contrast on the uniform interval between min and max test contrast
 pedNum = find(exptData.pedestals==pedestal);
 deltaContrast = testMinContrast(pedNum)+(testMaxContrast(pedNum)-testMinContrast(pedNum))*rand;
 exptData.deltaContrast(trialNum) = deltaContrast;
 %  deltaContrast = exptData.deltaContrast(trialNum);
 % get the response difference and sigma
 responseDiff = crf(pedestal+deltaContrast,exptData.crf)-crf(pedestal,exptData.crf);
 sigma = sqrt(crfsigma(pedestal+deltaContrast,exptData.noise,exptData.crffit)^2+crfsigma(pedestal,exptData.noise,exptData.crffit)^2);
 
 % compute the probability of getting the trial correct
 pcorrect(trialNum) = cdf('norm',responseDiff,0,sigma);
end

% now compute a randomized set of responses
exptData.correctIncorrect = rand(1,exptData.n)<pcorrect;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   modellike: to do noise fit   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function loglike = modellike(params,exptData,noisetype)

p = getNoiseParams(params,noisetype);

% for each trial compute the probability
% of being correct given the paramaters
% and the data
for iTrial = 1:exptData.n
 % get the pedestal contrast and the thresh contrast (i.e. pedestal+deltaC)
 % for this trial
 pedContrast = exptData.pedContrast(iTrial);
 threshContrast = pedContrast + exptData.deltaContrast(iTrial);
 
 % compute the response difference and standard deviation of response
 % responseDifference = mu and sigma is standard deviation
 responseDifference = crf(threshContrast,exptData.crf) - crf(pedContrast,exptData.crf);
 sigma = sqrt(crfsigma(threshContrast,p,exptData.crf)^2+crfsigma(pedContrast,p,exptData.crf)^2);
 
 % now compute probability of correct for this trial
 % NOTE this assumes that the difference of two response is
 % gaussian distributed. This is correct if the noise is
 % gaussian distributed. If the noise is *not* gaussian
 % distributed, the central limit theorem gets us close,
 % depending on what the actual distribution is. We need
 % to explore this issue later.
 pcorrect(iTrial) = cdf('norm',responseDifference,0,sigma);
end

% compute probabilty incorrect
pincorrect = 1-pcorrect;

% compute the probablity with which the model
% generates the observers correct and incorrect
% responses for each trial
pObserverCorrectTrials = pcorrect(find(exptData.correctIncorrect==1));
pObserverIncorrectTrials = pincorrect(find(exptData.correctIncorrect==0));

% now compute the loglikelihood of the model
% generating the observers responses
loglike = sum(log([pObserverCorrectTrials pObserverIncorrectTrials]));

% we are using fminsearch, so need to minimize
loglike = -loglike;


%%%%%%%%%%%%%%%%%%%%%%%%
%   getNoiseParamsStr  %
%%%%%%%%%%%%%%%%%%%%%%%%
function str = getNoiseParamsStr(noise)

% extract the params for the array
switch noise.noisetype
 case {'additive'}
  str = sprintf('additive; k=%0.4f',noise.k);
  
 case {'multiplicative'}
  str = sprintf('multiplicative; k=%0.4f',noise.k);
  
 case {'both'} % noise has additive and multiplicative components
  str = sprintf('both; k=%0.4f k1=%0.4f',noise.k,noise.k1);
  
 case {'square'} % noise is aditive and scales with response
  str = sprintf('square; k=%0.4f k1=%0.4f',noise.k,noise.k1);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   getNoiseParams: to do noise fit   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function p = getNoiseParams(params,noisetype)

% set the model type
p.noisetype = noisetype;

% extract the params for the array
switch noisetype
 case {'additive'}
  p.k = params(1);
  
 case {'multiplicative'}
  p.k = params(1);
  
 case {'both'} % noise has additive and multiplicative components
  p.k = params(1);
  p.k1 = params(2);
  
 case {'square'} % noise is aditive and scales with response
  p.k = params(1);
  p.k1 = params(2);
  p.k1 = params(3);
  
 otherwise
  %disp(sprintf('(dprimefit3_r2:getNoiseParams) Unknown model type %s',noisetype));
  keyboard
end


%%%%%%%%%%%%%%%%%%%%%%%
%   initNoiseParams   %
%%%%%%%%%%%%%%%%%%%%%%%
function initParams = initNoiseParams(noisetype,noise)

% init parameters of model
switch (noisetype)
 case {'additive' 'multiplicative'}
  initParams = .00001;
  
 case  {'both'}
  initParams = [noise.k-.001 noise.k1-.001];
  
 case {'square'}
  initParams = [noise.k-.001 noise.k1-.001 noise.k2-.001];
  
 otherwise
  %disp(sprintf('[dprimefit3_r2:computeNoiseValue] Unknown noisetype %s',noisetype));
  keyboard
end


%%%%%%%%%%%%%%%%%
%  fitNoiseDual %
%%%%%%%%%%%%%%%%%
function noise = fitNoiseDual(exptData,noisetype)

% dispay what we are doing
%disp(sprintf('[dprimefit3_r2] fitting noise - CRF: "%s:, noise type: "%s"',exptData.used_crf,noisetype))

% init the noise params
initParams = initNoiseParams(noisetype,exptData.noise);

% optimization parameters
maxIter = inf;
MaxFunEvals = inf;
optimParams = optimset('MaxIter',maxIter,'MaxFunEvals',MaxFunEvals,'Display','off','Diagnostics','off');

% use fminsearch to find function minimum
[fitParams maxLike exitflag fitOutput] = fminsearch(@modellike,initParams,optimParams,exptData,noisetype);

if exitflag
 noise = getNoiseParams(fitParams,noisetype);
else
 %disp(sprintf('[dprimefit3_r2:computeNoiseValue] Could not fit noise check ''modellike'''));
 noise = [];
 return
end

% display that we are done
%disp(sprintf('[dprimefit3_r2] noise fit DONE - CRF: "%s", noise type: %s value: "%s"',exptData.used_crf,noisetype,num2str(fitParams)))


%%%%%%%%%%%%%%%%%%%%%
%   fitCRFfromTvC   %
%%%%%%%%%%%%%%%%%%%%%
function noise = fitCRFfromTvC(exptData,noisetype)

disp('[dprimefit3_r2:fitCRFfromTvC] fitting noise: CRF from TvC');

% optimization parameters
maxIter = inf;
MaxFunEvals = inf;
optimParams = optimset('LevenbergMarquardt','on','MaxIter',maxIter,'MaxFunEvals',MaxFunEvals,'Display','off','Diagnostics','off');

% init the noise params
initParams = initNoiseParams(noisetype,[]);
minParams = zeros(1,length(initParams));
maxParams = exptData.maxk.*ones(1,length(initParams));

% add the offset param
initParams(end+1) = .2;
minParams(end+1) = exptData.minOffset;
maxParams(end+1) = inf;

% use lsqnonlin to find function minimum
contrast = exptData.pedestalsCRF;
response = exptData.use_crf;

[fitParams resnorm residual exitflag output lambda jacobian] = lsqnonlin(@fitCRFfromTvCResidual, ...
 initParams,minParams,maxParams,optimParams,contrast,response,noisetype,exptData);

% get the noise params and the offset from the params
noise = getNoiseParams(fitParams,noisetype);
noise.responseOffset = fitParams(end);

% if we have the best fit then keep it.
noise.bestFit.params = fitParams;

% compute R2 values and save them in the fit:
% (1) reshape the errors by condition and contrast
noise.bestFit.err = fitCRFfromTvCResidual(fitParams,contrast,response,noisetype,exptData);

% (2) run through conditions and compute r2
% for each one of them
thisError = noise.bestFit.err;
thisResponses = response;
noise.bestFit.r2 = 1-var(thisError)/var(thisResponses);

noise.bestFit.resnorm = resnorm;
noise.bestFit.residual = residual;
noise.bestFit.exitflag = exitflag;
noise.bestFit.output = output;
noise.bestFit.lambda = lambda;
noise.bestFit.jacobian = jacobian;

% computed a chi-suared test fo goodness of fit:
noise.bestFit.r2 = 1 - var(thisError)/var(thisResponses);
noise.bestFit.chi2.chi2value = sum((thisError./exptData.use_crf_ste).^2);
noise.bestFit.chi2.M = 2;
noise.bestFit.chi2.N = 8;
noise.bestFit.chi2.nparams = noise.bestFit.chi2.N - noise.bestFit.chi2.M;
noise.bestFit.chi2.dataVariability = exptData.use_crf_ste;
noise.bestFit.chi2.modelFitError = thisError;
noise.bestFit.chi2.p = 1 - gammainc(0.5*noise.bestFit.chi2.chi2value,0.5*noise.bestFit.chi2.nparams); 
% this last line is taken from numerical recepies

disp(sprintf('[dprimefit3_r2:fitCRFfromTvC] DONE fitting noise: CRF from TvC: R2 = [%s]',num2str(noise.bestFit.r2)));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   fitCRFfromTvCResiudal   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function residual = fitCRFfromTvCResidual(p,contrast,response,noisetype,exptData)

% get the noise params and the offset from the params
noise = getNoiseParams(p,noisetype);
responseOffset = p(end);

% disp(sprintf('Entered in fitCRFfromTvCResidual'))

% get the CRF
[modelContrast modelResponse] = makeModelCRFfromTvC(exptData,noise,responseOffset,0);

% now interpolate the model responses to the pedestal responses
pedC = (exptData.pedestalsCRF);
c = modelContrast;
r = modelResponse;

modelResponse = interp1(c,r,pedC,'pchip');

% compute residual
residual = response-modelResponse;

if any(isnan(modelResponse))
 keyboard
end


dispFit = 1;
if dispFit
 smartfig('dprimefit3_r_CRFfromTvC','reuse');
 cla
 semilogx(pedC,modelResponse,'ro-',pedC,response,'bo-')
 legend('CRF model', 'CRF data',4)
 ylabel('fMRI response')
 xlabel('Contrast')
 
end


%%%%%%%%%%%%%%%%%%%%
%  fitTvCfromCRF   %
%%%%%%%%%%%%%%%%%%%%
function noise = fitTvCfromCRF(exptData,noisetype)

% optimization parameters
optimParams = optimset('LevenbergMarquardt','on', ...
 'MaxIter',inf, ...
 'MaxFunEvals',inf, ...
 'Display','off', ...
 'Diagnostics','off', ...
 'TolFun',10^-10);

% init the noise params
initParams = initNoiseParams(noisetype,[]);
%minParams = zeros(1,length(initParams));
minParams = -inf(1,length(initParams));
maxParams = inf(1,length(initParams));

% use lsqnonlin to find function minimum
[fitParams resnorm residual exitflag output lambda jacobian] = lsqnonlin(@fitTvCfromCRFResidual,initParams,minParams,maxParams,optimParams,exptData.behavior.tvc,exptData.pedestals,exptData.crf,exptData.behavior.dprime,noisetype,exptData);

% get the noise params and the offset from the params
noise = getNoiseParams(fitParams,noisetype);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   fitTvCfromCRFResiudal   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function residual = fitTvCfromCRFResidual(p,threshold,contrast,crffit,dprime,noisetype,exptData)

% get the noise params and the offset from the params
noise = getNoiseParams(p,noisetype);

% get the CRF
[modelThreshold] = makeModelTvCfromCRF(contrast,noise,crffit,dprime,exptData,0);

% compute residual
residual = log(threshold)-log(modelThreshold);


%%%%%%%%%%%%%%%
%   skewFun   %
%%%%%%%%%%%%%%%
function y = skewFun(p,x)
% this functon fits the tvc
y = (p.peak*exp(-(log(x./p.mean)./(p.sd+p.skew.*log(x./p.mean))).^2)-exp(-1/p.skew^2));


%%%%%%%%%%%%%%%%%%%%
%   analyzeQuest   %
%%%%%%%%%%%%%%%%%%%%
function [results stimulus] = analyzeQuest(filename)
% function [results stimulus] = analyzeQuest(filename)
%
% this function analyzes quest in a data file
% ##change it to do the analysis on a pssed quest structure##

load(filename);

indexes = size(stimulus.quest.q);

if length(indexes)==2
 for hh = 1:indexes(1) % pedestal
  for jj = 1:indexes(2) % attentional condition
   if ~isempty(stimulus.quest.q{hh,jj})
    qthreshold(hh,jj) = QuestMean(stimulus.quest.q{hh,jj});
    qthSD(hh,jj) = QuestSd(stimulus.quest.q{hh,jj});
    bEstimate(hh,jj) = QuestBetaAnalysis(stimulus.quest.q{hh,jj});
   end
  end
 end
 % making errorbars:
 se = qthSD;
 s = 1;
else
 for hh = 1:indexes(1) % pedestal
  for jj = 1:indexes(2) % attentional condition
   for ii = 1:indexes(3) % number of quest
    if ~isempty(stimulus.quest.q{hh,jj,ii})
     qthreshold(hh,jj,ii) = QuestMean(stimulus.quest.q{hh,jj,ii});
     qthSD(hh,jj,ii) = QuestSd(stimulus.quest.q{hh,jj,ii});
     bEstimate(hh,jj,ii) = QuestBetaAnalysis(stimulus.quest.q{hh,jj,ii});
    end
   end
  end
 end
 % making errorbars:
 s = size(qthreshold);
 se = sum(qthSD.^2,3)./s(3);
 s = s(3);
end

results.meanThreshold = qthreshold;
results.thresholdSD = qthSD;
results.betaAnalysis = bEstimate;
results.pedestals = stimulus.pedestals;


%%%%%%%%%%%%%%%
%  fitexptest %
%%%%%%%%%%%%%%%
function bestfit = fitexptest(contrasts,responses, testType)

% check arguments
if (nargin ~= 3)
 help fitexptest
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
numConditions = [1 1]; % only one crf at the time

[initparams minfit maxfit] = switchInitialParamsEX(testType,numConditions);

% set optimization parameters
optimizationParams = optimset( ...
 'LevenbergMarquardt','on', ...
 'MaxIter',inf, ...
 'TolFun',10^-15, ...
 'MaxFunEvals',10000);


global numIters;numIters = 0;

% fit function using lsqnonlin in LevenbergMarquardt mode.
[fitparams resnorm residual exitflag output lambda jacobian] = ...
 lsqnonlin(@experr,initparams,minfit,maxfit,optimizationParams, ...
 contrasts,responses,testType,numConditions);

% Taken from Numerical Recipies,
% the leastsq function seems to return the transposed gradient
% instead of the jacobian...
hessian = jacobian'*jacobian;
reducedChiSquared = (residual'*residual)/(numel(responses)-length(initparams));
covar = (reducedChiSquared * inv(hessian));
M = max(diag(covar)); % get the maximum value on the
covar = covar/M; % normalized covariance matrix

% set params to the size of the conditions (e.g., 3*4*3):
[c50 n offset] = setExpParams(fitparams,testType,numConditions);

numIters = nan;
[bestfit.err bestfit.yfit] = experr(fitparams,contrasts,responses,testType,numConditions);

% compute R2 values and sasve them in the fit:
% (1) reshape the errors by condition and contrast
bestfit.err = reshape(bestfit.err,size(responses));

% (2) run through conditions and compute r2
% for each one of them
thisError = bestfit.err;
thisResponses = responses;
bestfit.r2 = 1-var(thisError)/var(thisResponses);

% if we have the best fit then keep it.
bestfit.params = fitparams;
bestfit.eachparam.c50 = c50;
bestfit.eachparam.n = n;
bestfit.eachparam.offset = offset;
bestfit.covar = covar;
bestfit.output = output;
bestfit.testType = testType;
bestfit.numConditions = numConditions;

% compute the function returned by the best fit:
x = .0001:.0001:1;
bestfit.fitx = zeros([numConditions length(x)]);
for i = 1:numConditions(1)
 for j = 1:numConditions(2)
  bestfit.fitx(i,j,:) = x;
 end
end

% [dummy bestfit.fity] = sigmoiderr(bestfit.params,bestfit.fitx,bestfit.fitx,testType);
% compute th current model estimate:
[bestfit.fity c50 n offset] = testexp(bestfit.params,bestfit.fitx,testType,numConditions);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function for fitting exponential to crf %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [err, fitfun] = experr(fitparams,contrasts,responses,testType,numConditions)
% computes the error from the model to the real data using the exponent
% function as model. not to be used anymore at riken.

% compute th current model estimate:
[fitfun c50 n offset] = testexp(fitparams,contrasts,testType,numConditions);

% calculate error:
err = responses-fitfun;
err = err(:);

% display the current fit (FIX: This doesn't yet work for displaying simultaneous fit)
global numIters
numIters = numIters+1;
displayAllFits = 0;

if displayAllFits
 % display the fit as we go along
 if ~isnan(numIters)
  % clear the figure:
  f = smartfig('dprimefit3_exp','reuse');
  h=get(f,'Child'); % handles to the axes
  for a = 1:length(h)
   cla(h(a))
  end
  % extract the current values:
  thisC = 100*squeeze(contrasts);
  thisR = squeeze(responses);
  thisFitFun = squeeze(fitfun);
  this_c50 = 100*c50;
  this_n = n;
  this_offset = offset;
  [sortcontrasts sortindex] = sort(thisC);
  
  semilogx(thisC,thisR,'ko');hold on
  semilogx(sortcontrasts,thisFitFun(sortindex),'r-');
  vline(this_c50);
  title(sprintf('c50=%0.2f,\n n=%0.2f,\n offset=%0.5f\n numIters=%i', ...
   this_c50,this_n,this_offset,numIters));
  axis([1 100 0 1.5])
  drawnow
 end
end


%%%%%%%%%%%
% testexp %
%%%%%%%%%%%
function [fitfun c50 n offset] = testexp(fitparams,contrasts,testType,numConditions)
% this functions uses an exponent instead of a naka-rushton to fit the crf.
% it should not be used anymore at riken.

% set params of the naka rushton depending on the type of test being
% performed:
[c50 n offset] = setExpParams(fitparams,testType,numConditions);

% calculate function
% retrieve parameters for this EXPONENTIAL
thisc50 = c50;
thisn = n;
thisoffset = offset;
thiscontrasts = squeeze(contrasts);
% now compute responses for this CRF
% fitfun = thisoffset+(thiscontrasts/thisc50).^thisn;
fitfun = thisoffset+(thiscontrasts.^thisn)/thisc50;


%%%%%%%%%%%%%%%%
% setExpParams %
%%%%%%%%%%%%%%%%
function [c50 n offset] = setExpParams(fitparams,testType,numConditions)
% this function sets the parameters for the exponential function
% it should not be used anymore.

switch testType
 case {'111'} % fit all params
  % and grab the values
  c50    = fitparams(1);
  n      = fitparams(2);
  offset = fitparams(3);
  
 otherwise
  keyboard
end


%%%%%%%%%%%%%%%%%%%%%%%%%
% switchInitialParamsEX %
%%%%%%%%%%%%%%%%%%%%%%%%%
function [initparams minfit maxfit] = switchInitialParamsEX(testType,numConditions)
% this function sets the initial parameters of the xponet that was used to
% fit the crf. not used anymore at riken.
% c50 = fitparams(1);
% n = fitparams(2);
% offset = fitparams(3);

switch testType
 case {'111'} % fit all params
  % c50
  initparams(1) = .1;
  minfit(1) = 0;
  maxfit(1) = inf;
  % n
  initparams(2) = .1;
  minfit(2) = 0;
  maxfit(2) = inf;
  % offset
  initparams(3) = .15;
  minfit(3) = 0;
  maxfit(3) = inf;
  
 otherwise
  keyboard
end


%%%%%%%%%%%%%%%%%%%
%  fitsigmoidtest %
%%%%%%%%%%%%%%%%%%%
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
numConditions = size(responses);    % first is adaptation second is attention
numConditions = numConditions(1:2); % the fit is run for adaptation and attention

[initparams minfit maxfit] = switchInitialParamsNK(testType,numConditions);

% set optimization parameters
optimizationParams = optimset(...
 'LevenbergMarquardt','on',   ...
 'MaxIter',inf,               ...
 'TolFun',10^-10,             ...
 'MaxFunEvals',10000);

global numIters;numIters = 0;

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

numIters = nan;
[bestfit.err bestfit.yfit] = sigmoiderr(fitparams,contrasts,responses,testType,numConditions);

% compute R2 values and save them in the fit:
% (1) reshape the errors by condition and contrast
bestfit.err = reshape(bestfit.err,size(responses));

% (2) run through conditions and compute r2
% for each one of them
thisError = bestfit.err;
thisResponses = responses;
bestfit.r2 = 1-var(thisError)/var(thisResponses);

% (3) compute the overal r2 of the whole model
thisError = bestfit.err(:);
thisResponses = responses(:);
bestfit.model_r2 = 1-var(thisError)/var(thisResponses);

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
x = .00001:.00001:1;
bestfit.fitx = x;

% compute th current model estimate:
[bestfit.fity Rmax c50 n offset] = testnaka(bestfit.params,bestfit.fitx,testType,numConditions);


%%%%%%%%%%%%%%
% sigmoiderr %
%%%%%%%%%%%%%%
function [err fitfun] = sigmoiderr(fitparams,contrasts,responses,testType,numConditions)

% compute th current model estimate:
[fitfun Rmax c50 n offset] = testnaka(fitparams,contrasts,testType,numConditions);

% calculate error:
err = responses-fitfun;
err = err(:);

% calculate residual sum of squares

% display the current fit (FIX: This doesn't yet work for displaying simultaneous fit)
global numIters 
displayAllFits = 1;
numIters = numIters+1;
if displayAllFits
 % display the fit as we go along
 if ~isnan(numIters)
  % extract the current values:
  thisC = 100*contrasts;
  thisR = responses;
  thisFitFun = fitfun;
  this_c50 = 100*c50;
  this_Rmax = Rmax;
  this_n = n;
  this_offset = offset;
  [sortcontrasts sortindex] = sort(thisC);

  smartfig('dprimefit3_r_nakarushton','reuse');
  cla
  semilogx(thisC,thisR,'bs');hold on
  semilogx(sortcontrasts,thisFitFun(sortindex),'r.-');
  vline(this_c50);
  title(sprintf('Rmax=%0.2f,\n c50=%0.2f,\n n=%0.2f,\n offset=%0.2f\n numIters=%i', ...
   this_Rmax,this_c50,this_n,this_offset,numIters));
  axis([1 100 0 1.5])
  drawnow
 end
end


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


%%%%%%%%%%%%%%%%%
% setNakaParams %
%%%%%%%%%%%%%%%%%
function [Rmax c50 n offset] = setNakaParams(fitparams,testType,numConditions)

switch testType
 case {'1111'} % fit all params
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
  % Rmax
  initparams(1) = .999;
  minfit(1) = .99;
  maxfit(1) = 1;
  % c50
  initparams(2) = .01;
  minfit(2) = 0;
  maxfit(2) = .2;
  % n
  initparams(3) = .7;
  minfit(3) = 0;
  maxfit(3) = 4;
  % offset
  initparams(4) = .15;
  minfit(4) = 0;
  maxfit(4) = .5;
  
 otherwise
  keyboard
end

