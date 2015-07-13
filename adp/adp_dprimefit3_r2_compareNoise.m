function [noise exptData] = adp_dprimefit3_r2_compareNoise(data,varargin)
% adp_dprimefit3_r2_compareNoise.m
%
%        $Id: adp_dprimefit3_r2_compareNoise.m
%      usage: adp_dprimefit3_r2_compareNoise(events,d1)
%         by: franco pestilli
%       date: 01/23/09
%    purpose: compute the contrast
%             threshold at a certain d' value for each adaptation condition
%
%  This function computes estimates of noise params
%  it works only for distributed target and attended.
%  the rest of the conditions do not have behavior.
%
% it is used to generate the plots of the test of the noise model:
% 1) get noise for distributed and use it to estimate the crf of
%    distributed
% 2) use the nosie for distributed to estimate the crf of attended
% 3) estimate noise for attended and use it to estimate the crf for attended
%
% all these tests are performed by makenoisestatisticsfigs2_r2_compareNoise.m
%
% the function can load either the individual observers' data or the
% avergae data.
%
% version 3_r uses the inv/deriv of a Naka-Rushton to do the fits
% _r stand for riken. this is done in the last stages of the paper.
%
% it is now called by: makenoisestatisticsfigs2_r2_compareNoise.m
%
% NB 2009.307.30 - skewfit and nk2 are implemented
%
% franco pestilli 2010.01.10

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% variables that define fit types.
crfFittype =[];    % Possible: {'linear' {'naka' 'Naka' 'NakaRushton' 'nk'}};
psychoFittype=[];  % Possible: {'mllfit' 'quest' 'average'}
dprimefittype=[];  % possible {'nakarushton' 'm'} {'poly' 'polynomial'} {'powerlaw' 'exp' 'exponent' 'power' 'e'} {'skewedgaussian' 'sg'}
fudgeIndex=[];     % this used to be used as follwos now it fixed to 1 -> events and d1 are assumed to be structures for each adaptation condition (1=adapt-0, 2=adapt-28, 3=adapt-100)
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
adaptation = [];   % This is the index of the adaptation condition fo rthe mean observer inside the behavioral data
conditionIndex=[]; % this is an index for plot info the first one is
minOffset=[];      % this value passed in sets the minimum offset value, when fitting distributed frst and
% attended after i use the distributed offset as min
% value for attended. i hope this solves the problems
% with negative noise reduction index values in fp-v3
maxk=[];           % same thing as min offset
doMean=[];         % if '1' sets up the average data, sets up individual observers otherwise.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set arguments
getArgs(varargin,{ ...
 'crfFittype=exp', ...
 'psychoFittype=mllfit', ...
 'dprimefittype=nk', ...
 'doNoise=1', ...
 'doBootstrap=0', ...
 'numBootstraps=1000', ...
 'dprime=1', ...
 'use0contrast=1', ...
 'whichCRF=crf_distributed_nontarget', ...
 'noisetype=additive', ...
 'noise=[]', ...
 'dispFit=0', ...
 'dataTitle=[]', ...
 'fudgeIndex=1', ...
 'VisualArea','v1', ...
 'adaptation',1, ...
 'figureInfo',[], ...
 'conditionIndex',[1 1], ...
 'minOffset',0, ...
 'maxk',0, ...
 'doMean',0},'verbose=0');

% check arguments
if nargin < 1
 help adp_dprimefit3_r2_compareNoise
 return
end

% do a large bootstrap only for one visual area:
if strcmpi(VisualArea,'v1')
 numBootstraps=numBootstraps;
else
 numBootstraps=2;
end

% switch between individual observers results or average results:
if ~doMean
 % setup exptData structure
 exptData = setDataStruct(data{fudgeIndex}.events,data{fudgeIndex}.d1,whichCRF, ...
  psychoFittype,simulateData,doBootstrap,numBootstraps,dprime,use0contrast, ...
  VisualArea,fudgeIndex,conditionIndex,minOffset,maxk);
 clear data;
 
 % (1) threshold is estimated from the data only for the individual observers,
 % not from the average
 % - estimate threshold contrast from data by using a maxLikelihood fit:
 exptData = computeXcontrast(exptData);
 
else
 exptData = setDataStructAVRG(data{fudgeIndex}.meanC,data{fudgeIndex}.d1,whichCRF,dprime,VisualArea,conditionIndex,minOffset,maxk,adaptation);
end

exptData.title = dataTitle;


exptData.crf = []; % NB this field must be set always because some of the calls of the noise model (e.g. multiplicative) require it

% (2) compute noise estimate given the crf and the TvC
switch doNoise
 case {'dual'}
  % not implemented in this function
  
 case {'tvc2crf'}
  % fit a tvc
  exptData = fitTVC(exptData,dprimefittype,dispFit>1);
  
  % fit the CRF from the TvC allowing noise parameters and offset to be adjusted
  if isempty(noise) % if noise is not empty it was passed and needs to be used to estimate the CRF
   noise = fitCRFfromTvC(exptData,noisetype);
   doRecomputeR2 = 0;

  else
   doRecomputeR2 = 1;
   disp(sprintf('[adp_dprimefit3_r2_compareNoise] noise was not empty, but passed through, we are not going to fit...'))
  end
  
  % display the fit
  [contrast response] = makeModelCRFfromTvC(exptData,noise,noise.responseOffset,dispFit,figureInfo);
   
  % substitute the bestFit in noise witht
  % the one given by the current model
  % this is done only for the noise test, when the noise is passed in from
  % distributed and used to predict attended
  if doRecomputeR2
   % (a) interpolate the responses to find responses at the pedestal contrast:
   model_response = interp1(contrast,response,exptData.pedestalsCRF,'pchip');
   
   
   % (b) the error of the mode; to the data
   % computed a chi-suared test fo goodness of fit:
   thisError = exptData.use_crf - model_response;
   thisResponses = exptData.use_crf;
   thisResponses_ste = exptData.use_crf_ste;
   
   noise.bestFit.responses_model = model_response;
   noise.bestFit.responses_data = response;
   
   SSreg = var(model_response);
   SStot = SSreg + var(thisError);
   
   noise.bestFit.resnorm  = [];
   noise.bestFit.residual = [];
   noise.bestFit.exitflag = [];
   noise.bestFit.output   = [];
   noise.bestFit.lambda   = [];
   noise.bestFit.jacobian = [];
   
   % computed a chi-suared test fo goodness of fit:
   noise.bestFit.r2 = SSreg/SStot;
   noise.bestFit.chi2.chi2value = sum((thisError./thisResponses_ste).^2);
   noise.bestFit.chi2.M = 2;
   noise.bestFit.chi2.N = 8;
   noise.bestFit.chi2.nparams = noise.bestFit.chi2.N - noise.bestFit.chi2.M;
   noise.bestFit.chi2.dataVariability = exptData.use_crf_ste;
   noise.bestFit.chi2.modelFitError = thisError;
   noise.bestFit.chi2.p = 1 - gammainc(0.5*noise.bestFit.chi2.chi2value,0.5*noise.bestFit.chi2.nparams); % this last line is taken from numerical recepies
   disp(sprintf('[adp_dprimefit3_r2_compareNoise] DONE computing the fit of the model given a passed noise: R2 = [%s]', ...
    num2str(noise.bestFit.r2)));
  end
  
  % save the output variables into:
  exptData.noise = noise; % save the whole noise structure here
  exptData.noise.response = response; % add the numerically reconstructed crf
  exptData.noise.contrast = contrast; % and contrast (the x-values for the crf)
  
  % distracters contrast inormation (this is used only for the decision model) 
  exptData.distractorContrast = {[0 .0175 .035] [0 .0175 .035] [.0175 .035 .07] [.035 .07 .14] ...
                                   [ .07 .14 .28] [ .14 .28 .56] [.28   .56  .84] [.28  .56 .84]};

  
 case {'crf2tvc'}
  % not implemented in this function
  keyboard
  
 otherwise
  keyboard
end

% %%%%%%%%%%%%%%%%%%%%% end main call %%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%
%   fitTVC   %
%%%%%%%%%%%%%%
function exptData = fitTVC(exptData,dprimefittype,dispFit)
% this function fits different type of TvC to the data:

% disp(sprintf('[adp_dprimefit3_r2_compareNoise:fitTVC] start, (%s)',dprimefittype));
% global figLabel;

% fit the crf with the appropriate function
switch dprimefittype
 case {'poly' 'polynonmial' 'pfit' 'p'}
  
  pedestals = exptData.pedestalsTvC;
  thresholds = exptData.behavior.tvc.thisTvC;
  
  % fit a polynomial:
  exptData.behavior.tvcfit.polyorder = 2;
  fitParams = fitpolynomial(pedestals,thresholds,exptData.behavior.tvcfit.polyorder);
  exptData.behavior.tvcfit.fitParams = fitParams;
  
 case {'skewedgaussian' 'sg' 'skew' 'gaussian' 'g'}
  % fit a a skewed gaussian on log axis:
  TVC = exptData.behavior.tvc.thisTvC(2:end);
  TVC = TVC(:);
  pedestals = exptData.pedestalsTvC(2:end);
  pedestals = pedestals(:);
  fitParams = skewfit(pedestals,TVC,dispFit);
  exptData.behavior.tvcfit.fitParams = fitParams;
  
  
 case {'nkn2' 'nkn' 'nakarushtonNumerical' 'nknum' 'NakaNum' 'nakanum'}
  
  % fit a naka-rushton sigmodal function:
  %             <n>    <thr>    <p>     <dr>    <scaleC50> <dipSlope>
  %   initParams = [1.6    0.0120    0.7    0.01      10       1];
  %   minParams  = [.2     0.00001   0.01   0.000001  1       .1];
  %   maxParams  = [3.5    0.85      1      0.2     100      3];
  
  initParams = [1.6    0.0120    0.9   0.001     50       1];
  minParams  = [.2    0.0001    0.01  0.000001   1       .1];
  maxParams  = [6      0.75      1       1       100      4];
  
  
  pedestals = exptData.pedestalsTvC;
  thresholds = exptData.behavior.tvc.thisTvC;
  thresholds_ste = exptData.behavior.tvc.thisTvCste;
  
  bestFit = fitnakarushtonTvCNumerical(pedestals,thresholds,thresholds_ste,initParams,minParams,maxParams);
  exptData.behavior.tvcfit.fitParams.fit = bestFit;
  
 otherwise
  %disp(sprintf('[adp_dprimefit3_r2_compareNoise:fitTVC] Unknown TvC psychoFittype: %s',dprimefittype));
  keyboard
end

% set the fittype
exptData.behavior.tvcfit.type = dprimefittype;

% make sure we can compute the thresholds properly
dispFit = 0;
if dispFit
 fitContrast = 0.008:0.001:.85;
 
 %disp(sprintf('[adp_dprimefit3_r2_compareNoise:fitTVC] showingthe fit: %s',dprimefittype));
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
disp(sprintf('[adp_dprimefit3_r2_compareNoise:fitTVC] done'));


%%%%%%%%%%%%%%%%%%%%%%%
%    fitpolynomial    %
%%%%%%%%%%%%%%%%%%%%%%%
function fitParams = fitpolynomial(x,y,polyorder)

% get the glm
x = x(:);y = y(:);
glm = getPolynomialGLM(x,polyorder);

% fit the polynomial
fitParams = pinv(glm)*y;

% display the fit
dispFit = 0;
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   getThresholdContrast   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function threshold = getThresholdContrast(pedestal,behavior)
% check that we have done the fit
if ~isfield(behavior,'tvcfit')
 %disp(sprintf('[adp_dprimefit3_r2_compareNoise] tvcfit has not been done'));
 keyboard
end

% return the threshold contrast for the given pedestal. Depends
% on which fit type we did.
switch (behavior.tvcfit.type)
 case {'poly' 'polynomial'}
  threshold = 0;
  for p = 0:behavior.tvcfit.polyorder
   threshold = threshold+behavior.tvcfit.fitParams(p+1)*pedestal^p;
  end
  
 case {'skewedgaussian' 'sg' 'skew' 'gaussian' 'g'}
  threshold = skewFun(behavior.tvcfit.fitParams.params,pedestal);
  
 case {'nkn2' 'nakarushtonNumerical' 'NakaNum' 'nknum' 'nakanum'}
  params = behavior.tvcfit.fitParams.fit.fitParams;
  threshold = nakarushtonTvCNumerical(pedestal,params);
 otherwise
  %disp(sprintf('[adp_dprimefit3_r2_compareNoise:getThresholdContrast) Unknown TvC fittype: %s',behavior.tvcfit.type));
  keyboard
end


%%%%%%%%%%%%%%%
%   skewFun   %
%%%%%%%%%%%%%%%
function y = skewFun(p,x)
% this functon fits the tvc
y = (p.peak*exp(-(log(x./p.mean)./(p.sd+p.skew.*log(x./p.mean))).^2)-exp(-1/p.skew^2));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    fitnakarushtonTvCNumerical    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
if size(t,1) > 1
 t = t';
end

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
bestFit.chi2.chi2value = sum((thisError./t_ste(:)').^2);
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%
% resNakarushtonNumerical %
%%%%%%%%%%%%%%%%%%%%%%%%%%%
function residual = resNakarushtonNumerical(params,pedestal,data)

model = nakarushtonTvCNumerical(pedestal,params);
residual = log10(data)-log10(model);


dispFit = 0;
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    makeModelCRFfromTvC    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [contrast response] = makeModelCRFfromTvC(exptData,noise,responseOffset,dispFit,figureInfo)

if nargin < 4,dispFit = 0;end
% disp(sprintf('[adp_dprimefit3_r2_compareNoise:makeModelCRFfromTvC] start'));

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
 
 %  disp(sprintf('[adp_dprimefit3_r2_compareNoise:makeModelCRFfromTvC] current DeltaC %s',num2str(deltaContrast)));
 
 % get the sigma
 sigma = sqrt(crfsigma(contrast(end),noise,exptData.crf)^2+crfsigma(contrast(end)+deltaContrast,noise,exptData.crf)^2);
 
 % now compute the response at the pedestal+delta contrast
 contrast(end+1) = contrast(end)+deltaContrast;
 response(end+1) = dprime*sigma+response(end);
 
end
% $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ END NUMERICAL $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ %%


%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    fitnakarushtonTvC    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
exptData.pedestalsCRF = exptData.pedestals + exptData.behavior.tvc.thisTvC/2;


%%%%%%%%%%%%%%%%%%%%%%%
%    setDataStruct    %
%%%%%%%%%%%%%%%%%%%%%%%
function exptData = setDataStruct(events,d1,whichCRF,psychoFittype,simulateData, doBootstrap,numBoostraps,dprime,use0contrast, VisualArea,a,conditionIndex,minOffset,maxk)

if strcmp(whichCRF,'crf_attended')
 indexAttention = 1;
  attCond = 'A';
elseif strcmp(whichCRF,'crf_distributed_target')
 indexAttention = 2;
 attCond = 'D';
elseif strcmp(whichCRF,'crf_distributed_nontarget')
 indexAttention = 2;
 attCond = [];
else
 keyboard
end


% figure formatting info:

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

exptData.plotInfo.XYColor        = [.4 .4 .4];
exptData.plotInfo.Fsize          = 8;
exptData.plotInfo.MarkerSize     = 10;
plotInfo.MarkerEdgeColor         = [.4 .4 .4];
exptData.plotInfo.LineWidth      = 1;
exptData.plotInfo.TickLength     = [0.025 .01];
exptData.plotInfo.PlotBoxAspectRatio = [1 1 1];
exptData.plotInfo.YLim           = [0 1.3];
exptData.plotInfo.Xlim           = [.01,1];
exptData.plotInfo.YTicksCRF      = [0 .25 .5 .75 1 1.25];
exptData.plotInfo.YTicksLabelCRF = [0 .25 .5 .75 1 1.25];
exptData.plotInfo.YTicksTvC      = [.00875 .0175 .035 .07 .14 .28];
exptData.plotInfo.YTicksLabelTvC = 100*[.00875 .0175 .035 .07 .14 .28];

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

% change the x-values of contrast for fits and plots
% depending on which crf we are using:
switch whichCRF
 case {'crf_distributed_target' 'crf_attended'}
  meanC = computedMeanC(events,.5);
  exptData.pedestalsCRF = meanC.(sprintf('quantile%scontrast',attCond))./100;
 case {'crf_distributed_nontarget' 'crf_unattended'}
  exptData.pedestalsCRF = exptData.pedestals;
end

exptData.fitcontrast = logspace(log10(.0001),log10(1),100);

% get contrast response function
exptData.crf_distributed_target    = d1.Dt_amplitude;
exptData.crf_attended              = d1.A_amplitude;
exptData.crf_distributed_nontarget = d1.Dnt_amplitude;
exptData.crf_unattended            = d1.U_amplitude;

% and standard error
exptData.crf_distributed_target_ste    = d1.Dt_ste;
exptData.crf_attended_ste              = d1.A_ste;
exptData.crf_distributed_nontarget_ste = d1.Dnt_ste;
exptData.crf_unattended_ste            = d1.U_ste;

% used crf and ste:
exptData.use_crf = exptData.(whichCRF);
exptData.use_crf_ste = exptData.([whichCRF,'_ste']);
exptData.used_crf = whichCRF;

% set the minimum value for noise offset during fit.
% for attended this value is set to be the value for distrbuted
exptData.minOffset = minOffset;
exptData.maxk = maxk;

exptData.distractorContrast = {[0 .0175 .035] [0 .0175 .035] [.0175 .035 .07] [.035 .07 .14] ...
                               [ .07 .14 .28] [ .14 .28 .56] [.28   .56  .84] [.28  .56 .84]};


%%%%%%%%%%%%%%%%%%%%%
% setDataStructAVRG %
%%%%%%%%%%%%%%%%%%%%%
function exptData = setDataStructAVRG(meanC,d,whichCRF,dprime,VisualArea,conditionIndex,minOffset,maxk,adaptation)
% this function set up the structure for the average:

if strcmp(whichCRF,'crf_attended')
 indexAttention = 1;
 attCond = 'A';
 icrf = 2;
elseif strcmp(whichCRF,'crf_distributed_nontarget') % this was originally fit to distributed target then we changeed 
                                                    % it because it makes sense to use the non-target 
                                                    % given that we are saying that distributed target has an extra 
                                                    % response due to memory
 indexAttention = 2;
 attCond = 'D';
 icrf = 1;
elseif strcmp(whichCRF,'crf_distributed_target')
 indexAttention = 2;
 attCond = 'D';
 icrf = 1;
else
 keyboard
end

% figure formatting info:

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
exptData.plotInfo.MarkerEdgeColor     = [.4 .4 .4];
exptData.plotInfo.LineWidth  = 1;
exptData.plotInfo.TickLength = [0.025 .01];
exptData.plotInfo.PlotBoxAspectRatio = [1 1 1];
exptData.plotInfo.YLim       = [0 1.3];
exptData.plotInfo.Xlim       = [.01,1];
exptData.plotInfo.YTicksCRF  = [0 .25 .5 .75 1 1.25];
exptData.plotInfo.YTicksLabelCRF = [0 .25 .5 .75 1 1.25];
exptData.plotInfo.YTicksTvC  = [.00875 .0175 .035 .07 .14 .28];
exptData.plotInfo.YTicksLabelTvC = 100*[.00875 .0175 .035 .07 .14 .28];

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
% load the file containing the average TvC:
 filename = '/dev/data/crfaa/adaptation_ms/data/12-Jul-2009averageTvC.mat';

% filename = '/Volumes/data/riken/crfaa/fmridata/crfaa/data_used_files/12-Jul-2009averageTvC.mat';
exptData.behavior.tvc = load(sprintf('%s',filename),'meanTvC', 'steTvC');
exptData.behavior.tvc.thisTvC    = squeeze(exptData.behavior.tvc.meanTvC(adaptation,icrf,:));
exptData.behavior.tvc.thisTvCste = squeeze(exptData.behavior.tvc.steTvC(adaptation,icrf,:));

exptData.behavior.dprime = dprime;% set the d-prime value for the threshold wanted:

% get pedestal contrasts
exptData.pedestals    = [0      .0175 .035 .07 .14 .28 .56 .84];
exptData.pedestalsTvC = [.00875 .0175 .035 .07 .14 .28 .56 .84]; % used to compute and plot the tvc

% change the x-values of contrast for fits and plots
% depending on which crf we are using:
switch whichCRF
 case {'crf_distributed_target' 'crf_attended'}
  exptData.pedestalsCRF = meanC.(sprintf('quantile%scontrast',attCond))./100;
 case {'crf_distributed_nontarget' 'crf_unattended'}
  exptData.pedestalsCRF = exptData.pedestals;
end

exptData.fitcontrast = logspace(log10(exptData.pedestalsTvC(1)),log10(exptData.pedestalsTvC(end)),100);

% get contrast response function
exptData.crf_distributed_target    = d.Dt_amplitude;
exptData.crf_attended              = d.A_amplitude;
exptData.crf_distributed_nontarget = d.Dnt_amplitude;
exptData.crf_unattended            = d.U_amplitude;

% and standard error
exptData.crf_distributed_target_ste    = d.Dt_ste;
exptData.crf_attended_ste              = d.A_ste;
exptData.crf_distributed_nontarget_ste = d.Dnt_ste;
exptData.crf_unattended_ste            = d.U_ste;

% used crf and ste:
exptData.use_crf = exptData.(whichCRF);
exptData.use_crf_ste = exptData.([whichCRF,'_ste']);
exptData.used_crf = whichCRF;

% set the minimum value for noise offset during fit.
% for attended this value is set to be the value for distrbuted
exptData.minOffset = minOffset;
exptData.maxk = maxk;

exptData.distractorContrast = {[0 .0175 .035] [0 .0175 .035] [.0175 .035 .07] [.035 .07 .14] ...
                               [ .07 .14 .28] [ .14 .28 .56] [.28   .56  .84] [.28  .56 .84]};


%%%%%%%%
% dip2 %
%%%%%%%%
function y = dip2(c,T,beta)
% this one fixes the decrease in t at high contrast
% parameters of the copressive function at high contrast
% T center (threshold) of the function
% beta slope of the function
y = .015+exp(-(c./T).^beta-.025);


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
  %disp(sprintf('[adp_dprimefit3_r2_compareNoise:crf] Unknown crf fittype: %s',crffit.fittype));
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
  %disp(sprintf('[adp_dprimefit3_r2_compareNoise:crfsigma] Unknown model type %s',p.noisetype));
end


%%%%%%%%%%%%%%%%%%
%   computeTvC   %
%%%%%%%%%%%%%%%%%%
function exptData = computeTvC(exptData,psychoFittype)
% this function computes a tvc
%disp(sprintf('[adp_dprimefit3_r2_compareNoise:computeTvC] start...'));

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
  keyboard
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
     disp(sprintf('[adp_dprimefit3_r2_compareNoise] Pedestal <%s> - percent bootstrap done <%s>',num2str(whichPedc),num2str(100*bootStrapN/numBoostraps)))
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
    prc = .995; % percentile used to compute errorbars (=/-95th percentile = two-tailed t-test)
    exptData.behavior.tvc.thisTvC(pedc) = th;
    exptData.behavior.tvc.thisTvCste(pedc) = 10.^(log10(th) + log10(prctile(btThreshold,prc)));% i think this is now fixed
    exptData.behavior.tvc.thisTvCste_percentileUsed = prc;
    exptData.behavior.tvc.thistvcSD = [];
    exptData.behavior.tvc.BootMedian(pedc) =  prctile(btThreshold,.5);
    exptData.behavior.tvc.BootStrapWBLparams{pedc} = bootstrpWBLfitp;
    exptData.behavior.tvc.fitresults{pedc}.wblparams = wblfitp;
    exptData.behavior.tvc.fitresults{pedc}.maxMLL = maxMLL;
    exptData.behavior.tvc.fitresults{pedc}.exitflag = exitflag;
    exptData.behavior.tvc.fitresults{pedc}.fitOutput = fitOutput;
   else
    %disp(sprintf('[adp_dprimefit3_r2_compareNoise:computeTvC] no bootstrap, pedestal %0.4f',whichPedc))
    % use fminsearch to find function minimum
    % init parameters for threshold fit:
    thresholdInit = median(deltaContrast); % mll weight the number of occurrences the median is a good start
    if isfield(exptData,'syntheticThresholdContrast')
     thresholdInit = exptData.syntheticThresholdContrast(pedc);
     % disp(sprintf('[adp_dprimefit3_r2_compareNoise] Using synthetic threshold of %0.4f as initial threshold',exptData.syntheticThresholdContrast(pedc)));
    else
     % disp(sprintf('[adp_dprimefit3_r2_compareNoise] Using %0.4f as initial threshold',thresholdInit));
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
% if (p(3) < 0) || (p(3) > .075) || (p(2) < 1 || p(2) > 5) || (p(1) < 0 || % NB ths is the old set of limits before the problem with JG p(1) > .575)
if (p(3) < 0) || (p(3) > .075) || (p(2) < 1 || p(2) > 6.5) || (p(1) < 0 || p(1) > .28575)
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
if length(yvals) > 1
x = interp1(yvals,xvals(i),y,'pchip');
else
 x = nan;
end
 

%%%%%%%%%%%%%%%%%%%
%    wbl2dprime   %
%%%%%%%%%%%%%%%%%%%
function dp = wbl2dprime(x,params)
% transform pcorrect retunred by
% a weibul function into dprime
% by dp = sqrt(2)*Z(PC);

dp = sqrt(2)*icdf('normal',weibullFun(x,params),0,1);


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
  %disp(sprintf('(adp_dprimefit3_r2_compareNoise:getNoiseParams) Unknown model type %s',noisetype));
  keyboard
end


%%%%%%%%%%%%%%%%%%%%%%%
%   initNoiseParams   %
%%%%%%%%%%%%%%%%%%%%%%%
function initParams = initNoiseParams(noisetype,noise)

% init parameters of model
switch (noisetype)
 case {'additive' 'multiplicative'}
  initParams = .0001;%noise.k-.001;
  
 case  {'both'}
  initParams = [noise.k-.001 noise.k1-.001];
  
 case {'square'}
  initParams = [noise.k-.001 noise.k1-.001 noise.k2-.001];
  
 otherwise
  %disp(sprintf('[adp_dprimefit3_r2_compareNoise:computeNoiseValue] Unknown noisetype %s',noisetype));
  keyboard
end


%%%%%%%%%%%%%%%%%%%%%
%   fitCRFfromTvC   %
%%%%%%%%%%%%%%%%%%%%%
function noise = fitCRFfromTvC(exptData,noisetype)
disp('[adp_dprimefit3_r2_compareNoise:fitCRFfromTvC] fitting noise: CRF from TvC');

% optimization parameters
maxIter = 1000000000000;
MaxFunEvals = 1000000000000000;
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
noise.bestFit.err = residual;%fitCRFfromTvCResidual(fitParams,contrast,response,noisetype,exptData);

% compute the model response
[model_contrast M_response] = makeModelCRFfromTvC(exptData,noise,noise.responseOffset);

% (2) run through conditions and compute r2
% for each one of them
thisError = noise.bestFit.err;
thisResponses = response;

noise.bestFit.responses_model = M_response;
noise.bestFit.responses_data = response;

noise.bestFit.resnorm = resnorm;
noise.bestFit.residual = residual;
noise.bestFit.exitflag = exitflag;
noise.bestFit.output = output;
noise.bestFit.lambda = lambda;
noise.bestFit.jacobian = jacobian;

% computed a chi-suared test fo goodness of fit:
% this was take from wikipedia.
% http://en.wikipedia.org/wiki/Coefficient_of_determination
SSreg = var(M_response);        % explained variance
SStot = SSreg + var(thisError); % total variance this was taken from wikipedia, 
noise.bestFit.r2 = SSreg/SStot;
% NB we are not using the chi2 values for the noise model anymore, these
% calculations might be slightly imprecise...
noise.bestFit.chi2.chi2value = sum((thisError./exptData.use_crf_ste).^2);
noise.bestFit.chi2.M = 2;
noise.bestFit.chi2.N = 8;
noise.bestFit.chi2.nparams = noise.bestFit.chi2.N - noise.bestFit.chi2.M;
noise.bestFit.chi2.dataVariability = exptData.use_crf_ste;
noise.bestFit.chi2.modelFitError = thisError;
noise.bestFit.chi2.p = 1 - gammainc(0.5*noise.bestFit.chi2.nparams,0.5*noise.bestFit.chi2.chi2value); % this last line is taken from numerical recepies

disp(sprintf('[adp_dprimefit3_r2_compareNoise:fitCRFfromTvC] DONE fitting noise: CRF from TvC: R2 = [%s]',num2str(noise.bestFit.r2)));


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


dispFit = 0;
if dispFit
 smartfig('dprimefit3_r_CRFfromTvC','reuse');
 cla
 semilogx(pedC,modelResponse,'ro-',pedC,response,'bo-')
 legend('CRF model', 'CRF data',4)
 ylabel('fMRI response')
 xlabel('Contrast')
 
end


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









