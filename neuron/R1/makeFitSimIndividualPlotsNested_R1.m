function r = makeFitSimIndividualPlotsNested_R1(data,meanC,varargin)
% function r = makeFitSimIndividualPlotsNested_R1(data,meanC,varargin)
%
% this function is called by runNestedHypothesisTest_r.m.
%
% it gets some data{adaptation} and contrasts and fits different models
% to the resulting contrast response functions.
% 
% - the 6 conditions version works to do a test on the correct/incorrect and 
% 1st/2nd interval anlysis
%
% varargin (Default):
% {'adaptation=0', 'testType=''1111''','dispFit=1', 'baselineParams=[]', ...
%  'crfFittype=naka','numBootstraps=200','use0contrast=0','dataTitle=[]', 'doBootstrap=0', ...
%  'subplotInfo',[1 1 1]}%
%
% franco pestilli 2009.07.24

% check arguments
if nargin < 2
 help makeFitSimIndividualPlotsNested_R1
 return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% variables that define fit types.
global displayAllFits
testType    = [];    % which params to test
crfFittype  = [];    % Possible: {'exp' 'naka' 'Naka' 'nk'};
visualArea  = [];    % which visual area are we wroking on (used as figure title)
dispFit     = [];    % 1 or 0, to display the fit (set to 2 for subsidary fits like TvC etc)
dataTitle   = [];    % arbitrary string to use as data title (should be subject data filename)
displayAllFits = []; % this is to display the fits as the procedure advances
saveFig        = []; % it saves the individual fits in a nice looking fashion
saveFitFig     = []; % it saves a figure with a panel fo each crf fit, this is used mostly to debug the fit procedure
observer       = []; % this is used to construct the filename for the individual figures
r2             = []; % r2 of the roi used
subtractionType= []; % type of subtraction for the 0-response
analysisType   = []; % passes in which analys we are using (standard or correct/incorrect or 1st/2nd interval, 
                     % so to created a folder name for the figures)
adaptation     = []; % 0, 28 or 100
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set arguments
[argName argValue] = getArgs(varargin,{ ...
 'testType','1111', ...
 'dispFit',0, ...
 'crfFittype','nk2', ...
 'dataTitle',[], ...
 'visualArea','v1', ...
 'displayAllFits',0, ...
 'saveFig',0, ...
 'saveFitFig', [], ...
 'observer','NS', ...
 'r2', '06',    ...
 'subtractionType', 'ORIG', ...
 'analysisType',analysisType, ...
 'adaptation',0});

% make the figure folder and data folder global to pass them along each
% function that needs them.
global figDir figDirEach defaultDataFolder;

if strcmpi(visualArea,'v1') && strcmpi(r2,'05')
disp('did not compute roi v1 at r2:0.05, exiting.');
return
end

%%%%%%%%%%%%%%%%%%%%%%%%
% set up save folders: %
%%%%%%%%%%%%%%%%%%%%%%%%
% all fits and results are saved in here:
% defaultDataFolder = '/data2/crfaa/crfaa/nestedHypothesisTest_R1/';
figDir = sprintf('figs_nestedFit_eachFit_%s_%s_%s',analysisType,r2,subtractionType);%['figs_nestedFit_eachFit_',analysisType]; % figures with individual function fits are saved here
figDirEach = sprintf('figs_nestedFit_individual_%s_%s_%s',analysisType,r2,subtractionType); % fits are saved here

% if this var was not passed then do not display the fit figure.
if isempty(saveFitFig)
 saveFitFig.save = 0;
 saveFitFig.observer   = 'nsO';
 saveFitFig.testType   = 'nsT';
 saveFitFig.visualarea = 'nsV';
end

if (length(testType) > 3) && (strcmp(crfFittype,'exp')) || ...
  (length(testType) < 4) && (strcmp(crfFittype,'naka'))
 disp(sprintf('[%s] number of parameters and test type do not match.',mfilename))
 keyboard
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ste up the data structure to pass around the local function %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
d = setDataStructure(data, meanC, visualArea,analysisType);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% do the simulataneous fit: %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
r = fitCRF(d,crfFittype, testType,saveFitFig);

%%%%%%%%%%%%%%%%%%
% save call info %
%%%%%%%%%%%%%%%%%%
r.callInfo.function.mfilename  = mfilename('fullpath');
r.callInfo.function.dbstack    = dbstack;
r.callInfo.function.varargin   = varargin;
r.callInfo.function.argsName   = argName;
r.callInfo.function.argsValue  = argValue;

r.callInfo.whichCRF       = 'all';     % 'dnt' 'dt' 'a' 'u';
r.callInfo.testType       = testType;  % e.g., '01000'
r.callInfo.crfFittype     = crfFittype;% Possible: {'nk2'};
r.callInfo.dispFit        = dispFit;   % 1 or 0, to display the fit (set to 2 for subsidary fits like TvC etc)
r.callInfo.dataTitle      = dataTitle; % arbitrary string to use as data title (should be subject data filename)
r.d = d;

%%%%%%%%%%%%%%%%%%%
% display the fit %
%%%%%%%%%%%%%%%%%%%
if dispFit
 p = dispCRF(r.d,r.crf,testType,saveFig,observer,visualArea,r2,subtractionType,adaptation);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        dispCRF         %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
function r = dispCRF(r,crf,testType,saveFig,observer,visualArea,r2,subtractionType,adaptation)
global defaultDataFolder figDirEach

if adaptation == 0
 i = 1;
elseif adaptation == 28
 i = 2;
elseif adaptation == 100
 i = 3;
else
 keyboard
end
 
% bring up a new figure for each adaptation condition and do the plot on it:
figName = sprintf('%scrf_simFit_%sAdapt%i_%s_%s_%s', ...
 upper(observer),testType,adaptation,visualArea,r2,subtractionType);

h = smartfig(figName,'reuse');

% to be clear get all the contrasts, crfs and crf_ste for the current
% adaptation condition:
pedestals  = squeeze(r.pedestals);
thisCRF    = squeeze(r.use_crf);
thisCRFste = squeeze(r.use_crf_ste);

% now run through each crf and plot it:
for j = 1:size(r.pedestals,1) % for each attention condition
 
 % the conditions without a tvc have 0 as first
 % pedestal, which is difficult to plot on log
 if pedestals(j,1) == 0
  pedestals(j,1) = r.plotInfo.Xlim(1);
 end
 
 % to be clear take out the current vectors to plot:
 plotTheseContrasts = squeeze(pedestals(j,:));
 plotTheseCRF = squeeze(thisCRF(j,:));
 plotTheseCRFste = squeeze(thisCRFste(j,:));
 
 % display the data
 myerrorbar(plotTheseContrasts,plotTheseCRF, ...
  'yError',plotTheseCRFste, ...
  'Symbol',r.plotInfo.thisSymbol{j}(1),...
  'MarkerFaceColor',r.plotInfo.MarkerFaceColor{i}{j}, ...
  'MarkerEdgeColor',r.plotInfo.MarkerEdgeColor,...
  'MarkerSize',r.plotInfo.MarkerSize, ...
  'Color',r.plotInfo.thisColor{i}{j});
 
 % display the smooth fit function:
 thisSmoothX = squeeze(crf.fit.fitx(j,:));
 index = find(thisSmoothX > r.plotInfo.plotPeds(1));
 thisSmoothX = thisSmoothX(index:end);
 thisSmoothY = squeeze(crf.fit.fity(j,(index:end)));
 hold on
 semilogy(thisSmoothX,thisSmoothY,...
  r.plotInfo.thisSymbol{j}(2:end), ...
  'Color',r.plotInfo.thisColor{i}{j}, ...
  'LineWidth', r.plotInfo.LineWidth);
end

% axis formatting:
set(gca,...
 'XScale','log', ...
 'YScale','lin', ...
 'XTick', [0.00109375 .0021875  .004375 .00875 .0175 .035 .07 .14 .28 .56 .84 1], ...
 'XTickLabel', [0 .0021875  .004375 .00875 .0175 .035 .07 .14 .28 .56 .84 1], ...
 'YTick', [0 .25 .5 .75 1 1.25] ,...
 'XLim',[0.00109375  1], ...
 'YLim',[0 1.3]);
axis('square')
myaxisCRF;

set(h,'PaperPosition',[.25 .25 8 10.5]);
set(h,'PaperOrientation','Portrait');
drawnow

if saveFig
% % % % % % %  defaultDataFolder = '/Users/frakkopesto/data/riken/crfaa/fmridata/crfaa/nestedHypothesisTest/';
% % % % % % %  figDirEach = 'figs_nestedFit_individual';
 figDir = savefig(h,'figName',figName,'defaultDataFolder',defaultDataFolder,'figDir',figDirEach,'verbose',1);
end


%%%%%%%%%%%%%%%
%   fitCRF   %%
%%%%%%%%%%%%%%%
function r = fitCRF(d,crfFittype, testType,saveFitFig)
disp(sprintf('[makeFitSimIndividualPlotsNested_R1:fitCRF] Using function: <%s>,  with test: <%s> start...',crfFittype,testType,saveFitFig.testType));

switch crfFittype
 case {'naka' 'NAKA' 'nk'}
  % fit a naka-rushton sigmodal function:
  % not implemented for the attention paper
  keyboard
  r.crf.fit = fitsigmoidtest(d.pedestals,d.use_crf, testType,saveFitFig);
  
 case {'naka2' 'NAKA2' 'nk2'}
  % fit a naka-rushton sigmodal function with an exponent at the numerator larger than that at the denominator:
  r.crf.fit = fitsigmoidtest2(d.pedestals,d.use_crf, d.use_crf_ste, testType,saveFitFig);
  
 case {'exp' 'EXP' 'e'}
  % this is the fuction i need to fit an exponential
  % not implemented for the attention paper
  keyboard
  r.crf.fit = fitexptest(d.pedestals,d.use_crf, testType,saveFitFig);
  
 otherwise
  disp(sprintf('[makeFitSimIndividualPlotsNested_R1:fitCRF] %s fit type not defined',crfFittype))
  keyboard
end

disp(sprintf('[makeFitSimIndividualPlotsNested_R1:fitCRF] Using function: <%s>,  with test: <%s> done...',crfFittype,saveFitFig.testType));


%%%%%%%%%%%%%%%%%%%
% fitsigmoidtest2 %
%%%%%%%%%%%%%%%%%%%
function bestfit = fitsigmoidtest2(contrasts,responses, ste, testType,saveFitFig)

% check arguments
if (nargin ~= 5)
 help fitsigmoidtest2
 bestfit = nan;
 keyboard
end

% make sure we have a column vector
if (size(contrasts,1) == 1)
 contrasts = contrasts';
end
if (size(responses,1) == 1)
 responses = responses';
end

% set initial params
numConditions = size(responses,1);

[initparams minfit maxfit] = switchInitialParamsNK2(testType,numConditions);

% set optimization parameters
optimizationParams = optimset( ...
 'LevenbergMarquardt','on', ...
 'MaxIter',10^3, ...
 'TolFun',10^-3, ...
 'MaxFunEvals', 10^3, ...
 'TolX',10^-3);


global numIters;numIters = 0;

originalSaveFitFig = saveFitFig.save;
saveFitFig.save = 0;
saveFitFig.displayAllFits=0;%fix fix change this later now only displays the last fit
% fit function using lsqnonlin in LevenbergMarquardt mode.
[fitparams resnorm residual exitflag output lambda jacobian] = ...
 lsqnonlin(@sigmoiderr2,initparams,minfit,maxfit,optimizationParams, ...
 contrasts,responses,ste,testType,numConditions,saveFitFig);

% Taken from Numerical Recipies,
% the leastsq function seems to return the transposed gradient
% instead of the jacobian...
hessian = jacobian'*jacobian;
reducedChiSquared = (residual'*residual)/(numel(responses)-length(initparams));
covar = (reducedChiSquared * inv(hessian));
M = max(diag(covar)); % get the maximum value on the
covar = covar/M; % normalized covariance matrix

% set params to the size of the conditions (e.g., 3*4*4):
[Rmax c50 n q offset] = setNakaParams2(fitparams,testType,numConditions);

% fix when one fo the ste is '0'
ste(ste==0) = .000001;

saveFitFig.save = originalSaveFitFig;
saveFitFig.displayAllFits=1;%fix fix change this later now only displays the last fit
[bestfit.err bestfit.yfit] = sigmoiderr2(fitparams,contrasts,responses,ste, testType,numConditions,saveFitFig);

% compute R2 values and save them in the fit:
% (1) reshape the errors by condition and contrast
bestfit.err = reshape(bestfit.err,size(responses));

% (2) run through conditions and compute r2
% for each one of them
for j = 1:size(bestfit.err,1) % attention conditions
 thisError = squeeze(bestfit.err(j,:));
 thisResponses = squeeze(responses(j,:));
 bestfit.r2(j) = 1-var(thisError)/var(thisResponses);
end

% (3) compute the overal r2 of the whole model
thisError = bestfit.err(:);
thisResponses = responses(:);
bestfit.model_r2 = 1-var(thisError)/var(thisResponses);

% if we have the best fit then keep it.
bestfit.params = fitparams;
bestfit.eachparam.Rmax = Rmax;
bestfit.eachparam.c50 = c50;
bestfit.eachparam.n = n;
bestfit.eachparam.q = q;
bestfit.eachparam.offset = offset;
bestfit.covar = covar;
bestfit.output = output;
bestfit.testType = testType;
bestfit.numConditions = numConditions;

% (4) compute  chi-squared test of goodness of fit.
% taken from numerical receipts, p. 
bestfit.chi2.chi2value = sum((thisError./ste(:)).^2);
bestfit.chi2.M = length(fitparams);
bestfit.chi2.N = numel(responses);
bestfit.chi2.nparams = bestfit.chi2.N - bestfit.chi2.M;
bestfit.chi2.dataVariability = ste(:);
bestfit.chi2.modelFitError = thisError;
bestfit.chi2.p = 1 - gammainc(0.5*bestfit.chi2.chi2value,0.5*bestfit.chi2.nparams);

% (5) compute the function returned by the best fit:
x = .0001:.0001:1;
bestfit.fitx = zeros([numConditions length(x)]);
for i = 1:numConditions
  bestfit.fitx(i,:) = x;
end

% compute th current model estimate:
[bestfit.fity Rmax c50 n q offset] = testnaka2(bestfit.params,bestfit.fitx,testType,numConditions);


%%%%%%%%%%%%%%%
% sigmoiderr2 %
%%%%%%%%%%%%%%%
function [err fitfun] = sigmoiderr2(fitparams,contrasts,responses,ste, testType,numConditions,saveFitFig)

% compute th current model estimate:
[fitfun Rmax c50 n q offset] = testnaka2(fitparams,contrasts,testType,numConditions);


% calculate error:

err = (responses-fitfun);
err = err(:);

% this would be a chi2 calculation (numerialc receips 14.1.5):
% err = (responses-fitfun)./ste;

% display the current fit (FIX: This doesn't yet work for displaying simultaneous fit)
global numIters figDir defaultDataFolder;

numIters = numIters+1;
if saveFitFig.displayAllFits
 figName = sprintf('show_all_fitNested_OBS%s_TEST%s_VISAR%s',saveFitFig.observer,saveFitFig.testType,saveFitFig.visualarea);
 f = smartfig(figName,'reuse');
 colors = {'r' 'b' 'm' 'g' 'y' 'c'};
 
 % display the fit as we go along
 if ~isnan(numIters)
  h = get(gcf,'Child'); % handles to the axes
  if numel(h) > 1
   for a = size(responses,1):-1:1
    cla(h(a))
   end
  end
  
  for j = 1:size(responses,1) % attention
   % extract the current values:
   thisC = 100*squeeze(contrasts(j,:));
   thisR = squeeze(responses(j,:));
   thisFitFun = squeeze(fitfun(j,:));
   this_c50 = 100*c50(j);
   this_Rmax = Rmax(j);
   this_n = n(j);
   this_q = q(j);
   this_offset = offset(j);
   [sortcontrasts sortindex] = sort(thisC);
   
   subplot(1,size(responses,1),j),
   semilogx(thisC,thisR,'ko');hold on
   semilogx(sortcontrasts,thisFitFun(sortindex),[colors{j},'-']);
   vline(this_c50);
   title(sprintf('Rmax=%0.2f,\n c50=%0.2f,\n n=%0.2f,\n q=%0.2f,\n offset=%0.2f\n numIters=%i', ...
    this_Rmax,this_c50,this_n,this_q,this_offset,numIters));
   axis([1 100 0 1.3])
   drawnow
  end
  
  % clear the figure:
  if saveFitFig.save
   set(f,'PaperPosition',[.25 .25 8 10.5]);
   set(f,'PaperOrientation','Portrait');
% figDir = 'figs_nestedFit_eachFit';
% defaultDataFolder = '/Users/frakkopesto/data/riken/crfaa/fmridata/crfaa/nestedHypothesisTest/';
   savefig(f,'defaultDataFolder',defaultDataFolder,'figName',figName,'figDir',figDir);
  end
 end
end


%%%%%%%%%%%%%
% testnaka2 %
%%%%%%%%%%%%%
function [fitfun Rmax c50 n q offset] = testnaka2(fitparams,contrasts,testType,numConditions)
% set params of the naka rushton depending on the type of test being
% performed:
[Rmax c50 n q offset] = setNakaParams2(fitparams,testType,numConditions);

% calculate function
for iAttend = 1:numConditions
 % retrieve parameters for this CRF
 thisRmax = Rmax(iAttend);
 thisc50 = c50(iAttend);
 thisn = n(iAttend);
 thisq = q(iAttend);
 thisoffset = offset(iAttend);
 thiscontrasts = squeeze(contrasts(iAttend,:));
 % now compute responses for this CRF
 fitfun(iAttend,:) = thisRmax*(thiscontrasts.^thisn)./(thiscontrasts.^(thisn*thisq)+thisc50.^(thisn*thisq))+thisoffset;
end


%%%%%%%%%%%%%%%%%%%
% setNakaParams2 %%
%%%%%%%%%%%%%%%%%%%
function [Rmax c50 n q offset] = setNakaParams2(fitparams,testType,numConditions)
% init values
Rmax   = zeros(1,numConditions);
c50    = zeros(1,numConditions);
n      = zeros(1,numConditions);
q      = zeros(1,numConditions);
offset = zeros(1,numConditions);

% this is for general debug of the params:
% disp(sprintf('[makeSimFitIndividaulPlots:setNakaParams2] Test: %s - Fitparams: [%s].',testType,num2str(fitparams)))

% testType works as follows. It is a string of four numbers for each of the four
% parametes in the naka-rushton equation (Rmax, c50, n and offset - in that order).
% (0) If you set the number to 0, it means that all curves gets the same value
% for that parameter.
% (1) If you set the number to 1, it means that every curve gets its own value
% for that parameter.
% (2) If you set the number to 2, it means that every adaptation condition gets its
% own value for that parameter but the value is the same across attention conditions
% (3) If you set the number to 3, it means that every attention condition gets its
% own value for that parameter but the value is the same across adaptation conditions
%
% note that parameters are specified as a nAdapt x nAttend x numParams array
switch testType
 case {'11111'} % fit all params
  % and grab the values
  Rmax   = fitparams(1:numConditions);
  c50    = fitparams(5:8);
  n      = fitparams(9:12);
  q      = fitparams(13:16);
  offset = fitparams(17:20);
  
 case {'11000'} % fit all params
  % Rmax
  for iAttend = 1:numConditions
   Rmax(iAttend) = fitparams(iAttend);
  end
  
  % c50
  startIndex = iAttend;
  for iAttend = 1:numConditions
   c50(iAttend) = fitparams(startIndex+iAttend);
  end
  
  % n
  nextIndex = startIndex + 1;
  n(:) = fitparams(nextIndex);
  
  % q (q is for the whole fit between 0 and 1, and it is multiplied by n only at the end of everything)
  nextIndex = nextIndex + 1;
  q(:) = fitparams(nextIndex);
  
  % offSet
  nextIndex = nextIndex+1;
  offset(:) = fitparams(nextIndex);
  
 case {'11001'}
  % Rmax
  for iAttend = 1:numConditions
   Rmax(iAttend) = fitparams(iAttend);
  end
  
  % c50
  startIndex = iAttend;
  for iAttend = 1:numConditions
   c50(iAttend) = fitparams(startIndex+iAttend);
  end
  
  % n
  nextIndex = startIndex+iAttend;
  n(:) = fitparams(nextIndex);
  
  % q (q is for the whole fit between 0 and 1, and it is multiplied by n only at the end of everything)
  nextIndex = nextIndex+1;
  q(:) = fitparams(nextIndex);
  
  % offSet
  startIndex = nextIndex;
  for iAttend = 1:numConditions
   offset(iAttend) = fitparams(startIndex+iAttend);
  end
  
 case {'01001'}
  % Rmax
  Rmax(:) = fitparams(1);
  
  % c50
  startIndex = 1;
  for iAttend = 1:numConditions
   c50(iAttend) = fitparams(startIndex+iAttend);
  end
  
  % n
  nextIndex = startIndex+1;
  n(:) = fitparams(nextIndex);
  
  % q (q is for the whole fit between 0 and 1, and it is multiplied by n only at the end of everything)
  nextIndex = nextIndex+1;
  q(:) = fitparams(nextIndex);

  % offSet
  startIndex = nextIndex;
  for iAttend = 1:numConditions
   offset(iAttend) = fitparams(startIndex+iAttend);
  end
  
 case {'10001'}
  % Rmax
  for iAttend = 1:numConditions
   Rmax(iAttend) = fitparams(iAttend);
  end
  
  % c50
  startIndex = iAttend+1;
  c50(:) = fitparams(startIndex);
  
  % n
  nextIndex = startIndex+1;
  n(:) = fitparams(nextIndex);
  
  % q (q is for the whole fit between 0 and 1, and it is multiplied by n only at the end of everything)
  nextIndex = nextIndex+1;
  q(:) = fitparams(nextIndex);

  % offSet
  startIndex = nextIndex;
  for iAttend = 1:numConditions
   offset(iAttend) = fitparams(startIndex+iAttend);
  end
  
 case {'10000'}
  % Rmax
  for iAttend = 1:numConditions
   Rmax(iAttend) = fitparams(iAttend);
  end
  
  % c50
  startIndex = iAttend+1;
  c50(:) = fitparams(startIndex);
  
  % n
  nextIndex = startIndex+1;
  n(:) = fitparams(nextIndex);
  
  % q (q is for the whole fit between 0 and 1, and it is multiplied by n only at the end of everything)
  nextIndex = nextIndex+1;
  q(:) = fitparams(nextIndex);
  
  % offSet
  startIndex = nextIndex;
  offset(:) = fitparams(startIndex);
  
 case {'01000'}
  % Rmax
  Rmax(:) = fitparams(1);
  
  % c50
  startIndex = 1;
  for iAttend = 1:numConditions
   c50(iAttend) = fitparams(startIndex+iAttend);
  end
  
  % n
  nextIndex = startIndex+iAttend+1;
  n(:) = fitparams(nextIndex);
  
  % q (q is for the whole fit between 0 and 1, and it is multiplied by n only at the end of everything)
  nextIndex = nextIndex+1;
  q(:) = fitparams(nextIndex);
  
  % offSet
  nextIndex = nextIndex+1;
  offset(:) = fitparams(nextIndex);
  
 case {'00001'}
  % Rmax
  Rmax(:) = fitparams(1);
  
  % c50
  nextIndex = 1+1;
  c50(:) = fitparams(nextIndex);
  
  % n
  nextIndex = nextIndex+1;
  n(:) = fitparams(nextIndex);
  
  % q (q is for the whole fit between 0 and 1, and it is multiplied by n only at the end of everything)
  nextIndex = nextIndex+1;
  q(:) = fitparams(nextIndex);
  
  % offSet
  startIndex = nextIndex;
  for iAttend = 1:numConditions
   offset(iAttend) = fitparams(startIndex+iAttend);
  end
  
 otherwise
  keyboard
end


%%%%%%%%%%%%%%%%%%%%%%%%%%
% switchInitialParamsNK %%
%%%%%%%%%%%%%%%%%%%%%%%%%%
function [initparams minfit maxfit] = switchInitialParamsNK2(testType,numConditions)
% Rmax = fitparams(1);
% c50 = fitparams(2);
% n = fitparams(3);
% offset = fitparams(4);

switch testType
 case {'11111'} % fit all params
  % init the params to 0
  initparams = zeros([numConditions 5]);
  minfit     = zeros([numConditions 5]);
  maxfit     = zeros([numConditions 5]);
  
  % Rmax
  initparams(:,1) = 1;
  minfit(:,1) = .5;
  maxfit(:,1) = inf;
  
  % c50
  initparams(:,2) = .1;
  minfit(:,2) = 0;
  maxfit(:,2) = 1;
  
  % n
  initparams(:,3) = .2;
  minfit(:,3) = 0;
  maxfit(:,3) = inf;
  
  % q (q is for the whole fit between 0 and 1, and it is multiplied by n only at the end of everything)
  initparams(:,4) = 0.8;
  minfit(:,4) = .1;
  maxfit(:,4) = 1;
  
  % offset
  initparams(:,5) = .05;
  minfit(:,5) = 0;
  maxfit(:,5) = inf;
  
  % make parameters into a linear array
  initparams = initparams(:);
  minfit = minfit(:);
  maxfit = maxfit(:);
  
  
 case {'11000'} % fit all params
  % Rmax
  startIndex = 1;
  endIndex = startIndex+numConditions;
  initparams(startIndex:endIndex) = 1;
  minfit(startIndex:endIndex) = .5;
  maxfit(startIndex:endIndex) = inf;
  
  % c50
  startIndex = endIndex+1;
  endIndex = startIndex+numConditions;
  initparams(startIndex:endIndex) = .1;
  minfit(startIndex:endIndex) = 0;
  maxfit(startIndex:endIndex) = inf;
  
  % n
  nextIndex = endIndex+1;
  initparams(nextIndex) = .2;
  minfit(nextIndex) = 0;
  maxfit(nextIndex) = inf;
  
  % q (q is for the whole fit between 0 and 1, and it is multiplied by n only at the end of everything)
  nextIndex = nextIndex+1;
  initparams(nextIndex) = .8;
  minfit(nextIndex) = 0;
  maxfit(nextIndex) = 1;
  
  % offSet
  nextIndex = nextIndex+1;
  initparams(nextIndex) = .1;
  minfit(nextIndex) = 0;
  maxfit(nextIndex) = inf;
  
 case {'01001'} % fit free c50 and offset but one n, q and rmax for all
  % Rmax
  nextIndex = 1;
  initparams(nextIndex) = 1;
  minfit(nextIndex) = .5;
  maxfit(nextIndex) = inf;
  
  % c50
  startIndex = nextIndex+1;
  endIndex = startIndex+numConditions-1;
  initparams(startIndex:endIndex) = .1;
  minfit(startIndex:endIndex) = 0;
  maxfit(startIndex:endIndex) = inf;
  
  % n
  nextIndex = endIndex+1;
  initparams(nextIndex) = .2;
  minfit(nextIndex) = 0;
  maxfit(nextIndex) = inf;
  
  % q (q is for the whole fit between 0 and 1, and it is multiplied by n only at the end of everything)
  nextIndex = nextIndex+1;
  initparams(nextIndex) = .8;
  minfit(nextIndex) = 0;
  maxfit(nextIndex) = 1;
  
  % offSet
  startIndex = nextIndex+1;
  endIndex = startIndex+numConditions-1;
  initparams(startIndex:endIndex) = .01;
  minfit(startIndex:endIndex) = 0;
  maxfit(startIndex:endIndex) = inf;

  
 case {'11001'} % fit free c50 and offset but one n, q and rmax for all
  % Rmax
  startIndex = 1;
  endIndex = numConditions;
  initparams(startIndex:endIndex) = 1;
  minfit(startIndex:endIndex)     = .5;
  maxfit(startIndex:endIndex)     = inf;
  
  % c50
  startIndex = endIndex+1;
  endIndex = startIndex+numConditions-1;
  initparams(startIndex:endIndex) = .1;
  minfit(startIndex:endIndex)     = 0;
  maxfit(startIndex:endIndex)     = inf;
  
  % n
  nextIndex = endIndex+1;
  initparams(nextIndex) = .2;
  minfit(nextIndex) = 0;
  maxfit(nextIndex) = inf;
  
  % q (q is for the whole fit between 0 and 1, and it is multiplied by n only at the end of everything)
  nextIndex = nextIndex+1;
  initparams(nextIndex) = .8;
  minfit(nextIndex) = 0;
  maxfit(nextIndex) = 1;
  
  % offSet
  startIndex = nextIndex+1;
  endIndex = startIndex+numConditions-1;
  initparams(startIndex:endIndex) = .01;
  minfit(startIndex:endIndex) = 0;
  maxfit(startIndex:endIndex) = inf;
  
 case {'10001'} % fit free c50 and offset but one n, q and rmax for all
  % Rmax
  startIndex = 1;
  endIndex = startIndex+numConditions-1;
  initparams(startIndex:endIndex) = 1;
  minfit(startIndex:endIndex)     = .5;
  maxfit(startIndex:endIndex)     = inf;
  
  % c50
  nextIndex = endIndex+1;
  initparams(nextIndex) = .1;
  minfit(nextIndex)     = 0;
  maxfit(nextIndex)     = inf;
  
  % n
  nextIndex = nextIndex+1;
  initparams(nextIndex) = .2;
  minfit(nextIndex) = 0;
  maxfit(nextIndex) = inf;
  
  % q (q is for the whole fit between 0 and 1, and it is multiplied by n only at the end of everything)
  nextIndex = nextIndex+1;
  initparams(nextIndex) = .8;
  minfit(nextIndex) = 0;
  maxfit(nextIndex) = 1;
  
  % offSet
  startIndex = nextIndex+1;
  endIndex = startIndex+numConditions-1;
  initparams(startIndex:endIndex) = .01;
  minfit(startIndex:endIndex) = 0;
  maxfit(startIndex:endIndex) = inf;
  
 case {'10000'} % fit free c50 and offset but one n, q and rmax for all
  % Rmax
  startIndex = 1;
  endIndex = startIndex+numConditions;
  initparams(startIndex:endIndex) = 1;
  minfit(startIndex:endIndex)     = .5;
  maxfit(startIndex:endIndex)     = inf;
  
  % c50
  nextIndex = endIndex+1;
  initparams(nextIndex) = .1;
  minfit(nextIndex)     = 0;
  maxfit(nextIndex)     = inf;
  
  % n
  nextIndex = endIndex+1;
  initparams(nextIndex) = .2;
  minfit(nextIndex) = 0;
  maxfit(nextIndex) = inf;
  
  % q (q is for the whole fit between 0 and 1, and it is multiplied by n only at the end of everything)
  nextIndex = nextIndex+1;
  initparams(nextIndex) = .8;
  minfit(nextIndex) = 0;
  maxfit(nextIndex) = 1;
  
  % offSet
  nextIndex = nextIndex+1;
  initparams(nextIndex) = .01;
  minfit(nextIndex) = 0;
  maxfit(nextIndex) = inf;
  
 case {'01000'} % fit free c50 and offset but one n, q and rmax for all
  % Rmax
  nextIndex = 1;
  initparams(nextIndex) = 1;
  minfit(nextIndex)     = .5;
  maxfit(nextIndex)     = inf;
  
  % c50
  startIndex = nextIndex + 1;
  endIndex = startIndex+numConditions-1;
  initparams(startIndex:endIndex) = .1;
  minfit(startIndex:endIndex)     = 0;
  maxfit(startIndex:endIndex)     = inf;
  
  % n
  nextIndex = endIndex+1;
  initparams(nextIndex) = .2;
  minfit(nextIndex) = 0;
  maxfit(nextIndex) = inf;
  
  % q (q is for the whole fit between 0 and 1, and it is multiplied by n only at the end of everything)
  nextIndex = nextIndex+1;
  initparams(nextIndex) = .8;
  minfit(nextIndex) = 0;
  maxfit(nextIndex) = 1;
  
  % offSet
  nextIndex = nextIndex+1;
  initparams(nextIndex) = .01;
  minfit(nextIndex) = 0;
  maxfit(nextIndex) = inf;
  
 case {'00001'} % fit free c50 and offset but one n, q and rmax for all
  % Rmax
  nextIndex = 1;
  initparams(nextIndex) = 1;
  minfit(nextIndex)     = .5;
  maxfit(nextIndex)     = inf;
  
  % c50
  nextIndex = nextIndex + 1;
  initparams(nextIndex) = .1;
  minfit(nextIndex)     = 0;
  maxfit(nextIndex)     = inf;
  
  % n
  nextIndex = nextIndex + 1;
  initparams(nextIndex) = .2;
  minfit(nextIndex) = 0;
  maxfit(nextIndex) = inf;
  
  % q (q is for the whole fit between 0 and 1, and it is multiplied by n only at the end of everything)
  nextIndex = nextIndex+1;
  initparams(nextIndex) = .8;
  minfit(nextIndex) = 0;
  maxfit(nextIndex) = 1;
  
  % offSet
  startIndex = nextIndex + 1;
  endIndex = startIndex + numConditions-1;
  initparams(startIndex:endIndex) = .01;
  minfit(startIndex:endIndex) = 0;
  maxfit(startIndex:endIndex) = inf;
  
 otherwise
  keyboard
end

% check that we are setting the right number of parameters:
freeParams = testType =='1';
numFreeParams = sum(freeParams * numConditions);
numSingleParams = sum(~freeParams);
wantedParams = numFreeParams + numSingleParams;

if ~(length(initparams) == wantedParams)
 disp(sprintf('[makeFitSimIndividualPlotsNested_R1] NOT the right number of parameters - TestType: %s, Number of wanted parameters: %s, Number of parameters: %s', ...
  testType,num2str(length(initparams)),num2str(wantedParams)));
 keyboard
end


%%%%%%%%%%%%%%%%%%%%
% setDataStructure %
%%%%%%%%%%%%%%%%%%%%
function r = setDataStructure(data,meanC, visualArea,analysisType)

% % chose adaptation, contrast and attention condtion:
% % check how many conditions the data set has (4=main analysis, 6=correct/incorrect or 1st/2nd)
% if size(fieldnames(data),1)/2 == 6
%  analysisType = '6conditions';
% elseif size(fieldnames(data),1)/2 == 4
%  analysisType = 'standard4conditions';
% else
%  keyboard
% end
 
[thisCRF  thisSTE thisContrast plotInfo] = chooseCondition(visualArea,analysisType);

% make arrays for pedestal contrasts and crfs:
for j = 1:length(thisContrast)
 r.pedestals(j,:)   = meanC.(thisContrast{j})./100;
  
 r.use_crf(j,:)     = data.(thisCRF{j});
 r.use_crf_ste(j,:) = data.(thisSTE{j});
end

r.used_crf = 'all';
r.plotInfo = plotInfo;
r.plotInfo.figureTitle = visualArea;


%%%%%%%%%%%%%%%%%%%%%%%%
%   chooseCondition    %
%%%%%%%%%%%%%%%%%%%%%%%%
function [thisCRF  thisSTE thisContrast plotInfo] = chooseCondition(visualArea,analysisType)

switch analysisType
 case {'standard' '4conditions' '4c' '4' 's'}
  % NB these fields need to be match to the ones set in
  % runNestedHypothesis_r.m by the local function <loadData>
  thisCRF = {'A_amplitude' 'Dt_amplitude' 'Dnt_amplitude' 'U_amplitude'};
  thisContrast = {'quantileAcontrast' 'quantileDcontrast' 'pedC' 'pedC'};
  thisSTE = {'A_ste' 'Dt_ste' 'Dnt_ste' 'U_ste'};
 
 case {'CorrectIncorrect' 'corrincorr' 'CorrIncorr' 'ci' '1st2ndInterval' '1st2nd' '12'} 
  % this is the correct/incorrect or 1st/2nd interval set up
  % NB these fields need to be match to the ones set in
  % runNestedHypothesis_r.m by the local function <loadData>
  thisCRF = {'A1_amplitude' 'Dt1_amplitude' 'Dnt_amplitude' 'U_amplitude'  'A2_amplitude' 'Dt2_amplitude' };
  thisContrast = {'quantileAcontrast' 'quantileDcontrast' 'pedC' 'pedC' 'quantileAcontrast' 'quantileDcontrast'};
  thisSTE = {'A1_ste' 'Dt1_ste' 'Dnt_ste' 'U_ste' 'A2_ste' 'Dt2_ste' };
  
 otherwise
  keyboard
end

% figure formatting info:
plotInfo.MarkerEdgeColor = [0 0 0];
plotInfo.thisSymbol = {'o-' 'o-' 'o-' 'o-' 'o-' 'o-'};

% generate some nice colors:
numC =60;
for i = 1:numC
 c{i} = getSmoothColor(i,numC,'hsv');
%   plot(i,1,'ko','MarkerFaceColor',c{i},'MarkerSize',26);hold on
%   text(i,1,sprintf('%i',i),'HorizontalAlignment','Center','Color',[0 0 0]);
end


% choose the ones i like:
reds    = {c{2} c{5} c{7}};
greens  = {c{23} c{17} c{15}};
blues   = {c{40} c{37} c{34}};
purples = {c{46} c{49} c{54}};
yellows = {c{9} c{10} c{11}}; % NB these last two conditions are used for the 1st/2nd intraval and correct/incorrect analysis
cyans   = {c{31} c{30} c{28}};% they DO NOT change the color coding of the figure in general

for i = 1:3 % adaptation
 plotInfo.thisColor{i}       = {reds{i} blues{i} purples{i} greens{i} yellows{i} cyans{i}};
 plotInfo.MarkerFaceColor{i} = plotInfo.thisColor{i};
end


% plot basic set up:
plotInfo.XYColor = [0 0 0];
plotInfo.Fsize = 6;
plotInfo.MarkerSize = 10;
plotInfo.LineWidth = 1.5;
plotInfo.TickLength = [0.025 .01];
plotInfo.PlotBoxAspectRatio = [1 1 1];
plotInfo.YLim = [0 1.45];
plotInfo.Xlim = [0.00109375 1];
plotInfo.plotPeds = [plotInfo.Xlim(1) 0.0175 0.035 0.07 0.14 0.28 0.56 0.84];


% set up subplotinfo
switch visualArea
 case {'v1' 'V1'}
  plotInfo.subplotInfo{1} = [4,1,1];
 case {'v2' 'V2'}
  plotInfo.subplotInfo{1} = [4,1,1+3];
 case {'v3' 'V3'}
  plotInfo.subplotInfo{1} = [4,1,1+6];
 case {'v4' 'V4'}
  plotInfo.subplotInfo{1} = [4,1,1+9];
 otherwise
end


