function r = makeFitSim(data,meanC,varargin)
% function r = makeFitSim(data,meanC,varargin)
%
% this function is called by runStatTestAttentionSim or
% runStatTestAdaptation.
%
% it gets some data{adaptation} and contrasts and fits different models
% to the resulting contrast response functions.
%
% varargin (Default):
% {'adaptation=0', 'testType=''1111''','dispFit=1', 'baselineParams=[]', ...
%  'crfFittype=naka','numBootstraps=200','use0contrast=0','dataTitle=[]', 'doBootstrap=0', ...
%  'subplotInfo',[1 1 1]}%
%
% franco pestilli 2009.01.31

% check arguments
if nargin < 2
 help makeFitSim
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
subplotInfo = [];    % the number of rows and columns for the subplot and the index to the plot
displayAllFits = []; % this is to display the fits as the procedure advances
saveFig     = [];
saveFitFig     = []; % it saves a figure with a panel fo each crf fit, this is used mostly to debug the fit procedure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set arguments
[argName argValue] = getArgs(varargin,{ ...
 'testType','1111', ...
 'dispFit',0, ...
 'crfFittype','naka', ...
 'dataTitle',[], ...
 'subplotInfo',[1 1 1], ...
 'visualArea','v1', ...
 'saveFig',0, ...
 'saveFitFig', [], ...
 'displayAllFits',0});


if (length(testType) > 3) && (strcmp(crfFittype,'exp')) || ...
  (length(testType) < 4) && (strcmp(crfFittype,'naka'))
 disp(sprintf('[%s] number of parameters and test type do not match.',mfilename))
 keyboard
end

d = setDataStructure(data, meanC, subplotInfo, visualArea);


% do the simulataneous fit:
r = fitCRF(d,crfFittype, testType,saveFitFig);

% save call info
r.callInfo.function.mfilename  = mfilename('fullpath');
r.callInfo.function.dbstack    = dbstack;
r.callInfo.function.varargin   = varargin;
r.callInfo.function.argsName   = argName;
r.callInfo.function.argsValue  = argValue;

r.callInfo.whichCRF       = 'all';     % 'dnt' 'dt' 'a' 'u';
r.callInfo.testType       = testType;
r.callInfo.crfFittype     = crfFittype;% Possible: {'linear' {'naka' 'Naka' 'NakaRushton' 'nk'}};
r.callInfo.dispFit        = dispFit;   % 1 or 0, to display the fit (set to 2 for subsidary fits like TvC etc)
r.callInfo.dataTitle      = dataTitle; % arbitrary string to use as data title (should be subject data filename)
r.d = d;

if dispFit
 d = dispCRF(r.d,r.crf,testType,saveFig);
 r.plotInfo.figureHandle = d.plotInfo.figureHandle;
end


%%%%%%%%%%%%%%%%%
%%   dispCRF   %%
%%%%%%%%%%%%%%%%%
function r = dispCRF(r,crf,testType,saveFig)

% bring up one figure for each adaptation condition and do the plot on it:
r.plotInfo.figureHandle = ...
 smartfig(sprintf('crf_simFit_%s',testType),'reuse');

for i = 1:size(r.pedestals,1) % each adaptation condition

 % subplot layout:
 row = r.plotInfo.subplotInfo{i}(1);
 col = r.plotInfo.subplotInfo{i}(2);
 num = r.plotInfo.subplotInfo{i}(3);

 % to be clear get all the contrasts, crfs and crf_ste for the current
 % adaptation condition:
 pedestals = squeeze(r.pedestals(i,:,:));
 thisCRF = squeeze(r.use_crf(i,:,:));
 thisCRFste = squeeze(r.use_crf_ste(i,:,:));

 % now run through each crf and plot it:
 for j = 1:size(r.pedestals,2) % for each attention condition
  hold on

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
  subplot(row,col,num),
  myerrorbar(plotTheseContrasts,plotTheseCRF, ...
   'yError',plotTheseCRFste, ...
   'Symbol',r.plotInfo.thisSymbol{j}(1),...
   'MarkerFaceColor',r.plotInfo.MarkerFaceColor{i}{j}, ...
   'MarkerEdgeColor',r.plotInfo.MarkerEdgeColor,...
   'MarkerSize',r.plotInfo.MarkerSize, ...
   'Color',r.plotInfo.thisColor{i}{j});
  
   % display the smooth fit function:
  thisSmoothX = squeeze(crf.fit.fitx(i,j,:));
  thisSmoothY = squeeze(crf.fit.fity(i,j,:));

  subplot(row,col,num),
  semilogy(thisSmoothX,thisSmoothY,...
   r.plotInfo.thisSymbol{j}(2:end), ...
   'Color',r.plotInfo.thisColor{i}{j});

 end

 %% now plot the fit params in the title:
 %  this_c50    = num2str(round(100*crf.fit.eachparam.c50(i,:)));
 %  this_n      = num2str(crf.fit.eachparam.n(i,:));
 %  this_offset = num2str(crf.fit.eachparam.offset(i,:));
 %  thisTile    = sprintf('%s\nc50: [%s]\nn: [%s]\noffset: [%s]\n - ', ...
 %              r.plotInfo.figureTitle,this_c50,this_n, this_offset);
 %  title(thisTile,'FontSize',1.5*r.plotInfo.Fsize)

 % figure formatting:
 if sum([row col num]) == 1
  xlabel('Contrast');
  ylabel('fMRI response (% signal change)');
 end
 axis('square')

 % axis formatting:
 set(gca,...
  'FontName','Helvetica','FontSize',r.plotInfo.Fsize, ...
  'PlotBoxAspectRatio',r.plotInfo.PlotBoxAspectRatio, ...
  'XLim', r.plotInfo.Xlim,'YLim', r.plotInfo.YLim,...
  'LineWidth',r.plotInfo.LineWidth,'TickLength',r.plotInfo.TickLength, ...
  'XScale','log', ...
  'XColor',r.plotInfo.XYColor,'YColor',r.plotInfo.XYColor,'Box','off');
%   'XTick', [r.plotInfo.plotPeds 1], ...
%   'XTickLabel', 100*[0 r.plotInfo.plotPeds(2:end) 1] ,...
 
  set(gca,'XScale','log','Xlim',[.001 1],'Ylim',r.plotInfo.YLim);

  xticks.major = [r.plotInfo.plotPeds 1];
  yticks.major = [0 .5 1];
  yticks.minor = [0 .25 .5 .75 1 1.25];
  myaxisfp(xticks,yticks);

 
 set(r.plotInfo.figureHandle,'PaperPosition',[.25 .25 8 10.5]);
 set(r.plotInfo.figureHandle,'PaperOrientation','Portrait');
end
drawnow



%%%%%%%%%%%%%%%%
%%   fitCRF   %%
%%%%%%%%%%%%%%%%
function r = fitCRF(d,crfFittype, testType,saveFitFig)
disp(sprintf('[makeFitSim:fitCRF] Using function: <%s>,  with test: <%s> start...',crfFittype,testType));

switch crfFittype
 case {'naka' 'NAKA' 'nk'}
  % fit a naka-rushton sigmodal function:
  r.crf.fit = fitsigmoidtest(d.pedestals,d.use_crf, testType,saveFitFig);
 
 case {'naka2' 'NAKA2' 'nk2'}
  % fit a naka-rushton sigmodal function:
  r.crf.fit = fitsigmoidtest2(d.pedestals,d.use_crf, testType,saveFitFig);
  
 case {'exp' 'EXP' 'e'}
  %% this is the fuction i need to fit an exponential
  r.crf.fit = fitexptest(d.pedestals,d.use_crf, testType,saveFitFig);

 otherwise
  disp(sprintf('[makeFitSim:fitCRF] %s fit type not defined',crfFittype))
  keyboard
end


disp(sprintf('[makeFitSimIndividualPlots:fitCRF] Using function: <%s>,  with test: <%s> done...',crfFittype,testType));

%% new stuff starts here:
%%%%%%%%%%%%%%%%%%%%%%%%%%%
function bestfit = fitsigmoidtest2(contrasts,responses, testType,saveFitFig)

% check arguments
if (nargin ~= 4)
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
numConditions = size(responses);    % first is adaptation second is attention
numConditions = numConditions(1:2); % the fit is run for adaptation and attention

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
 contrasts,responses,testType,numConditions,saveFitFig);

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

saveFitFig.save = originalSaveFitFig;
saveFitFig.displayAllFits=1;%fix fix change this later now only displays the last fit
[bestfit.err bestfit.yfit] = sigmoiderr2(fitparams,contrasts,responses,testType,numConditions,saveFitFig);

% compute R2 values and save them in the fit:
% (1) reshape the errors by condition and contrast
bestfit.err = reshape(bestfit.err,size(responses));

% (2) run through conditions and compute r2
% for each one of them
for i = 1:size(bestfit.err,1)  % adaptation conditions
 for j = 1:size(bestfit.err,2) % attention conditions
  thisError = squeeze(bestfit.err(i,j,:));
  thisResponses = squeeze(responses(i,j,:));
  bestfit.r2(i,j) = 1-var(thisError)/var(thisResponses);
 end
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

% THIS IS NOT DOABLE NOW UNLESS I FIND A WAY TO ESTIMATE THE OVERAL VARIABILITY IN THE DATA:
% (4) compute  chi-squared test of goodness of fit.
% bestFit.chi2.chi2value = sum((thisError./thisResponses).^2);
% bestFit.chi2.M = length(fitparams);
% bestFit.chi2.N = numel(responses);
% bestFit.chi2.nparams = bestFit.chi2.N - bestFit.chi2.M;
% bestFit.chi2.dataVariability = **ESTIMATE ME**; % maybe the sum of the variances for each datapoint divided the number of datapoints? 
% bestFit.chi2.modelFitError = thisError;
% bestFit.chi2.p = 1 - gammainc(0.5*bestFit.chi2.chi2value,0.5*bestFit.chi2.nparams); 

% (5) compute the function returned by the best fit:
x = .0001:.0001:1;
bestfit.fitx = zeros([numConditions length(x)]);
for i = 1:numConditions(1)
 for j = 1:numConditions(2)
  bestfit.fitx(i,j,:) = x;
 end
end

% compute th current model estimate:
[bestfit.fity Rmax c50 n q offset] = testnaka2(bestfit.params,bestfit.fitx,testType,numConditions);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [err fitfun] = sigmoiderr2(fitparams,contrasts,responses,testType,numConditions,saveFitFig)

% compute th current model estimate:
[fitfun Rmax c50 n q offset] = testnaka2(fitparams,contrasts,testType,numConditions);

% calculate error:
err = responses-fitfun;
err = err(:);

% display the current fit (FIX: This doesn't yet work for displaying simultaneous fit)
global numIters displayAllFits

numIters = numIters+1;
if saveFitFig.displayAllFits
 figName = sprintf('show_all_fit_mfs_OBS%s_TEST%s_VISAR%s',saveFitFig.observer,saveFitFig.testType,saveFitFig.visualarea);
 f = smartfig(figName,'reuse');
 colors = {'r' 'b' 'm' 'g'};
 
 % display the fit as we go along
 if ~isnan(numIters)
  h = get(gcf,'Child'); % handles to the axes
  if numel(h) > 1
   for a = prod([size(responses,1) size(responses,2)]):-1:1
    cla(h(a))
   end
  end
  
  c=1;
  for i = 1:size(responses,1)  % adaptation
   for j = 1:size(responses,2) % attention
    % extract the current values:
    thisC = 100*squeeze(contrasts(i,j,:));
    thisR = squeeze(responses(i,j,:));
    thisFitFun = squeeze(fitfun(i,j,:));
    this_c50 = 100*c50(i,j);
    this_Rmax = Rmax(i,j);
    this_n = n(i,j);
    this_q = q(i,j);
    this_offset = offset(i,j);
    [sortcontrasts sortindex] = sort(thisC);
    
    subplot(size(responses,1),size(responses,2),c),
    semilogx(thisC,thisR,'ko');hold on
    semilogx(sortcontrasts,thisFitFun(sortindex),[colors{j},'-']);
    vline(this_c50);
    title(sprintf('Rmax=%0.2f,\n c50=%0.2f,\n n=%0.2f,\n q=%0.2f,\n offset=%0.2f\n numIters=%i', ...
     this_Rmax,this_c50,this_n,this_q,this_offset,numIters));
    axis([1 100 0 1.3])
    drawnow
    c=c+1;
   end
  end
  
  % clear the figure:
  if saveFitFig.save
   set(f,'PaperPosition',[.25 .25 8 10.5]);
   set(f,'PaperOrientation','Portrait');
   savefig(f,'figName',figName,'figDir','testnaka2_bestfit_figs_freeC50');
  end
 end
end

%%%%%%%%$%%%%%%%%%%%%%%%%%%%
function [fitfun Rmax c50 n q offset] = testnaka2(fitparams,contrasts,testType,numConditions)

% set params of the naka rushton depending on the type of test being
% performed:
[Rmax c50 n q offset] = setNakaParams2(fitparams,testType,numConditions);

% calculate function
for iAdapt = 1:numConditions(1)
 for iAttend = 1:numConditions(2)
  % retrieve parameters for this CRF
  thisRmax = Rmax(iAdapt,iAttend);
  thisc50 = c50(iAdapt,iAttend);
  thisn = n(iAdapt,iAttend);
  thisq = q(iAdapt,iAttend);
  thisoffset = offset(iAdapt,iAttend);
  thiscontrasts = squeeze(contrasts(iAdapt,iAttend,:));
  % now compute responses for this CRF
  fitfun(iAdapt,iAttend,:) = thisRmax*(thiscontrasts.^thisn)./(thiscontrasts.^(thisn*thisq)+thisc50.^(thisn*thisq))+thisoffset;
 end
end


%%%%%%%%%%%%%%%%%%%%
%% setNakaParams2 %%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Rmax c50 n q offset] = setNakaParams2(fitparams,testType,numConditions)
% init values
Rmax   = zeros(numConditions);
c50    = zeros(numConditions);
n      = zeros(numConditions);
q      = zeros(numConditions);
offset = zeros(numConditions);

% this is for general debug of the params:
% disp(sprintf('[makeSim:setNakaParams2] Test: %s - Fitparams: [%s].',testType,num2str(fitparams)))

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
  % make fitparams back into a nAttention x nAdaptation x nParams matrix;
  fitparams = reshape(fitparams,[numConditions 5]);
  % and grab the values
  Rmax   = fitparams(:,:,1);
  c50    = fitparams(:,:,2);
  n      = fitparams(:,:,3);
  q      = fitparams(:,:,4);
  offset = fitparams(:,:,5);
  
 case {'11000'} % fit all params
  % Rmax
  startIndex = 1;
  c = 0;
  for iAdapt = 1:numConditions(1)
   for iAttend = 1:numConditions(2)
    c = c+1;
    Rmax(iAdapt,iAttend) = fitparams(startIndex+c);
   end
  end
  
  % c50
  % fix fix change
  startIndex = c;
  c = 0;
  for iAdapt = 1:numConditions(1)
   for iAttend = 1:numConditions(2)
    c = c+1;
    c50(iAdapt,iAttend) = fitparams(startIndex+c);
   end
  end
  
  % n
  nextIndex = startIndex + c;
  n(:) = fitparams(nextIndex);
  
  
  % q (q is for the whole fit between 0 and 1, and it is multiplied by n only at the end of everything)
  nextIndex = nextIndex+1;
  q(:) = fitparams(nextIndex);
  
  % offSet
  nextIndex = nextIndex+1;
  offset(:) = fitparams(nextIndex);
  
 case {'01001'} % fit free c50 and offset but one n, q and rmax for all
  % Rmax
  Rmax(:) = fitparams(1);
  
  % c50
  startIndex = 1;
  c = 0;
  for iAdapt = 1:numConditions(1)
   for iAttend = 1:numConditions(2)
    c = c + 1;
    c50(iAdapt,iAttend) = fitparams(c+startIndex);
   end
  end
  
  % n
  nextIndex = c+startIndex+1;
  n(:) = fitparams(nextIndex);
  
  % q (q is for the whole fit between 0 and 1, and it is multiplied by n only at the end of everything)
  nextIndex = nextIndex+1;
  q(:) = fitparams(nextIndex);
  
  % offSet
  startIndex = nextIndex;
  c = 0;
  for iAdapt = 1:numConditions(1)
   for iAttend = 1:numConditions(2)
    c = c + 1;
    offset(iAdapt,iAttend) = fitparams(c+startIndex);
   end
  end
  
 case {'02003'} % allow c50 to change across adaptation and offset across attention conditions
  % Rmax, n and q are common across conditions:
  % now we interpret the fitparams as the following
  % e.g. [Rmax c50Adapt1 c50Adapt2 c50Adapt3 n q offsetAttend1 offsetAttend2 offsetAttend3 offsetAttend4]
  
  % Rmax
  Rmax(:) = fitparams(1);
  
  % c50
  startIndex = 1;
  for iAdapt = 1:numConditions(1)
   c50(iAdapt,:) = fitparams(startIndex+iAdapt);
  end
  
  % n
  nextIndex = startIndex + 1;
  n(:) = fitparams(nextIndex);
  
  % q (q is for the whole fit between 0 and 1, and it is multiplied by n only at the end of everything)
  nextIndex = nextIndex + 1;
  q(:) = fitparams(nextIndex);
  
  % offset
  startIndex = nextIndex;
  for iAttend = 1:numConditions(2)
   offset(:,iAttend) = fitparams(startIndex+iAttend);
  end
  
  
 case {'03002'} % allow offset to change across adaptation conditions and c50 across attention
  % Rmax, n and q are common across conditions:
  % now we interpret the fitparams as the following
  % e.g. [Rmax c50att1 c50att2 c50att3 c50att4 n offsetAdp1 offsetAdp2 offsetAdp3]
  
  % Rmax
  nextIndex = 1;
  Rmax(:) = fitparams(nextIndex);
  
  % c50
  startIndex = nextIndex;
  for iAttend = 1:numConditions(2)
   c50(:,iAttend) = fitparams(startIndex+iAttend);
  end
  
  % n
  nextIndex = startIndex + iAttend;
  n(:) = fitparams(nextIndex);
  
  % q
  nextIndex = nextIndex + iAttend;
  q(:) = fitparams(nextIndex);
  
  % offset
  startIndex = nextIndex;
  for iAdapt = 1:numConditions(1)
   offset(iAdapt,:) = fitparams(startIndex+iAdapt);
  end
  
  
 case {'23000'} % allow Rmax to change across attention and c50 across adaptation
  % now we interpret the fitparams as the following
  % e.g. [Rmax Rmax Rmax Rmax c50 c50 c50 n q offset]
  
  % Rmax
  startIndex = 0;
  for iAttend = 1:numConditions(2)
   Rmax(:,iAttend)  = fitparams(startIndex + iAttend);
  end
  
  % c50
  startIndex = startIndex+iAttend;
  for iAdapt = 1:numConditions(1)
   c50(iAdapt,:)   = fitparams(startIndex+iAdapt);
  end
  
  % n
  nextIndex = startIndex+iAdapt+1;
  n(:)      = fitparams(nextIndex);
  
  % q
  nextIndex = nextIndex+1;
  q(:)      = fitparams(nextIndex);
  
  % offset
  nextIndex = nextIndex+1;
  offset(:) = fitparams(nextIndex);
  
  
 case {'23111'} % allow Rmax to change across attention and c50 across adaptation
  % now we interpret the fitparams as the following
  % e.g. [Rmax Rmax Rmax Rmax c50 c50 c50 n q offset]
  
  % Rmax
  for iAttend = 1:numConditions(2)
   Rmax(:,iAttend)  = fitparams(iAttend);
  end
  
  % c50
  startIndex = iAttend;
  for iAdapt = 1:numConditions(1)
   c50(iAdapt,:)   = fitparams(startIndex+iAdapt);
  end
  
  % n
  startIndex = startIndex+iAdapt;
  c = 0;
  for iAdapt = 1:numConditions(1)
   for iAttend = 1:numConditions(2)
    c = c+1;
    n(iAdapt,iAttend) = fitparams(startIndex+c);
   end
  end
  
  % q
  startIndex = c+startIndex;
  c = 0;
  for iAdapt = 1:numConditions(1)
   for iAttend = 1:numConditions(2)
    c = c+1;
    q(iAdapt,iAttend) = fitparams(startIndex+c);
   end
  end
  
  % offset
  startIndex = c+startIndex;
  c = 0;
  for iAdapt = 1:numConditions(1)
   for iAttend = 1:numConditions(2)
    c = c+1;
    offset(iAdapt,iAttend) = fitparams(startIndex+c);
   end
  end
  
  
 case {'12113'} % allow c50 to change across adaptation and offset across attention conditions
  % Rmax, n and q are free to vary across conditions:
  % now we interpret the fitparams as the following
  % e.g. [Rmax Rmax Rmax Rmax Rmax Rmax Rmax Rmax Rmax Rmax Rmax Rmax
  %       c50Adapt1 c50Adapt2 c50Adapt3
  %       n n n n n n n n n n n n
  %       q q q q q q q q q q q q
  %       offsetAttend1 offsetAttend2 offsetAttend3 offsetAttend4]
  
  % Rmax
  startIndex = 0;
  c = 0;
  for iAdapt = 1:numConditions(1)
   for iAttend = 1:numConditions(2)
    c = c+1;
    Rmax(iAdapt,iAttend) = fitparams(startIndex+c);
   end
  end
  
  % c50
  startIndex = c + startIndex;
  for iAdapt = 1:numConditions(1)
   c50(iAdapt,:) = fitparams(startIndex+iAdapt);
  end
  
  % n
  startIndex = startIndex+iAdapt;
  c = 0;
  for iAdapt = 1:numConditions(1)
   for iAttend = 1:numConditions(2)
    c = c+1;
    n(iAdapt,iAttend) = fitparams(startIndex+c);
   end
  end
  
  % q
  startIndex = c+startIndex;
  c = 0;
  for iAdapt = 1:numConditions(1)
   for iAttend = 1:numConditions(2)
    c = c+1;
    q(iAdapt,iAttend) = fitparams(startIndex+c);
   end
  end
  
  % offset
  startIndex = startIndex + c;
  for iAttend = 1:numConditions(2)
   offset(:,iAttend) = fitparams(startIndex+iAttend);
  end
  
  
 case {'13112'} % allow offset to change across adaptation conditions and c50 across attention
  % Rmax, n and q are free to vary across conditions:
  % now we interpret the fitparams as the following
  % e.g. [Rmax Rmax Rmax Rmax Rmax Rmax Rmax Rmax Rmax Rmax Rmax Rmax
  %       c50Attend1 c50Attend2 c50Attend3 c50Attend4
  %       n n n n n n n n n n n n
  %       q q q q q q q q q q q q
  %       offsetAdapt1 offsetAdapt2 offsetAdapt3]
  
  % Rmax
  startIndex = 0;
  c = 0;
  for iAdapt = 1:numConditions(1)
   for iAttend = 1:numConditions(2)
    c = c+1;
    Rmax(iAdapt,iAttend) = fitparams(startIndex+c);
   end
  end
  
  % c50
  startIndex = c + startIndex;
  for iAttend = 1:numConditions(2)
   c50(:,iAttend) = fitparams(startIndex+iAttend);
  end
  
  % n
  startIndex = startIndex+iAdapt;
  c = 0;
  for iAdapt = 1:numConditions(1)
   for iAttend = 1:numConditions(2)
    c = c+1;
    n(iAdapt,iAttend) = fitparams(startIndex+c);
   end
  end
  
  % q
  startIndex = c+startIndex;
  c = 0;
  for iAdapt = 1:numConditions(1)
   for iAttend = 1:numConditions(2)
    c = c+1;
    q(iAdapt,iAttend) = fitparams(startIndex+c);
   end
  end
  
  % offset
  startIndex = startIndex+c;
  for iAdapt = 1:numConditions(1)
   offset(iAdapt,:) = fitparams(startIndex+iAdapt);
  end
  
 otherwise
  keyboard
end

% this is for general debug of the params:
% disp(sprintf('[makeSim:setNakaParams2] Test: %s - Fitparams: [%s].',testType,num2str(fitparams)))



%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% switchInitialParamsNK %%

%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
  initparams(:,:,1) = 1;
  minfit(:,:,1) = .5;
  maxfit(:,:,1) = inf;
  
  % c50
  initparams(:,:,2) = .1;
  minfit(:,:,2) = 0;
  maxfit(:,:,2) = inf;
  
  % n
  initparams(:,:,3) = .2;
  minfit(:,:,3) = 0;
  maxfit(:,:,3) = inf;
  
  % q (q is for the whole fit between 0 and 1, and it is multiplied by n only at the end of everything)
  initparams(:,:,4) = 0.8;
  minfit(:,:,4) = .1;
  maxfit(:,:,4) = 1;
  
  % offset
  initparams(:,:,5) = .15;
  minfit(:,:,5) = 0;
  maxfit(:,:,5) = inf;
  
  % make parameters into a linear array
  initparams = initparams(:);
  minfit = minfit(:);
  maxfit = maxfit(:);
 
  
 case {'02003'} % allow c50 to change across adaptation and offset across attention conditions
  % Rmax, n and q are common across conditions:
  % now we interpret the fitparams as the following
  % e.g. [Rmax c50Adapt1 c50Adapt2 c50Adapt3 n q
  %       offsetAttend1 offsetAttend2 offsetAttend3 offsetAttend4]
  
  % Rmax
  initparams(1) = 1;
  minfit(1) = .5;
  maxfit(1) = inf;
  
  % c50
  startIndex = 2;
  endIndex = startIndex+numConditions(1);
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
  
  % offsSet
  startIndex = nextIndex+1;
  endIndex = startIndex+numConditions(2);
  initparams(startIndex:endIndex) = .1;
  minfit(startIndex:endIndex) = 0;
  maxfit(startIndex:endIndex) = inf;
  
  
 case {'03002'} % allow offset to change across adaptation conditions and c50 across attention
  % Rmax, n and q are common across conditions:
  % now we interpret the fitparams as the following
  % e.g. [Rmax c50att1 c50att2 c50att3 c50att4 n q offsetAdp1 offsetAdp2 offsetAdp3]
  
  % Rmax
  initparams(1) = 1;
  minfit(1)     = .5;
  maxfit(1)     = inf;
  
  % c50
  startIndex = 2;
  endIndex = startIndex+numConditions(2);
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
  
  % offset
  startIndex = nextIndex+1;
  endIndex = startIndex+numConditions(1);
  initparams(startIndex:endIndex) = .1;
  minfit(startIndex:endIndex) = 0;
  maxfit(startIndex:endIndex) = inf;
  
 case {'11000'} % fit all params
  % Rmax
  startIndex = 1;
  endIndex = startIndex+prod(numConditions);
  initparams(startIndex:endIndex) = 1;
  minfit(startIndex:endIndex) = .5;
  maxfit(startIndex:endIndex) = inf;
  
  % c50
  startIndex = endIndex+1;
  endIndex = startIndex+prod(numConditions);
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
  endIndex = startIndex+prod(numConditions);
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
  endIndex = startIndex+prod(numConditions);
  initparams(startIndex:endIndex) = .1;
  minfit(startIndex:endIndex) = 0;
  maxfit(startIndex:endIndex) = inf;
  
  
 case {'12113'} % allow c50 to change across adaptation and offset across attention conditions
  % Rmax, n and q are FREE to very across conditions:
  % now we interpret the fitparams as the following
  % e.g. [Rmax Rmax Rmax Rmax Rmax Rmax Rmax Rmax Rmax Rmax Rmax Rmax
  %       c50Attend1 c50Attend2 c50Attend3 c50Attend4
  %       n n n n n n n n n n n n
  %       q q q q q q q q q q q q
  %       offsetAdapt1 offsetAdapt2 offsetAdapt3]
  
  % Rmax
  startIndex = 1;
  endIndex = startIndex+prod(numConditions);
  initparams(startIndex:endIndex) = 1;
  minfit(startIndex:endIndex) = .5;
  maxfit(startIndex:endIndex) = inf;
  
  % c50
  startIndex = endIndex+1;
  endIndex = startIndex +numConditions(2);
  initparams(startIndex:endIndex) = .1;
  minfit(startIndex:endIndex)     = 0;
  maxfit(startIndex:endIndex)     = inf;
  
  % n
  startIndex = endIndex+1;
  endIndex = startIndex+prod(numConditions);
  initparams(startIndex:endIndex) = .2;
  minfit(startIndex:endIndex)     = 0;
  maxfit(startIndex:endIndex)     = inf;
  
  % q (q is for the whole fit between 0 and 1, and it is multiplied by n only at the end of everything)
  startIndex = endIndex+1;
  endIndex = startIndex+prod(numConditions);
  initparams(startIndex:endIndex) = .8;
  minfit(startIndex:endIndex)     = 0;
  maxfit(startIndex:endIndex)     = 1;
  
  % offsSet
  startIndex = endIndex+1;
  endIndex   = startIndex +numConditions(1);
  initparams(startIndex:endIndex) = .1;
  minfit(startIndex:endIndex)     = 0;
  maxfit(startIndex:endIndex)     = inf;
  
  
 case {'13112'} % allow offset to change across adaptation conditions and c50 across attention
  % Rmax, n and q are FREE to vary across conditions:
  % now we interpret the fitparams as the following
  % e.g. [Rmax Rmax Rmax Rmax Rmax Rmax Rmax Rmax Rmax Rmax Rmax Rmax
  %       c50Attend1 c50Attend2 c50Attend3 c50Attend4
  %       n n n n n n n n n n n n
  %       q q q q q q q q q q q q
  %       offsetAdapt1 offsetAdapt2 offsetAdapt3]
  
  % Rmax
  startIndex = 1;
  endIndex = startIndex+prod(numConditions);
  initparams(startIndex:endIndex) = 1;
  minfit(startIndex:endIndex) = .5;
  maxfit(startIndex:endIndex) = inf;
  
  % c50
  startIndex = endIndex+1;
  endIndex = startIndex +numConditions(1);
  initparams(startIndex:endIndex) = .1;
  minfit(startIndex:endIndex)     = 0;
  maxfit(startIndex:endIndex)     = inf;
  
  % n
  startIndex = endIndex+1;
  endIndex = startIndex+prod(numConditions);
  initparams(startIndex:endIndex) = .2;
  minfit(startIndex:endIndex)     = 0;
  maxfit(startIndex:endIndex)     = inf;
  
  % q (q is for the whole fit between 0 and 1, and it is multiplied by n only at the end of everything)
  startIndex = endIndex+1;
  endIndex = startIndex+prod(numConditions);
  initparams(startIndex:endIndex) = .8;
  minfit(startIndex:endIndex)     = 0;
  maxfit(startIndex:endIndex)     = 1;
  
  % offsSet
  startIndex = endIndex+1;
  endIndex   = startIndex +numConditions(2);
  initparams(startIndex:endIndex) = .1;
  minfit(startIndex:endIndex)     = 0;
  maxfit(startIndex:endIndex)     = inf;
  
  
 case {'23000'} % allow Rmax to change across attention and c50 across adaptation
  % now we interpret the fitparams as the following
  % e.g. [Rmax Rmax Rmax Rmax c50 c50 c50 n q offset]
  
  % Rmax
  startIndex = 1;
  endIndex   = numConditions(2);
  initparams(startIndex:endIndex)  = 1;
  minfit(startIndex:endIndex) = .5;
  maxfit(startIndex:endIndex) = inf;
  
  % c50
  startIndex = endIndex+1;
  endIndex   = startIndex+numConditions(1);
  initparams(startIndex:endIndex)  = .1;
  minfit(startIndex:endIndex) = 0;
  maxfit(startIndex:endIndex) = inf;
  
  % n
  nextIndex = endIndex+1;
  initparams(nextIndex)  = 1;
  minfit(nextIndex) = 0;
  maxfit(nextIndex) = inf;
  
  % q
  nextIndex = nextIndex+1;
  initparams(nextIndex)  = .8;
  minfit(nextIndex) = .1;
  maxfit(nextIndex) = 1;
  
  % offset
  nextIndex = nextIndex+1;
  initparams(nextIndex)  = .15;
  minfit(nextIndex) = 0;
  maxfit(nextIndex) = inf;
  
  
 case {'23111'} % allow Rmax to change across attention and c50 across adaptation
  % now we interpret the fitparams as the following
  % e.g. [Rmax Rmax Rmax Rmax c50 c50 c50 n q offset]
  
  % Rmax
  startIndex = 1;
  endIndex   = startIndex+numConditions(2);
  initparams(startIndex:endIndex)  = 1;
  minfit(startIndex:endIndex) = .5;
  maxfit(startIndex:endIndex) = inf;
  
  % c50
  startIndex = endIndex+1;
  endIndex   = startIndex+numConditions(1);
  initparams(startIndex:endIndex)  = .1;
  minfit(startIndex:endIndex) = 0;
  maxfit(startIndex:endIndex) = inf;
  
  % n
  startIndex = endIndex+1;
  endIndex = endIndex + prod(numConditions);
  initparams(startIndex:endIndex)  = 1;
  minfit(startIndex:endIndex) = 0;
  maxfit(startIndex:endIndex) = inf;
  
  % q
  startIndex = endIndex+1;
  endIndex = endIndex + prod(numConditions);
  initparams(startIndex:endIndex)  = .8;
  minfit(startIndex:endIndex) = .1;
  maxfit(startIndex:endIndex) = 1;
  
  % offset
  startIndex = endIndex+1;
  endIndex = endIndex + prod(numConditions);
  initparams(startIndex:endIndex)  = .15;
  minfit(startIndex:endIndex) = 0;
  maxfit(startIndex:endIndex) = inf;
  
 otherwise
  keyboard
end


% this is for general debug of the params:
disp(sprintf('[makeSimFitIndividaulPlots] Test: %s - InitParamsLen: %s (minfitLen %s, maxfitLen %s) initparams: [%s].', ...
 testType,num2str(length(initparams)),num2str(length(minfit)),num2str(length(maxfit)),num2str(initparams)))



%% NEW STUFF ENDS HERE:


%%%%%%%%%%%%%%%%%%%%%
%%  fitsigmoidtest %%

%%%%%%%%%%%%%%%%%%%%%
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
optimizationParams = optimset( ...
 'LevenbergMarquardt','on', ...
 'MaxIter',inf, ...
 'TolFun',10^-10);

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
for i = 1:size(bestfit.err,1)  % adaptation conditions
 for j = 1:size(bestfit.err,2) % attention conditions
  thisError = squeeze(bestfit.err(i,j,:));
  thisResponses = squeeze(responses(i,j,:));
  bestfit.r2(i,j) = 1-var(thisError)/var(thisResponses);
 end
end

% (3) compute the overal r2 of the whole model
thisError = bestfit.err(:);
thisResponses = responses(:);
bestfit.model_r2 = 1-var(thisError)/var(thisResponses);

% if we have the best fit then keep it.
bestfit.params = reshape(fitparams,[numConditions 4]);
bestfit.eachparam.Rmax = Rmax;
bestfit.eachparam.c50 = c50;
bestfit.eachparam.n = n;
bestfit.eachparam.offset = offset;
bestfit.covar = covar;
bestfit.output = output;
bestfit.testType = testType;
bestfit.numConditions = numConditions;

% compute the function returned by the best fit:
x = .001:.00001:1;
bestfit.fitx = zeros([numConditions length(x)]);
for i = 1:numConditions(1)
 for j = 1:numConditions(2)
  bestfit.fitx(i,j,:) = x;
 end
end

% compute th current model estimate:
[bestfit.fity Rmax c50 n offset] = testnaka(bestfit.params,bestfit.fitx,testType,numConditions);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% function for fitting naka-rushton to crf %%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [err rss fitfun] = sigmoiderr(fitparams,contrasts,responses,testType,numConditions)

% compute th current model estimate:
[fitfun Rmax c50 n offset] = testnaka(fitparams,contrasts,testType,numConditions);

% calculate error:
err = responses-fitfun;
err = err(:);

% calculate residual sum of squares

% display the current fit (FIX: This doesn't yet work for displaying simultaneous fit)
global numIters displayAllFits
numIters = numIters+1;
if displayAllFits
 % display the fit as we go along
 if ~isnan(numIters)
  c=1;
  for i = 1:size(responses,1)  % adaptation
   for j = 1:size(responses,2) % attention
    
    % extract the current values:
    thisC = 100*squeeze(contrasts(i,j,:));
    thisR = squeeze(responses(i,j,:));
    thisFitFun = squeeze(fitfun(i,j,:));
    this_c50 = 100*c50(i,j);
    this_Rmax = Rmax(i,j);
    this_n = n(i,j);
    this_offset = offset(i,j);
    [sortcontrasts sortindex] = sort(thisC);
    
    subplot(3,4,c),
    semilogx(thisC,thisR,'ko');hold on
    semilogx(sortcontrasts,thisFitFun(sortindex),'r-');
    vline(this_c50);
    title(sprintf('Rmax=%0.2f,\n c50=%0.2f,\n n=%0.2f,\n offset=%0.2f\n numIters=%i', ...
     this_Rmax,this_c50,this_n,this_offset,numIters));
    axis([1 100 0 1.5])
    drawnow
    c=c+1;
   end
  end
  % clear the figure:
  h=get(gcf,'Child'); % handles to the axes
  for a = 12:-1:1
   cla(h(a))
  end
 end
end



%%%%%%%%%%%%%%%%%
%%  fitexptest %%

%%%%%%%%%%%%%%%%%
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
numConditions = size(responses);    % first is adaptation second is attention
numConditions = numConditions(1:2); % the fit is run for adaptation and attention

[initparams minfit maxfit] = switchInitialParamsEX(testType,numConditions);

% set optimization parameters
optimizationParams = optimset( ...
 'LevenbergMarquardt','on', ...
 'MaxIter',inf, ...
 'TolFun',10^-10);

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
for i = 1:size(bestfit.err,1)  % adaptation conditions
 for j = 1:size(bestfit.err,2) % attention conditions
  thisError = squeeze(bestfit.err(i,j,:));
  thisResponses = squeeze(responses(i,j,:));
  bestfit.r2(i,j) = 1-var(thisError)/var(thisResponses);
 end
end

% (3) compute the overal r2 of the whole model
thisError = bestfit.err(:);
thisResponses = responses(:);
bestfit.model_r2 = 1-var(thisError)/var(thisResponses);

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
x = .001:.0001:1;
bestfit.fitx = zeros([numConditions length(x)]);
for i = 1:numConditions(1)
 for j = 1:numConditions(2)
  bestfit.fitx(i,j,:) = x;
 end
end

% [dummy bestfit.fity] = sigmoiderr(bestfit.params,bestfit.fitx,bestfit.fitx,testType);
% compute th current model estimate:
[bestfit.fity c50 n offset] = testexp(bestfit.params,bestfit.fitx,testType,numConditions);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% function for fitting exponential to crf %%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [err, fitfun] = experr(fitparams,contrasts,responses,testType,numConditions)

% compute th current model estimate:
[fitfun c50 n offset] = testexp(fitparams,contrasts,testType,numConditions);

% calculate error:
err = responses-fitfun;
err = err(:);

% display the current fit (FIX: This doesn't yet work for displaying simultaneous fit)
global numIters displayAllFits
numIters = numIters+1;
if displayAllFits
 % display the fit as we go along
 if ~isnan(numIters)
  c=1;
  % clear the figure:
  f = smartfig('SimFitIterations','reuse');
  h=get(f,'Child'); % handles to the axes
  for a = 1:length(h)
   cla(h(a))
  end
  for i = 1:size(responses,1)  % adaptation
   for j = 1:size(responses,2) % attention

    % extract the current values:
    thisC = 100*squeeze(contrasts(i,j,:));
    thisR = squeeze(responses(i,j,:));
    thisFitFun = squeeze(fitfun(i,j,:));
    this_c50 = 100*c50(i,j);
    this_n = n(i,j);
    this_offset = offset(i,j);
    [sortcontrasts sortindex] = sort(thisC);

    subplot(3,4,c),
    semilogx(thisC,thisR,'ko');hold on
    semilogx(sortcontrasts,thisFitFun(sortindex),'r-');
    vline(this_c50);
    title(sprintf('c50=%0.2f,\n n=%0.2f,\n offset=%0.2f\n numIters=%i', ...
     this_c50,this_n,this_offset,numIters));
    axis([1 100 0 1.5])
    drawnow
    c=c+1;
   end
  end
 end
end


%%%%%%%%%%%%%
%% testexp %%

%%%%%%%%%%%%%
function [fitfun c50 n offset] = testexp(fitparams,contrasts,testType,numConditions)

% set params of the naka rushton depending on the type of test being
% performed:
[c50 n offset] = setExpParams(fitparams,testType,numConditions);

% calculate function
for iAdapt = 1:numConditions(1)
 for iAttend = 1:numConditions(2)
  % retrieve parameters for this EXPONENTIAL
  thisc50 = c50(iAdapt,iAttend);
  thisn = n(iAdapt,iAttend);
  thisoffset = offset(iAdapt,iAttend);
  thiscontrasts = squeeze(contrasts(iAdapt,iAttend,:));
  % now compute responses for this CRF
  fitfun(iAdapt,iAttend,:) = thisoffset+(thiscontrasts.^thisn)/thisc50;
 end
end


%%%%%%%%%%%%%%%%%%
%% setExpParams %%

%%%%%%%%%%%%%%%%%%
function [c50 n offset] = setExpParams(fitparams,testType,numConditions)

% init values
c50    = zeros(numConditions);
n      = zeros(numConditions);
offset = zeros(numConditions);

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
 case {'111'} % fit all params

  % make fitparams back into a nAttention x nAdaptation x nParams matrix;
  fitparams = reshape(fitparams,[numConditions 3]);
  % and grab the values
  c50    = fitparams(:,:,1);
  n      = fitparams(:,:,2);
  offset = fitparams(:,:,3);

 case {'101'} % fit all free c50 and ofset but one n
  c = 0;
  for iAdapt = 1:numConditions(1)
   for iAttend = 1:numConditions(2)
    c = c+1;
    c50(iAdapt,iAttend) = fitparams(c);
   end
  end
  
  n      = fitparams(13).*ones(numConditions);
  
  c = 0;
  for iAdapt = 1:numConditions(1)
   for iAttend = 1:numConditions(2)
    c = c+1;
    offset(iAdapt,iAttend) = fitparams(14+c-1);
   end
  end
  
 case {'003'} % allow offset to change across attention conditions
  % now we interpret the fitparams as the following
  % e.g. [c50 n offsetAttend1 offsetAttend2 offsetAttend3 offsetAttend4]

  % init values
  c50    = fitparams(1)*ones(numConditions);
  n      = fitparams(2)*ones(numConditions);
  
  for iAdapt = 1:numConditions(1)
   for iAttend = 1:numConditions(2)
    offset(iAdapt,iAttend) = fitparams(3+iAttend-1);
   end
  end

 case {'203'} % allow offset to change across attention conditions
  % and c50 across adaptation
  % now we interpret the fitparams as the following
  % e.g. [c50adapt1 c50adapt2 c50adapt3 n offsetAttend1 offsetAttend2 offsetAttend3 offsetAttend4]
 
  for iAdapt = 1:numConditions(1)
   for iAttend = 1:numConditions(2)
    c50(iAdapt,iAttend) = fitparams(1+iAdapt-1);
   end
  end

  n     = fitparams(4)*ones(numConditions);
  
  for iAdapt = 1:numConditions(1)
   for iAttend = 1:numConditions(2)
    offset(iAdapt,iAttend) = fitparams(4+iAttend);
   end
  end
  
  case {'213'} 
   % allow offset to change across attention conditions 
   % c50 across adaptation and n across al conditions:
   
  % now we interpret the fitparams as the following
  % e.g. [c50adapt1 c50adapt2 c50adapt3 
  %       nAdp1Att1 nAdp1Att2 nAdp1Att3 nAdp1Att4
  %       nAdp2Att1 nAdp2Att2 nAdp2Att3 nAdp2Att4
  %       nAdp3Att1 nAdp3Att2 nAdp3Att3 nAdp3Att4
  %       offsetAttend1 offsetAttend2 offsetAttend3 offsetAttend4]
 
  % c50
  for iAdapt = 1:numConditions(1)
   for iAttend = 1:numConditions(2)
    c50(iAdapt,iAttend) = fitparams(1+iAdapt-1);
   end
  end

  % n
  c=0;
  for iAdapt = 1:numConditions(1)
   for iAttend = 1:numConditions(2)
    c= c+1;
    n(iAdapt,iAttend) = fitparams(3+c);
   end
  end
  
  for iAdapt = 1:numConditions(1)
   for iAttend = 1:numConditions(2)
    offset(iAdapt,iAttend) = fitparams(15+iAttend);
   end
  end

 case {'302'}
  % allow offset to change across adaptation conditions
  % and c50 across attention n is fixed
  % now we interpret the fitparams as the following
  % e.g. [c50att1 c50att2 c50att3 c50att4 n offsetAdp1 offsetAdp2 offsetAdp3]

  for iAdapt = 1:numConditions(1)
   for iAttend = 1:numConditions(2)
    c50(iAdapt,iAttend) = fitparams(iAttend);
   end
  end

  n     = fitparams(5)*ones(numConditions);

  for iAdapt = 1:numConditions(1)
   for iAttend = 1:numConditions(2)
    offset(iAdapt,iAttend) = fitparams(5+iAdapt);
   end
  end

 case {'312'}
  % allow offset to change across adaptation conditions
  % and c50 across attention n is free to vary

  % now we interpret the fitparams as the following
  % e.g. [c50adapt1 c50adapt2 c50adapt3 c50adapt4
  %       nAdp1Att1 nAdp1Att2 nAdp1Att3 nAdp1Att4
  %       nAdp2Att1 nAdp2Att2 nAdp2Att3 nAdp2Att4
  %       nAdp3Att1 nAdp3Att2 nAdp3Att3 nAdp3Att4
  %       offsetAttend1 offsetAttend2 offsetAttend3]
  initparams = zeros(1,8);
  minfit = initparams;
  maxfit = initparams;

   % c50
   for iAdapt = 1:numConditions(1)
   for iAttend = 1:numConditions(2)
    c50(iAdapt,iAttend) = fitparams(iAttend);
   end
  end
  
  % n
  c=0;
  for iAdapt = 1:numConditions(1)
   for iAttend = 1:numConditions(2)
    c= c+1;
    n(iAdapt,iAttend) = fitparams(4+c);
   end
  end
  
  for iAdapt = 1:numConditions(1)
   for iAttend = 1:numConditions(2)
    offset(iAdapt,iAttend) = fitparams(16+iAdapt);
   end
  end

 otherwise
  keyboard
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% switchInitialParamsEX %%

%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [initparams minfit maxfit] = switchInitialParamsEX(testType,numConditions)
% c50 = fitparams(1);
% n = fitparams(2);
% offset = fitparams(3);

switch testType
 case {'111'} % fit all params
  % init the params to 0
  initparams = zeros([numConditions 3]);
  minfit = zeros([numConditions 3]);
  maxfit = zeros([numConditions 3]);
  
  % c50
  initparams(:,:,1) = .1;
  minfit(:,:,1) = 0;
  maxfit(:,:,1) = inf;
  % n
  initparams(:,:,2) = .1;
  minfit(:,:,2) = 0;
  maxfit(:,:,2) = inf;
  % offset
  initparams(:,:,3) = .15;
  minfit(:,:,3) = 0;
  maxfit(:,:,3) = inf;
  
  % make parameters into a linear array
  initparams = initparams(:);
  minfit = minfit(:);
  maxfit = maxfit(:);
  
 case {'101'} % fit all free c50 and ofset but one n

  initparams(1:12) = .1;
  minfit(1:12) = 0;
  maxfit(1:12) = inf;
  
  initparams(13) = .2;
  minfit(13) = 0;
  maxfit(13) = inf;
  
  initparams(14:26) = .1;
  minfit(14:26) = 0;
  maxfit(14:26) = inf;

 case {'003'} % allow offset to change across attention conditions
  % init the params to 0, there will be 2 values to fit c50 and n
  % for all curves and then 4 values for the offset to change across attention conditions
  % e.g. initparams = [c50 n offsetAttend1 offsetAttend2 offsetAttend3 offsetAttend4]
  
  % c50
  initparams(1) = .1; % one far all conditions
  minfit(1) = 0;
  maxfit(1) = inf;
  % n
  initparams(2) = .2;  % one far all conditions
  minfit(2) = 0;
  maxfit(2) = inf;
  % offset
  initparams(3:3+numConditions(2)-1) = .15;
  minfit(3:3+numConditions(2)-1) = 0;
  maxfit(3:3+numConditions(2)-1) = inf;

 case {'203'} % allow c50 to change across adaptation and offset across attention conditions
  % now we interpret the fitparams as the following
  % e.g. [c50Adapt1 c50Adapt2 c50Adapt3 n offsetAttend1 offsetAttend2 offsetAttend3 offsetAttend4]
  initparams = zeros(1,8);
  minfit = initparams;
  maxfit = initparams;
  
  % c50
  initparams(1:3) = .1;
  minfit(1:3)  = 0;
  maxfit(1:3)  = inf;

  % n
  initparams(4) = .2;
  minfit(4)  = 0;
  maxfit(4)  = inf;
  
  % offset
  initparams(5:8) = .15;
  minfit(5:8)  = 0;
  maxfit(5:8)  = inf;
  
  case {'302'} 
   % allow offset to change across adaptation conditions
  % and c50 across attention
  % now we interpret the fitparams as the following
  % e.g. [c50att1 c50att2 c50att3 c50att4 n offsetAdp1 offsetAdp2 offsetAdp3]
   initparams = zeros(1,8);
   minfit = initparams;
   maxfit = initparams;

   % c50
   initparams(1:4) = .1;
   minfit(1:4)  = 0;
   maxfit(1:4)  = inf;

   % n
   initparams(5) = .2;
   minfit(5)  = 0;
   maxfit(5)  = inf;

   % offset
   initparams(6:8) = .15;
   minfit(6:8)  = 0;
   maxfit(6:8)  = inf;

 case {'213'}
  % allow offset to change across attention conditions 
  % c50 across adaptation and n across al conditions:
   
  % now we interpret the fitparams as the following
  % e.g. [c50adapt1 c50adapt2 c50adapt3 
  %       nAdp1Att1 nAdp1Att2 nAdp1Att3 nAdp1Att4
  %       nAdp2Att1 nAdp2Att2 nAdp2Att3 nAdp2Att4
  %       nAdp3Att1 nAdp3Att2 nAdp3Att3 nAdp3Att4
  %       offsetAttend1 offsetAttend2 offsetAttend3 offsetAttend4]
  initparams = zeros(1,19);
  minfit = initparams;
  maxfit = initparams;

  % c50
  initparams(1:3) = .1;
  minfit(1:3)  = 0;
  maxfit(1:3)  = inf;

  % n
  initparams(4:15) = .2;
  minfit(4:15)  = 0;
  maxfit(4:15)  = inf;

  % offset
  initparams(16:19) = .15;
  minfit(16:19)  = 0;
  maxfit(16:19)  = inf;
  
  
  case {'312'} 
   % allow offset to change across adaptation conditions
   % and c50 across attention n is free to vary
   
   % now we interpret the fitparams as the following
   % e.g. [c50adapt1 c50adapt2 c50adapt3 c50adapt4
   %       nAdp1Att1 nAdp1Att2 nAdp1Att3 nAdp1Att4
   %       nAdp2Att1 nAdp2Att2 nAdp2Att3 nAdp2Att4
   %       nAdp3Att1 nAdp3Att2 nAdp3Att3 nAdp3Att4
   %       offsetAttend1 offsetAttend2 offsetAttend3]
   initparams = zeros(1,8);
   minfit = initparams;
   maxfit = initparams;

   % c50
   initparams(1:4) = .1;
   minfit(1:4)     = 0;
   maxfit(1:4)     = inf;

   % n
   initparams(5:16) = .2;
   minfit(5:16)     = 0;
   maxfit(5:16)     = inf;

   % offset
   initparams(17:19) = .15;
   minfit(17:19)     = 0;
   maxfit(17:19)     = inf;
   
 otherwise
  keyboard
end



%%%%%%%%%%%%%%
%% testnaka %%

%%%%%%%%%%%%%%
function [fitfun Rmax c50 n offset] = testnaka(fitparams,contrasts,testType,numConditions)

% set params of the naka rushton depending on the type of test being
% performed:
[Rmax c50 n offset] = setNakaParams(fitparams,testType,numConditions);

% calculate function
for iAdapt = 1:numConditions(1)
 for iAttend = 1:numConditions(2)
  % retrieve parameters for this CRF
  thisRmax = Rmax(iAdapt,iAttend);
  thisc50 = c50(iAdapt,iAttend);
  thisn = n(iAdapt,iAttend);
  thisoffset = offset(iAdapt,iAttend);
  thiscontrasts = squeeze(contrasts(iAdapt,iAttend,:));
  % now compute responses for this CRF
  fitfun(iAdapt,iAttend,:) = thisRmax*(thiscontrasts.^thisn)./(thiscontrasts.^thisn+thisc50.^thisn)+thisoffset;
 end
end



%%%%%%%%%%%%%%%%%%%
%% setNakaParams %%

%%%%%%%%%%%%%%%%%%%
function [Rmax c50 n offset] = setNakaParams(fitparams,testType,numConditions)
% init values
Rmax   = zeros(numConditions);
c50    = zeros(numConditions);
n      = zeros(numConditions);
offset = zeros(numConditions);

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
 case {'1111'} % fit all params
  % make fitparams back into a nAttention x nAdaptation x nParams matrix;
  fitparams = reshape(fitparams,[numConditions 4]);
  % and grab the values
  Rmax   = fitparams(:,:,1);
  c50    = fitparams(:,:,2);
  n      = fitparams(:,:,3);
  offset = fitparams(:,:,4);
 case {'0003'} % allow offset to change across attention conditions
  % now we interpret the fitparams as the following
  % e.g. [Rmax c50 n offsetAttend1 offsetAttend2 offsetAttend3 offsetAttend4]
  Rmax(:)  = fitparams(1);
  c50(:)   = fitparams(2);
  n(:)     = fitparams(3);
  for iAttend = 1:numConditions(2)
   offset(:,iAttend) = fitparams(3+iAttend);
  end
 case {'0203'} % allow offset to change across attention conditions
  % now we interpret the fitparams as the following
  % e.g. [Rmax c50Adapt1 c50Adapt2 c50Adapt3 n offsetAttend1 offsetAttend2 offsetAttend3 offsetAttend4]
  Rmax(:)  = fitparams(1);
  keyboard
  for iAdapt = 1:numConditions(1)
   c50(iAdapt,:) = fitparams(2+iAdapt-1);
  end
  
  n(:)     = fitparams(5);
  
  for iAttend = 1:numConditions(2)
   offset(:,iAttend) = fitparams(6+iAttend-1);
  end
 otherwise
  keyboard
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% switchInitialParamsNK %%

%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [initparams minfit maxfit] = switchInitialParamsNK(testType,numConditions)
% Rmax = fitparams(1);
% c50 = fitparams(2);
% n = fitparams(3);
% offset = fitparams(4);

switch testType
 case {'1111'} % fit all params
  % init the params to 0
  initparams = zeros([numConditions 4]);
  minfit = zeros([numConditions 4]);
  maxfit = zeros([numConditions 4]);
  % Rmax
  initparams(:,:,1) = 1;
  minfit(:,:,1) = .5;
  maxfit(:,:,1) = inf;
  % c50
  initparams(:,:,2) = .5;
  minfit(:,:,2) = 0;
  maxfit(:,:,2) = 1;
  % n
  initparams(:,:,3) = .2;
  minfit(:,:,3) = 0;
  maxfit(:,:,3) = inf;
  % offset
  initparams(:,:,4) = .15;
  minfit(:,:,4) = 0;
  maxfit(:,:,4) = inf;
  % make parameters into a linear array
  initparams = initparams(:);
  minfit = minfit(:);
  maxfit = maxfit(:);

 case {'0003'} % allow offset to change across attention conditions
  % init the params to 0, there will be 3 values to fit Rmax, c50 and n
  % for all curves and then 4 values for the offset to change across attention conditions
  % e.g. [Rmax c50 n offsetAttend1 offsetAttend2 offsetAttend3 offsetAttend4]
  initparams = [3 numConditions(2)];
  minparams = [3 numConditions(2)];
  maxparams = [3 numConditions(2)];
  % Rmax
  initparams(1) = 1;
  minfit(1) = .5;
  maxfit(1) = inf;
  % c50
  initparams(2) = .5;
  minfit(2) = 0;
  maxfit(2) = 1;
  % n
  initparams(3) = .2;
  minfit(3) = 0;
  maxfit(3) = inf;
  % offset
  initparams(4:4+numConditions(2)-1) = .15;
  minfit(4:4+numConditions(2)-1) = 0;
  maxfit(4:4+numConditions(2)-1) = inf;

 otherwise
  keyboard
end



%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   setDataStructure   %%

%%%%%%%%%%%%%%%%%%%%%%%%%%
function r = setDataStructure(data,meanC, subplotInfo,visualArea)

% chose adaptation, contrast and attention condtion:
[a thisCRF  thisSTE thisContrast plotInfo] = chooseCondition(subplotInfo,visualArea);

% make arrays for pedestal contrasts and crfs:
for i = 1:length(data)
 for j = 1:length(thisContrast)
  r.pedestals(i,j,:) = meanC{a(i)}.(thisContrast{j})./100;
  r.use_crf(i,j,:) = data{a(i)}.(thisCRF{j});
  r.use_crf_ste(i,j,:) = data{a(i)}.(thisSTE{j});
 end
end

r.used_crf = 'all';
r.plotInfo = plotInfo;
r.plotInfo.figureTitle = visualArea;


%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   chooseCondition    %%

%%%%%%%%%%%%%%%%%%%%%%%%%%
function [a thisCRF  thisSTE thisContrast plotInfo] = chooseCondition(subplotInfo,visualArea)

% set up adaptation indexes
% change this when turning this code for individual observers
a = [1 2 3];

thisCRF = {'A_amplitude' 'Dt_amplitude' 'Dnt_amplitude' 'U_amplitude'};
thisContrast = {'quantileAcontrast' 'quantileDcontrast' 'pedC' 'pedC'};
thisSTE = {'A_ste' 'Dt_ste' 'Dnt_ste' 'U_ste'};

% figure formatting info:
plotInfo.MarkerEdgeColor = [.4 .4 .4];
plotInfo.thisSymbol = {'o-' 'o-' 'o-' 'o-'};

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
plotInfo.XYColor = [.4 .4 .4];
plotInfo.Fsize = 6;
plotInfo.MarkerSize = 5;
plotInfo.LineWidth = 1;
plotInfo.TickLength = [0.025 .01];
plotInfo.PlotBoxAspectRatio = [1 1 1];
plotInfo.YLim = [0 1.3];
plotInfo.Xlim = [.001,1];
plotInfo.plotPeds = [plotInfo.Xlim(1) 0.0175 0.035 0.07 0.14 0.28 0.56 0.84];

% set up subplotinfo
for i = 1:length(a)
 switch visualArea
  case {'v1' 'V1'}
   plotInfo.subplotInfo{i} = [4,3,i];
  case {'v2' 'V2'}
   plotInfo.subplotInfo{i} = [4,3,i+3];
  case {'v3' 'V3'}
   plotInfo.subplotInfo{i} = [4,3,i+6];
  case {'v4' 'V4'}
   plotInfo.subplotInfo{i} = [4,3,i+9];
  otherwise
 end
end



