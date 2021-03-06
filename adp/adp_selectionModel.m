% adp_selectionModel.m
%
%        $Id:$ 
%      usage: adp_selectionModel(CRF,CDF,<>)
%         by: justin gardner
%       date: 04/22/11
%    purpose: Fit the selection model code.
%             Pestilli, Carrasco, Heeger & Gardner Neuron 2011
%
% INPUTS:
%  CRF - a structure containign the data and the conditions for the contrast
%        brain response, the contrast response function (CRF). 
%
%        CRF is a structure that has two different contrast response functions (cued and uncued)
%        for each condition you are testing. The names of the conditions are arbitrary but
%        each condition must contain a cued and uncued CRF. So, for example if you have a 
%        focal cue condition then you would have:         
%        - CRF.focal.cued   (target) 
%        - CRF.focal.uncued (non-target)
%                
%        There should also be a field to tell whether the cued and uncued functions are the
%        same (i.e. for a distributed cue condition or different):
%        - CRF.focal.sameCuedUncued       = false;  
%        - CRF.distributed.sameCuedUncued = true; % in the distributed cue condition target 
%                                                 % and non-target are the same.
%        For each of the cued/uncued fields above there should be fields for data and fit
%        which contain the contrasts and responses. For example, something like:
%        - CRF.focal.cued.data.c   = [0 0.1 0.5 1]; % Contrast of the stimulus on the 
%                                                   % x-axis of the CRF)
%        - CRF.focal.uncued.data.r = [0.01 0.3 0.7 1.10]; % The actual CRF BOLD responses
%  
%  CDF - a structure containing the data and the condition for the
%        behavior, the contrast discrimination function (CDF)
%  
%        CDF is a similar struncture containing contrast discrimination data. So, continuing
%        from the exmaple above, if you have a focal cue condition you would have fields for
%        both data and fit that contain contrast discrimination data:
%        - CDF.focal.data.c     = [0.01 0.5 0.7];  % Pedestal contrast 
%        - CDF.focal.data.t     = [0.04 0.08 0.12];% Threshold, delta-c
%        - CDF.focal.data.tste  = [0.01 0.02 0.01];% ste of threshold
%        - CDF.focal.data.d     = {[0.0 0.01 0.02], ...
%                                  [0.01 0.5 0.7],  ...
%                                  [0.01 0.5 0.7]};
%        - CDF.focal.data.nCued = 1;
% 
%        Where c are the pedestal contrasts, t are the measured thresholds and d is the
%        distractor contrasts in a cell array. The cued field contains how many targets
%        were cued. The tste field contains the standard error of the CDF data (for display)
%
%  dispCDF  - Whether to dispaly (1) or not (0) the contrast discrimination functions
%  dispCRF  - Whether to display (1) of not (0) the contrast responst functions
%  dispCond - Whoch condtion to display: 
%             (a) dispCRF(s,m,CRF);
%             (b) dispCDF(s,m,CDF);
%  fitCDF   - Whether to fit or no tthe contrast discrimination function
%  fitCRFCDF- Whether to fit both contrast response function and contrast
%             discrimination function.
%  verbose  - Show progress, 1,0
%  dispFig  - Make  figure?
%  fitCDF   - fit the CDF
%  dispIter - display each  iteraton 1/0
%  fitCRFCDF- fit the CRF and the CDF (1/0)
%  recompute- recompute model?
%  saveName - file name to save
%  crossValidate - cross valudate results
%  saveIter - save each iteration 1/0
%  dispSavedIter - display the saved iterations 0/1
%  reslog   - 
%  fixedParams - 
%  initParams  -
%  reuseFig    - 
%  noFig       -
%  savePathDir=[]'
%
% OUTPUTS:
%    retval - what is this justin?
%        
% NOTEs:
%   (1) You can skip the data field if you are not simultaneously fitting
%   CRF and CDF. Also, the fit contrast and response should have fairly finely
%   spaced contrast/response values since the code will interpolate from these
%   values to get the response for any particular contrast. 
%   (2) For display purposes you can also supply the standard errors of the 
%   response by adding a rste field: CRF.focal.uncued.data.rste
%


function retval = adp_selectionModel(CRF,CDF,varargin)

retval = [];
% check arguments
if nargin < 2
  help adp_selectionModel
  return
end

% set up modelParams depending on input arguments
[s vars] = setSystemParams(varargin);
m = setModelParams(s,vars);

% validate the input arguments
[tf m CRF] = checkDataStructures(CRF,CDF,s,m);
if ~tf,return,end

% run various things depending on todo list
if todo(s,'dispSavedIter'), dispSavedIter(s,m);end
if todo(s,'dispCRF'),       dispCRF(s,m,CRF);end
if todo(s,'dispCDF'),       dispCDF(s,m,CDF);end
if todo(s,'dispCond'),      dispCond(s,m,CDF);end
if todo(s,'fitCDF'),        retval = fitCDF('CDF',s,m,CRF,CDF);end
if todo(s,'fitCRFCDF'),     retval = fitCDF('CRFCDF',s,m,CRF,CDF);end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Model fits
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%
%    fitCDF    %
%%%%%%%%%%%%%%%%
function retval = fitCDF(fitType,s,m,CRF,CDF)

% need to split the data for cross validation purposes
if s.crossValidate
  % keep original CDF
  originalCDF = CDF;
  % get the even crossval
  [m CDF] = getCDFCrossVal(m,originalCDF,'even');
end

% set parameter to show that we are fitting the CRF or not
if strcmp(fitType,'CRFCDF')
  m.fitCRF = true;
else
  % this is a CDF fit only
  m.fitCRF = false;
  % see if we need to compute a fit
  if ~isfield(CRF.(m.conditions{1}).cued,'fit')
    disp(sprintf('(adp_selectionModel) No preexisting fit found. Recomputing fit'));
    [CRF CRFParams m] = fitCRF(s,m,CRF,1);
  end
  % set the interp function to use the fit
  m.CRFInterpFun = 'fit';
end
  
% show save filename
disp(sprintf('Filename: %s',getSaveName(fitType,s,m)));
%dispHeader(sprintf('Fit CDF: %s',fitType));

% first decide if we are actually computing the fit, or using fixed/preloaded params
if (s.recompute || ~isComputed(fitType,s,m)) && isempty(s.fixedParams)
  if strcmp(fitType,'CRFCDF')
    % fit the CRF functions
    disppercent(-inf,'(adp_selectionModel:fitCRFCDF) Getting initial fit params for CRF');
    % fit the CRF (if we have init parameters, then we setup the
    % structure to have the correct paramsIndex, but don't actually do fit
    [CRF CRFParams m] = fitCRF(s,m,CRF,isempty(s.initParams) && ~m.emptyCRF);
    disppercent(inf);
    % display one line header
    %dispHeader;
  else
    CRFParams = [];
  end

  % get the initial K
  if strcmp(m.poolingRule,'sensitivity')
    initK = [];
  else
    initK = 10;
  end

  % fixed k-value
  if m.fixedK,initK = [];end
  
  % initial sigma value
  initSigma = 0.017;

  % set up params array
  if m.sigmaForEachCond
    initParams = [repmat(initSigma,1,m.numConditions) initK CRFParams];
  elseif m.sigmaCuedUncued
    initParams = [initSigma initSigma initK CRFParams];
  else
    initParams = [initSigma initK CRFParams];
  end

  % override, if we have initParams
  if ~isempty(s.initParams),initParams = s.initParams;end

  % if we are using an output nonlinearity, compute the range (how to normalize BOLD responses).
  if ~any(m.outputNonlin==[0 1]) && isempty(m.outputNonlinRange)
    m.outputNonlinRange = getRangeOfCRF(s,m,CRF);
  end
  
  % if outputNonlinearity is a fit parameter
  if m.outputNonlinFit
    % init outputNonlin in params
    initParams(end+1) = m.outputNonlin;
    m.outputNonlinFitParamNum = length(initParams);
  end
  
  % use Nelder-Mead method to do fit
  [fitParams fval exitflag output] = fminsearch(@computeModelResidualUsingThresholds,initParams,s.optimParams,s,m,CRF,CDF);

  % get the model thresholds
  t = computeModelThreshold(fitParams,s,m,CRF,CDF);

  % save 
  saveState(fitType,s,m,CRF,CDF,fitParams,t);

% passed in fixed params
elseif ~isempty(s.fixedParams)
  % use the passed in params
  fitParams = s.fixedParams;
  % compute the fit on the CRF, first setup the CRF with the
  % right paramaterIndexes by calling fitCRF with doFIT = 0
  if strcmp(fitType,'CRFCDF')
    [CRF CRFParams m] = fitCRF(s,m,CRF,0);
  end
  % compute thresholds
  t = computeModelThreshold(fitParams,s,m,CRF,CDF);
else
  % load precomputed values
  [sold m CRF CDF fitParams t] = loadState(fitType,s,m);
  % get the model thresholds
  if isempty(t)
    t = computeModelThreshold(fitParams,s,m,CRF,CDF);
  end
end

% get cross val data if we are doing this
if s.crossValidate
  % get the odd crossval
  [m CDF] = getCDFCrossVal(m,originalCDF,'odd');
  t = computeModelThreshold(fitParams,s,m,CRF,CDF);
end

% display model fit.
dispModelThresholds(s,m,CRF,CDF,t,fitParams,'');
dispCRF(s,m,CRF,1,fitParams);

% calculate statistics
stats = calcStats(s,m,CRF,CDF,t,fitParams);

% compute smooth fit of CRF for passing back the fit
% to the calling function.
CRF = computeSmoothFit(s,m,CRF,fitParams);

% package up return arguments
retval.params = fitParams;
retval.stats = stats;
retval.CRF = CRF;
%YUKO we want these params too - when refitting for smoothCDFs, we want to see the output thresholds.
retval.s = s;
retval.m = m;
retval.CDF = CDF;
retval.t = t;
retval.aux = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    computeModelResidualUsingThresholds    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function residual = computeModelResidualUsingThresholds(params,s,m,CRF,CDF)

% so we can stop by hitting shift and ctrl at the same time and check out / modify params
%if all(mglGetKeys([57 60])),keyboard,end

% get data thresholds
for iCond = 1:m.numConditions
  data(iCond,:) = CDF.(m.conditions{iCond}).data.t;
  c(iCond,:) = CDF.(m.conditions{iCond}).data.c;
end

% get the model thresholds
t = computeModelThreshold(params,s,m,CRF,CDF);

% compute difference in log units as residual
residual = logb(data(:),s.reslog) - logb(t(:),s.reslog);
residual = sqrt(sum(residual.^2));

if m.normResByVar %YUKO 1/3
% residual is normalized by how much variance there is in the data
  residual = var(residual)/var(logb(data(:),s.reslog));
end

% if we are fitting CRF, then compute that residual as well
if m.fitCRF && ~m.emptyCRF
  iCRF = 1;
  for iCond = 1:m.numConditions
    for iCued = {'cued','uncued'}
      if CRF.(m.conditions{iCond}).sameCuedUncued && strcmp(iCued,'uncued')
        %skip
      else
        thisCRF = CRF.(m.conditions{iCond}).(iCued{1});
        % get actual response
        r{iCRF} = thisCRF.data.r;
        % get fit response
        fitResponse{iCRF} = getResponse(s,m,params,thisCRF,thisCRF.data.c)';
        iCRF = iCRF+1;
      end
    end
  end
  % now convert to an array of all measured responses against all fitted responses
  r = cell2mat(r);
  fitResponse = cell2mat(fitResponse);
  % and add to the residual
  res2 = 2*sqrt(sum((r-fitResponse).^2));
if m.normResByVar %YUKO 2/3
  %  residual is normalized by how much variance there is in the data
  res2 = var(r-fitResponse)/var(r);
end
  resString = sprintf('%s+%s',mynum2str(residual,'sigfigs=3','compact=1'),mynum2str(res2,'sigfigs=3','compact=1'));
% YUKO 3/3: weight residuals to bias fits toward CRF/CDF (these weights are usually set to 1)
  residual = m.CDFres*residual + m.CRFres*res2; 
  resString = sprintf('%s=%s',mynum2str(residual,'sigfigs=3','compact=1'),resString);
else
  resString = sprintf('%s',mynum2str(residual,'sigfigs=3','compact=1'));
end

% display current parameter settings
if s.verbose>0
  disp(sprintf('(adp_selectionModel:computeModelResidualUsingThresholds)\nParameters: %s (res:%s)',mynum2str(params,'sigfigs=3','compact=1'),resString));
end

% save each iteration of the loop
if s.saveIter
  global selectionModelSaveIter;
  % save the state
  saveState(sprintf('CRFCDF_iter%05i',selectionModelSaveIter),s,m,CRF,CDF,params,t,resString);
  % update the iter count
  selectionModelSaveIter = selectionModelSaveIter+1;
end

% display figure of model fit
if s.dispIter
  dispModelThresholds(s,m,CRF,CDF,t,params,resString);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    computeModelThreshold    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% computes the thresholds for the model parameters.
% Note that this uses an fminsearch, since simulating trials
% returns percent correct for the parameters given a threshold increment
% when what we need is the threshold increment that gives
% a certain percent correct (i.e. 76.05% for dprime=1)
function t = computeModelThreshold(params,s,m,CRF,CDF)

% This function computes the model thresholds that give approxiametly
% d'=1 performance (p=0.7605)

% cycle over experimental conditions (e.g. focal, distributed)
if s.verbose>0,disppercent(-inf,sprintf('(adp_selectionModel:computeModelThresholds) Computing model thresholds for parameters:\n %s',mynum2str(params,'compact=1')));end
for iCond = 1:m.numConditions
  % get condition name and relevant info
  condName = m.conditions{iCond};
  thisCRF = CRF.(condName);
  thisCDF = CDF.(condName).data;
  nPed = length(thisCDF.c);
  % cycle over pedestal contrasts computing the thresholds from the model
  % parameters. We do this in two for loops quickly precomputing some values
  % for the simulation and then running the fminsearch as a separate
  % for loop - so that this slow loop can be run as a parfor
  for iPed = 1:nPed
    % get target and distractor contrasts for this pedestal
    distractorContrasts{iPed} = thisCDF.d{iPed};
    nDistractors{iPed} = length(distractorContrasts{iPed});
    targetContrast{iPed} = thisCDF.c(iPed);
    % get all the contrasts for this pedestal
    contrasts{iPed} = [targetContrast{iPed} distractorContrasts{iPed}];
    % now put the necessary number of cued CRF functions into CRFs
    for iTarg = 1:thisCDF.nCued
      CRFs{iPed}{iTarg} = thisCRF.cued;
    end
    % and fill the rest with the uncued CRFs
    for iTarg = thisCDF.nCued+1:nDistractors{iPed}+1
      CRFs{iPed}{iTarg} = thisCRF.uncued;
    end
  end
  % create vector containing whether each location is cued or not
  cuedUncued = zeros(1,m.nStimuli);
  for iTarg = 1:thisCDF.nCued
    cuedUncued(iTarg) = 1;
  end
  % now search for the best threshold contrast that gives a d' use parfor loop
  % to spped this slow step up.
  parfor iPed = 1:nPed
    % now do a non-linear fit to find the expected thresholds
    [t(iCond,iPed) fval exitflag output] = fminsearch(@computeModelThresholdResidual,thisCDF.t(iPed),s.optimParams,s,m,params,iCond,contrasts{iPed},CRFs{iPed},cuedUncued);
  end
  % check for any failed fits
  failedFit = find(thisCDF.t(:)' == t(iCond,:));
  if ~isempty(failedFit)
    % go through each failed fit
    for iPed = failedFit
      % some parameters
      maxTries = 100;thisTry = 1;startThreshold = thisCDF.t(iPed);
      % as long as we have not exceeded the maximum number of tries, if the fminsearch
      % has returned the same threshold as the startThreshold than it has failed, so continue trying
      while (startThreshold == t(iCond,iPed)) && (thisTry < maxTries)
	% set the starting threshold to some random value
	startThreshold = rand;
	% now do a non-linear fit to find the expected thresholds
	[t(iCond,iPed) fval exitflag output] = fminsearch(@computeModelThresholdResidual,startThreshold,s.optimParams,s,m,params,iCond,contrasts{iPed},CRFs{iPed},cuedUncued);
	% update number of tries
	thisTry = thisTry + 1;
      end
      if thisTry >= maxTries
	disp(sprintf('(adp_selectionModel:computeModelThreshold) Failed threshold fit on %s:%i',condName,iPed));
      end
    end
  end
  % percent done
  if s.verbose>0,disppercent(iCond/m.numConditions);end
end
if s.verbose>0,disppercent(inf);end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    computeModelThresholdResidual    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function residual = computeModelThresholdResidual(threshold,s,m,params,condNum,contrasts,CRFs,cuedUncued)

% if threshold is less than 0, than give a residual of inf
if threshold < 0,residual = inf;return;end

% simulate the trial with these conditions
p = simTrial(s,m,params,condNum,contrasts,CRFs,threshold,cuedUncued);

% if p returns 0, then the parameters are way out of whack - set to nan
if p == 0
  disp(sprintf('(computeModelThresholdResidual) !!! simTrial returned p==0 params out of range? !!!'));
end

% compute residual relative to a dprime of 1 (76.06%)
residual = sqrt((p-0.7605).^2);

%%%%%%%%%%%%%%%%%%
%    simTrial    %
%%%%%%%%%%%%%%%%%%
function p = simTrial(s,m,params,condNum,contrasts,CRFs,deltaC,cuedUncued)

% returns the percent correct of the model for a trial condition given model parameters, params
% contrasts are the contrasts of each stimulus location. e.g. [0.1 0.2 0.3 0.2];
% CRFs is the contrast response functions for each of those locations {CRF1 CRF2 CRF3 CRF4}
% params are the params for the model: e.g. [sigma k]
% m.r is a set of random distributions for each locataion and for each temporal interval
% i.e. it should be like {{Rloc1int1 Rloc2int1 Rloc3int1 Rloc4int1}{Rloc1int2 Rloc2int2 Rloc3int2 Rloc4int2}}
% the distribution are fixed throughout the simulation so that we can use nonlinear search algorithms
% which depend on getting the same answer each time.

% compute response distributions for target location
if ~m.randDistractors
  for i = 1:length(contrasts)
    % Get the contrast of the stimulus
    thisContrast = contrasts(i);
    % get mean responses, for each contrast in the stimulus from its corresponding CRF
    r(1,i) = getResponse(s,m,params,CRFs{i},thisContrast,true);
    if i == 1
      % threshold contrast added at target location  
      r(2,i) = getResponse(s,m,params,CRFs{i},thisContrast+deltaC,true);
    else
      % no threshold contrast added for distractors
      r(2,i) = r(1,i);
    end
    
    % Get sigma that we are using
    sigma = getSigma(m,params,condNum,cuedUncued(i));

    % and compute response distributons
    rdist{1,i} = m.r{1,i}*sigma+r(1,i);
    rdist{2,i} = m.r{2,i}*sigma+r(2,i);
  end
else
  % deal with the case in which distractors are randomized trial by trial
  % the first contrast is the target, so it is always the fixed contrast
  % and has a deltaC added to the second interval
  r(1,1,1:m.n) = getResponse(s,m,params,CRFs{1},contrasts(1),true);
  r(1,2,1:m.n) = getResponse(s,m,params,CRFs{1},contrasts(1)+deltaC,true);
  for i = 2:length(contrasts)
    % get a random contrast for one of the distractors
    r(i,1,:) = getResponse(s,m,params,CRFs{i},contrasts(m.randDistractorsDist(i-1,:)),true)';
    % the second interval is the same as the first interval
    r(i,2,:) = r(i,1,:);
  end
  for i = 1:m.nStimuli
    % Get sigma that we are using
    sigma = getSigma(m,params,condNum,cuedUncued(i));
    % compute the distributions
    rdist{1,i} = m.r{1,i}*sigma+squeeze(r(i,1,:))';
    rdist{2,i} = m.r{2,i}*sigma+squeeze(r(i,2,:))';
  end
end

% now pump through the pooling rule
pooled = poolResponses(s,m,params,rdist);

% display response distributions if asked for
if s.dispFig,dispResponseDist(s,m,params,rdist,pooled);end

% compute the performance from the pooled distribution
p = computeAreaUnderROC(pooled);

%%%%%%%%%%%%%%%%%%%%%
%    getResponse    %
%%%%%%%%%%%%%%%%%%%%%
function response = getResponse(s,m,params,CRF,contrast,outputNonlin)

% only use the output nonlin when simulating trials
if nargin < 6
  outputNonlin = false;
end

% compute response from contrast, depends on whether the model calls
% for fitting a function to the contrast-response data or using
% a precomputed fit. 
switch (m.CRFInterpFun)
  case {'fit'}
   response = interpolateResponse(CRF.fit.c,CRF.fit.r,contrast);
  case {'loglinear'}
   % get parameters
   p = getCRFParams(m,params);
   % now get the params for this CRF - this index is set in
   % the fitCRFpoly function above.
   p = p(CRF.fit.paramsIndex);
   p = p(:);
   % and compute the response for the contrast.
   response = makePolyMatrix(contrast,2,1) * p;
  case {'linear'}
   % get parameters
   p = getCRFParams(m,params);
   % now get the params for this CRF - this index is set in
   % the fitCRFpoly function above.
   p = p(CRF.fit.paramsIndex);
   p = p(:);
   % and compute the response for the contrast.
   response = makePolyMatrix(contrast,2) * p;
  case {'poly'}
   % get parameters
   p = getCRFParams(m,params);
   % now get the params for this CRF - this index is set in
   % the fitCRFpoly function above.
   p = p(CRF.fit.paramsIndex);
   p = p(:);
   % and compute the response for the contrast. Note that
   % the order of the polynomial is the first CRFInterpFunParam
   response = makePolyMatrix(contrast,m.CRFInterpFunParams{1}) * p;
  case {'nakarushton'}
   % get parameters
   p = getCRFParams(m,params);
   % now get the params for this CRF - this index is set in
   % the fitCRFpoly function above.
   p = p(CRF.fit.paramsIndex);
   p = p(:);
   %account for fixedN if necessary
   if isfield(m,'fixedN')
     if any(m.fixedN)
       p = [p(1:2); m.fixedN; p(3)];
     end
   end
   % and compute the response for the contrast. Note that
   % the order of the polynomial is the first CRFInterpFunParam
   response = nakaRushton(contrast(:),p);
  case {'exponential'}
   % get parameters
   p = getCRFParams(m,params);
   % now get the params for this CRF - this index is set in
   % the fitCRFpoly function above.
   p = p(CRF.fit.paramsIndex);
   p = p(:);
   % and compute the response for the contrast. Note that
   % the order of the polynomial is the first CRFInterpFunParam
   response = exponential(contrast(:),p);
 otherwise
  disp(sprintf('(adp_selectionModel:getResponse) Unknown InterpFun type: %s',m.CRFInterpFun));
  keyboard
end

% apply output nonlinearity if applicable
if outputNonlin
%  disp(sprintf('(adp_selectionModel:getResponse) Applying output nonlinearity: %f',m.outputNonlin));
  % normalize response range
  response = applyOutputNonlin(m,response,params);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    applyOutputNonlin    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%
function response = applyOutputNonlin(m,response,params)

% get the output nonlin from the parameters
outputNonlin = getOutputNonlin(m,params);

% apply the output nonlinearity if needed
if ~isnan(outputNonlin)
  response = (response-m.outputNonlinRange(1))/diff(m.outputNonlinRange);
  response(response<0) = 0;
  response = response.^m.outputNonlin;
end

%%%%%%%%%%%%%%%%%%%%%%%
%    poolResponses    %
%%%%%%%%%%%%%%%%%%%%%%%
function pooled = poolResponses(s,m,params,rdist)

switch (m.poolingRule)
  case {'exponential'}
   k = getK(m,params);

   % pool over 1st interval
   pooled{1} = exponentSum(cell2mat({rdist{1,:}}'),k);
   % pool over 2st interval
   pooled{2} = exponentSum(cell2mat({rdist{2,:}}'),k);
 case {'sensitivity'}
  pooled{1} = rdist{1,1};
  pooled{2} = rdist{2,1};
 otherwise
  disp(sprintf('(adp_selectionModel:poolResponses) Unknown pooling rule: %s',m.poolingRule));
  keyboard
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    computeAreaUnderROC    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function areaUnderROC = computeAreaUnderROC(r)

nSamples = size(r{1}(:),1);

explicitComputation = false;
if explicitComputation
  % number of criterion points to compute the ROC curve from -- the more the
  % finer the ROC curve will be sampled
  nCriterion = 10000; % make sure this is an even number

  % now compute min/max and step size for stepping through criterion
  minResponse = min([r{1}(:)' r{2}(:)']);
  maxResponse = max([r{1}(:)' r{2}(:)']);
  responseStep = (maxResponse-minResponse)/(nCriterion-1);
  criterions = minResponse:responseStep:maxResponse;

  % Now get hits and false alarams for ROC function.
  % i.e. compute the proportion of hist and false alarms as a function
  % of the criterion (i.e. each element in hits/false alarms is for
  % a criterion as specified in the criterions array). 
  hits = max(0,1-cumsum(histc(r{2},criterions,2)/nSamples,2));
  falseAlarms = max(0,1-cumsum(histc(r{1},criterions,2)/nSamples,2));

  % compute area under ROC
  areaUnderROC = abs(trapz(falseAlarms(:),hits(:)));
else
  % performance can be comupted as a difference between trials instead
  % of area under the ROC
  areaUnderROC = sum((r{2}-r{1})>0,2)/nSamples;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Cross-val functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%
%    getCDFCrossVal    %
%%%%%%%%%%%%%%%%%%%%%%%%
function [m crossvalCDF] = getCDFCrossVal(m,CDF,whichOne)

for iCond = 1:m.numConditions
  % select which peds to keep for this crossval
  switch (whichOne)
   case {'even','e',1}
    whichPeds = 2:2:length(CDF.(m.conditions{iCond}).data.c);
   case {'odd','o',2}
    whichPeds = 1:2:length(CDF.(m.conditions{iCond}).data.c);
   otherwise
    disp(sprintf('(adp_selectionModel:getCDFCrossVal) Unknown crossvalidation set'));
    keyboard
  end
  % get all fields in data
  dataFields = fieldnames(CDF.(m.conditions{iCond}).data);
  crossvalFields = {'c','d','t','tste'};
  % copy fields, but for any of the crossval fileds only copy part of the data
  for iField = 1:length(dataFields)
    if any(strcmp(dataFields{iField},crossvalFields))
      % copy only crossval part
      crossvalCDF.(m.conditions{iCond}).data.(dataFields{iField}) = CDF.(m.conditions{iCond}).data.(dataFields{iField})(whichPeds);
    else
      % copy whole structure
      crossvalCDF.(m.conditions{iCond}).data.(dataFields{iField}) = CDF.(m.conditions{iCond}).data.(dataFields{iField});
    end
  end
  % now copy rest of fields
  dataFields = fieldnames(CDF.(m.conditions{iCond}));
  for iField = 1:length(dataFields)
    if ~strcmp(dataFields{iField},'data')
      crossvalCDF.(m.conditions{iCond}).(dataFields{iField}) = CDF.(m.conditions{iCond}).(dataFields{iField});
    end
  end
  % update number of pedestals
  m.numPeds(iCond) = length(whichPeds);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    CRF fitting functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%
%    fitCRF    %
%%%%%%%%%%%%%%%%
function [CRF CRFParams m] = fitCRF(s,m,CRF,doFit)

if nargin < 4, doFit = 1;end    
switch (m.CRFInterpFun)
 case {'loglinear'}
  % linear is just a polynomial fit of order 2
  [CRF CRFParams] = fitCRFPoly(s,m,CRF,2,m.CRFInterpFunParams{1},1);
 case {'linear'}
  % linear is just a polynomial fit of order 2
  [CRF CRFParams] = fitCRFPoly(s,m,CRF,2,m.CRFInterpFunParams{1});
 case {'poly'}
  % these parameters are passed to the poly fit. The second
  % specifies which CRFs have the same parameters (all or each) and
  % the first is the order of poly that is being fit. Offset
  % is always floating
  m.CRFInterpFunParams{end+1} = 4;
  % get initial parameters for the fit to the CRF. Here we are fitting
  % with a polynomial of order specified in the 1st parameter and
  % whether we are doing it across each or all CRFs as the 2nd param
  [CRF CRFParams] = fitCRFPoly(s,m,CRF,m.CRFInterpFunParams{1},m.CRFInterpFunParams{2});
 case {'nakarushton'}
  [CRF CRFParams] = fitCRFNakaRushton(s,m,CRF,m.CRFInterpFunParams{1},doFit);  
 case {'exponential'}
  [CRF CRFParams] = fitCRFExponential(s,m,CRF,m.CRFInterpFunParams{1},doFit);  
 otherwise
  disp(sprintf('(adp_selectionModel:fitCRFCDF) Unknown CRF interp type: %s',m.CRFInterpFun));
  keyboard
end

%%%%%%%%%%%%%%%%%%%%
%    fitCRFPoly    %
%%%%%%%%%%%%%%%%%%%%
function [CRF params] = fitCRFPoly(s,m,CRF,polyOrder,fitType,uselog)

% default
if ieNotDefined('fitType'),fitType = 'each';end
if ieNotDefined('uselog'),uselog = false;end
  
switch fitType
  %%%%%%%%%%%%%%% 
  % Fit all CRFs with same params except offset
  %%%%%%%%%%%%%%% 
 case {'all'}
  % number of CRFs we are fitting
  n = m.numConditions*2-1; %should not include cue4uncued, so subtract 1
  thisCRFnum = 1;
  A = [];
  R = [];
  for iCond = 1:m.numConditions
    for jCued = {'cued','uncued'}
      thisCRF = CRF.(m.conditions{iCond}).(jCued{1}).data;
      % create columns of polynomial matrix (i.e. first column is x, second column is x2...)
      % this is the same for all responses
      thisA = makePolyMatrix(thisCRF.c,polyOrder,uselog);
      % strip of ones column
      thisA = thisA(:,2:end);
      % add a zero column for each CRF's offset parameter
      thisA(:,polyOrder:polyOrder+n-1) = 0;
      % and set this CRF's offset column to 1
      thisA(:,polyOrder+thisCRFnum-1) = 1;
      % if we have the same CRF for both cued and uncued then, this
      % is not an independent CRF so assign the same number as in
      % the last (cued) CRF
      if ~(CRF.(m.conditions{iCond}).sameCuedUncued && strcmp(jCued{1},'uncued'))
        % concatenate together this matrix with all of the rest
        A = [A;thisA];
        % and add the responses for this CRF on to all of the responses
        R = [R;thisCRF.r'];
        thisCRFnum = thisCRFnum+1;
      end
    end
  end
  % now get least-squares fit params, by taking the pseudo inverse
  % of the A matrix computed above and multiplying by y (the responses).
  params = pinv(A)*R;
  % now go back throuch each condition, sort out which parameters belong to which CRF
  thisCRFnum = 0;
  for iCond = 1:m.numConditions
    for jCued = {'cued','uncued'}
      % update fit number as long as this is not a case in which it should
      % have the same fit as the cued condition
      if ~(CRF.(m.conditions{iCond}).sameCuedUncued && strcmp(jCued{1},'uncued'))
	thisCRFnum = thisCRFnum+1;
      end
      % get where in the parameter array the parameters for this function live. 
      % this is used in other parts of the code like getResponse which needs
      % to know where in the params array the params for that function should live
      CRF.(m.conditions{iCond}).(jCued{1}).fit.paramsIndex = [polyOrder+thisCRFnum-1 1:polyOrder-1];
    end
  end
  %%%%%%%%%%%%
  % Fit each CRF independently
  %%%%%%%%%%%%
 case {'each'}
  params = [];
  for iCond = 1:m.numConditions
    for jCued = {'cued','uncued'}
      if ~(CRF.(m.conditions{iCond}).sameCuedUncued && strcmp(jCued{1},'uncued'))
	thisCRF = CRF.(m.conditions{iCond}).(jCued{1}).data;
	% create columns of polynomial matrix (i.e. first column is 1, second column
	% is x, third column is x.^2, etc.
	A = makePolyMatrix(thisCRF.c,polyOrder,uselog);
	% now get least-squares fit params, by taking the pseudo inverse
	% of this matrix and multiplying by y (the responses).
	thisParams = pinv(A)*thisCRF.r';
	% compute the paramsIndex from the params
	CRF.(m.conditions{iCond}).(jCued{1}).fit.paramsIndex = [length(params)+1:length(params)+length(thisParams)];
	% keep params 
	params(end+1:end+length(thisParams)) = thisParams;
      else
	% this is sameCuedUncued, set the paramsIndex same as lat one, no need to compute an independent fit
	CRF.(m.conditions{iCond}).(jCued{1}).fit.paramsIndex = [length(params)+1-length(thisParams):length(params)];
      end
    end
  end
end

% fit a smooth version of the function with these parameters and return params
params = params(:)';
m.CRFParamsAreRaw = true;
CRF = computeSmoothFit(s,m,CRF,params);
m.CRFParamsAreRaw = false;

%%%%%%%%%%%%%%%%%%%%%%%%%%
%    computeSmoothFit    %
%%%%%%%%%%%%%%%%%%%%%%%%%%
function CRF = computeSmoothFit(s,m,CRF,params)

% if interp fun is fit then we are using a precomputed fit,
% no need to recompute smooth fit
if strcmp(m.CRFInterpFun,'fit'),return,end

% this function takes the params and computes a smooth fit of each CRF

% The contrast, arbitrairly chosen, over which to compute the fit
c = 0:0.01:1;

% cycle through CRFs
for iCond = 1:m.numConditions
  for iCued = {'cued','uncued'}
    if CRF.(m.conditions{iCond}).sameCuedUncued && strcmp(iCued,'uncued')
      %skip
    else
    % set the contrasts for the smooth fit
    CRF.(m.conditions{iCond}).(iCued{1}).fit.c = c;
    % get fit response over these contrasts
    CRF.(m.conditions{iCond}).(iCued{1}).fit.r = getResponse(s,m,params,CRF.(m.conditions{iCond}).(iCued{1}),c)';
    end
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    computeCRFFitResidual   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [fit] = computeCRFFitResidual(s,m,CRF,params)

fit.r = [];
fit.rste = [];
fit.fitR = [];
if m.emptyCRF,return,end

% if interp fun is fit then we are using a precomputed fit,
% no need to recompute smooth fit
if strcmp(m.CRFInterpFun,'fit'),return,end

% this function takes the params and computes a smooth fit of each CRF

% cycle through CRFs
for iCond = 1:m.numConditions
  for iCued = {'cued','uncued'}
    if CRF.(m.conditions{iCond}).sameCuedUncued && strcmp(iCued,'uncued')
      %skip
    else    % get the contrasts/responses of data
      c = CRF.(m.conditions{iCond}).(iCued{1}).data.c;
      fit.r(end+1,:) = CRF.(m.conditions{iCond}).(iCued{1}).data.r;
      fit.rste(end+1,:) = CRF.(m.conditions{iCond}).(iCued{1}).data.rste;
      % compute the responses by the fit
      fit.fitR(end+1,:) = getResponse(s,m,params,CRF.(m.conditions{iCond}).(iCued{1}),c)';
    end
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%
%    makePolyMatrix    %
%%%%%%%%%%%%%%%%%%%%%%%%
function A = makePolyMatrix(x,polyOrder,uselog)

if nargin < 3, uselog = false;,end

x = x(:);

if uselog, x = log(x);end

% this function makes a poly matrix for the passed in array x. That
% is the first column is 1, the second column is x, the third column
% is x^2 up until polyOrder. This is useful for fitting and computing
% polynomial functions
A = repmat(x,1,polyOrder).^repmat(0:polyOrder-1,size(x,1),1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    fitCRFNakaRushton    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [CRF params] = fitCRFNakaRushton(s,m,CRF,fitType,doFit)

% in some cases we just want to set up the CRF structure with the correct
% parameter indexes and not actually do the fit (because we have passed in param values)
if nargin < 5, doFit = 1;end

% get the CRF data, each row will be a different CRF fit
c = [];r = [];nCRF = 0;
for iCond = 1:m.numConditions
  for jCued = {'cued','uncued'}
    if CRF.(m.conditions{iCond}).sameCuedUncued && strcmp(jCued,'uncued') 
      % skip cue4uncued condition
      	CRF.(m.conditions{iCond}).(jCued{1}).fit.paramsIndex = CRF.(m.conditions{iCond}).cued.fit.paramsIndex;
    else
      thisCRF = CRF.(m.conditions{iCond}).(jCued{1}).data;
      % get all the contrast and response values
      c = [c;thisCRF.c];
      r = [r;thisCRF.r];
      % assign parameter index for each CRF, that is a number that
      % is used in decoding which set of parameters are used
      % to fit the CRF
      % first get which CRF this is
      nCRF = nCRF+1;
      % if we have the same CRF for both cued and uncued then, this
      % is not an independent CRF so assign the same number as in
      % the last (cued) CRF
% % % %       this was removed because the if sameCuedUncued statement takes care of this problem
% % % %       if CRF.(m.conditions{iCond}).sameCuedUncued && strcmp(jCued{1},'uncued')
% % % %         nCRF = nCRF-1;
% % % %       end
      % set the number
      CRF.(m.conditions{iCond}).(jCued{1}).parameterIndex = nCRF;
    end
  end
end


switch fitType
 case {'all'}
  % params for fitting
  maxiter = inf;
  % params are Rmax c50 n and offset for each CRF
  initParams = [1 0.07 3 zeros(1,nCRF)];
  minParams =  [0 0 -inf -inf(1,nCRF)];
  maxParams =  [inf inf inf inf(1,nCRF)];
 
  % YUKOnaka: make initParams from data. Rmax = mean(r@50%) c50 = mean(r@25%) n=3 offset = mean(r@12.5%)
  % fixed n-value
  offsetsStartHere = 3;
  if m.fixedN
    initN = []; 
    maxN=[]; 
    minN=[]; 
    offsetsStartHere = offsetsStartHere-1;
  else
    initN=2.1;
    maxN=3;
    minN=2;
  end
               %Rmax                    c50     n     offsets x 5
  initParams = [mean(r(1:nCRF,3))       0.25    initN r(1:nCRF,1)'];
  minParams =  [min(r(1:nCRF,3))-0.05   0       minN  r(1:nCRF,1)'-0.05];
  maxParams =  [max(r(1:nCRF,3))+0.05   1       maxN  r(1:nCRF,1)'+0.05];

  
%   optimParams = optimset('LevenbergMarquardt','on','MaxIter',maxiter); 
  optimParams = optimset('Algorithm','levenberg-marquardt','MaxIter',maxiter);
  % set the parameter indexes for each condition,
  % also save this as a structure so we can decode the
  % parameters that are passed in
  parameterIndex = [];
  for iCond = 1:m.numConditions
    for jCued = {'cued','uncued'}
      if CRF.(m.conditions{iCond}).sameCuedUncued && strcmp(jCued,'uncued') 
        % skip cue4uncued condition
      else
        % set the parameter indexes so getResponse knows how to get the parameters
        parameterIndex(end+1,:) = [1:offsetsStartHere CRF.(m.conditions{iCond}).(jCued{1}).parameterIndex+offsetsStartHere];
        CRF.(m.conditions{iCond}).(jCued{1}).fit.paramsIndex = parameterIndex(end,:);
      end
    end
  end
  
  % now go fit
  [params resnorm residual exitflag output lambda jacobian] = lsqnonlin(@nakaRushtonResidual,initParams,minParams,maxParams,optimParams,c,r,parameterIndex,m);
 case {'each'}
  % params for fitting
  maxiter = inf;
  % params are Rmax c50 n and offset for each CRF
  initParams = [ones(1,nCRF) ones(1,nCRF)*0.07 ones(1,nCRF)*3 zeros(1,nCRF)];
  minParams = [zeros(1,nCRF) zeros(1,nCRF) -inf(1,nCRF) -inf(1,nCRF)];
  maxParams = [inf(1,nCRF) ones(1,nCRF) inf(1,nCRF) inf(1,nCRF)];
  optimParams = optimset('LevenbergMarquardt','on','MaxIter',maxiter);
  optimParams.MaxFunEvals = 3000;
  
  % set the parameter indexes for each condition,
  % also save this as a structure so we can decode the
  % parameters that are passed in
  parameterIndex = [];
  for iCond = 1:m.numConditions
    for jCued = {'cued','uncued'}
      % set the parameter indexes so getResponse knows how to get the parameters
      parameterIndex(end+1,:) = CRF.(m.conditions{iCond}).(jCued{1}).parameterIndex+[1 2 3 4]*nCRF-nCRF;
      CRF.(m.conditions{iCond}).(jCued{1}).fit.paramsIndex = parameterIndex(end,:);
    end
  end
  % now go fit
  if doFit
    [params resnorm residual exitflag output lambda jacobian] = lsqnonlin(@nakaRushtonResidual,initParams,minParams,maxParams,optimParams,c,r,parameterIndex,m);
  else
    params = initParams;
  end
 otherwise
  disp(sprintf('(adp_selectionModel:fitCRFNakaRushoton) Unknown fit type: %s',fitType));
end

% fit a smooth version of the function with these parameters and return params
params = params(:)';
m.CRFParamsAreRaw = true;
CRF = computeSmoothFit(s,m,CRF,params);
m.CRFParamsAreRaw = false;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    nakaRushtonResidual    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function residual = nakaRushtonResidual(params,c,r,parameterIndex,m)

dispFit = 0;
 
% calculate residual for each crf
for i = 1:size(c,1)
  thisC = c(i,:);
  if ~m.fixedN
    % decode parameters
    Rmax = params(parameterIndex(i,1));
    c50 = params(parameterIndex(i,2));
    n = params(parameterIndex(i,3));
    offset = params(parameterIndex(i,4));
    % calculate naka-rushton, note that we use the log of c
    fitR(i,:) = nakaRushton(thisC,[Rmax c50 n offset]); 
  else %fixedN
    % decode parameters
    Rmax = params(parameterIndex(i,1));
    c50 = params(parameterIndex(i,2));
    n = m.fixedN;
    offset = params(parameterIndex(i,3));
    % calculate naka-rushton, note that we use the log of c
    fitR(i,:) = nakaRushton(thisC,[Rmax c50 n offset]); 
  end
end

% display fit if called for
if dispFit
  f = smartfig('selectionModel_nakaRushtonResidual','reuse');  
  clf;
  for i = 1:size(c,1)
    subplot(1,size(c,1),i)
    semilogx(c(i,:),r(i,:),'ko');
    hold on
    semilogx(c(i,:),fitR(i,:),'k-')
    if ~m.fixedN
      titleStr = sprintf('Rmax: %0.3f c50: %0.2f n: %0.3f\n',params(parameterIndex(i,1)),params(parameterIndex(i,2)),params(parameterIndex(i,3)));
      title(sprintf('%s offset: %f',titleStr,params(parameterIndex(i,4))));
    else
      titleStr = sprintf('Rmax: %0.3f c50: %0.2f n: %0.3f\n',params(parameterIndex(i,1)),params(parameterIndex(i,2)),m.fixedN);
      title(sprintf('%s offset: %f',titleStr,params(parameterIndex(i,3))));
    end
  end
  makeEqualYaxis(1,size(c,1));
  drawnow
end

residual = r-fitR;
residual = residual(:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    nakaRushtonTransformContrast    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function c = nakaRushtonTransformContrast(c,doInverse)

% function to transform contrast so that we can handle
% taking the log
if nargin == 1,doInverse = 0;end
if ~doInverse
  c = log10(c*100+10);
else
  c = ((10.^c)-10)/100;
end

%%%%%%%%%%%%%%%%%%%%%
%    nakaRushton    %
%%%%%%%%%%%%%%%%%%%%%
function response = nakaRushton(c,params)

% params are [Rmax c50 n offset]
c = nakaRushtonTransformContrast(c);
c50 = nakaRushtonTransformContrast(params(2));
response = params(1) * ((c.^params(3)) ./ ((c.^params(3)) + c50.^params(3))) + params(4);

%%%%%%%%%%%%%%%%%%%%%%%%%%
%    fitCRFExponential   %
%%%%%%%%%%%%%%%%%%%%%%%%%%
function [CRF params] = fitCRFExponential(s,m,CRF,fitType,dofit)

% in some cases we just want to set up the CRF structure with the correct
% parameter indexes and not actually do the fit (because we have passed in param values)
if nargin < 5, doFit = 1;end

% get the CRF data, each row will be a different CRF fit
c = [];r = [];nCRF = 0;
for iCond = 1:m.numConditions
  for jCued = {'cued','uncued'}
    thisCRF = CRF.(m.conditions{iCond}).(jCued{1}).data;
    % get all the contrast and response values
    c = [c;thisCRF.c];
    r = [r;thisCRF.r];
    % assign parameter index for each CRF, that is a number that
    % is used in decoding which set of parameters are used
    % to fit the CRF
    % first get which CRF this is
    nCRF = nCRF+1;
    % if we have the same CRF for both cued and uncued then, this
    % is not an independent CRF so assign the same number as in
    % the last (cued) CRF
    if CRF.(m.conditions{iCond}).sameCuedUncued && strcmp(jCued{1},'uncued')
      nCRF = nCRF-1;
    end
    % set the number
    CRF.(m.conditions{iCond}).(jCued{1}).parameterIndex = nCRF;
  end
end


switch fitType
 case {'all'}
  % params for fitting
  maxiter = inf;
  % params are scale, tau, offset
  initParams = [0.787 .573 0.3*ones(1,nCRF)];
  minParams = [-inf -inf -inf(1,nCRF)];
  maxParams = [inf inf inf(1,nCRF)];
  optimParams = optimset('LevenbergMarquardt','on','MaxIter',maxiter);
  % set the parameter indexes for each condition,
  % also save this as a structure so we can decode the
  % parameters that are passed in
  parameterIndex = [];
  for iCond = 1:m.numConditions
    for jCued = {'cued','uncued'}
      % set the parameter indexes so getResponse knows how to get the parameters
      parameterIndex(end+1,:) = [1:2 CRF.(m.conditions{iCond}).(jCued{1}).parameterIndex+2];
      CRF.(m.conditions{iCond}).(jCued{1}).fit.paramsIndex = parameterIndex(end,:);
    end
  end
  % now go fit
  [params resnorm residual exitflag output lambda jacobian] = lsqnonlin(@exponentialResidual,initParams,minParams,maxParams,optimParams,c,r,parameterIndex);
 case {'each'}
  % params for fitting
  maxiter = inf;
  % params are scale, center, tau, offset
  initParams = [0.787*ones(1,nCRF) 0.572*ones(1,nCRF) 0.3*ones(1,nCRF)];
  minParams = [-inf(1,nCRF) -inf(1,nCRF) -inf(1,nCRF)];
  maxParams = [inf(1,nCRF) inf(1,nCRF) inf(1,nCRF)];
  optimParams = optimset('LevenbergMarquardt','on','MaxIter',maxiter);
  % set the parameter indexes for each condition,
  % also save this as a structure so we can decode the
  % parameters that are passed in
  parameterIndex = [];
  for iCond = 1:m.numConditions
    for jCued = {'cued','uncued'}
      % set the parameter indexes so getResponse knows how to get the parameters
      parameterIndex(end+1,:) = CRF.(m.conditions{iCond}).(jCued{1}).parameterIndex+[1 2 3]*nCRF-nCRF;
      CRF.(m.conditions{iCond}).(jCued{1}).fit.paramsIndex = parameterIndex(end,:);
    end
  end
  % now go fit
  if doFit
    [params resnorm residual exitflag output lambda jacobian] = lsqnonlin(@exponentialResidual,initParams,minParams,maxParams,optimParams,c,r,parameterIndex);
  else
    params = initParams;
  end
 otherwise
  disp(sprintf('(adp_selectionModel:fitCRFExponential) Unknown fit type: %s',fitType));
end

% fit a smooth version of the function with these parameters and return params
params = params(:)';
m.CRFParamsAreRaw = true;
CRF = computeSmoothFit(s,m,CRF,params);
m.CRFParamsAreRaw = false;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    exponentialResidual    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function residual = exponentialResidual(params,c,r,parameterIndex)

dispFit = 0;

% calculate residual for each crf
for i = 1:size(c,1)
  thisC = c(i,:);
  % decode parameters
  scale = params(parameterIndex(i,1));
  tau = params(parameterIndex(i,2));
  offset = params(parameterIndex(i,3));
  % calculate exponential
  fitR(i,:) = exponential(thisC,[scale tau offset]); 
end

% display fit if called for
if dispFit
  f = smartfig('selectionModel_exponentialResidual','reuse');  
  clf;
  for i = 1:size(c,1)
    subplot(3,size(c,1),i)
    plot(c(i,:),r(i,:),'ko');
    hold on
    plot(c(i,:),fitR(i,:),'k-')
    titleStr = sprintf('scale: %0.3f tau: %0.3f\n',params(parameterIndex(i,1)),params(parameterIndex(i,2)));
    title(sprintf('%s offset: %f',titleStr,params(parameterIndex(i,3))));
    subplot(3,size(c,1),i+size(c,1))
    semilogx(c(i,:),r(i,:),'ko');
    hold on
    semilogx(c(i,:),fitR(i,:),'k-')
    title('semilog axis');
    subplot(3,size(c,1),i+2*size(c,1))
    plot(r(i,:)-fitR(i,:),'ko');
    title(sprintf('Residual: %s',mynum2str(r(i,:)-fitR(i,:),'sigfigs=2')));
  end
  makeEqualYaxis(3,size(c,1),1:size(c,1)*2);
  drawnow
end

%residual = log(r)-log(fitR);
residual = r-fitR;
residual = residual(:);

%%%%%%%%%%%%%%%%%%%%%
%    exponential    %
%%%%%%%%%%%%%%%%%%%%%
function response = exponential(c,params)

% params are [scale tau offset]
response = params(1) * (1-exp(-c/params(2)))+params(3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Parameter decoding
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%
%    getOutputNonlin    %
%%%%%%%%%%%%%%%%%%%%%%%%%
function outputNonlin = getOutputNonlin(m,params)

% if not a fit parameter, then just get from m
if (nargin==1) || ~m.outputNonlinFit
  outputNonlin = m.outputNonlin;
  % outputNonlin setting of 0 or 1 means non outputNonlin
  if any(outputNonlin==[0 1])
    outputNonlin = nan;
  end
else
  % get it from params
  outputNonlin = params(m.outputNonlinFitParamNum);
end


%%%%%%%%%%%%%%%%%%
%    getSigma    %
%%%%%%%%%%%%%%%%%%
function sigma = getSigma(m,params,condNum,cuedUncued)

% if condNum is -1 then we return a printable string
if condNum == -1
  cuedUncued = [1 0];
  for iCuedUncued = 1:2
    for i = 1:m.numConditions
      % get the sigma for each cue conditon and for cued or uncued
      sigma(i,iCuedUncued) = getSigma(m,params,i,cuedUncued(iCuedUncued));
    end
  end
  % format
  if ~m.sigmaCuedUncued && isempty(m.sigmaCuedUncuedRatio)
    sigma = mynum2str(sigma(:,1),'sigfigs=3','compact=1');
  else
    % have one sigma for cued and one for uncued
    sigmaStr = '';
    for i = 1:size(sigma,1)
      sigmaStr = sprintf('%s%s, ',sigmaStr,mynum2str(sigma(i,:),'sigfigs=4','compact=1'));
    end
    sigma = sigmaStr(1:end-2);
  end
  return
end

% get the sigma for the condition number. 
if m.sigmaForEachCond
  sigma = params(condNum);
elseif m.sigmaCuedUncued
  % first sigma is cued, second is uncued
  if cuedUncued
    sigma = params(1);
  else
    sigma = params(2);
  end
%  disp(sprintf('(adp_selectionModel:getSigma) sigma:%f condNum:%i cuedUncued:%i',sigma,condNum,cuedUncued));
else
  % Only one sigma, always first parameter
  sigma = params(1);
  % if we have given a fixed sigma ratio then apply it here
  if ~isempty(m.sigmaRatio)
    sigma = sigma*m.sigmaRatio(condNum);
  end
  % if we have given a fixed sigma ratio for cued versus uncued then apply it here
  if ~isempty(m.sigmaCuedUncuedRatio)
    % if the location is cued, then multiply by the cued uncued ratio
    if cuedUncued
      sigma = sigma*m.sigmaCuedUncuedRatio;
    end
  end
end


%%%%%%%%%%%%%%
%    getK    %
%%%%%%%%%%%%%%
function k = getK(m,params)

% if this is the sensitivity model then k doesn't matter, just return 1
if strcmp(m.poolingRule,'sensitivity')
  k = 1;
  return
end

if m.fixedK ~= 0
  k = m.fixedK;
  return
end

% get the pooling parameter
if m.sigmaForEachCond
  % sigma for each condition, k is after the sigmas
  k = params(m.numConditions+1);
elseif m.sigmaCuedUncued
  % sigma for cued and uncued, k is after 2 sigmas
  k = params(3);
else
  % single sigma k is second parameters
  k = params(2);
end

%%%%%%%%%%%%%%%%%%%%%%
%    getCRFParams    %
%%%%%%%%%%%%%%%%%%%%%%
function p = getCRFParams(m,params);

% just return the end of the params set (the first two are sigma and k).
if m.CRFParamsAreRaw
  p = params;
else
  % account for the sigma parameter(s)
  if m.sigmaForEachCond
    CRFParamsStartAt = m.numConditions+1;
  elseif m.sigmaCuedUncued
    CRFParamsStartAt = 3;
  else
    CRFParamsStartAt = 2;
  end
  % account for the k parameter
  if ~strcmp(m.poolingRule,'sensitivity')  && ~m.fixedK
    CRFParamsStartAt = CRFParamsStartAt+1;
  end
  p = params(CRFParamsStartAt:end);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Display functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%
%    dispSavedIter    %
%%%%%%%%%%%%%%%%%%%%%%%
function dispSavedIter(s,m)

% initialize counter
iterCount = 1;

% fit name
fitName = 'CRFCDF';
% get fitIterName
fitIterName{1} = sprintf('%s_iter%05i',fitName,iterCount);
while isComputed(fitIterName{end},s,m)
  % update counter
  iterCount = iterCount+1;
  % get next fitIterName
  fitIterName{end+1} = sprintf('%s_iter%05i',fitName,iterCount);
end

% no iterations found
if iterCount == 1
  disp(sprintf('(adp_selectionModel:dispSavedIter) No saved iterations found for %s',fitName));
  %dispHeader
  return
end

% display how many iterations found
maxIter = iterCount-1;
%dispHeader(sprintf('(adp_selectionModel:dispSavedIter) Found %i iterations of %s',maxIter,fitName));

% load the fitst state
iterNum = 1;
[s m CRF CDF params t resString] = loadState(sprintf('%s_iter%05i',fitName,iterNum),s,m);
% display it
dispModelThresholds(s,m,CRF,CDF,t,params,resString);

% some key codes
leftArrowKey = 124;
rightArrowKey = 125;
enterKey = 37;
escKey = 54;

% check if the state is saved
%dispHeader('Hit left or right arrow or ESC to exit');
while (1)
%  keyDown = first(find(mglGetKeys([leftArrowKey rightArrowKey enterKey escKey])));
%  if ~isempty(keyDown)
    % break on esc
%    if keyDown==4,break,end
    % increment or decrement iterNum
%    if any(keyDown==[2 3]),iterNum = min(iterNum+1,maxIter);end
%    if any(keyDown==[1]),iterNum = max(iterNum-1,1);end
iterNum  = min(iterNum+1,maxIter);
    disp(sprintf('(adp_selectionModel:dispSavedIter) Displaying iteration: %i',iterNum));
    % load the state
    [sold m CRF CDF params t resString] = loadState(sprintf('%s_iter%05i',fitName,iterNum),s,m);
    % display it
    dispModelThresholds(s,m,CRF,CDF,t,params,resString);
    drawnow;
%  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    dispModelThresholds    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dispModelThresholds(s,m,CRF,CDF,t,params,dispStr)

% don't do anything if noFig is srt
if isfield(s,'noFig') && s.noFig,return,end

if nargin < 6, dispStr = '';end

% get data thresholds
for iCond = 1:m.numConditions
  data(iCond,:) = CDF.(m.conditions{iCond}).data.t;
  c(iCond,:) = CDF.(m.conditions{iCond}).data.c;
  % get ste if available
  if isfield(CDF.(m.conditions{iCond}).data,'tste')
    dataste(iCond,:) = CDF.(m.conditions{iCond}).data.tste;
  else
    dataste(iCond,1:size(data,2)) = nan;
  end
end

% break up into two rows if the contrasts are non-monotonic
% this is just so the plot below will not plot lines between
% data taken with different conditions (like low and hi distracter contrasts)
nonMonotonic  = find(diff(c(1,:))<0);
if ~isempty(nonMonotonic)
  % make sure we are splitting at midpoint of matrix
  splitPoint = size(c,2)/2;
  if nonMonotonic == splitPoint
    % split data at this location
    c = [c(:,1:splitPoint) ; c(:,splitPoint+1:end)];
    data = [data(:,1:splitPoint) ; data(:,splitPoint+1:end)];
    dataste = [dataste(:,1:splitPoint) ; dataste(:,splitPoint+1:end)];
    t = [t(:,1:splitPoint) ; t(:,splitPoint+1:end)];
  end
end

dispFit = 0;
if dispFit
  % display fit to CDF
  f = smartfig('dispModelThresholds',s.reuse);
  a = gca(f);
  cla(a);h = [];
  for i = 1:size(c,1)
    thish = myloglog(c(i,:),data(i,:),'a',a,'Symbol','o','yError',dataste(i,:),'Color',getcolor(i+1),'MarkerFaceColor',getcolor(i+1),'MarkerEdgeColor=w');
    h(i) = thish{1};
    hold on
  end
  for i = 1:size(c,1)
    myloglog(c(i,:),t(i,:),'a',a,'Symbol','s-','Color',getcolor(i+1),'MarkerFaceColor','w','MarkerEdgeColor',getcolor(i+1));
  end
  legend(h,m.conditions,'Location','SouthEast');
  
  % display fit to CRF if we are also doing that
  if m.fitCRF
    h = title(sprintf('%s Parameters: sigma=%s k=%0.4f\nCRF:%s %s',s.saveName,getSigma(m,params,-1),getK(m,params),mynum2str(params(3:end),'compact=1','sigfigs=3'),dispStr));
    set(h,'Interpreter','none');
    % recompute smooth fit
    CRF = computeSmoothFit(s,m,CRF,params);
    % and display
    dispCRF(s,m,CRF,1,params);
  else
    dispCRF(s,m,CRF,1,params);
    title(sprintf('Parameters: sigma=%s k=%0.4f %s',getSigma(m,params,-1),getK(m,params),dispStr));
  end
  drawnow
end


%%%%%%%%%%%%%%%%%%
%    dispCond    %
%%%%%%%%%%%%%%%%%%
function dispCond(s,m,CDF)

% open figure
smartfig('adp_selectionModel:dispCond','reuse');clf;

for i = 1:m.numConditions
  for j = 1:m.numPeds(i)
    % set the plot number
    subplot(m.numConditions,max(m.numPeds),(i-1)*max(m.numPeds)+j);
    % get the colors to show bars in
    colorMap = repmat([0 0 0],length(CDF.(m.conditions{i}).data.d{j}),1);
    colorMap = [1 0.3 0.4;colorMap];
    % display the bar graph
    mybar(100*[CDF.(m.conditions{i}).data.c(j) CDF.(m.conditions{i}).data.d{j}],'colorMap','hsv','drawAxis=0');
    % place a title
    if j == round(m.numPeds(i)/2),title(m.conditions{i});end
  end
  % set y labels
  subplot(m.numConditions,max(m.numPeds),(i-1)*max(m.numPeds)+1);
  ylabel('Contrast (%)');
end
makeEqualYaxis(m.numConditions,max(m.numPeds));

%%%%%%%%%%%%%%%%%
%    dispCRF    %
%%%%%%%%%%%%%%%%%
function dispCRF(s,m,CRF,plotSemilog,params)

% don't do anything if noFig is srt
if isfield(s,'noFig') && s.noFig,return,end

if nargin < 4, plotSemilog = 1; end
if nargin < 5, params = [];end

% open figure
smartfig('adp_selectionModel:dispCRF','reuse');clf;

if plotSemilog,numRows = 2;else numRows = 1;end
if ~any(m.outputNonlin == [0 1])
  plotOutputNonlin = true;
  numRows = 3;
else
  plotOutputNonlin = false;
end

for rowNum = 1:numRows
  % if we are on the third row, it means to plot the outputNonlin
  if rowNum ==3
    doOutputNonlin = true;
  else
    doOutputNonlin = false;
  end
  % plot the curves
  for i = 1:m.numConditions
    subplot(numRows,m.numConditions+1,i+(m.numConditions+1)*(rowNum-1));
    c1 = getSmoothColor(i*2-1,2*m.numConditions,'hsv');
    plotCRF(CRF.(m.conditions{i}).cued,c1,doOutputNonlin,m,params);
    
    c2 = getSmoothColor(i*2,2*m.numConditions,'hsv');
    if ~CRF.(m.conditions{i}).sameCuedUncued %if these are not replicated points
      plotCRF(CRF.(m.conditions{i}).uncued,c2,doOutputNonlin,m,params);
    end
    % label
    xlabel('Contrast');
    ylabel('Response');
    titleStr = m.conditions{i};
    
    %YUKO: title fixed. paramsIndex must be offset if fixedK
    if m.fixedK
      accountForFixedK = 1;
    else
      accountForFixedK = 0;
    end    
	% 'each' : [sigma k     offsetA slopeA offsetB slopeB offsetC slopeC offsetD slopeD offsetE slopeE ...]
	% 'all'  : [sigma k     slopeALL offsetA offsetB offsetC offsetD offsetE ...]
    CRFParamsStartHere = 2;

    if ~isempty(params) && (isfield(m,'fitCRF') && m.fitCRF) && isfield(CRF.(m.conditions{i}).cued,'fit') && isfield(CRF.(m.conditions{i}).cued.fit,'paramsIndex')
      if ~isfield(CRF.(m.conditions{i}).uncued,'fit') %sometimes, uncued will not have a fit value if it was not calculated. (cue4uncued)
        CRF.(m.conditions{i}).uncued.fit = CRF.(m.conditions{i}).cued.fit; %just make it equal to cued for display purposes
      end
      titleStr = sprintf('%s\nc:%s u:%s',titleStr,mynum2str(params(CRF.(m.conditions{i}).cued.fit.paramsIndex+CRFParamsStartHere-accountForFixedK)),mynum2str(params(CRF.(m.conditions{i}).uncued.fit.paramsIndex+CRFParamsStartHere-accountForFixedK)));
    end
    title(titleStr);
    %mylegend({'cued','uncued'},{{'.' c1} {'.' c2}},2);
    if rowNum==2,set(gca,'XScale','log');end
  end

  subplot(numRows,m.numConditions+1,m.numConditions+1+(m.numConditions+1)*(rowNum-1));
  for i = 1:m.numConditions
    c1 = getSmoothColor(i*2-1,2*m.numConditions,'hsv');
    plotCRF(CRF.(m.conditions{i}).cued,c1,doOutputNonlin,m,params);
    
    c2 = getSmoothColor(i*2,2*m.numConditions,'hsv');
    if ~CRF.(m.conditions{i}).sameCuedUncued %if these are not replicated points
      plotCRF(CRF.(m.conditions{i}).uncued,c2,doOutputNonlin,m,params);
    end
  end
  if rowNum==2,set(gca,'XScale','log');end
  % label
  xlabel('Contrast');
  ylabel('Response');
end

% equal y-axis
makeEqualYaxis(numRows,m.numConditions+1);
drawnow;

%%%%%%%%%%%%%%%%%%
%    plotCRF     %
%%%%%%%%%%%%%%%%%%
function plotCRF(CRF,c,doOutputNonlin,m,params)

% only pass in m if you want to apply output nonlinearity
if nargin < 3
  doOutputNonlin = false;
end

% plot the data
if isfield(CRF,'data')
  r = CRF.data.r;
  if doOutputNonlin
    r = applyOutputNonlin(m,r,params);
  end
  plot(CRF.data.c,r,'.','Color',c);
  hold on
end  
% plot the fit
if isfield(CRF,'fit')
  r = CRF.fit.r;
  if doOutputNonlin
    r = applyOutputNonlin(m,r,params);
  end
  plot(CRF.fit.c,r,'-','Color',c);
  hold on
end

%%%%%%%%%%%%%%%%%
%    dispCDF    %
%%%%%%%%%%%%%%%%%
function dispCDF(s,m,CDF)

% open figure
smartfig('adp_selectionModel:dispCDF','reuse');clf;

for i = 1:m.numConditions
  c = getSmoothColor(i,m.numConditions,'hsv');
  plotCDF(CDF.(m.conditions{i}),c);
end
% label axis
xlabel('Contrast');
ylabel('Threshold');

%%%%%%%%%%%%%%%%%%
%    plotCDF     %
%%%%%%%%%%%%%%%%%%
function plotCDF(CDF,c)

% plot the data
if isfield(CDF,'data')
  myloglog(CDF.data.c,CDF.data.t,'Symbol','.','Color',c);
  hold on
end  
% plot the fit
if isfield(CDF,'fit')
  myloglog(CDF.fit.c,CDF.fit.t,'Symbol','-','Color',c);
  hold on
end

%%%%%%%%%%%%%%%%%%%%%%%%%%
%    dispResponseDist    %
%%%%%%%%%%%%%%%%%%%%%%%%%%
function dispResponseDist(s,m,params,rdist,pooled)

% display response distributions
smartfig('adp_selectionModel:dispResponseDist','reuse');clf;
nLocs = size(rdist,2);
% get number of rows
if ieNotDefined('pooled'),nRows = nLocs;else nRows = nLocs+1;end

for i = 1:nLocs
  subplot(nRows,1,i);
  myhist(rdist{2,i},40,'r');
  myhist(rdist{1,i},40,'k');
  title(sprintf('Loc %i: %0.3f %0.3f',i,mean(rdist{1,i}),mean(rdist{2,i})));
  xlabel('Response');ylabel('n');
end

% now draw pooled
if ~ieNotDefined('pooled')
  subplot(nRows,1,nRows);
  myhist(pooled{2},40,'r');
  myhist(pooled{1},40,'k');
  title(sprintf('Pooled: k=%0.4f area under ROC=%0.4f',getK(m,params),computeAreaUnderROC(pooled)));
  xlabel('Response');ylabel('n');
end  

makeEqualXaxis(nRows,1);
drawnow;

%%%%%%%%%%%%%%%%%%%%%%%%%%
%    displModelFitCDF    %
%%%%%%%%%%%%%%%%%%%%%%%%%%
function dispModelFitCDF(s,m,params,CDF,fitThreshold)

smartfig('adp_selectionModel:dispModelFit','reuse');clf;
for iCond = 1:m.numConditions
  % get condition name
  condName = m.conditions{iCond};
  % and color
  c = getSmoothColor(iCond,m.numConditions,'hsv');
  % and plot data
  myloglog(CDF.(condName).data.c,CDF.(condName).data.t,'Symbol','.','Color',c);
  hold on;
  % plot model
  myloglog(CDF.(condName).data.c,fitThreshold(iCond,:),'Symbol','-','Color',c);
end
xlabel('Pedestal Contrast');
ylabel('Threshold');
title(sprintf('sigma: %0.4f k: %0.4f',getSigma(m,params,1),getK(m,params)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Interpolation functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

x = x(:);
y = y(:);
% fit the line
A = ones(size(x,1),2);
A(:,1) = x;
b = pinv(A)*y;

A = ones(size(xi,1),2);
A(:,1) = xi;
yi = A*b;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    pooling rules
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%
%   softmax   %%
%%%%%%%%%%%%%%%%%
function p = softmax(q,tau)

% default tau setting
if nargin == 1,tau = 1;end

% make sure to pass in a *column* array e.g. [3 4 1 5]'

% softmax operates column-wise, that is, it computes the softmax operator 
% across each column -- converting each element into a weight such that
% the weights across each column sum to 1. With tau set small (say 0.1), this operates
% like max weighting (i.e. the entry in the row with the largest value approach 1 and
% all others approach 0). With tau set large (say 10), this acts like average weighting
% where each element in a row approaches the same value.
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


%%%%%%%%%%%%%%%%%%%
%    calcStats    %
%%%%%%%%%%%%%%%%%%%
function stats = calcStats(s,m,CRF,CDF,t,fitParams)

% get data thresholds and ste
for iCond = 1:m.numConditions
  data(iCond,:) = CDF.(m.conditions{iCond}).data.t;
  dataste(iCond,:) = CDF.(m.conditions{iCond}).data.tste;
end

% check data size
if ~isequal(size(data),size(t))
  disp(sprintf('(adp_selectionModel:calcStats) Size of data and model thresholds are mismatched'));
  keyboard
end

% get the CRF fit residual
CRFFit = computeCRFFitResidual(s,m,CRF,fitParams);
CRFss = sum((CRFFit.r(:) - CRFFit.fitR(:)).^2);

% calculate sum-of-squares difference between data and model
CDFss = sum((data(:)-t(:)).^2);

% calculate chi-squared (similar to above, but normalized by ste)
CDFchi = sum(((data(:)-t(:))./dataste(:)).^2);
CRFchi = sum(((CRFFit.r(:)-CRFFit.fitR(:))./CRFFit.rste(:)).^2);

% k is number of parameters
k = length(fitParams);

% n is number of data points
CDFn = length(data(:));
CRFn = length(CRFFit.r(:));

% get combined stats
if m.fitCRF
  ss = CDFss+CRFss;
  chi = CDFchi+CRFchi;
  n = CDFn+CRFn;
else
  ss = CDFss;
  chi = CDFchi;
  n = CDFn;
end

% calc AIC and BIC
BIC = n * log (ss/n) + k*log(n);
AIC = 2*k + n * log (ss/n) + 2*k*(k+1)/(n-k-1);

% calculate degrees of freedom
v = n-k;
% calculate chi statistic, this number should be big
q = gammainc(0.5*v,0.5*chi);

% display
disp(sprintf('(adp_selectionModel:calcStats) Params: %s',mlrnum2str(fitParams,'sigfigs=4')));
disp(sprintf('(adp_selectionModel:calcStats) AIC: %0.3f BIC: %0.3f (ss=%0.3f k=%i n=%i)',AIC,BIC,ss,k,n))

% calculate variance accounted for
stats.CDFr2 = 1-sum((log(data(:))-log(t(:))).^2)/sum((log(data(:))-mean(log(data(:)))).^2);
stats.CRFr2 = 1-sum((CRFFit.r(:)-CRFFit.fitR(:)).^2)/sum((CRFFit.r(:)-mean(CRFFit.r(:))).^2);
disp(sprintf('(adp_selectionModel:calcStats) CDF r^2: %s CRF r^2: %s',mynum2str(stats.CDFr2,'sigfigs=2','compact=1'),mynum2str(stats.CRFr2,'sigfigs=2','compact=1')));

% package up to return
stats.k = k;
stats.n = n;
stats.ss = ss;
stats.AIC = AIC;
stats.BIC = BIC;

% compute chi-squared for CDF and CRF


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Save / load functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%
%%   saveState   %%
%%%%%%%%%%%%%%%%%%%
function saveState(fitType,s,m,CRF,CDF,params,t,aux)

% some other variables that we can save optionally
if nargin < 7, t = [];end
if nargin < 8, aux = [];end

% remove randomization distributions
m = rmfield(m,'r');
m = rmfield(m,'randDistractorsDist');

% get filename
filename = getSaveName(fitType,s,m);

% display what we are doing
disp(sprintf('(adp_selectionModel) Saving state: %s',filename));

% save the file
save(filename,'s','m','CRF','CDF','params','t','aux');

%%%%%%%%%%%%%%%%%%%
%%   loadState   %%
%%%%%%%%%%%%%%%%%%%
function [s m CRF CDF params t aux] = loadState(fitType,s,m)

% auxillary variables that may not be saved
t = [];aux = [];

% get save filename
filename = getSaveName(fitType,s,m);

% load or give error
if isfile(filename)
  disp(sprintf('(adp_selectionModel) Loading state: %s',filename));
  load(filename);
else
  disp(sprintf('(adp_selectionModel:loadState) Could not find file: %s',filename));
  keyboard
end

% make fixed random distributions
m = makeFixedRandDistributions(m);
m = checkModelFields(m);

%%%%%%%%%%%%%%%%%%%%
%%   isComputed   %%
%%%%%%%%%%%%%%%%%%%%
function tf = isComputed(fitType,s,m)

filename = getSaveName(fitType,s,m);
tf = isfile(filename);

%%%%%%%%%%%%%%%%%%%%%
%%   getSaveName   %%
%%%%%%%%%%%%%%%%%%%%%
function filename = getSaveName(fitType,s,m)

% get fit type
%fitType = '';
%if isfield(s.todo,'fitCRFCDF'),fitType = 'CRFCDF';end
%if isfield(s.todo,'fitCRF'),fitType = 'CRF';end
%if isfield(s.todo,'fitCDF'),fitType = 'CDF';end

% get CRF interp function
CRFInterpFun = '';
if isfield(m,'CRFInterpFun')
  CRFInterpFun = m.CRFInterpFun;
end

% get each or all
eachOrAll = '';
if isfield(m,'CRFInterpFunParams')
  for i = 1:length(m.CRFInterpFunParams)
    % look for each/all
    if ~isempty(find(strcmp(m.CRFInterpFunParams{i},{'any','each'})))
      eachOrAll = m.CRFInterpFunParams{i};
    end
  end
end

% get pooling rule
poolingRule = m.poolingRule;

% cross-validation
crossVal = '';
if s.crossValidate,crossVal = 'xval';end

% sigma for each condition
sigmaForEachCond = '';
if m.sigmaForEachCond
  sigmaForEachCond = 'sigmaForEachCond';
end

% get sigma ratio
sigmaRatio = '';
if ~isempty(m.sigmaRatio)
  sigmaRatio = mynum2str(m.sigmaRatio,'sigfigs=1','compact=1');
end

% sigma cued uncued
sigmaCuedUncued = '';
if m.sigmaCuedUncued
  sigmaCuedUncued = 'sigmaCuedUncued';
end

% get sigma cued uncuedratio
sigmaCuedUncuedRatio = '';
if ~isempty(m.sigmaCuedUncuedRatio)
  sigmaCuedUncuedRatio = mynum2str(m.sigmaCuedUncuedRatio,'sigfigs=1','compact=1');
end

% non zero or 1 outputNonlin
outputNonlin='';
if ~any(m.outputNonlin == [0 1])
  outputNonlin = sprintf('outnonlin%s',mynum2str(m.outputNonlin,'sigfigs=2','compact=1'));
end
if m.outputNonlinFit
  outputNonlin = 'outnonlinFit';
end

% non 2 reslog
reslog = '';
if s.reslog ~= 2
  reslog = sprintf('reslog%s',mynum2str(s.reslog,'sigfigs=2','compact=1'));
end

%YUKO added fixedKStr and weightedResStr and fixedN
% fixedKStr
if m.fixedK
  fixedKStr = sprintf('fixedK%s',fixBadChars(num2str(m.fixedK),{'.' 'p'}));
else
  fixedKStr = [];
end

% weightedResStr
weightedResStr = [];
if m.CDFres~=1
  weightedResStr = [weightedResStr sprintf('CDFres%i',m.CDFres)];
end
if m.CRFres~=1
  weightedResStr = [weightedResStr sprintf('CRFres%i',m.CRFres)];
end
if m.normResByVar==1
  weightedResStr = [weightedResStr sprintf('normResByVar')];
end

if m.fixedN
  fixedNStr = sprintf('fixedN%s',fixBadChars(num2str(m.fixedN),{'.' 'p'}));
else
  fixedNStr = [];
end
  

% combine together variables to get filename - if the variable is not
% empty then add an _ before it to distinguish it from other parts of the filename
filenameVars = {'fitType','poolingRule','CRFInterpFun','eachOrAll','sigmaRatio','crossVal','sigmaForEachCond','reslog','fixedKStr','fixedNStr','weightedResStr','sigmaCuedUncued','sigmaCuedUncuedRatio','outputNonlin'};
if ~isempty(s.saveName)
  filename = sprintf('%s_%s',s.saveName,fitType);
else
  filename = sprintf('%s',fitType);
end  
for i = 2:length(filenameVars)
  if ~isempty(eval(filenameVars{i}))
    filename = sprintf('%s_%s',filename,eval(filenameVars{i}));
  end
end
filename = fixBadChars(filename,{{'.','p'},{' ','_'}},[],255);

if ~isempty(m.filename) %YUKO: if filename is specificed by user
  m.filename = stripext(m.filename);
  filename = m.filename;
end

% add on path and extension
filename = setext(fullfile(s.savePath,filename),'mat');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Obsolete - for reference only
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%
%    fitCDF    %
%%%%%%%%%%%%%%%%
function fitCDF_old(s,m,CRF,CDF)

% set parameter to show that we are *not* fitting the CRF
m.fitCRF = false;

% set parameter info
%initParams = [0.005 20];
%initParams = [0.011278 1.946632];
initParams = [0.022069 10];
minParams = [0 1];
maxParams = [inf inf];

% for this, we are using an already precomputed CRF
m.CRFInterpFun = 'fit';

% use Levenberg-Marquardt to find function minimum
%[fitParams resnorm residual exitflag output lambda jacobian] = lsqnonlin(@computeModelResidualUsingPercentCorrect,initParams,minParams,maxParams,s.optimParams,s,m,CRF,CDF);

% use Nelder-Mead method
[fitParams fval exitflag output] = fminsearch(@computeModelResidualUsingThresholds,initParams,s.optimParams,s,m,CRF,CDF);

% get the model thresholds
t = computeModelThreshold(fitParams,s,m,CRF,CDF);

% display model fit.
dispModelFitCDF(s,m,fitParams,CDF,t);

% calculate statistics
stats = calcStats(s,m,CRF,CDF,t);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    computeModelResidualUsingPrecentCorrect    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function is not used any more - it computes an error based on
% percent correct which does not do as good a job of fitting
% as fitting based on log thresholds
function residual = computeModelResidualUsingPercentCorrect(params,s,m,CRF,CDF)

% This function computes the model residual (i.e. the difference between the
% expected percent correct for d'=1 (0.7605) and what the model gives with
% the current settings

% cycle over experimental conditions (e.g. focal, distributed)
for iCond = 1:m.numConditions
  % get condition name and relevant info
  condName = m.conditions{iCond};
  thisCRF = CRF.(condName);
  thisCDF = CDF.(condName).data;
  nPed = length(thisCDF.c);
  % cycle over pedestal contrasts
  for iPed = 1:nPed
    % get target and distractor contrasts for this pedestal
    distractorContrasts = thisCDF.d{iPed};
    nDistractors = length(distractorContrasts);
    targetContrast = thisCDF.c(iPed);
    % get all the contrasts for this pedestal
    contrasts = [targetContrast distractorContrasts];
    % now put the necessary number of cued CRF functions into CRFs
    for iTarg = 1:thisCDF.nCued
      CRFs{iTarg} = thisCRF.cued;
    end
    % and fill the rest with the uncued CRFs
    for iTarg = thisCDF.nCued+1:nDistractors+1
      CRFs{iTarg} = thisCRF.uncued;
    end
    % now get the probability correct form the model
    p(iCond,iPed) = simTrial(s,m,params,iCond,contrasts,CRFs,thisCDF.t(iPed));
  end
end

residual = p(:)-0.7605;
residual = sqrt(sum(residual.^2));

% display current parameter settings
if s.verbose>0
  disp(sprintf('(adp_selectionModel:computeModelResidual) Parameters: %s (r2=%0.6f)',mynum2str(params),sqrt(sum(residual.^2))));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Helper Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%
%    getRangeOfCRF    %
%%%%%%%%%%%%%%%%%%%%%%%
function range = getRangeOfCRF(s,m,CRF)

% this code gets the range of response values across conditions. Note
% that it takes the lowest value as the lowest value measured. The highest
% value is the highest value that would have been obtained for a 100% contrast
range(1) = inf;
range(2) = -inf;
measuredResponseMax = -inf;
cuedUncued = {'cued','uncued'};
for iCond = 1:m.numConditions
  for jCued = length(cuedUncued)
    if isfield(CRF.(m.conditions{iCond}).(cuedUncued{jCued}).data,'r')
      % if the response fit is there, then set the range to the min and max values found
      range(1) = min(range(1),min(CRF.(m.conditions{iCond}).(cuedUncued{jCued}).data.r));
      measuredResponseMax = max(measuredResponseMax,max(CRF.(m.conditions{iCond}).(cuedUncued{jCued}).data.r));
    end
    % if there is a fit, then set the range to go out to a maximum for the 
    % inferred response to the highest contrast (1)
    if isfield(CRF.(m.conditions{iCond}).(cuedUncued{jCued}).fit,'r')
      range(2) = max(range(2),max(CRF.(m.conditions{iCond}).(cuedUncued{jCued}).fit.r));
    end
  end
end

% display what we have found
if s.verbose
  disp(sprintf('(adp_selectionModel:getRangeOfCRF) Output nonlin range: [%0.3f %0.3f]. Measured responses minmax: [%0.3f %0.3f] Max response inferred by fit for full contrast: %0.3f',range(1),range(2),range(1),measuredResponseMax,range(2)));
end



%%%%%%%%%%%%%%
%    logb    %
%%%%%%%%%%%%%%
function r = logb(x,b)

r = log(x)./log(b);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    checkDataSturctures    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [tf m CRF] = checkDataStructures(CRF,CDF,s,m)

tf = 1;

if isempty(CDF)
  disp(sprintf('(adp_selectionModel:CDF) Empty CDF'));
  tf = 0;
  return;
end

% first get the experiment name
m.conditions = fieldnames(CDF);
m.numConditions = length(m.conditions);

% display message
if s.verbose>0
  %dispHeader('Checking data structures')
  for i = 1:m.numConditions
    disp(sprintf('Found condition: %s',m.conditions{i}));
  end
end

% Check on the CDF that there are fields that specify the experiment type
for i = 1:m.numConditions
  if ~isfield(CDF.(m.conditions{i}),'data')
    disp(sprintf('(adp_selectionModel:checkDataStructures) Missing field CDF.%s.data - mixed up CRF/CDF arguments?',m.conditions{i}));
    tf = 0;
    return
  end
  [tf missingFieldName] = checkForFields(CDF.(m.conditions{i}).data,{'d','nCued'});
  if ~tf
    disp(sprintf('(adp_selectionModel:checkDataStructures) Missing field CDF.%s.data.%s',m.conditions{i},missingFieldName));
    return
  end
  % make sure that the number of pedestal contrasts in CDF match the number of distractors
  if length(CDF.(m.conditions{i}).data.d) ~= length(CDF.(m.conditions{i}).data.c)
    disp(sprintf('(adp_selectionModel:checkDataStructures) In CDF.%s.data, the number of pedestals: %i must match the number of distractor arrays: %i',m.conditions{i},length(CDF.(m.conditions{i}).data.d),length(CDF.(m.conditions{i}).data.c)));
  end
  % get the number of pedestals tested
  m.numPeds(i) = length(CDF.(m.conditions{i}).data.d);
end

% check for empty CRF. No error if empty, since that means that
% we are fitting with an uncosntrained CRF
if isempty(CRF)
  disp(sprintf('(adp_selectionModel:CRF) Empty CRF, building structure'));
  m.emptyCRF = true;
  for iCond = 1:m.numConditions
    
    % set the cued and uncued data to empty
    CRF.(m.conditions{iCond}).cued.data.c   = [];
    CRF.(m.conditions{iCond}).cued.data.r   = [];
    CRF.(m.conditions{iCond}).uncued.data.c = [];
    CRF.(m.conditions{iCond}).uncued.data.r = [];

    % if 4 locations are cued then that means all
    % locations are cued, so set sameCuedUncued
    if CDF.(m.conditions{iCond}).data.nCued   == 4
      CRF.(m.conditions{iCond}).sameCuedUncued = 1;
    else
      CRF.(m.conditions{iCond}).sameCuedUncued = 0;
    end
  end
end

% check that there is a matching field on CRF;
[tf missingFieldName] = checkForFields(CRF,m.conditions);
if ~tf
  disp(sprintf('(adp_selectionModel:checkDataStructures) Missing field %s if CRF',missingFieldName));
  return
end

% now for each condition name, go check that there is a cued and non-cued CRF
for i = 1:m.numConditions
  [tf missingFieldName] = checkForFields(CRF.(m.conditions{i}),{'cued','uncued','sameCuedUncued'});
  if ~tf
    disp(sprintf('(adp_selectionModel:checkDataStructures) Missing field CRF.%s.%s',m.conditions{i},missingFieldName));
    return
  end
end

% check for sigma ratio to match correctly
if ~isempty(m.sigmaRatio)
  if length(m.sigmaRatio) ~= m.numConditions
    disp(sprintf('(adp_selectionModel:checkDataStructures) Must have a sigmaRatio number for each condition'));
    keyboard
  end
end

%dispHeader;

%%%%%%%%%%%%%%%%%%%%%%%%
%    checkForFields    %
%%%%%%%%%%%%%%%%%%%%%%%%
function [tf missingFieldName]  = checkForFields(s,fieldsToCheck)

% defaults
tf = 1;
missingFieldName = [];

% if missing any field, spit out error and return
for i = 1:length(fieldsToCheck)
  if ~isfield(s,fieldsToCheck{i})
    missingFieldName = fieldsToCheck{i};
    tf = 0;
    return
  end
end

%%%%%%%%%%%%%%
%    todo    %
%%%%%%%%%%%%%%
function tf = todo(s,todoName)

tf = 0;
if isfield(s.todo,todoName) && s.todo.(todoName)
  tf = 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    makeFixedrandDistributions    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function m = makeFixedRandDistributions(m)

% older versions did not have this field
if ~isfield(m,'nStimuli'),m.nStimuli = 4;end

% set the random distributions for the distractors if necessary
m.randDistractorsDist = [];
if m.randDistractors
  m.randDistractorsDist = ceil(rand(m.nStimuli-1,m.n)*(m.nStimuli-1))+1;
end

% create fixed random distributions
% one for each interval
for i = 1:2
  % one for each location
  for j = 1:m.nStimuli
    m.r{i,j} = randgauss(0,1,m.n);
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%
%    setSystemParams    %
%%%%%%%%%%%%%%%%%%%%%%%%%
function [s args] = setSystemParams(args)

% returns s for the system parameters based on input arguments to function
% s should contain high level parameters, like which model is being run
% whether to display figures and verbose output etc.

% parse input arguments
verbose = [];dispFig = [];dispCRF = [];dispCDF = [];dispCond = [];fitCDF = [];dispIter = [];fitCRFCDF = [];recompute = [];saveName = [];crossValidate = [];saveIter = [];dispSavedIter = [];reslog = [];fixedParams = [];initParams = [];reuseFig = [];noFig = 0; savePathDir=[];
[argNames argValues args] = getArgs(args,{'verbose=0',...
  'dispFig=0',  'dispCRF=0',  'dispCDF=1',...
  'fitCDF=1',   'dispIter=0', 'dispCond=0',...
  'fitCRFCDF=0','recompute=0','saveName=[]',...
  'crossValidate=0','saveIter=0','dispSavedIter=1',...
  'reslog=2','fixedParams=[]','initParams=[]',...
  'reuseFig=1','noFig=0','savePathDir=[]'});

% make a todo list
s.todo = [];
if dispSavedIter, s.todo.dispSavedIter = true;end
if dispCRF,       s.todo.dispCRF       = true;end
if dispCDF,       s.todo.dispCDF       = true;end
if dispCond,      s.todo.dispCond      = true;end
if fitCDF,        s.todo.fitCDF        = true;end
if fitCRFCDF,     s.todo.fitCRFCDF     = true;end

% Whether to crossvalidate or not
s.crossValidate = crossValidate;

% set in structure some variables that get run
s.verbose = verbose;

s.dispIter = dispIter;

% dispFig is for showing more detailed figs or not
s.dispFig = dispFig;
% no fig is for displaying no figs at all.
s.noFig = noFig;

% maximum number of iterations of nonlinear serach 
s.maxIter = inf;

% optimization parameters
s.optimParams = optimset('MaxIter',s.maxIter);

% create save path
if isempty(savePathDir)
  s.savePath = '/home/frk/Dropbox/adaptation/adp_selectionModel';
else %%YUKO added to allow users to save fitStates in different dirs (useful for separating insideScanner/outsideScanner fits)
  s.savePath = savePathDir;
end
if ~isdir(s.savePath),mkdir(s.savePath);end
s.saveName = saveName;
fprintf('[%s] SaveName: %s\nn',mfilename,saveName)
fprintf('[%s] SaveDir : %s\nn',mfilename,savePathDir)

% whether to load data or recompute
s.recompute = recompute;

% if we were passed in fixed params (in which case we just use the params instead of fitting)
s.fixedParams = fixedParams;

% if we want to use specific initParams
s.initParams = initParams;

% whether to save every iteration of the optimization (useful for debugging the fits)
s.saveIter = saveIter;
global selectionModelSaveIter;selectionModelSaveIter = 1;

% save the log base to use for residuals
s.reslog = reslog;

% whether to reuse figures
if reuseFig, s.reuse = 'reuse';else s.reuse = '';end

if s.verbose >= 0
  %%dispHeader('System parameters');
  disp(sprintf('Recompute: %i',s.recompute));
  disp(sprintf('Verbose: %i',s.verbose));
  disp(sprintf('Cross-validate: %i',s.crossValidate));
  disp(sprintf('Disp iterations: %i',s.dispIter));
  disp(sprintf('Save iter: %i',s.saveIter));
  disp(sprintf('reslog: %0.2f',s.reslog));
  disp(sprintf('fixedParams: %s',mynum2str(s.fixedParams,'sigfigs=3','compact=1')));
  disp(sprintf('initParams: %s',mynum2str(s.initParams,'sigfigs=3','compact=1')));
end

%%%%%%%%%%%%%%%%%%%%%%%%
%    setModelParams    %
%%%%%%%%%%%%%%%%%%%%%%%%
function m = setModelParams(s,args)

m = [];
% parse arguments
poolingRule = [];sigmaRatio = [];sigmaForEachCond=[];eachOrAll = [];fixedK = [];fixedN = [];CRFInterpFun = [];randDistractors = [];CDFres=[];CRFres=[];normResByVar=[]; filename=[];
getArgs(args,{'poolingRule=exponential','sigmaRatio=[]','sigmaForEachCond=0', ...
              'eachOrAll=each','fixedK=0','fixedN=0','CRFInterpFun=nakarushton', ...
              'randDistractors=0','CDFres=1','CRFres=1','normResByVar=0','filename=[]', ...
              'sigmaCuedUncuedRatio=[]','sigmaCuedUncued=0','outputNonlin=0','outputNonlinFit=0'});

% number of repeats 
m.n = 100000;

% number of stimuli
m.nStimuli = 4;

% default is to use fit
m.CRFInterpFun = 'fit';

% just a flag for getCRFParams to know
% when params arrays have other params in them or not
m.CRFParamsAreRaw = false;

% set pooling rule
m.poolingRule = poolingRule;

% if we have a single sigma and want to make a fixed ratio of sigma
% between different conditions (for instance to model affect of noise
% reduction with attention, then this should be set as an array
% with each element corresponding to each condition: e.g. [1 1.5] would
% give 50% more noise to the second condition
m.sigmaRatio = sigmaRatio;

% If we want to have separate sigmas for each condition (i.e. let them float
% in the model like we do for the sensitivity model set this to true
m.sigmaForEachCond = sigmaForEachCond;

% This is to have the model set different sigmas for the cued and
% uncued stimuli. Note that this does not set different sigmas for
% different cueing conditons (i.e. cue-one, cue-two and cue-four, which
% is set above with sigmaRatio or sigmaForEachCond). Since these
% are mutually exclusive items, we check for conflicts
m.sigmaCuedUncuedRatio = sigmaCuedUncuedRatio;
m.sigmaCuedUncued = sigmaCuedUncued;

% check for conflicting conditions for sigma
if ~isempty(m.sigmaRatio) && (m.sigmaForEachCond==1)
  disp(sprintf('(adp_selectionModel) Either set sigmaRatio (which sets the ratio between different sigmas for each condition) or sigmaEachCond (which allows the sigma to float for each condition)'));
  keyboard
end
if ~isempty(m.sigmaCuedUncuedRatio) && (m.sigmaCuedUncued==1)
  disp(sprintf('(adp_selectionModel) Either set sigmaCuedUncuedRatio (which sets the ratio between different sigmas for cued and uncued targets) or sigmaCuedUncued (which allows the sigma to float between cued and uncued)'));
  keyboard
end
if (~isempty(m.sigmaCuedUncuedRatio) || (m.sigmaCuedUncued==1)) && (~isempty(m.sigmaRatio) || (m.sigmaForEachCond==1))
  disp(sprintf('(adp_selectionModel) Either allow sigma to change between cued/uncued locations (sigmaCuedUncued/sigmaCuedUncuedRatio) or sigma to change between different cueing conditons (sigmaRatio/sigmaForEachCond)'));
  keyboard
end

% if we want to fix the k-vale
m.fixedK = fixedK;

% for this, we are going to use a polynomial fit to the data
m.CRFInterpFun = CRFInterpFun;

% if we want to fix the n-value in nakarushton
m.fixedN = fixedN;
if ~strcmp(m.CRFInterpFun,'nakarushton')
  m.fixedN = 0;
end

% sets whether to fit all CRFs with the same parameters (except offset) or with different
% paramaters for each
m.CRFInterpFunParams = {eachOrAll};

% If the distractor contrasts were chosen at random with replacement from a set of distractor contrasts
% on each trial. 
m.randDistractors = randDistractors;

% create fixed random distributions
m = makeFixedRandDistributions(m);

% default to assuming that there is CRF, when the data structure gets checked
% later this can get set to true in which case CDFCRF fits are not constrained
% by any data
m.emptyCRF = false;

% weight residuals as 1 unless there is an input
m.CDFres = CDFres;
m.CRFres = CRFres;
m.normResByVar = normResByVar;
if m.normResByVar
  m.CDFres=1;
  m.CRFres=1;
end

% set the output nonlinearity assumed between firing rates and BOLD responses
% Note that what we do is after reading off the BOLD response we pass it
% through an exponent set to whatever value this is set to. Thus if ooutputNonlin
% is set to 1, this assumes a linear function between firing rates and BOLD.
% However, if we assume something like a threshold non-linearity between firing-rate
% and BOLD (see Schummers, Yu and Sur), then to get back something like firing rate
% from BOLD responses we need an expansive non-linearity, so this value should be
% set to < 1.
m.outputNonlin = outputNonlin;
% The range parameter is important. It sets what range to consider the range
% of 0-1. That is the power law fits ans exponent to the range 0-1 and in some
% papers (like Devor) for instance they have normalized the BOLD response to
% a min of 0 max of 1 - which depends on your experiment (the smaller the
% range that your experiment covers, the more likely to be linear. Here
% we start by setting the range to empty. This will make the code later
% set the range between the minimum and maximum of the values wer measured.
m.outputNonlinRange = [];
% outputNonlin can be a fit parameter
m.outputNonlinFit = outputNonlinFit;

% if user specified a file to load
m.filename = stripext(filename);

if s.verbose >= 0
  %%dispHeader('Model parameters');
  disp(sprintf('Pooling rule: %s',m.poolingRule));
  disp(sprintf('Sigma ratio: %s',mynum2str(m.sigmaRatio,'sigfigs=2','compact=1')));
  disp(sprintf('Sigma for each condition: %i',m.sigmaForEachCond));
  disp(sprintf('Sigma cued/uncued ratio: %s',mynum2str(m.sigmaCuedUncuedRatio,'sigfigs=2','compact=1')));
  disp(sprintf('Sigma cued/uncued: %i',m.sigmaCuedUncued));
  disp(sprintf('Fixed k: %0.2f',m.fixedK));
  disp(sprintf('Fixed n: %0.2f',m.fixedN));
  disp(sprintf('CRF fit: %s',eachOrAll));
  disp(sprintf('randDistractors: %i',randDistractors));
  disp(sprintf('output nonlinearity: %f',m.outputNonlin))
  disp(sprintf('output nonlinearity fit: %i',m.outputNonlinFit))
end

%%%%%%%%%%%%%%%%%%%%%%%%%%
%    checkModelFields    %
%%%%%%%%%%%%%%%%%%%%%%%%%%
function m = checkModelFields(m)

% add fields that might be missing (since the loaded simulation was done before
% these fields were added to structure).
if ~isfield(m,'sigmaCuedUncued')
  m.sigmaCuedUncued = false;
end
if ~isfield(m,'sigmaCuedUncuedRatio')
  m.sigmaCuedUncuedRatio = [];
end

% fix a typo in parameter
if isfield(m,'ouptutNonlinFitParamNum')
  m.outputNonlinFitParamNum = m.ouptutNonlinFitParamNum;
  m = rmfield(m,'ouptutNonlinFitParamNum');
end

  









