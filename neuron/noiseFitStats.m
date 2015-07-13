function p = noiseFitStats(varargin)
% function p = noiseFitStats(varargin)
%
% this function computes some statistics on the Noise fit as returned by
% the the noise comparison (i.e., by the function
% makenoisestatisticsfigs2_r2_compareNoise.m)
%
% the data file used is:
% /Volumes/data/riken/crfaa/fmridata/crfaa/03-Aug-2009CompareNoiseFit_nkn2.mat
% /Volumes/data/riken/crfaa/fmridata/crfaa/30-Jul-2009CompareNoiseFit_sg.mat
%
% results are not saved. But a screen print out is given.
%
% NB the program canbe called with 'fitType' set to 'nk2' to do the test on
%    the naka-rushton fit, and with 'sg' to do the test on the skewed
%    gaussian fit.
%
% franco pestilli 2009.07.31

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
observer  =[];
filename  =[];
adaptation=[];
fitType   =[];     % load the naka-rushton fit or the skewed gaussian fit
savePlots =[];     % 1 saves the plot in the current dierectory
dispFit   =[];     % if 1 the plot is displayed
numBoots  =[];     % number of bootstrap samples
saveFig   =[];     % if 1 saves a figure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set default variables:
getArgs(varargin,{...
 'observer',{'fm' 'fp' 'jg'}, ...  
 'adaptation',[0 28 100], ...
 'dispFit',0, ...
 'saveFig',1, ...
 'fitType', 'nk2', ... % 'sg' or 'nk2'
 'filename',[date,'allNoise.mat']});

defaultDataFolder = '/Volumes/data/riken/crfaa/fmridata/crfaa';

% load the proper file.
[sigma offset r2] = loadData(defaultDataFolder,adaptation(1),fitType);

meanR2 = mean(squeeze(r2),1);
steR2  = std(squeeze(r2));

% make one t-test for each visual area:
for v = 1:size(sigma,2)
 thisSigma = squeeze(sigma(:,v,:));
 diffSigma = diff(thisSigma,[],2);
 [H0_sigma(v) p_sigma(v)] = ttest(diffSigma);
 
 thisOffset = squeeze(offset(:,v,:));
 diffOffset = diff(thisOffset,[],2);
 [H0_offset(v) p_offset(v)] = ttest(diffOffset);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (A) compute the average sigma across observers %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Mv1234_sigma = squeeze(mean(squeeze(sigma),2));


% make one test across visual areas and observers:
% sigma:
thisSigma_a = squeeze(sigma(2:end,:,1));
thisSigma_a = thisSigma_a(:);
thisSigma_d = squeeze(sigma(2:end,:,2));
thisSigma_d = thisSigma_d(:);
diffSigma_all = diff([thisSigma_d,thisSigma_a],[],2);
[H0_sigmaAll p_sigmaAll] = ttest(diffSigma_all);
disp(sprintf('Mean Sigma Distributed %1.3f.',mean(thisSigma_a)))
disp(sprintf('ste  Sigma Distributed %1.3f.',std(thisSigma_a)/sqrt(3)))
disp(sprintf('Mean Sigma Attended    %1.3f.',mean(thisSigma_d)))
disp(sprintf('ste  Sigma Attended    %1.3f.',std(thisSigma_d)/sqrt(3)))
disp('%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%')


% offset
thisOffset_a = squeeze(offset(2:end,:,1));
thisOffset_a = thisOffset_a(:);
thisOffset_d = squeeze(offset(2:end,:,2));
thisOffset_d = thisOffset_d(:);
disp(sprintf('Mean Offset Distributed %1.3f.',mean(thisOffset_a)))
disp(sprintf('ste  Offset Distributed %1.3f.',std(thisOffset_a)/sqrt(3)))
disp(sprintf('Mean Offset Attended    %1.3f.',mean(thisOffset_d)))
disp(sprintf('ste  Offset Attended    %1.3f.',std(thisOffset_d)/sqrt(3)))
disp('%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%')

diffOffset_all = diff([thisOffset_d,thisOffset_a],[],2);
[H0_offsetAll p_offsetAll] = ttest(diffOffset_all);


disp('%% <T-Test Across Visual areas and observers> %%')
disp(sprintf('%%    Offset p = %1.12f    %%',p_offsetAll))
disp(sprintf('%%    Sigma  p = %1.12f    %%',p_sigmaAll))
disp('%% <T-Test Across Visual areas and observers> %%')


% nested anova test:
doAnova = 0;
if doAnova
 s_factors = {{[sigma(:,1,1)]' [sigma(:,1,2)]'} ...
  {[sigma(:,2,1)]' [sigma(:,2,2)]'} ...
  {[sigma(:,3,1)]' [sigma(:,3,2)]'} ...
  {[sigma(:,4,1)]' [sigma(:,4,2)]'}};
 nestedanova(s_factors);
 
 offset_factors = {{[offset(:,1,1)]' [offset(:,1,2)]'} ...
  {[offset(:,2,1)]' [offset(:,2,2)]'} ...
  {[offset(:,3,1)]' [offset(:,3,2)]'} ...
  {[offset(:,4,1)]' [offset(:,4,2)]'}};
 nestedanova(offset_factors);
end

keyboard

% plot a bar graph of the r2 values:


%%%%%%%%%%%%
% loadData %
%%%%%%%%%%%%
function [sigma  offset  r2 filename] = loadData(defaultDataFolder,adaptation,fitType)

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

% load the file
[filename vars2load] = makeFileName(defaultDataFolder,fitType);
data = load(sprintf('%s',filename),vars2load{1});


% collect all noise estimates in on array
for o = 1:size(data.noise,2) % for each observer
 for v = 1:size(data.noise,3) % for each visual area
  r2(v,o,1) = data.noise{1,o,v}.bestFit.r2; % these are the r2 of Dnoise-2-Dcrf
  r2(v,o,2) = data.noise{2,o,v}.bestFit.r2; % these are the r2 of Dnoise-2-Acrf
  r2(v,o,3) = data.noise{3,o,v}.bestFit.r2; % these are the r2 of Dnoise-2-Acrf

  sigma(v,o,1) = data.noise{1,o,v}.k;
  sigma(v,o,2) = data.noise{3,o,v}.k;
  
  offset(v,o,1) = data.noise{1,o,v}.responseOffset;
  offset(v,o,2) = data.noise{3,o,v}.responseOffset;
 end
end


%%%%%%%%%%%%%%%%
% makeFileName %
%%%%%%%%%%%%%%%%
function [filename vars2load] = makeFileName(defaultDataFolder,fitType)

switch fitType
 case {'nk2' 'nakarushton'}
filename = [defaultDataFolder,'/03-Aug-2009CompareNoiseFit_nkn2.mat'];
 case {'sg' 'skewedgaussin'}
filename = [defaultDataFolder,'/30-Jul-2009CompareNoiseFit_sg.mat'];
otherwise
keyboard
end
vars2load = {'noise'};



