function p = tvcParamsStats(varargin)
% function p = tvcParamsStats(varargin)
%
% thi sfunction computes some statistics on the TvC fit as returned by
% the the noise comparison (i.e., by the function
% makenoisestatisticsfigs2_r2_compareNoise.m)
% 
% the data file used is:
% /Volumes/data/riken/crfaa/fmridata/crfaa/31-Jul-2009CompareNoiseFit_nkn2.mat
% 
% results are saved in:
% defaultDataFolder = '/Users/frakkopesto/data/riken/crfaa/fmridata/crfaa/';
%
% and figures are saved in:
% figDir = 'fig_secondDip';
% defaultDataFolder = '/Users/frakkopesto/data/riken/crfaa/fmridata/crfaa/';
%
% these is the way the TvC parameters are coded:
% rmax   = 1; - fixed to '1';
% n      = params(1);
% c50    = params(2);
% p      = params(3);
% dr     = params(4); % differentiation factor
%
% - the CRF has this form in the code:
% r = rmax.*((c.^n)./((c.^(n*p)) + (c50.^(n*p))));
% 
% - but this form in the text:
% r = rmax.*((c.^p+q)./((c.^(q)) + (c50.^(q))));
% 
% so  the real p as in the paper is: n - (n*p)
% and the real q is: n*p
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
 'observer',{ 'average' 'fm' 'fp' 'jg'}, ...
 'adaptation',[0 28 100], ...
 'dispFit',0, ...
 'saveFig',1, ...
 'numBoots', 10, ...
 'filename',[date,'allNoise.mat']});

defaultDataFolder = '/Volumes/data/riken/crfaa/fmridata/crfaa';

% NB % for the moment the datat files were created only for adapt-0
iAdapt = 1; 

for iObserver = 1:length(observer)
 disp(sprintf('[tvcParamsStats] OBS: %s',observer{iObserver}))
 % load the proper file.
 tvcfit = loadData(observer{iObserver},'v3',defaultDataFolder,adaptation(iAdapt));
 for a = 1:2 % run throught each cue condition
  if ~isempty(tvcfit)
   disp(sprintf('[tvcParamsStats] Second dip Rmax <%s> Threshold <%s> SlopeP <%s> and SlopeQ <%s>.', ...
    num2str(tvcfit(a, 1)),num2str(tvcfit(a, 2)), ...
    num2str(tvcfit(a, 3)),num2str(tvcfit(a, 4))));
   
   p.tvcRmax(iObserver,a) = ones(size(tvcfit(a, 2)));
   p.tvcThreshold(iObserver,a) = tvcfit(a, 2);
   p.tvcP(iObserver,a) = tvcfit(a, 1) - (tvcfit(a, 1)*tvcfit(a, 3));
   p.tvcQ(iObserver,a) = tvcfit(a, 1).*tvcfit(a, 3);
   p.tvcDr(iObserver,a) = tvcfit(a, 4);
  
   % all params
   p.tvcParams(iObserver,a,1:5) = [p.tvcRmax(iObserver,a) p.tvcThreshold(iObserver,a) ...
                                   p.tvcP(iObserver,a) p.tvcQ(iObserver,a) p.tvcDr(iObserver,a)]; % extract the second dip parameters

  end
 end
end

% remove 0's from the matrices
p.tvcDr(p.tvcRmax == 0) = nan;
p.tvcThreshold(p.tvcThreshold == 0) = nan;
p.tvcP(p.tvcP == 0) = nan;
p.tvcQ(p.tvcQ == 0) = nan;
p.tvcDr(p.tvcRmax == 0) = nan;

% make average across observers
dr_mean = nanmean(p.tvcDr,1);
dr_ste = nanstd(p.tvcDr,1)./sqrt(length(observer)-1);

threshold_mean = log10(nanmean(10.^p.tvcThreshold,1));
threshold_ste = nanstd(10.^p.tvcThreshold,1)./sqrt(length(observer)-1);

rmax_mean = nanmean(p.tvcRmax,1);
rmax_ste = nanstd(p.tvcRmax,1)./sqrt(length(observer)-1);

tvcP_mean = nanmean(p.tvcP,1);
tvcP_ste = nanstd(p.tvcThreshold,1)./sqrt(length(observer)-1);

tvcQ_mean = nanmean(p.tvcQ,1);
tvcQ_ste = nanstd(p.tvcQ,1)./sqrt(length(observer)-1);

allparams = [threshold_mean;rmax_mean;tvcP_mean;tvcQ_mean;dr_mean];
allparams_ste = [threshold_ste;rmax_ste;tvcP_ste;tvcQ_ste;dr_ste];

% make a bra graph fo the average parameters:
plotBarTvCparams(allparams,allparams_ste,saveFig)


% make a t-test for significance:
[h_dr,p_dr] = ttest([10.^p.tvcDr(:,1)-10.^p.tvcDr(:,2)])
[h_thr,p_thr] = ttest([10.^p.tvcThreshold(:,1)-10.^p.tvcThreshold(:,2)])
[h_p,p_p] = ttest([p.tvcP(:,1)-p.tvcP(:,2)])
[h_q,p_q] = ttest([p.tvcQ(:,1)-p.tvcQ(:,2)])
keyboard



%%%%%%%%%%%%%%%%%%%%
% plotBarTvCparams %
%%%%%%%%%%%%%%%%%%%%
function plotBarTvCparams(allparams,allparams_ste,saveFig)

figurename = 'tvcparams';
h = smartfig(figurename,'reuse');
subplot(2,2,1)
mybar(allparams(1,:), ...
 'yError',allparams_ste(1,:), ...
 'groupLabels',{'Threshold'}, ...
 'withinGroupLabels',{'Distributed cue, target' 'Focal cue, target'}, ...
 'yAxisMin',0, ...
 'yAxisMax',.04, ...
 'yLabelText',{'Contrast (%)'},'dispValues',1);

subplot(2,2,2)
mybar(allparams(3,:), ...
 'yError',allparams_ste(3,:), ...
 'groupLabels',{'p'}, ...
 'withinGroupLabels',{'Distributed cue, target' 'Focal cue, target'}, ...
 'yAxisMin',0, ...
 'yAxisMax',.5, ...
 'yLabelText',{'-na-'},'dispValues',1);

subplot(2,2,3)
mybar(allparams(4,:), ...
 'yError',allparams_ste(4,:), ...
 'groupLabels',{'q'}, ...
 'withinGroupLabels',{'Distributed cue, target' 'Focal cue, target'}, ...
 'yAxisMin',0, ...
 'yAxisMax',3, ...
 'yLabelText',{'-na-'},'dispValues',1);

subplot(2,2,4)
mybar(allparams(5,:), ...
 'yError',allparams_ste(5,:), ...
 'groupLabels',{'Dr'}, ...
 'withinGroupLabels',{'Distributed cue, target' 'Focal cue, target'}, ...
 'yAxisMin',0, ...
 'yAxisMax',.25, ...
 'yLabelText',{sprintf('fMRI response\n(%s signal change)','%')},'dispValues',1);

if saveFig
 set(h,'PaperPosition',[.25 .25 8 10.5]);
 set(h,'PaperOrientation','Portrait');
 
 disp(sprintf('Saving figure %s',figurename));
 figDir = 'fig_tvc_params';
 defaultDataFolder = '/Users/frakkopesto/data/riken/crfaa/fmridata/crfaa/';
 figDir = savefig(h,'figName',figurename,'defaultDataFolder', defaultDataFolder,'figDir',figDir);
end

%%%%%%%%%%%%
% loadData %
%%%%%%%%%%%%
function [tvcfit filename] = loadData(observer,visualAreas,defaultDataFolder,adaptation)

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
% this is from compareNoise the fucntio who made the files being loaded
% 'observer',{'avrg','jg' 'fm' 'fp'}, ...
if strcmpi(observer,'fm')
 o = 3;
elseif  strcmpi(observer,'fp')
 o = 4;
elseif strcmpi(observer,'jg')
 o = 2;
elseif strcmpi(observer,'average')
 o = 1; 
end

if ~(strcmpi('fm',observer) && adaptation == 100)
 [filename vars2load] = makeFileName(observer,adaptation,visualAreas,defaultDataFolder);
 data = load(sprintf('%s',filename),vars2load{1});
 tvcfit(1,1:6) = data.results{o,adp,1}.behavior.tvcfit.fitParams.fit.fitParams;
 tvcfit(2,1:6) = data.results{o,adp,2}.behavior.tvcfit.fitParams.fit.fitParams;
 
else
 tvcfit = [];
end


%%%%%%%%%%%%%%%%
% makeFileName %
%%%%%%%%%%%%%%%%
function [filename vars2load] = makeFileName(observer,adaptation,visualArea,defaultDataFolder)  
filename = [defaultDataFolder,'/31-Jul-2009CompareNoiseFit_nkn2.mat'];
vars2load = {'results'};
