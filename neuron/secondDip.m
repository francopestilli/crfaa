function p = secondDip(varargin)
% function p = secondDip(varargin)
%
% thi sfunction computes some statistics on the second dip as returned by
% the TvC fits done during the the noise comparison (i.e., by the function
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
 'numBoots', 10000, ...
 'filename',[date,'allNoise.mat']});

defaultDataFolder = '/Volumes/data/riken/crfaa/fmridata/crfaa';
if isdir(defaultDataFolder)
 
else
 defaultDataFolder = '/Volumes/homosacer_data/riken/crfaa/fmridata/crfaa/data_used_files';
end

for iAdapt = 1:length(adaptation)
 disp(sprintf('[secondDip] ADP: %i',adaptation(iAdapt))) 
 for iObserver = 1:length(observer)
  disp(sprintf('[secondDip] OBS: %s',observer{iObserver}))
  % load the proper file.
  tvcfit = loadData(observer{iObserver},'v3',defaultDataFolder,adaptation(iAdapt));
  for iAttend = 1:2
   if ~isempty(tvcfit{1})
    disp(sprintf('[secondDip] Second dip Threshold <%s> and slope <%s>.', ...
     num2str(tvcfit{iAttend}.fitParams(5)),num2str(tvcfit{iAttend}.fitParams(6))));
    
    p.dipParams(iObserver,iAdapt,iAttend,1:2) = tvcfit{iAttend}.fitParams(5:6); % extract the second dip parameters
    p.dipThreshold(iObserver,iAdapt,iAttend) = p.dipParams(iObserver,iAdapt,iAttend,1).*tvcfit{iAttend}.fitParams(2);
    p.dipXc50(iObserver,iAdapt,iAttend) = p.dipParams(iObserver,iAdapt,iAttend,1);
    p.dipSlope(iObserver,iAdapt,iAttend) = p.dipParams(iObserver,iAdapt,iAttend,2);
   end
  end
 end
 p.dipSlope(p.dipSlope == 0) = nan;
 p.dipThreshold(p.dipThreshold == 0) = nan; 
end
keyboard
% make average across observers
slope_mean = squeeze(nanmean(p.dipSlope,1));
slope_ste = squeeze(nanstd(p.dipSlope,1))./sqrt(length(observer));

threshold_mean = squeeze(10.^(nanmean(log10(p.dipThreshold),1)));
threshold_ste = squeeze(nanstd(log10(p.dipThreshold),1))./sqrt(length(observer));

c50factor_mean = squeeze(nanmean(p.dipXc50,1));
c50factor_ste = squeeze(nanstd(p.dipXc50,1))./sqrt(length(observer));

cue = {'Distributed' 'Focal      '};

disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('%%                        Adapt-0   %%   Adapt-28   %%  Adapt-100 %%')
for c  = 1:length(cue)
 disp(sprintf('<%s> Mean Slope: %s - %s - %1.3f ',cue{c},slope_mean(:,c)))
end
disp('%%                        Adapt-0   %%   Adapt-28   %%  Adapt-100 %%')
for c  = 1:length(cue)
 disp(sprintf('<%s> Ste  Slope: %s - %s - %1.3f ',cue{c},slope_ste(:,c)))
end

disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('%%                            Adapt-0   %%   Adapt-28   %%  Adapt-100 %%')
for c  = 1:length(cue)
 disp(sprintf('<%s> Mean Threshold: %s - %s - %1.3f ',cue{c},threshold_mean(:,c)))
end
disp('%%                            Adapt-0   %%   Adapt-28   %%  Adapt-100 %%')
for c  = 1:length(cue)
 disp(sprintf('<%s> Ste  Threshold: %s - %s - %1.3f ',cue{c},threshold_ste(:,c)))
end
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')


keyboard

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
 tvcfit{1} = data.results{o,adp,1}.behavior.tvcfit.fitParams.fit;
 tvcfit{2} = data.results{o,adp,2}.behavior.tvcfit.fitParams.fit;
 
else
 tvcfit{1} = [];
end



%%%%%%%%%%%%%%%%
% makeFileName %
%%%%%%%%%%%%%%%%
function [filename vars2load] = makeFileName(observer,adaptation,visualArea,defaultDataFolder)  
filename = [defaultDataFolder,'/31-Jul-2009CompareNoiseFit_nkn2.mat'];
if isfile(filename)
else
filename = [defaultDataFolder,'/03-Sep-2009CompareNoiseFit_nkn2.mat']; 
end
vars2load = {'results'};
