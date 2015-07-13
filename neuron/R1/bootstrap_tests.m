function [p d] = bootstrap_tests(varargin)
% function p = bootstrap_tests(varargin)
%
% computes the statistical tests 
% using the results of the bootstrap fits 
% for both the sensory and the selection model
%
%
% franco pestilli 2010/09/02

% initialize input vars
model = [];
dataPath = [];
getArgs(varargin,{ ...
        'model', {'sensory' 'selection'}, ...
        'dataPath', '/data2/crfaa/crfaa/data_used_files_R1/'});

d = cell(1,2); p = d;
% (A) run the tests
for t = 1:length(model)
 switch model{t}
  case {'sensory', 'se', 'early'}
   [p{1} d{1}] = doSensory(dataPath);

  case {'selection' 'sm', 'late'}
   [p{2}  d{2}] = doSelection(dataPath);
  
  otherwise
   keyboard
 end
end


%%%%%%%%%%%%%
% doSensory %
%%%%%%%%%%%%%
function [p d] = doSensory(dataPath)
% (1) load the data 
% sensoryNoise dimensions = (iBoot,iFit,iObserver, iArea)
dataFile = '2010-50-AugCompareNoiseFit_Dnt_bootNum100_nkn2_070_BOOT.mat';

d = load(fullfile(dataPath,dataFile),'sensoryNoise','responseOffset');

% sigma
 p.sensory = computeProbability(d.sensoryNoise,'sensory');
 
% OFFSET
 p.baseline = computeProbability(d.responseOffset,'baseline');


%%%%%%%%%%%%%%%
% doSelection %
%%%%%%%%%%%%%%%
function [p d] = doSelection(dataPath)
% (1) load the data 

d1.sensitivity = nan(100,3,3,4);
d1.selection  = d1.sensitivity;
for bt = 1:100
 disp(sprintf('[%s] Computing selection model sigma for bootstrap#%i.',mfilename,bt))                 
 sigma = crfaaMaxPool(sprintf('dataDir=%scrfaaMaxPoolBootstrap100',dataPath), ...
                      'bootstrapNum',bt,'doSigma=1','dispFig=0');
                     
 d1.sensitivity(bt,:,:,:) = sigma.sensitivity;
 d1.selection(bt,:,:,:)   = sigma.selection;
end

% compute the p-values
[p.sensory     d.sensory]  = computeProbability(d1.sensitivity,'sensory');
[p.selection   d.selection]= computeProbability(d1.selection,'selection');


%%%%%%%%%%%%%%%%%%%%%%
% computeProbability %
%%%%%%%%%%%%%%%%%%%%%%
function [p d] = computeProbability(value,title)
% this function computes the p-value 
% from the results fo the bootstraps
% can be used for sigma and for the baseline either 
% returned by makenoisestatistis_boot
% or
% by crfaamaxpool.m
count = 0;
n_diff = nan(size(value,1),size(value,3),size(value,4));
p = ones(size(value,3),size(value,4));
for o = 1:size(value,3)
 for va = 1:size(value,4)
  count = count + 1;
  n_focal(:,o,va) = squeeze(value(:,3,o,va));
  n_distributed(:,o,va) = squeeze(value(:,1,o,va));
  switch title
   case {'baseline'}
    n_diff(:,o,va) = n_focal(:,o,va) - n_distributed(:,o,va);
    H0 = 0;
   case {'sensory' 'selection' }
    n_diff(:,o,va) = n_distributed(:,o,va) ./ n_focal(:,o,va);
    H0 = 1;
   otherwise
    keyboard
  end
  
  % computing probability of H1
  if min(n_diff(:,o,va)) <= H0
   temp1 = sort(n_diff(:,o,va));
   temp2 = temp1 <= H0;
   p(o,va) = max(find(temp2))/length(temp1);
  else
   p(o,va) = 1/size(value,1);
  end
  % plot a figure
  if count == 1
   smartfig(title,'reuse');
  end
  subplot(size(value,3),size(value,4),count);
  hist(n_diff(:,o,va)); 
  if ~ieNotDefined('temp2') 
   hold on, vline(temp1(max(find(temp2))));
  end
  drawnow
 end
end

d.ratio       = n_diff;
d.focal       = n_focal;
d.distributed = n_distributed;

