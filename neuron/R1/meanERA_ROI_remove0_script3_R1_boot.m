function [d d1 d2 ts events] = meanERA_ROI_remove0_script3_R1_boot(groupNum,scanNum,r2thr,saveTag,target_roiFileNames,recomputedDeconv, session)
% function [d d1 d2 ts events] = meanERA_ROI_remove0_script3_R1_boot(groupNum,scanNum,r2thr,saveTag,target_roiFileNames, recomputedDeconv, session)
%
% this file calls a bunch of functions to do several types of projection
% out
% it has a cutomized analysis for each subject and adaptation condition
%
%
% fp 2008.07.29

tic

% check if exist a mrSession.mat
if ~(exist('mrSession.mat','file') == 2)
 disp(sprintf('[meanERA_ROI_remove0_script3_R1_boot] mrSession.mat not found in this directory. please cd to a mrLoadret dir'));
 return
end

if ieNotDefined('saveTag')
 saveTag = num2str(100*r2thr);
end

target_roiFileNames = switchVisualArea(target_roiFileNames);
saveTag = ['2010-20-July','_',target_roiFileNames{1}(1:2),'_',saveTag];
stdGroups = {1:32}; % one grop for each hrf
tr =.8;

% choose whethe to recompute the deconvolution or load a precomputed file.
if recomputedDeconv
 disp(sprintf('[meanERA_ROI_remove0_script3_R1_boot] DO Deconv ANALYSIS ON group %i, scan %i roi: %s',groupNum,scanNum,target_roiFileNames{1}(1:2)));
 
 % this is the original analysis
 [d  ts  events] = DOerainroi0remove2(target_roiFileNames,groupNum,scanNum, r2thr, sprintf('crfXloc4conditions_R1'));
 
 d.r2thr = r2thr;
 d.original_ts = [];
 
 % use fitTimesries:
 disp('[meanERA_ROI_remove0_script3_R1_boot] Running the GLM fit. Using fitTimecourse.m');

 d1 = fitTimecourse(ts.subtracted_mean_roi_ts,ts.stimvol,tr,'concatInfo',ts.concatInfo,'displayFit=1', ...
  'fitType=glm','option=std','stdGroups',stdGroups);
 clear ts d;
else
 disp(sprintf('[meanERA_ROI_remove0_script3_R1_boot] Loading  file %s to recompute bootstrap',saveTag));
 if isfile([saveTag,'.mat'])
  loadFile = [saveTag,'.mat'];
 elseif ~isempty(session)
  switch session
   case {'/data2/crfaa/crfaa/FP_A0_testStream/fp20071019/'}
    loadFile = sprintf('/data2/crfaa/crfaa/data_used_files_R1/FP_A0_2010-20-July_%s_era_R1_subtrBOOT_%s.mat', ...
               target_roiFileNames{1}(1:2),target_roiFileNames{1}(1:2));
    
   case {'/data2/crfaa/crfaa/FM_A0_testStream/fm20080209/'}
    loadFile = sprintf('/data2/crfaa/crfaa/data_used_files_R1/FM_A0_2010-20-July_%s_era_R1_subtrBOOT_%s.mat', ...
               target_roiFileNames{1}(1:2),target_roiFileNames{1}(1:2));
    
   case {'/data2/crfaa/crfaa/JG_A0_testStream/jg20070919/'}
    loadFile = sprintf('/data2/crfaa/crfaa/data_used_files_R1/JG_A0_2010-20-July_%s_era_R1_subtrBOOT_%s.mat', ...
               target_roiFileNames{1}(1:2),target_roiFileNames{1}(1:2));
    
   otherwise
    keyboard
  end
  
 else
  keyboard
 end
 data = load(loadFile);
 d1     = data.d1;
 events = data.events;
 clear data;
end

% do bootstrap
d2.nBoots = 10000;
d2.bootType = 'glm';
disp(sprintf('[meanERA_ROI_remove0_script3_R1_boot]finding positive responses to compute %i bootstrap samples of CRFs',d2.nBoots));

% chose how to do the bootstrap
% either using the single-trial responses from the deconvolution
% or the
switch d2.bootType
 case {'deconv' 'd'} % take the individual trials responses from the glm
  % remove any response < 0
  for i = 1: length(stdGroups{1})
   if recomputedDeconv
    d1.numAmplitudes{i} = length(d1.allAmplitudes{i});
   end
   % remove responses below 0
   d2.minval.cutoff = -1;
   d2.minval.accept{i} = d1.allAmplitudes{i} > d2.minval.cutoff;
   
   % remove responses more than a certain z-score
   d2.zs.cutoff = 2.5;
   d2.zs.values{i} = abs(zscore(d1.allAmplitudes{i}));
   d2.zs.accept{i} = d2.zs.values{i} < d2.zs.cutoff;
   d2.accepted{i} = and(d2.zs.accept{i},d2.minval.accept{i});
   
   d2.allAmplitudes{i} = d1.allAmplitudes{i}(d2.accepted{i});
   disp(sprintf('number of deleted amplitudes: %i/%i',length(find(~d2.accepted{i})),d1.numAmplitudes{i}));
   
  end
  
  % create boot-samples
  disppercent(-inf,sprintf('[meanERA_ROI_remove0_script3_R1_boot] Computing %i bootstrap samples of CRFs',d2.nBoots));
  for bt = 1:d2.nBoots
   disppercent(bt/nBoots,sprintf('[meanERA_ROI_remove0_script3_R1_boot] Computing %i bootstrap samples of CRFs',d2.nBoots));
   for i = 1:length(stdGroups{1})
    sample = randsample(d2.allAmplitudes{i},d1.numAmplitudes{i},true);
    d2.amplitude(i,bt) = mean(sample);
    d2.amplitudeSTE(i,bt) = std(sample)/sqrt(d1.numAmplitudes{i});
   end
  end
  
 case {'glm' 'GLM'}
  % create boot-samples
  disppercent(-inf,sprintf('[meanERA_ROI_remove0_script3_R1_boot] Computing %i bootstrap samples of CRFs',d2.nBoots));
  d2.amplitude    = nan(length(stdGroups{1}),d2.nBoots);
  for i = 1: length(stdGroups{1})
   d2.amplitude(i,:)    = d1.amplitude(i) + d1.amplitudeSTE(i).*randn(d2.nBoots,1);
  end
  
 otherwise
  keyboard
end

disppercent(inf,sprintf('[meanERA_ROI_remove0_script3_R1_boot] Done computing %i bootstrap samples of CRFs',d2.nBoots));

% plot the bootstrapped crfs
plot(d2.amplitude)
title('Contrast responses after cut-off')
ylabel(sprintf('BOLD response\n(% signal change)'))
drawnow

disp(sprintf('[meanERA_ROI_remove0_script3_R1_boot] Saving d1, d2, and events as %s',saveTag));
eval(sprintf('save(''%s'',''d2'', ''d1'',''events'')', saveTag))

% plotting and saving figures:
disp(sprintf('[meanERA_ROI_remove0_script3_R1_boot] DONE ANALYSIS ON group %i, scan %i roi: %s',groupNum,scanNum,target_roiFileNames{1}(1:2)));
disp(sprintf('[meanERA_ROI_remove0_script3_R1_boot] pwd = [%s]',pwd));

close all
toc


%%%%%%%%%%%%%%%%%%%%
% switchVisualArea %
%%%%%%%%%%%%%%%%%%%%
function target_roiFileNames = switchVisualArea(target_roiFileNames)

switch target_roiFileNames
 case 'v1'
  target_roiFileNames = {'v1_lhv_ruq_loc1_0_7', 'v1_lhd_rdq_loc2_0_7', 'v1_rhd_ldq_loc3_0_7', 'v1_rhv_luq_loc4_0_7'};
 case 'v1r030'
  target_roiFileNames = {'v1_lhv_ruq_loc1_r2_03r2', 'v1_lhd_rdq_loc2_r2_03r2', 'v1_rhd_ldq_loc3_r2_03r2', 'v1_rhv_luq_loc4_r2_03r2'};
 case 'v1r050'
  target_roiFileNames = {'v1_lhv_ruq_loc1_r2_05r2',  'v1_lhd_rdq_loc2_r2_05r2',  'v1_rhd_ldq_loc3_r2_05r2',  'v1_rhv_luq_loc4_r2_05r2'};
 case 'v1r060'
  target_roiFileNames = {'v1_lhv_ruq_loc1_r2_06r2',  'v1_lhd_rdq_loc2_r2_06r2',  'v1_rhd_ldq_loc3_r2_06r2',  'v1_rhv_luq_loc4_r2_06r2'};
 case 'v1r070'
  target_roiFileNames = {'v1_lhv_ruq_loc1_r2_07r2',  'v1_lhd_rdq_loc2_r2_07r2',  'v1_rhd_ldq_loc3_r2_07r2',  'v1_rhv_luq_loc4_r2_07r2'};
 case 'v1r075'
  target_roiFileNames = {'v1_lhv_ruq_loc1_r2_075r2', 'v1_lhd_rdq_loc2_r2_075r2', 'v1_rhd_ldq_loc3_r2_075r2', 'v1_rhv_luq_loc4_r2_075r2'};
 case 'v1r085'
  target_roiFileNames = {'v1_lhv_ruq_loc1_r2_085r2', 'v1_lhd_rdq_loc2_r2_085r2', 'v1_rhd_ldq_loc3_r2_085r2', 'v1_rhv_luq_loc4_r2_085r2'};
 case 'v1r090'
  target_roiFileNames = {'v1_lhv_ruq_loc1_r2_09r2', 'v1_lhd_rdq_loc2_r2_09r2', 'v1_rhd_ldq_loc3_r2_09r2', 'v1_rhv_luq_loc4_r2_09r2'};
  
 case 'v2'
  target_roiFileNames = {'v2_lhv_ruq_loc1_0_5', 'v2_lhd_rdq_loc2_0_5', 'v2_rhd_ldq_loc3_0_5', 'v2_rhv_luq_loc4_0_5'};
 case 'v2r030'
  target_roiFileNames = {'v2_lhv_ruq_loc1_r2_03r2', 'v2_lhd_rdq_loc2_r2_03r2', 'v2_rhd_ldq_loc3_r2_03r2', 'v2_rhv_luq_loc4_r2_03r2'};
 case 'v2r050'
  target_roiFileNames = {'v2_lhv_ruq_loc1_r2_05r2',  'v2_lhd_rdq_loc2_r2_05r2',  'v2_rhd_ldq_loc3_r2_05r2',  'v2_rhv_luq_loc4_r2_05r2'};
 case 'v2r060'
  target_roiFileNames = {'v2_lhv_ruq_loc1_r2_06r2',  'v2_lhd_rdq_loc2_r2_06r2',  'v2_rhd_ldq_loc3_r2_06r2',  'v2_rhv_luq_loc4_r2_06r2'};
 case 'v2r070'
  target_roiFileNames = {'v2_lhv_ruq_loc1_r2_07r2',  'v2_lhd_rdq_loc2_r2_07r2',  'v2_rhd_ldq_loc3_r2_07r2',  'v2_rhv_luq_loc4_r2_07r2'};
 case 'v2r075'
  target_roiFileNames = {'v2_lhv_ruq_loc1_r2_075r2', 'v2_lhd_rdq_loc2_r2_075r2', 'v2_rhd_ldq_loc3_r2_075r2', 'v2_rhv_luq_loc4_r2_075r2'};
 case 'v2r085'
  target_roiFileNames = {'v2_lhv_ruq_loc1_r2_085r2', 'v2_lhd_rdq_loc2_r2_085r2', 'v2_rhd_ldq_loc3_r2_085r2', 'v2_rhv_luq_loc4_r2_085r2'};
 case 'v2r090'
  target_roiFileNames = {'v2_lhv_ruq_loc1_r2_09r2', 'v2_lhd_rdq_loc2_r2_09r2', 'v2_rhd_ldq_loc3_r2_09r2', 'v2_rhv_luq_loc4_r2_09r2'};
  
 case 'v3'
  target_roiFileNames = {'v3_lhv_ruq_loc1_0_5', 'v3_lhd_rdq_loc2_0_5', 'v3_rhd_ldq_loc3_0_5', 'v3_rhv_luq_loc4_0_5'};
 case 'v3r030'
  target_roiFileNames = {'v3_lhv_ruq_loc1_r2_03r2', 'v3_lhd_rdq_loc2_r2_03r2', 'v3_rhd_ldq_loc3_r2_03r2', 'v3_rhv_luq_loc4_r2_03r2'};
 case 'v3r070'
  target_roiFileNames = {'v3_lhv_ruq_loc1_r2_07r2',  'v3_lhd_rdq_loc2_r2_07r2',  'v3_rhd_ldq_loc3_r2_07r2',  'v3_rhv_luq_loc4_r2_07r2'};
  
 case 'v4'
  target_roiFileNames = {'v4_lhv_ruq_loc1_0_5', 'v4_lhv_rdq_loc2_0_5', 'v4_rhv_ldq_loc3_0_5', 'v4_rhv_luq_loc4_0_5'};
 case 'v4r030'
  target_roiFileNames = {'v4_lhv_ruq_loc1_r2_03r2', 'v4_lhv_rdq_loc2_r2_03r2', 'v4_rhv_ldq_loc3_r2_03r2', 'v4_rhv_luq_loc4_r2_03r2'};
 case 'v4r070'
  target_roiFileNames = {'v4_lhv_ruq_loc1_r2_07r2',  'v4_lhv_rdq_loc2_r2_07r2',  'v4_rhv_ldq_loc3_r2_07r2',  'v4_rhv_luq_loc4_r2_07r2'};
  
 otherwise
  keyboard
end

