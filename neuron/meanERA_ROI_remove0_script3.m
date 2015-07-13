function [d d1 ts events] = meanERA_ROI_remove0_script3(groupNum,scanNum,r2thr,saveTag,target_roiFileNames)
% function [d d1 ts events] = meanERA_ROI_remove0_script3(groupNum,scanNum,r2thr,saveTag,target_roiFileNames)
%
% this file calls a bunch of functions to do several types of projection
% out
% it has a cutomized analysis for each subject and adaptation condition
%
% fp 2008.07.29

clc
tic

% check if exist a mrSession.mat
if ~(exist('mrSession.mat','file') == 2)
 disp(sprintf('[meanERA_ROI_remove0_script3] mrSession.mat not found in this directory. please cd to a mrLoadret dir'));
 return
end

if ieNotDefined('saveTag')
 saveTag = num2str(100*r2thr);
end

if nargin == 3 || nargin == 4
 % the new rois after freeSurfer:
 target_roiFileNames = switchVisualArea('v1');
 saveTag = ['2010-30-May','_',target_roiFileNames{1}(1:2),'_',saveTag];

elseif nargin == 5
 disp(sprintf('[meanERA_ROI_remove0_script3] all paramteters passed by users'));
 target_roiFileNames = switchVisualArea(target_roiFileNames);
 saveTag = ['2010-30-May','_',target_roiFileNames{1}(1:2),'_',saveTag];

else
 d  = [];d1  = [];ts  = [];events = [];
 disp(sprintf('[meanERA_ROI_remove0_script3] ERROR ... nargin must be 3, 4 or 6 ...'));
 return
end
disp(sprintf('[meanERA_ROI_remove0_script3] DO ANALYSIS ON group %i, scan %i roi: %s',groupNum,scanNum,target_roiFileNames{1}(1:2)));

% this is the original analysis
[d  ts  events] = DOerainroi0remove2(target_roiFileNames,groupNum,scanNum, r2thr, sprintf('crfXloc4conditions'));

% this is the analysis for correct and incorrect responses 
% [d  ts  events] = DOerainroi0remove2(target_roiFileNames,groupNum,scanNum, r2thr, sprintf('crfXloc4conditionsCorrIncor'));

% this is the analysis for 1st and 2nd interval to look at the effect of transient attention to generate the distributed target/non-target difference 
% [d  ts  events] = DOerainroi0remove2(target_roiFileNames,groupNum,scanNum, r2thr, sprintf('crfXloc4conditions1st2ndInterv'));

d.r2thr = r2thr;
ts.original_ts = [];

% plot time series to check that everything is good:
% figure(100)
% plot(1:length(ts.original_ts),ts.original_ts,'r', 1:length(ts.original_ts),.9+ts.estimated_ts,'g',1:length(ts.original_ts),.1+ts.subtracted_mean_roi_ts,'b')
% legend({'original','estimated','subtracted'});
% title('concatenated tSeries')
% drawnow;
% disp(sprintf('Saving figure %s',[saveTag,'_ts']));
% eval(sprintf('print(gcf,''-depsc2'',''-tiff'',''-r300'', ''%s'')', [saveTag,'_ts']));

% use fitTimesries:
tr =.8;
% d1 = fitTimecourse(ts.subtracted_mean_roi_ts,ts.stimvol,tr,'concatInfo',ts.concatInfo,'displayFit=1','fitType=glm');
% d1 = fitTimecourse(ts.subtracted_mean_roi_ts,ts.stimvol,tr,'concatInfo',ts.concatInfo,'displayFit=1','fitType=glm','option=std','stdGroups',{1:8 9:16 17:24 25:32});
% d1 = fitTimecourse(ts.subtracted_mean_roi_ts,ts.stimvol,tr,'concatInfo',ts.concatInfo,'displayFit=1','fitType=glm','option=std');

% try the analysis on the time series before subtraction (riken, 2009.06.29)
% d1 = fitTimecourse(ts.original_ts,ts.stimvol,tr,'concatInfo',ts.concatInfo,'displayFit=1','fitType=glm');

% this is the original analysis of the response and of the standard deviation
d1 = fitTimecourse(ts.subtracted_mean_roi_ts,ts.stimvol,tr,'concatInfo',ts.concatInfo,'displayFit=1','fitType=glm','option=std','stdGroups',{1:8 9:16 17:24 25:32 });

% this is the analysis for correct and incorrect responses or for frist and
% second intervals
% d1 = fitTimecourse(ts.subtracted_mean_roi_ts,ts.stimvol,tr,'concatInfo',ts.concatInfo,'displayFit=1','fitType=glm','option=std','stdGroups',{1:8 9:16 17:24 25:32 33:40 41:48});


disp(sprintf('[meanERA_ROI_remove0_script3] Saving d, d1, ts and events as %s',saveTag));
eval(sprintf('save(''%s'',''d'', ''d1'',''ts'',''events'')', saveTag))

% % plot and save figures:roiTest
% disp(sprintf('[meanERA_ROI_remove0_script3] Plotting results of group %i, scan %i roi: %s',groupNum,scanNum,target_roiFileNames{1}(1:2)));
% % % plot and save figures:roiTest
c = plotCRFtimecourse(d1,events,saveTag);

% plotting and saving figures:
disp(sprintf('[meanERA_ROI_remove0_script3] DONE ANALYSIS ON group %i, scan %i roi: %s',groupNum,scanNum,target_roiFileNames{1}(1:2)));
disp(sprintf('[meanERA_ROI_remove0_script3] pwd = [%s]',pwd));

close all
toc

function target_roiFileNames = switchVisualArea(target_roiFileNames)

switch target_roiFileNames
 case 'v1r030'
  target_roiFileNames = {'v1_lhv_ruq_loc1_r2_03r2', 'v1_lhd_rdq_loc2_r2_03r2', 'v1_rhd_ldq_loc3_r2_03r2', 'v1_rhv_luq_loc4_r2_03r2'};
 case 'v1r050'
  target_roiFileNames = {'v1_lhv_ruq_loc1_r2_05r2',  'v1_lhd_rdq_loc2_r2_05r2',  'v1_rhd_ldq_loc3_r2_05r2',  'v1_rhv_luq_loc4_r2_05r2'};
 case 'v1r060'
  target_roiFileNames = {'v1_lhv_ruq_loc1_r2_06r2',  'v1_lhd_rdq_loc2_r2_06r2',  'v1_rhd_ldq_loc3_r2_06r2',  'v1_rhv_luq_loc4_r2_06r2'};
 case 'v1'
  target_roiFileNames = {'v1_lhv_ruq_loc1_0_7', 'v1_lhd_rdq_loc2_0_7', 'v1_rhd_ldq_loc3_0_7', 'v1_rhv_luq_loc4_0_7'};
 case 'v1r070'
  target_roiFileNames = {'v1_lhv_ruq_loc1_r2_07r2',  'v1_lhd_rdq_loc2_r2_07r2',  'v1_rhd_ldq_loc3_r2_07r2',  'v1_rhv_luq_loc4_r2_07r2'};
 case 'v1r075'
  target_roiFileNames = {'v1_lhv_ruq_loc1_r2_075r2', 'v1_lhd_rdq_loc2_r2_075r2', 'v1_rhd_ldq_loc3_r2_075r2', 'v1_rhv_luq_loc4_r2_075r2'};
 case 'v1r085'
  target_roiFileNames = {'v1_lhv_ruq_loc1_r2_085r2', 'v1_lhd_rdq_loc2_r2_085r2', 'v1_rhd_ldq_loc3_r2_085r2', 'v1_rhv_luq_loc4_r2_085r2'};
 case 'v1r090'
  target_roiFileNames = {'v1_lhv_ruq_loc1_r2_09r2', 'v1_lhd_rdq_loc2_r2_09r2', 'v1_rhd_ldq_loc3_r2_09r2', 'v1_rhv_luq_loc4_r2_09r2'};
  
 case 'v2r030'
  target_roiFileNames = {'v2_lhv_ruq_loc1_r2_03r2', 'v2_lhd_rdq_loc2_r2_03r2', 'v2_rhd_ldq_loc3_r2_03r2', 'v2_rhv_luq_loc4_r2_03r2'};
 case 'v2'
  target_roiFileNames = {'v2_lhv_ruq_loc1_0_5', 'v2_lhd_rdq_loc2_0_5', 'v2_rhd_ldq_loc3_0_5', 'v2_rhv_luq_loc4_0_5'};
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
  
 case 'v3r030'
  target_roiFileNames = {'v3_lhv_ruq_loc1_r2_03r2', 'v3_lhd_rdq_loc2_r2_03r2', 'v3_rhd_ldq_loc3_r2_03r2', 'v3_rhv_luq_loc4_r2_03r2'};
 case 'v3'
  target_roiFileNames = {'v3_lhv_ruq_loc1_0_5', 'v3_lhd_rdq_loc2_0_5', 'v3_rhd_ldq_loc3_0_5', 'v3_rhv_luq_loc4_0_5'};
 case 'v3r070'
  target_roiFileNames = {'v3_lhv_ruq_loc1_r2_07r2',  'v3_lhd_rdq_loc2_r2_07r2',  'v3_rhd_ldq_loc3_r2_07r2',  'v3_rhv_luq_loc4_r2_07r2'};
  
 case 'v4r030'
  target_roiFileNames = {'v4_lhv_ruq_loc1_r2_03r2', 'v4_lhv_rdq_loc2_r2_03r2', 'v4_rhv_ldq_loc3_r2_03r2', 'v4_rhv_luq_loc4_r2_03r2'};
 case 'v4'
  target_roiFileNames = {'v4_lhv_ruq_loc1_0_5', 'v4_lhv_rdq_loc2_0_5', 'v4_rhv_ldq_loc3_0_5', 'v4_rhv_luq_loc4_0_5'};
 case 'v4r070'
  target_roiFileNames = {'v4_lhv_ruq_loc1_r2_07r2',  'v4_lhv_rdq_loc2_r2_07r2',  'v4_rhv_ldq_loc3_r2_07r2',  'v4_rhv_luq_loc4_r2_07r2'};
  
 otherwise
  keyboard
end


%  case 'v1Noise'
%   target_roiFileNames = {'v1_lhv_ruq_loc1_noise', 'v1_lhd_rdq_loc2_noise', 'v1_rhd_ldq_loc3_noise', 'v1_rhv_luq_loc4_noise'};
%  case 'v1Cerebellum'
%   target_roiFileNames = {'v1_lhv_ruq_loc1_cerebellum', 'v1_lhd_rdq_loc2_cerebellum', 'v1_rhd_ldq_loc3_cerebellum', 'v1_rhv_luq_loc4_cerebellum'};
%  case 'v1Skull'
%   target_roiFileNames = {'v1_lhv_ruq_loc1_skull', 'v1_lhd_rdq_loc2_skull', 'v1_rhd_ldq_loc3_skull', 'v1_rhv_luq_loc4_skull'};
%  case 'v1CSF'
%   target_roiFileNames = {'v1_lhv_ruq_loc1_csf', 'v1_lhd_rdq_loc2_csf', 'v1_rhd_ldq_loc3_csf', 'v1_rhv_luq_loc4_csf'};
%  case 'v2'
%   target_roiFileNames = {'v2_lhv_ruq_loc1_0_5', 'v2_lhd_rdq_loc2_0_5', 'v2_rhd_ldq_loc3_0_5', 'v2_rhv_luq_loc4_0_5'};
%  case 'v3'
%   target_roiFileNames = {'v3_lhv_ruq_loc1_0_5', 'v3_lhd_rdq_loc2_0_5', 'v3_rhd_ldq_loc3_0_5', 'v3_rhv_luq_loc4_0_5'};
%  case 'v4'
%   target_roiFileNames = {'v4_lhv_ruq_loc1_0_5', 'v4_lhv_rdq_loc2_0_5', 'v4_rhv_ldq_loc3_0_5', 'v4_rhv_luq_loc4_0_5'};
%  
%  case 'combined'
%   target_roiFileNames = {'v1_3_lhv_ruq_loc1', 'v1_3_lhd_rdq_loc2', 'v1_3_rhd_ldq_loc3', 'v1_3_rhv_luq_loc4'};
%  case 'combinedR' % restircted to a stricter r2
%   target_roiFileNames = {'v1_3_lhv_ruq_loc1_r', 'v1_3_lhd_rdq_loc2_r', 'v1_3_rhd_ldq_loc3_r', 'v1_3_rhv_luq_loc4_r'};
%  case 'combinedR2' % restircted to a stricter r2
%   target_roiFileNames = {'v1_3_lhv_ruq_loc1_r2', 'v1_3_lhd_rdq_loc2_r2', 'v1_3_rhd_ldq_loc3_r2', 'v1_3_rhv_luq_loc4_r2'};
%  case 'combinedR203' % restircted to a stricter r2
%   target_roiFileNames = {'v1_3_lhv_ruq_loc1_r2_03r2', 'v1_3_lhd_rdq_loc2_r2_03r2', 'v1_3_rhd_ldq_loc3_r2_03r2', 'v1_3_rhv_luq_loc4_r2_03r2'};
%  case 'combinedR205' % restircted to a stricter r2
%   target_roiFileNames = {'v1_3_lhv_ruq_loc1_r2_05r2', 'v1_3_lhd_rdq_loc2_r2_05r2', 'v1_3_rhd_ldq_loc3_r2_05r2', 'v1_3_rhv_luq_loc4_r2_05r2'};
%  case 'combinedR206' % restircted to a stricter r2
%   target_roiFileNames = {'v1_3_lhv_ruq_loc1_r2_06r2', 'v1_3_lhd_rdq_loc2_r2_06r2', 'v1_3_rhd_ldq_loc3_r2_06r2', 'v1_3_rhv_luq_loc4_r2_06r2'};
%  case 'combinedR207' % restircted to a stricter r2
%   target_roiFileNames = {'v1_3_lhv_ruq_loc1_r2_07r2', 'v1_3_lhd_rdq_loc2_r2_07r2', 'v1_3_rhd_ldq_loc3_r2_07r2', 'v1_3_rhv_luq_loc4_r2_07r2'};
%  case 'combinedR2075' % restircted to a stricter r2
%   target_roiFileNames = {'v1_3_lhv_ruq_loc1_r2_075r2', 'v1_3_lhd_rdq_loc2_r2_075r2', 'v1_3_rhd_ldq_loc3_r2_075r2', 'v1_3_rhv_luq_loc4_r2_075r2'};
%  case 'combinedR2085' % restircted to a stricter r2
%   target_roiFileNames = {'v1_3_lhv_ruq_loc1_r2_085r2', 'v1_3_lhd_rdq_loc2_r2_085r2', 'v1_3_rhd_ldq_loc3_r2_085r2', 'v1_3_rhv_luq_loc4_r2_085r2'};
%  case 'combinedR209' % restircted to a stricter r2
%   target_roiFileNames = {'v1_3_lhv_ruq_loc1_r2_09r2', 'v1_3_lhd_rdq_loc2_r2_09r2', 'v1_3_rhd_ldq_loc3_r2_09r2', 'v1_3_rhv_luq_loc4_r2_09r2'};
%  
%  case 'v1_3r030'
%   target_roiFileNames = {'v1_3_lhv_ruq_loc1_r2_03r2', 'v1_3_lhd_rdq_loc2_r2_03r2', 'v1_3_rhd_ldq_loc3_r2_03r2', 'v1_3_rhv_luq_loc4_r2_03r2'};
%  case 'v1_3r050'
%   target_roiFileNames = {'v1_3_lhv_ruq_loc1_r2_05r2',  'v1_3_lhd_rdq_loc2_r2_05r2',  'v1_3_rhd_ldq_loc3_r2_05r2',  'v1_3_rhv_luq_loc4_r2_05r2'};
%  case 'v1_3r060'
%   target_roiFileNames = {'v1_3_lhv_ruq_loc1_r2_06r2',  'v1_3_lhd_rdq_loc2_r2_06r2',  'v1_3_rhd_ldq_loc3_r2_06r2',  'v1_3_rhv_luq_loc4_r2_06r2'};
%  case 'v1_3r070'
%   target_roiFileNames = {'v1_3_lhv_ruq_loc1_r2_07r2',  'v1_3_lhd_rdq_loc2_r2_07r2',  'v1_3_rhd_ldq_loc3_r2_07r2',  'v1_3_rhv_luq_loc4_r2_07r2'};
%  case 'v1_3r075'
%   target_roiFileNames = {'v1_3_lhv_ruq_loc1_r2_075r2', 'v1_3_lhd_rdq_loc2_r2_075r2', 'v1_3_rhd_ldq_loc3_r2_075r2', 'v1_3_rhv_luq_loc4_r2_075r2'};
%  case 'v1_3r085'
%   target_roiFileNames = {'v1_3_lhv_ruq_loc1_r2_085r2', 'v1_3_lhd_rdq_loc2_r2_085r2', 'v1_3_rhd_ldq_loc3_r2_085r2', 'v1_3_rhv_luq_loc4_r2_085r2'};
%  case 'v1_3r090'
%   target_roiFileNames = {'v1_3_lhv_ruq_loc1_r2_09r2', 'v1_3_lhd_rdq_loc2_r2_09r2', 'v1_3_rhd_ldq_loc3_r2_09r2', 'v1_3_rhv_luq_loc4_r2_09r2'};
  