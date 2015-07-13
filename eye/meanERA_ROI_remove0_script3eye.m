function [d d1 ts events] = meanERA_ROI_remove0_script3eye(groupNum,scanNum,r2thr,saveTag,target_roiFileNames)
% this file calls a bunch of functions to do several types of projection
% out
% it has a cutomized analysis for each subject and adaptation condition
%
%
% fp 2009.10.03

clc
tic

% check if exist a mrSession.mat
if ~(exist('mrSession.mat','file') == 2)
 disp(sprintf('[meanERA_ROI_remove0_script3eye] mrSession.mat not found in this directory. please cd to a mrLoadret dir'));
 return
end

if ieNotDefined('saveTag')
 saveTag = num2str(100*r2thr);
end

if nargin == 3 | nargin == 4
 % the new rois after freeSurfer:
 target_roiFileNames = switchVisualArea('v1');
 saveTag = [date,'_',target_roiFileNames{1}(1:2),'_',saveTag];

elseif nargin == 5
 disp(sprintf('[meanERA_ROI_remove0_script3eye] all paramteters passed by users'));
 target_roiFileNames = switchVisualArea(target_roiFileNames);
 saveTag = [date,'_',target_roiFileNames{1}(1:2),'_',saveTag];

else
 d  = [];d1  = [];ts  = [];events = [];
 disp(sprintf('[meanERA_ROI_remove0_script3eye] ERROR ... nargin must be 3, 4 or 6 ...'));
 return
end
disp(sprintf('[meanERA_ROI_remove0_script3eye] DO ANALYSIS ON group %i, scan %i roi: %s',groupNum,scanNum,target_roiFileNames{1}(1:2)));

% this is the original analysis
[d ts events] = DOerainroi0remove2eye(target_roiFileNames,groupNum,scanNum, r2thr, sprintf('cueXloc'));
d1 = [];

disp(sprintf('[meanERA_ROI_remove0_script3eye] Saving d, d1, ts and events as %s',saveTag));
eval(sprintf('save(''%s'',''ts'')', saveTag))

% plotting and saving figures:
disp(sprintf('[meanERA_ROI_remove0_script3eye] DONE ANALYSIS ON group %i, scan %i roi: %s',groupNum,scanNum,target_roiFileNames{1}(1:2)));
disp(sprintf('[meanERA_ROI_remove0_script3eye] pwd = [%s]',pwd));

close all
toc

% switchVisualArea
function target_roiFileNames = switchVisualArea(target_roiFileNames)

switch target_roiFileNames
 case 'v1'
  target_roiFileNames = {'v1_lhv_ruq_loc2_0_7', 'v1_lhd_rdq_loc2_0_7', 'v1_rhd_ldq_loc3_0_7', 'v1_rhv_luq_loc4_0_7'};
  
 otherwise
  keyboard
end
