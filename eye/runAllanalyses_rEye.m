function runAllanalyses_rEye(roinum)
% function runAllanalyses_rEye(roinum)
%
% this function runs all the analyses to do generate
% the stimfiles for the event realted analysis of the eye movement
% for the supplementary material
%
% it only runs for one session, a session with one eyepos.mat file.
% the session must have a concatenation group which is the concatenation
% of all thes scans
%
% the fucntion saves a .mat file (e.g., testEyeLoc1.mat) in the dir of the session 
%
% to visualize the eye position load the file and use 
% eraEypos.m to plot the results.
%
% franco pestilli 2009.10.05

saveFigure = 1;
SESSIONS = {'/Volumes/homosacer_data/riken/crfaa/fmridata/crfaa/JG_A0_eye/jg20071023'};
%           '/Volumes/homosacer_data/riken/crfaa/fmridata/crfaa/FP_A0_eye/fp20071023' ...
%           '/Volumes/homosacer_data/riken/crfaa/fmridata/crfaa/FP_A0_eye/fp20071027' ...
%           '/Volumes/homosacer_data/riken/crfaa/fmridata/crfaa/JG_A0_eye/jg20070919',...
%           '/Volumes/homosacer_data/riken/crfaa/fmridata/crfaa/JG_A0_eye/jg20071019'};

for s = 1:length(SESSIONS)
 for roi = 1:length(roinum)
 disp('%  START ')
 disp('%  START ')
 disp(sprintf('###################### %s #############################',SESSIONS{s}))
 cd(SESSIONS{s})
 filename = sprintf('testEyeLoc%i',roinum(roi));

 close all; gs = getGS;
 [d d1 ts events] = meanERA_ROI_remove0_script3eye(gs(1),gs(2),roinum(roi),filename,'v1');
 clear global MLR;
 clear d d1 events;
 
 % visualize the eye position:
 eraEyepos(gs,'scanInfo',ts)
 end
end
disp(sprintf('###################### %s #############################',SESSIONS{s}))
disp('%  END ')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% meanERA_ROI_remove0_script3eye %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [d d1 ts events] = meanERA_ROI_remove0_script3eye(groupNum,scanNum,roinum,saveTag,target_roiFileNames)
% this file calls a bunch of functions to do several types of projection
% out
% it has a cutomized analysis for each subject and adaptation condition
%
%
% fp 2009.10.03
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
[d ts events] = DOerainroi0remove2eye(target_roiFileNames,groupNum,scanNum, roinum, sprintf('cueXloc'));
d1 = [];

% svae useful info in the results
ts.filename = saveTag;
ts.mfilename = mfilename('fullpath');
ts.currentdir = pwd;

disp(sprintf('[meanERA_ROI_remove0_script3eye] Saving d, d1, ts and events as %s',saveTag));
eval(sprintf('save(''%s'',''ts'')', saveTag))

disp(sprintf('[meanERA_ROI_remove0_script3eye] DONE ANALYSIS ON group %i, scan %i roi: %s',groupNum,scanNum,target_roiFileNames{1}(1:2)));
disp(sprintf('[meanERA_ROI_remove0_script3eye] pwd = [%s]',pwd));


%%%%%%%%%%%%%%%%%%%%%%%%%
% DOerainroi0remove2eye %
%%%%%%%%%%%%%%%%%%%%%%%%%
function [d timeSeries events] = DOerainroi0remove2eye(roinames,groupNum,scanNum,roinum,preprocessFun)
% function [d timeSeries events] = DOerainroi0remove2eye(roinames,groupNum,scanNum,roinum,preprocessFun)
%      usage: DOerainroi0remove2eye()
%         by: franco pestilli
%       date: 08/06/10
%    purpose: 0. does analysis concatenating 4 quadrant rois.
%             1. uses doerainroi to compute a single timecourse er in an roi.
%             2. it multiply the scm with the response to 0% contrast
%             3. it subtracts the obtained tseries from the original tseries.
%
% if you don't know how to change something look up how it is done in the
% GUI by calling: guide('mrLoadRetGUI.fig') click on the 'Menu Editor icon'and then edit
% the GUI function (edit mrLoadRetGUI) to look for the call.
%
% franco pestilli

% check arguments
if ~any(nargin == 5)
 help DOerainroi0remove2eye
 mean_hdr = [];
 timeSeries = [];
 return
end

% 1. get responses in each roi and scm:
% hdrLengthInSeconds = 19;
% open new view:
v = newView;
groupName = viewGet(v,'groupName',groupNum);
disp(sprintf('[DOerainroi0remove2eye] Working on ALL ROIs, in group %s, scan %i.',groupName,scanNum));

% load scan info:
d = loadScan(v,scanNum,groupNum,0);

% stack together timeseries and scm for all rois
scm = [];

% select the right volumes/conditions:
disp(sprintf('[DOerainroi0remove2eye] computing stimVols with d=%s, on roi: %s',preprocessFun,roinames{roinum}));
preprocessType = sprintf('run=%s(d,%s)',preprocessFun,num2str(roinum));
d = eventRelatedPreProcess(d,preprocessType);

% save concatenation info for each roi (needed for fitTimecourse):
timeSeries.concatInfo{1} = d.concatInfo;
timeSeries.stimvol{1} = d.stimvol;

disp(sprintf('[DOerainroi0remove2eye] stimVols computed with d=%s, on roi: %s',preprocessFun,roinames{roinum}));

% getting all info out of analysis:
tr = d.tr;
hdrlen = 10;% from seconds (19) to trs (24)
hdrNum = d.hdrNum;

events = d.e;

% make a stimulus convolution matrix
% (the '1' option applies a filter to it) :
disp(sprintf('[DOerainroi0remove2eye] make scm for roi %s and stack it with the rest of the scm',roinames{roinum}));
d = makescm(d,hdrlen,1);
timeSeries.stimvol{1} = d.stimvol;
d.scm = [];

% % save the concatInfo:
timeSeries = concatInfoMultyROI(timeSeries);
timeSeries = concatStimvol(timeSeries,1);
timeSeries.dim = [d.dim(1:3) size(scm,1)];
timeSeries.location = roinum;

% clear d and view
clear roi_ts d;deleteView(v);

d = [];

%%%%%%%%%%%%%%%%%%%%%%
% concatInfoMultyROI %
%%%%%%%%%%%%%%%%%%%%%%
function d = concatInfoMultyROI(d)
% function d = concatInfoMultyROI(d)
%
% function to concatenate the concatenateInfo for cases in which several
% time series are concatenated
%
% fp 2008.07.29

if ~isfield(d,'concatInfo')
 disp('[ConcatInfoMultyROI] ''d'' structure has no concatInfo field');
 return
end
if length(d.concatInfo) > 1
 % check that all concatenation info is the same:
 for ii = 1:length(d.concatInfo)-1
  fieldNames1 = fieldnames(d.concatInfo{ii});
  fieldNames2 = fieldnames(d.concatInfo{ii+1});
  ok = strcmp(fieldNames1,fieldNames2);
  if ~any(ok)
   disp('[ConcatInfoMultyROI] ''d'' structure has no concatInfo field');
   return
  end
 end

 % grab the fields out of one concatenation:
 fieldNames = fieldnames(d.concatInfo{1});

 % inititalizing field names to []:
 for jj = 1:length(fieldNames)
  concatInfo.(fieldNames{jj}) = [];
 end

 % find field names that do not necessitate special treatment:
 fieldNames = setdiff(fieldNames,{'n','whichScan','runTransition'});

 % run through all the concatenated info and stack them:
 for ii = 1:length(d.concatInfo)
  % run through the field names except the one in fieldNames:
  for jj = 1:length(fieldNames)
   concatInfo.(fieldNames{jj}) = [concatInfo.(fieldNames{jj}) d.concatInfo{ii}.(fieldNames{jj})];
  end
  % these fields get a special treatment
  if ii == 1
   concatInfo.runTransition = d.concatInfo{ii}.runTransition;
   concatInfo.n = d.concatInfo{ii}.n;
   concatInfo.whichScan = d.concatInfo{ii}.whichScan;
  else
   concatInfo.runTransition = [concatInfo.runTransition; max(max(concatInfo.runTransition))+d.concatInfo{ii}.runTransition];
   concatInfo.whichScan = [concatInfo.whichScan d.concatInfo{ii}.whichScan+concatInfo.n];
   concatInfo.n = concatInfo.n + d.concatInfo{ii}.n;
  end
 end
 d.concatInfo = concatInfo;

else
 d.concatInfo = d.concatInfo{1};
end

%%%%%%%%%%%%%%%%%
% concatStimvol %
%%%%%%%%%%%%%%%%%
function d = concatStimvol(d,roinum)
% function d = concatStimvol(d,roinum)
%
% function to concatenate the stimvol for cases in which several
% time series are concatenated
%
% fp 2008.07.29

if roinum > 1
 if ~isfield(d,'stimvol')
  disp('[ConcatInfoMultyROI] ''d'' structure has no concatInfo field');
  return
 end

 % initializing d2 for each stimulus condition:
 for jj = 1:length(d.stimvol{1})
  d2.stimvol{jj} = [];
 end

 timeSeriesLength = length(d.concatInfo.whichVolume)/roinum;
 keyboard

 % run through the rois:
 for ii = 1:roinum
  % run through the conditions
  for jj = 1:length(d.stimvol{ii})
   if ii == 1
    d2.stimvol{jj} = [d2.stimvol{jj} d.stimvol{ii}{jj}];
   else
    %             d2.stimvol{jj} = [d2.stimvol{jj} max(d2.stimvol{jj})+d.stimvol{ii}{jj}];
    d2.stimvol{jj} = [d2.stimvol{jj} d.stimvol{ii}{jj}+(ii-1)*timeSeriesLength];
   end
  end
 end
 d.stimvol = d2.stimvol;

end



%%%%%%%%%%%%%%%%%%%%
% switchVisualArea %
%%%%%%%%%%%%%%%%%%%%
function target_roiFileNames = switchVisualArea(target_roiFileNames)

switch target_roiFileNames
 case 'v1'
  target_roiFileNames = {'v1_lhv_ruq_loc2_0_7', 'v1_lhd_rdq_loc2_0_7', 'v1_rhd_ldq_loc3_0_7', 'v1_rhv_luq_loc4_0_7'};
  
 otherwise
  keyboard
end


