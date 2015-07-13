function [d timeSeries events] = DOerainroi0remove2eye(roinames,groupNum,scanNum,minR2,preprocessFun)
% function [d timeSeries events] = DOerainroi0remove2eye(roinames,groupNum,scanNum,minR2,preprocessFun)
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

% check if exist a mrSession.mat
if ~(exist('mrSession.mat','file') == 2)
 disp(sprintf('[DOerainroi0remove2eye] mrSession.mat not found in this directory.'));
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
roinum = 2;

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

%  if roinum == 1
%   % save target contrasts, reaction, and response:
events = d.e;
%  end

% make a stimulus convolution matrix
% (the '1' option applies a filter to it) :
disp(sprintf('[DOerainroi0remove2eye] make scm for roi %s and stack it with the rest of the scm',roinames{roinum}));
d = makescm(d,hdrlen,1);
% save the scm:
%  scm = [scm; d.scm];

timeSeries.stimvol{1} = d.stimvol;

d.scm = [];

% % save the concatInfo across roi:
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