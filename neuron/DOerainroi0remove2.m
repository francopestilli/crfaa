function [d timeSeries events] = DOerainroi0remove2(roinames,groupNum,scanNum,minR2,preprocessFun)
% function [d timeSeries events] = DOerainroi0remove2(roinames,groupNum,scanNum,minR2,preprocessFun)
%      usage: DOerainroi0remove2()
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
    help DOerainroi0remove2
    mean_hdr = [];
    timeSeries = [];
    return
end

% check if exist a mrSession.mat
if ~(exist('mrSession.mat','file') == 2)
    disp(sprintf('[DOerainroi0remove2] mrSession.mat not found in this directory.'));
    return
end

% 1. get responses in each roi and scm:
hdrLengthInSeconds = 19;
% open new view:
v = newView;
groupName = viewGet(v,'groupName',groupNum);
disp(sprintf('[DOerainroi0remove2] Working on ALL ROIs, in group %s, scan %i.',groupName,scanNum));

% load scan info:
d = loadScan(v,scanNum,groupNum,0);

% stack together timeseries and scm for all rois
scm = []; mean_roi_ts = [];
for roinum = 1:length(roinames)
    % select the right volumes/conditions:
    disp(sprintf('[DOerainroi0remove2] computing stimVols with d=%s, on roi: %s',preprocessFun,roinames{roinum}));
    preprocessType = sprintf('run=%s(d,%s)',preprocessFun,num2str(roinum));
    d = eventRelatedPreProcess(d,preprocessType);
    
    % save concatenation info for each roi (needed for fitTimecourse):
    timeSeries.concatInfo{roinum} = d.concatInfo;
    timeSeries.stimvol{roinum} = d.stimvol;
    
    disp(sprintf('[DOerainroi0remove2] stimVols computed with d=%s, on roi: %s',preprocessFun,roinames{roinum}));

    % getting all info out of analysis:
    tr = d.tr;
    hdrlen = ceil(hdrLengthInSeconds/tr);% from seconds (19) to trs (24)
    hdrNum = d.hdrNum;
    
    if roinum == 1
        % save target contrasts, reaction, and response:
        events = d.e;
    end
 
    % make a stimulus convolution matrix
    % (the '1' option applies a filter to it) :
    disp(sprintf('[DOerainroi0remove2] make scm for roi %s and stack it with the rest of the scm',roinames{roinum}));
    d = makescm(d,hdrlen,1);
    % save the scm:
    scm = [scm; d.scm];
    
    timeSeries.stimvol{roinum} = d.stimvol;
   
    d.scm = [];

    % load time series from the right ROI:
    disp(sprintf('[DOerainroi0remove2] load timeseries for roi %s and stack it with the rest of the timeseries',roinames{roinum}));
    roi_ts = loadROITSeries(v,roinames{roinum},scanNum,groupNum);
    
    % change any ZERO in the tseries to nan to work properly with nanmean:
    roi_ts.tSeries(roi_ts.tSeries==0) = nan;
    
    % stack with the rest of the time series
    mean_roi_ts = [mean_roi_ts squeeze(nanmean(roi_ts.tSeries,1))];

    disp(sprintf('[DOerainroi0remove2] timeseries loaded for roi %s and stacked with the rest of scm (size tSeries=[%s]) ',roinames{roinum},num2str(size(mean_roi_ts))));
end

% % save the concatInfo across roi:
timeSeries = concatInfoMultyROI(timeSeries);
timeSeries = concatStimvol(timeSeries);
timeSeries.dim = [d.dim(1:3) size(scm,1)];

% clear d and view
clear roi_ts d;deleteView(v);

% computing hdr on average response in the roi:
disp('[DOerainroi0remove2] compute mean response on thresholded mean roi time series');
d.BeforeSubtraction_d  = getr2timecourse(mean_roi_ts,hdrNum,hdrlen,scm,tr);

% empting the data and scm fields to save memory
mean_hdr.data = [];
mean_hdr.scm = [];


% 2. select the hdr for 0 contrast unattended:
disp('[DOerainroi0remove2] now removing the response to the 0%-contrast unattended stimulus');
if hdrNum == 48
 unattended0index = 17;
elseif hdrNum == 32
 unattended0index = 9;
else
 disp('[DOerainroi0remove2] canno''t find the unattended response to the 0-pedestal contrast wth this number of HRF, please check the following lines.');
 keyboard
end
zeroContrastResponse = (d.BeforeSubtraction_d.ehdr(unattended0index,:)/100); % here i divide by 100 because getr2timecourse returns results in % [0-100] i want [0-1]

% 3. convolve 0-contrast response with event model
% % (this is equivalent to multipling for the scm, which is actually what i do here):
m = repmat(zeroContrastResponse,1,hdrNum)';
estimatedTseries = scm*m;

clear m;

% 4. subtract estimated 0-contrast response out of the original timeseries:
% % do subtraction: on each voxel
disp('[DOerainroi0remove2] now subtracting the estimated time series from the 0%-contrast unattended stimulus');
disp('[DOerainroi0remove2] off the time course of the voxels abve the r2 threshold');
subtracted_mean_roi_ts = mean_roi_ts - estimatedTseries';

% up to here launched
% 5. re-do era on the new tseries:
% computing average hdr in roi:
disp('[DOerainroi0remove2] now compute mean response on subtracted mean roi time series');
d.AfterSubtraction_d  = getr2timecourse(subtracted_mean_roi_ts,hdrNum,hdrlen,scm,tr);
clear scm;

% empting the data and scm fields to save memory
mean_hdr.data = [];
mean_hdr.scm = [];

% saving time series for output:
timeSeries.original_ts = mean_roi_ts;
timeSeries.estimated_ts = estimatedTseries;
timeSeries.subtracted_mean_roi_ts = subtracted_mean_roi_ts;
timeSeries.zeroContrastResponse = zeroContrastResponse;

disp(sprintf('[DOerainroi0remove2] Worked on ALL ROIs, in group %s, scan %i.',groupName,scanNum));


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


%%%%%%%%%%%%%%%%%
% concatStimvol %
%%%%%%%%%%%%%%%%%
function d = concatStimvol(d)
% function d = concatStimvol(d)
%
% function to concatenate the stimvol for cases in which several
% time series are concatenated
%
% fp 2008.07.29

if ~isfield(d,'stimvol')
    disp('[ConcatInfoMultyROI] ''d'' structure has no concatInfo field');
    return
end

% initializing d2 for each stimulus condition:
for jj = 1:length(d.stimvol{1})
    d2.stimvol{jj} = [];
end

timeSeriesLength = length(d.concatInfo.whichVolume)/length(d.stimvol);

% run through the rois:
for ii = 1:length(d.stimvol)
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

