function d = focalVSdistributed(d)
% function d = focalVSdistributed(d)
%
% gets a d structure selects the stimvols for each
% cue condition and each location on the screen.
%
% this function gets called as preprocessing:
% run = crfXlocation(d,1) [in the preprocessing field of the GUI]
d.stimvol = {};
phaseNum = 2;

% set the number of HRF (= 2 cue conditions)
hdrNum = 2;

% initialize matrix for each one of the 3 attentional conditions
% (Attended, Distributed, Unattended)
allstimvol = cell(1,hdrNum);
for runNum = 1:length(d.stimfile)
    % getting the right parameters for this run
    myscreen = d.stimfile{runNum}.myscreen;
    task     = d.stimfile{runNum}.task;
    e{runNum}        = getTaskParameters(myscreen,task);
    
    % initialize variables for the new run (scan):
    stimvol = cell(1,hdrNum);

    % start saving the volumes in the right slot:
    for tNum = 1:e{runNum}(phaseNum).nTrials
        % getting all the variables out of the structure:
        cueCondition = e{runNum}(phaseNum).parameter.cueCondition(tNum);
            
        % save volume
        stimvol{cueCondition} = [stimvol{cueCondition} e{runNum}(phaseNum).trialVolume(tNum)];
    end
    % get starting vol
    for attCond = 1:2
        stimvol{attCond} = [stimvol{attCond} + d.concatInfo.runTransition(runNum,1) - 1];
        % add the stimvols for this current scan into the rest from all the
        % scans:
        allstimvol{attCond} = [allstimvol{attCond} stimvol{attCond}];
    end
    disp(sprintf('[focalVSdistributed] volumes sorted out for SCAN%i',runNum));
end

% save the collected stimVols into the d structure for makescm to use it:
d.stimvol = allstimvol;
d.hdrNum = hdrNum;
