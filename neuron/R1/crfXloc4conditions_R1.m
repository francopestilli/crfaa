function d = crfXloc4conditions_R1(d,location)
% function d = crfXloc4conditions_R1(d,location)
%
% gets a d structure selects the stimvols for each
% stimulus condition and each location on the screen.
%
% does the analysis for the quadrant location selected (1,2,3,4)
%
% this function gets called as preprocessing:
% run = crfXlocation(d,1) [in the preprocessing field of the GUI]
%
% Conditions are:
% {Attended, Distributed, Unattended}x{PedContrast}x{Location1, Location2, Location3, Location4}
%
% franco pestilli 2007.09.28
%
%
% 1. get number of pedestal form the d
% 2. collect CRF for:
%    cue=1 location=target    (Attended)
%    cue=1 location=notTarget (Unattended)
%    cue=2 location=target    (Distributed Target)
%    cue=2 location=notTarget (Distributed not Target)
%

d.stimvol = {};
phaseNum = 2;
c = 0;

% get contrast used in the experiment
pedestals = d.stimfile{1}.stimulus.pedestals;
numPedestals = length(pedestals);
disp(sprintf('[crfXloc4conditions_R1: LOC%s] Pedestal contrasts: %s',num2str(location),num2str(pedestals)));

% set the number of HRF (length(pedestals) X 4 attentional conditions)
hdrNum = length(pedestals)*4;

% initialize matrix for each one of the 3 attentional conditions
% (Attended, Distributed, Unattended)
allstimvol = cell(1,hdrNum);

% setting to 0 the initial counter for checking number of volumes in each
% condition and the number of fixed volumes (not found and then set to vol=1)
for pedNum = 1:hdrNum;
    volXcond(pedNum) = 0;
    numFixed(pedNum) = 0;
end

for runNum = 1:length(d.stimfile)

    % getting the right parameters for this run
    myscreen{runNum} = d.stimfile{runNum}.myscreen;
    task{runNum} = d.stimfile{runNum}.task;
    e{runNum} = getTaskParameters(myscreen{runNum},task{runNum});
    
    % initialize variables for the new run (scan):
    stimvol = cell(1,hdrNum);

    % start saving the volumes in the right slot:
    for tNum = 1:e{runNum}(phaseNum).nTrials
        
        % getting all the variables out of the structure:
        cueCondition = e{runNum}(phaseNum).parameter.cueCondition(tNum);
        targetInterval = e{runNum}(phaseNum).randVars.targetInterval(tNum);
        targetLocation = e{runNum}(phaseNum).randVars.targetLocation(tNum);
        targetContrast = e{runNum}(phaseNum).traces(targetLocation,tNum);
        deltaC = e{runNum}(phaseNum).traces(5,tNum);
        
        %
        %%
        %%%
        %%%% SAVE ALL EVENTS
        c = c+1; % this counts over runs and trials
        % saving the deltaC contrast to plot crf shifted by 1/2 this contrast
        % saving the response as correct ('1') or incorrect ('0')
        % saving reaction time on each trial for each pedestal
        d.e.labels = {'1=attentionCondition' '2=targetInterval' '3=targetLocation'  '4=correctIncorrect' ...
            '5=buttonPressed' '6=reactionTime' '7=targetContrast' '8=deltaContrast'  '9=pedestalContrast' '10=runNum'};
        d.e.allEvents(1,c)  = e{runNum}(phaseNum).parameter.cueCondition(tNum);           % attentional condition 1 = attended, 2 = distributed
        d.e.allEvents(2,c)  = e{runNum}(phaseNum).randVars.targetInterval(tNum);          % target interval
        d.e.allEvents(3,c)  = e{runNum}(phaseNum).randVars.targetLocation(tNum);          % targetl location
        d.e.allEvents(4,c)  = (d.e.allEvents(2,c) == e{runNum}(phaseNum).response(tNum)); % correct (1) incorrect (0) response
        d.e.allEvents(5,c)  = e{runNum}(phaseNum).response(tNum);                         % button pressed
        d.e.allEvents(6,c)  = e{runNum}(phaseNum).reactionTime(tNum);                     % reaction time
        d.e.allEvents(7,c)  = e{runNum}(phaseNum).traces(d.e.allEvents(3,c),tNum);        % targetContrast
        d.e.allEvents(8,c)  = e{runNum}(phaseNum).traces(5,tNum);                         % deltaContrast
        d.e.allEvents(9,c)  = pedestals(e{runNum}(phaseNum).parameter.PedContrast(tNum)); % pedestal contrast
        d.e.allEvents(10,c) = runNum;                                                     % run number (scan number)
        d.e.allEvents(11,c) = e{runNum}(phaseNum).traces(1,tNum);                         % contrast at location 1
        d.e.allEvents(12,c) = e{runNum}(phaseNum).traces(2,tNum);                         % contrast at location 2
        d.e.allEvents(13,c) = e{runNum}(phaseNum).traces(3,tNum);                         % contrast at location 3
        d.e.allEvents(14,c) = e{runNum}(phaseNum).traces(4,tNum);                         % contrast at location 4
        %%%%
        %%%
        %%
        %
        
        %%% This volume is attended
        if (cueCondition == 1) && (targetLocation == location)
            % if attended then pedNum is PedContrast:
            pedNum  = e{runNum}(phaseNum).parameter.PedContrast(tNum);        
            conditionNum = pedNum;
            thisContrast = e{runNum}(phaseNum).traces(location,tNum)-deltaC;
%             disp(sprintf('[crfXloc4conditions_R1: LOC%s] Attended       - Scan: %s Trial: %s Contrast: %0.4f deltaC: %0.4f PedNum: %s hdrNum: %s',num2str(location),num2str(runNum),num2str(tNum),thisContrast,deltaC,num2str(pedNum),num2str(conditionNum)));

        %%% This volume is unattended
        elseif (cueCondition == 1) && (targetLocation ~= location)
            % if unattended then pedNum is the distracter contrast at this location is
            thisContrast = e{runNum}(phaseNum).traces(location,tNum);
            pedNum = find(thisContrast == pedestals);
            conditionNum = numPedestals+pedNum;

%           disp(sprintf('[crfXloc4conditions_R1: LOC%s] Unattended     - Scan: %s Trial: %s Contrast: %0.4f deltaC: %0.4f PedNum: %s hdrNum: %s',num2str(location),num2str(runNum),num2str(tNum),thisContrast,deltaC,num2str(pedNum),num2str(conditionNum)));

            %%% This volume is distributed and it is a target
        elseif (cueCondition == 2) && (targetLocation == location)

            thisContrast = e{runNum}(phaseNum).traces(location,tNum)-deltaC;
            pedNum = find(thisContrast == pedestals);            
            conditionNum = 2*numPedestals+pedNum;
%           disp(sprintf('[crfXloc4conditions_R1: LOC%s] Distributed  T - Scan: %s Trial: %s Contrast: %0.4f deltaC: %0.4f PedNum: %s hdrNum: %s',num2str(location),num2str(runNum),num2str(tNum),thisContrast,deltaC,num2str(pedNum),num2str(conditionNum)));

            %%% This volume is distributed BUT it is NOT a target
        elseif (cueCondition == 2) && (targetLocation ~= location)
            thisContrast = e{runNum}(phaseNum).traces(location,tNum);
            pedNum = find(thisContrast == pedestals);
            conditionNum = 3*numPedestals+pedNum;
%           disp(sprintf('[crfXloc4conditions_R1: LOC%s] Distributed NT - Scan: %s Trial: %s Contrast: %0.4f deltaC: %0.4f PedNum: %s hdrNum: %s',num2str(location),num2str(runNum),num2str(tNum),thisContrast,deltaC,num2str(pedNum),num2str(conditionNum)));
        else
            thisContrast = e{runNum}(phaseNum).traces(location,tNum);
            disp(sprintf('[crfXloc4conditions_R1: LOC%s] ERROR NO MATCH FOR THIS VOLUME! -Scan: %s Trial: %s Contrast: %0.4f deltaC: %0.4f PedNum: %s hdrNum: %s',num2str(location),num2str(runNum),num2str(tNum),thisContrast,deltaC,num2str(pedNum),num2str(conditionNum)));
            keyboard
        end
        
        % save volume
        stimvol{conditionNum} = [stimvol{conditionNum} e{runNum}(phaseNum).trialVolume(tNum)];
    end
    % get starting vol
    for attCond = 1:hdrNum
        stimvol{attCond} = [stimvol{attCond} + d.concatInfo.runTransition(runNum,1) - 1];
        % add the stimvols for this current scan into the rest from all the
        % scans:
        allstimvol{attCond} = [allstimvol{attCond} stimvol{attCond}];
    end
    disp(sprintf('[crfXloc4conditions_R1: LOC%i] volumes sorted out for SCAN%i',location,runNum));
end

% save the collected stimVols into the d structure for makescm to use it:
d.stimvol = allstimvol;
d.hdrNum = hdrNum;
