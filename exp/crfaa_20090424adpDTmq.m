function crfaa_20090424adpDTmq(observer,adapter)
%  crfaa_20090424adpDTmq.m
%
%        $Id: crfaa_MMDDYY.m
%      usage: crfaa_YYYYMMDDadpDTmq(observer,adapter)
%         by: franco pestilli
%       date: 02/06/2007
%  copyright: (c) 2006 Justin Gardner (GPL see mgl/COPYING)
%    purpose: measure attended & unattended TvC functions at
%             different adaptation levels seglen = [10 30];
%             1.  does multiple quests.
%             2.  creates folder directory
%             3.  creates observer directoy
%             4.  DT does only detection test for adaptation
%      input: observer is a string
%             are created in
%             the current directory: currentdir/datafoldername/observer
%             (only two datafoldernames are alowed: data_psychophysics & data_scanner)
%             FP: i run this from the CRFAA folder
%
% ---- help end ----

global stimulus;
global MGL;

% check arguments
if ~any(nargin == [2])
 help  crfaa_20090424adpDTmq
 return
end

% check which on which computer this is running:
[status,whichcomputer] = system('hostname');
whichcomputer = deblank(whichcomputer(end-19:end));

% check if running psychophysiscs or scanner exp and on the right computers (jackson OR amacbookpro-4.local):
if (~status && (strcmp(whichcomputer,'jackson.cns.nyu.edu') || strcmp(whichcomputer,'amacbookpro-4.local')))

 if strcmp(whichcomputer,'jackson.cns.nyu.edu')
  datafoldername = 'data_psychophysics';
  psychoORscanner = 1;
  disp('[crfaa_20090424adpDTmq] ... running psychophysics on jackson.');
  disp('[crfaa_20090424adpDTmq] cd to /private/var/automount/i/7.2/frakkopesto/fmritools/ATT_ADPT_FILES/CRFAA/')
  cd('/private/var/automount/i/7.2/frakkopesto/fmritools/ATT_ADPT_FILES/CRFAA/');

 elseif strcmp(whichcomputer,'amacbookpro-4.local')
  datafoldername = 'data_psychophysics';
  psychoORscanner = 1;
  disp('[crfaa_20090424adpDTmq] ... running psychophysics ON amacbookpro-4.local.');
  disp('[crfaa_20090424adpDTmq] cd to /Users/frakkopesto/data/nyu/attadp')
  cd('/Users/frakkopesto/data/nyu/attadp/')

 else
  % only two data folder names allowed:
  disp(sprintf('[crfaa_20090424adpDTmq] Please choose datafoldername for psychophysics (data_psychophysics) OR for the scanner (data_scanner). %s', datafoldername));
  help  crfaa_20090424adpDTmq
  return
 end

else
 disp(sprintf('[crfaa_20090424adpDTmq] command "system" returned status "%i" on the call: "[status,whichcomputer] = system(''hostname'');"',status));
 disp(sprintf('[crfaa_20090424adpDTmq] command "system" returned computer "%s" on the call: "[status,whichcomputer] = system(''hostname'');"',whichcomputer));
 disp('[crfaa_20090424adpDTmq] exiting if status ~= 0 (system not checked) OR if whichcomputer not "jackson.cns.nyu.edu" OR "stimulus-g5"');
 return
end

% make a data directory if necessary
datadirname = fullfile(deblank(pwd),datafoldername);
if ~isdir(datadirname)
 disp(sprintf('[crfaa_20090424adpDTmq] Data dir "%s" not found, making data directory.', datadirname));
 mkdir(datadirname);
else
 disp(sprintf('[crfaa_20090424adpDTmq] Found data directory "%s".', datadirname));
end

% make an observer directory if necessary
datadirname = fullfile(datadirname,observer);
if ~isdir(datadirname);
 disp(sprintf('[crfaa_20090424adpDTmq] Observer dir "%s" not found, making observer directory.', datadirname));
 mkdir(datadirname);
else
 disp(sprintf('[crfaa_20090424adpDTmq] Found observer directory "%s".', datadirname));

end

% look for old stim files
stimfileindatadir = dir(fullfile(datadirname,'*stim*.mat'));
if ~isempty(stimfileindatadir) % there is a stimfile
 % load last file in directory
 load(fullfile(datadirname,stimfileindatadir(end).name));
 disp(sprintf('[crfaa_20090424adpDTmq] Previous quest data found, loading quest parameters from: %s',stimfileindatadir(end).name));
 disp(sprintf('[crfaa_20090424adpDTmq] ... saving data in: %s',datadirname));
else % there is no file yet
 disp('[crfaa_20090424adpDTmq] No quest data found, will initialize a new quest.');
 disp(sprintf('[crfaa_20090424adpDTmq] ... saving data in: %s',datadirname));
end

% ereasing everything in stimulus except the quest structure
if isfield(stimulus,'quest')
 quest = stimulus.quest;
 stimulus = [];
 stimulus.quest = quest;
else
 stimulus = [];
end

% clearing old variables:
clear task myscreen;

% initalize the screen
if psychoORscanner == 1
 if strcmp(whichcomputer,'jackson.cns.nyu.edu')
  myscreen.calibFilename = '~/fmritools/ATT_ADPT_FILES/CRFAA/mgl/task/displays/0002_jackson_080205.mat';
  if isfile(myscreen.calibFilename)
   myscreen.screenParams{1} = {'jackson',[],2,1024,768,57,[40 30],75,1,0,1.8,myscreen.calibFilename,[0 0]};
   disp(sprintf('[crfaa_20090424adpDTmq] found suggested calibration file (%s)',myscreen.calibFilename));
  else
   disp(sprintf('[crfaa_20090424adpDTmq] suggested calibration file NOT found (%s)',myscreen.calibFilename));
   return
  end
 elseif strcmp(whichcomputer,'amacbookpro-4.local')
  myscreen.calibFilename = '~/code/mgl/task/displays/0002_amacbook-3_081216';
  if isfile([myscreen.calibFilename,'.mat'])
   disp(sprintf('[crfaa_20090424adpDTmq] found suggested calibration file (%s)',myscreen.calibFilename));
  else
   disp(sprintf('[crfaa_20090424adpDTmq] suggested calibration file NOT found (%s)',myscreen.calibFilename));
   return
  end
 end
else
 disp(sprintf('[crfaa_20090424adpDTmq] suggested calibration file NOT found (%s)',myscreen.calibFilename));
 return
end

% initalize the screen
myscreen = initScreen(myscreen);

%% saving some useul stuff:
myscreen.allowpause = 0;
myscreen.saveData = 1;
myscreen.datadir = datadirname;
myscreen.whichComputer = whichcomputer;
myscreen.whichFunction = which('crfaa_20090424adpDTmq');
myscreen.stimfileindatadir = datadirname;
% myscreen.calibrationFile = myscreen.calibFilename;

%% initialize tasks strucure:
% task 1 = initial blank delay + adaptation period
task{1}.numTrials = 1;
task{1}.waitForBacktick = 1;
task{1}.seglen = [2 70];

% real experiment:
% blanks between adp and stim:
task{2}.numTrials = 80;
task{2}.segmin =      [1 .1 .6 .2 .6 .1 1.3 1];%[2 .4 .6 .2 .6 .4 1.3 4.5];
task{2}.segmax =      [1 .1 .6 .2 .6 .1 1.3 3];%[2 .4 .6 .2 .6 .4 1.3 4.5];
task{2}.segquant =    [0  0  0  0  0  0  0  1];
% wait for back tick
if psychoORscanner
 % do not wait for backtick in psychophysiscs
 task{2}.synchToVol =  [ 0  0  0  0  0  0  0   0];
else
 % wait for backtick on the scanner
 task{2}.synchToVol =  [ 0  0  0  0  0  0  0   1];
end
task{2}.getResponse = [ 0  0  0  0  0  0  1   0];
task{2}.parameter.PedContrast = 1;
task{2}.parameter.cueCondition = 1:2;
task{2}.randVars.uniform.targetInterval = 1:2;
task{2}.randVars.uniform.targetLocation = 1:4;
task{2}.random = 1;

% initialize the task
[task{1} myscreen] = initTask(task{1},myscreen,@adaptationStartSegmentCallback,@adaptationDrawStimulusCallback);
[task{2} myscreen] = initTask(task{2},myscreen,@stimStartSegmentCallback,@stimDrawStimulusCallback,@responseCallback);

%% global initialization of the stimulus
myscreen = initStimulus('stimulus',myscreen);

%% MY stimulus inititialization:
stimulus.observer  = observer;

% adaptation test:
stimulus.pedestals = [0];
% stimulus.distractorContrast = {[0 0 0]}; % NB change distracor range also... 
stimulus.distractorContrast = {[.75 1 1.25]}; % ration around the delta-c
% contrast
stimulus.distractorCrange = 0.3;               % normally distributed variability around the distractor 
                                              % contrast that will be generated at each trial
                                              
stimulus.minContrast = 1/508;
stimulus.maxContrast = [.2];
stimulus.adaptationContrast = adapter; % TEST ADAptation C: .0001 .0175 .625 .125 .25 .5 OR 1
stimulus.flickerFrequency = 5; % in Hz

% sound effects (intervals and feedbacks)
stimulus.IntervalSound = find(strcmp(MGL.soundNames,'Morse'));
stimulus.rightResponseSound = find(strcmp(MGL.soundNames,'Purr'));
stimulus.wrongResponseSound = find(strcmp(MGL.soundNames,'Basso'));
stimulus.tooLateResponseSound = find(strcmp(MGL.soundNames,'Pong'));
stimulus.tmp = [];

%%%%% if the QUEST was not taken froma previous file then set its parameters:
if ~isfield(stimulus,'quest')
 stimulus.quest.length = 40;

 stimulus.quest.tGuess     = log10([.05; .05]);
 stimulus.quest.tGuessSd   =   3;
 stimulus.quest.pThreshold = .75;
 stimulus.quest.beta       =   3;
 stimulus.quest.delta      = .05;
 stimulus.quest.gamma      =  .5;

 for jj = 1:2
  for ii = 1:length(stimulus.pedestals)
   stimulus.quest.questnum(ii,jj) = 1;
   stimulus.quest.q{ii,jj,1} = QuestCreate(stimulus.quest.tGuess(jj,ii),stimulus.quest.tGuessSd, ...
    stimulus.quest.pThreshold,stimulus.quest.beta,stimulus.quest.delta,stimulus.quest.gamma);
  end
 end
end
myInitStimulus;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% run the eye calibration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% myscreen = eyeCalibDisp(myscreen);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main display loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
phaseNum = 1;
while (phaseNum <= length(task)) && ~myscreen.userHitEsc
 % update the task
 [task myscreen phaseNum] = updateTask(task,myscreen,phaseNum);
 % flip screen
 [myscreen task] = tickScreen(myscreen,task);
end
clear stimulus.tmp

% if we got here, we are at the end of the experiment
myscreen = endTask(myscreen,task);

% show the results:
results = crfaa_analquest(observer)


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TASK 1: function that gets called at the start of each segment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = adaptationStartSegmentCallback(task, myscreen)
global stimulus;
mglClearScreen(stimulus.grayColor);
setGammaTable(stimulus.adaptationContrast);


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TASK 1: function that gets called to draw the stimulus each frame
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = adaptationDrawStimulusCallback(task, myscreen)
global stimulus;
mglClearScreen(stimulus.grayColor);

% putting up the frame:
for ii = 1:size(stimulus.locations,2)
 mglFillOval(stimulus.eccentricity*stimulus.locations{ii}(1), ...
  stimulus.eccentricity*stimulus.locations{ii}(2),[stimulus.height+stimulus.frameThick stimulus.width+stimulus.frameThick],stimulus.black);
 mglFillOval(stimulus.eccentricity*stimulus.locations{ii}(1), ...
  stimulus.eccentricity*stimulus.locations{ii}(2),[stimulus.height stimulus.width],128);
end

if (task.thistrial.thisseg == 1)
 % blank screen at the beginning of each scan
 mglFixationCross(1,1,stimulus.black);

elseif (task.thistrial.thisseg == 2)
 % adaptation stimulation
 phaseNum = mod(floor(2*mglGetSecs(task.thistrial.trialstart)/stimulus.flickerPeriod),2)+1;
 drawGabor(stimulus.adaptationContrast,stimulus.eccentricity*[stimulus.locations{1}],phaseNum);
 drawGabor(stimulus.adaptationContrast,stimulus.eccentricity*[stimulus.locations{2}],phaseNum);
 drawGabor(stimulus.adaptationContrast,stimulus.eccentricity*[stimulus.locations{3}],phaseNum);
 drawGabor(stimulus.adaptationContrast,stimulus.eccentricity*[stimulus.locations{4}],phaseNum);
 mglFixationCross(1,1,stimulus.black);
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TASK 2: function that gets called at the start of each segment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = stimStartSegmentCallback(task, myscreen)

global stimulus;
mglClearScreen(stimulus.grayColor);
if (task.thistrial.thisseg == 1) % precueing

 %     disp(sprintf('PRECUE'));
 mglClearScreen(stimulus.grayColor);
 %set gamma table to adaptation contrast
 setGammaTable(stimulus.adaptationContrast);
 mglFixationCross(1,1,stimulus.black);

 % order stimuli locations (target is first in vector):
 stimulus.tmp.stimLocationindx = (1:4) == task.thistrial.targetLocation;
 stimulus.tmp.stimLocationindx = [task.thistrial.targetLocation find(~stimulus.tmp.stimLocationindx)];

 % set precue and response cue colors:
 stimulus.tmp.precueColor = {stimulus.background stimulus.background stimulus.background stimulus.background};
 stimulus.tmp.respcueColor = stimulus.tmp.precueColor;
 stimulus.tmp.respcueColor{task.thistrial.targetLocation} = stimulus.green;

 % set up target for first or second interval:
 stimulus.tmp.targetIntindx = (1:2) == task.thistrial.targetInterval;
 stimulus.tmp.targetIntindx = [task.thistrial.targetInterval find(~stimulus.tmp.targetIntindx)];

 % set pedestal and pedestal+delta contrast (target):
 % choose contrast with quest here:
 stimulus.tmp.deltaC = 10^QuestQuantile(stimulus.quest.q{task.thistrial.PedContrast, task.thistrial.cueCondition, ...
  stimulus.quest.questnum(task.thistrial.PedContrast, task.thistrial.cueCondition)});

 % SET CONTRAST TO DISPLAY MIN INCREMENT RESOLUTION (1/512-reserved_colors )
 stimulus.tmp.deltaC = max(stimulus.tmp.deltaC,stimulus.minContrast);

 stimulus.tmp.PedCdeltaC = stimulus.pedestals(task.thistrial.PedContrast)+stimulus.tmp.deltaC;

 % SET CONTRAST TO DISPLAY MAX: max contrast is given as parameter at the beginning for each pedestal contrast
 stimulus.tmp.PedCdeltaC = min(stimulus.tmp.PedCdeltaC,stimulus.maxContrast(task.thistrial.PedContrast));

 % check which contrast was displayed at the stimulus:
 stimulus.tmp.deltaC = stimulus.tmp.PedCdeltaC-stimulus.pedestals(task.thistrial.PedContrast);

 stimulus.tmp.targetContrasts(stimulus.tmp.targetIntindx(1)) = stimulus.tmp.PedCdeltaC;
 stimulus.tmp.targetContrasts(stimulus.tmp.targetIntindx(2)) = stimulus.pedestals(task.thistrial.PedContrast);

 % put the target in the right position in the vector of all contrasts:
 stimulus.tmp.allContrasts(task.thistrial.targetLocation) = stimulus.tmp.targetContrasts(stimulus.tmp.targetIntindx(1));

 % The next few lines choose distracter contrast
 % (1) take the current c+delta-c contrast
 % (2) have one distracter higher that that but less then macContrast
 % (3) have one distracter nearly half that contrast OR zero
 % (4) have the third distracter lower than delta c
 stimulus.tmp.distracterindx = setdiff(1:4,task.thistrial.targetLocation);

 % shuffle distacters indexes around to move them across locations:
 randDindexes =  randperm(length(stimulus.distractorContrast{task.thistrial.PedContrast}));

 for i = 1:3
  dindex = randDindexes(i);
 % this line creates distracters contrasts as a proportion (above and below) of deltaC
 % they clip the distracter contrast either at 0 or at max  
  stimulus.tmp.allContrasts(stimulus.tmp.distracterindx(i)) = ...
   max([0, ...
   min([stimulus.maxContrast(task.thistrial.PedContrast), ...
   stimulus.tmp.PedCdeltaC*stimulus.distractorContrast{task.thistrial.PedContrast}(dindex) + randn * stimulus.distractorCrange])]);
 end
 

 % set distracter contrast for display:
 stimulus.tmp.dispDistrContrast = stimulus.tmp.allContrasts(stimulus.tmp.distracterindx);
 % disp(sprintf('Distracters contras per location %0.2f %0.2f %0.2f', stimulus.tmp.dispDistrContrast));

 % setting the maximum of the clut contrast to the max of the contrast
 % to display:
 stimulus.tmp.clutcontrast = max(stimulus.tmp.allContrasts);

 % save a trace for each of the 4 possible distracters locations (-1 is for the target location)
 for ii = 1:4
  myscreen = writeTrace(stimulus.tmp.allContrasts(ii),myscreen.stimtrace+ii-1,myscreen,1);
 end
 myscreen = writeTrace(stimulus.tmp.deltaC,myscreen.stimtrace+4,myscreen,1);
 myscreen = writeTrace(task.thistrial.targetLocation,myscreen.stimtrace+5,myscreen,1);
 myscreen = writeTrace(task.thistrial.targetInterval,myscreen.stimtrace+6,myscreen,1);

 % choose attentional condition for precueing:
 if (task.thistrial.cueCondition == 1) % attended
  stimulus.tmp.precueColor{task.thistrial.targetLocation} = stimulus.white;

 elseif (task.thistrial.cueCondition == 2) % distributed
  stimulus.tmp.precueColor={stimulus.white stimulus.white stimulus.white stimulus.white};
 end
 stimulus.tmp.responseCueColor = stimulus.black;

 disp(sprintf('%i: CueCond: %i targetLoc: %i targetInt: %i ped#: %i pedC: %0.3f; targetC: %0.4f; DistC: %0.3f; %0.3f; %0.3f;',...
  task.trialnumTotal,...
  task.thistrial.cueCondition,...
  task.thistrial.targetLocation, ...
  task.thistrial.targetInterval,...
  task.thistrial.PedContrast,...
  stimulus.tmp.targetContrasts(stimulus.tmp.targetIntindx(2)),...
  stimulus.tmp.targetContrasts(stimulus.tmp.targetIntindx(1)),...
  stimulus.tmp.dispDistrContrast));

elseif (task.thistrial.thisseg == 2) % delay1
 mglClearScreen(stimulus.grayColor);
 setGammaTable(.000001);
 mglFixationCross(1,1,stimulus.black);

elseif (task.thistrial.thisseg == 3) % stimulus presentation 1st interval
 %     disp(sprintf('STIMULUS1'));
 mglClearScreen(stimulus.grayColor);
 setGammaTable(stimulus.tmp.clutcontrast);
 mglPlaySound(stimulus.IntervalSound);
 mglFixationCross(1,1,stimulus.black);

elseif (task.thistrial.thisseg == 4) % ISI
 %     disp(sprintf('ISI'));
 mglClearScreen(stimulus.grayColor);
 mglFixationCross(1,1,stimulus.black);

elseif (task.thistrial.thisseg == 5) % stimulus presentation 2nd interval
 %     disp(sprintf('STIMULUS2'));
 mglClearScreen(stimulus.grayColor);
 setGammaTable(stimulus.tmp.clutcontrast);
 mglPlaySound(stimulus.IntervalSound);

elseif (task.thistrial.thisseg == 6) % delay2
 %     disp(sprintf('DELAY'));
 mglClearScreen(stimulus.grayColor);
 setGammaTable(0.0001);
 mglFixationCross(1,1,stimulus.black);

elseif (task.thistrial.thisseg == 7) % response interval
 %     disp(sprintf('RESPONSE'));
 mglClearScreen(stimulus.grayColor);
 setGammaTable(stimulus.adaptationContrast);
 mglFixationCross(1,1,stimulus.black);

elseif (task.thistrial.thisseg == 8) % top-up (ITI)
 %     disp(sprintf('TOPUP'));
 %set gamma table to adaptation contrast
 mglClearScreen(stimulus.grayColor);
 setGammaTable(stimulus.adaptationContrast);
 mglFixationCross(1,1,stimulus.black);
 % disp(sprintf('Contrast = %i',task.thistrial.contrast));
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TASK 2: function that gets called to draw the stimulus each frame
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = stimDrawStimulusCallback(task, myscreen)

global stimulus
mglClearScreen(stimulus.grayColor);

% putting up the frame:
for ii = 1:size(stimulus.locations,2)
 mglFillOval(stimulus.eccentricity*stimulus.locations{ii}(1), ...
  stimulus.eccentricity*stimulus.locations{ii}(2),[((stimulus.height))+stimulus.frameThick ((stimulus.width))+stimulus.frameThick],stimulus.black);
 mglFillOval(stimulus.eccentricity*stimulus.locations{ii}(1), ...
  stimulus.eccentricity*stimulus.locations{ii}(2),[((stimulus.height)) ((stimulus.width))],128);
end

% sequence of events in a trial:
if (task.thistrial.thisseg == 1) % precueing
 % top-up stimulation
 phaseNum = mod(floor(2*mglGetSecs(task.thistrial.trialstart)/stimulus.flickerPeriod),2)+1;
 drawGabor(stimulus.adaptationContrast,stimulus.eccentricity*[stimulus.locations{1}],phaseNum);
 drawGabor(stimulus.adaptationContrast,stimulus.eccentricity*[stimulus.locations{2}],phaseNum);
 drawGabor(stimulus.adaptationContrast,stimulus.eccentricity*[stimulus.locations{3}],phaseNum);
 drawGabor(stimulus.adaptationContrast,stimulus.eccentricity*[stimulus.locations{4}],phaseNum);

 % draw cues
 mglQuad(stimulus.arrowLocationX{1},stimulus.arrowLocationY{1}, stimulus.tmp.precueColor{1}', 0);
 mglQuad(stimulus.arrowLocationX{2},stimulus.arrowLocationY{2}, stimulus.tmp.precueColor{2}', 0);
 mglQuad(stimulus.arrowLocationX{3},stimulus.arrowLocationY{3}, stimulus.tmp.precueColor{3}', 0);
 mglQuad(stimulus.arrowLocationX{4},stimulus.arrowLocationY{4}, stimulus.tmp.precueColor{4}', 0);
 mglFixationCross(1,1,stimulus.black);

elseif (task.thistrial.thisseg == 2) % Delay1 (NO top-up)
 mglQuad(stimulus.arrowLocationX{1},stimulus.arrowLocationY{1}, stimulus.tmp.precueColor{1}', 0);
 mglQuad(stimulus.arrowLocationX{2},stimulus.arrowLocationY{2}, stimulus.tmp.precueColor{2}', 0);
 mglQuad(stimulus.arrowLocationX{3},stimulus.arrowLocationY{3}, stimulus.tmp.precueColor{3}', 0);
 mglQuad(stimulus.arrowLocationX{4},stimulus.arrowLocationY{4}, stimulus.tmp.precueColor{4}', 0);
 mglFixationCross(1,1,stimulus.black);

elseif (task.thistrial.thisseg == 3) % first  interval
 % draw cues
 mglQuad(stimulus.arrowLocationX{1},stimulus.arrowLocationY{1}, stimulus.tmp.precueColor{1}', 0);
 mglQuad(stimulus.arrowLocationX{2},stimulus.arrowLocationY{2}, stimulus.tmp.precueColor{2}', 0);
 mglQuad(stimulus.arrowLocationX{3},stimulus.arrowLocationY{3}, stimulus.tmp.precueColor{3}', 0);
 mglQuad(stimulus.arrowLocationX{4},stimulus.arrowLocationY{4}, stimulus.tmp.precueColor{4}', 0);

 % Display stimulus and change gabor's phase
 phaseNum = mod(floor(2*mglGetSecs(task.thistrial.trialstart)/stimulus.flickerPeriod),2)+1;
 drawGabor(stimulus.tmp.targetContrasts(1),stimulus.eccentricity*[stimulus.locations{stimulus.tmp.stimLocationindx(1)}],phaseNum);
 drawGabor(stimulus.tmp.dispDistrContrast(1),stimulus.eccentricity*[stimulus.locations{stimulus.tmp.stimLocationindx(2)}],phaseNum);
 drawGabor(stimulus.tmp.dispDistrContrast(2),stimulus.eccentricity*[stimulus.locations{stimulus.tmp.stimLocationindx(3)}],phaseNum);
 drawGabor(stimulus.tmp.dispDistrContrast(3),stimulus.eccentricity*[stimulus.locations{stimulus.tmp.stimLocationindx(4)}],phaseNum);
 mglFixationCross(1,1,stimulus.white);

elseif (task.thistrial.thisseg == 4) % inter-stimulus-interval
 % draw cues
 mglQuad(stimulus.arrowLocationX{1},stimulus.arrowLocationY{1}, stimulus.tmp.precueColor{1}', 0);
 mglQuad(stimulus.arrowLocationX{2},stimulus.arrowLocationY{2}, stimulus.tmp.precueColor{2}', 0);
 mglQuad(stimulus.arrowLocationX{3},stimulus.arrowLocationY{3}, stimulus.tmp.precueColor{3}', 0);
 mglQuad(stimulus.arrowLocationX{4},stimulus.arrowLocationY{4}, stimulus.tmp.precueColor{4}', 0);
 mglFixationCross(1,1,stimulus.black);

elseif (task.thistrial.thisseg == 5) % second interval
 % draw cues
 mglQuad(stimulus.arrowLocationX{1},stimulus.arrowLocationY{1}, stimulus.tmp.precueColor{1}', 0);
 mglQuad(stimulus.arrowLocationX{2},stimulus.arrowLocationY{2}, stimulus.tmp.precueColor{2}', 0);
 mglQuad(stimulus.arrowLocationX{3},stimulus.arrowLocationY{3}, stimulus.tmp.precueColor{3}', 0);
 mglQuad(stimulus.arrowLocationX{4},stimulus.arrowLocationY{4}, stimulus.tmp.precueColor{4}', 0);

 % display stimulus and change gabor's phase
 phaseNum = mod(floor(2*mglGetSecs(task.thistrial.trialstart)/stimulus.flickerPeriod),2)+1;
 drawGabor(stimulus.tmp.targetContrasts(2),stimulus.eccentricity*[stimulus.locations{stimulus.tmp.stimLocationindx(1)}],phaseNum);
 drawGabor(stimulus.tmp.dispDistrContrast(1),stimulus.eccentricity*[stimulus.locations{stimulus.tmp.stimLocationindx(2)}],phaseNum);
 drawGabor(stimulus.tmp.dispDistrContrast(2),stimulus.eccentricity*[stimulus.locations{stimulus.tmp.stimLocationindx(3)}],phaseNum);
 drawGabor(stimulus.tmp.dispDistrContrast(3),stimulus.eccentricity*[stimulus.locations{stimulus.tmp.stimLocationindx(4)}],phaseNum);
 mglFixationCross(1,1,stimulus.white);

elseif (task.thistrial.thisseg == 6) % Delay2 (NO top-up)
 mglQuad(stimulus.arrowLocationX{1},stimulus.arrowLocationY{1}, stimulus.tmp.precueColor{1}', 0);
 mglQuad(stimulus.arrowLocationX{2},stimulus.arrowLocationY{2}, stimulus.tmp.precueColor{2}', 0);
 mglQuad(stimulus.arrowLocationX{3},stimulus.arrowLocationY{3}, stimulus.tmp.precueColor{3}', 0);
 mglQuad(stimulus.arrowLocationX{4},stimulus.arrowLocationY{4}, stimulus.tmp.precueColor{4}', 0);
 mglFixationCross(1,1,stimulus.black);

elseif (task.thistrial.thisseg == 7) % get response
 % top-up stimulation
 phaseNum = mod(floor(2*mglGetSecs(task.thistrial.trialstart)/stimulus.flickerPeriod),2)+1;
 drawGabor(stimulus.adaptationContrast,stimulus.eccentricity*[stimulus.locations{1}],phaseNum);
 drawGabor(stimulus.adaptationContrast,stimulus.eccentricity*[stimulus.locations{2}],phaseNum);
 drawGabor(stimulus.adaptationContrast,stimulus.eccentricity*[stimulus.locations{3}],phaseNum);
 drawGabor(stimulus.adaptationContrast,stimulus.eccentricity*[stimulus.locations{4}],phaseNum);
 % draw cues
 mglQuad(stimulus.arrowLocationX{1},stimulus.arrowLocationY{1}, stimulus.tmp.respcueColor{1}', 0);
 mglQuad(stimulus.arrowLocationX{2},stimulus.arrowLocationY{2}, stimulus.tmp.respcueColor{2}', 0);
 mglQuad(stimulus.arrowLocationX{3},stimulus.arrowLocationY{3}, stimulus.tmp.respcueColor{3}', 0);
 mglQuad(stimulus.arrowLocationX{4},stimulus.arrowLocationY{4}, stimulus.tmp.respcueColor{4}', 0);
 mglFixationCross(1,1,stimulus.tmp.responseCueColor);

elseif (task.thistrial.thisseg == 8) % inter-trial-interval (top-up)
 % top-up stimulation
 phaseNum = mod(floor(2*mglGetSecs(task.thistrial.trialstart)/stimulus.flickerPeriod),2)+1;
 drawGabor(stimulus.adaptationContrast,stimulus.eccentricity*[stimulus.locations{1}],phaseNum);
 drawGabor(stimulus.adaptationContrast,stimulus.eccentricity*[stimulus.locations{2}],phaseNum);
 drawGabor(stimulus.adaptationContrast,stimulus.eccentricity*[stimulus.locations{3}],phaseNum);
 drawGabor(stimulus.adaptationContrast,stimulus.eccentricity*[stimulus.locations{4}],phaseNum);
 mglFixationCross(1,1,stimulus.black);

end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to deal with subject responses
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = responseCallback(task, myscreen)

global stimulus
mglClearScreen(stimulus.grayColor);

whichButton = task.thistrial.whichButton;
if task.thistrial.gotResponse == 0

 % check response correct or not
 stimulus.tmp.response = whichButton == task.thistrial.targetInterval;

 % get which quest run we are on
 questnum = stimulus.quest.questnum(task.thistrial.PedContrast, task.thistrial.cueCondition);

 % update quest here after response
 stimulus.quest.q{task.thistrial.PedContrast, task.thistrial.cueCondition, questnum} = ...
  QuestUpdate(stimulus.quest.q{task.thistrial.PedContrast, task.thistrial.cueCondition, questnum}, ...
  log10(stimulus.tmp.deltaC),stimulus.tmp.response);
 %     disp(sprintf('Subject pressed %i',whichButton));

 % change fixaton point color to give feeback:
 if stimulus.tmp.response
  stimulus.tmp.responseCueColor = stimulus.green;
  mglPlaySound(stimulus.rightResponseSound);
  %         disp(sprintf('Response is correct'));
 else
  stimulus.tmp.responseCueColor = stimulus.red;
  mglPlaySound(stimulus.wrongResponseSound);
  %         disp(sprintf('Response is incorrect'));
 end
 %length(stimulus.quest{end}.q{task.thistrial.PedContrast, task.thistrial.cueCondition}.intensity) > stimulus.quest{1}.length

 % restart quest if number of trials exceeded:
 %     disp(sprintf('this Quest length %i',length(stimulus.quest.q{task.thistrial.PedContrast, task.thistrial.cueCondition, questnum}.intensity)));
 if length(stimulus.quest.q{task.thistrial.PedContrast, task.thistrial.cueCondition, questnum}.intensity) > stimulus.quest.length

  % increment the questnum
  questnum = questnum+1;
  % change it for this type
  stimulus.quest.questnum(task.thistrial.PedContrast, task.thistrial.cueCondition) = questnum;
  % and reinit the quest
  %         starttime = mglGetSecs;

  q = QuestRecompute(stimulus.quest.q{task.thistrial.PedContrast, task.thistrial.cueCondition,questnum-1});

  stimulus.quest.q{task.thistrial.PedContrast, task.thistrial.cueCondition,questnum} = ...
   QuestCreate(stimulus.quest.tGuess(task.thistrial.cueCondition,task.thistrial.PedContrast),stimulus.quest.tGuessSd, ...
   q.pThreshold,stimulus.quest.beta,stimulus.quest.delta,stimulus.quest.gamma);

  %         lentime = 1000*mglGetSecs(starttime);
  %         disp(sprintf('QuestCreate took %0.1f (ms)',lentime));
 end
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to init the stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function myInitStimulus

global stimulus;

% let's get the linearized gamma table
stimulus.linearizedGammaTable = mglGetGammaTable;

if ~isfield(stimulus,'init') || 1
 % stimulus parameters:
 stimulus.init = 1;
 stimulus.width = 8;          % in deg
 stimulus.height = 8;         % in deg
 stimulus.sf = 2;             % in cpd
 stimulus.orientation = 90;   % in deg
 stimulus.phase = [0 180];    % in deg
 stimulus.eccentricity = 5.2; % in deg

 % frames and locations:
 % frames and locations:
 % frames and locations:
 stimulus.locations = {[1 1], [1 -1], [-1 -1], [-1 1]};
 stimulus.frameThick = .15;

 % cues:
 % cues:
 % cues:
 stimulus.arrowSize = .4;       % in deg (the arrow is sqrt(3) longer then this value)
 stimulus.cueCondition = [1 2]; % attended (1) & distributed (2)
 stimulus.arrowLocationX = {    % cue X coorinates
  stimulus.arrowSize*[ 1;  sqrt(3); 0;   1], ...
  stimulus.arrowSize*[ 1;  sqrt(3); 0;   1], ...
  stimulus.arrowSize*[-1; -sqrt(3); 0;  -1], ...
  stimulus.arrowSize*[-1; -sqrt(3); 0;  -1]};

 stimulus.arrowLocationY = {    % cue Y coorinates
  stimulus.arrowSize*[0;  sqrt(3);  1;  1],...
  stimulus.arrowSize*[0; -sqrt(3); -1; -1],...
  stimulus.arrowSize*[0; -sqrt(3); -1; -1],...
  stimulus.arrowSize*[0;  sqrt(3);  1;  1]};

 % stimulus size
 stimulus.gaussSdx = stimulus.width/7;
 stimulus.gaussSdy = stimulus.height/7;

 % reserved colors (feedback and such)
 stimulus.reservedColors = [0 0 0; 1 1 1; 0 .65 0; 1 0 0];
end

%% some computations:
% check that maxContrast is effectively the max contrast:
% for ii = 1:length(stimulus.pedestals)
%  for jj = 1:size(stimulus.distractorContrast{ii},2)
%   if (stimulus.distractorContrast{ii}(jj) > stimulus.maxContrast(ii))
%    disp(sprintf('Distracter contrast %i (%.2f) for pedestal %i (%.2f) is greater then maxContrast (%.2f) for the same pedestal', ...
%     jj,stimulus.distractorContrast{ii}(jj),ii,stimulus.pedestals(ii),stimulus.maxContrast(ii)));
%    keyboard;
%   end
%  end
% end
% for ii = 1:length(stimulus.pedestals)
%  if (stimulus.pedestals(ii) > stimulus.maxContrast(ii))
%   disp(sprintf('Pedestal contrast %i (%.2f) is greater then maxContrast (%.2f) for the same pedestal', ...
%    ii,stimulus.pedestals(ii),stimulus.maxContrast(ii)));
%   keyboard;
%  end
% end

% for stimulus flicker frequency
stimulus.flickerPeriod = 1/stimulus.flickerFrequency;
stimulus.tmp.targetContrasts=zeros(1,2); % blank vect initialization to decide first or second interval later

% for contrast
stimulus.nPedestals=length(stimulus.pedestals);

stimulus.nReservedColors=size(stimulus.reservedColors,1);
stimulus.nGratingColors = 256-(2*floor(stimulus.nReservedColors/2)+1);
stimulus.minGratingColors = 2*floor(stimulus.nReservedColors/2)+1;
stimulus.midGratingColors = stimulus.minGratingColors+floor(stimulus.nGratingColors/2);
stimulus.maxGratingColors = 255;
stimulus.deltaGratingColors = floor(stimulus.nGratingColors/2);

% to set up color values
stimulus.black = [0 0 0];
stimulus.white = [1/255 1/255 1/255];
stimulus.background = [ceil((255-stimulus.nReservedColors)/2)/255 ceil((255-stimulus.nReservedColors)/2)/255 ceil((255-stimulus.nReservedColors)/2)/255];
stimulus.green = [2/255 2/255 2/255];
stimulus.red = [3/255 3/255 3/255];

% calculate a grating and gaussian envelope (gaussian is in the alpha channel)
for thisPhase = 1:length(stimulus.phase)
 gratingMatrix{thisPhase} = makeGrating(stimulus.width,stimulus.height,stimulus.sf,stimulus.orientation,stimulus.phase(thisPhase));
 grating{thisPhase}(:,:,4) = 255*(makeGaussian(stimulus.width,stimulus.height,stimulus.gaussSdx,stimulus.gaussSdy) > exp(-14/2));
end

% making the texture for all the stimuli:
disppercent(-inf,'Calculating gratings');
for thisPhase = 1:length(stimulus.phase)
 for thisContrast = 0:stimulus.deltaGratingColors
  %% stimulus.texture
  grating{thisPhase}(:,:,1) = stimulus.midGratingColors+gratingMatrix{thisPhase}*thisContrast;
  grating{thisPhase}(:,:,2) = grating{thisPhase}(:,:,1);
  grating{thisPhase}(:,:,3) = grating{thisPhase}(:,:,1);
  stimulus.tex{thisPhase}(thisContrast+1) = mglCreateTexture(grating{thisPhase});
  disppercent(thisContrast/stimulus.deltaGratingColors);
 end
end
disppercent(inf);
stimulus.nDisplayContrasts = stimulus.deltaGratingColors;

% calculate gray color
stimulus.grayColor = stimulus.midGratingColors/255;
% set background to gray for eye calibration:
mglClearScreen(stimulus.grayColor);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to draw the gabor stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function drawGabor(desiredContrast,position,phase)
% drawGabor
%
%        $Id: drawGabor.m, v 1 2007/01/18 19:40:56 ?? ?? $
%      usage: drawGabor(desiredContrast,position,pedestalNum,phase)
%         by: justin gardner, modified by franco pestilli
%       date: 01/11/07 - 01/18/07
%  copyright: (c) 2007 Justin Gardner
%    purpose: draw a gabor stimulus on the screen with a specified contrast
%             within a given clut range (it finds the closest contrast to
%             the requested one within the available clut range)
%
% it makes two operations:
% displayContrastNum = min(round(stimulus.nDisplayContrasts*desiredContrast/stimulus.currentMaxContrast)...
%                      ,stimulus.nDisplayContrasts); % it finds the closest
%                      match to the requesed contrast given the current
%                      clut
% mglBltTexture(stimulus.tex{phase}(displayContrastNum+1),position); % it plots the pregenrated texture with the requested phase

global stimulus;

% now find closest matching contrast we can display with this gamma table
displayContrastNum = min(round(stimulus.nDisplayContrasts*desiredContrast/stimulus.currentMaxContrast),stimulus.nDisplayContrasts);
% disp(sprintf('Desired contrast: %0.4f Actual contrast: %0.4f Difference: %0.4f',desiredContrast,stimulus.currentMaxContrast*(displayContrastNum/stimulus.nDisplayContrasts),desiredContrast-stimulus.currentMaxContrast*(displayContrastNum/stimulus.nDisplayContrasts)));
if round(stimulus.nDisplayContrasts*desiredContrast/stimulus.currentMaxContrast)>stimulus.nDisplayContrasts
 disp(sprintf('(drawGabor) Desired contrast out of range %0.2f > %0.2f',desiredContrast,stimulus.currentMaxContrast));
 keyboard
end

mglBltTexture(stimulus.tex{phase}(displayContrastNum+1),position);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to create a gamma table %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function setGammaTable(maxContrast)

global stimulus;

% set the reserved colors
gammaTable(1:size(stimulus.reservedColors,1),1:size(stimulus.reservedColors,2))=stimulus.reservedColors;

% create the gamma table
cmax = 0.5+maxContrast/2;cmin = 0.5-maxContrast/2;
luminanceVals = cmin:((cmax-cmin)/(stimulus.nGratingColors-1)):cmax;

% no get the linearized range
redLinearized = interp1(0:1/255:1,stimulus.linearizedGammaTable.redTable,luminanceVals,'linear');
greenLinearized = interp1(0:1/255:1,stimulus.linearizedGammaTable.greenTable,luminanceVals,'linear');
blueLinearized = interp1(0:1/255:1,stimulus.linearizedGammaTable.blueTable,luminanceVals,'linear');

% check the table
% plot(stimulus.linearizedGammaTable.redTable,'k');
% hold on
% plot(256*(0.25:0.5/250:0.75),redLinearized,'r');

% set the gamma table
gammaTable((stimulus.minGratingColors+1):256,:)=[redLinearized;greenLinearized;blueLinearized]';

% set the gamma table
mglSetGammaTable(gammaTable);

% remember what the current maximum contrast is that we can display
stimulus.currentMaxContrast = maxContrast;


%% show results:
function results = crfaa_analquest(filename)
global stimulus;

indexes = size(stimulus.quest.q);

if length(indexes)==2
 for hh = 1:indexes(1) % pedestal
  for jj = 1:indexes(2) % attentional condition
   if ~isempty(stimulus.quest.q{hh,jj})
    qthreshold(hh,jj) = QuestMean(stimulus.quest.q{hh,jj});
    qthSD(hh,jj) = QuestSd(stimulus.quest.q{hh,jj});
    bEstimate(hh,jj) = QuestBetaAnalysis(stimulus.quest.q{hh,jj});
   end
  end
 end
 % making errorbars:
 se = qthSD;
 eBarsLow = (squeeze(qthreshold(:,1))-se(:,1));
 eBarsHi = (squeeze(qthreshold(:,1))+se(:,1));
 eBars{1} = [eBarsLow; flipud(eBarsHi)];

 eBarsLow = squeeze(qthreshold(:,2))-se(:,2);
 eBarsHi = squeeze(qthreshold(:,2))+se(:,2);
 eBars{2} = [eBarsLow; flipud(eBarsHi)];
 s = 1;
else
 for hh = 1:indexes(1) % pedestal
  for jj = 1:indexes(2) % attentional condition
   for ii = 1:indexes(3) % number of quest
    if ~isempty(stimulus.quest.q{hh,jj,ii})
     qthreshold(hh,jj,ii) = QuestMean(stimulus.quest.q{hh,jj,ii});
     qthSD(hh,jj,ii) = QuestSd(stimulus.quest.q{hh,jj,ii});
     bEstimate(hh,jj,ii) = QuestBetaAnalysis(stimulus.quest.q{hh,jj,ii});
    end
   end
  end
 end
 % making errorbars:
 s = size(qthreshold);
 %     se = sum(qthSD.^2,3)./s(3);
 %     eBarsLow = [(squeeze(qthreshold(:,1,:))-repmat(se(:,1),s(3),1))];
 %     eBarsHi = [(squeeze(qthreshold(:,1,:))+repmat(se(:,1),s(3),1))];
 %     eBars{1} = [eBarsLow; flipud(eBarsHi)];
 %
 %     eBarsLow = [(squeeze(qthreshold(:,2,:))-repmat(se(:,2),s(3),1))];
 %     eBarsHi = [(squeeze(qthreshold(:,2,:))+repmat(se(:,2),s(3),1))];
 %     eBars{2} = [eBarsLow; flipud(eBarsHi)];
 s = s(3);
end

if stimulus.pedestals(1) == 0;
 stimulus.pedestals(1) = .01;
end

results.meanThreshold = qthreshold;
results.thresholdSD = qthSD;
results.betaAnalysis = bEstimate;
results.pedestals = stimulus.pedestals;
resutlts.quests = stimulus.quest;

%%  plot:

fv = figure(1);clf;
plotColor = {'r.', 'k.'};
plotColorMedian = {'ro-', 'ks--'};
lightRed = [.99 .4 .1];
gray = [.75 .75 .75];

loglog(stimulus.pedestals(1),10.^median(squeeze(qthreshold(1,1,1)),2),'k.');
hold on
% C = [stimulus.pedestals(:)' stimulus.pedestals(:)'];
% for ii = 1:s
%     pt2(ii) = patch([C,fliplr(C)], ...
%         10.^squeeze(eBars{2}(:,ii)) ,gray);
% end
% for ii = 1:s
%     pt1(ii) = patch([C,fliplr(C)], ...
%         10.^squeeze(eBars{1}(:,ii)) ,lightRed);
% end

loglog(stimulus.pedestals(:),10.^median(squeeze(qthreshold(:,1,:)),2),plotColorMedian{1}, ...
 stimulus.pedestals(:),10.^median(squeeze(qthreshold(:,2,:)),2),plotColorMedian{2},...
 stimulus.pedestals(:),10.^squeeze(qthreshold(:,1,:)),plotColor{1}, ...
 stimulus.pedestals(:),10.^squeeze(qthreshold(:,2,:)),plotColor{2});

axis([.001, .1 .0005 1]);
axis('square');

%%%%%%%%%%%%%%%%%%%%
% figure formatting:
% colors:
myRED = [.7 .1 .1];
myBLACK = [0 0 0 ];
myGREEN = [.1 .5 .1];
myGRAY = [.7 .7 .7];
myLightRED = [.97 .79 .69];
myLightGREEN = [.69 .97 .79];

lineThick = 3;
thinLine = 1;
markerSize = 16;
Fsize = 30;
XYColor = [.2 .2 .2];

% set(fv,'NumberTitle','off','MenuBar','none');

xlabel('Pedestal contrast','FontName','Arial','FontSize',Fsize);
ylabel('Delta contrast','FontName','Arial','FontSize',Fsize);
title([filename, ' - TvC (Number of quests = ', num2str(s),')'],'FontName','Arial','FontSize',Fsize)

% set(pt1(:),'FaceAlpha', .1,'EdgeColor', lightRed, 'EdgeAlpha', .1);
% set(pt2(:),'FaceAlpha', .1,'EdgeColor', gray, 'EdgeAlpha', .1);

scrsz = get(0,'ScreenSize');
set(fv,'Position',[1 scrsz(4) scrsz(3)/2 scrsz(4)]);
ah = get(fv,'Child'); % getting axis
lh = get(ah,'Child'); % getting points and lines
set(ah,'FontName','Arial','FontSize',Fsize,'PlotBoxAspectRatio',[1 1 1], ...
 'LineWidth',thinLine,'TickLength',[0.025 .01], ...
 'XColor',XYColor,'XTick',[.0001 0.001 0.01 0.1 1],'XTickLabel',[0.0001 0.001 0.01 0.1 1], ...
 'YColor',XYColor,'Box','off', ...
 'YTick',[0.0001 0.001 0.01 0.1 1],'YTickLabel',[0.0001 0.001 0.01 0.1 1]);

% attended
foh = findobj(lh,'Marker','o');
set(foh,'MarkerSize',markerSize,'MarkerFaceColor',myRED,'MarkerEdgeColor',myRED)

%% neutral
foh = findobj(lh,'Marker','s');
set(foh,'MarkerSize',markerSize,'MarkerFaceColor','k','MarkerEdgeColor','k')

