% crfaa.m
%
%        $Id: crfaa.m,v 1.1 2007/01/04 19:40:56 justin Exp $
%      usage: crfaa
%         by: justin gardner modified by franco pestilli
%       date: 09/07/06 - 01/11/07
%  copyright: (c) 2006 Justin Gardner (GPL see mgl/COPYING)
%    purpose: example program to show how to use the task structure
%
function myscreen = crfaa_localizer(observer)

% check arguments
if ~any(nargin == [1])
    help crfaa_localizer
    return
end

thisdir = pwd;
% make a data directory if necessary
if ~isdir(fullfile(thisdir,'data'))
    disp('Making data directory');
    mkdir('data');
end

% make an observer directory if necessary
datadirname = fullfile(thisdir,'data',[observer,'Localizer']);
if ~isdir(datadirname);
    disp(sprintf('Making observer directory %s',datadirname));
    mkdir(datadirname);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initalize the screen
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mglScreenCoordinates;
myscreen.datadir = datadirname;
myscreen.allowpause = 1;
myscreen.saveData = 1;
myscreen.displayName = 'projector';
% initalize the screen
myscreen.screenParams{1} = {'Mockscanner',[],2,1024,768,70,[32 25],60,1,0,1.8,'/Users/frakkopesto/fmritools/mgl/task/displays/0001_network-administrators-power-mac-g4_070228.mat'};
myscreen = initScreen(myscreen);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialize the task
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

task{1}.waitForBacktick = 1;
task{1}.seglen   =   6*[1 1 1 1];
[task{1} myscreen] = initTask(task{1},myscreen,@localizerStartSegmentCallback,@localizerDrawStimulusCallback);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialize the stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global stimulus;
myscreen = initStimulus('stimulus',myscreen);
stimulus = myInitStimulus(stimulus,myscreen);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% run the eye calibration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
myscreen = eyeCalibDisp(myscreen);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main display loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
phaseNum = 1;
while (phaseNum <= length(task)) && ~myscreen.userHitEsc
    % update the task
    [task myscreen phaseNum] = updateTask(task,myscreen,phaseNum);
    % flip screen
    myscreen = tickScreen(myscreen,task);
end
clear stimulus.tmp

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% end of the experiment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
myscreen = endTask(myscreen,task);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TASK 1: function that gets called at the start of each segment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = adaptationStartSegmentCallback(task, myscreen)
global stimulus;
setGammaTable(stimulus.adaptationContrast);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TASK 1: function that gets called to draw the stimulus each frame
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = localizerDrawStimulusCallback(task, myscreen)
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
    % NE SE
    phaseNum = mod(floor(2*mglGetSecs(task.thistrial.trialstart)/stimulus.flickerPeriod),2)+1;
    drawGabor(stimulus.adaptationContrast,stimulus.eccentricity*[stimulus.locations{1}],phaseNum);
    drawGabor(stimulus.adaptationContrast,stimulus.eccentricity*[stimulus.locations{2}],phaseNum);
    mglFixationCross(1,1,stimulus.black);

elseif (task.thistrial.thisseg == 2)
    % SW SE
    phaseNum = mod(floor(2*mglGetSecs(task.thistrial.trialstart)/stimulus.flickerPeriod),2)+1;
    drawGabor(stimulus.adaptationContrast,stimulus.eccentricity*[stimulus.locations{2}],phaseNum);
    drawGabor(stimulus.adaptationContrast,stimulus.eccentricity*[stimulus.locations{3}],phaseNum);
    mglFixationCross(1,1,stimulus.black);

elseif (task.thistrial.thisseg == 3)
    % SW SE
    phaseNum = mod(floor(2*mglGetSecs(task.thistrial.trialstart)/stimulus.flickerPeriod),2)+1;
    drawGabor(stimulus.adaptationContrast,stimulus.eccentricity*[stimulus.locations{3}],phaseNum);
    drawGabor(stimulus.adaptationContrast,stimulus.eccentricity*[stimulus.locations{4}],phaseNum);
    mglFixationCross(1,1,stimulus.black);

elseif (task.thistrial.thisseg == 4)
    % NE NW
    phaseNum = mod(floor(2*mglGetSecs(task.thistrial.trialstart)/stimulus.flickerPeriod),2)+1;
    drawGabor(stimulus.adaptationContrast,stimulus.eccentricity*[stimulus.locations{1}],phaseNum);
    drawGabor(stimulus.adaptationContrast,stimulus.eccentricity*[stimulus.locations{4}],phaseNum);
    mglFixationCross(1,1,stimulus.black);

end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TASK 1: function that gets called at the start of each segment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = localizerStartSegmentCallback(task, myscreen)
global stimulus;

setGammaTable(stimulus.adaptationContrast);
mglFixationCross(1,1,stimulus.black);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to init the stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stimulus = myInitStimulus(stimulus,myscreen)
global stimulus;

% let's get the linearized gamma table
stimulus.linearizedGammaTable = mglGetGammaTable;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% stimulus parameters:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
stimulus.init = 1;
stimulus.width = 8;          % in deg
stimulus.height = 8;         % in deg
stimulus.sf = 2;             % in cpd
stimulus.orientation = 90;   % in deg
stimulus.phase = [0 180];    % in deg
stimulus.eccentricity = 5.2;   % in deg


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% frames and locations:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
stimulus.locations = {[1 1], [1 -1], [-1 -1], [-1 1]};
stimulus.frameThick = .15;

% gabors
% gabors
% gabors
stimulus.gaussSdx = stimulus.width/7;
stimulus.gaussSdy = stimulus.height/7;
stimulus.adaptationContrast = 1;
stimulus.reservedColors = [0 0 0; 1 1 1; 0 .65 0; 1 0 0];
stimulus.flickerFrequency = 2; % in Hz

%% some computations:

% for stimulus flicker frequency
stimulus.flickerPeriod = 1/stimulus.flickerFrequency;

% for contrast
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
disppercent(-inf,'Calculating gabors');
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
mglClearScreen(stimulus.grayColor);


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to draw the gabor stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function drawGabor(desiredContrast,position,phase)
% drawGaborPedCdeltaC
%
%        $Id: drawGabor.m, v 1 2007/01/18 19:40:56 ?? ?? $
%      usage: drawGabor(desiredContrast,position,pedestalNum,phase)
%         by: justin gardner, modified by franco pestilli
%       date: 01/11/07 - 01/18/07
%  copyright: (c) 2007 Justin Gardner
%    purpose: draw a gabor stimulus on the screen with a specified contrast
%             within a given clut range (it finds the closest contrast to
%             the requested one within the available clut range)

global stimulus;

% now find closest matching contrast we can display with this gamma table
displayContrastNum = min(round(stimulus.nDisplayContrasts*desiredContrast/stimulus.currentMaxContrast),stimulus.nDisplayContrasts);;
% disp(sprintf('Desired contrast: %0.4f Actual contrast: %0.4f Difference: %0.4f',desiredContrast,stimulus.currentMaxContrast*(displayContrastNum/stimulus.nDisplayContrasts),desiredContrast-stimulus.currentMaxContrast*(displayContrastNum/stimulus.nDisplayContrasts)));
if round(stimulus.nDisplayContrasts*desiredContrast/stimulus.currentMaxContrast)>stimulus.nDisplayContrasts
    disp(sprintf('(drawGabor) Desired contrast out of range %0.2f > %0.2f',desiredContrast,stimulus.currentMaxContrast));
    keyboard
end

mglBltTexture(stimulus.tex{phase}(displayContrastNum+1),position);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to create a gamma table
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
%plot(stimulus.linearizedGammaTable.redTable,'k');
%hold on
%plot(256*(0.25:0.5/250:0.75),redLinearized,'r');

% set the gamma table
gammaTable((stimulus.minGratingColors+1):256,:)=[redLinearized;greenLinearized;blueLinearized]';

% set the gamma table
mglSetGammaTable(gammaTable);

% remember what the current maximum contrast is that we can display
stimulus.currentMaxContrast = maxContrast;
