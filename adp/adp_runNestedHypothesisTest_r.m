function results = adp_runNestedHypothesisTest_r(varargin)
% function adp_runNestedHypothesisTest_r(varargin)
%
% version _r is the one used at riken.
% this function the nested hypothesis test.
% 
% it tests whether the full model is better than a model with some reduced
% parameters.
%
% e.g., it tests whether allowing all parameters of the contrast response
% function to vary freely across attention (code 2) or adaptation (code 3)
% is better than allowing only c50 to vary across attention (code 2)
% given the number fo parameters of the full model and that of the reduced
% model.
%
% varargin (Default):
% {'adaptation',{0 28 100},'visualAreas',{'v1','v2','v3','v4'}, ...
%  'doTestOnly',0,'testType','0101'};
% -
% data files are loaded from: 
%  - the default data folder for the files used by this function at riken is:
%   /Users/frakkopesto/data/riken/crfaa/fmridata/crfaa/data_used_files/
%
% Figures are saved in:
%  - /Users/frakkopesto/data/riken/crfaa/fmridata/crfaa/adp_nestedHypothesisTest/figs_nestedFits_all;
%
% Results are saved in:
%  - /Users/frakkopesto/data/riken/crfaa/fmridata/crfaa/adp_nestedHypothesisTest/allOBS_results;
% 
% It calls the fucntions:
% - adp_makeFitSimIndividualPlotsNested.m
% - ...
%
% This code:
% 1. it is written for the adaptation paper, and runs across
% attention only.
% 2. it only uses the naka-rushton function with two exponents
%
% Test codes:
% - '0' there is one of the parameters set to 0 to fit all the functions
% - '1' means that the parameter is free to vary across conditions (e.g., acrss adaptation)
%
% Implemented test types:
% 1) 11001 rmax, c50 and offset free, n and q common
% 2) 01001 c50 and offset free, rmax, n and q common
% 3) 10001 rmax and offset free, c50, n and q common
% 4) 10000 rmax free, c50, offset, n and q common
% 5) 01000 c50 free, rmax, offset, n and q common
% 6) 00001 offset free, rmax, c50, n and q common
%
%
% franco pestilli 2010.10.07


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
visualAreas = []; % selects which visual area to use (to load the file)
whichCRF    = []; % selects which contrast response function in that visual area
doTestOnly  = []; % only the statistical test without plots
analysisType= []; % choose which analysis files to load, standard is the main results, 1st/2nd and corr/incorrect are the other options
testType    = []; % sets the test type
crfFittype  = []; % sets the type of crf-model to use (only nk2 implemented at the moment)
doIndividualPlots = []; % plots individual graphs formatted for pubblication, one for each visual area and test type
displayAllFits    = []; % it displays a figure with each individual crf best fit for the current model, it saves the figure
saveData          = []; % if 1, saves the results
observer          = []; % selects the observers to load the data (average also)
adaptation        = []; % which adaptation contiion for the current observer, default for the attention paper is '0'
saveFig           = []; % saves the figure just plotted (check functionality it has been inherited)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

getArgs(varargin,{...
 'visualAreas',{'v1','v2','v3','v4'}, ... 
 'doTestOnly',0, ...
 'analysisType', 'standard',...'standard', ... % 12 or ci
 'crfFittype','nk2', ...
 'testType',{'00001'}, ... '01000' '10000' '11001', '01001', '10001', '10000', '01000', '00001'}, ...
 'displayAllFits',0, ...
 'saveData',1, ...
 'observer',{'avrg' 'fp' 'jg' 'fm'}, ... 'fm' 'fp' 'jg'
 'doIndividualPlots',1, ...
 'adaptation',0, ...
 'saveFig',1});
%  'whichCRF',{'dnt' 'dt' 'a' 'u'}, ...
global defaultDataFolder;

% run trhough the requested test types, check that they are all strings:
for t = 1:length(testType)
 if ~ischar(testType{t})
  disp(sprintf('(%s) testType must be a string!',mfilename));
  return
 end
end

% set figure filder (to save) and data folder (to save):
figDir = ['figs_nestedAttTests_',date];
defaultDataFolder = '/Users/frakkopesto/data/riken/crfaa/fmridata/crfaa/data_used_files/';

% now run through the requested observers, by default the first observer is
% the average data across observers
for o = 1:length(observer)      % do the whole thing for each observer
 disp(sprintf('[adp_runNestedHypothesisTest_r] Running fit on Observer: %s',observer{o}))
 disp(sprintf('[adp_runNestedHypothesisTest_r] Running fit on Observer: %s',observer{o}))
 disp(sprintf('[adp_runNestedHypothesisTest_r] Running fit on Observer: %s',observer{o}))
 
 % check that we are not trying to run the correct/incorrect or 1st/2nd
 % interval analysis on the averge responses
 if strcmpi(observer{o},'avrg')
  switch analysisType
   case {'CorrectIncorrect' 'corrincorr' 'CorrIncorr' 'ci' '1st2ndInterval' '1st2nd' '12'}
    disp('[adp_runNestedHypothesisTest_r] ERROR! only the 4 conditions analysis was computed for the average.')
    keyboard
  end
 end

 % perform the requested test: first runs the full model
 % then runs each reduced model (NB this is different that previous versions of test files)
 for t = 1:length(testType)
  % do the tests by visual area. each visual area independently tested.
  for v = 1:length(visualAreas) % v1-v4
   
   % load data:
   [d meanC] = loadData(observer{o},visualAreas{v},defaultDataFolder,adaptation,analysisType);
   if isempty(d) || isempty(meanC)
    disp(sprintf('[adp_runNestedHypothesisTest_r] WARNING! Subject fm does not have the 100%% adaptation, exiting.\n[adp_runNestedHypothesisTest_r] To incur in no problems please run fit on observer fm as last in a series of observers.'));
    return
   end

   % set up test parameters:
   [fullTest nParamsFull nParamsRedux nDataPoints] = setUpTest(testType{t},observer{o});
   
   % NB this is to debug the fits.
   % set parameters for saving intermadiate fit figure, this is useful for
   % looking at each fit individually by attention and adaptation
   % set .save to 0 if not wanted.
   % it saves a figure in a folder called 'testnake2_bestfit_figs'
   % in the current folder:
   saveFitFig.save = 1;
   saveFitFig.observer   = observer{o};
   saveFitFig.testType   = testType{t};
   saveFitFig.visualarea = visualAreas{v};
   saveFitFig.displayAllFits= 0; % this displays the idividual curves being fitted
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % (1) compute the reduced model %
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   results{t}{o,v} =  makeFitSimIndividualPlotsNested(d,meanC, ...
    'dispFit', 1, ...
    'crfFittype',crfFittype, ...
    'testType',testType{t}, ...
    'visualArea',visualAreas{v}, ...
    'displayAllFits',displayAllFits, ...
    'saveFig',1, ...
    'saveFitFig',saveFitFig, ...
    'observer',observer{o}, ...
    'analysisType',analysisType);
   
   results{t}{o,v}.mfilename = mfilename;
   results{t}{o,v}.observer = observer{o};
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % (2) now compute the full model to do test %
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % - do not display the fit this time
   % - also do it only once at the begining and keep the same value
   % for the rest of the tests
   % makeFitSimIndividualPlotsNested
   saveFitFig.testType   = fullTest;
   fullResults{t}{o,v} =  makeFitSimIndividualPlotsNested(d,meanC, ...
    'dispFit', 0,...
    'crfFittype',crfFittype, ...
    'testType',fullTest, ...
    'visualArea',visualAreas{v}, ...
    'saveFitFig',saveFitFig, ...
    'displayAllFits',displayAllFits, ...
    'analysisType',analysisType);
   
   % store full and reduced test R2 for F-test:
   R2full = fullResults{t}{o,v}.crf.fit.model_r2;
   R2redux = results{t}{o,v}.crf.fit.model_r2;
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%
   % (3) now compute F-test: %
   %%%%%%%%%%%%%%%%%%%%%%%%%%%
   test{t}{o,v} = fTestR2(R2full,R2redux,nParamsFull,nParamsRedux,nDataPoints);
   test{t}{o,v}.redux = testType{t};
   test{t}{o,v}.full  = fullTest;
   
   % save data
   thisDate = [datestr(now,10),datestr(now,7),datestr(now,5)];
   fileName = sprintf('%s_%s_nestedFTest',observer{o},deblank(thisDate));
   figureName = [fileName,'_',analysisType];
   
   % save each figure:
   if ~doIndividualPlots
    if saveFig
     defaultDataFolder = '/Users/frakkopesto/data/riken/crfaa/fmridata/crfaa/adp_nestedHypothesisTest/';
     figDir = 'figs_nestedFits_all';
     savefig(results{end}{o,v}.plotInfo.figureHandle,'figName',[fileName,testType{t}],'defaultDataFolder',defaultDataFolder,'figDir',figDir);
    end
   end
  end
  
  %%%%%%%%%%%%%%%%%
  % save the data %
  %%%%%%%%%%%%%%%%%
  if saveData
   saveDataFolder = '/Users/frakkopesto/data/riken/crfaa/fmridata/crfaa/adp_nestedHypothesisTest/';
   fiDir = 'allOBS_results';
   fileName = [date,'nestsedTests','_',analysisType];
   saveresults(results, fullResults, test,saveDataFolder,fiDir,fileName)
  else
   disp(sprintf('%s: Results not saved',mfilename));
  end
  
 end
%  close all, drawnow
 disp(sprintf('[adp_runNestedHypothesisTest_r] Running fit on Observer: %s',observer{o}))
 disp(sprintf('[adp_runNestedHypothesisTest_r] Running fit on Observer: %s',observer{o}))
 disp(sprintf('[adp_runNestedHypothesisTest_r] Running fit on Observer: %s',observer{o}))
end


%%%%%%%%%%%%%%%%%%%%%%%%%%
% (4) plot test results: %
%%%%%%%%%%%%%%%%%%%%%%%%%%
saveFig = 1;
doOneGraph = 1;
whichAdapt = [1:3];
% plotFtest(test,figureName,saveFig);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (5) plot the scatter plot of the fit parameters %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot a scatter plot of the parameters obtained in the full model:
% simFitIndScatterPlot(fullResults,saveFig,doOneGraph,observer)


% %%%%%%%%%%%%% END main call %%%%%%%%%%%%%%% %



%%%%%%%%%%%%
% loadData %
%%%%%%%%%%%%
function [d meanC] = loadData(observer,visualAreas,defaultDataFolder,adaptation,analysisType)

% find adaptation index:
if adaptation == 0
 adp = 1;
elseif adaptation == 28
 adp = 2;
elseif adaptation == 100
 adp = 3;
else
 keyboard;
end
d = [];
meanC = [];

% load each file:
if strcmpi(observer,'avrg')
 defaultDataFolder = '/Users/frakkopesto/data/riken/crfaa/fmridata/crfaa/data_used_files/';
 
 % (0) load data file
 [filename vars2load] = makeFileName(observer,[],visualAreas,defaultDataFolder, analysisType);
 [data] = load(sprintf('%s',filename),vars2load{1}, vars2load{2});
 d = data.d{adp};
 % NB for some reason v1 data of the average has a different threshold at the first
 % pedestal contrast, i am using the v2 mean contrast instead
 if strcmpi(visualAreas,'v1')
  [filename vars2load] = makeFileName(observer,[],'v2',defaultDataFolder,analysisType);
  [data] = load(sprintf('%s',filename),vars2load{1}, vars2load{2});
  meanC = data.meanC{adp};
 else
  meanC = data.meanC{adp};
 end
 clear data;
 
else % individual observers
 if ~(strcmpi('fm',observer) && adaptation == 100)
  defaultDataFolder = '/Users/frakkopesto/data/riken/crfaa/fmridata/crfaa/data_used_files/';
  
  % (0) load data file
  [filename vars2load] = makeFileName(observer,adaptation,visualAreas,defaultDataFolder, analysisType);
  data = load(sprintf('%s',filename),vars2load{1}, vars2load{2});
  
  % now make the average contrast
  meanC = computedMeanC(data.events,.5);
  meanC.quantileAcontrast = 100*meanC.quantileAcontrast';
  meanC.quantileDcontrast = 100*meanC.quantileDcontrast';
  meanC.pedC = 100*meanC.pedC;
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % now sort the conditions depending on the type of analysis loaded %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  switch analysisType
   case {'standard' '4conditions' '4c' '4' 's'}
    % mean
    d.A_amplitude   = data.d1.amplitude(1:8);
    d.U_amplitude   = data.d1.amplitude(9:16);
    d.Dt_amplitude  = data.d1.amplitude(17:24);
    d.Dnt_amplitude = data.d1.amplitude(25:32);
    
    % ste
    d.A_ste         = data.d1.amplitudeSTE(1:8);
    d.U_ste         = data.d1.amplitudeSTE(9:16);
    d.Dt_ste        = data.d1.amplitudeSTE(17:24);
    d.Dnt_ste       = data.d1.amplitudeSTE(25:32);
   
   case {'CorrectIncorrect' 'corrincorr' 'CorrIncorr' 'ci' '1st2ndInterval' '1st2nd' '12'}        
    % mean
    d.A1_amplitude   = data.d1.amplitude(1:8); % t attended correct
    d.A2_amplitude   = data.d1.amplitude(9:16);% t attende incorrect

    d.U_amplitude   = data.d1.amplitude(17:24);
    
    d.Dt1_amplitude  = data.d1.amplitude(25:32);
    d.Dt2_amplitude  = data.d1.amplitude(33:40);
    
    d.Dnt_amplitude = data.d1.amplitude(41:48);
    
    % ste
    d.A1_ste         = data.d1.amplitudeSTE(9:16);
    d.A2_ste         = data.d1.amplitudeSTE(1:8);
    
    d.U_ste         = data.d1.amplitudeSTE(17:24);
    d.Dt1_ste        = data.d1.amplitudeSTE(25:32);
    d.Dt2_ste        = data.d1.amplitudeSTE(33:40);
    
    d.Dnt_ste       = data.d1.amplitudeSTE(41:48);
     
   otherwise
    keyboard
  end
 end
end

if isempty(d) || isempty(meanC)
 return
end

%%%%%%%%%%%%%
% setUpTest %
%%%%%%%%%%%%%
function [fullTest nParamsFull nParamsRedux nDataPoints] = setUpTest(testType,observer)
% set up F-test parameters

% depending on the type of test being performed:
switch testType
 case { '10000' '01000' '00001'}
  fullTest = '11001';
  nParamsFull  = (4*3+2);  % 4 attention X 3 parameters free + 2 common
  nParamsRedux = (4+4);  % 4 attention X 2 parameters free + 3 common
  
 case {'01001'  '10001'}
  fullTest = '11001';
  nParamsFull  = (4*3+2);  % 4 attention X 3 parameters free + 2 common
  nParamsRedux = (2*4+3);  % 4 attention X 2 parameters free + 3 common
  
 case {'11001'}
  fullTest = '11111';
  nParamsFull  = (4*5);  % 4 attention X 3 parameters free + 2 common
  nParamsRedux = (3*4+2);  % 4 attention X 2 parameters free + 3 common
  
 otherwise
  keyboard
end

nDataPoints  = (8*4);    % 8 contrasts X 4 attention


%%%%%%%%%%%%%%%%%%%%%%%%
% simFitIndScatterPlot %
%%%%%%%%%%%%%%%%%%%%%%%%
function simFitIndScatterPlot(fullResults,saveFig,doOneGraph,observer)
% function noiseScatterPlot(noise)
% *symbols* are observers
% *colors* are visual areas
% *hue* is adaptation

if ieNotDefined('doOneGraph')
 doOneGraph = 1;
end

% attention figure
if doOneGraph
 figurename = 'ATTsimfitScatterPlot_c50';
 h = smartfig(figurename,'reuse');
 set(h,'Name',figurename);
end
% set up basic plot info:
plotInfo = setUpPlotInfo;
plotInfo.thisSymbol = {'s' 'o' '^' 'd'}; % one symbol per observer
plotInfo.title = 'all values';
plotInfo.MarkerSize = 12;
plotInfo.YLim = [.001 1];
plotInfo.Xlim = [.001 1];
visareas = {'v1' 'v2' 'v3' 'v4'};

% Attention plot:
for v = 1:size(fullResults,2)   % visual areas
 % it will make a graph for each visual area
 if ~doOneGraph
  figurename = sprintf('ATTsimfitScatterPlot_c50_%s',visareas{v});
  h = smartfig(figurename,'reuse');
  set(h,'Name',figurename);
 end
 for o = 1:size(fullResults,1) % observers
  % fm does not have the 100% adaptation
  if strcmp('fm',observer{o})
   adaptations = 1:2;
  else
   adaptations = 1:3;
  end
  
  for j = 1:adaptations % adaptation conditions
   
   % c50:
   c50    = fullResults{o,v}.crf.fit.eachparam.c50(j,:); % first dimension is adaptation the second attention
   
   loglog([1,1;5,5],[1,1;5,5],'k-','Color',plotInfo.XYColor)
   hold on;
   % plotting atteded vs. unattended:
   loglog(c50(4),c50(1),plotInfo.thisSymbol{o}, ...
    'LineWidth',plotInfo.LineWidth/2, ...
    'Color',plotInfo.thisColor{v}{j}, ...
    'MarkerFaceColor',plotInfo.MarkerFaceColor{v}{j}, ...
    'MarkerEdgeColor',plotInfo.MarkerEdgeColor, ...
    'MarkerSize', plotInfo.MarkerSize);
   
   %    title(sprintf('%s - c50',plotInfo.title))
   %    if v == 4
   %     xlabel('Distributed','FontSize',plotInfo.Fsize);
   %     ylabel(sprintf('Attended'),'FontSize',plotInfo.Fsize);
   %    end
   
   set(gca,...
    'FontName','Helvetica','FontSize',plotInfo.Fsize, ...
    'PlotBoxAspectRatio',plotInfo.PlotBoxAspectRatio, ...
    'XLim', [1 4], ...
    'YLim', [1 4],...
    'LineWidth',plotInfo.LineWidth,'TickLength',plotInfo.TickLength, ...
    'XTick', [1 2 4], 'XTickLabel', [1 2 4] ,...
    'YTick', [1 2 4], 'YTickLabel', [1 2 4] ,...
    'XColor',plotInfo.XYColor,'YColor',plotInfo.XYColor,'Box','off');
   
  end
 end
 if ~doOneGraph
  myaxisTvC(0,1);
  if saveFig
   set(h,'PaperPosition',[.25 .25 8 10.5]);
   set(h,'PaperOrientation','Portrait');
   
   disp(sprintf('Saving figure %s',figurename));
   figDir = savefig(h,'figName',figurename,'figDir','fig_noise_scatter');
  end
 end
end

if doOneGraph
 myaxisTvC(0,1);
 if saveFig
  set(h,'PaperPosition',[.25 .25 8 10.5]);
  set(h,'PaperOrientation','Portrait');
  
  disp(sprintf('Saving figure %s',figurename));
  figDir = savefig(h,'figName',figurename,'figDir','fig_noise_scatter');
 end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%
% Attention OFFSET plot: %
%%%%%%%%%%%%%%%%%%%%%%%%%%
if doOneGraph
 figurename = 'ATTsimfitScatterPlot_offset';
 h = smartfig(figurename,'reuse');
 set(h,'Name',figurename);
end

plotInfo.XTicks = [0 .25 .5 .75 1];
plotInfo.YTicks = [0 .25 .5 .75 1];

for v = 1:size(fullResults,2)   % visual areas
 % it will make a graph for each visual area
 if ~doOneGraph
  figurename = sprintf('ATTsimfitScatterPlot_offset_%s',visareas{v});
  h = smartfig(figurename,'reuse');
  set(h,'Name',figurename);
 end
 
 for o = 1:size(fullResults,1) % observers
  % fm does not have the 100% adaptation
  if strcmp('fm',observer{o})
   adaptations = 1:2;
  else
   adaptations = 1:3;
  end
  
  for j = 1:adaptations % adaptation conditions
   
   offset = fullResults{o,v}.crf.fit.eachparam.offset(j,:); % first dimension is adaptation the second attention
   
   % response off-set
   plot([0,0;1,1],[0,0;1,1],'k-','Color',plotInfo.XYColor)
   hold on;
   plot(offset(4),offset(1), ...
    plotInfo.thisSymbol{o}, ...
    'LineWidth',plotInfo.LineWidth/2, ...
    'Color',plotInfo.thisColor{v}{j},...
    'MarkerFaceColor',plotInfo.MarkerFaceColor{v}{j}, ...
    'MarkerEdgeColor',plotInfo.MarkerEdgeColor, ...
    'MarkerSize', plotInfo.MarkerSize);
   
   %    title(sprintf('%s - offset',plotInfo.title))
   %    if v == 4
   %     xlabel('Distributed','FontSize',plotInfo.Fsize);
   %     ylabel('Attended','FontSize',plotInfo.Fsize);
   %    end
   set(gca,...
    'FontName','Helvetica','FontSize',plotInfo.Fsize, ...
    'PlotBoxAspectRatio',plotInfo.PlotBoxAspectRatio, ...
    'XLim', [0 1],'YLim', [0 1],...
    'LineWidth',plotInfo.LineWidth,'TickLength',plotInfo.TickLength, ...
    'XTick', [0 .5 1], 'XTickLabel', [0 .5 1] ,...
    'YTick', [0 .5 1], 'YTickLabel', [0 .5 1] ,...
    'XColor',plotInfo.XYColor,'YColor',plotInfo.XYColor,'Box','off');
  end
 end
 if ~doOneGraph
  myaxis;
  if saveFig
   set(h,'PaperPosition',[.25 .25 8 10.5]);
   set(h,'PaperOrientation','Portrait');
   
   disp(sprintf('Saving figure %s',figurename));
   figDir = savefig(h,'figName',figurename,'figDir','fig_noise_scatter');
  end
 end
end

if doOneGraph
 myaxis;
 if saveFig
  set(h,'PaperPosition',[.25 .25 8 10.5]);
  set(h,'PaperOrientation','Portrait');
  
  disp(sprintf('Saving figure %s',figurename));
  figDir = savefig(h,'figName',figurename,'figDir','fig_noise_scatter');
 end
end

%%%%%%%%%%%%%%%%%%%%
% adaptation plot: %
%%%%%%%%%%%%%%%%%%%%
% k parameter
if doOneGraph
 figurename = 'ADPsimfitScatterPlot_c50';
 h = smartfig(figurename,'reuse');
 set(h,'Name',figurename);
end

for v = 1:size(fullResults,2) % visual area
 % it will make a graph for each visual area
 if ~doOneGraph
  figurename = sprintf('ADPsimfitScatterPlot_c50_%s',visareas{v});
  h = smartfig(figurename,'reuse');
  set(h,'Name',figurename);
 end
 
 for j = 1:4 % attention
  for o = 1:size(fullResults,1) % observers
   
   c50 = fullResults{o,v}.crf.fit.eachparam.c50(:,j); % first dimension is adaptation the second attention
   
   loglog([.5,.5;5,5],[.5,.5;5,5],'k-','Color',plotInfo.XYColor)
   hold on;
   % adapt-0 and adapt-28
   loglog(c50(1),c50(2),plotInfo.thisSymbol{o}, ...
    'LineWidth',plotInfo.LineWidth/2, ...
    'Color',plotInfo.thisColor{v}{j}, ...
    'MarkerFaceColor',plotInfo.MarkerFaceColor{v}{j}, ...
    'MarkerEdgeColor',plotInfo.MarkerEdgeColor, ...
    'MarkerSize', plotInfo.MarkerSize);
   
   % fm does not have 100%
   if ~strcmp('fm',observer{o})
    loglog(c50(1),c50(3),plotInfo.thisSymbol{o}, ...
     'LineWidth',plotInfo.LineWidth/2, ...
     'Color',plotInfo.thisColor{v}{j}, ...
     'MarkerFaceColor',plotInfo.MarkerFaceColor{v}{j}, ...
     'MarkerEdgeColor',plotInfo.MarkerEdgeColor, ...
     'MarkerSize', plotInfo.MarkerSize);
   end
   %    title(sprintf('%s - k',plotInfo.title))
   %    if v == 4
   %     xlabel('Adapt 0% contrast','FontSize',plotInfo.Fsize);
   %     ylabel('Adapt 28 or 100% contrast','FontSize',plotInfo.Fsize);
   %    end
   set(gca,...
    'FontName','Helvetica','FontSize',plotInfo.Fsize, ...
    'PlotBoxAspectRatio',plotInfo.PlotBoxAspectRatio, ...
    'XLim', [1 4], ...
    'YLim', [1 4],...
    'LineWidth',plotInfo.LineWidth,'TickLength',plotInfo.TickLength, ...
    'XTick', [1 2 4], 'XTickLabel', [1 2 4] ,...
    'YTick', [1 2 4], 'YTickLabel', [1 2 4] ,...
    'XColor',plotInfo.XYColor,'YColor',plotInfo.XYColor,'Box','off');
  end
 end
 if ~doOneGraph
  myaxisTvC(0,1);
  if saveFig
   set(h,'PaperPosition',[.25 .25 8 10.5]);
   set(h,'PaperOrientation','Portrait');
   
   disp(sprintf('Saving figure %s',figurename));
   figDir = savefig(h,'figName',figurename,'figDir','fig_noise_scatter');
  end
 end
end

if doOneGraph
 myaxisTvC(0,1);
 if saveFig
  set(h,'PaperPosition',[.25 .25 8 10.5]);
  set(h,'PaperOrientation','Portrait');
  
  disp(sprintf('Saving figure %s',figurename));
  figDir = savefig(h,'figName',figurename,'figDir','fig_noise_scatter');
 end
end


% offset plot:
if doOneGraph
 figurename = 'ADPsimfitScatterPlot_offset';
 h = smartfig(figurename,'reuse');
 set(h,'Name',figurename);
end

for v = 1:size(fullResults,2) % visual area
 % it will make a graph for each visual area
 if ~doOneGraph
  figurename = sprintf('ADPsimfitScatterPlot_offset_%s',visareas{v});
  h = smartfig(figurename,'reuse');
  set(h,'Name',figurename);
 end
 
 for j = 1:4 % attention
  for o = 1:size(fullResults,1) % observers
   
   offset = fullResults{o,v}.crf.fit.eachparam.offset(:,j); % first dimension is adaptation the second attention
   
   % response off-set
   % adapt-0 and adapt-28
   plot([0,0;.5,.5],[0,0;.5,.5],'k-','Color',plotInfo.XYColor)
   hold on;
   plot(offset(1),offset(2), ...
    plotInfo.thisSymbol{o}, ...
    'LineWidth',plotInfo.LineWidth/2,...
    'Color',plotInfo.thisColor{v}{j}, ...
    'MarkerFaceColor',plotInfo.MarkerFaceColor{v}{j}, ...
    'MarkerEdgeColor',plotInfo.MarkerEdgeColor, ...
    'MarkerSize', plotInfo.MarkerSize);
   
   % fm does not have 100%
   if ~strcmp('fm',observer{o})
    plot(offset(1),offset(3), ...
     plotInfo.thisSymbol{o}, ...
     'LineWidth',plotInfo.LineWidth/2, ...
     'Color',plotInfo.thisColor{v}{j}, ...
     'MarkerFaceColor',plotInfo.MarkerFaceColor{v}{j}, ...
     'MarkerEdgeColor',plotInfo.MarkerEdgeColor, ...
     'MarkerSize', plotInfo.MarkerSize);
   end
   %    title(sprintf('%s - offset',plotInfo.title))
   %    if v == 4
   %     xlabel('Adapt 0% contrast','FontSize',plotInfo.Fsize);
   %     ylabel('Adapt 28 or 100% contrast','FontSize',plotInfo.Fsize);
   %    end
   set(gca,...
    'FontName','Helvetica','FontSize',plotInfo.Fsize, ...
    'PlotBoxAspectRatio',plotInfo.PlotBoxAspectRatio, ...
    'XLim', [0 .5], ...
    'YLim', [0 .5],...
    'LineWidth',plotInfo.LineWidth,'TickLength',plotInfo.TickLength, ...
    'XTick', [0 .25 .5], 'XTickLabel', [0 .25 .5] ,...
    'YTick', [0 .25 .5], 'YTickLabel', [0 .25 .5] ,...
    'XColor',plotInfo.XYColor,'YColor',plotInfo.XYColor,'Box','off');
  end
 end
 
 if ~doOneGraph
  myaxis;
  if saveFig
   set(h,'PaperPosition',[.25 .25 8 10.5]);
   set(h,'PaperOrientation','Portrait');
   
   disp(sprintf('Saving figure %s',figurename));
   figDir = savefig(h,'figName',figurename,'figDir','fig_noise_scatter');
  end
 end
end

if doOneGraph
 myaxis;
 if saveFig
  set(h,'PaperPosition',[.25 .25 8 10.5]);
  set(h,'PaperOrientation','Portrait');
  
  disp(sprintf('Saving figure %s',figurename));
  figDir = savefig(h,'figName',figurename,'figDir','fig_noise_scatter');
 end
end


%%%%%%%%%%%%%%%%%
% setUpPlotInfo %
%%%%%%%%%%%%%%%%%
function plotInfo = setUpPlotInfo
% figure formatting info:
plotInfo.MarkerEdgeColor = [.4 .4 .4];
plotInfo.thisSymbol = {'o-' 'o-' 'o-' 'o-'};
plotInfo.plotPeds = [0 1.75 3.5 7 14 28 56 84];

% generate some nice colors:
numC =60;
for i = 1:numC
 c{i} = getSmoothColor(i,numC,'hsv');
 %  plot(i,1,'ko','MarkerFaceColor',c{i},'MarkerSize',26);hold on
 %  text(i,1,sprintf('%i',i),'HorizontalAlignment','Center','Color',[0 0 0]);
end

% choose the ones i like:
reds = {c{2} c{5} c{7} c{11}};
greens = {c{23} c{17} c{15} c{13}};
blues = {c{40} c{37} c{34} c{30}};
purples = {c{46} c{49} c{54} c{50}};

for i = 1:4 % visual area
 plotInfo.thisColor{i} = {reds{i} blues{i} purples{i} greens{i}};
 plotInfo.MarkerFaceColor{i} = plotInfo.thisColor{i};
end

% plot basic set up:
plotInfo.XYColor = [0 0 0];
plotInfo.Fsize = 6;
plotInfo.MarkerSize = 10;
plotInfo.LineWidth = 2;
plotInfo.TickLength = [0.025 .01];
plotInfo.PlotBoxAspectRatio = [1 1 1];
plotInfo.YLim = [0 1.3];
plotInfo.Xlim = [.001,1];
plotInfo.plotPeds = [plotInfo.Xlim(1) 0.0175 0.035 0.07 0.14 0.28 0.56 0.84];


%%%%%%%%%%%%%%%%
% makeFileName %
%%%%%%%%%%%%%%%%
function [filename vars2load] = makeFileName(observer,adaptation,visualArea,defaultDataFolder,analysisType)
% NB whereas the main analysis files are saved into the special folder
% <data_used_files> so that changes to that folder do not screw up the
% files...
% the files for the other analyses are instad in each folder/session. so
% for construtcting the filename we use getsession.

if strcmpi(observer,'avrg')
 filename = sprintf('%sADAPTATION-EFFECT-%s-19-Jan-2009.mat',defaultDataFolder,visualArea);
 vars2load = {'d' 'meanC'};
else
 % choose the type of analysis files to load
 switch analysisType
  case {'standard' '4conditions' '4c' '4' 's'} % this is the original analysis for the main results
   switch visualArea
    case 'v1'
     filename = sprintf('%s%s_A%s_29-Dec-2008_%s_0_7roi.mat',defaultDataFolder,upper(observer),num2str(adaptation),visualArea);
    case {'v2' 'v3' 'v4'}
     filename = sprintf('%s%s_A%s_29-Dec-2008_%s_0_5roi.mat',defaultDataFolder,upper(observer),num2str(adaptation),visualArea);
    otherwise
     keyboard
   end
   
  case {'CorrectIncorrect' 'corrincorr' 'CorrIncorr' 'ci'} % this is the correct incorrect analysis
   session = getsession(observer,adaptation);
   filename = sprintf('%s/28-Jul-2009_%s_correctIncorrectAnalysis_STD_%s.mat',session,visualArea,visualArea);

  case {'1st2ndInterval' '1st2nd' '12'} % This is the first and second interval analysis
   session = getsession(observer,adaptation);
    filename = sprintf('%s/28-Jul-2009_%s_firstSecondIntervalAnalysis_STD_%s.mat',session,visualArea,visualArea);
  otherwise
   keyboard
 end
 vars2load = {'d1' 'events'};
end

%%%%%%%%%%%%%%%%
%  getsession  %
%%%%%%%%%%%%%%%%
function session = getsession(observer,adaptation)

% check which computer are we using:
[notok c] = system('hostname');
if ~(~notok && strcmp('riken.jp',deblank(c(end-8:end))))
 % check if there is a drive connected which one is it ...
 if isdir('/Volumes/homosacer_data')
  hd = '/Volumes/homosacer_data';
 elseif isdir('/Volumes/homosacer')
  hd = '/Volumes/homosacer';
 elseif isdir('/Volumes/homosacer_backup')
  hd = '/Volumes/homosacer_backup';
 else
  disp('No data hard drive found.')
  keyboard
 end
elseif (~notok && strcmp('riken.jp',deblank(c(end-8:end))))
 hd = '/Users/frakkopesto/data/riken/crfaa';
end

switch lower(observer)
 case {'fm'}
  switch adaptation
   case {0}
    mrsession = 'fm20080209';
   case 28
    mrsession = 'fm20080406';
   otherwise
    disp(sprintf('This condition is NOT set (%s, adp-%i).',observer,adaptation))
    mrsession = [];
  end
  
 case {'fp'}
  switch adaptation
   case {0}
    mrsession = 'fp20071019';
   case {28}
    mrsession = 'fp20080402';
   case {100}
    mrsession = 'fp20080415';
   otherwise
    disp(sprintf('This condition is NOT set (%s, adp-%i).',observer,adaptation))
    keyboard
  end
  
 case {'jg'}
  switch adaptation
   case {0}
    mrsession = 'jg20070919';
   case {28}
    mrsession = 'jg20080402';
   case {100}
    mrsession = 'jg20080414';
   otherwise
    disp(sprintf('This condition is NOT set (%s, adp-%i).',observer,adaptation))
    keyboard
  end
  
 otherwise
  disp(sprintf('Observer (%s) NOT found.',observer))
  keyboard
end

if ~isempty(mrsession)
 session = deblank(sprintf('%s/fmridata/crfaa/%s_A%i_testStream/%s',hd,upper(observer),adaptation,mrsession));
else
 session = '';
end
disp(sprintf('Loading: %s',session))


%%%%%%%%%%%%%
% plotFtest %
%%%%%%%%%%%%%
function plotFtest(test,figurename,saveFig)
% test as to be a struct as returned by makeFtest
global defaultDataFolder;

for t = 1:length(test)
 figName = deblank(char([figurename{1} figurename{t+1}]));
 h = smartfig(figName,'reuse');
 
 % combine al F and p values
 for v = 1:length(test{t})
  toplot(v,1) = test{t}{v}.F;
  toplot(v,2) = test{t}{v}.p;
 end
 
 % num of params:
 df1 = num2str(test{t}{1}.df1);
 df2 = num2str(test{t}{1}.df2);
 mybar(toplot, ...
  'groupLabels',{'v1','v2', 'v3' 'v4'}, ...
  'withinGroupLabels',{sprintf('F(%s,%s)',df1,df2),'p'}, ...
  'yAxisMin=0', ...
  sprintf('yAxisMax=%s',num2str(round(max(toplot(:))))),...
  'xLabelText', ' ', ...
  'yLabelText',sprintf('F value and\n probability'), ...
  'dispValues',1, ...
  'hline',0.05);
 
 if saveFig
  figDir = savefig(h,'figName',figName,'defaultDataFolder',defaultDataFolder,'figDir','testFigs');
 end
end

% now make a summary plot of all type of tests:
figName = deblank(char([figurename{1} figurename{2},'_',figurename{3}]));
h = smartfig(figName,'reuse');

% combine al F and p values
for t = 1:length(test)
 for v = 1:length(test{t})
  toplot(v,t) = test{t}{v}.p;
 end
end

% num of params:
mybar(toplot, ...
 'groupLabels',{'v1','v2', 'v3' 'v4'}, ...
 'withinGroupLabels', ...
 {sprintf('c50 - adaptation\noffset - attention'),...
 sprintf('c50 - attention\noffset - adaptation')}, ...
 'yAxisMin=0', ...
 'yAxisMax=1',...
 'xLabelText', ' ', ...
 'yLabelText',sprintf('Probability of the full-model\nbeing better by chance.'), ...
 'dispValues',1, ...
 'hline',0.05);

if saveFig
 figDir = savefig(h,'figName',figName,'defaultDataFolder',defaultDataFolder,'figDir','testFigs');
end


%%%%%%%%%%%%%%%
% saveresults %
%%%%%%%%%%%%%%%
function fiDir = saveresults(results, fullResults, test, defaultDataFolder,fiDir,filename)
% save the figure if requested
fiDir = sprintf('%s%s',defaultDataFolder,fiDir);
currentDir = pwd;
if isdir(fiDir)
 cd(fiDir);
else
 mkdir(fiDir);
 cd(fiDir);
end

disp(sprintf('[%s] saving data %s/%s',mfilename,fiDir,filename));
eval(sprintf('save(''%s.mat'', ''%s'', ''%s'', ''%s'')',filename,'results','fullResults', 'test'));
cd(currentDir);


%%%%%%%%%%%%
%  fTestR2 %
%%%%%%%%%%%%
function test = fTestR2(R2full,R2redux,nParamsFull,nParamsRedux,nDataPoints)
% test best model given number of paramters and r-squares:
% this is testing the null hypothesis that the additional
% parameter is not needed.
%
%   The F-test tests the probability that the increase in r2 given by
%   the additional parameters of the full model is due to chance.
%   if F is big and p is small the full model wins over the reduced model
%
% franco pestilli 2009.01.23

% save intput to output:
test.R2full = R2full;
test.R2redux = R2redux;
test.nParamsFull = nParamsFull;
test.nParamsRedux = nParamsRedux;
test.nDataPoints = nDataPoints;

% do the test
test.df1 = nParamsFull-nParamsRedux;
test.df2 = nDataPoints - nParamsFull-1;
test.F = ((R2full-R2redux)/test.df1)./((1-R2full)/test.df2);
test.p = fpdf(test.F,test.df1,test.df2);
test.p(test.p>1) = .999;





