function adp_runSelectionModelProclus(runType)
%
% Batch function that allows to run seelciton model fits on seevral
% condtions. Visual areas and observers.
%
% adp_runSelectionModelProclus(runType)
%
% INPUTS:
%  runType - is a number between 1 and 32 identifying different
%            combinations of adaptation, visual area and observers.
%
% Franco (c) Stanford University 2012
%
%  dispCDF  - Set to 0. Whether to dispaly (1) or not (0) the contrast discrimination functions
%  dispCRF  - Set to 0. Whether to display (1) of not (0) the contrast responst functions
%  dispCond - Set to 0. Whoch condtion to display: 
%             (a) dispCRF(s,m,CRF);
%             (b) dispCDF(s,m,CDF);
%  fitCDF   - Set this to 1. Fits the CRF first then uses the fit params to fit the CDF using the pooling model. 
%             This is the the Neuron fit.
%             Whether to fit or no the contrast discrimination function - that is do the model fit
%             fitting the CRF first, fixing those parameters and then only allowing sigma/k to change
%             to fit CDF.
%  fitCRFCDF- This is adjusts the CRF to best fit the CDF also at once so the parames of the CRF cna change depending on the CDF. 
%             Whether to fit both contrast response function and contrast discrimination function simultaneously. 
%             That is, run the model so that the CRF and the
%             CDF parameters all are allowed to change to fit both the CRF and CDF data at the same time.
%  verbose  - Show progress, 1
%  dispFig  - Display figure. Set to 0.
%  dispIter - Set to 0. display each  iteraton 1/0
%  recompute- recompute model - after running once, the fit parameters will be saved and running
%             with the same model parameters will load the fit parameters and display them w/out refitting.
%             If you want to force a refit, then set here.
%  savePathDir- Directory where to save fit results. 
%  saveName - file name to save the model parameters as. This is the stem of the filename - a unique identifier
%             for the model settings should be concatenated on to this to differentiate different settings.
%  crossValidate - cross valudate results
%  saveIter      - save each iteration 1/0
%  dispSavedIter - display the saved iterations 0/1
%  reslog        - compute residuals in log units - this set to the base of the log (e.g. 2)
%  fixedParams   - Run the model with the input parameters - does not fit.
%  initParams    - Initial parameters for the model.
%  reuseFig      - Always use the same figure to display results in instead of creating a new figure. 
%  noFig         - Do not display figure. This will overwrites all the
%                  previous disp fig parameters. If set to 0 will not display any figure.
%

% Detrmine the decies for the conditions to use
[adaptation, area, observer] = getConditions(runType);

% Input files and directories
dataPath = '/home/frk/Dropbox/adaptation/selectionModel';
dataFile = fullfile(dataPath,'08-Nov-2012_seletionModelData_Adp0_28_100.mat');

% Output files and directories
saveDir  = sprintf('savePathDir=''/home/frk/Dropbox/adaptation/selectionModel''');
saveFile = sprintf('saveName=smA%i_V%i_OBS%i',adaptation,area,observer);

% Load the data and run the process
load(dataFile);

% Run the fit
selectionModel(CRF(observer,area),CDF(observer,area),...
  'fitCDF=1',   ...
  'dispCDF',0,  ...
  'dispCRF',0,  ...
  'dispCond',0, ...
  'dispFig',0,  ...
  'dispIter',0, ...
   saveFile,saveDir)

end


%-------------------------------------%
function [adaptation, area, observer] = getConditions(runType)
switch runType
  %% ADP0 * 4 obsr * 4 visual areas
  % OBS = 1  
  case 1
    adaptation = 0;area = 1;observer = 1;
  case 2
    adaptation = 0;area = 2;observer = 1;
  case 3
    adaptation = 0;area = 3;observer = 1;
  case 4
    adaptation = 0;area = 4;observer = 1;
    % OBS =2
  case 5
    adaptation = 0;area = 1;observer = 2;
  case 6
    adaptation = 0;area = 2;observer = 2;
  case 7
    adaptation = 0;area = 3;observer = 2;
  case 8
    adaptation = 0;area = 4;observer = 2;
    % OBS =2
  case 9
    adaptation = 0;area = 1;observer = 3;
  case 10
    adaptation = 0;area = 2;observer = 3;
  case 11
    adaptation = 0;area = 3;observer = 3;
  case 12
    adaptation = 0;area = 4;observer = 3;
  case 13
    adaptation = 0;area = 1;observer = 4;
  case 14
    adaptation = 0;area = 2;observer = 4;
  case 15
    adaptation = 0;area = 3;observer = 4;
  case 16
    adaptation = 0;area = 4;observer = 4;
    
    
   
  %% ADP28 * 4 obsr * 4 visual areas
  % OBS = 1
  case 1+16
    adaptation = 28;area = 1;observer = 1;
  case 2+16
    adaptation = 28;area = 2;observer = 1;
  case 3+16
    adaptation = 28;area = 3;observer = 1;
  case 4+16
    adaptation = 28;area = 4;observer = 1;
    % OBS =2
  case 5+16
    adaptation = 28;area = 1;observer = 2;
  case 6+16
    adaptation = 28;area = 2;observer = 2;
  case 7+16
    adaptation = 28;area = 3;observer = 2;
  case 8+16
    adaptation = 28;area = 4;observer = 2;
    % OBS =3
  case 9+16
    adaptation = 28;area = 1;observer = 3;
  case 10+16
    adaptation = 28;area = 2;observer = 3;
  case 11+16
    adaptation = 28;area = 3;observer = 3;
  case 12+16
    adaptation = 28;area = 4;observer = 3;
     % OBS =4
 case 13+16
    adaptation = 28;area = 1;observer = 4;
  case 14+16
    adaptation = 28;area = 2;observer = 4;
  case 15+16
    adaptation = 28;area = 3;observer = 4;
  case 16+16
    adaptation = 28;area = 4;observer = 4;
    
    
  %% ADP28 * 3 obsr * 4 visual areas
  % OBS = 1
  case 1+16*2
    adaptation = 100;area = 1;observer = 1;
  case 2+16*2
    adaptation = 100;area = 2;observer = 1;
  case 3+16*2
    adaptation = 100;area = 3;observer = 1;
  case 4+16*2
    adaptation = 100;area = 4;observer = 1;
    % OBS =2
  case 5+16*2
    adaptation = 100;area = 1;observer = 2;
  case 6+16*2
    adaptation = 100;area = 2;observer = 2;
  case 7+16*2
    adaptation = 100;area = 3;observer = 2;
  case 8+16*2
    adaptation = 100;area = 4;observer = 2;
    % OBS =3
  case 9+16*2
    adaptation = 100;area = 1;observer = 3;
  case 10+16*2
    adaptation = 100;area = 2;observer = 3;
  case 11+16*2
    adaptation = 100;area = 3;observer = 3;
  case 12+16*2
    adaptation = 100;area = 4;observer = 3;
    % OBS = 1  
  
  case 1+12+16*2
    adaptation = 128;area = 1;observer = 1;
  case 2+12+16*2
    adaptation = 128;area = 2;observer = 1;
  case 3+12+16*2
    adaptation = 128;area = 3;observer = 1;
  case 4+12+16*2
    adaptation = 128;area = 4;observer = 1;
    % OBS =2
  case 5+12+16*2
    adaptation = 128;area = 1;observer = 2;
  case 6+12+16*2
    adaptation = 128;area = 2;observer = 2;
  case 7+12+16*2
    adaptation = 128;area = 3;observer = 2;
  case 8+12+16*2
    adaptation = 128;area = 4;observer = 2;
    % OBS =2
  case 9+12+16*2
    adaptation = 128;area = 1;observer = 3;
  case 10+12+16*2
    adaptation = 128;area = 2;observer = 3;
  case 11+12+16*2
    adaptation = 128;area = 3;observer = 3;
  case 12+12+16*2
    adaptation = 128;area = 4;observer = 3;
  case 13+12+16*2
    adaptation =128;area = 1;observer = 4;
  case 14+12+16*2
    adaptation = 128;area = 2;observer = 4;
  case 15+12+16*2
    adaptation = 128;area = 3;observer = 4;
  case 16+12+16*2
    adaptation = 128;area = 4;observer = 4;

  otherwise
    keyboard
end
  end