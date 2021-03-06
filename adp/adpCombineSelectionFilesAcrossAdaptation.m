function [CRF, CDF] = adpCombineSelectionFilesAcrossAdaptation(selectionFileName, saveFile)
%
% (1) Loads up selectionModel files (see adpSensoryNoiseFileTOSelectionFile.m
%     and selectionModel.m)
% (2) Creates a copy of the first adaptation file (adap-0) and adds fields
%     for adapt 28 and adap 100.
%
%   [CRF, CDF] = adpCombineSelectionFilesAcrossAdaptation(snFname,selectionFname)
%
% INPUTS:
%     selectionFileName - full path to a three files generated by the function:
%                         adpSensoryNoiseFileTOSelectionFile.m. One file
%                         per adaptation condition, 0, 28 and 100
%     saveFile             - 1      = save a file in the same folder of the
%                                     sensoryNoiseFileName.
%                            string = fullpath to a file name to use for saving the CRF
%                                     and CDF.
%                            0      = do not save a file.
%
% %% ** selectionFile structure spec ** %% 
%  CRF - a structure containign the data and the conditions for the contrast
%        brain response, the contrast response function (CRF). 
%
%        CRF is a structure that has two different contrast response functions (cued and uncued)
%        for each condition you are testing. The names of the conditions are arbitrary but
%        each condition must contain a cued and uncued CRF. So, for example if you have a 
%        focal cue condition then you would have:         
%        - CRF.focal.cued   (target) 
%        - CRF.focal.uncued (non-target)
%                
%        There should also be a field to tell whether the cued and uncued functions are the
%        same (i.e. for a distributed cue condition or different):
%        - CRF.focal.sameCuedUncued       = false;  
%        - CRF.distributed.sameCuedUncued = true; % in the distributed cue condition target 
%                                                 % and non-target are the same.
%        For each of the cued/uncued fields above there should be fields for data and fit
%        which contain the contrasts and responses. For example, something like:
%        - CRF.focal.cued.data.c   = [0 0.1 0.5 1]; % Contrast of the stimulus on the 
%                                                   % x-axis of the CRF)
%        - CRF.focal.uncued.data.r = [0.01 0.3 0.7 1.10]; % The actual CRF BOLD responses
%  
%  CDF - a structure containing the data and the condition for the
%        behavior, the contrast discrimination function (CDF)
%  
%        CDF is a similar struncture containing contrast discrimination data. So, continuing
%        from the exmaple above, if you have a focal cue condition you would have fields for
%        both data and fit that contain contrast discrimination data:
%        - CDF.focal.data.c     = [0.01 0.5 0.7];  % Pedestal contrast 
%        - CDF.focal.data.t     = [0.04 0.08 0.12];% Threshold, delta-c
%        - CDF.focal.data.tste  = [0.01 0.02 0.01];% ste of threshold
%        - CDF.focal.data.d     = {[0.0 0.01 0.02], ...
%                                  [0.01 0.5 0.7],  ...
%                                  [0.01 0.5 0.7]};
%        - CDF.focal.data.nCued = 1;
% 
%        Where c are the pedestal contrasts, t are the measured thresholds and d is the
%        distractor contrasts in a cell array. The cued field contains how many targets
%        were cued. The tste field contains the standard error of the CDF data (for display)
%
% EXAMPLE:
%   sensoryNoiseFileName = fullpath(adp_defaultDataFolder,'04-Nov-2012CompareNoiseFit_nkn2_A0.mat');
%   
%
% Franco Pestilli (c) Stanford University 2012

if ieNotDefined('saveFile'), saveFile = 1;end
if ieNotDefined('selectionFileName'), selectionFileName = {       ...
    '~/Dropbox/adaptation/selectionModel/08-Nov-2012CompareNoiseFit_nkn2_A0_seletionModelData.mat',  ...
    '~/Dropbox/adaptation/selectionModel/08-Nov-2012CompareNoiseFit_nkn2_A28_seletionModelData.mat', ...
    '~/Dropbox/adaptation/selectionModel/08-Nov-2012CompareNoiseFit_nkn2_A100_seletionModelData.mat'};
end


% Load the sensoryFile
a0   = load(selectionFileName{1});
a28  = load(selectionFileName{2});
a100 = load(selectionFileName{3});

%% CRF - contrast response function
for iArea = 1:size(a0.CRF,2)
  for iObs  = 1:size(a0.CRF,1)
    %-% Focal cue condition
    CRF(iObs,iArea).focal0   =   a0.CRF(iObs,iArea).focal;
    CRF(iObs,iArea).focal28  =  a28.CRF(iObs,iArea).focal;
    if iObs <=3
      CRF(iObs,iArea).focal100 = a100.CRF(iObs,iArea).focal;
    end
    
    %-% Distributed cue condition
    CRF(iObs,iArea).distributed0     =   a0.CRF(iObs,iArea).distributed;
    CRF(iObs,iArea).distributed28    =  a28.CRF(iObs,iArea).distributed;
    if iObs <=3
      CRF(iObs,iArea).distributed100 = a100.CRF(iObs,iArea).distributed;
    end
  end
end

%% CDF - contrast discrimination function
for iArea = 1:size(a0.CDF,2)
  for iObs  = 1:size(a0.CDF,1)
      %-% Focal cue condition
    CDF(iObs,iArea).focal0     =   a0.CDF(iObs,iArea).focal;
    CDF(iObs,iArea).focal28    =  a28.CDF(iObs,iArea).focal;
    if iObs <=3
      CDF(iObs,iArea).focal100 = a100.CDF(iObs,iArea).focal;
    end
    
    %-% Distributed cue condition
    CDF(iObs,iArea).distributed0     =   a0.CDF(iObs,iArea).distributed;
    CDF(iObs,iArea).distributed28    =  a28.CDF(iObs,iArea).distributed;
    if iObs <=3
      CDF(iObs,iArea).distributed100 = a100.CDF(iObs,iArea).distributed;
    end
  end
end

% Save file to disk
if ischar( saveFile ) % A full path to a file was passed-in
  fprintf('[%s] Saving file:\n %s\n',mfilename,saveFile);
  save(saveFile,'CRF','CDF');
  
elseif (saveFile == 1) % The file will be saved in the save folder of the sensoryNoiseFile
  [p,n,~] = fileparts(selectionFileName{1});
  
  fprintf('[%s] Saving file:\n %s\n',mfilename,fullfile(p,sprintf('%s_%s',n,'seletionModelDataAdp0_28_100.mat')))
  save(fullfile(p,sprintf('%s_%s',n,'seletionModelDataAdp0_28_100.mat')),'CRF','CDF');
end

return

