function runAllanalyses_r
% function runAllanalyses_r
%
% this function runs all the analyses to do do subtractiona nd fitTimecourse.m
% it is now set up to run the stream at riken. 
%
%
% franco pestilli 2009.06.29

SESSIONS = { ...
 '/Volumes/data/riken/crfaa/fmridata/crfaa/FP_A0_testStream/fp20071019/', ...
 '/Volumes/data/riken/crfaa/fmridata/crfaa/FP_A28_testStream/fp20080402/', ...
 '/Volumes/data/riken/crfaa/fmridata/crfaa/FP_A100_testStream/fp20080415/', ...
 '/Volumes/data/riken/crfaa/fmridata/crfaa/JG_A0_testStream/jg20070919/', ...
 '/Volumes/data/riken/crfaa/fmridata/crfaa/JG_A28_testStream/jg20080402/', ...
 '/Volumes/data/riken/crfaa/fmridata/crfaa/JG_A100_testStream/jg20080414/', ...
 '/Volumes/data/riken/crfaa/fmridata/crfaa/FM_A28_testStream/fm20080406/', ...
 '/Volumes/data/riken/crfaa/fmridata/crfaa/FM_A0_testStream/fm20080209/', ...
 };

% these are the rois used to estimate the STD of responses by combining
% across visual areas (v1-v4):
% roiType = {'combinedR206' 'combinedR207'};
%'combinedR2025' 'combinedR205' 'combinedR2075' 'combinedR2085' 'combinedR2088'};

% these are the original rois used to estimate the responses in each visual
% area:
roiType = {'v1' 'v2' 'v3' 'v4'};

for s = 1:length(SESSIONS) 
 disp('%  START ')
 disp('%  START ')
 disp('##########################################################################################################################')
 disp('##########################################################################################################################')
 disp(sprintf('###################### %s #############################',SESSIONS{s}))
 disp('##########################################################################################################################')
 disp('##########################################################################################################################')
 cd(SESSIONS{s})
 for roi = 1:length(roiType)
  % this is the original call to do the STD analysis used for the plot in
  % the paper.
    filename = sprintf('subTSnoiseSTD_%s',roiType{roi});
  
  % this is the call to do the correct/incorrect analysis
%   filename = sprintf('correctIncorrectAnalysis_STD_%s',roiType{roi});

  % this is the call to do the 1st/2nd interval analysis
%   filename = sprintf('firstSecondIntervalAnalysis_STD_%s',roiType{roi});
  
  close all; drawnow, gs = getGS;
  [d d1 ts events] = meanERA_ROI_remove0_script3(gs(1),gs(2),1,filename,roiType{roi});
 
  drawnow; 
  clear d d1 ts events gs; clear global MLR;
 end
 disp('##########################################################################################################################')
 disp('##########################################################################################################################')
 disp(sprintf('###################### %s #############################',SESSIONS{s}))
 disp('##########################################################################################################################')
 disp('##########################################################################################################################')
 disp('%  END ')
 disp('%  END ')

end