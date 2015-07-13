function runAllanalyses_r_R1(subtractionType)
% function runAllanalyses_r_R1(subtractionType)
%
% this function runs all the analyses to do do subtractiona nd fitTimecourse.m
% it is now set up to run the stream at riken. 
%
% version R1 is the one used to reply to the reviewers in Nn
% after the first submission
%
% subtractionType determines which type of 0-response to subtract.
% 1 = original 0-contrast focal non-target response
% 2 = reviewers' suggested 0-contrast response subtracted by trial type 
%     (e.g., non-target focal from target focal, and non-target distributed 
%      from target distributed)
%
% franco pestilli 2010.04.04

SESSIONS = { ...
  '/data2/crfaa/crfaa/FP_A0_testStream/fp20071019/'     ...
  '/data2/crfaa/crfaa/FM_A0_testStream/fm20080209/'     ...
  '/data2/crfaa/crfaa/JG_A0_testStream/jg20070919/'}; % ... 
%  '/data2/crfaa/crfaa/FP_A0_testStream/fp20071019/'      };
%  '/data2/crfaa/crfaa/JG_A100_testStream/jg20080414/', ...
%   '/data2/crfaa/crfaa/JG_A28_testStream/jg20080402/',  ...
%  '/data2/crfaa/crfaa/FM_A28_testStream/fm20080406/',  ...
%  '/data2/crfaa/crfaa/FP_A28_testStream/fp20080402/'  };
%  '/data2/crfaa/crfaa/FP_A100_testStream/fp20080415/'  };
% these are the rois used to estimate the STD of responses by combining
% across visual areas (v1-v4):
% roiType = {'combinedR206' 'combinedR207'};
%'combinedR2025' 'combinedR205' 'combinedR2075' 'combinedR2085' 'combinedR2088'};

% these are the original rois used 
% to estimate the responses in each visual
roiType = {  'v1' 'v2' 'v3' 'v4'};%...
%            'v1r030' 'v1r050' 'v1' 'v1r090' ...
%            'v2r030' 'v2' 'v2r070' ...
%            'v3r030' 'v3' 'v3r070' ...
%            'v4r030' 'v4' 'v4r070'};
 
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
  close all; drawnow; gs = getGS();
  
  switch subtractionType
   case {1}
    subtraction = 'ORIG';
   case {2}
    subtraction = 'BYTRIAL';
   case {3}
    subtraction = 'BOOT';
    recomputedDeconv = 0; % 1=re-do the deconvolution and fit time course stuff.
   otherwise
    keyboard
  end
  
  filename = sprintf('era_R1_subtr%s_%s',subtraction,roiType{roi});
  
  % chose which function to call to make the requested subtraction
  switch subtractionType
   case {1}
    [d d1 ts events] = meanERA_ROI_remove0_script3(gs(1),gs(2),1,filename,roiType{roi});
   case {2}
    [d d1 ts events] = meanERA_ROI_remove0_script3_R1(gs(1),gs(2),1,filename,roiType{roi});
   case {3}
    meanERA_ROI_remove0_script3_R1_boot(gs(1),gs(2),1,filename,roiType{roi},recomputedDeconv,SESSIONS{s});
   otherwise
    keyboard
  end
  
  drawnow; 
  clear d d1 ts events gs; 
  clear global MLR;
 end
 disp('##########################################################################################################################')
 disp('##########################################################################################################################')
 disp(sprintf('###################### %s #############################',SESSIONS{s}))
 disp('##########################################################################################################################')
 disp('##########################################################################################################################')
 disp('%  END ')
 disp('%  END ')

end