function p = crfFitStats(varargin)
% function p = crfFitStats(varargin)
%
% this function computes some statistics on the Noise fit as returned by
% the the test comparison (i.e., by the function
% makenoisestatisticsfigs2_r2_compareNoise.m)
%
% the data file used is:
% /Volumes/data/riken/crfaa/fmridata/crfaa/03-Aug-2009CompareNoiseFit_nkn2.mat
% /Volumes/data/riken/crfaa/fmridata/crfaa/30-Jul-2009CompareNoiseFit_sg.mat
%
% results are not saved. But a screen print out is given.
%
% NB the program canbe called with 'fitType' set to 'nk2' to do the test on
%    the naka-rushton fit, and with 'sg' to do the test on the skewed
%    gaussian fit.
%
% franco pestilli 2009.07.31
%
%
% pls not indexes are computed this way for the 1st/2nd interval and
% correct/incorrect analyses during the GLM analysis:
%
% Correct/incorrect:
% conditions types in this order:
%
% (1) focal cue, target - correct
% (2) focal cue, target - incorrect
% (3) focal cue, non-target
% (4) distributed cue, target - correct
% (5) distributed cue, target - incorrect
% (6) ditributed cue, non-target
%
%
% 1st and 2nd interval:
% conditions types in this order:
%
% (1) focal cue, target - 1st interval
% (2) focal cue, target - 2nd interval
% (3) focal cue, non-target
% (4) distributed cue, target - 1st interval
% (5) distributed cue, target - 2nd interval
% (6) ditributed cue, non-target
%
% !BUT!
% the same indexes are changed in the following way during the fitting procedure
% so to have consistency in color coding between the 4 conditions analyses (standard)
% and the 6 condition analyses (ci, 12)
%
% Correct/incorrect:
% conditions types in this order:
%
% (1) focal cue, target - correct
% (2) focal cue, target - incorrect
% (3) focal cue, non-target
% (4) distributed cue, target - correct
% (5) distributed cue, target - incorrect
% (6) ditributed cue, non-target
%
%
% 1st and 2nd interval:
% conditions types in this order:
%
% (1) focal cue, target - 1st interval
% (2) focal cue, target - 2nd interval
% (3) focal cue, non-target
% (4) distributed cue, target - 1st interval
% (5) distributed cue, target - 2nd interval
% (6) ditributed cue, non-target
%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
observer  =[];
filename  =[];
adaptation=[];
fitType   =[];     % load the naka-rushton fit or the skewed gaussian fit
savePlots =[];     % 1 saves the plot in the current dierectory
dispFit   =[];     % if 1 the plot is displayed
numBoots  =[];     % number of bootstrap samples
saveFig   =[];     % if 1 saves a figure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set default variables:
getArgs(varargin,{...
 'observer',{'fm' 'fp' 'jg'}, ...  
 'dispFit',0, ...
 'saveFig',1, ...
 'fitType', '12', ... % 's' or 'ci' or '12'
 'filename',[date,'allNoise.mat']});

defaultDataFolder = '/Volumes/data/riken/crfaa/fmridata/crfaa/nestedHypothesisTest/allOBS_results/';

% load the proper file.
[test offsets filename] = loadData(defaultDataFolder,fitType);
[test_s offsets_s filename] = loadData(defaultDataFolder,'s');

visareas = {'v1' 'v2' 'v3' 'v4'};

for t = 1:size(test,3)  % test type
 for o = 1%:size(test,1)
  for v = 1:size(test,2)  
           disp('%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%')
   disp(sprintf('%% Visual Area <%s> --  Test Type: <%s>  %%',visareas{v},test(o,v,t).test))
   disp(sprintf('%% Visual Area <%s> --  Observer#  <%i>      %%',visareas{v},o))
   disp(sprintf('%% R2 full    model                         %1.3f.',test(o,v,t).r2f))
   disp(sprintf('%% R2 reduced model                         %1.3f.',test(o,v,t).r2r))
   disp(sprintf('%% Probability of rejecting the full model: %1.3f.',test(o,v,t).p))
   disp(sprintf('%% Value of the F statistic               : %1.3f.',test(o,v,t).F))
   if test(o,v,t).p >= .05
           disp('%% !! The  winner  is  the  reduced  model :) !!')
   else
           disp('%% :( The  winner  is  the  full  model :(')    
   end
           disp('%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%')
  end
           disp('%%                                           %%')

 end
end


c = 0;
for v = 1:size(offsets,1)    % across visal areas
 for o = 1:size(offsets,2) % across observers
  c = c + 1;
  offset(c,1:size(offsets,3)) = offsets(v,o,:);
  firstCondition(c,1) = offsets(v,o,3);
  secondCondition(c,1) = offsets(v,o,6);
 end
end

% computed paired t-test for 1st/2nd intervals
a_ttest(log(firstCondition),log(secondCondition),1)


keyboard

% plot a bar graph of the test values:


%%%%%%%%%%%%
% loadData %
%%%%%%%%%%%%
function [test offsets filename] = loadData(defaultDataFolder,fitType)

% load the file
[filename vars2load] = makeFileName(defaultDataFolder,fitType);
data = load(sprintf('%s',filename),vars2load{1},vars2load{2});

% collect all test estimates in on array
for o = 1:size(data.test{1},1) % for each observer
 for v = 1:size(data.test{1},2) % for each visual area
  for t = 1:length(data.test)
   test(v,o,t).p = data.test{t}{o,v}.p;
   test(v,o,t).F = data.test{t}{o,v}.F;
   test(v,o,t).r2f = data.test{t}{o,v}.R2full;
   test(v,o,t).r2r = data.test{t}{o,v}.R2redux;
   test(v,o,t).test = data.test{t}{o,v}.redux;
   if t == 3 % collect only the offset test
    last = length(data.results{t}{o,v}.crf.fit.params(5:end));
    offsets(v,o,1:last) = data.results{t}{o,v}.crf.fit.params(5:end);
   end
  end
 end
end


%%%%%%%%%%%%%%%%
% makeFileName %
%%%%%%%%%%%%%%%%
function [filename vars2load] = makeFileName(defaultDataFolder,fitType)

switch fitType
 case {'standard' 's'}
  filename = [defaultDataFolder,'03-Aug-2009nestsedTests_standard.mat'];
 case {'12' 'firstsecond' 'fs'}
  filename = [defaultDataFolder,'05-Aug-2009nestsedTests_12.mat'];
 case {'correctincorrect' 'ci'}
  filename = [defaultDataFolder,'05-Aug-2009nestsedTests_CorrIncorr.mat']; 
 otherwise
  keyboard
end
vars2load = {'test' 'results'};



