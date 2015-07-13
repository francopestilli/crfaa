function crfeffect(varargin)
% function crfeffect
%
% this function:
% (1) loads the results of the crf fit,
% (2) strips out the correct test results
% (3) normalizes them (see code)
% (4) makes a plot of U,D,A by offset
%
%
% franco pestilli 2009/12/22

% (0) set up the arguments here
normalize = [];
getArgs(varargin,{'normalize=su','whichCue',1});

% (1) load the file, normalize it
[b crf] = loadData(normalize);

% (2) plot b
plot_b(b)

% (3) plot crf
plot_crf(crf,whichCue)




%%%%%%%%%%
% plot_b %
%%%%%%%%%%
function plot_b(b)
plotIndex = [4 2 1];

%          V1    V2    V3    V4
colors = {'ro-' 'bo-' 'mo-' 'go-'};
smartfig('baseline','reuse');
for o = 1%:size(b,1)
 for v = 1:size(b,2)
  plot([0 .26 .52],b{o,v}(plotIndex),colors{v})
  hold on
 end
end
plot([0 .26 .52],[0 .26 .52],'k--')
axis([0 .52 0 .52])
legend({'V1' 'V2' 'V3' 'V4'},2)
ylabel('b (fMRI response offset)')
xlabel('Cue conditions')
set(gca,'XTick', [0 .26 .52],'XTickLabel', {'focal cue, non-target'; 'distributed cue, target'; 'focal cue, target'})


%%%%%%%%%%%%
% plot_crf %
%%%%%%%%%%%%
function plot_crf(crf,whichCue)

%          V1    V2    V3    V4
colors = {'ro-' 'bo-' 'mo-' 'go-'};
cueconditions = {'Focal, target', 'Distributed, target' 'Distributed, non-target' 'Focal, non-target'};
smartfig('crfs','reuse');

% plot V1 vs V2-4
for o = 1%:size(crf,1)
 for v = 2:4
  plot(crf{o,1}(whichCue,:),crf{o,v}(whichCue,:),colors{v})
  hold on
 end
end
plot([0 1.3],[0 1.3],'k--')
axis([0 1.3 0 1.3]),
axis('square')
legend({'V2' 'V3' 'V4'},2)
ylabel('fMRI response V2,3,4')
xlabel('fMRI response V1')
title(sprintf('cue: %s',cueconditions{whichCue}))


%%%%%%%%%%%%
% loadData %
%%%%%%%%%%%%
function [b crf] = loadData(normalize)

% (0) load the data
d = load('~/data/riken/crfaa/fmridata/crfaa/data_used_files/nestedHypothesisTest/allOBS_results/03-Aug-2009nestsedTests_standard.mat');

% (1) strip out the beseline values
offsetFitIndex = 3; % index to the fit of the offset only, this is the statistically winning model, reported in the attention paper
d = d.results{offsetFitIndex};

for o = 1:size(d,1)
 for v = 1:size(d,2)
  b{o,v} = d{o,v}.crf.fit.eachparam.offset;
 end
end

% (2) normalize the baseline values
switch normalize
 case {'su' 'subtract_unattended'}
  for o = 1:size(b,1)
   for v = 1:size(b,2)
    b{o,v} = b{o,v} - b{o,v}(4);
   end
  end
 case {'no' 'nonormalization'}
  % do nothing
 otherwise
  keyboard
end

% (3) strip out the crfs
for o = 1:size(d,1)
 for v = 1:size(d,2)
  crf{o,v} = d{o,v}.d.use_crf;
 end
end

