function c = computedMeanC(events,myQuantiles)
% function c = computedMeanC(events,myQuantiles)
% compute average contrast presented for attended and distributed-target
% conditions:
%
%
% franco pestilli 2008.12.20

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute average contrast presented for attended and distributed-target conditions: %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. get number of totla trials and the pedestal contrast:
trialNum = size(events.allEvents,2);
c.pedC = unique(events.allEvents(9,:));

% 2. atteneded and distributed indexes:
attendedTrials = events.allEvents(1,:) == ones(1,trialNum);
distributedTrials = ~attendedTrials;

% 3. sort attended and unattended target contrasts into vectors
for pnum = 1:length(c.pedC)
 c.deltaCattended{pnum} = [];
 c.deltaCdistributed{pnum} = [];
end

skippedCount = 0;
aCount = 0;
dtCount = 0;
% on each trial
for tnum = 1:trialNum
 % check the location
 for tLocation = unique(events.allEvents(3,:))
  % attended target: if this is an attended trial and the right target
  % location save the delta contrast
  if (events.allEvents(1,tnum) == 1 && events.allEvents(3,tnum) == tLocation)
   aCount = aCount +1;
   for pnum = 1:length(c.pedC)
    if events.allEvents(9,tnum) == c.pedC(pnum)
%      %disp(sprintf('[computedMeanC] trial %i attended target, pedNum %i',aCount, pnum));
     c.deltaCattended{pnum} = [c.deltaCattended{pnum} events.allEvents(8,tnum)];
    end
   end
  % distributed target: if this is a distributed trial and the
  % right target location , save the delta contrast
  elseif (events.allEvents(1,tnum) == 2  && events.allEvents(3,tnum) == tLocation)
   dtCount = dtCount +1;
   for pnum = 1:length(c.pedC)
    if events.allEvents(9,tnum) == c.pedC(pnum)
%      %disp(sprintf('[computedMeanC] trial %i distributed target, pedNum %i',dtCount,pnum));
     c.deltaCdistributed{pnum} = [c.deltaCdistributed{pnum} events.allEvents(8,tnum)];
    end
   end
  else
   skippedCount = skippedCount +1;
%    %disp(sprintf('[computedMeanC] trial %i skipped',skippedCount));
  end
 end
end

% 4. computed mean and sd for attended and unattended tcontrast:
if ieNotDefined('myQuantiles')
 myQuantiles = [.33 .5 .66];
end
c.quantiles = myQuantiles;
c.quantileAcontrast = zeros(length(c.pedC),length(c.quantiles));
c.quantileDcontrast = c.quantileAcontrast;
c.quantileAdeltaCs = c.quantileAcontrast;
c.quantileDdeltaCs = c.quantileAcontrast;
for pnum = 1:length(c.pedC)
 % Contrast presented on average at the target location:
 c.quantileAcontrast(pnum,:)   = c.pedC(pnum)+quantile(c.deltaCattended{pnum}/2,myQuantiles);
 c.quantileDcontrast(pnum,:)   = c.pedC(pnum)+quantile(c.deltaCdistributed{pnum}/2,myQuantiles);

 % DeltaC presented on average at the attended target location:
 c.quantileAdeltaCs(pnum,:)   = quantile(c.deltaCattended{pnum},myQuantiles);
 c.quantileDdeltaCs(pnum,:)   = quantile(c.deltaCdistributed{pnum},myQuantiles);

 % actual contrast presented on average at the distributed target
 % location:
 c.numAtrials(pnum,:) = length(c.deltaCattended{pnum});
 c.numDtrials(pnum,:) = length(c.deltaCdistributed{pnum});
end


