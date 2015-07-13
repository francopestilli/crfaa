
numContrasts = 8;
titleString = 'Figure 2A: large att, small stim';
stimWidth = 3; 
AxWidth = 7.5;
Atheta = 0;
AthetaWidth = 360;
baselineMod = 5e-7; 
baselineMod = 5e-7; 
baselineUnmod = 0;
cRange = [1e-5 1];
Apeak = 4;

IthetaWidth = 360;
% Sampling of space and orientation
x = [-200:200];
theta = [-180:180]';

% Make stimuli
stimCenter1 = 100;
stimOrientation1 = 0;
stimCenter2 = -100;
stimOrientation2 = 0;
stim1 = makeGaussian(theta,stimOrientation1,1,1) * makeGaussian(x,stimCenter1,stimWidth,1); 
stim2 = makeGaussian(theta,stimOrientation2,1,1) * makeGaussian(x,stimCenter2,stimWidth,1);

% Pick contrasts
logCRange = log10(cRange);
logContrasts = linspace(logCRange(1),logCRange(2),numContrasts);
contrasts = 10.^logContrasts;

% Pick neuron to record
j = find(theta==stimOrientation1);
i = find(x==stimCenter1);

attCRF = zeros(size(contrasts));
unattCRF = zeros(size(contrasts));
for c = 1:numContrasts
  stim = contrasts(c) * stim1 + contrasts(c) * stim2;
  if c == numContrasts
    showActivityMaps = 1;
    showModelParameters = 1;
  else
    showActivityMaps = 0;
    showModelParameters = 0;
  end
  % Population response when attending stim 1
  R1 = attentionModel(x,theta,stim,'Ax',stimCenter1,'AxWidth',AxWidth,'AthetaWidth',AthetaWidth,'Atheta',Atheta,'showActivityMaps',showActivityMaps,'showModelParameters',showModelParameters,'baselineMod',baselineMod,'baselineUnmod',baselineUnmod,'IthetaWidth',IthetaWidth,'Apeak',Apeak);
  % Population response when attending stim 2
  R2 = attentionModel(x,theta,stim,'Ax',stimCenter2,'AxWidth',AxWidth,'AthetaWidth',AthetaWidth,'Atheta',Atheta,'baselineMod',baselineMod,'baselineUnmod',baselineUnmod,'IthetaWidth',IthetaWidth,'Apeak',Apeak);
  attCRF(c) = R1(j,i);
  unattCRF(c) = R2(j,i);
  R1c(c,:,:) = R1;
  R2c(c,:,:) = R2;
end

smartfig('attModelTest3'); clf;
subplot(1,2,1)
semilogx(contrasts,unattCRF,contrasts,attCRF);
ylim([0 30]);
xlim(cRange);
legend('Att Away','Att RF');
ylabel('Normalized response');
xlabel('Log contrast');
title(titleString);
subplot(1,2,2);
semilogx(contrasts,100*(attCRF-unattCRF)./unattCRF);
ylim([0 100]);
xlim(cRange);
ylabel('Attentional modulation (%)');
xlabel('Log contrast');
drawnow

spaceSample = -5:1:5;
orientSample = -90:15:90;
nOrient = length(orientSample);
nSpace = length(spaceSample);
f = smartfig('attModelTest','reuse');clf;
R1sum = [];R2sum = [];spaceIndex = [];orientIndex = [];
for iSpace = 1:nSpace
  for iOrient = 1:nOrient
    subplot(nOrient,nSpace,iSpace + (iOrient-1)*nSpace)
    orientIndex(iOrient) = find((orientSample(iOrient)+stimOrientation1) == theta);
    spaceIndex(iSpace) = find((spaceSample(iSpace)+stimCenter1) == x);
    semilogx(contrasts,R2c(:,orientIndex(iOrient),spaceIndex(iSpace)),contrasts,R1c(:,orientIndex(iOrient),spaceIndex(iSpace)));
    if (spaceSample(iSpace) == 0) && (orientSample(iOrient) == 0)
      set(gca,'Color',[0.8 0.8 0.8]);
    end
    axis off
    hold on
    drawnow
    if (iSpace == round(nSpace/2)) && (iOrient == nOrient)
      xlabel('Contrast');
    end
    if isempty(R1sum)
      R1sum = R1c(:,orientIndex(iOrient),spaceIndex(iSpace));
      R2sum = R2c(:,orientIndex(iOrient),spaceIndex(iSpace));
    else
      R1sum = R1sum + R1c(:,orientIndex(iOrient),spaceIndex(iSpace));
      R2sum = R2sum + R2c(:,orientIndex(iOrient),spaceIndex(iSpace));
    end
      
  end
end
setSubplotSameRange(f);

% get where we are going to sum over to compute fMRI response
allSpace = min(spaceIndex):max(spaceIndex);
allOrient = 1:size(R2c,2);

% now compute sum of response across all neurons in the chosen range
fMRI1 = sum(sum(R1c(:,allOrient,allSpace),3),2);
fMRI2 = sum(sum(R2c(:,allOrient,allSpace),3),2);

% normalize to max of 1
m = max([fMRI1;fMRI2]);
fMRI1 = fMRI1 / m;
fMRI2 = fMRI2 / m;

% plot
smartfig('attModelTest2','reuse');clf;
semilogx(contrasts,fMRI2,contrasts,fMRI1);
%semilogx(contrasts,fMRI2+mean(fMRI1-fMRI2),'r-')
xlabel('Contrast');ylabel('Response');
yaxis([],1.05);

smartfig('attModelTest4','reuse');clf;
subplot(1,2,1);
semilogx(contrasts,R2c(:,j,i),contrasts,R1c(:,j,i));
xlabel('Contrast');ylabel('Response');
subplot(1,2,2);
semilogx(contrasts,R2c(:,j,spaceIndex(1)),contrasts,R1c(:,j,spaceIndex(1)));
xlabel('Contrast');ylabel('Response');



% testing ballon model
pedestals = [0 1.75 3.5 7 14 28 56 84];
seglen = [1 0.6 0.2 0.6 0.4 1.2];
deltaCdistributed = [1.76 3 6 8 15 28 24 20 5];
deltaCattended = [0.88 0.44 0.88 1.7 2 3.5 3.5 5 4];

triallen = sum(seglen);

% get time 
t = 0:0.01:triallen;

% now set neural response
n = zeros(size(t));

thisPedestal = 4;
c1 = pedestals(thisPedestal);
c2 = pedestals(thisPedestal)+deltaCdistributed(thisPedestal);

n((t>sum(seglen(1))) & (t<=sum(seglen(1:2)))) = c1;
n((t>sum(seglen(1:3))) & (t<=sum(seglen(1:4)))) = c2;


% now set flow response
f = zeros(size(t));
f((t<=sum(seglen(1:2)))) = 4;

d = balloonmodel(t,n,'flow=f');

