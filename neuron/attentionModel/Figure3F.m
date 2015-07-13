
titleString = 'Figure 3F (Williford & Maunsell, 2007)';
stimWidth = 7;
AxWidth = 7;
baselineMod = 5e-7;
baselineUnmod = 0;
cRange = [1e-5 1];

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
  % Population response when attending stim 1
  R1 = attentionModel(x,theta,stim,'Ax',stimCenter1,'AxWidth',AxWidth,...
    'baselineMod',baselineMod,'baselineUnmod',baselineUnmod);
  % Population response when attending stim 2
  R2 = attentionModel(x,theta,stim,'Ax',stimCenter2,'AxWidth',AxWidth,...
    'baselineMod',baselineMod,'baselineUnmod',baselineUnmod);
  attCRF(c) = R1(j,i);
  unattCRF(c) = R2(j,i);
end

figure; clf;
subplot(1,2,1)
semilogx(contrasts,unattCRF,contrasts,attCRF);
ylim([0,25]);
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
