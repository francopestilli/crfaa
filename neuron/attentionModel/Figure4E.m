
titleString = 'Figure 4E (yoked contrast experiment)';
stimWidth = 5;
AxWidth = 5;
AthetaWidth = 20;
Apeak = 5;
cRange = [1e-4 0.1];

% Sampling of space and orientation
x = [-200:200];
theta = [-180:180]';

% Make stimuli
stimCenter1 = 90;
stimOrientation1 = 0;
stimCenter2 = 110;
stimOrientation2 = 180;
stimCenter3 = -90;
stimOrientation3 = 0;
stimCenter4 = -110;
stimOrientation4 = 180;

% Choose neuron with RF centered at midpoint between the two stim intended
% to be in RF:
RF_center = round(mean([stimCenter1, stimCenter2]));

% Stim 1 and 2 in RF
stim1 = makeGaussian(theta,stimOrientation1,1,1) * makeGaussian(x,stimCenter1,stimWidth,1);
stim2 = makeGaussian(theta,stimOrientation2,1,1) * makeGaussian(x,stimCenter2,stimWidth,1);

% Stim 3 and 4 contralateral to RF
stim3 = makeGaussian(theta,stimOrientation1,1,1) * makeGaussian(x,stimCenter3,stimWidth,1);
stim4 = makeGaussian(theta,stimOrientation2,1,1) * makeGaussian(x,stimCenter4,stimWidth,1);

% Pick contrasts
logCRange = log10(cRange);
logContrasts = linspace(logCRange(1),logCRange(2),numContrasts);
contrasts = 10.^logContrasts;

% We are interested in a neuron that prefers orientation 1. Orientation 2
% is its null stimulus.
j = find(theta==stimOrientation1);
i = find(x==RF_center);

attCRF = zeros(size(contrasts));
unattCRF = zeros(size(contrasts));
for c = 1:numContrasts
  stim = contrasts(c) * stim1 + contrasts(c) * stim2 + contrasts(c) * stim3 + contrasts(c) * stim4;
  % Population response when attending preferred stim in RF:
  R1 = attentionModel(x,theta,stim,'Apeak',Apeak,...
    'Ax',stimCenter1,'AxWidth',AxWidth,...
    'Atheta',stimOrientation1,'AthetaWidth',AthetaWidth);
  % Population response when attending null stim in RF:
  R2 = attentionModel(x,theta,stim,'Apeak',Apeak,...
    'Ax',stimCenter2,'AxWidth',AxWidth,...
    'Atheta',stimOrientation2,'AthetaWidth',AthetaWidth);
  attCRF(c) = R1(j,i);
  unattCRF(c) = R2(j,i);
end

figure; clf;
subplot(1,2,1)
semilogx(contrasts,unattCRF,contrasts,attCRF);
ylim([0 7]);
xlim(cRange);
legend('Att non-pref','Att pref');
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
