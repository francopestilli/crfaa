
titleString = 'Figure 7C (Treue & Martinez-Trujillo, 1999)';
stimWidth = 5; 
AxWidth = 5;
AthetaWidth = 45;
Apeak = 5;

% Sampling of space and orientation
x = [-200:200];
theta = [-180:180]';

% Make stimuli 
stimCenter1 = 93;
stimCenter2 = 107;
att_away_loc = -100;
RF_center = round(mean([stimCenter1,stimCenter2]));

pair_att_vars  = [];
pair_att_nulls = [];
pair_att_aways = [];
Var_att_vars    = [];
Null_att_nulls  = [];
Var_att_aways   = [];

% Set contrast to 1
contrast = 1;  

% Pick neuron to record
j = find(theta==0);
i = find(x==RF_center);

orientations = linspace(-180,180,numOrientations);
for stimOrientation1 = orientations

  stimOrientation2 = 180;

  stim1 = contrast * makeGaussian(theta,stimOrientation1,1,1) * makeGaussian(x,stimCenter1,stimWidth,1);
  stim2 = contrast * makeGaussian(theta,stimOrientation2,1,1) * makeGaussian(x,stimCenter2,stimWidth,1);
  pair = stim1 + stim2;

  % Population response when attending stim 1 (varying stim)
  Pair_resp_att_var_pop = attentionModel(x,theta,pair,'Apeak',Apeak,...
    'Ax',stimCenter1,'AxWidth',AxWidth,...
    'Atheta',stimOrientation1,'AthetaWidth',AthetaWidth);

  % Population response when attending stim 2 (null stim)
  Pair_resp_att_null_pop = attentionModel(x,theta,pair,'Apeak',Apeak,...
    'Ax',stimCenter2,'AxWidth',AxWidth,...
    'Atheta',stimOrientation2,'AthetaWidth',AthetaWidth);

  % Population response when attending to fixation point
  Pair_resp_att_away_pop = attentionModel(x,theta,pair,'Apeak',Apeak,...
    'Ax',att_away_loc,'AxWidth',AxWidth,...
    'Atheta',NaN);

  % Population response, attention to var presented alone
  Var_att_var_pop = attentionModel(x,theta,stim1,'Apeak',Apeak,...
    'Ax',stimCenter1,'AxWidth',AxWidth,...
    'Atheta',stimOrientation1,'AthetaWidth',AthetaWidth);

  % Population response, attention to null presented alone
  Null_att_null_pop = attentionModel(x,theta,stim2,'Apeak',Apeak,...
    'Ax',stimCenter2,'AxWidth',AxWidth,...
    'Atheta',stimOrientation2,'AthetaWidth',AthetaWidth);

  % Population response, attention away, var alone
  Var_att_away_pop = attentionModel(x,theta,stim1,'Apeak',Apeak,...
    'Ax',att_away_loc,'AxWidth',AxWidth,...
    'Atheta',NaN);

  pair_att_var    = Pair_resp_att_var_pop(j,i);
  pair_att_null   = Pair_resp_att_null_pop(j,i);
  pair_att_away   = Pair_resp_att_away_pop(j,i);
  Var_att_var     = Var_att_var_pop(j,i);
  Null_att_null   = Null_att_null_pop(j,i);
  Var_att_away    = Var_att_away_pop(j,i);

  pair_att_vars   = [pair_att_vars, pair_att_var];
  pair_att_nulls  = [pair_att_nulls, pair_att_null];
  pair_att_aways  = [pair_att_aways, pair_att_away];
  Var_att_vars    = [Var_att_vars, Var_att_var];
  Null_att_nulls  = [Null_att_nulls, Null_att_null];
  Var_att_aways   = [Var_att_aways, Var_att_away];

end

figure; clf;

plot(orientations,pair_att_vars,orientations,pair_att_nulls,orientations,pair_att_aways,orientations,Var_att_vars,orientations,Null_att_nulls,orientations,Var_att_aways);
xlim([-180 180]);
legend('Pair Var','Pair Null','Pair Away','Var var','Null Null','Var away');
title(titleString);
drawnow
