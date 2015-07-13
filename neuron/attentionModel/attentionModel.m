function R = attentionModel(x,theta,stimulus,varargin)
% 
% R = attentionModel(x,theta,stimulus,[param1],[value1],[param2],[value2],...,[paramN],[valueN])
%
% Required arguments
% x is a row vector of spatial coordinates
% theta is a column vector of feature/orientation coordinates
% stimulus is NxM where N is the length of theta and M is the length of x
% 
% Optional parameters are passed as string/value pairs. If any of them are
% not specified then default values are used. Valid parameters are as
% follows.
%
% ExWidth: spatial spread of stimulation field
% EthetaWidth: feature/orientation tuning width of stimulation field
% IxWidth: spatial spread of suppressive field
% IthetaWidth: feature/orientation extent/width of suppressive field
% Ax: spatial center of attention field
% Atheta: feature/orientation center of attention field
% AxWidth: spatial extent/width  attention field
% AthetaWidth: feature/orientation extent/width  attention field
% Apeak: peak amplitude of attention field
% Abase: baseline of attention field for unattended locations/features
% Ashape: either 'oval' or 'cross'
% sigma: constant that determines the semi-saturation contrast
% baselineMod: amount of baseline added to stimulus drive
% baselineUnmod: amount of baseline added after normalization
% showActivityMaps: if non-zero, then display activity maps
% showModelParameters: if non-zero, then display stimulus, stimulation
%    field, suppressive field, and attention field.
%
% If Ax or Atheta are NaN or not specified then attention is spread evenly
% such that attnGain = 1 (a constant) for all spatial positions or
% features/orientations, respectively.
%
% Returns the population response (R), same size as stimulus, for neurons
% with receptive fields centered at each spatial position and tuned to each
% feature/orientation.

% Parse varargin to get parameters and values
for index = 1:2:length(varargin)
  field = varargin{index};
  val = varargin{index+1};
  switch field
    case 'ExWidth'
      ExWidth = val;
    case 'EthetaWidth'
      EthetaWidth = val;
    case 'IxWidth'
      IxWidth = val;
    case 'IthetaWidth'
      IthetaWidth = val;
    case 'Ax'
      Ax = val;
    case 'Atheta'
      Atheta = val;
    case 'AxWidth'
      AxWidth = val;
    case 'AthetaWidth'
      AthetaWidth = val;
    case 'Apeak'
      Apeak = val;
    case 'Abase'
      Abase = val; 
    case 'Ashape'
      Ashape = val;
    case 'sigma'
      sigma = val;
    case 'baselineMod'
      baselineMod = val;
    case 'baselineUnmod'
      baselineUnmod = val;
    case 'showActivityMaps'
      showActivityMaps = val;
    case 'showModelParameters'
      showModelParameters = val;
    otherwise
      warning(['attentionModel: invalid parameter: ',field]);
  end
end

% Choose default values for unspecified parameters
if notDefined('ExWidth')
  ExWidth = 5;
end
if notDefined('EthetaWidth')
  EthetaWidth = 60;
end
if notDefined('IxWidth')
  IxWidth = 20;
end
if notDefined('IthetaWidth')
  IthetaWidth = 360;
end
if notDefined('Ax')
  Ax = NaN;
end
if notDefined('Atheta')
  Atheta = NaN;
end
if notDefined('AxWidth')
  AxWidth = ExWidth;
end
if notDefined('AthetaWidth')
  AthetaWidth = EthetaWidth;
end
if notDefined('Apeak')
  Apeak =  2;
end
if notDefined('Abase')
  Abase = 1;
end
if notDefined('Ashape')
  Ashape = 'oval';
end
if notDefined('sigma')
  sigma = 1e-6;
end
if notDefined('baselineMod')
  baselineMod = 0;
end
if notDefined('baselineUnmod')
  baselineUnmod = 0;
end
if notDefined('showActivityMaps')
  showActivityMaps = 0;
end
if notDefined('showModelParameters')
  showModelParameters = 0;
end

% Stimulation field and suppressive field
ExKernel = makeGaussian(x,0,ExWidth);  
IxKernel = makeGaussian(x,0,IxWidth);  
EthetaKernel = makeGaussian(theta,0,EthetaWidth);
IthetaKernel = makeGaussian(theta,0,IthetaWidth);

% Attention field
if isnan(Ax) & isnan(Atheta)
  attnGain = ones(size(stimulus));
else
  if isnan(Ax)
    attnGainX = ones(size(x));
  else
    attnGainX = makeGaussian(x,Ax,AxWidth,1);
    if strcmp(Ashape,'cross')
      attnGainX =  (Apeak-Abase)*attnGainX + Abase;
    end
  end
  if isnan(Atheta)
    attnGainTheta = ones(size(theta));
    Atheta = 0;
  else
    attnGainTheta = makeGaussian(theta,0,AthetaWidth,1);
    if strcmp(Ashape,'cross')
      attnGainTheta = (Apeak-Abase)*attnGainTheta + Abase;
    end
  end
  impulse = (theta == Atheta);
  tmp = impulse * attnGainX;
  attnGain = conv2sepYcirc(tmp,[1],attnGainTheta);
  attnGain = (Apeak-Abase)*attnGain + Abase;
end

% Stimuulus drive
Eraw = conv2sepYcirc(stimulus,ExKernel,EthetaKernel) + baselineMod;
Emax = max(Eraw(:));
E = attnGain .* Eraw;

% Suppressive drive
I = conv2sepYcirc(E,IxKernel,IthetaKernel);
Imax = max(I(:));

% Normalization
R = E ./ (I + sigma) + baselineUnmod;

Rmax = max(R(:));

if showModelParameters == 1
  figure(1); clf;
  subplot(2,2,1); plot(x,ExKernel,x,-IxKernel);
  grid on;
  subplot(2,2,2); plot(theta,EthetaKernel,theta,-IthetaKernel);
  grid on;
  subplot(2,2,3); plot(x,attnGainX);
  grid on;
  subplot(2,2,4); plot(theta,attnGainTheta);
  grid on;
  drawnow
end

if showActivityMaps == 1
  figure(2); clf;
  subplot(3,2,1);
  imshow(stimulus,[0,1]);
  xlabel('Space');
  ylabel('Orientation');
  title ('Stimulus');
  subplot(3,2,2);
  imshow(Eraw,[0,Emax]);
  xlabel('Receptive field center');
  ylabel('Orientation preference');
  title ('Stimulus drive');
  subplot(3,2,3);
  imshow(attnGain,[0,max(max(attnGain))]);
  xlabel('Receptive field center');
  ylabel('Orientation preference');
  title ('Attention field');
  subplot(3,2,4);
  imshow(I,[0,Imax]);
  xlabel('Receptive field center');
  ylabel('Orientation preference');
  title ('Suppressive drive');
  subplot(3,2,5);
  imshow(R,[0,Rmax]);
  xlabel('Receptive field center');
  ylabel('Orientation preference');
  title ('Population response');
  drawnow
end

return


% Test/debug

stimWidth = 5; 
AxWidth = 30;
Atheta = 0;
AthetaWidth = 60;

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
stim = stim1 + stim2;
% Attending stim 1
R1 = attentionModel(x,theta,stim,'Ax',stimCenter1,'AxWidth',AxWidth,...
  'showActivityMaps',1,'showModelParameters',1);
% Attend orientation
R2 = attentionModel(x,theta,stim,'Atheta',Atheta,'AthetaWidth',AthetaWidth,...
  'showActivityMaps',1,'showModelParameters',1);
% Attending stim 1 and orientation, oval
R3 = attentionModel(x,theta,stim,...
  'Ax',stimCenter1,'AxWidth',AxWidth,...
  'Atheta',Atheta,'AthetaWidth',AthetaWidth,...
  'showActivityMaps',1,'showModelParameters',1);
% Attending stim 1 and orientation, cross
R4 = attentionModel(x,theta,stim,'Ashape','cross',...
  'Ax',stimCenter1,'AxWidth',AxWidth,...
  'Atheta',Atheta,'AthetaWidth',AthetaWidth,...
  'showActivityMaps',1,'showModelParameters',1);


