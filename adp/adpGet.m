function val = adpGetD(d,param,varargin)
%
% Getaway function to request values and calculations from the dtaa
% structures used for the adaptation adn attention project.
%
% INPUTS:
%  d - The data structure:
%       d1 = 
%        amplitude: [1x32 double]
%     amplitudeSTE: [1x32 double]
%    amplitudeType: 'max'
%            covar: [32x32 double]
%           deconv: [1x1 struct]
%             ehdr: [32x26 double]
%          fitType: 'glm'
%               r2: 0.0502
%             time: [1x26 double]
%
%  param    - is a string indicating the type of output requested. See
%             below. 
%  varargin - it is possble to pass other parameters that can be
%             used to compute returns using the dat ain d1.
% 
% OUTPUS:
%   val - the value requested by param
%
% EXAMPLE:
%   p = adp_defaultDataFolder;
%   load(fullfile(p,'FM_A0_29-Dec-2008_v1_0_7roi.mat'))
%   amplitudes = adpGet(d1,'crfamp')
%   semilogx(amplitudes)
%
% Franco Pestilli (c) Stanford University 2012

% Format the input parameters.
param = lower(strrep(param,' ',''));
val = [];

switch param
  case {'crfamplitude','crfamp'}
    val = d.amplitude;
    val = reshape(val,8,4);
    
  case {'crfamplitudeste'}
    val = d.amplitudeSTE;
    loval = reshape(val,8,4);

  case {'hrfsfit','ehdrfit'}
    val = d.ehdr;
   
  case {'hrfstime','ehdrtime'}
    val = d.time;
 
  case {'hrfsdata','ehdrdata'}
    val = d.deconv.ehdr;
    
  case {'hrfsdataste','ehdrdataste'}
    val = d.deconv.ehdrste;
    
  case {'amplitudetype'}
    val = d.amplitudeType;
    
  case {'pedestals'}
    val = [0 0.0175 0.035 .07 .14 .28 .56 .84];
    
  otherwise
    help('adpGet')
    fprintf('[adpGet] Unknown parameter << %s >>...\n',param);
    return
end