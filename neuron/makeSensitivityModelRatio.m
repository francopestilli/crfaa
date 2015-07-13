function makeSensitivityModelRatio(varargin)
% function makeSensitivityModelRatio(varargin)
%
% this function makes a ration of 
% sigma attended/sigma distributed
% as returned by the sensiivity model.
% it plots a bar graph to be used for fig 6.
% 
% default data folder is:
% ~/data/riken/crfaa/
%
% default data name:
% 03-Sep-2009CompareNoiseFit_nkn2.mat
%
% data is the result of the noise model comparison fit
%  makenoisestatisticsfigs2_r2_compareNoise.m 
%
% noise{iFit,iObserver, iArea}.k;
%
% where 'k' is sigma. 
% iFit = 1 --> D2D
% iFit = 2 --> D2A
% iFit = 3 --> A2A
% iObserver = 1 --> average
% 
% franco pestilli 2010/01/30

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dataFolder = [];
dataFile   = [];
noise      = [];
figname    = [];
figFolder  = [];
getArgs(varargin,{'noise',[],...
        'dataFolder','~/data/riken/crfaa/', ...
        'dataFile','03-Sep-2009CompareNoiseFit_nkn2.mat', ...
        'figname','ratios', ...
        'figFolder', 'figures_modelRatios'});
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
% load the noise file
file2load = fullfile(dataFolder,dataFile);
noise = load(file2load,'noise');
c = 0; % dummy counter for the cue condition
for iFit = [1 3] % fit 1 is D2D, fit 2 is D2A, fit 3 is A2A
 c = c + 1;
 for iObserver = 1:size(noise.noise,2)
  for iArea = 1:size(noise.noise,3)
   sigma(c, iObserver, iArea) = noise.noise{iFit,iObserver, iArea}.k;
   offset(c, iObserver, iArea)= noise.noise{iFit,iObserver, iArea}.responseOffset;
  end
 end
end

% NOW: size(sigma,1) = 2, D & A
%      size(sigma,2) = 4, average, obs1, obs2 and obs3
%      size(sigma,3) = 4, V1, V2, V3 and V4 
%      

distributed = 1;
focal       = 2;

%%%%%%%%%%%%%%%%%%%%%%%%%
% (1) do sigma ratio plot

% get the sensitivity model sigma ratios:
sigmaFocal       = squeeze(sigma(focal,:,:));
sigmaDistributed = squeeze(sigma(distributed,:,:));
sigmaRatio       = sigmaDistributed./sigmaFocal;
% NOW: size(sigmaRation,1) = 4 --> observers
%      size(sigmaRation,2) = 4 --> visual areas

% some things we need for doing small n correction
averageOver = 1; % observers
n = size(sigmaRatio,averageOver);
smallNcorrection = sqrt(2/(n-1))*gamma(n/2)/gamma((n-1)/2);

% now calculate mean and standard error
ratioMean = [median(sigmaRatio,averageOver)];
meanOfRatios = mean(ratioMean);

% calculate ratio STE using small n corection
ratioSTE = [std(sigmaRatio,[],averageOver)/smallNcorrection]/sqrt(n);

% display the bar graph
fignameS = [figname,'_sigma'];
h = smartfig(fignameS,'reuse');
mybar(ratioMean,'yError',ratioSTE,'groupLabels',{'Sensitivity'}, ...
                'withinGroupLabels',{'V1','V2','V3','V4'},       ...
                'yAxisMin=0','yAxisMax=8',                       ...
                'yLabelText=Ratio of sigmaDistributed to sigmaFocal',...
                'xLabelText=Model','dispValues=1',sprintf('hline=%2.2f',meanOfRatios),'hLineStyle=--');

% save the figure
savefig(h,'figName',fignameS,'figDir',figFolder,'defaultDataFolder',dataFolder(1:end-1),'verbose=1');


%%%%%%%%%%%%%%%%%%%%%%%%%%
% (2) do offset difference plot

% get the sensitivity model sigma ratios:
offsetFocal       = squeeze(offset(focal,:,:));
offsetDistributed = squeeze(offset(distributed,:,:));
offsetDiff       = offsetFocal - offsetDistributed;

% some things we need for doing small n correction
averageOver = 1; % observers
n = size(offsetDiff,averageOver);
smallNcorrection = sqrt(2/(n-1))*gamma(n/2)/gamma((n-1)/2);

% now calculate mean and standard error
diffMean = [median(offsetDiff,averageOver)];
meanOfDiff = mean(diffMean);

% calculate ratio STE using small n corection
ratioSTE = [std(offsetDiff,[],averageOver)/smallNcorrection]/sqrt(n);

% display the bar graph
fignameO = [figname,'_offset'];
h = smartfig(fignameO,'reuse');
mybar(diffMean,'yError',ratioSTE,'groupLabels',{'Sensitivity'}, ...
                'withinGroupLabels',{'V1','V2','V3','V4'},       ...
                'yAxisMin=0','yAxisMax=0.25',                       ...
                'yLabelText=Difference btween offset_f to offset_d',...
                'xLabelText=Model','dispValues=1',sprintf('hline=%2.2f',meanOfDiff),'hLineStyle=--');

% save the figure
savefig(h,'figName',fignameO,'figDir',figFolder,'defaultDataFolder',dataFolder(1:end-1),'verbose=1');
keyboard
