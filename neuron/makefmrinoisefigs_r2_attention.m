function makefmrinoisefigs_r2_attention(varargin)
%
% function makefmrinoisefigs_r2_attention({plotType},{observer},{visualarea},{adaptation})
%
% this version has been developed while at riken. it uses one roi that is the
% combination of data across v1, v2, v3.
%
% this function loads the results of fitTimecourse.m
% run with the 'std' option and the 'group' option.
% it generates a bar graph of the noise vlues computed over the
% fmri responses.
%
% e.g.,
% makefmrinoisefigs('plotType',[1 2],'observer',{'jg' 'fp' 'fm'},'visualarea',{'v1'},'adaptation',{0})
%
% version r2 calls the data computed on several rois with different r2
% values. it does not load the older files anymore.
% use makefmrinoisefigs_r.m for the.
%
%
% franco pestilli 2008.11.12

% global figureHandle;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
observer=[];       %                                                         %
visualarea=[];     %                                                         %
adaptation=[];     %                                                         %
doMeanObserver=[]; % average across observers or not                         %
iscombined=[];     % if rois are the combined result of v1-v3                %
averageD=[];       % averages the Dnt and Dt STDs                            %
normalize=[];      % nomralizes the STD by the STD of the attended condition %
savePlots=[];      % 1 saves the plot in the current dierectory              %
dispFit=[];        % if 1 the plot is displayed                              %
saveFig=[];        % if 1 saves a figure                                     %
day=[];            % day of analysis                                         %
doOnePlot=[];      % do one plot with attention and adaptation results       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set default variables:
getArgs(varargin,{...
 'observer',{'fm' 'fp' 'jg'}, ...
 'visualarea',{'v1' 'v2' 'v3' 'v4'}, ...{'05' '06' '07' '075'}, OR  {'v1' 'v2' 'v3' 'v4'}
 'adaptation',{0},   ...
 'iscombined', 0,    ...
 'averageD',1,       ...
 'normalize',0,      ...
 'doMeanObserver',1, ...
 'day', '15-Feb-2009'... combined roi date: '05-Jul-2009', individual roi date: '15-Feb-2009', ...
 'savePlots',0,      ...
 'dispFit',0,        ...
 'saveFig',1,        ...
 'doOnePlot',1});

if iscombined
 day = '05-Jul-2009';
end

if averageD
 allSTD = nan.*ones(length(adaptation),length(observer),3);
else
 keyboard
 % not implemented now for the average graphs
 allSTD = nan.*ones(size(allfiles,1),size(allfiles,2),4);
end

% set up geenral variable:
defaultDataFolder = '/Users/frakkopesto/data/riken/crfaa/fmridata/crfaa/';

for v = 1:length(visualarea) % loop through each r2 roi
 % loop thorugh each adaptation condition:
 for adpt = 1:length(adaptation)
  for o = 1:length(observer)
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%
   % get the right file name %
   %%%%%%%%%%%%%%%%%%%%%%%%%%%
   if ((adaptation{adpt} == 100) && strcmpi(observer{o},'fm'))
    disp(sprintf('NO adapt-%i condition for observer %s',adaptation{adpt},observer{o}));
    allfiles{v,adpt,o} = 'nofile';
   else
    allfiles{v,adpt,o} = getAllFileNames(observer{o},adaptation{adpt},day,visualarea{v},iscombined);
   end

   % prepare the arrays
   if ((adaptation{adpt} == 100) && strcmpi(observer{o},'fm'))
    disp(sprintf('NO adapt-%i condition for observer %s',adaptation{adpt},observer{o}));
    if averageD
     allSTD(adpt,o,:) = [1 1 1];
    else
     allSTD(adpt,o,:) = [1 1 1 1];
    end
    
   else
    
    %%%%%%%%%%%%%%%%%%%%%%%
    % load the data here: %
    %%%%%%%%%%%%%%%%%%%%%%%
    d1{v,adpt,o} = load(deblank(sprintf('%s',allfiles{v,adpt,o})),'d1');
    
    %%%%%%%%%%%%%%%%%%%%%%
    % average Dnt and Dt %
    %%%%%%%%%%%%%%%%%%%%%%
    if logical(averageD)
     allSTD(adpt,o,1) = d1{v,adpt,o}.d1.groupSTD.amplitudeSTD(1);
     allSTD(adpt,o,2) = mean(d1{v,adpt,o}.d1.groupSTD.amplitudeSTD([2 3]));%sqrt(sum(d1.d1.groupSTD.amplitudeSTD([2 3]).^2))./2;
     allSTD(adpt,o,3) = d1{v,adpt,o}.d1.groupSTD.amplitudeSTD(3);
    else
     allSTD(adpt,o,:) = d1{v,adpt,o}.d1.groupSTD.amplitudeSTD;
    end
   end
  end
 end
 
 
 %%%%%%%%%%%%%%%%%%%%%%%
 % make attention plot %
 %%%%%%%%%%%%%%%%%%%%%%%
 figTag = sprintf('Att_%s',visualarea{v});

 % average across observers and adaptation conditions.
 attNorm_allSTD = zeros(size(allSTD));
 if normalize
  for a = 1:size(allSTD,1) % adaptation
   for o = 1:size(allSTD,2) % observer  
    attNorm_allSTD(a,o,:) = allSTD(a,o,:)./allSTD(a,o,1); % normalizing by attended
   end
  end
 else
  % NOT normalizing by the attended value
  attNorm_allSTD = allSTD;
 end

 if doMeanObserver
  if length(adaptation) == 1
   % not sure about this
   meanNoiseAtt = squeeze(nanmean(attNorm_allSTD,2)); % mean across observers
   steNoiseAtt =  squeeze(nanstd(attNorm_allSTD,[],2))/sqrt(3); % ste acrss observers
  else
   meanNoise = squeeze(nanmedian(attNorm_allSTD,2)); % mean across observers
   steNoiseAtt =  nanstd(meanNoise,[],1)/sqrt(3); % ste acrss observers
   meanNoiseAtt = nanmean(meanNoise,1);
  end
  noiseATTIndexPlot(meanNoiseAtt,steNoiseAtt,saveFig,figTag,averageD,normalize);
 
 else % do not average across observers
   meanNoiseAtt = squeeze(nanmean(attNorm_allSTD,1)); % mean across adaptation
   steNoiseAtt =  squeeze(nanstd(attNorm_allSTD,[],1))./sqrt(3); % ste acrss adaptation

  % plot
  for o = 1:size(allSTD,2)
   figTag1 = sprintf('%s%s%s',figTag,'_O',num2str(o));
   noiseATTIndexPlot(meanNoiseAtt(o,:),steNoiseAtt(o,:),saveFig,figTag1,averageD,normalize);
   figTag1 = [];
  end
 end
end % end across r2 rois

% %%%%%%%%%%%%%%%% END MAIN CALL %%%%%%%%%%%%%%%% %


%%%%%%%%%%%%%%%%%%%%%
% noiseATTIndexPlot %
%%%%%%%%%%%%%%%%%%%%%
function noiseATTIndexPlot(STD,STDste,saveFig,figTag,averageD,normalize)
% function noiseScatterPlot(noise)
% this function is able to plot bar-graphs across visual areas bu tit is
% currently not used.
%
figurename = sprintf('fmriNoise_%s',figTag);
h = smartfig(figurename,'reuse');
set(h,'Name',figurename);

plotInfo = plotSetUp;
if strcmp(figTag(1:3),'Att')

 if strcmpi(figTag(5:end),'v4')
  ymax = 1.75;
 else
  ymax = 1.75;
 end
 
 if averageD % average Dn and Dnt
  label = {'A','D','U'};
  labelG= {'Attention'};
 else
  label = {'A','Dnt','Dt','U'};
 end
 if normalize
  ylabel = sprintf('STD of fMRI response\n (normalized to neutral STD)');
 else
  ylabel = sprintf('STD of fMRI response\n (%s signal change)','%');
 end
else
 keyboard
end


if size(STD,1) > 1
 STD = STD';
end
if size(STDste,1) > 1
 STDste = STDste';
end

mybar(STD, 'yError',STDste, ...
 'withinGroupLabels',label, ...
 'yAxisMin=0.5',sprintf('yAxisMax=%i',ymax), ...
 'groupLabels',labelG, ...
 'yLabelText',ylabel,'dispValues',1)


if saveFig
 set(h,'PaperPosition',[.25 .25 8 10.5]);
 set(h,'PaperOrientation','Portrait');
 
 disp(sprintf('Saving figure %s',figurename));
 figDir = 'fig_attention_fMRInoise';
 defaultDataFolder = '/Users/frakkopesto/data/riken/crfaa/fmridata/crfaa/';
 figDir = savefig(h,'figName',figurename,'defaultDataFolder', defaultDataFolder,'figDir',figDir);
end


%%%%%%%%%%%%%%%%%%%%%
%  getAllFileNames  %
%%%%%%%%%%%%%%%%%%%%%
function filename = getAllFileNames(observer,adaptation,day,visualarea,iscombined)

% check which computer are we using:
[notok c] = system('hostname');
if ~(~notok && strcmp('riken.jp',deblank(c(end-8:end))))
 % check if there is a drive connected which one is it ...
 if isdir('/Volumes/homosacer_data')
  hd = '/Volumes/homosacer_data';
 elseif isdir('/Volumes/homosacer')
  hd = '/Volumes/homosacer';
 elseif isdir('/Volumes/homosacer_backup')
  hd = '/Volumes/homosacer_backup';
 else
  disp('No data hard drive found.')
  keyboard
 end
elseif (~notok && strcmp('riken.jp',deblank(c(end-8:end))))
 hd = '/Volumes/data/riken/crfaa';
end

% display the info regarding the HD being used.
disp(sprintf('Using data in HD: <%s>',hd));

switch lower(observer)
 case {'fm'}
  switch adaptation
   case {0}
    session = 'fm20080209';
   case 28
    session = 'fm20080406';
   otherwise
    disp(sprintf('This condition is NOT set (%s, adp-%i).',observer,adaptation))
    session = [];
  end
  
 case {'fp'}
  switch adaptation
   case {0}
    session = 'fp20071019';
   case {28}
    session = 'fp20080402';
   case {100}
    session = 'fp20080415';
   otherwise
    disp(sprintf('This condition is NOT set (%s, adp-%i).',observer,adaptation))
    keyboard
  end
  
 case {'jg'}
  switch adaptation
   case {0}
    session = 'jg20070919';
   case {28}
    session = 'jg20080402';
   case {100}
    session = 'jg20080414';
   otherwise
    disp(sprintf('This condition is NOT set (%s, adp-%i).',observer,adaptation))
    keyboard
  end
  
 otherwise
  disp(sprintf('Observer (%s) NOT found.',observer))
  keyboard
end

if ~isempty(session)
 thispath = deblank(sprintf('%s/fmridata/crfaa/%s_A%i_testStream/%s',hd,upper(observer),adaptation,session));
 if iscombined
  filename = deblank(sprintf('%s/%s_v1_subTSnoiseSTD_combinedR2%s.mat',thispath,day,lower(visualarea)));
 else
  filename = deblank(sprintf('%s/%sNoiseAnalysis%s.mat',thispath,day,lower(visualarea)));
 end
else
 filename = [];
end
disp(sprintf('Loading: %s',filename))


%%%%%%%%%%%%%
% plotSetUp %
%%%%%%%%%%%%%
function plotInfo = plotSetUp

% basic plot set up:
col = 1; % slightly change the color depending on adaptation
myRed = [col 0 0];
myBlack = [0 0 0];
myBlue = [0 0 1];
myGreen = [0 col 0];
myCyan = [col col 0];

plotInfo.contrast = [.001 .0175 .035 .07 .14 .28 .56 .84];
plotInfo.MarkerFaceColor{1} = myRed;
plotInfo.MarkerFaceColor{2} = myBlack;
plotInfo.MarkerFaceColor{3} = myBlack;
plotInfo.MarkerFaceColor{4} = myGreen;
plotInfo.thisColor = myBlack;
plotInfo.thisSymbol = {'o-' 's-' 's--' '^-'};
plotInfo.XYColor = [.4 .4 .4];
plotInfo.Fsize = 8;
plotInfo.MarkerSize = 10;
plotInfo.MarkerEdgeColor = [0 0 0];
plotInfo.LineWidth = 1;
plotInfo.TickLength = [0.025 .01];
plotInfo.PlotBoxAspectRatio = [1 1 1];
plotInfo.YLim = [0 1];
plotInfo.Xlim = [0 1];
plotInfo.XTicks = [0 .25 .5 .75 1];
plotInfo.YTicks = [0 .25 .5 .75 1];
plotInfo.title = {'v1' 'v2' 'v3' 'v4'};


