function results = crfaa_analDetection(observer, tc , adaptconditions)
% results = crfaa_analDetection(observer, tc , adaptconditions)
%
% makes a plot of the results.
% tc is 0 < tc < 1 choses the max sd for a threshold to be used.
%
% franco pestili may, 2008


% d to data dir:
datadir = '/Users/frakkopesto/data/nyu/attadp/data_psychophysics/';
cd(datadir);



for ii = 1:length(adaptconditions)
 % check if the current observers dir exists,
 datadirs{ii} = [observer,'_A',adaptconditions{ii}];%,'_gaps'];

 if ~isdir(datadirs{ii})
  disp(sprintf('[crfaa_analDetection] Could not find observer dir : %s',datadirs{ii}));
  return;
 end

 % get all files in each directory:
 ac{ii} = dir([datadir,datadirs{ii}]);
 analyzeFile{ii} = ac{ii}(end).name;
end

keyboard

% construct the set of file names we would like to analyze (the last one in
% each adapt condition):
removeme = find(strcmp(analyzeFile,'..'));
if ~isempty(removeme)
 useme = ~([1:length(adaptconditions) == removeme]);
 for ii = find(useme)
  filename{ii} = analyzeFile{ii};
 end
else
 filename = analyzeFile;
 useme = [1:length(adaptconditions)] == [1:length(adaptconditions)];
end

% get the quest data out of the files:
for ff = 1:length(filename)
 filetoload{ff} = deblank([datadirs{ff},'/',filename{ff}]);
 dataVars{ff} = load(filetoload{ff});
 % disregarding adaptation conditions with no stim files:
 quest{ff} = squeeze(dataVars{ff}.stimulus.quest.q);
end

% accept only quest with more than 39 trials:
for ac = 1:length(quest) % adaptation condition
 for aa = 1:size(quest{ac},1) % attentional condition
  for qn = 1:size(quest{ac},2) % quest number
   if ~isempty(quest{ac}{aa,qn})
    if length(quest{ac}{aa,qn}.intensity) >= 39
     lastQuest(ac, aa) = qn;
    end
   end
  end
 end
end

% analyze the quest:
qthreshold = nan.*ones(length(quest),size(quest{1},1),max(max(lastQuest)));
qthSD = qthreshold; bEstimate = qthreshold;
for ac = 1:length(quest) % adaptation condition
 for aa = 1:size(quest{1},1) % attentional condition
  for qn = 1:lastQuest(ac,aa) % quest number
   qthreshold(ac,aa,qn) = QuestMean(quest{ac}{aa,qn});
   qthSD(ac,aa,qn) = QuestSd(quest{ac}{aa,qn});
%    bEstimate(ac,aa,qn) = QuestBetaAnalysis(quest{ac}{aa,qn});
  end
 end
end

results.thrZScriterion = tc;
results.thresholds = qthreshold;
results.thresholdSD = qthSD;
results.zscores = zscore(results.thresholds);

% select only the threshold that are within the accepted range:
results.usedThresholds = zeros(size(results.thresholds));

for ac = 1:length(quest) % adaptation condition
 for aa = 1:size(quest{1},1) % attentional condition
  for qn = 1:lastQuest(ac,aa) % quest number
   if abs(results.zscores(ac,aa,qn)) <= results.thrZScriterion
    results.usedThresholds(ac,aa,qn) = 1;
   end
  end
 end
end
results.usedThresholds = logical(results.usedThresholds);

% average only across the accepted thresholds:
results.meanThreshold = nan.*ones(size(results.usedThresholds));
for ac = 1:length(quest) % adaptation condition
 for aa = 1:size(quest{1},1) % attentional condition
  results.meanThreshold(ac,aa,:) = nanmean(qthreshold(ac,aa,results.usedThresholds(ac,aa,:)),3);
  results.sdThreshold(ac,aa,:) = nanstd(qthreshold(ac,aa,results.usedThresholds(ac,aa,:)),[],3);
 end
end

results.betaAnalysis = bEstimate;
results.pedestal = dataVars{1}.stimulus.pedestals;
results.quests = quest;
results.usedQuests = lastQuest;
results.usedAdaptation = useme;
results.usedFiles = filetoload;

%% %%%%%%%%%%%%%%%%
%  plot:

fv = figure(1);clf;
plotColor = {'r.', 'k.'};
plotColorMedian = {'ro-', 'ks--'};
lightRed = [.99 .4 .1];
gray = [.75 .75 .75];

x = 1:length(quest);
at = 100.*(10.^results.meanThreshold(:,1));
aLowerSD = 100.*10.^(results.meanThreshold(:,1) - squeeze(results.sdThreshold(:,1,:)));
aUpperSD = 100.*10.^(results.meanThreshold(:,1) + squeeze(results.sdThreshold(:,1,:)));

ut = 100.*10.^results.meanThreshold(:,2);
uLowerSD = 100.*10.^(results.meanThreshold(:,2) - squeeze(results.sdThreshold(:,2,:)));
uUpperSD = 100.*10.^(results.meanThreshold(:,2) + squeeze(results.sdThreshold(:,2,:)));

% do the actual plot
loglog(x,at,plotColorMedian{1}, x,ut,plotColorMedian{2})%,...
hold on;
% patch([x fliplr(x)], [aLowerSD;aUpperSD],lightRed);
% patch([x fliplr(x)], [uLowerSD;uUpperSD],gray);
loglog(x,at,plotColorMedian{1}, x,ut,plotColorMedian{2});
loglog([repmat(x,2)'],[aLowerSD;aUpperSD],plotColor{1});
loglog([repmat(x,2)'],[uLowerSD;uUpperSD],plotColor{2});

% loglog(1:length(quest),100.*(10.^results.meanThreshold(:,1)),plotColorMedian{1}, ...
%     1:length(quest),100.*(10.^results.meanThreshold(:,2)),plotColorMedian{2},...
%     [repmat(1:length(quest),size(squeeze(qthreshold(:,2,:)),2),1)],100.*(10.^squeeze(qthreshold(:,1,:)))',plotColor{1}, ...
%     [repmat(1:length(quest),size(squeeze(qthreshold(:,2,:)),2),1)],100.*(10.^squeeze(qthreshold(:,2,:)))',plotColor{2});

axis([0, length(quest)+1 .1 100]);
axis('square');

%%%%%%%%%%%%%%%%%%%%
% figure formatting:
% colors:
myRED = [.7 .1 .1];
myBLACK = [0 0 0 ];
myGREEN = [.1 .5 .1];
myGRAY = [.7 .7 .7];
myLightRED = [.97 .79 .69];
myLightGREEN = [.69 .97 .79];

lineThick = 3;
thinLine = 2;
markerSize = 20;
Fsize = 40;
XYColor = [.2 .2 .2];

xlabel('Adapter contrast (%)','FontName','Arial','FontSize',Fsize);
ylabel('Detection threshold (%)','FontName','Arial','FontSize',Fsize);
title([observer, ' - TvC'],'FontName','Arial','FontSize',Fsize)

scrsz = get(0,'ScreenSize');
set(fv,'Position',[1 scrsz(4) scrsz(3)/2 scrsz(4)]);
ah = get(fv,'Child'); % getting axis
lh = get(ah,'Child'); % getting points and lines
for ii = find(useme)
 adpLabel(ii) = str2double(adaptconditions{ii});
end

set(ah,'FontName','Arial','FontSize',Fsize,'PlotBoxAspectRatio',[1 1 1], ...
 'LineWidth',thinLine,'TickLength',[0.025 .01], ...
 'XColor',XYColor,'XTick',[1:length(quest)],'XTickLabel',adpLabel, ...
 'YColor',XYColor,'Box','off', ...
 'YTick',[.001 .01 .1 1 10],'YTickLabel',[.001 .01 .1 1 10]);

% attended
foh = findobj(lh,'Marker','o');
set(foh,'MarkerSize',markerSize,'MarkerFaceColor',myRED,'MarkerEdgeColor',myRED)

%% neutral
foh = findobj(lh,'Marker','s');
set(foh,'MarkerSize',markerSize,'MarkerFaceColor','k','MarkerEdgeColor','k')

