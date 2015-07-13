function results = crfaa_analDiscrimination(observer,[scanDir])

% d to data dir:
if exist(scanDir)
    datadir = scanDir;
else
    datadir = '~/fmritools/ATT_ADPT_FILES/CRFAA/data_scanner/';
end
cd(datadir);

% get all files in observer directory:
ac = dir([datadir,observer]);
filename = ac(end).name;
keyboard

% get the quest data out of the files:
dataVars = load(filename);
% disregarding adaptation conditions with no stim files:
quest = squeeze(dataVars.stimulus.quest.q);

% accept only quest with more than 39 trials:
for pedC = 1:size(quest,1) % pedC
    for aa = 1:size(quest,2) % attentional condition
        for qn = 1:size(quest,3)
            if ~isempty(quest{pedC,aa,qn})
                if length(quest{pedC,aa,qn}.intensity) >= 36
                    lastQuest(pedC,aa) = qn;
                end
            end
        end
    end
end

if any(lastQuest<3)
   disp('[crfaa_analDiscrimination] One of the pedestals has less then 3 quests.') 
end

% analyze the quest:
qthreshold = nan.*ones(size(quest,1),size(quest,2),max(max(lastQuest)));
qthSD = qthreshold; bEstimate = qthreshold;
for pedC = 1:size(quest,1) % pedC
    for aa = 1:size(quest,2) % attentional condition
        for qn = 1:lastQuest(pedC) % quest number
            qthreshold(pedC,aa,qn) = QuestMean(quest{pedC,aa,qn});
            qthSD(pedC,aa,qn) = QuestSd(quest{pedC,aa,qn});
            bEstimate(pedC,aa,qn) = QuestBetaAnalysis(quest{pedC,aa,qn});
        end
    end
end


results.meanThreshold = squeeze(nanmean(qthreshold,3));
results.thresholds = qthreshold;
results.thresholdSD = qthSD;
results.betaAnalysis = bEstimate;
results.pedestal = dataVars.stimulus.pedestals;
results.quests = quest;
results.usedQuests = lastQuest;
results.usedFiles = filename;

%%  plot:

fv = figure;clf;
plotColor = {'r.', 'k.'};
plotColorMedian = {'ro-', 'ks--'};
lightRed = [.99 .4 .1];
gray = [.75 .75 .75];

loglog(1:length(quest),100.*(10.^results.meanThreshold(:,1)),plotColorMedian{1}, ...
    1:length(quest),100.*(10.^results.meanThreshold(:,2)),plotColorMedian{2},...
    [repmat(1:length(quest),size(squeeze(qthreshold(:,2,:)),2),1)],100.*(10.^squeeze(qthreshold(:,1,:)))',plotColor{1}, ...
    [repmat(1:length(quest),size(squeeze(qthreshold(:,2,:)),2),1)],100.*(10.^squeeze(qthreshold(:,2,:)))',plotColor{2});

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
thinLine = 1;
markerSize = 16;
Fsize = 30;
XYColor = [.2 .2 .2];

xlabel('Adapter contrast (%)','FontName','Arial','FontSize',Fsize);
ylabel('Detection threshold (%)','FontName','Arial','FontSize',Fsize);
title([observer, ' - TvC'],'FontName','Arial','FontSize',Fsize)

scrsz = get(0,'ScreenSize');
set(fv,'Position',[1 scrsz(4) scrsz(3)/2 scrsz(4)]);
ah = get(fv,'Child'); % getting axis
lh = get(ah,'Child'); % getting points and lines

set(ah,'FontName','Arial','FontSize',Fsize,'PlotBoxAspectRatio',[1 1 1], ...
    'LineWidth',thinLine,'TickLength',[0.025 .01], ...
    'XColor',XYColor,'XTick',[1:length(results.pedestal)],'XTickLabel',results.pedestal, ...
    'YColor',XYColor,'Box','off', ...
    'YTick',[.01 .1 1 10 100],'YTickLabel',[.01 .1 1 10 100]);

% attended
foh = findobj(lh,'Marker','o');
set(foh,'MarkerSize',markerSize,'MarkerFaceColor',myRED,'MarkerEdgeColor',myRED)

%% neutral
foh = findobj(lh,'Marker','s');
set(foh,'MarkerSize',markerSize,'MarkerFaceColor','k','MarkerEdgeColor','k')

