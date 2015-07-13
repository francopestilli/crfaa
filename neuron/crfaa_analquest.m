function results = crfaa_analquest(filename)

load(filename);

indexes = size(stimulus.quest.q);

if length(indexes)==2
    for hh = 1:indexes(1) % pedestal
        for jj = 1:indexes(2) % attentional condition
            if ~isempty(stimulus.quest.q{hh,jj})
                qthreshold(hh,jj) = QuestMean(stimulus.quest.q{hh,jj});
                qthSD(hh,jj) = QuestSd(stimulus.quest.q{hh,jj});
                bEstimate(hh,jj) = QuestBetaAnalysis(stimulus.quest.q{hh,jj});
            end
        end
    end
    % making errorbars:
    se = qthSD;
%     eBarsLow = (squeeze(qthreshold(:,1))-se(:,1));
%     eBarsHi = (squeeze(qthreshold(:,1))+se(:,1));
%     eBars{1} = [eBarsLow; flipud(eBarsHi)];
% 
%     eBarsLow = squeeze(qthreshold(:,2))-se(:,2);
%     eBarsHi = squeeze(qthreshold(:,2))+se(:,2);
%     eBars{2} = [eBarsLow; flipud(eBarsHi)];
    s = 1;
else
    for hh = 1:indexes(1) % pedestal
        for jj = 1:indexes(2) % attentional condition
            for ii = 1:indexes(3) % number of quest
                if ~isempty(stimulus.quest.q{hh,jj,ii})
                    qthreshold(hh,jj,ii) = QuestMean(stimulus.quest.q{hh,jj,ii});
                    qthSD(hh,jj,ii) = QuestSd(stimulus.quest.q{hh,jj,ii});
                    bEstimate(hh,jj,ii) = QuestBetaAnalysis(stimulus.quest.q{hh,jj,ii});
                end
            end
        end
    end
    % making errorbars:
    s = size(qthreshold);
    se = sum(qthSD.^2,3)./s(3);
%     eBarsLow = [(squeeze(qthreshold(:,1,:))-repmat(se(:,1),1,s(3)))];
%     eBarsHi = [(squeeze(qthreshold(:,1,:))+repmat(se(:,1),1,s(3)))];
%     eBars{1} = [eBarsLow; flipud(eBarsHi)];
% 
%     eBarsLow = [(squeeze(qthreshold(:,2,:))-repmat(se(:,2),1,s(3)))];
%     eBarsHi = [(squeeze(qthreshold(:,2,:))+repmat(se(:,2),1,s(3)))];
%     eBars{2} = [eBarsLow; flipud(eBarsHi)];
    s = s(3);
end

if stimulus.pedestals(1) == 0;
    stimulus.pedestals(1) = .002;
end

results.meanThreshold = qthreshold;
results.thresholdSD = qthSD;
results.betaAnalysis = bEstimate;
results.pedestals = stimulus.pedestals;
resutlts.quests = stimulus.quest;

%%  plot:

fv = figure(1);clf;
plotColor = {'r.', 'k.'};
plotColorMedian = {'ro-', 'ks--'};
lightRed = [.99 .4 .1];
gray = [.75 .75 .75];

% loglog(stimulus.pedestals(1),10.^median(squeeze(qthreshold(1,1,1)),2),plotColorMedian{1});
% hold on
% for ii = 1:s
%     pt2(ii) = patch([stimulus.pedestals(:)',fliplr(stimulus.pedestals(:)')], ...
%         10.^squeeze(eBars{2}(:,ii)) ,gray);
% end
% for ii = 1:s
%     pt1(ii) = patch([stimulus.pedestals(:)',fliplr(stimulus.pedestals(:)')], ...
%         10.^squeeze(eBars{1}(:,ii)) ,lightRed);
% end

loglog(stimulus.pedestals(:),10.^median(squeeze(qthreshold(:,1,:)),2),plotColorMedian{1}, ...
    stimulus.pedestals(:),10.^median(squeeze(qthreshold(:,2,:)),2),plotColorMedian{2},...
    stimulus.pedestals(:),10.^squeeze(qthreshold(:,1,:)),plotColor{1}, ...
    stimulus.pedestals(:),10.^squeeze(qthreshold(:,2,:)),plotColor{2});

axis([min(stimulus.pedestals)-.0001, 1 .0005 1]);
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

set(fv,'NumberTitle','off','MenuBar','none');

xlabel('Pedestal contrast','FontName','Arial','FontSize',Fsize);
ylabel('Delta contrast','FontName','Arial','FontSize',Fsize);
title([filename, ' - TvC (Number of quests = ', num2str(s),')'],'FontName','Arial','FontSize',Fsize)

% set(pt1(:),'FaceAlpha', .1,'EdgeColor', lightRed, 'EdgeAlpha', .1);
% set(pt2(:),'FaceAlpha', .1,'EdgeColor', gray, 'EdgeAlpha', .1);

scrsz = get(0,'ScreenSize');
set(fv,'Position',[1 scrsz(4) scrsz(3)/2 scrsz(4)]);
ah = get(fv,'Child'); % getting axis
lh = get(ah,'Child'); % getting points and lines
set(ah,'FontName','Arial','FontSize',Fsize,'PlotBoxAspectRatio',[1 1 1], ...
    'LineWidth',thinLine,'TickLength',[0.025 .01], ...
    'XColor',XYColor,'XTick',[.0001 0.001 0.01 0.1 1],'XTickLabel',[0.0001 0.001 0.01 0.1 1], ...
    'YColor',XYColor,'Box','off', ...
    'YTick',[0.0001 0.001 0.01 0.1 1],'YTickLabel',[0.0001 0.001 0.01 0.1 1]);

% attended
foh = findobj(lh,'Marker','o');
set(foh,'MarkerSize',markerSize,'MarkerFaceColor',myRED,'MarkerEdgeColor',myRED)

%% neutral
foh = findobj(lh,'Marker','s');
set(foh,'MarkerSize',markerSize,'MarkerFaceColor','k','MarkerEdgeColor','k')

