function plotRT(rt, observer)
% function plotRT(rt,observer)
%
% observer can be 'FM' or 'FP' or 'JG'
%
%
%

% check arguments
if ~any(nargin == 2)
    help plotRT
    return
end

if ~(strcmp(observer,'FM') | strcmp(observer,'FP') | strcmp(observer,'JG'))
    help plotRT
    return
end

if ~iscell(rt)
    rt1{1} = rt;
    clear rt;
    rt = rt1;
    clear rt1;
end

% plot them:
% plot's layout parameters:
XYColor = [.4 .4 .4];
Fsize = 16;
lineThickness = 3;
myRed = [.6 0 0];
myBlack = [0 0 0];
myGreen = [0 .4 0];
minX = .01;
adaptation = {'0' '28' '100'};

fv(1) = figure;
% turn off menu/title etc.
thetitle =[observer, '_RT'];
set(fv(1),'NumberTitle','off','Name',thetitle);
scrsz = get(0,'ScreenSize');
set(fv(1),'Position',[1 scrsz(4) scrsz(3)/2 scrsz(4)]);

numadapt = length(rt);
for i = 1:numadapt

    rt{i}.pedC = 100*rt{i}.pedC;
    rt{i}.pedC(1) = minX;

    subplot(1,numadapt,i),semilogx(rt{i}.pedC,rt{i}.quantileArt(:,2),'r^-','MarkerSize',8,'Color',myRed,'LineWidth',lineThickness);
    hold on
    subplot(1,numadapt,i),semilogx(rt{i}.pedC,rt{i}.quantileDrt(:,2),'ks-','MarkerSize',8,'Color',myBlack,'LineWidth',lineThickness);
    subplot(1,numadapt,i),semilogx(rt{i}.pedC,rt{i}.quantileArt(:,[1 3]),'r-','Color',myRed,'LineWidth',.5);
    subplot(1,numadapt,i),semilogx(rt{i}.pedC,rt{i}.quantileDrt(:,[1 3]),'k-','Color',myBlack,'LineWidth',.5);

    title(sprintf('Reaction time - %s',adaptation{i}),'FontSize',16)
    ylabel(' Reaction time (s) ','FontSize',Fsize);
    xlabel(' Contrast (%) ','FontSize',Fsize);
    set(gca,...
        'FontName','Helvetica','FontSize',Fsize,'PlotBoxAspectRatio',[1 1 1], 'XLim', [minX 100],...
        'LineWidth',1,'TickLength',[0.045 .015], 'XTick', [minX rt{i}.pedC(2:end)], ...
        'XTickLabel', [0 rt{i}.pedC(2:end)],'YLim', [0 1.1], ...
        'XColor',XYColor,'YColor',XYColor,'Box','off');
    drawnow
end
h = legend({'Attended','Distributed'},2);
set(h, 'interpreter','none','Box','off');

disp(sprintf('Saving figure %s',thetitle));
eval(sprintf('print(gcf,''-depsc2'',''-tiff'',''-r300'', ''%s'')', thetitle));

