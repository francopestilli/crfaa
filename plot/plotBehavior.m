function plotBehavior(c, observer)
% function plotBehavior(c, observer)
%
%
% observer can be 'FM' or 'FP' or 'JG'
%
% 2008.10.12 franco pestilli

% check arguments
if ~any(nargin == 2)
    help plotRT
    return
end

if ~(strcmp(observer,'FM') | strcmp(observer,'FP') | strcmp(observer,'JG'))
    help plotRT
    return
end

if ~iscell(c)
    c1{1} = c;
    clear c;
    c = c1;
    clear c1;
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
thetitle =[observer, '_TVC'];
set(fv(1),'NumberTitle','off','Name',thetitle);
scrsz = get(0,'ScreenSize');
set(fv(1),'Position',[1 scrsz(4) scrsz(3)/2 scrsz(4)]);
numadapt = length(c);

for i = 1:numadapt
%% strip out contrast quantiles:
for ii = 1:length(c.CI_ATcontrast)
    % attended contrast
    ac_min(ii)    = 100.*c.CI_ATcontrast{ii}(1);
    ac_fiftyQ(ii) = 100.*c.CI_ATcontrast{ii}(2);
    ac_max(ii)    = 100.*c.CI_ATcontrast{ii}(3);

    % attended delta c
    adc_min(ii)    = 100.*c.CI_ATdeltaC{ii}(1);
    adc_fiftyQ(ii) = 100.*c.CI_ATdeltaC{ii}(2);
    adc_max(ii)    = 100.*c.CI_ATdeltaC{ii}(3);
    
    % distributed contrast
    dc_min(ii)    = 100.*c.CI_DTcontrast{ii}(1);
    dc_fiftyQ(ii) = 100.*c.CI_DTcontrast{ii}(2);
    dc_max(ii)    = 100.*c.CI_DTcontrast{ii}(3);

    % distributed delta c
    ddc_min(ii)    = 100.*c.CI_DTdeltaC{ii}(1);
    ddc_fiftyQ(ii) = 100.*c.CI_DTdeltaC{ii}(2);
    ddc_max(ii)    = 100.*c.CI_DTdeltaC{ii}(3);    
end
c.pedC = 100.*c.pedC;

minC = 1;
c.pedC(1) = minC;
f = pwd;
save_tag = f(end-26:end-22);

%% plot mean contrast presented at each pedc and sd:
fv(7) = figure;
thetitle =['Behavior'];
set(fv(7),'NumberTitle','off','Name',[save_tag,'_',thetitle]);

for ii = 1:length(c.CI_ATcontrast)
    a = 100.*c.aTcontrasts{ii};
    a = a((a>adc_min(ii) & a<adc_max(ii)));
    ax = c.pedC(ii)+10.^randn(size(a))/40;pwd
    
    d = 100.*c.dTcontrasts{ii};
    d = d((d>ddc_min(ii) & d<ddc_max(ii)));
    dx = c.pedC(ii)+10.^randn(size(d))/40;
    
    loglog(ax,a,'rs','MarkerSize',1,'Color',myRed,'LineWidth',lineThickness);%, ...
    hold on
    loglog(dx,d,'ks','MarkerSize',1,'Color',myBlack,'LineWidth',lineThickness);
end
loglog(c.pedC,adc_fiftyQ,'ro-','MarkerSize',10,'Color',myRed,'LineWidth',lineThickness);%, ...
loglog(c.pedC,ddc_fiftyQ,'ks-','MarkerSize',10,'Color',myBlack,'LineWidth',lineThickness)

xlabel('Pedestal contrast (%)','FontSize',Fsize);
ylabel('Delta contrast (%)','FontSize',Fsize);
title('Discrimination performance','FontSize',Fsize)
set(gca,...
    'FontName','Helvetica','FontSize',Fsize,'PlotBoxAspectRatio',[1 1 1], 'XLim', [.9 100],...
    'LineWidth',1,'TickLength',[0.025 .01], 'XTick', [.9 c.pedC 100], ...
    'XTickLabel', [0 0 c.pedC(2:end) 0],'YTickLabel', [.1 1 10 100],...
    'XColor',XYColor,'YColor',XYColor,'Box','off');
drawnow