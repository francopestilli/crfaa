function c = plotCRFtimecourse(d,events,save_tag)
% c = plotCRFtimecourse(d,events,save_tag)
%
% - d is a struct as returned by fiTimecourse:
% d = 
% 
%         amplitude: [1x32 double]
%     amplitudeMEAN: [1x32 double]
%      amplitudeSTD: [1x32 double]
%      amplitudeSTE: [1x32 double]
%     amplitudeType: 'max'
%             covar: [32x32 double]
%            deconv: [1x1 struct]
%              ehdr: [32x26 double]
%           fitType: 'glm'
%                r2: 0.0851
%              time: [1x26 double]
% 
%
%
%
% franco pestilli 2008.6.02


if nargin == 2
    save_tag = [date,'_'];

elseif nargin < 2
    help plotCRFtimecourse;
    return

end

%% compute average contrast presented for attended and distributed-target conditions:

disp('[plotCRFtimecourse] not computing contrast.');
c = computedMeanC(events);

% attentional conditions:
A_conditions = 1:8;
U_conditions = 8+1:2*8;
Dt_conditions = 2*8+1:3*8;
Dnt_conditions = 3*8+1:4*8;

% plot basic set up:
XYColor = [.4 .4 .4];
Fsize = 16;
lineThickness = 3;
myRed = [.6 0 0];
myBlack = [0 0 0];
myGreen = [0 .4 0];
Ymax = 1.5;

% fix contrast values
c.quantileAcontrast = 100.*c.quantileAcontrast;
c.quantileDcontrast = 100.*c.quantileDcontrast;
c.quantileAdeltaCs = 100.*c.quantileAdeltaCs;
c.quantileDdeltaCs = 100.*c.quantileDdeltaCs;
c.pedC = 100.*c.pedC;
c.pedC(1) = .1;

fv(1) = figure;
thetitle =['amplitude_CRF'];
set(fv(1),'NumberTitle','off','Name',[save_tag,'_',thetitle]);

semilogx(c.quantileAcontrast(:,2),d.amplitude(A_conditions),'r^-','MarkerSize',8,'Color',myRed,'LineWidth',lineThickness);
hold on;
semilogx(c.quantileDcontrast(:,2),d.amplitude(Dt_conditions),'ro-','MarkerSize',8,'Color',myBlack,'LineWidth',lineThickness);
semilogx(c.pedC,d.amplitude(Dnt_conditions),'ro--','MarkerSize',8,'Color',myBlack,'LineWidth',lineThickness);
semilogx(c.pedC,d.amplitude(U_conditions),'rv-','MarkerSize',8,'Color',myGreen,'LineWidth',lineThickness);

h = legend({'Attended','Dist_target','Dist_noTarget','Unattended'},2); set(h, 'interpreter','none', 'box','off');

% plot error  bars:
semilogx([c.quantileAcontrast(:,2), c.quantileAcontrast(:,2)]', ... 
    [d.amplitude(A_conditions)-d.amplitudeSTE(A_conditions); d.amplitude(A_conditions)+d.amplitudeSTE(A_conditions)],'r-','Color',myRed,'LineWidth',.5);
semilogx([c.quantileDcontrast(:,2),c.quantileDcontrast(:,2)]', ...
    [d.amplitude(Dt_conditions)-d.amplitudeSTE(Dt_conditions); d.amplitude(Dt_conditions)+d.amplitudeSTE(Dt_conditions)],'r-','Color',myBlack,'LineWidth',.5);
semilogx([c.pedC;c.pedC], ...
    [d.amplitude(Dnt_conditions)-d.amplitudeSTE(Dnt_conditions); d.amplitude(Dnt_conditions)+d.amplitudeSTE(Dnt_conditions)],'r-','Color',myBlack,'LineWidth',.5);
semilogx([c.pedC;c.pedC], ...
    [d.amplitude(U_conditions)-d.amplitudeSTE(U_conditions); d.amplitude(U_conditions)+d.amplitudeSTE(U_conditions)],'r-','Color',myGreen,'LineWidth',.5);

xlabel('% Contrast','FontSize',Fsize);
ylabel('fMRI response (% Signal change)','FontSize',Fsize);
title('Contrast response','FontSize',16)
set(gca,...
    'FontName','Helvetica','FontSize',Fsize,'PlotBoxAspectRatio',[1 1 1], 'XLim', [min(c.pedC) 100],...
    'YLim', [-0.05 Ymax],...
    'LineWidth',1,'TickLength',[0.025 .01], 'XTick', [c.pedC], ...
    'XTickLabel', [0 c.pedC(2:end)],...
    'XColor',XYColor,'YColor',XYColor,'Box','off');

disp(sprintf('Saving figure %s',[save_tag,'_',thetitle]));
eval(sprintf('print(gcf,''-depsc2'',''-tiff'',''-r300'', ''%s'')', [save_tag,'_',thetitle]));

% % plot crf obtained by nlinfit
% fv(2) = figure;
% thetitle =['nlinfit_CRF'];
% set(fv(1),'NumberTitle','off','Name',[save_tag,'_',thetitle]);
% 
% semilogx(c.quantileAcontrast(:,2),d.amplitude(A_conditions),'r^-','MarkerSize',8,'Color',myRed,'LineWidth',lineThickness);
% hold on;
% semilogx(c.quantileDcontrast(:,2),d.amplitude(Dt_conditions),'ro-','MarkerSize',8,'Color',myBlack,'LineWidth',lineThickness);
% semilogx(c.quantileAcontrast(:,2),d.amplitude(Dnt_conditions),'ro--','MarkerSize',8,'Color',myBlack,'LineWidth',lineThickness);
% semilogx(c.quantileAcontrast(:,2),d.amplitude(U_conditions),'rv-','MarkerSize',8,'Color',myGreen,'LineWidth',lineThickness);
% 
% h = legend({'Attended','Dist_target','Dist_noTarget','Unattended'},2); set(h, 'interpreter','none');
% 
% % plot error  bars:
% semilogx([c.quantileAcontrast(:,2), c.quantileAcontrast(:,2)]', ... 
%     [d.amplitude(A_conditions)-d.amplitudeSTD(A_conditions); d.amplitude(A_conditions)+d.amplitudeSTD(A_conditions)],'r-','Color',myRed,'LineWidth',.5);
% semilogx([c.quantileDcontrast(:,2),c.quantileDcontrast(:,2)]', ...
%     [d.amplitude(Dt_conditions)-d.amplitudeSTD(Dt_conditions); d.amplitude(Dt_conditions)+d.amplitudeSTD(Dt_conditions)],'r-','Color',myBlack,'LineWidth',.5);
% semilogx([c.quantileAcontrast(:,2),c.quantileAcontrast(:,2)]', ...
%     [d.amplitude(Dnt_conditions)-d.amplitudeSTD(Dnt_conditions); d.amplitude(Dnt_conditions)+d.amplitudeSTD(Dnt_conditions)],'r-','Color',myBlack,'LineWidth',.5);
% semilogx([c.quantileAcontrast(:,2),c.quantileAcontrast(:,2)]', ...
%     [d.amplitude(U_conditions)-d.amplitudeSTD(U_conditions); d.amplitude(U_conditions)+d.amplitudeSTD(U_conditions)],'r-','Color',myGreen,'LineWidth',.5);
% 
% 
% xlabel('% Contrast','FontSize',Fsize);
% ylabel('% Signal change','FontSize',Fsize);
% title('Contrast response','FontSize',16)
% set(gca,...
%     'FontName','Helvetica','FontSize',Fsize,'PlotBoxAspectRatio',[1 1 1], 'XLim', [min(c.pedC) 100],...
%     'LineWidth',1,'TickLength',[0.025 .01], 'XTick', [c.pedC], ...
%     'XTickLabel', [0 c.pedC(2:end)],...
%     'XColor',XYColor,'YColor',XYColor,'Box','off');
% 
% disp(sprintf('Saving figure %s',[save_tag,'_',thetitle]));
% eval(sprintf('print(gcf,''-depsc2'',''-tiff'',''-r300'', ''%s'')', [save_tag,'_',thetitle]));
% 
% % plot mean aplitude vs std
% fv(3) = figure;
% thetitle =['mean_VS_variance'];
% set(fv(1),'NumberTitle','off','Name',[save_tag,'_',thetitle]);
% 
% plot(d.amplitude(A_conditions),d.amplitudeSTD(A_conditions).^2,'r^','MarkerSize',8,'Color',myRed,'LineWidth',lineThickness)
% hold on
% plot(d.amplitude(Dt_conditions),d.amplitudeSTD(Dt_conditions).^2,'ko','MarkerSize',8,'Color',myBlack,'LineWidth',lineThickness)
% plot(d.amplitude(Dnt_conditions),d.amplitudeSTD(Dnt_conditions).^2,'ko','MarkerSize',8,'Color',myBlack,'LineWidth',lineThickness)
% plot(d.amplitude(U_conditions),d.amplitudeSTD(U_conditions).^2,'gv','MarkerSize',8,'Color',myGreen,'LineWidth',lineThickness)
% 
% xlabel('Mean amplitude','FontSize',Fsize);
% ylabel('Variance','FontSize',Fsize);
% title('Correlation of mean and variance','FontSize',16)
% set(gca,...
%     'FontName','Helvetica','FontSize',Fsize,'PlotBoxAspectRatio',[1 1 1],'XLim', [0 max(d.amplitudeSTD.^2)], 'YLim', [0 max(d.amplitudeSTD.^2)], ...
%     'LineWidth',1,'TickLength',[0.025 .01], ...
%     'XColor',XYColor,'YColor',XYColor,'Box','off');
% 
% disp(sprintf('Saving figure %s',[save_tag,'_',thetitle]));
% eval(sprintf('print(gcf,''-depsc2'',''-tiff'',''-r300'', ''%s'')', [save_tag,'_',thetitle]));


