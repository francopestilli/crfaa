function makeAverageTvC
% function makeAverageTvC
%
% this function makes an average TvC using the 'results'
% structure returned by makenoisestatistics2_r.m
%
% it is called by itself and it saves a file in the current dir
% the saved file is used by dprimefit3_rAVRG.m to compute noise estimates
% on the average across observers.
%
% franco pestilli 2009/07/12


% load the file with the individual results (this is now hard coded):
r = load('/Volumes/data/riken/crfaa/fmridata/crfaa/10-Jul-2009AllNoise.mat','results');
r = squeeze(r.results);

% filed meaning:
% {iObserver, iAdapt, iArea, icrf}
% icrf=1 -> distributed
% icrf=2 -> attended

h = smartfig('makeAverageTvC','reuse');

% make average TvC:
meanTvC   = zeros(size(r,2),size(r,4),8);
steTvC    = meanTvC;
model_tvc = meanTvC;
for iAdapt = 1:size(r,2)
 tvc = nan.*ones(size(r,1),8);
 for icrf = 1:size(r,4)
  for iObserver = 1:size(r,1)
   if ~(iAdapt==3 && iObserver ==3)
    tvc(iObserver,:) = r{iObserver,iAdapt,1,icrf}.behavior.tvc;
   else
    disp('observer ''3'' (fm), does not have the last adaptation condition')
   end
  end
  meanTvC(iAdapt,icrf,:) = 10.^nanmean(log10(tvc));
  steTvC(iAdapt,icrf,:) = nanstd(tvc)./sqrt(3);
  
  % plot the results:
  symbol= 'o';
  pedestals = [.00875 .0175 .035 .07 .14 .28 .56 .84];
  plotTvC(h,pedestals,squeeze(meanTvC(iAdapt,icrf,:)),squeeze(steTvC(iAdapt,icrf,:))',iAdapt,icrf,symbol)
  
  % fit a tvc
  %            slope  thr        d         T    beta
  initParams = [1.4    0.1       0.065    10    2];
  minParams  = [1      0.0001    0         5     0];
  maxParams  = [3.5    0.25      1        inf   inf];
  
  % fit the function
  pedestals = [.00875 .0175 .035 .07 .14 .28 .56 .84];
  t = squeeze(meanTvC(iAdapt,icrf,:))';
  [fitParams(iAdapt,icrf,:) model_tvc] = fitnakarushtonTvC(pedestals(1:end),t(1:end),initParams,minParams,maxParams);
  
  % plot the model:
  symbol = '-';
  plotTvC(h,model_tvc.c,model_tvc.t,ones(size(model_tvc.c)),iAdapt,icrf,symbol)
  
 end
end


% save results:
filename = '/Volumes/data/riken/crfaa/fmridata/crfaa/12-Jul-2009averageTvC.mat';
functionName = mfilename;
save(sprintf('%s',filename),'meanTvC', 'steTvC','model_tvc','functionName');
disp(sprintf('Saving file: <%s>',filename))

keyboard



% the code that follwos fits a function to this average data.
%%%%%%%%%%%%%%%%%%%%%%%
%% fitnakarushtonTvC %%
%%%%%%%%%%%%%%%%%%%%%%%
function [fitParams model_tvc] = fitnakarushtonTvC(c,t,initParams,minParams,maxParams)

% optimization parameters
maxIter = inf;
MaxFunEvals = 10^4;
optimParams = optimset('LevenbergMarquardt','on', ...
 'MaxIter',maxIter,         ...
 'MaxFunEvals',MaxFunEvals, ...
 'Display','off',           ...
 'Diagnostics','off',       ...
 'TolFun',10^-4);


% use lsqnonlin to find function minimum
[fitParams resnorm residual exitflag output lambda jacobian] = lsqnonlin(@resNakarushton,initParams,minParams,maxParams,optimParams,c,t);

% return the fitted function
model_tvc.c = c(1):.0001:c(end);
model_tvc.t = nakarushtonTvC(model_tvc.c,fitParams)';


%%%%%%%%%%%%%%%%%%%%%%%%%
%%    resNakarushton    %%
%%%%%%%%%%%%%%%%%%%%%%%%%
function residual = resNakarushton(params,pedestal,data)

model = nakarushtonTvC(pedestal,params);

residual = log(data)-log(model);


dispFit = 1;
if dispFit
 smartfig('dprimefit_nakaTvC','reuse');
 cla;
 loglog(pedestal,data,'bs');hold on
 loglog(pedestal,model,'r.-');
 axis('square')
 
 
 rmax  = '1';%num2str(params(1)); % it is now fixed to 1
 n     = num2str(params(1));
 c50   = num2str(params(2));
 d     = num2str(params(3));
 T   = num2str(params(2));
 beta     = num2str(params(3));
 
 txt = sprintf('rmax: %s, n: %s, c50: %s,\nd: %s\nT: %s beta: %s.',rmax,n,c50,d, T, beta); 
 title(txt)
end


%%%%%%%%%%%%%%%%%%%%%%%%%
%%    nakarushtonTvC   %%
%%%%%%%%%%%%%%%%%%%%%%%%%
function deltac = nakarushtonTvC(pedestals,params)
rmax  = 1;%params(1);
n     = params(1);
c50   = params(2);
d     = params(3); % differentiation factor
T     = params(4).*c50;
beta  = params(5);

c     = pedestals;
y = dip2(c,T,beta);

% breaking down the function into
% components for simplicity (kinda)
v1 = c.^n;
v2 = c50^n;
v3 = v1+v2;
v4 = v1./v3;
v5 = d/rmax;
v6 = (1./(v5+v4))-1;
v7 = v2./v6;
v8 = v7.^(1/n);
v9 = v8 - c;
deltac = v9.*y;


%%%%%%%%%%%
%% adapt %%
%%%%%%%%%%%
function y = dip2(c,T,beta)
% this one fixes the decrease in t at high contrast
% parameters of the copressive function at high contrast
% T center (threshold) of the function
% beta slope of the function
% y = .015+exp(-(c./T).^beta-.015);
y = exp(-(c./T).^beta);

%%%%%%%%%%%%%%
%% plotTvC  %%
%%%%%%%%%%%%%%
function plotTvC(h,pedestals,meanTvC,steTvC,iAdapt,icrf,symbol)

figure(h)
if icrf == 1
 color = 'b';
else
 color = 'r';
end

subplot(1,3,iAdapt),
hold on
myerrorbar(pedestals,meanTvC, ...
 'yError',steTvC, ...
 'Symbol',symbol,'MarkerFaceColor',color, ...
 'MarkerEdgeColor',[0 0 0],'MarkerSize',8, ...
 'Color',color);

ticks = [.00875 .0175 .035 .07 .14 .28 .56 .84];

set(gca,...
 'FontName','Helvetica','FontSize',12, ...
 'PlotBoxAspectRatio',[1 1 1], ...
 'XLim', [.0035 1],'YLim', [.0035 .45],...
 'LineWidth',1, ...
 'yScale','log','xScale','log','XTick', [], 'XTickLabel', 100*[0 ticks(2:end)] ,...
 'YTick', ticks, 'YTickLabel', 100*[0 ticks(2:end)], ...
 'XColor',[0 0 0],'YColor',[0 0 0],'Box','off');
