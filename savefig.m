% function figDir = savefig(h,varargin)
% saves a figure for publishing purpose.
%
% <options>:
% figName            % the name of the figure file (default is YYYYMMDDThhmmss_figure)
% defaultDataFolder  % where all the data and results for this project are (default is current dir)
% figDir             % name of the subfolder where to save this figure (default is figures)
% dopng              % to make a dokuwiki compatible file
% verbose            % display on screen what it is going on
% 
% it calls this command:
% print(handle, '-cmyk', '-painters','-depsc2','-tiff','-r500', '-noui', 'figureName')
%
% example call:
% savefig(figureHandle,'verbose',1,'figName',figure1,'figDir','thisfigures','defaultDataFolder','~/allfigures');
%
% Franco Pestilli 2009.03.04

function figDir = savefig(h,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set up default variables:
figName           = []; % the name of the figure file
defaultDataFolder = []; % where all the data and results for this project are (default is current dir)
figDir            = []; % the subfolder where to save this figure
dopng             = []; % to make an easy to read and dokuwiki compatible file
verbose           = []; % display on screen what it is going on

getArgs(varargin,{ ...
        'figName',sprintf('%s_figure',datestr(now,30)), ...
        'defaultDataFolder','.', ...
        'figDir','figures', ...
        'dopng',0, ...
        'verbose',0});
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% make the figure dir if it does not exist:
figDir = sprintf('%s/%s',defaultDataFolder,figDir);
currentDir = pwd;
if isdir(figDir)
 cd(figDir);
else
 mkdir(figDir);
 cd(figDir);
end

% create a print command that will save the figure
printCommand = sprintf('print(%s, ''-cmyk'', ''-painters'',''-depsc2'',''-tiff'',''-r500'' , ''-noui'', ''%s'')', num2str(h),figName);
if verbose
 disp(sprintf('[%s] saving eps figure %s/%s\nUsing command: %s',mfilename,figDir,figName,printCommand));
end

% do the printing here:
eval(printCommand);

% save a dokuwiki-compatible file
if dopng
 disp(sprintf('[%s] saving png figure %s',mfilename,figName));
 eval(sprintf('print(%s, ''-painters'',''-dpng'', ''-noui'', ''%s'')',  ...
  num2str(h),figName));
end
cd(currentDir);
