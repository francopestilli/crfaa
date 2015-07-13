function make_crf_structures(varargin)
% function make_crf_structures(varargin)
%
% this function loads up analyses files created by 
% runAllanalyses_r_R1.m and does three things:
%
% (1) for each roi r2 level, generates an average file similar to 
%     ADAPTATION-EFFECT-v<1-4>-19-Jan-2009.mat, used to run the 
%     crf fit and noise/decision models. so that the analyses can be 
%     repeated at each r2 level. 
%
% 
% (2) saves all the new files into the folder: crfaa/data_used_files_R1
%
% franco pestilli 2010/05/18

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
observer  =[];  % the observer
visualarea=[];  % the visual area
adaptation=[];  % the adaptation condition
subtractionType=[]; % type of subtraction done on the time series ORIG or BYTRIAL
r2=[];              % the r2 values that define the file name/roi used
projectPath=[];     % path to the folder where all the data for crfaa are
makeAverageFile=[];        % choose whether to recreate the average file
combineContrastResponse=[];% choose whether to re-create a file that combines crf with different r2 values 
deleteOriginalFiles=[];    % choose whether to delete the files being imported in data_used_files_R1 from their original location
dataFolder=[];             % data_used_files_R1
originalDataFolder=[];     % data_used_files
getArgs(varargin,{...
 'makeAverageFile',1,                 ...
 'combineContrastResponse',1,         ...
 'deleteOriginalFiles',0,             ...
 'subtractionType', {'ORIG' 'BYTRIAL'},                   ...
 'dataFolder','data_used_files_R1',   ...
 'originalDataFolder','data_used_files', ...
 'observer',{'fp' 'fm' 'jg'},         ...
 'visualarea',{'v1' 'v2' 'v3' 'v4'},  ...
 'adaptation',{0},                    ...
 'r2',{'030' '050' '070'},...
 'projectPath','/data2/crfaa/crfaa/'});
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% (0) make file names for files to load for each observer, area, r2  and
% adaptation
filename = makeFiles2load(observer,visualarea,adaptation,subtractionType,r2,dataFolder);

% (1) make average files with data for each adaptation condition and each roi/r2 
makeAverage(filename,originalDataFolder,visualarea,adaptation,projectPath,observer, subtractionType,dataFolder,r2)




%%%%%%%%%%%%%%%
% makeAverage %
%%%%%%%%%%%%%%%
function makeAverage(filename,originalDataFolder, visualarea,adaptation, ...
                     projectPath,observer, subtractionType, dataFolder,r2)

d     = cell(1,3);
amplitudes = nan.*ones(3,32);
% for v = 1:length(visualarea)
%  file2load = fullfile(projectPath,originalDataFolder,sprintf('ADAPTATION-EFFECT-%s-19-Jan-2009.mat',visualarea{v}));
%  disp(sprintf('[%s] loading: %s',mfilename,file2load))
%  old = load(file2load,'d','meanC');
 


 for s = 1:length(subtractionType)
  for v = 1:length(visualarea)
   for r = 1:length(r2)
    for a = 1:length(adaptation)
     for o = 1:length(observer)
      if o == 1 && r == 1 && s == 1
       %load the old file to extract the mean contrast computed from behavior
       % this is not changed and will need to be saved back into the new file.
       old_file2load = fullfile(projectPath,originalDataFolder, ...
        sprintf('ADAPTATION-EFFECT-%s-19-Jan-2009.mat',visualarea{v}));
       disp(sprintf('[%s] loading old file to extract meanC: %s',mfilename,old_file2load))
       old = load(old_file2load,'d','meanC');
       
       % save meanC
       meanC = old.meanC;
       clear old old_file2load;
      end
      
     % this is the new file that should be loaded
     file2load = filename{o,s,v,r,a};
     if ~isempty(file2load)
      if o == 1
       disp(sprintf('[%s] extracting amplitudes from:\n%s',mfilename,file2load));
      else
       disp(sprintf('%s',file2load));
      end
      new = load(file2load);
      amplitudes(o,:) = new.d1.amplitude;
     else
      keyboard
     end
    end
    d{a}.Dt_amplitude  = nanmean(amplitudes(:,17:24));
    d{a}.Dnt_amplitude = nanmean(amplitudes(:,25:32));
    d{a}.A_amplitude   = nanmean(amplitudes(:,1:8));
    d{a}.U_amplitude   = nanmean(amplitudes(:,9:16));
    d{a}.Dt_ste        = nanstd(amplitudes(:,17:24))/3;
    d{a}.Dnt_ste       = nanstd(amplitudes(:,25:32))/3;
    d{a}.A_ste         = nanstd(amplitudes(:,1:8))/3;
    d{a}.U_ste         = nanstd(amplitudes(:,9:16))/3;
   end
   
   % save the file
   if ~isempty(file2load)
    file2save = fullfile(projectPath,dataFolder, ...
     sprintf('ADAPTATION-EFFECT-%s-%s-%s.mat', visualarea{v},r2{r},subtractionType{s}));
    disp(sprintf('[%s] saving: %s',mfilename,file2save))
    save(file2save,'d','meanC');
   end
  end
 end
end


%%%%%%%%%%%%%%%%%%
% makeFiles2load %
%%%%%%%%%%%%%%%%%%
function filename = makeFiles2load(observer,visualarea,adaptation,subtractionType,r2,defaultDataFolder)
temp_filename = '';filename = [];

for o = 1:length(observer)
 for s = 1:length(subtractionType)
  for v = 1:length(visualarea)
   for r = 1:length(r2)
    for a = 1:length(adaptation)
     temp_filename = makeFileName(observer{o},adaptation{a},visualarea{v},defaultDataFolder,r2{r}, subtractionType{s});
     if (isfile(temp_filename))
      fn = temp_filename;
      filename{o,s,v,r,a} = fn;
     else
      keyboard
     end
    end
   end
  end
 end
end




%%%%%%%%%%%%%%%%
% makeFileName %
%%%%%%%%%%%%%%%%
function [filename vars2load] = makeFileName(observer,adaptation,visualArea,defaultDataFolder,r2, subtractionType)
% NB whereas the main analysis files are saved into the special folder
% <data_used_files> so that changes to that folder do not screw up the
% files...
% the files for the other analyses are instad in each folder/session. so
% for construtcting the filename we use getsession.

if strcmpi(observer,'avrg')
 if strcmp(r2,'none')
  filename =  fullfile(defaultDataFolder,sprintf('ADAPTATION-EFFECT-%s-19-Jan-2009.mat', visualArea));
 else
  filename =  fullfile(defaultDataFolder,sprintf('ADAPTATION-EFFECT-%s-%s-%s.mat', visualArea,r2,subtractionType));
 end
 
else
 
 if strcmp(r2,'070') && strcmp(visualArea,'v1')
  filename = fullfile('/data2/crfaa/crfaa/data_used_files/', ...
   sprintf('%s_A%i_29-Dec-2008_%s_0_7roi.mat',upper(observer), adaptation,visualArea));
 elseif strcmp(r2,'050') && ~strcmp(visualArea,'v1')
  filename = fullfile('/data2/crfaa/crfaa/data_used_files/', ...
   sprintf('%s_A%i_29-Dec-2008_%s_0_5roi.mat',upper(observer), adaptation,visualArea));
 else
  filename = fullfile('/data2/crfaa/crfaa/data_used_files_R1/', ...
   sprintf('%s_A%i_2010-30-May_%s_era_R1_subtr%s_%sr%s.mat',upper(observer), adaptation,visualArea,subtractionType, visualArea,r2));
 end
end

% %%%%%%%%%%%%%%%
% % getFileName %
% %%%%%%%%%%%%%%%
% function file = getFileName(visualarea,subtractionType,filedate,r2)
% keyboard
% if strcmp(r2,'050') && ~strcmp(visualarea,'v1')
%  file = sprintf('%s_%s_era_R1_subtr%s_%s.mat',filedate,visualarea,subtractionType,visualarea);
% elseif strcmp(r2,'070') && strcmp(visualarea,'v1')
%  file = sprintf('%s_%s_era_R1_subtr%s_%s.mat',filedate,visualarea,subtractionType,visualarea);
% elseif  strcmp(r2,'090') && ~strcmp(visualarea,'v1')
%  file = [];
% else
%  file = sprintf('%s_%s_era_R1_subtr%s_%sr%s.mat',filedate,visualarea,subtractionType,visualarea,r2);
% end
%  
%  
% %%%%%%%%%%%%%%%%%
% % getFoldername %
% %%%%%%%%%%%%%%%%%
% function observerFolder = getFolderName(observer,adaptation)
% observer = upper(observer);
% 
% % choose the corrrect mrTools session name 
% obsadapt = [observer,num2str(adaptation)];
% switch obsadapt
%  case {'FP0'}
%   session = 'fp20071019';
%  case {'FP28'}
%   session = 'fp20080402';
%  case {'FP100'}
%   session = 'fp20080415';
%  case {'JG0'}
%   session = 'jg20070919';
%  case {'JG28'}
%   session = 'jg20080402';
%  case {'JG100'}
%   session = 'jg20080414';
%  case {'FM0'}
%   session = 'fm20080209';
%  case {'FM28'}
%   session = 'fm20080406';
%  case {'FM100'}
%   session = '';
%  otherwise
%   keyboard
% end
% 
% observerFolder = fullfile(sprintf('%s_A%i_testStream',observer,adaptation),session);
