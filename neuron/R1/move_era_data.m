function move_era_data(varargin)
% function move_era_data(varargin)
%
% this function moves the era files in a session folder
% to the era files in data_used-files_R1
%
% new files are saved with observer initials and adaptation level
%
% these files are to be used to run the noise model  
%
%
% franco pestilli 2010/06/01

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
visualAreas     = [];
saveData        = [];
observer        = [];
adaptation      = [];
r2              = [];
subtractionType = [];
testonly        = [];
analysisType    = [];

getArgs(varargin,{...
 'visualAreas',{'v1','v2','v3','v4'}, ... 
 'saveData',1, ...
 'observer',{'fp' 'jg' 'fm'}, ...
 'adaptation',0, ...
 'r2', '090',    ...
 'subtractionType', 'BYTRIAL', ...
 'testonly',0, ...
 'analysisType', 'standard'});
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% folder where all the project is
defaultDataFolder = '/data2/crfaa/crfaa/';

% folder where the files will be moved to:
save_data_folder = 'data_used_files_R1';

for v = 1:length(visualAreas)
 for o = 1:length(observer)
  % getting the original file name
  [filename this_file] = makeFileName(observer{o},adaptation,visualAreas{v},defaultDataFolder, analysisType,r2, subtractionType);
  
  % check if the file exist
  if ~isfile(filename)
   keyboard
  end
  
  % creating the new file name with observer and adaptation info
  newfilename = sprintf('%s_A%s_%s',upper(observer{o}),num2str(adaptation),this_file);
  
  % saving the file with the new name
  file2save = fullfile(defaultDataFolder, save_data_folder,newfilename);
  if testonly
   eval(sprintf('!echo ''mv -v %s %s''',filename,file2save));
  else
   eval(sprintf('!mv -v %s %s',filename,file2save));
  end
 end
end



%%%%%%%%%%%%%%%%
% makeFileName %
%%%%%%%%%%%%%%%%
function [filename this_file session] = makeFileName(observer,adaptation,visualArea,defaultDataFolder, analysisType,r2, subtractionType)
                                                     
if strcmpi(observer,'avrg')
  %  this_file = sprintf('ADAPTATION-EFFECT-%s-%s-%s.mat', visualArea,r2,subtractionType);
  %  filename =  fullfile(defaultDataFolder,this_file);
   return
else
 % choose the type of analysis files to load
 [session this_date]  = getsession(observer,adaptation);
 switch analysisType
  case {'standard' '4conditions' '4c' '4' 's'} % this is the original analysis for the main results
   if (~strcmp(visualArea,'v1') && strcmp(r2,'050')) || (strcmp(visualArea,'v1') && strcmp(r2,'070'))
    this_file = sprintf('%s_%s_era_R1_subtr%s_%s.mat',this_date,visualArea, subtractionType,visualArea);
    filename = fullfile(session,this_file);
    
   else
    this_file = sprintf('%s_%s_era_R1_subtr%s_%sr%s.mat',this_date,visualArea,subtractionType,visualArea,r2);
    filename = fullfile(session,this_file);
   end

  otherwise
   keyboard
 end
end


%%%%%%%%%%%%%%%%
%  getsession  %
%%%%%%%%%%%%%%%%
function [session this_date] = getsession(observer,adaptation)

% check which computer are we using:
% check if there is a drive connected which one is it ...
if isdir('/data2/crfaa/crfaa')
 hd = '/data2/crfaa/crfaa';
else
 disp('No data hard drive found.')
 keyboard
end

switch lower(observer)
 case {'fm'}
  switch adaptation
   case {0}
    mrsession = 'fm20080209';
    this_date = '2010-30-May';
   case 28
    mrsession = 'fm20080406';
    this_date = '2010-30-May';
   otherwise
    disp(sprintf('This condition is NOT set (%s, adp-%i).',observer,adaptation))
    mrsession = [];
    this_date = '';
  end
  
 case {'fp'}
  switch adaptation
   case {0}
    mrsession = 'fp20071019';
    this_date = '2010-30-May';
   case {28}
    mrsession = 'fp20080402';
    this_date = '2010-30-May';
   case {100}
    mrsession = 'fp20080415';
    this_date = '2010-30-May';
   otherwise
    disp(sprintf('This condition is NOT set (%s, adp-%i).',observer,adaptation))
    keyboard
  end
  
 case {'jg'}
  switch adaptation
   case {0}
    mrsession = 'jg20070919';
    this_date = '2010-30-May';
   case {28}
    mrsession = 'jg20080402';
    this_date = '2010-30-May';
   case {100}
    mrsession = 'jg20080414';
    this_date = '2010-30-May';
   otherwise
    disp(sprintf('This condition is NOT set (%s, adp-%i).',observer,adaptation))
    keyboard
  end
  
 otherwise
  disp(sprintf('Observer (%s) NOT found.',observer))
  keyboard
end

if ~isempty(mrsession)
 session = fullfile(hd,sprintf('%s_A%i_testStream',upper(observer),adaptation),mrsession);
else
 session = '';
end

