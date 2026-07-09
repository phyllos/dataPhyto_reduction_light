% script_read_stream_files_sacla2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script read the stream file "merged_detwinned.stream" and extract 
% the following information:
%  nRun:  run numbers
%  nEvent: event numbers
%  nRefl: num of reflections
%  DRL:   diffraction resolution limits
% Notice-1: The above information for the snapshots that have not been 
% indexed must be removed, and this is done at the end of this code.
% Notice-2: The extracted DRL are in 1/A (i.e., it is q NOT d=1/q) here. 
% In RS data files DRL numbers are in [A] (i.e., d=1/q are extracted).
% Ahmad H., Jan-12-2018
%--------------------------------------------------------------------------
% updated on Feb-2022 for SACLA data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear 

tic
%##########################################################################
% % read the strem file
%NAME = 'upto1ps';
NAME = '3ps';

%NAME = 'dark';
%NAME = 'dark-clean';
%NAME = 'puredark-clean';  % Newest dark file from Uppsala

fid  = fopen([NAME,'.stream'],'r');
text = textscan(fid,'%s','Delimiter','');
text = text{1};
fclose(fid);


%##########################################################################
% extract num of reflections
nRefl = regexp(text,'num_reflections[\s\.=]+(\d+)','tokens');
nRefl = [nRefl{:}]; nRefl = str2double([nRefl{:}]).';


%##########################################################################
% % extract event numbers
nEvent = regexp(text,'Event[:\s\.tag-]+(\d+)','tokens');
nEvent = [nEvent{:}]; nEvent = str2double([nEvent{:}]).';

%
%##########################################################################
% % extract diffraction resolution limits
DRL = regexp(text,'diffraction_resolution_limit[\s\.=]+(\d+\.+\d*)','tokens');
DRL = [DRL{:}]; DRL = str2double([DRL{:}]).';


%#########################################################################################
% % extract run numbers
Key = 'run';
row_idx = find(~cellfun('isempty',strfind(text,Key)));
nRun = cell(length(row_idx),1);
for k=1:length(row_idx)
 Str = text{row_idx(k)};
 Index = strfind(Str, Key);
 if ~isempty(Index) 
  nRun{k} = sscanf(Str(Index(1) + length(Key):end), '%f', 1);
 end
end
nRun = cell2mat(nRun);


%##########################################################################
% grab snapshots with "indexed_by=none" 
KEY = 'indexed_by';
row_index = find(~cellfun('isempty',strfind(text,KEY)));
VAL = cell(length(row_index),1);
for ii=1:length(row_index)
 Str = text{row_index(ii)};
 Ind = strfind(Str, KEY);
 if ~isempty(Ind) 
  Str(strfind(Str, '=')) = [];
  VAL{ii} = sscanf(Str(Ind(1) + length(KEY):end), '%s', 1);
 end
end
idx_none = strcmp(VAL,'none');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% remove snapshots (run & event numbers) assigned to "indexed_by=none"
nRun(idx_none) = [];
nEvent(idx_none) = [];

if ~isequal(length(nRun), length(nEvent), length(nRefl), length(DRL))
   error('the lengths of vectors are different!')
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% save the extracted information
save(['merged_snapshotInfo_phyto_',NAME,'.mat'],'nRun','nEvent','DRL','nRefl')

% % save the same thing as a *.dat file
% M = [nRun nEvent nRefl DRL];
% dlmwrite(['merged_snapshotInfo_phytochrome_',NAME,'.dat'],M,'delimiter','\t','precision',16)

myData = [nRun, nEvent, nRefl, DRL]';
fileID = fopen(['merged_snapshotInfo_phyto_',NAME,'.dat'],'w');
%fprintf(fileID,'%12s %16s %8s %8s\n','nRun','nEvent','nRefl','DRL');
fprintf(fileID,'%12d %16d %8d %8.2f\n',myData);
fclose(fileID);


toc
% EOF
