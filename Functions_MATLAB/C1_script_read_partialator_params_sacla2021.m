% script_read_partialator_params
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script read the partialator parameters files (*.params files) and 
% extract the following information:
%  nRun:  run numbers
%  nEvent: event numbers
%  OSF: Waller-Debye factor-1
%  relB: Waller-Debye factor-2 
% Ahmad H., Feb-2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear 

format long g
%NAME = 'upto1ps';
NAME = '3ps';

%NAME = 'dark';
%NAME = 'puredark-clean';  % Newest dark file from Uppsala

fid  = fopen([NAME,'.params'],'r');
text = textscan(fid,'%s %s %s %s %s %s %s');
fclose(fid);
whos text

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Extract Partialator parameters OSF and relB:
column2 = sprintf('%s ', text{2}{2:end}); 
column3 = sprintf('%s ', text{3}{2:end}); 

OSF  = sscanf(column2, '%f');
relB = sscanf(column3, '%f');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Extract run numbers:
Key = 'run';
row_idx = find(~cellfun('isempty',strfind(text{6},Key)));
nRun = cell(length(row_idx),1);
for k=1:length(row_idx)
 Str = text{6}{row_idx(k)};
 Index = strfind(Str, Key);
 if ~isempty(Index) 
  nRun{k} = sscanf(Str(Index(1) + length(Key):end), '%f', 1);
 end
end
nRun = cell2mat(nRun);


%##########################################################################
% % extract event numbers
Key = 'tag-';
row_idx = find(~cellfun('isempty',strfind(text{7},Key)));
nEvent = cell(length(row_idx),1);
for k=1:length(row_idx)
 Str = text{7}{row_idx(k)};
 Index = strfind(Str, Key);
 if ~isempty(Index) 
  nEvent{k} = sscanf(Str(Index(1) + length(Key):end), '%f', 1);
 end
end
nEvent = cell2mat(nEvent);


%##########################################################################
% % double-check before saving
%{
id=isnan(nRun); find(id==1)
id=isnan(nEvent); find(id==1)
id=isnan(relB); find(id==1)
id=isnan(OSF); find(id==1)
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% save the extracted information
save(['partialator_params_phyto_',NAME,'.mat'],'nRun','nEvent','OSF','relB')


% EOF
