%script_pack_phyto_allInfo_light
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script reads all required information from two DAT files and save 
% them together in a single DAT file.
% Ahmad H., Feb-2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear

% % read (runID,eventID,num_reflection,delay,qx,qy,qz) information
% % (qx,qy,qz) is the vector from (0,0,0) to center of spherical fit to the Ewald's sphere

file1 = 'merged_snapshotInfo_phyto_3ps.mat';
load(file1,'DRL','nEvent','nRun','nRefl');
nR1 = nRun;
nE1 = nEvent;

file2 = 'partialator_params_phyto_3ps.mat';
load(file2,'nRun','nEvent','OSF','relB')

% to save:
fileSave1 = 'merged_snapshotInfo_phyto_3ps_allInfo.mat';
fileSave2 = 'merged_snapshotInfo_phyto_3ps_allInfo.dat';

% double-check the two files have the same runID and eventID:
isequal(nR1,nRun)
isequal(nE1,nEvent)

file3 = 'delays_original.mat';
load(file3,'delay_3ps','tag_3ps')         % AHZ

% put info from two files together
info_ = [nRun, nEvent, nRefl, DRL, OSF, relB];

% finding intersection of info_ and delay file:
[C, IDX_info_, IDX_delay] = intersect(nEvent, tag_3ps);

delay   = delay_3ps(IDX_delay);

runID   = nRun(IDX_info_);
eventID = nEvent(IDX_info_);
OSF     = OSF(IDX_info_);
relB    = relB(IDX_info_);
nBrg    = nRefl(IDX_info_);
DRL     = DRL(IDX_info_);

% check the un-sorted delay times:
figure(1); subplot(121); plot(delay); ylabel('delay'); xlabel('delay #')
figure(1); subplot(122); histogram(delay, 600)

%%%%%%%%%%%%%%%%% SORTING WILL BE DONE IN THE NEXT SCRIPT %%%%%%%%%%%%%%%%%
% % sorting based on delay times
% [~, idSort] = sort(delay,'ascend');
% 
% delay = delay(idSort);
% runID = runID(idSort);
% eventID = eventID(idSort);
% OSF = OSF(idSort);
% relB = relB(idSort);
% nBrg = nBrg(idSort);
% DRL = DRL(idSort);
% 
% % check the sorted delay times:
% figure(2); subplot(121); plot(delay); ylabel('delay'); xlabel('delay #')
% figure(2); subplot(122); histogram(delay, 100)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nS = length(runID) 
% % put all info together
B = inf(nS,7);
B(:,1) = runID;         % run numbers
B(:,2) = eventID;       % event numbers
B(:,3) = nBrg;          % Bragg spots number
B(:,4) = DRL;           % diffraction resolution limit
B(:,5) = delay;         % delay time
B(:,6) = OSF;           % OSF param
B(:,7) = relB;          % relB param

myData = [runID, eventID, nBrg, DRL, delay, OSF, relB]';


notice='parameters NOT sorted based on delay yet.';
% % save the information
save(fileSave1,'runID','eventID',...
     'nBrg', 'DRL', 'delay', 'OSF', 'relB', 'notice', '-v7.3');
%
fileID = fopen(fileSave2 ,'w');
%fprintf(fileID,'%12s %16s %8s %8s\n','nRun','nEvent','nRefl','DRL');
fprintf(fileID,'%2d %14d %8d %8.2f %10.2f %8.2f %8.2f\n',myData);
fclose(fileID);


% EOF



