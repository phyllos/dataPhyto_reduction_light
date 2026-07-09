% script_pack_phyto_light_lpf_drl
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script reads PYP data (i.e., h, k, l, T, M, qmax), calculates vector "q(h,k,l)" for 
% Bragg spots in each snapshot, then remove the "q(h,k,l)" below qmax for each snapshot.
% Then finds the average and std and save them into an HKL file for next steps
% Ahmad H., Dec-22-2017, March 2018, May 2018
% update for phytochromes on Feb-2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear

tic
dataName = 'upto1ps'; %'3ps'
dataForm = 'int'; %'int' vs. 'amp';  % amplitude vs. intensity

% %if dataName == '3ps'
%  nBrg = 61784;  % number of active Bragg spots
%  nS   = 35366;  % number of snapshots
%if dataName == 'upto1ps'
 nBrg = 62530;  % number of active Bragg spots
 nS   = 68084;%24726;%8206; %(uniformed) % 135041 (original); %number of snapshots
 
remove_negative_pixels = false;
    
% lattice sizes for group P212121:
a = 54.95; %54.22;  % in Angstrom
b = 116.5; %115.78; % in Angstrom
c = 117.7; %117.08; % in Angstrom

% if dataName == '3ps'
% fileData = 'dataPhyto_3ps_int_sortdelay_LPF_nS35366_nBrg61784.mat'; % to load
% %fileMask_DRL = sprintf('maskDRL_%s_%s_nS%d_nBrg%d.mat',dataName, dataForm,nS,nBrg); % to load
% fileData_DRL = 'dataPhyto_3ps_int_sortdelay_LPF_DRL_nS35366_nBrg61784.mat';  % to save
% if dataName == '1ps'
fileData = sprintf('dataPhyto_upto1ps_int_sortdelay_unifdelay_LPF_nS%d_nBrg%d.mat',nS,nBrg); % to load
% %fileMask_DRL = sprintf('maskDRL_%s_%s_nS%d_nBrg%d.mat',dataName, dataForm,nS,nBrg); % to load
fileData_DRL = sprintf('dataPhyto_upto1ps_int_sortdelay_unifdelay_LPF_DRL_nS%d_nBrg%d.mat',nS,nBrg); % to save

fileDir = ['./DATA_',dataName];  % the directory to save the output

%-----------------------------------------------------------------------------------------
FIL1 = fullfile(fileDir, fileData);
load(FIL1,'T_lpf','M','DRL','LPF','OSF','relB','miller_h','miller_k',...
   'runID','eventID','delay','miller_l','sort_notice','notice_negative_pix');

% Applying DRL mask:
%FIL2 = fullfile(fileDir, fileMask_DRL);
%load(FIL2,'maskDRL');
maskDRL = ones(nS,nBrg);    % i.e., aplying NO mask. test by AHZ, April-2022
M_drl = maskDRL.*M;
T_lpf_drl = M_drl.*T_lpf;

%noticeDRL = 'Bragg reflections beyond DRLs are removed, and qmax is q at DRL.';
noticeDRL = 'Bragg reflections beyond DRL NOT removed: maskDRL = ones(nS,nBrg)';
FIL3 = fullfile(fileDir, fileData_DRL);
save(FIL3, 'T_lpf_drl', 'M_drl','miller_h','miller_k','miller_l',...
           'delay','noticeDRL','notice_negative_pix','DRL','LPF',...
           'runID','eventID','sort_notice','OSF','relB','-v7.3');  
toc
%EOF
