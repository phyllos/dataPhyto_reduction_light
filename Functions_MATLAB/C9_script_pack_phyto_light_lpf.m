% script_pack_phyto_light_lpf
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script reads data (h, k, l, I), calculates Lorentz-Polarisation
% factor for each snapshot and divide the corresponding Bragg intensity by 
% that and save them into an HKL file for next steps.
% Ahmad H., for phytochromes on April-2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

lambda = 1.77;    % wavelength in Angstrom (E=7000 eV)
remove_negative_pixels = false;
    
% %lattice sizes for group P212121:
a = 54.95; %54.22;  % in Angstrom
b = 116.5; %115.78; % in Angstrom
c = 117.7; %117.08; % in Angstrom

%if dataName == 'upto1ps'
fileData = sprintf('dataPhyto_upto1ps_int_sortdelay_unifdelay_nS%d_nBrg%d.mat',nS,nBrg);
fileData_LPF = sprintf('dataPhyto_upto1ps_int_sortdelay_unifdelay_LPF_nS%d_nBrg%d.mat',nS,nBrg);
% %if dataName == '3ps'
% fileData = 'dataPhyto_3ps_int_sortdelay_nS35366_nBrg61784.mat'; % to load
% fileData_LPF = 'dataPhyto_3ps_int_sortdelay_LPF_nS35366_nBrg61784.mat'; % to save

fileDir = ['./DATA_',dataName];  % the directory to save the output
%-----------------------------------------------------------------------------------------
FIL1 = fullfile(fileDir, fileData);
load(FIL1,'T','M','DRL','OSF','relB','miller_h','miller_k','delay',...
        'runID','eventID','miller_l','sort_notice','notice_negative_pix');

          
% %Computing Loorentz-Polarization factor for each Bragg reflection          
qvec = [miller_h./a, miller_k./b, miller_l./c];
q2   = qvec(:,1).^2 + qvec(:,2).^2 + qvec(:,3).^2; 
q    = sqrt(q2);
d    = 1./q;
find(isinf(d))
find(isnan(d))

% % Some notes on calcullation polarization factor:
% % n*lambda = 2*d*sin(theta)
% % Energy ~ 7000 eV => lambda=1.77 A
% % sin(theta) = (n*lambda)/(2*d)
% % cos(2*theta) = 1 - 2*sin(theta).^2

%lambda = 1.77;  % in Angstrom
sin_theta = lambda./(2*d);  % Bragg's law
find(isinf(sin_theta))
find(isnan(sin_theta))
cos_2theta = 1 - 2*sin_theta.^2;

LPF = (1 + cos_2theta.^2 )/2;

T_lpf = T./LPF';

% % % double-check:
% nS = 35366;
% kr = randperm(nS);
% k = kr(1);
% idx = find(T(k,:)~=0);
% nr = randperm(length(idx));
% n = nr(1);
% T(k,idx(n));
% T(k,idx(n))/P(idx(n))
% T_lpf(k,idx(n))
% % % % % % % % % % 

noticeLPF = 'Bragg intensities divided by Polarization factor';
FIL2 = fullfile(fileDir, fileData_LPF);
save(FIL2, 'T_lpf', 'M','miller_h','miller_k','miller_l','noticeLPF',...
                   'delay','DRL','notice_negative_pix','sort_notice',...
                   'LPF','OSF','relB','runID','eventID','-v7.3');  

%EOF
