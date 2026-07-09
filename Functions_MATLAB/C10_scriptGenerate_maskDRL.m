% scriptGenerate_maskDRL
% Feb-2022

clear 

dataName = 'upto1ps'; %'3ps'; %'puredark-clean'

fileDir = ['./DATA_',dataName];  % the directory to save the output

filename_inp = 'dataPhyto_upto1ps_int_sortdelay_unifdelay_LPF_nS68084_nBrg62530.mat'; % to load info
fileMask = 'maskDRL_upto1ps_int_nS68084_nBrg62530.mat';  % to save the mask file

% filename_inp = 'dataPhyto_3ps_int_sortdelay_LPF_nS35366_nBrg61784.mat'; % to load info
% fileMask = 'maskDRL_3ps_int_nS35366_nBrg61784.mat';  % to save the mask file

% filename_inp = 'dataPhyto_puredark-clean_int_sortEvent_LPF_nS70395_nBrg62123.mat'; % to load info
% fileMask = 'maskDRL_puredark-clean_nS70395_nBrg62123.mat';  % to save the mask file

FIL1 = fullfile(fileDir, filename_inp);

note = 'this mask removes Bragg points beyond DRL. Also, qmax=1/DRL';

load(FIL1,'DRL','miller_h','miller_k','miller_l');
DRL(DRL == 0) = mean(DRL); % to get rid of possible NaN values  
%

disp('--- Careful if data have been sorted or not! ---');

fileParams = matfile(FIL1);
[nS , nBrg] = size(fileParams,'T_lpf');

maskDRL = ones(nS,nBrg);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % removing spots beyond DRL in reconst data
% lattice sizes for group P212121:
a = 54.95; %54.22;  % in Angstrom
b = 116.5; %115.78; % in Angstrom
c = 117.7; %117.08; % in Angstrom

%DRL = 1.5*ones(nS,1); %test on April-07-2022 for a constant DRL for all
qmax = 1./DRL; % notice that DRL is in Angstrom

% q-vector in reciprocal space:
%qvec = [miller_h./a, miller_h./(sqrt(3)*a) + 2*miller_k./(sqrt(3)*b),miller_l./c]; PYP 

qvec = [miller_h./a, miller_k./b, miller_l./c];
q2   = qvec(:,1).^2 + qvec(:,2).^2 + qvec(:,3).^2; 
q    = sqrt(q2);

% find reflections with better resolution than DRL for each snapshot
tic
for ii = 1:nS
    IND = find(q > qmax(ii));
    maskDRL(ii,IND) = 0;      % AH
end
toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
maskDRL = sparse(maskDRL);
FIL2 = fullfile(fileDir, fileMask);
save(FIL2,'maskDRL','DRL','qmax','note','-v7.3');








