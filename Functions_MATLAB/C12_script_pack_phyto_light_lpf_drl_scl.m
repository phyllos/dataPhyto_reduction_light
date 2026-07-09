% script_pack_phyto_light_lpf_drl_scl
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script reads Partialator parameters of crystFEL, i.e., OSF & relB, 
% for PYP data and then multiply them onto the snapshots, which are already
% passed the DRL cut step.
% Ahmad H., Dec-29-2017, May 2018, July-03-2018
% update for phytochromes on Feb-2022, June 2022
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
 
generate_hkl_avg = true;
remove_negative_pixels = false;  
STD_factor = 10; % RS took 10

% lattice sizes for group P6(3):
a = 54.95; %54.22;  % in Angstrom
b = 116.5; %115.78; % in Angstrom
c = 117.7; %117.08; % in Angstrom

% % if dataName == '3ps'
% file2 = 'dataPhyto_3ps_int_sortdelay_LPF_DRL_nS35366_nBrg61784.mat'; %to load             
% % 
% file3 = 'dataPhyto_3ps_int_sortdelay_LPF_DRL_SCL_nS35366_nBrg61784.mat'; %to save
% fileHKL = 'dataPhyto_3ps_int_sortdelay_LPF_DRL_SCL_AVG_nS35366_nBrg61784.hkl'; %to save

% % % if dataName == '1ps'
% %file2 = sprintf('merged_snapshotInfo_phyto_%s_allInfo.dat',dataName);  %to load             
 file2 = sprintf('dataPhyto_upto1ps_int_sortdelay_unifdelay_LPF_DRL_nS%d_nBrg%d.mat',nS,nBrg); %to load          
 file3 = sprintf('dataPhyto_upto1ps_int_sortdelay_unifdelay_LPF_DRL_SCL_nS%d_nBrg%d.mat',nS,nBrg); %to save
 fileHKL = sprintf('dataPhyto_upto1ps_int_sortdelay_unifdelay_LPF_DRL_SCL_AVG_nS%d_nBrg%d.hkl',nS,nBrg); %to save  
        
fileDir = ['./DATA_',dataName];  % the directory to save the output

%--------------------------------------------------------------------------
% load snapshots with resolution cuts beyond DRL 
FIL2 = fullfile(fileDir, file2);
load(FIL2,'T_lpf_drl','M_drl','miller_h','miller_k','miller_l','DRL','LPF',...
     'OSF','relB','delay','noticeDRL','notice_negative_pix', 'sort_notice');

%--------------------------------------------------------------------------
% % strategy suggested by RS:
% blocked on Feb-2022
% OSF(OSF>5)=5;OSF(OSF<-2.5)=-2.5;     %numbers selected based on plot(relB)
% relB(relB>50)=50;relB(relB<-50)=-50; %numbers selected based on plot(relB)

% blocked on Feb-2022:
% relB(relB>20)=20;relB(relB<-60)=-60; %numbers selected based on plot(relB)

qvec = [miller_h./a, miller_k./b, miller_l./c];
q2   = qvec(:,1).^2 + qvec(:,2).^2 + qvec(:,3).^2; 
%notice: sinTheta_lambda2 = q2/4 = 1/(4*d2) 

T_lpf_drl_scl = sparse(size(T_lpf_drl,1),size(T_lpf_drl,2));
blocksize = 10000;  % for large matrix multiplication
N = blocksize;

for i = 1: ceil(size(T_lpf_drl,1)/N)
    i
  if i*N > size(T_lpf_drl,1)
    T_lpf_drl_scl(1+(i-1)*N:end,:) = T_lpf_drl(1+(i-1)*N:end,:) ...
    .*repmat(exp(-OSF(1+(i-1)*N:end)),1,nBrg).* exp(relB(1+(i-1)*N:end)*q2'/4); 
  else
    T_lpf_drl_scl(1+(i-1)*N:i*N,:) = T_lpf_drl(1+(i-1)*N:i*N,:)...
    .*repmat(exp(-OSF(1+(i-1)*N:i*N)),1,nBrg).* exp(relB(1+(i-1)*N:i*N)*q2'/4); 
  end  
end
noticeSCL = 'Data scaled by partialator parameters';

FIL3 = fullfile(fileDir, file3);
save(FIL3,'T_lpf_drl_scl','M_drl','LPF','miller_h','miller_k','miller_l', ...
'delay','noticeDRL','noticeSCL','sort_notice','notice_negative_pix','-v7.3'); 

% %________________________________________________________________________
% % Finding average and std and save in a *hkl file 
if generate_hkl_avg
T_avg = full(sum(T_lpf_drl_scl)./sum(M_drl));
T_std = zeros(1,nBrg);
for i = 1:nBrg
    T_std(1,i) = std(nonzeros(T_lpf_drl_scl(:,i)));  % Not sure yet!!
end

% % set NaN, Inf and negative values to zero
T_avg(isnan(T_avg)) = 0;
T_std(isnan(T_std)) = 0;
T_avg(isinf(T_avg)) = 0;
T_std(isinf(T_std)) = 0;
if remove_negative_pixels
 T_avg(T_avg<0) = 0;
 T_std(T_avg == 0) =0;
end

T_avg(T_avg<0) = 0;     % added by AHZ, March-2022
T_STD_new = sqrt(T_avg);

filHKL = fullfile(fileDir, fileHKL);
fileID = fopen(filHKL,'w');  
%fprintf(fileID,'%3d %3d %3d %6.6f %6.6f\n',[miller_h,miller_k,miller_l,T_avg',T_std'/STD_factor]');
fprintf(fileID,'%6d %6d %6d %12.2f %12.2f\n',[miller_h,miller_k,miller_l,T_avg',T_STD_new']');
fclose(fileID);
end
%EOF
