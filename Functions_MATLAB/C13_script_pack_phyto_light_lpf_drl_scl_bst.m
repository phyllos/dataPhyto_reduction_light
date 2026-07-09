% script_pack_phyto_light_lpf_drl_scl_wgt_bst
% Correct for multiplicities
% April-2022, June-2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear

dataName = 'upto1ps'; %'3ps'
nBrg = 62530; %3ps: 61784;  % number of active Bragg spots
nS   = 107404;%24726;%8206;%68084; %(uniformed) %

fileDir = ['./DATA_',dataName];  % the directory to save the output

file1 = sprintf('--dataPhyto_upto1ps_int_sortdelay_LPF_DRL_SCL_nS%d_nBrg%d.mat',nS,nBrg);
file2 = sprintf('--dataPhyto_upto1ps_int_sortdelay_LPF_DRL_SCL_BST_nS%d_nBrg%d.mat',nS,nBrg);

FIL1 = fullfile(fileDir, file1);
load(FIL1)

[nS, nBrg]=size(T_lpf_drl_scl);

T_lpf_drl_scl_bst = T_lpf_drl_scl;

% Boosting, RS version:
for k=1:nBrg
    co(k) = sum(M_drl(:,k));
    if co(k)
        T_lpf_drl_scl_bst(:,k) = T_lpf_drl_scl(:,k)./co(k);
    end
end

noticeBST = 'boosting applied to T_lpf_drl_scl';


FIL2 = fullfile(fileDir, file2);
save(FIL2,'T_lpf_drl_scl_bst','M_drl','miller_k','miller_h','miller_l',...
          'delay','sort_notice','noticeBST','notice_negative_pix','-v7.3')

% % lattice sizes for group P212121:
% a = 54.22;  % in Angstrom
% b = 115.78; % in Angstrom
% c = 117.08; % in Angstrom
% qvec = [miller_h./a, miller_k./b, miller_l./c];
% q2   = qvec(:,1).^2 + qvec(:,2).^2 + qvec(:,3).^2; 
% q = sqrt(q2);
% [qs, ind]=sort(q);
% figure; plot(qs,co(ind));

%EOF
