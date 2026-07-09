% find_hkl_redundancy_sacla2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%snapshotInfo_file = {'merged_snapshotInfo_common.dat'};
snapshotInfo_file = {'3ps_snapshotInfo.dat'};
redundancy_file_name = 'redundancy_3ps.mat';
num_snapshotInfo_file = length(snapshotInfo_file);

h_max = 30;  % % 1ps:29, 3ps:30, dark_All:31, dark_pure:30, dark-clean: 29
k_max = 62;  % % 1ps:63, 3ps:62, dark_All:62, dark_pure:62, dark-clean: 61
l_max = 64;  % % 1ps:67, 3ps:64, dark_All:66, dark_pure:64, dark-clean: 62
% % NOTICE: 
% % For dark data whether "All" or "pure", we use above hkl_max from "pure" 

redundancy = zeros((h_max+1),(k_max+1),(l_max+1));
for nn=1:num_snapshotInfo_file
  nn;
  fid_snapshotInfo = fopen(snapshotInfo_file{nn},'r');
  snapshotInfo = fscanf(fid_snapshotInfo,'%d');
  snapshotInfo = reshape(snapshotInfo,3,[])';
  num_snapshot = size(snapshotInfo,1);
  for jj=1:num_snapshot
    runID = snapshotInfo(jj,1);
    eventID = snapshotInfo(jj,2);
    num_reflection = snapshotInfo(jj,3);
    %filename = ['cxig0915_hklI/cxig0915_' num2str(runID) '_' num2str(eventID) '_hklI.dat'];
    % % for dark data whether All or pure, we use "sacla2021_dark_All_hklI"
    % directory, but "snapshotInfo_file" should be specified either for All
    % or pure (not the same for both!) as done in the beginning of this script 
    filename = ['./sacla2021_3ps_hklI/sacla2021_run' num2str(runID) '_tag' num2str(eventID) '_hklI.dat'];
    fid_hklI = fopen(filename,'r');
    hklI = fscanf(fid_hklI,'%f');
    hklI = reshape(hklI,4,num_reflection)';
    for kk=1:num_reflection
      h = hklI(kk,1);
      k = hklI(kk,2);
      l = hklI(kk,3);
      redundancy(h+1,k+1,l+1) = redundancy(h+1,k+1,l+1)+1;
    end
    fclose(fid_hklI);
  end
  fclose(fid_snapshotInfo);
end
save(redundancy_file_name,'redundancy')

return

% addpath('../PYP-work-space/export_fig/'); 
% redundancy(redundancy<1)=nan;
% for l=0:10:l_max
%   h = figure(1);
%   set(h,'color','w','resize','off')
%   pos = get(h,'pos');
%   set(h,'pos',[pos(1) pos(2) pos(3) 1.2*pos(4)])
%   hsp = subplot(1,1,1);
%   pcolor([0:h_max],[0:k_max],squeeze(redundancy(:,:,l+1))')
%   axis equal
%   set(hsp,'xlim',[0 h_max],'ylim',[0 k_max])
%   set(hsp,'fontSize',15,'lineWidth',2)
%   xlabel('h','fontSize',15)
%   ylabel('k','fontSize',15)
%   title(['Redundancy,l=' num2str(l)],'fontSize',20)
%   colorbar('fontSize',15)
%   export_fig('-jpeg','-r200',['redundancy_l_' num2str(l) '.jpg'])
%   close(1)
% end
