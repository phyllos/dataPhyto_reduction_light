% find_hkl_max_sacla2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%snapshotInfo_file = {'merged_snapshotInfo_common.dat'};
binName = '3ps';
snapshotInfo_file = {'3ps_snapshotInfo.dat'};
num_snapshotInfo_file = length(snapshotInfo_file);

h_min = inf;
h_max = -inf;
k_min = inf;
k_max = -inf;
l_min = inf;
l_max = -inf;

for nn = 1:num_snapshotInfo_file
  fid_snapshotInfo = fopen(snapshotInfo_file{nn},'r');
  snapshotInfo = fscanf(fid_snapshotInfo,'%d');
  snapshotInfo = reshape(snapshotInfo,3,[])';
  num_snapshot = size(snapshotInfo,1);
  for jj=1:num_snapshot
    runID = snapshotInfo(jj,1);
    eventID = snapshotInfo(jj,2);
    num_reflection = snapshotInfo(jj,3);
    %filename = ['cxig0915_hklI/cxig0915_' num2str(runID) '_' num2str(eventID) '_hklI.dat'];
    %filename = ['sacla2021_dark-clean_hklI/sacla2021_run' num2str(runID) '_tag' num2str(eventID) '_hklI.dat'];
    filename = ['sacla2021_',binName,'_hklI/sacla2021_run' num2str(runID) '_tag' num2str(eventID) '_hklI.dat'];

    fid_hklI = fopen(filename,'r');
    hklI = fscanf(fid_hklI,'%f');
    hklI = reshape(hklI,4,num_reflection)';
    for kk=1:num_reflection
      h = hklI(kk,1);
      k = hklI(kk,2);
      l = hklI(kk,3);
      if (h<h_min), h_min = h; end
      if (h>h_max), h_max = h; end
      if (k<k_min), k_min = k; end
      if (k>k_max), k_max = k; end
      if (l<l_min), l_min = l; end
      if (l>l_max), l_max = l; end
    end
    fclose(fid_hklI);
  end
  fclose(fid_snapshotInfo);
end

h_min
h_max
k_min
k_max
l_min
l_max
