% script_pack_phyto_final_light
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script reads data information (runID, eventID, number of Bragg spots) 
% from DAT files and read the redundancy file for the data as well. Then it 
% packs the data in a sparse matrix and generates a mask too. Results are 
% saved in a MAT file for analysis. 
% Ahmad H., Dec-07-17, updated on May-24-18 for fs data. updated July-03-18
% Updated on Feb-2022 for Phytochrome data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clear

tic
% % data ID
dataName = '3ps'; % '3ps' vs. '1ps'; % data name
dataForm = 'int';     % 'int' vs. 'amp'; (intensity vs. amplitude) 
sort_delay = true;    % sort based on delay time
remove_negative_pixels = false;  % decide to set negative pixels to zero or not
fileDir = ['./DATA_',dataName];  % the directory to save the output
  if ~exist(fileDir,'dir')
    system(['mkdir ' fileDir]);
  end
  
fid = sprintf('merged_snapshotInfo_phyto_%s_allInfo.dat',dataName);  % to load
filRedundancy = ['redundancy_',dataName,'.mat']; 

% "miller_h", "miller_k", "miller_l" are Miller indices
% they are intentionally created to have the same packing convention as "redundancy"
% if dataName == '3ps'
    h_max = 30;  
    k_max = 62;  
    l_max = 64; 
%  elseif dataName == 'upto1ps'
%    h_max = 29; 
%    k_max = 63; 
%    l_max = 67; 
% end
    
%                      
if sort_delay
   sort_notice = 'data sorted based on delay time';
else
   sort_notice = 'data NOT sorted'; 
end

fid = fopen(fid,'r');
snapshotInfo = fscanf(fid,'%f');
fclose(fid);

snapshotInfo = reshape(snapshotInfo,7,[])'; 
num_snapshot = size(snapshotInfo,1);

if ~exist('debug','var')
  debug = false;  % AHZ: was true
end
if debug
  num_snapshot = 10;
  debugfile = 'debug.dat';
  delete(debugfile);
  fileA = 'fileA.dat';
  fileB = 'fileB.dat';
end

% % sort snapshotInfo by delay, which is the 5-th column
[snapshotInfo, indSort_delay] = sortrows(snapshotInfo,5); 

runID   = snapshotInfo(:,1);   % run numbers
eventID = snapshotInfo(:,2);   % event numbers
nRefl   = snapshotInfo(:,3);   % num of Bragg reflections per snapshot
DRL     = snapshotInfo(:,4);   % diffraction resolution limit 
delay   = snapshotInfo(:,5);   % delay times
OSF     = snapshotInfo(:,6);   % Waller-Debye factor-1
relB    = snapshotInfo(:,7);   % Waller-Debye factor-2

% % check the sorted delay times:
% figure; plot(delay); ylabel('delay');

num_h = h_max+1;
num_k = k_max+1;
num_l = l_max+1;
miller_h = zeros(num_h,num_k,num_l);
miller_k = zeros(num_h,num_k,num_l);
miller_l = zeros(num_h,num_k,num_l);
num_miller = num_h*num_k*num_l;
for hh=0:h_max
  for kk=0:k_max
    for ll=0:l_max
      miller_h(hh+1,kk+1,ll+1) = hh;
      miller_k(hh+1,kk+1,ll+1) = kk;
      miller_l(hh+1,kk+1,ll+1) = ll;
    end
  end
end

% "redundany" has the same packing convention as "miller_h", "miller_k" and "miller_l"
load(filRedundancy,'redundancy')      

% % Specific for PYP (lines 94-102), RF Jan 2024:
% %%%%%%%%%%%%%%%%%%%%
% %%%kk=0;  not used in these loops!, AHZ
% for hh=1:h_max
%   for ll=0:l_max
%     redundancy_h0l = redundancy(hh+1, 0+1,ll+1); % the (h,0,l) and (0,h,l) reflections are symmetrically equivalent.                 
%     redundancy_0hl = redundancy( 0+1,hh+1,ll+1); % Thus, here we add them together to give the new redundancy for (h,0,l) 
%     redundancy(hh+1, 0+1,ll+1) = redundancy_h0l+redundancy_0hl;
%     redundancy( 0+1,hh+1,ll+1) = 0;  % and then the redundancy of (0,h,l) is set to zero here. 
%   end
% end
% %%%%%%%%%%%%%%%%%%%%

active_reflections = find(redundancy(:)>0);         % size: (61124,1)
num_unique_reflections = numel(active_reflections);  % i.e., active lattice vertices (Millers) that have
% Brag spots for "all" snapshots (notice, each snapshot partially contains a subset of active spots)

% follow the Glownia convention
T = zeros(num_snapshot,num_unique_reflections);
M = zeros(num_snapshot,num_unique_reflections);

%X_AHZ = cell(num_snapshot,1);
tic
for jj=1:num_snapshot
% "temp_vol" and "count" also follow the same packing convention
  temp_vol = zeros(num_h,num_k,num_l);   % size: (30,64,68) like redundancy file
  count = zeros(num_h,num_k,num_l);      % size: (30,64,68) like redundancy file
  temp_vec = nan(num_miller,1);          % size: (30x64x68,1) 
  runID_ = snapshotInfo(jj,1);            % run number for each snapshot
  eventID_ = snapshotInfo(jj,2);          % event number for each snapshot
  num_reflection = snapshotInfo(jj,3);   % Bragg' spot number for each snapshot
  filename = ['./sacla2021_',dataName,'_hklI/sacla2021_run'...
              num2str(runID_) '_tag' num2str(eventID_) '_hklI.dat'];
  fid_hklI = fopen(filename,'r');        
  hklI = fscanf(fid_hklI,'%f');          % read (h,k,l,I) for each snapshot from data file
  fclose(fid_hklI);
  hklI = reshape(hklI,4,num_reflection)';% make it as 2D matrix for all Bragg's points
  %K_AHZ = 0;
  for nn=1:num_reflection      % for each Brag spot in a typical snapshot:
    hh = hklI(nn,1);           % read h from 1st column
    kk = hklI(nn,2);           % read k from 2nd column
    ll = hklI(nn,3);           % read l from 3rd column
    intensity = hklI(nn,4);    % read I from 4th column

% % Specific for PYP (lines 136-139), RF Jan 2024:  
% %%%%%%%%%%%%%%%%%%%%
%     if ((hh==0) && (kk>0))     % if h=0 and k>0, set h=k and set k=0   
%       hh = kk;
%       kk = 0;
%       %K_AHZ = K_AHZ+1;
%     end
% %%%%%%%%%%%%%%%%%%%%
%again the same packing convention as redundancy. Like creating redundancy file,
%here we make a similar 3D lattice and add intensities of Bragg's points to
%vertices (h,k,l), and in fact we are doing this for all snapshots, so one 3D lattice 
%for intensities of all snapshots 
    temp_vol(hh+1,kk+1,ll+1) = temp_vol(hh+1,kk+1,ll+1)+intensity; 
% we count how many times a certain Bragg's point (lattice vertex) has been repeated when
% we add intensities on top of each other (needed for intensity normalization of vertices) 
    count(hh+1,kk+1,ll+1) = count(hh+1,kk+1,ll+1)+1;
  end % end-loop for a snapshot jj
  %X_AHZ{jj} = K_AHZ;
  
  count = count(:);              % size: (48*48*38,1)  
  M(jj,:) = (count(active_reflections)>0); % make a mask for snapshot jj by counting non-zero vertices 
  count(count==0) = 1;           % take care of Inf because of devision in next line
  temp_vec = temp_vol(:)./count; % point-wise intensity normalization for each vertex (Bragg's spot) in snapshot jj
  T(jj,:) = temp_vec(active_reflections);  % take the above normalized 3D lattice of intensities for snapshot jj
                                          % and only consider the total active spots as a snapshot TT(jj,:)                                                               
end
toc

miller_h = miller_h(:);
miller_h = miller_h(active_reflections);
miller_k = miller_k(:);
miller_k = miller_k(active_reflections);
miller_l = miller_l(:);
miller_l = miller_l(active_reflections);
T = sparse(T);
M = sparse(M);

if remove_negative_pixels
 T(T<0) = 0;   % AHZ, July-2018
 M(T<0) = 0;   % AHZ, July-2018
 notice_negative_pix = 'negative pixels set to zero';
else
 notice_negative_pix = 'negative pixels servived';
end
 
if strcmp(dataForm,'amp')
    T = sqrt(T);  % convert to amplitudes
    disp('-- data format: Amplitude --');
elseif strcmp(dataForm,'int')
    disp('-- data format: Intensity --')
end

% to save
fileName = sprintf('dataPhyto_%s_%s_sortdelay_nS%d_nBrg%d.mat',dataName, ...
                    dataForm, num_snapshot,num_unique_reflections);         
FIL = fullfile(fileDir, fileName);
save(FIL,'T','M','miller_h','miller_k','miller_l','delay','DRL',...
              'OSF','relB','sort_notice','runID','eventID',...
              'indSort_delay','notice_negative_pix', '-v7.3')      
toc
%EOF
