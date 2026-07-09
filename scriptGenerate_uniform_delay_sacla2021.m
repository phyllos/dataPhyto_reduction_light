% scriptGenerate_uniform_delay_sacla2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Uniformizing delay-time distribution by distributing them into smaller bins.
% Adapted from R. Sepehr, UWM, 2016-2022 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear

set(0,'DefaultAxesFontSize',16);
set(0,'DefaultLineLineWidth', 2);
set(0,'defaultfigurecolor','w');

%% ------------------------------ Load data -------------------------------
dataName = 'upto1ps'; 
fileDir = ['./DATA_',dataName];  % the directory to save the output

fileData = 'dataPhyto_upto1ps_int_sortdelay_nS135041_nBrg62530.mat';
FIL1 = fullfile(fileDir, fileData);
load(FIL1)

%% --------------------- Original delay distribution ----------------------
D0 = unique(delay);
d0 = [];
for ii = 1:length(D0)
    bin0{ii} = delay(delay ==D0(ii));  
    d0 = [d0;bin0{ii}(:,1)];                 % double check
    l0(ii) = size(bin0{ii}(:,1),1);          % number of samples per bin
end

%% ------------- diagnostic plots: Original delay histogram ---------------
figure; bar(D0,l0,'FaceColor','b','EdgeColor','k','LineWidth',0.8)
title(['Phytochromes: histogram of unique values ',...
       'of delay times (',int2str(length(l0)),' unique delays)']);
xlabel('Delay time (fs)'); ylabel('# of snapshots');

%% ---------- Time window selection and taking data within that ----------- 
nS0 = length(delay);

idx = find(delay<550 & delay>-84);  % generates nS=107598 =>68,084(unif)
%idx = find(delay<947 & delay>647); % generates 8206, region around 800 fs

delay = delay(idx);
DRL = DRL(idx);
runID = runID(idx);
eventID = eventID(idx);
indSort_delay = indSort_delay(idx);
OSF = OSF(idx);
relB = relB(idx);
M = M(idx,:);
T = T(idx,:);

%% ---------------------- Group by unique delay values --------------------
ind_delay_uniform = 1:nS0;
ind_delay_uniform = ind_delay_uniform(idx);

D = unique(delay);
d = [];
for ii = 1:length(D)
    bin{ii} = delay(delay == D(ii));
    d = [d; bin{ii}(:,1)];          
    l(ii) = size(bin{ii}(:,1), 1);  % number of samples per bin
end

%% ------------------ diagnostic plots: Unique delays ---------------------
figure; plot(D,'.');
title('Unique delay values');
xlabel('# of snapshots');
ylabel('Delay time (fs)');

%% ----------------- Uniform subsampling per delay bin --------------------
ss = 155; %(for 68084); %for -84<t<550 % older version: floor(mean(l));
%ss = 65  %(for 8206)   %for the region aorund 800fs 

id_unif = [];
delay_uniform = [];
sum_lii = 0;
for ii = 1:length(l) 
    l(ii) = size(bin{ii}(:,1),1);  % bin size
    ind = (randperm(l(ii)));       % random index selection
    maxx = min(ss,l(ii));          % new bin size
    ind = sort(ind(1:maxx));       % sorting random indexes
    bin2{ii}(:,1) = bin{ii}(ind,1);% removing excessive delays in each bin
    delay_uniform = [delay_uniform;bin2{ii}(:,1)];  % reduced delay indexes
    id_unif = [id_unif;ind'+sum_lii]; % indexes that will be used 
    sum_lii = sum_lii+l(ii);
    l2(ii) =  size(bin2{ii}(:,1),1);  % bin size
end

%% --------------- diagnostic plots: Compare distributions ----------------
figure;
subplot(3,1,1:2);
bar(D,l,'FaceColor','b','EdgeColor','b','LineWidth',0.8)

ylabel('# of snapshots');
title(['Non-uniform delays: bar graph based on the number of unique ',...
       'delay times (',int2str(length(l)),' unique delays)']);
subplot(313);
%bar(l2,'b','LineWidth',2)
bar(D,l2,'BarWidth', 1,'FaceColor','b','EdgeColor','b','LineWidth',0.8)
ylim([0, ss+5])
title(['Flattened delays (',int2str(length(l2)),' unique delays)']);
xlabel('Delay time (fs)'); 

%% ---------------------- Final dataset assignment ------------------------
nS_unif = length(id_unif);

T = T(id_unif, :);
M = M(id_unif, :);
OSF = OSF(id_unif);
relB = relB(id_unif);
DRL = DRL(id_unif);
delay = delay(id_unif);
runID = runID(id_unif);
eventID = eventID(id_unif);
indSort_delay = indSort_delay(id_unif);
ind_delay_uniform = ind_delay_uniform(id_unif);

%% ---------------------------- Save output -------------------------------
file2 = ['dataPhyto_upto1ps_int_sortdelay_unifdelay_nS', int2str(nS_unif),...
         '_nBrg62530.mat'];         
FIL2 = fullfile(fileDir, file2);
save(FIL2,'T', 'M', 'miller_h', 'miller_k', 'miller_l', 'delay','idx',...
           'OSF','relB','DRL','indSort_delay','runID','eventID',...
          'sort_notice','ind_delay_uniform','notice_negative_pix', '-v7.3')

%% -------------------------------- End -----------------------------------
%EOF