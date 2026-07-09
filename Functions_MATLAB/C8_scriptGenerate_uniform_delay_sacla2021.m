% scriptGenerate_uniform_delay_sacla2021
% Adapted from RS code named as "Rey2016_1_PYPprePro_delay".
% A. Hosseini, May-2018, Feb & Nov 2022,Jan 2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % distributing the delay time into smallest bins.
clear

set(0,'DefaultAxesFontSize',20);
set(0,'DefaultLineLineWidth', 2);
set(0,'defaultfigurecolor','w');
%
dataName = 'upto1ps'; 
fileDir = ['./DATA_',dataName];  % the directory to save the output

fileData = 'dataPhyto_upto1ps_int_sortdelay_nS135041_nBrg62530.mat';
FIL1 = fullfile(fileDir, fileData);
load(FIL1)
%
% figure;
% histogram(delay,1000);
% title(['Phytochromes: histogram of delay times'])
% xlabel('Delay time (fs)'); ylabel('# of snapshots');
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ 
% % Added Nov-11 for new plotting based on unique delays:
D0 = unique(delay);
d0 = [];
for ii = 1:length(D0)
    bin0{ii} = delay(find(delay ==D0(ii)));  % changed bin to bin0 (Jan 2026)
    d0 = [d0;bin0{ii}(:,1)];          %%%% double check
    l0(ii) = size(bin0{ii}(:,1),1);  %%%% number of samples per bin
end
figure; bar(D0,l0,'FaceColor','b','EdgeColor','b','LineWidth',0.8)
%xticks(D0(1):100:D0(end))
%xlim([D0(1)-20 D0(end)+20])
title(['Phytochromes: bar graph based on the number of unique values',...
       'of delay times (',int2str(length(l0)),' unique delays)']);
xlabel('Delay time (fs)'); ylabel('# of snapshots');
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ 
% % taking subset time-intervals:

% %idx = find(delay<440 & delay>-80);    % generates nS=97817 =>59,988(unif)
idx = find(delay<550 & delay>-84);      % generates nS=107598 =>68,084(unif)

%idx = find(delay<947 & delay>647); % generates 8206, region around 800 fs
%idx = find(delay<940 & delay>-84);  % generates 24726, full range (Oct-22)

nS = size(T,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% % take the data in the above delay time interval
delay = delay(idx);
DRL = DRL(idx);
runID = runID(idx);
eventID = eventID(idx);
indSort_delay = indSort_delay(idx);
OSF = OSF(idx);
relB = relB(idx);
M = M(idx,:);
T = T(idx,:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ind_delay_uniform = 1:nS;
ind_delay_uniform = ind_delay_uniform(idx);
%
%
D = unique(delay);
d = [];
for ii = 1:length(D)
    bin{ii} = delay(find(delay ==D(ii)));
    d = [d;bin{ii}(:,1)];          %%%% double check
    l(ii) = size(bin{ii}(:,1),1);  %%%% number of samples per bin
end
% figure;plot(D,'.');
% figure;plot(delay);hold on;plot(d,'-.r');

ss = 155 %(for 68084); %the region of -84<t<550 %an old one floor(mean(l));
%ss = 65  %(for 8206)   %the region aorund 800fs 
%ss = 32   %(for 24726)  %the full region of -84<t<940 

% Reassigning the delay times to make a more uniform distribusion
id_unif = [];
delay_uniform = [];
sum_lii = 0;
for ii = 1:length(l) %619
    l(ii) = size(bin{ii}(:,1),1);  % bin size
    ind = (randperm(l(ii)));       % random index selection
    maxx = min(ss,l(ii));          % new bin size
    ind = sort(ind(1:maxx));       % sorting random indexes
    bin2{ii}(:,1) = bin{ii}(ind,1);% removing excessive delay times in each bin
    delay_uniform = [delay_uniform;bin2{ii}(:,1)];  %%% new reduced delay index
    id_unif = [id_unif;ind'+sum_lii]; % indexes that have been used 
    sum_lii = sum_lii+l(ii);
    l2(ii) =  size(bin2{ii}(:,1),1);  % bin size
end
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ 
% % Older style of plotting
% figure;
% plot(delay,'-b'); 
% hold on;plot(delay_uniform,'-r'); title('Original vs. flattened delays');
% legend('Non-uniform delays','Flattened delays','Location','northwest')
% %
% figure;
% subplot(3,1,1:2);
% hist(delay,2000);  
% ylabel('# of snapshots'); title('Original delays')    
% subplot(313);
% hist(delay_uniform,2500);
% xlabel('Delay time (fs)'); title('Flattened delays');
% %
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ 
% % Added Nov-11 for new plotting based on unique delays:
figure;
subplot(3,1,1:2);
%bar(l,'b'); 
bar(D,l,'FaceColor','b','EdgeColor','b','LineWidth',0.8)

ylabel('# of snapshots');
title(['Non-uniform delays: bar graph based on the number of unique ',...
       'delay times (',int2str(length(l)),' unique delays)']);
subplot(313);
%bar(l2,'b'); 
bar(D,l2,'FaceColor','b','EdgeColor','b','LineWidth',0.8)
ylim([0, ss+5])
title(['Flattened delays (',int2str(length(l2)),' unique delays)']);
xlabel('Delay time (fs)'); 
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ 
%
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% to save
file2 = ['dataPhyto_upto1ps_int_sortdelay_unifdelay_nS', int2str(nS_unif),...
         '_nBrg62530.mat'];         
FIL2 = fullfile(fileDir, file2);
save(FIL2,'T', 'M', 'miller_h', 'miller_k', 'miller_l', 'delay','idx',...
           'OSF','relB','DRL','indSort_delay','runID','eventID',...
          'sort_notice','ind_delay_uniform','notice_negative_pix', '-v7.3')
