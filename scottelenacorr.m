recalculate = 0;
sr = 1000;

params.tapers = [2,5]; params.Fs = 1000; params.err = [2, 0.05]; params.trialave = 1;

recname = '1240';
foldername = 'ID1240-1241';

% if ~exist(['C:\Users\Julia\work\data\others\elenascottLFP\' recname '\filtered.mat']) || recalculate
% %     load(['C:\Users\Julia\work\data\others\elenascottLFP\' recname '\' foldername '.phy'], '-mat');
% %     load(['C:\Users\Julia\work\data\others\elenascottLFP\' recname '\' foldername '.ttl'], '-mat');
%     load(['C:\Users\Julia\work\data\others\elenascottLFP\' recname '\' foldername '.mat']);
%     for tr = 1:length(MCdata)
%         for chan = 1:size(MCdata{1},2)
%             lfptrial = eegfilt(MCdata{tr}(:,chan)',30000,0,200)';
%             lfp(chan,tr,:) = resample(lfptrial,1,30);
%             gamma(chan,tr,:) = eegfilt(lfp(chan,tr,:),sr,30,90);
%             h = hilbert(gamma(chan,tr,:)); gpow(chan,tr,:) = abs(h);
%         end
%     end
%     save(['C:\Users\Julia\work\data\others\elenascottLFP\' recname '\filtered.mat'],'lfp','gamma','gpow','stimsequence');
% else
%     load(['C:\Users\Julia\work\data\others\elenascottLFP\' recname '\filtered.mat']);
% end
respwin = 1501:2250;
load('C:\Users\Julia\work\data\others\elenascottLFP\SillehSom.mat')

if ~exist(['C:\Users\Julia\work\data\others\elenascottLFP\sillehlfp.mat']) || recalculate
    for tr = 1:length(LFPz)
        for chan = 1:size(LFPz{1},2)-1
            lfptrial = eegfilt(LFPz{tr}(:,chan)',30000,0,200)';
            lfp(chan,tr,:) = resample(lfptrial,1,30);
            [lfpspect(chan,tr,:),fax] = mtspectrumc(squeeze(lfp(chan,tr,respwin))',params);
        end
    end
    save(['C:\Users\Julia\work\data\others\elenascottLFP\sillehlfp.mat'],'lfp','lfpspect','stim');
else
    load(['C:\Users\Julia\work\data\others\elenascottLFP\sillehlfp.mat']);
end

stims = unique(stim);
for chan = 1:32
    for i = 1:2 % no light, light
        for j = 1:2 % no touch, touch
            if i == 1 & j == 1
                thisinds = find(stim == 9);
            elseif i == 1 & j == 2
                thisinds = find(stim == 4);
            elseif i == 2 & j == 1
                thisinds = find(stim == 18);
            else
                thisinds = find(stim == 13);
            end
            condlfpresp(chan,i,j,:) = nanmean(lfp(chan,thisinds,:),2);
            condlfpspect(chan,i,j,:) = nanmean(lfpspect(chan,thisinds,:),2);
        end
    end
end

for i = [8,20,3,21,11,28,16,22,4,29,10,30,12,25,2,23,5,17,9,24,13,31,1,32,15,18,7,26,14,19,6,27]
    figure
    semilogy(fax,squeeze(condlfpspect(i,1,1,:)),'linewidth',2)
    hold on
    semilogy(fax,squeeze(condlfpspect(i,1,2,:)),'c','linewidth',2)
    semilogy(fax,squeeze(condlfpspect(i,2,1,:)),'m','linewidth',2)
    semilogy(fax,squeeze(condlfpspect(i,2,2,:)),'r','linewidth',2)
    legend('L0 T0','L0 T1','L1 T0','L1 T1')
end
        

%normal trimming or schmenetics
ntrials = length(stimsequence)/2;
stimsequence(ntrials+1:end) = stimsequence(ntrials+1:end)+9;
l1 = find(stimsequence>9); 
l0 = find(stimsequence<=9);

% sham trimming
ntrials = length(stimsequence)/3;
stimsequence(ntrials+1:2*ntrials) = stimsequence(ntrials+1:2*ntrials)+9;
stimsequence(ntrials*2+1:end) = stimsequence(ntrials*2+1:end)+18;
l1 = find(stimsequence>18); 
l0 = find(stimsequence<=18&stimsequence>9);


% plot(squeeze(mean(gpow(9,l0(find(l0>67)),:),2)))
% hold on
% plot(squeeze(mean(gpow(9,l1(find(l1>67)),:),2)),'r')

plot(squeeze(mean(gpow(32,l0,:),2)))
hold on
plot(squeeze(mean(gpow(32,l1,:),2)),'r')


% l0p = 6; l1p = 15;
% chan = 10;
% l0inds = find(stimsequence == l0p); l0inds(find(l0inds<67)) = [];
% l1inds = find(stimsequence == l1p); l1inds(find(l1inds<67)) = [];
chan = 32;

% figure
% plot(squeeze(mean(gpow(chan,l0inds,:),2)));
% hold on
% plot(squeeze(mean(gpow(chan,l1inds,:),2)),'r');

% window = 1700:2211;
window = 1500:2400;
trialnfft = 2^nextpow2(500);
trialfax = sr/2*linspace(0,1,trialnfft/2+1);
for chan = 1:size(lfp,1)-1
%     for i = 68:size(lfp,2)
    for i = 1:size(lfp,2)
        y = fft(lfp(chan,i,window),trialnfft);
        lfpspect(chan,i,:) = abs(y(1:trialnfft/2+1));
    end    
end

l0p = 11; l1p = 20;
l0inds = find(stimsequence == l0p); %l0inds(find(l0inds<67)) = [];
l1inds = find(stimsequence == l1p); %l1inds(find(l1inds<67)) = [];
% l0inds = find(stimsequence == 2 | stimsequence == 11);
% l1inds = find(stimsequence == 9 | stimsequence == 18); 
% for i = [6,11,3,14,1,16,2,15,5,12,4,13,7,10,8,9] % 16 channel
for i = [8,20,3,21,11,28,16,22,4,29,10,30,12,25,2,23,5,17,9,24,13,31,1,32,15,18,7,26,14,19,6,27]
% for i = [5]
    figure
    semilogy(trialfax,squeeze(mean(lfpspect(i,l1inds,:),2)),'r','linewidth',2)
    hold on
    semilogy(trialfax,squeeze(mean(lfpspect(i,l0inds,:),2)),'linewidth',2)
    semilogy(trialfax,squeeze(mean(lfpspect(i,l1inds,:),2))-(squeeze(std(lfpspect(i,l1inds,:),1,2))./sqrt(length(l1inds))),'r')
    semilogy(trialfax,squeeze(mean(lfpspect(i,l1inds,:),2))+(squeeze(std(lfpspect(i,l1inds,:),1,2))./sqrt(length(l1inds))),'r')
    semilogy(trialfax,squeeze(mean(lfpspect(i,l0inds,:),2))-(squeeze(std(lfpspect(i,l0inds,:),1,2))./sqrt(length(l0inds))))
    semilogy(trialfax,squeeze(mean(lfpspect(i,l0inds,:),2))+(squeeze(std(lfpspect(i,l0inds,:),1,2))./sqrt(length(l0inds))))
    ax = axis;
    axis([0,120,ax(3),ax(4)])
end