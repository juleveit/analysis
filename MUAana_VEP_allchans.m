function MUAana_VEP_allchans

animalid = '180619';
block = 2;
basepath = ['C:\Users\Julia\work\data\' animalid '\'];
channels = 1:16;
recdepth = 600;
spacing = 25;
penangle = 25;

recalculate = 0;

prestim = 3000;
poststim = 500;
trialdur = 5000;
respwin = 1:5000; %501:1500; % after stimulus onset
respwin = respwin+prestim;
muaresultfile = [basepath, 'block' int2str(block) 'MUAresult.mat'];

depth = recdepth.*(cosd(penangle)*cosd(22));
spacing = spacing.*(cosd(penangle)*cosd(22));
evaldepth = 400;
for i = 1:length(depth)
    for j = 1:length(channels(1):channels(end))
        dm(i,j) = depth(i)-((j-1)*spacing(i)); % depth matrix
    end
    [c,di(i)] = min(abs(dm(i,:)-evaldepth)); % get depth index, least distance to evaldepth
end     

if ~exist(muaresultfile) || recalculate
    result = MUAdataprepareNOPTB(basepath,animalid,block,channels);
    save(muaresultfile,'result')
else
    load(muaresultfile);
end

wo = 6.5/(1000/2); bw = wo/1.5;
[b,a] = iirnotch(wo,bw);
binwidth = 50;
msstamps = result.msstamps;
chan = zeros(size(result.lfp));

for ch = 1:16    
    msstimes = round(result.msStimes{ch});
    chan(ch,msstimes) = 1;
    filtlfp = filtfilt(b,a,result.lfp(ch,:));
    for i = 1:length(msstamps)
        lfpresp(ch,i,:) = result.lfp(ch,msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim);
        filtlfpresp(ch,i,:) = filtlfp(msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim);
        resp(ch,i,:) = chan(ch,msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim);
        [binresp(ch,i,:),bta] = binit(squeeze(resp(ch,i,:)),binwidth);
    end
end
ta = -prestim+1:trialdur+poststim;
binresp = binresp.*(1000/binwidth);
bta = bta-prestim;

ch = 2;
figure
plot(ta,squeeze(lfpresp(ch,:,:))');
hold on
plot(ta,squeeze(mean(lfpresp(ch,:,:),2)),'k','linewidth',2)

figure
plot(ta,squeeze(filtlfpresp(ch,:,:))');
hold on
plot(ta,squeeze(mean(filtlfpresp(ch,:,:),2)),'k','linewidth',2)

figure
plot(bta,squeeze(nanmean(binresp(ch,:,:),2)),'k','linewidth',2)
xlabel('time (ms)')
ylabel('firing rate (Hz')

figure
plot(ta,squeeze(nanmean(filtlfpresp(ch,:,:),2)),'k','linewidth',2)
xlabel('time (ms)')

% imagescs
figure
imagesc(bta,dm,squeeze(nanmean(binresp,2)))
xlabel('time (ms)')
ylabel('cortical depth (um)')

figure
imagesc(ta,dm,squeeze(nanmean(filtlfpresp,2)))
xlabel('time (ms)')
ylabel('cortical depth (um)')

figure
hold on
for i = 1:16
plot(bta,squeeze(nanmean(binresp(i,:,:),2))+(i-1)*20,'k')
end
xlabel('time (ms)')

figure
hold on
for i = 1:16
plot(ta,squeeze(nanmean(filtlfpresp(i,:,:),2))+(i-1)*70,'k')
end
xlabel('time (ms)')
% figure
% errorbar(ta,squeeze(mean(lfpresp(ch,:,:),2)),squeeze(std(lfpresp(ch,:,:),1,2))./sqrt(size(lfpresp,2)));
% 
% figure
% errorbar(bta,nanmean(binresp),nanstd(binresp)./sqrt(length(msstamps)))

 disp('');

