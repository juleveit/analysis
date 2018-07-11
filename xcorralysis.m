function xcorralysis

animalid = '140618';
block = 3;
tetrode = 6;

cell1 = 6;
cell2 = 2;

freqbinwidth = 5;
sr = 1000;

supath = ['C:\Users\Julia\work\data\' animalid '\singleunits\'];
basename = [animalid '_block' int2str(block) '_tet' int2str(tetrode)];

files = dir([supath, basename, '*.mat']);

for i = 1:length(files)
    if ~isempty(regexp(files(i).name, [basename '_cell' int2str(cell1) '_']))
        in(1) = i;
    elseif ~isempty(regexp(files(i).name, [basename '_cell' int2str(cell2) '_']))
        in(2) = i;
    end
end

prestim = 0;
poststim = 0;
respwin = 501:1500; % after stimulus onset
respwin = respwin+prestim;
  
for cell = 1:2
    load([supath, files(in(cell)).name]);

    % get spiketimes
    msStimes = round(result.spikes);
    if ~isempty(msStimes) & msStimes(1) == 0, msStimes(1) = 1; end
    chan = zeros(1,length(result.lfp));
    chan(msStimes) = 1;     
    trialdur = result.stimduration*1000;
    msstamps = result.msstamps;

    wvchan = find(var(result.waveforms) == max(var(result.waveforms)));
    lfp = result.lfp(:,wvchan)';
    for i = 1:100/freqbinwidth
        filtmat(i,:) = eegfilt(lfp,sr,(i-1)*freqbinwidth+1,i*freqbinwidth);
        h = hilbert(filtmat(i,:));
        powmat(i,:) = abs(h); phasmat(i,:) = angle(h);
    end

    for i = 1:length(msstamps)
        resp(i,:) = chan(msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim);
        lfpresp(i,:) = result.lfp(msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim,wvchan);
        [lfpspect(i,:),fax] = pmtm(lfpresp(i,respwin),2,[],1000);

        for j = 1:size(phasmat,1)
            allfiltresp(j,i,:) = filtmat(j, msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim);
            allphaseresp(j,i,:) = phasmat(j, msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim);
            allpowresp(j,i,:) = powmat(j, msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim);
        end
    end
    
    l0 = find(result.light == 0);
    long(cell,:) = reshape(resp(l0,:)',1,trialdur*length(l0));

    tmpi = zeros(size(allphaseresp));
    for i = 1:size(allphaseresp,1)
        tmp = zeros(size(resp));
        tmp(find(resp)) = allphaseresp(i,find(resp));
        tmpi(i,:,:) = tmp;
    end
    phasemat = tmpi(:,l0,:);
    for b = 1:size(allphaseresp,1)
        condpowresp(b,:) = squeeze(mean(allpowresp(b,l0,:),2));
        condpow(b) = squeeze(mean(mean(allpowresp(b,l0,respwin),2),3));
        allphases = squeeze(phasemat(b,find(squeeze(phasemat(b,:,:)))));
        condr(cell,b) = circ_r(allphases');
        condcmean(cell,b) = circ_mean(allphases');
    end
end

figure
subplot(2,2,1)
plot(2:5:98,condr(1,:),'b');
hold on
plot(2:5:98,condr(2,:),'r');
legend({'RS','FS'})
title('strength of locking')

subplot(2,2,2)
plot(2:5:98,rad2deg(uncircle(condcmean(1,:))),'b')
hold on
plot(2:5:98,rad2deg(uncircle(condcmean(2,:))),'r')
title('mean phase of spike')

[xc,lags] = xcorr(long(1,:),long(2,:),200);
subplot(2,2,3)
plot(lags,xc)

subplot(2,2,4)
[xc,lags] = xcorr(binit(long(1,:),5),binit(long(2,:),5),200);
plot(lags,xc)
title('binned')

end


function out = uncircle(in)
    in(find(in<0)) = in(find(in<0))+2*pi;
    out = in;
end

  