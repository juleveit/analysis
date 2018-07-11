function MUAana_VEP

animalid = '180619';
block = 2;
basepath = ['C:\Users\Julia\work\data\' animalid '\'];
channels = 1:16;
recdepth = 600;
spacing = 25;
penangle = 25;

recalculate = 0;

prestim = 3000;
poststim = 3000;
trialdur = 5000;
respwin = 1:5000; %501:1500; % after stimulus onset
respwin = respwin+prestim;

depth = recdepth.*(cosd(penangle)*cosd(22));
spacing = spacing.*(cosd(penangle)*cosd(22));
evaldepth = 375;
for i = 1:length(depth)
    for j = 1:length(channels(1):channels(end))
        dm(i,j) = depth(i)-((j-1)*spacing(i)); % depth matrix
    end
    [c,di(i)] = min(abs(dm(i,:)-evaldepth)); % get depth index, least distance to evaldepth
end       

muaresultfile = [basepath, 'block' int2str(block) 'chan_' int2str(di) 'MUAresult.mat'];
if ~exist(muaresultfile) || recalculate
    recfiles = dir([basepath, '*.rec']);
    for i = 1:length(recfiles)
        thisname = recfiles(i).name;
        nu(i) = strread(thisname); % gets the blocknumber
    end
    blind = find(nu == block);

    analogdata = readTrodesFileDigitalChannels([basepath recfiles(blind).name]);
    % get the timestamps
    stimons = get_timestamps(analogdata.channelData(1).data);
    result.msstamps = round(stimons/30); % divide by 30 to go from 30kHz sampling to 1ms resolution
    chans = readTrodesFileContinuous([basepath recfiles(blind).name],[di,ones(1,1)]);

    lfp = eegfilt(chans.channelData',30000,0,200)';
    result.lfp = resample(lfp,1,30);

    depth = dm(di);

    Wp = [ 700 8000] * 2 / chans.samplingRate; % pass band for filtering
    Ws = [ 500 10000] * 2 / chans.samplingRate; % transition zone
    [N,Wn] = buttord( Wp, Ws, 3, 20); % determine filter parameters
    [B,A] = butter(N,Wn); % builds filter
    fltdata{1} = filtfilt( B, A, chans.channelData ); % runs filter

    spikes = ss_default_params(chans.samplingRate);
    spikes = ss_detect(fltdata,spikes);
    msstimes = round(spikes.spiketimes'.*1000);
    result.chan = zeros(size(result.lfp));
    result.chan(msstimes) = 1;

    save(muaresultfile,'result');
else
    load(muaresultfile);
end

wo = 6.5/(1000/2); bw = wo/1.5;
[b,a] = iirnotch(wo,bw);
filtlfp = filtfilt(b,a,result.lfp);
binwidth = 50;
msstamps = result.msstamps;
for i = 1:length(msstamps)-1
    lfpresp(i,:) = result.lfp(msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim);
    filtlfpresp(i,:) = filtlfp(msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim);
    resp(i,:) = result.chan(msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim);
    [binresp(i,:),bta] = binit(squeeze(resp(i,:)),binwidth);

end
ta = -prestim+1:trialdur+poststim;
binresp = binresp.*(1000/binwidth);
bta = bta-prestim;

figure
plot(result.lfp);
hold on
for i = 1:length(msstimes)
    line([msstimes(i),msstimes(i)],[0,100],'color','r')
end
for i = 1:length(msstamps)
    line([msstamps(i),msstamps(i)],[300,400],'color','g')
end

figure
plot(ta,lfpresp');
hold on
plot(ta,mean(lfpresp),'k','linewidth',2)

figure
plot(ta,filtlfpresp');
hold on
plot(ta,mean(filtlfpresp),'k','linewidth',2)

figure
errorbar(ta,mean(lfpresp),std(lfpresp)./sqrt(100));

figure
errorbar(bta,nanmean(binresp),nanstd(binresp)./sqrt(length(msstamps)))

 disp('');