function result = MUAdataprepareNOPTB(basepath,animalid,block,channels)

recfiles = dir([basepath, '*.rec']);
for i = 1:length(recfiles)
    thisname = recfiles(i).name;
    nu(i) = strread(thisname); % gets the blocknumber
end
blno = find(nu == block);
analogdata = readTrodesFileDigitalChannels([basepath recfiles(blno).name]);

% get the timestamps
stimons = get_timestamps(analogdata.channelData(1).data);
result.msstamps = round(stimons/30);

% get licks
licks = get_timestamps(analogdata.channelData(2).data);
result.licks = round(licks/30);

% get rewards
rewards = get_timestamps(analogdata.channelData(4).data);
result.rewards = round(rewards/30);

%get wheelspeed
[speed, ta] = get_wheelspeed(analogdata.channelData(3).data);

x_t = cumsum(speed);
[smooth_win, FWHM] = MakeGaussWindow(round(1000),23.5/2, 1000);
sw_len = length(smooth_win);

x_t(end+1:end+sw_len) = x_t(end); %pad with last value for length of smoothing kernel
d_smooth_win = [0;diff(smooth_win)]/(1/1000);
dx_dt = conv(x_t,d_smooth_win,'same');
dx_dt(end-sw_len+1:end) = []; %remove values produced by convolving kernel with padded values
mouserad = 6; %cm
scalefact = 2*pi*mouserad/360;
result.runspeed = dx_dt.*scalefact;

% chans = readTrodesFileChannels([basepath recfiles(blno).name],channels);
chans = readTrodesFileContinuous([basepath recfiles(blno).name],[channels',ones(length(channels),1)]);

Wp = [ 700 8000] * 2 / chans.samplingRate; % pass band for filtering
Ws = [ 500 10000] * 2 / chans.samplingRate; % transition zone
[N,Wn] = buttord( Wp, Ws, 3, 20); % determine filter parameters
[B,A] = butter(N,Wn); % builds filter
for i = 1:length(channels)

    fltdata{1} = filtfilt( B, A, chans.channelData(:,i) ); % runs filter
    muac = fltdata{1};

    spikes = ss_default_params(chans.samplingRate);
    spikes = ss_detect(fltdata,spikes);

    thisStimesms = spikes.spiketimes'.*1000;
    result.msStimes{i} = thisStimesms;
    
    lfp = eegfilt(chans.channelData(:,i)',30000,0,200)';
    result.lfp(i,:) = resample(lfp,1,30);
    result.muac(i,:) = resample(muac,1,30);
end
