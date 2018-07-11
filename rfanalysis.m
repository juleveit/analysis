function rfanalysis

normalize_plot = 0;

filename = 'C:\Users\Julia\work\data\140513\sn3.rec';
channels = [1:16];
refchan = 6;
for channo = channels
    chans = readTrodesFileChannels(filename,channo);
    refchannel = readTrodesFileChannels(filename,refchan);
    chans.channelData = chans.channelData-refchannel.channelData;
    data = readTrodesFileDigitalChannels(filename);

    % get the spikes
    [matclustdata, waves, spiketimes, filtData] = findSpikes(chans,-4);

    % get the timestamps
    stimons = get_timestamps(data.channelData(1).data);
    
    % get wheelspeed
    [speed, ta] = get_wheelspeed(data.channelData(3).data);    
    
    % put to 1ms resolution
    filtData = resample(filtData,1,30);
    msStimes = round(spiketimes/30);
    msstamps = round(stimons/30);

    % get the stimulus information
    load('C:\Users\Julia\work\data\140513\140513_block3');
    normimg = result.stimulus;
    normimg(find(normimg == 0)) = -1;
    normimg(find(normimg == 128)) = 0;
    normimg(find(normimg == 255)) = 1;
    trialdur = (1000/result.frameRate)*result.FramesperStim;

    chan = zeros(1,size(filtData,1));
    chan(msStimes) = 1;

    sr = 1000;
    %window for calculating spikerates
    srwinoffs = 20;
    srwinsize = 79;
    trialstart = -0;
    trialend = 2*round(trialdur);
    triallen = trialend-trialstart;
    trialnsamps = round((trialdur*sr)/1000);

    for trial = 1:length(msstamps)
        resp(trial,:) = chan(msstamps(trial):msstamps(trial)+triallen-1);
    end

    nspikes = sum(resp(:,srwinoffs:srwinoffs+srwinsize),2).*(1000/srwinsize+1);

    nstims = size(result.stimulus,3);
    mpT1 = zeros(nstims,result.repetitions);
    for rep = 1:result.repetitions
        for i = 1:nstims
            mp(result.stimulusIndex((rep-1)*nstims+i),rep) = nspikes((rep-1)*nstims+i);
        end
    end

    weightedimg = zeros(size(result.stimulus));
    for rep = 1:result.repetitions
        for i = 1:size(result.stimulus,1)
            for j = 1:size(result.stimulus,2)
                weightedimg(i,j,:) = squeeze(weightedimg(i,j,:)) + squeeze(normimg(i,j,:)).*mp(:,rep);
            end
        end
    end
    weightedimg = weightedimg./result.repetitions;
    pos = weightedimg; neg = weightedimg;
    pos(find(weightedimg<0)) = 0;
    neg(find(weightedimg>0)) = 0;

    trunc = result.sizePixel-1; % truncate by sizePixel-1 on the borders (undersampled)
    a = trunc+1;
    b = size(pos,1)-trunc;
    degperElem = (result.shownSize*result.sizePixel)/size(result.stimulus,1);
    nelem = result.stimulusElements;
    degImg = result.stimulusSize;
    ax = linspace(1,degImg,nelem);
    degperGridPos = degImg/nelem;
    off = abs(squeeze(sum(neg(a:b,a:b,:),3)));
    on = abs(squeeze(sum(pos(a:b,a:b,:),3)));
    off = off./(result.sizePixel)^2; %not yet compensated for!!
    on = on./(result.sizePixel)^2;
    if normalize_plot
        off = off - min(min(off));
        off = off./max(max(off));
        on = on-min(min(on));
        on = on./max(max(on));
    end

    a = ax*degperGridPos;
    w = result.stimulusSize/2;
    xax = linspace(result.position(1)-w,result.position(1)+w,result.stimulusElements);
    yax = linspace(result.position(2)-w,result.position(2)+w,result.stimulusElements);
    
    figure
    subplot(1,2,1)
    imagesc(xax,yax,on);
    axis square
    colorbar
    title(['contact ' int2str(channo) '  ON field'])

    subplot(1,2,2)
    imagesc(xax,yax,off);
    axis('square')
    colorbar;
    title('OFF field')
 
    disp('');

    onfields(:,:,channo) = on;
    offfields(:,:,channo) = off;
    respmats(channo,:,:) = resp;
end

disp('');
