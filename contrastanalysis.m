function contrastanalysis

filename = 'C:\Users\Julia\work\data\141126b\10fcgrating.rec';
channels = [10];
refchan = 8;
for channo = channels
    chans = readTrodesFileChannels(filename,channo);
    refchannel = readTrodesFileChannels(filename,refchan);
    chans.channelData = chans.channelData-refchannel.channelData;
    
    data = readTrodesFileDigitalChannels(filename);

    % get the spikes
    [matclustdata, waves, spiketimes, filtData] = findSpikes(chans,-4);

    % get the timestamps
    stimons = get_timestamps(data.channelData(1).data);

    % put to 1ms resolution
    filtData = resample(filtData,1,30);
    msStimes = round(spiketimes/30);
    msstamps = round(stimons/30);

    % get the stimulus information
    load('C:\Users\Julia\work\data\141126b\141126b_block10');
    
    chan = zeros(1,size(filtData,1));
    chan(msStimes) = 1;
    
    prestim = 300;
    trialdur = result.stimduration*1000;

    for i = 1:length(msstamps)
        resp(i,:) = chan(msstamps(i) - prestim+1:msstamps(i) + trialdur); 
    end
    
    ta = linspace(-prestim,trialdur,size(resp,2));
    
    lightresp = resp(find(result.light),:);
    nolightresp = resp(find(result.light == 0),:);
    
    binwidth = 20;
    [binnedlight,bta] = binit(mean(lightresp),binwidth);
    [binnednolight,bta] = binit(mean(nolightresp),binwidth);
    
    figure
    ta = bta-prestim;
    plot(ta,binnednolight./binwidth,'k','linewidth',2);
    hold on
    plot(ta,binnedlight./binwidth,'b','linewidth',2);
    mx = max(binnedlight/binwidth);
    axis([-prestim,trialdur,0,mx]);
    line([0,0],[0,mx],'color','k','linewidth',2);
    line([500,500],[0,mx],'color','r','linewidth',2)
    line([1500,1500],[0,mx],'color','r','linewidth',2);
    legend({'LED OFF','LED ON'})
    xlabel('time [ms]')
    ylabel('firing prob.')

    frs = sum(resp(:,prestim:prestim+trialdur-1),2);
    
    clevels = unique(result.gratingInfo.Contrast);
    for co = 1:length(clevels)
        thiscos = result.gratingInfo.Contrast == clevels(co);
        condresp(1,co,:) = mean(resp(find(thiscos&~result.light),:));
        condresp(2,co,:) = mean(resp(find(thiscos&result.light),:));
        condfr(1,co) = mean(frs(find(thiscos&~result.light)));
        conderr(1,co) =std(frs(find(thiscos&~result.light)))./sqrt(length(frs)/2);
        condfr(2,co) = mean(frs(find(thiscos&result.light)));
        conderr(2,co) =std(frs(find(thiscos&result.light)))./sqrt(length(frs)/2);
    end
    
    figure
    errorbar(clevels,condfr(2,:),conderr(2,:),'bo-','markersize',8,'linewidth',2)
    hold on
    errorbar(clevels,condfr(1,:),conderr(1,:),'ko-','markersize',8,'linewidth',2)
    xlabel('shown patch size [vd]')
    ylabel('Firing rate [Hz]')
    legend({'LED ON','LED OFF'})    
    set(gca,'xtick',clevels)

    
    oris = unique(result.gratingInfo.Orientation);
    for ori = 1:length(oris)
        thisori = result.gratingInfo.Orientation == oris(ori);
        oriresp(1,ori,:) = mean(resp(find(thisori&~result.light),:));
        oriresp(2,ori,:) = mean(resp(find(thisori&result.light),:));
        orifr(1,ori) = mean(frs(find(thisori&~result.light)));
        orierr(1,ori) = std(frs(find(thisori&~result.light)))./sqrt(length(frs)/2);
        orifr(2,ori) = mean(frs(find(thisori&result.light)));
        orierr(2,ori) = std(frs(find(thisori&result.light)))./sqrt(length(frs)/2);
    end
    orifr = sum(oriresp(:,:,prestim:prestim+trialdur-1),3)./result.stimduration;
    
    figure
    errorbar(oris,orifr(2,:),orierr(2,:),'bo-','markersize',8,'linewidth',2)
    hold on
    errorbar(oris,orifr(1,:),orierr(1,:),'ko-','markersize',8,'linewidth',2)
    xlabel('shown orientation')
    ylabel('Firing rate [Hz]')
    set(gca,'xtick',oris)
    legend({'LED ON','LED OFF'})
    
    disp('');
    
end