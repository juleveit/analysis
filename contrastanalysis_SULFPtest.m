function contrastanalysis_SULFPtest

filename = ['C:\Users\Julia\work\data\140520\5contrast.rec'];
path = ['C:\Users\Julia\work\Matlab\sorting\clustering\'];
basename = '0520_5_13-16_waveforms_';
channels = 13:16;

tfiles = dir([path, basename, '*.t']);

%read all t-files and put to struct clust
for i = 1:length(tfiles)
    fid = fopen([path tfiles(i).name],'rb','b');
    clust(i).stamps = fread(fid,inf,'double');
    fclose(fid);
end

 prestim = 300;
 poststim = 300;
 respwin = 500:1500; % after stimulus onset

for cell = 1:length(clust)
   
    data = readTrodesFileDigitalChannels(filename);
    % get the timestamps
    stimons = get_timestamps(data.channelData(1).data);
    
    chans = readTrodesFileChannels(filename,channels);
    lfp = chans.channelData(:,4)-chans.channelData(:,2);
    lfp = eegfilt(lfp',chans.samplingRate,0,100);
    lfp = resample(lfp,1,30);
    lplfp = eegfilt(lfp,1000,0,55);
    hplfp = eegfilt(lfp,1000,65,150);
     
    msStimes = round(clust(cell).stamps);
    msstamps = round(stimons/30);

    % get the stimulus information
    load('C:\Users\Julia\work\data\140520\140520_block5');
    
    chan = zeros(1,round(size(data.timestamps,1)/30));
    chan(msStimes) = 1;
    
    trialdur = result.stimduration*1000;
    respwin = respwin+prestim;

    for i = 1:length(msstamps)
        resp(i,:) = chan(msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim); 
        lfpresp(i,:) = lfp(msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim);
        lpresp(i,:) = lplfp(msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim);
        hpresp(i,:) = hplfp(msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim);
    end
    
    ta = linspace(-prestim,trialdur+poststim,size(resp,2));
    
    lightresp = resp(find(result.light),:);
    nolightresp = resp(find(result.light == 0),:);
    
    binwidth = 10;
    [binnedlight,bta] = binit(mean(lightresp),binwidth);
    [binnednolight,bta] = binit(mean(nolightresp),binwidth);
    
    figure
    ta = bta-prestim;
    plot(ta,binnednolight.*(1000/binwidth),'k','linewidth',2);
    hold on
    plot(ta,binnedlight.*(1000/binwidth),'r','linewidth',2);
    mx = max(binnedlight.*(1000/binwidth));
    axis([-prestim,trialdur+poststim,0,mx]);
    line([0,0],[0,mx],'color','k','linewidth',2);
    line([2000,2000],[0,mx],'color','k','linewidth',2);
    line([500,500],[0,mx],'color','b','linewidth',2)
    line([1500,1500],[0,mx],'color','b','linewidth',2);
    legend({'Light OFF','Light ON'})
    xlabel('time [ms]')
    ylabel('firingrate [Hz]')

    frs = sum(resp(:,respwin),2)./(length(respwin)/1000);
    bl = sum(resp(:,1:prestim),2)./(prestim/1000);
    
    nolbl = mean(bl(find(~result.light)));
    nolblerr = std(bl(find(~result.light)))./(sqrt(length(find(~result.light))));
    lbl = mean(bl(find(result.light)));
    lblerr = std(bl(find(result.light)))./(sqrt(length(find(result.light))));
    
    oris = unique(result.gratingInfo.Orientation);
    clevels = unique(result.gratingInfo.Contrast);
    for l = 1:2
        for ori = 1:length(oris)
            for co = 1:length(clevels)
                thisinds = find(result.gratingInfo.Orientation == oris(ori) &...
                    result.gratingInfo.Contrast == clevels(co) & ...
                    result.light == l-1);
                condresp(l,ori,co,:) = mean(resp(thisinds,:),1);
                condfr(l,ori,co) = mean(frs(thisinds));        
                conderr(l,ori,co) =std(frs(thisinds))./sqrt(length(thisinds));
                condlfp(l,ori,co,:) = mean(lfpresp(thisinds,:),1);
                condhp(l,ori,co,:) = mean(hpresp(thisinds,:),1);
                condlp(l,ori,co,:) = mean(lpresp(thisinds,:),1);
                
            end
        end
    end
     
    figure
    errorbar([0,clevels],[lbl, squeeze(mean(condfr(2,:,:),2))'],...
        [lblerr, squeeze(mean(conderr(2,:,:),2))'],'ro-','markersize',8,'linewidth',2)
    hold on
    errorbar([0,clevels],[nolbl, squeeze(mean(condfr(1,:,:),2))'],...
        [nolblerr, squeeze(mean(conderr(1,:,:),2))'],'ko-','markersize',8,'linewidth',2)
    xlabel('shown contrast')
    ylabel('Firing rate [Hz]')
    legend({'Light ON','Light OFF'})    
    set(gca,'xtick',[0,clevels])
    title('contrast response all orientations');

    fc = find(clevels == 1);
    figure
    errorbar(oris,squeeze(condfr(2,:,fc)),squeeze(conderr(2,:,fc)),'ro-','markersize',8,'linewidth',2)
    hold on
    errorbar(oris,squeeze(condfr(1,:,fc)),squeeze(conderr(1,:,fc)),'ko-','markersize',8,'linewidth',2)
    xlabel('shown orientation')
    ylabel('Firing rate [Hz]')
    set(gca,'xtick',oris)
    legend({'Light ON','Light OFF'})
    
    oneorifr = mean(reshape(condfr(:,:,fc),2,4,2),3);
    prefori = find(oneorifr(1,:) == max(oneorifr(1,:)),1);
    ortho = mod(prefori+2,length(oris));
    
    figure
    errorbar([0,clevels],[lbl,squeeze(mean(condfr(2,[prefori,prefori+(length(oris)/2)],:),2))'],...
        [lblerr, squeeze(mean(conderr(2,[prefori,prefori+(length(oris)/2)],:),2))'],'ro-','markersize',8,'linewidth',2)
    hold on
    errorbar([0,clevels],[nolbl,squeeze(mean(condfr(1,[prefori,prefori+(length(oris)/2)],:),2))'],...
        [nolblerr, squeeze(mean(conderr(1,[prefori,prefori+(length(oris)/2)],:),2))'],'ko-','markersize',8,'linewidth',2)
    errorbar([0,clevels],[lbl,squeeze(mean(condfr(2,[ortho,ortho+(length(oris)/2)],:),2))'],...
        [lblerr, squeeze(mean(conderr(2,[ortho,ortho+(length(oris)/2)],:),2))'],'ro--','markersize',8,'linewidth',1)
    hold on
    errorbar([0,clevels],[nolbl,squeeze(mean(condfr(1,[ortho,ortho+(length(oris)/2)],:),2))'],...
        [nolblerr, squeeze(mean(conderr(1,[ortho,ortho+(length(oris)/2)],:),2))'],'ko--','markersize',8,'linewidth',1)
    xlabel('shown contrast')
    ylabel('Firing rate [Hz]')
    legend({'Light ON preferred','Light OFF preferred', 'Light ON orthogonal', 'Light OFF orthogonal'})    
    set(gca,'xtick',[0,clevels])
    title('contrast response preferred orientations');
    
    disp('');
    
end