function sizeanalysis_SU

filename = ['C:\Users\Julia\work\data\140530\3size.rec'];
path = ['C:\Users\Julia\work\data\140530\sorting\'];
basename = '3size_1-4_waveforms_';
lcol = 'g'; %laser color for printing

tfiles = dir([path, basename, '*.t']);
cd(path);
S = LoadSpikes({tfiles.name})

for i = 1:length(S)
    clust(i).stamps = S{i}.T;
end

% %read all t-files and put to struct clust
% for i = 1:length(tfiles)
%     fid = fopen([path tfiles(i).name],'rb','b');
%     clust(i).stamps = fread(fid,inf,'uint32');
%     fclose(fid);
% end

for cell = 1:length(clust)
   
    data = readTrodesFileDigitalChannels(filename);
    % get the timestamps
    stimons = get_timestamps(data.channelData(1).data);

    msStimes = round(clust(cell).stamps*1000);
    msstamps = round(stimons/30);

    % get the stimulus information
    load('C:\Users\Julia\work\data\140530\140530_block3');
    
    chan = zeros(1,round(size(data.timestamps,1)/30));
    chan(msStimes) = 1;
    
    prestim = 300;
    poststim = 300;
    trialdur = result.stimduration*1000;

    for i = 1:length(msstamps)
        resp(i,:) = chan(msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim); 
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
    plot(ta,binnedlight.*(1000/binwidth),lcol,'linewidth',2);
    mx = max(binnedlight.*(1000/binwidth));
    axis([-prestim,trialdur+poststim,0,mx]);
    line([0,0],[0,mx],'color','k','linewidth',2);
    line([2000,2000],[0,mx],'color','k','linewidth',2);
    line([500,500],[0,mx],'color','b','linewidth',2)
    line([1500,1500],[0,mx],'color','b','linewidth',2);
    legend({'Light OFF','Light ON'})
    xlabel('time [ms]')
    ylabel('firing rate [Hz]')
    title(['cell ' int2str(cell)])

    frs = sum(resp(:,prestim:prestim+trialdur-1),2)./(trialdur/1000);
    
    sizes = unique(result.gratingInfo.size);    
    oris = unique(result.gratingInfo.Orientation);
    for l = 1:2
        for ori = 1:length(oris)
            for sz = 1:length(sizes)
                thisinds = find(result.gratingInfo.Orientation == oris(ori) &...
                    result.gratingInfo.size == sizes(sz) & ...
                    result.light == l-1);
                condresp(l,ori,sz,:) = mean(resp(thisinds,:),1);
                condfr(l,ori,sz) = mean(frs(thisinds));        
                conderr(l,ori,sz) =std(frs(thisinds))./sqrt(length(thisinds));
            end
        end
    end
    
    figure
    errorbar(sizes,squeeze(mean(condfr(2,:,:),2)),squeeze(mean(conderr(2,:,:),2)),'o-','color',lcol,'markersize',8,'linewidth',2)
    hold on
    errorbar(sizes,squeeze(mean(condfr(1,:,:),2)),squeeze(mean(conderr(1,:,:),2)),'ko-','markersize',8,'linewidth',2)
    xlabel('shown patch size [vd]')
    ylabel('Firing rate [Hz]')
    legend({'Light ON','Light OFF'})    
    set(gca,'xtick',sizes)  
    title(['cell ' int2str(cell) ' average all orientations'])
    
    prefsize = find(mean(condfr(1,:,:),2) == max(mean(condfr(1,:,:),2)));
    
    [nlprefratio, nlprefori, nlmeanori, nlosi, nlmeandir, nldsi] = getOSI(squeeze(condfr(1,:,prefsize)),oris);
    [lprefratio, lprefori, lmeanori, losi, lmeandir, ldsi] = getOSI(squeeze(condfr(2,:,prefsize)),oris);
    
    figure
    errorbar(oris,squeeze(condfr(2,:,prefsize)),squeeze(conderr(2,:,prefsize)),'o-','color',lcol,'markersize',8,'linewidth',2)
    hold on
    errorbar(oris,squeeze(condfr(1,:,prefsize)),squeeze(conderr(1,:,prefsize)),'ko-','markersize',8,'linewidth',2)
    xlabel('shown orientation')
    ylabel('Firing rate [Hz]')
    set(gca,'xtick',oris)
    legend({'Light ON','Light OFF'})
    title(['cell ' int2str(cell)])
    
    disp('');
    
end