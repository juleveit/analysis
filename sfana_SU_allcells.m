function sfana_SU_allcells

animalid = '160129';
block = 6;
lcol = 'r'; %lasercolor

supath = ['C:\Users\Julia\work\data\' animalid '\singleunits\'];
basename = [animalid '_block' int2str(block) '_tet'];

files = dir([supath, basename, '*.mat']);

prestim = 300;
poststim = 300;
respwin = 501:1500; % after stimulus onset
respwin = respwin+prestim;

% chronux parameters
params.tapers = [5,9]; params.Fs = 1000; params.err = [2, 0.05]; params.trialave = 1;

% get it to cm/s speed
screendeg = 62;
screencm = 28;
degpercm = screendeg/screencm;

% for fi = 1:length(vismods) 
for fi = 1:length(files)
   
    cell = fi;
%     cell = vismods(fi);

    load([supath, files(cell).name]);
%     disp(['now analyzing file: ' files(cell).name]);
    
   i = strfind(files(cell).name, 'tet');
    if strcmp(files(cell).name(i+4),'_')
        tetno = strread(files(cell).name(i+3)); % single character number
    else        
        tetno = strread(files(cell).name(i+3:i+4)); % number >10
    end
    if tetno>8
        v1(cell) = logical(0); v2(cell) = logical(1); cellstr = 'V2';
    else
        v1(cell) = logical(1); v2(cell) = logical(0); cellstr = 'V1';
    end
    
    wvchan = find(var(result.waveforms) == max(var(result.waveforms)));
    spike = result.waveforms(:,wvchan);
    interpspike = spline(1:32,spike,1:.1:32);
    [adiff(cell),swidth(cell),ptr(cell),eslope(cell)] = spikequant(interpspike);
    depth(fi) = result.depth;
    
    lfp = result.lfp(:,wvchan);
    
    
    msStimes = round(result.spikes);
    if ~isempty(msStimes) & msStimes(1) == 0, msStimes(1) = 1; end
    
    chan = zeros(1,length(result.lfp));
    chan(msStimes) = 1;
    
    trialdur = result.stimduration*1000;
    msstamps = result.msstamps;
        
    if length(msstamps)~=length(result.light)
        %             disp('');
        %             msstamps(16) = []; % for 150414 block 10
                    msstamps(239) = []; % for 160129 block 6
                    result.msstamps = msstamps;
                    save([supath, files(cell).name],'result');
%         pause;
    end
    
    
    for i = 1:length(msstamps)
        resp(i,:) = chan(msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim); 
        lfpresp(i,:) = result.lfp(msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim, wvchan);
        allS(i,:)=mtspectrumc(squeeze(lfpresp(i,1001:1800))',params);
    end
    
    ta = linspace(-prestim,trialdur+poststim,size(resp,2));
    
    %%%%%%%%%%%%%%%%%%%%%%%%
    %find gamma peaks for this animal
    beta = [15,40];
    large = find(result.gratingInfo.spFreq == 0.04 & result.gratingInfo.tFreq == 2 & result.light == 0);
    
    for i = 1:length(large)
        [pl(i,:),f] = mtspectrumc(lfp(result.msstamps(large(i))+700:result.msstamps(large(i))+1500),params);
    end
    b1 = find(f>beta(1),1); b2 = find(f>beta(2),1);
    bsig = nanmean(pl(:,b1:b2));
    if isempty(find(diff(bsig)>0)) % there is no clear beta peak
        bpi = round((b1+b2)/2);
    else
        peaks = find(diff(bsig)>0)+1;
        pvs = bsig(peaks);
        bpi = peaks(pvs == max(pvs));
        bpi = bpi+b1-1;
    end
    gammaind(cell) = bpi;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    lightresp = resp(find(result.light),:);
    nolightresp = resp(find(result.light == 0),:);
    
    lightlfpresp = lfpresp(find(result.light),:);
    nolightlfpresp = lfpresp(find(~result.light),:);
    
    frs = sum(resp(:,respwin),2)./(length(respwin)/1000);
    bl = sum(resp(:,1:prestim),2)./(prestim/1000);
    
    lfr(fi) = mean(frs(find(result.light)));
    nlfr(fi) = mean(frs(find(result.light == 0)));
    
     %determine if cell is visually modulated
    blfr = sum(resp(:,1:prestim),2);
    vrfr = sum(resp(:,prestim+40:2*prestim+40),2);
    vismod(fi) = ttest2(blfr,vrfr);
    
    %determine if cell is modulated by light
    lightmod(fi) = ttest2(frs(find(result.light)),frs(find(result.light == 0)));
    
    binwidth = 20;
    [binnedlight,bta] = binit(mean(lightresp),binwidth);
    binnedlight = binnedlight.*(1000/binwidth);
    [binnednolight,bta] = binit(mean(nolightresp),binwidth);
    binnednolight = binnednolight.*(1000/binwidth);
    
    printname = files(cell).name;
    printname(find(printname=='_')) = ' ';
    
    spFreqs = unique(result.gratingInfo.spFreq);   %delete control condition
    tFreqs = unique(result.gratingInfo.tFreq);
    spFreqscm = 1./(spFreqs.*degpercm); % in cm/cycle (* .5 = stripe thickness)
    for i = 1:length(spFreqs)
        tFreqscms(:,i) = spFreqscm(i).*tFreqs;
    end
    for l = 1:2
        for tf = 1:length(tFreqs)
            for sf = 1:length(spFreqs)
                thisinds = find(result.gratingInfo.tFreq == tFreqs(tf) &...
                    result.gratingInfo.spFreq == spFreqs(sf) & ...
                    result.light == l-1);
                condresp(l,tf,sf,:) = mean(resp(thisinds,:),1);
                condfr(l,tf,sf) = mean(frs(thisinds));%-mean(bl);        
                conderr(l,tf,sf) =std(frs(thisinds))./sqrt(length(thisinds));                
                condresperr(l,tf,sf,:) = nanstd(resp(thisinds,:),1,1)./sqrt(length(thisinds));
                [S,chf,Serr]=mtspectrumc(squeeze(lfpresp(thisinds,1001:1800))',params);
                condS(cell,l,tf,sf,:) = S(1:150); condSerr(cell,l,tf,sf,:,:) = Serr(:,1:150);
                condallS(cell,l,tf,sf,:) = nanmean(allS(thisinds,:),1);
                condallSerr(cell,l,tf,sf,:) = nanstd(allS(thisinds,:),1,1)./sqrt(length(thisinds));

                if ~isnan(condresp(l,tf,sf,:))
                    [bincondresp(l,tf,sf,:),bta] = binit(condresp(l,tf,sf,:),binwidth);
                else
                    bincondresp(l,tf,sf,:) = binit(condresp(l,tf,sf,:),binwidth);
                end
                binconderr(l,tf,sf,:) = binit(condresperr(l,tf,sf,:),binwidth);
                
                condlfpresp(l,tf,sf,:) = mean(lfpresp(thisinds,:),1);
            end
        end
    end
    
    bincondresp = bincondresp.*(1000/binwidth);
    bta = bta-prestim;
%     psthcondplot(bincondresp,binconderr,bta);
  
%     figure
%     errorbar(sizes,squeeze(nanmean(condfr(2,:,:),2)),squeeze(nanmean(conderr(2,:,:),2)),'o-','color',lcol,'markersize',8,'linewidth',2)
%     hold on
%     errorbar(sizes,squeeze(nanmean(condfr(1,:,:),2)),squeeze(nanmean(conderr(1,:,:),2)),'ko-','markersize',8,'linewidth',2)
%     xlabel('shown patch size [vd]')
%     ylabel('Firing rate [Hz]')
%     legend({'Light ON','Light OFF'})    
%     set(gca,'xtick',sizes)  
%     title(['cell ' int2str(cell) ' average all orientations'])
    
    prefspfreq = find(mean(condfr(1,:,:),2) == max(mean(condfr(1,:,:),2)),1);
    preftfreq = find(mean(condfr(1,:,:),3) == max(mean(condfr(1,:,:),3)),1);
    
    psf(fi) = spFreqs(prefspfreq);
    ptf(fi) = tFreqs(preftfreq);
    
    figure
    subplot(2,2,1)
    ta = bta;
    plot(ta,binnednolight,'k','linewidth',2);
    hold on
    plot(ta,binnedlight,lcol,'linewidth',2);
    mx = max([max(binnednolight),max(binnedlight),0.1]);
    axis([-prestim,trialdur+poststim,0,mx]);
    line([0,0],[0,mx],'color','k','linewidth',2);
    line([2000,2000],[0,mx],'color','k','linewidth',2);
%     line([500,500],[0,mx],'color','b','linewidth',2)
%     line([1500,1500],[0,mx],'color','b','linewidth',2);
%     legend({'Light OFF','Light ON'})
    xlabel('time [ms]')
    ylabel('firing rate [Hz]')
    title(['cell ' int2str(cell) ' depth: ' int2str(result.depth), 'cell ' printname ])
    
    subplot(2,2,2)
    errorbar(tFreqs,squeeze(condfr(2,:,prefspfreq)),squeeze(conderr(2,:,prefspfreq)),'o-','color',lcol,'markersize',8,'linewidth',2)
    hold on
    errorbar(tFreqs,squeeze(condfr(1,:,prefspfreq)),squeeze(conderr(1,:,prefspfreq)),'ko-','markersize',8,'linewidth',2)
    xlabel('shown temporal frequency')
    ylabel('Firing rate [Hz]')
    set(gca,'xtick',tFreqs)
%     legend({'Light ON','Light OFF'})
    title([' TF Tuning at preferred SF '])
    
    subplot(2,2,3)
    errorbar(spFreqs,squeeze(condfr(2,preftfreq,:)),...
        squeeze(conderr(2,preftfreq,:)),'o-','color',lcol,'markersize',8,'linewidth',2);
    hold on
    errorbar(spFreqs,squeeze(condfr(1,preftfreq,:)),...
        squeeze(conderr(1,preftfreq,:)),'ko-','markersize',8,'linewidth',2);
    xlabel('shown spfreq')
    ylabel('Firing rate [Hz]')
%     legend({'Light ON','Light OFF'})    
    set(gca,'xtick',spFreqs)  
    title(['SF Tuning at preferred TF'])
    
    subplot(2,4,7)
    plot(spike)
    axis([0,40,-100,100])
    legend(['width: ' int2str(swidth(fi)) ' adiff: ' num2str(adiff(fi))])

%     printpath = ['C:\Users\Julia\work\data\V2\SFTuning\' animalid '\'];
%     figSize = [30 21];
%     set(gcf,'paperunits','centimeters','papersize',figSize,'paperposition',[0 0 figSize])
%     if cell<10, printi = ['0', int2str(cell)]; else printi = int2str(cell); end
%     print([printpath ,  printi '__' files(fi).name '.pdf'],'-dpdf')
    
    disp('');
    
end

%spike classification
kmeansind = kmeans([eslope',ptr',swidth',adiff'],2);
if mean(swidth(find(kmeansind==1)))<mean(swidth(find(kmeansind==2)))  %1 is FS
    pfs = find(kmeansind==1); prs = find(kmeansind==2); pfsv = kmeansind==1; prsv = kmeansind==2;
else
    pfs = find(kmeansind==2); prs = find(kmeansind==1); pfsv = kmeansind==2; prsv = kmeansind==1;
end

figure
plot(swidth(prs),adiff(prs),'b.')
xlabel('spike width')
ylabel('amplitude diff')
hold on
plot(swidth(pfs),adiff(pfs),'r.')
axis([5.5,22.5,-.9,.7])

depth(v1) = depth(v1).*cosd(22);
depth(v2) = depth(v2).*cosd(45); % guess until histology

figure
plot(ptf,depth,'k.')
axis([0,5,250,1000])
axis ij
xlabel('preferred temporal frequency')
ylabel('depth')

figure
plot(psf,depth,'k.')
axis([0,.2,250,1000])
axis ij
xlabel('preferred spatial frequency')
ylabel('depth')

disp('');
