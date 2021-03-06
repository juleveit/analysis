function contrastanalysis_SU_nolight

animalid = '150225';
block = 19;
lcol = 'r';

onlymod = 0;
sfc = 1;

basepath = ['C:\Users\Julia\work\data\' animalid '\'];
supath = [basepath 'singleunits\'];
basename = [animalid '_block' int2str(block) '_tet'];

files = dir([supath, basename, '*.mat']);

prestim = 300;
poststim = 300;
respwin = 501:1500; % after stimulus onset
respwin = respwin+prestim;
freqbinwidth = 5;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %quick hack only for 06/03 because timestamps still missing
% recording = '8contrast.rec';
% filedepth = 950;
% filename = ['C:\Users\Julia\work\data\' animalid '\' recording];
% 
% analogdata = readTrodesFileDigitalChannels(filename);
% % get the timestamps
% stimons = get_timestamps(analogdata.channelData(1).data);
% msstamps = round(stimons/30);
% 
% %get wheelspeed
% [speed, ta] = get_wheelspeed(analogdata.channelData(3).data);
% speed = speed;
% speedta = ta;
% contactspacing = 25; % in microns
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% vlmods = [18 19 24 25 26 28 29 30 31 32 33 34 35 36 37 38 39 41 42];

cell = 1;
for fi = 1:length(files)
    
    if strfind(files(fi).name, 'MU')
        continue;
    end
    
    load([supath, files(fi).name]);
    
    % calc spiking to see if includable
    msStimes = round(result.spikes);
    if ~isempty(msStimes) & msStimes(1) == 0, msStimes(1) = 1; end
    
    chan = zeros(1,length(result.lfp));
    chan(msStimes) = 1;
    
    wvchan = find(var(result.waveforms) == max(var(result.waveforms)));
    
    % LFP GAMMA
%     lfp = lfp-mean(lfp(1,1:beg),2);
%     stdbl = std(lfp(:,1:beg)');
%     lfp = lfp./stdbl;
    
    sr = 1000;
    lfp = result.lfp(:,wvchan)';
    nfft = 2^nextpow2(length(lfp));
    fax = sr/2*linspace(0,1,nfft/2+1);
    y = fft(lfp,nfft);
    lfpspectrum = abs(y(1:nfft/2+1));

%     plot(fax,lfpspectrum)
    gamma = eegfilt(lfp,sr,30,90);
    h = hilbert(gamma); gpow = abs(h); gphas = angle(h);
    % %
    if sfc
        for i = 1:100/freqbinwidth
            filtmat(i,:) = eegfilt(lfp,sr,(i-1)*freqbinwidth+1,i*freqbinwidth);
            h = hilbert(filtmat(i,:));
            powmat(i,:) = abs(h); phasmat(i,:) = angle(h);
        end
    end
        

    trialdur = result.stimduration*1000;
    msstamps = result.msstamps;
    if length(msstamps)~=length(result.light)
%         msstamps([54,201,239,316]) = []; % for 140807 block 7
%         msstamps(200) = []; % for 140815 block 7
%         msstamps(123) = []; % for 140815 block 7
%         msstamps(385) = []; % for 140703 block 5
%         msstamps(37) = []; % for 150113 block 7
%         msstamps(39) = []; % for 150113 block 11
%         msstamps(37) = []; % for 150119 block 9
%         msstamps(25) = []; % for 150119 block 18
%         newmsstamps = zeros(80,1);  % for 150129 block 2 missing stamps
%         newmsstamps(1:34) = msstamps(1:34);
%         newmsstamps(37:80) = msstamps(35:78);
%         newmsstamps(35) = newmsstamps(34)+mean(diff(msstamps(1:34)));
%         newmsstamps(36) = newmsstamps(35)+mean(diff(msstamps(1:34)));
%         msstamps = newmsstamps;
%         result.msstamps = msstamps;
%         save([supath, files(fi).name],'result');
        pause;
    end
%     trialnfft = 2^nextpow2(800);
%     trialfax = sr/2*linspace(0,1,trialnfft/2+1);
    for i = 1:length(msstamps)
        resp(i,:) = chan(msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim); 
        lfpresp(i,:) = lfp(msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim);
%         y = fft(lfpresp(i,1001:1800),trialnfft);
%         lfpspect(i,:) = abs(y(1:trialnfft/2+1));
        [lfpspect(i,:),trialfax] = pmtm(lfpresp(i,1001:1800),3,[],sr);
        gammaresp(i,:) = gpow(msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim);
        gphaseresp(i,:) = gphas(msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim);
        
        if sfc
            for j = 1:size(phasmat,1)
                allphaseresp(j,i,:) = phasmat(j, msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim);
                allpowresp(j,i,:) = powmat(j, msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim);
            end
        end
        
        speed(i,:) = result.runspeed(msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim);     
    end
    
    
    % figure out sufficiently high and nonvariable runspeed trials
    meanspeed = mean(speed(:,respwin),2);
    stdspeed = std(speed(:,respwin),1,2);
    notstill = find(meanspeed>1);
    okspeed = find(meanspeed>( mean(meanspeed(notstill))-(1.5*std(meanspeed(notstill))) ) );
    okvar = find(stdspeed<( mean(stdspeed(notstill))+(1.5*std(stdspeed(notstill)))) & stdspeed>.5);
    oktrials = intersect(okspeed,okvar);
    nonoktrials = 1:size(resp,1); nonoktrials(oktrials) = [];
    stilltrials = 1:size(resp,1); stilltrials(notstill) = [];

    
    frs = sum(resp(:,respwin),2)./(length(respwin)/1000);
    bl = sum(resp(:,1:prestim),2)./(prestim/1000);
    
    %determine if cell is visually modulated
    blfr = sum(resp(:,1:prestim),2);
    vrfr = sum(resp(:,prestim+40:2*prestim+40-1),2);
    vismod(fi) = ttest2(blfr,vrfr);
    visdriven(fi) = mean(vrfr)>=mean(blfr)+2; % average firing rate is increase at least 2Hz above baseline
    
    if onlymod & ~vismod(fi) %~(lightmod(fi) & vismod(fi))
        continue;
    end
    
    cellname{cell} = files(cell).name;
    
    i = strfind(files(cell).name, 'tet');
    if strcmp(files(cell).name(i+4),'_')
        tetno = strread(files(cell).name(i+3)); % single character number
    else        
        tetno = strread(files(cell).name(i+3:i+4)); % number >10
    end
    tetnos(cell) = tetno;
    if tetno>8
        v1(cell) = logical(0); v2(cell) = logical(1); cellstr = 'V2';
    else
        v1(cell) = logical(1); v2(cell) = logical(0); cellstr = 'V1';
    end
    
    spike = result.waveforms(:,wvchan); 
    interpspike = spline(1:32,spike,1:.1:32);
    [adiff(cell),swidth(cell)] = spikequant(interpspike);
    
    % phases
    tmp = zeros(size(gphaseresp));
    tmp(find(resp)) = gphaseresp(find(resp));
    phases{cell} = tmp(find(tmp));
      
    phaser(cell) = circ_r(phases{cell});
    cmean(cell) = circ_mean(phases{cell});
    
    if sfc
        tmpi = zeros(size(allphaseresp));
        for i = 1:size(allphaseresp,1)
            tmp = zeros(size(gphaseresp));
            tmp(find(resp)) = allphaseresp(i,find(resp));
            tmpi(i,:,:) = tmp;
        end
        allphasemat = tmpi;
        for i = 1:size(allphaseresp,1)
            allphases{i} = allphasemat(i,find(squeeze(allphasemat(i,:,:)))); 
            allphaser(cell,i) = circ_r(allphases{i}');
            allcmean(cell,i) = circ_mean(allphases{i}');
        end
    end
    
    depth(cell) = result.depth;
    cellresp(cell,:,:) = resp;
    celllfpresp(cell,:,:) = lfpresp;
    cellchan(cell,:) = chan;
    
    figure
    subplot(2,2,1)
    semilogy(trialfax,mean(lfpspect),'linewidth',2);
    semilogy(trialfax,mean(lfpspect)-(std(lfpspect)./sqrt(size(lfpspect,1))));
    semilogy(trialfax,mean(lfpspect)+(std(lfpspect)./sqrt(size(lfpspect,1))));
    axis([0,120,...
        min([min(squeeze(mean(lfpspect(:,1:125)))),min(squeeze(mean(lfpspect(:,1:125))))]),...
        max([max(squeeze(mean(lfpspect(:,1:125)))),max(squeeze(mean(lfpspect(:,1:125))))])])
    xlabel('frequency [Hz]')
    ylabel('spectral power')
    title('LFP spectrum during light on vs off')
    
    subplot(2,2,2)
    plot(mean(gammaresp))
    title(['gamma power in time depth: ' int2str(depth(cell))])
    
    subplot(2,2,3)
    [to,ro] = rose(phases{cell});
    polar(to,ro,'b')
    title(['gamma phase locking of unit ' int2str(cell) ' spikewidth: ' int2str(swidth(cell))])
    
    subplot(2,2,4)
    plot(mean(lfpresp));
    
    msta = linspace(-prestim,trialdur+poststim,size(resp,2));
    
    baseline = mean(bl);
    baselineerr = std(bl)./(sqrt(size(bl,1)));
    
    binwidth = 50;
    oris = unique(result.gratingInfo.Orientation); oris(find(oris == -1)) = [];
    clevels = unique(result.gratingInfo.Contrast); clevels(find(clevels == 0)) = []; %delete control condition
    
    for ori = 1:length(oris)
        for co = 1:length(clevels)
            thisinds = find(result.gratingInfo.Orientation == oris(ori) &...
                result.gratingInfo.Contrast == clevels(co));
            condn(ori,co) = length(thisinds);
            condresp(ori,co,:) = nanmean(resp(thisinds,:),1);
            condresperr(ori,co,:) = nanstd(resp(thisinds,:),1,1)./sqrt(length(thisinds));
            if ~isnan(condresp(ori,co,:))
                [bincondresp(ori,co,:),bta] = binit(condresp(ori,co,:),binwidth);
            else
                bincondresp(ori,co,:) = binit(condresp(ori,co,:),binwidth);
            end
            binconderr(ori,co,:) = binit(condresperr(ori,co,:),binwidth);
            condfr(ori,co) = nanmean(frs(thisinds));
            conderr(ori,co) =nanstd(frs(thisinds))./sqrt(length(thisinds));
        end
    end

    bincondresp = bincondresp.*(1000/binwidth);
    bta = bta-prestim;
    
%     contindsnl = find(result.gratingInfo.Contrast == 0);
%     controlresp = mean(resp(contindsnl,:),1);
%     controlfr(cell) = mean(frs(contindsnl));
%     controlerr = std(frs(contindsnl))./sqrt(length(contindsnl));
%     
%     [binctrl,bta] = binit(controlresp,binwidth); binctrl = binctrl.*(1000/binwidth);
    
    fc = find(clevels == 1);
    pc = find(squeeze(nanmedian(condfr(:,:),1)) == max(squeeze(nanmedian(condfr(:,:),1))),1,'last');
    
    [oriprefratio, dirprefratio, prefori, meanori, osi(cell), meandir, dsi(cell)] = getOSI(squeeze(condfr(:,fc))',oris);
    
    blscondfr = condfr-mean(bl);
    oneorifr = mean(reshape(blscondfr(:,pc),4,2),2);
    prefori = find(oneorifr(1,:) == max(oneorifr(1,:)),1);
    [preffr(cell,:), prefdir] = max(condfr(:,pc),[],1); prefdir = prefdir(1);
    ortho = mod(prefori+2,length(oris)/2); if ortho == 0, ortho = length(oris)/2; end
    
    [binpref,bta] = binit(squeeze(condresp(prefdir,pc,:)),binwidth); binpref = binpref.*(1000/binwidth);
    
    [binavg,bta] = binit(mean(resp),binwidth); binavg = binavg.*(1000/binwidth); % in Hz
    binstd = binit(std(resp)./sqrt(size(resp,1)),binwidth);
    ta = bta-prestim;

    
    % bin and fit a sine of correct temporal frequency
    binpref = binpref-mean(bl);  % subtract baseline
    stwin = find(ta>700&ta<1500);  % only take part of response after transient
    tx = ta(stwin);        % time axis for afer transient
    tempfreq = unique(result.gratingInfo.tFreq);
    if ~isnan(binpref)
        spsig = binpref(stwin); 
        spp = fit_fixedfsin(tx,spsig,tempfreq,sr); % the fit parameters ( Amplitude, Phase and Offset)
    else
        spp = [NaN,NaN,NaN];
    end
    
    f1f0(cell) = (spp(1))/spp(3); % Amplitude = F1, Offset = F0
    f1(cell) = spp(1); f0(cell) = spp(3);
    
%     figure
%     plot(ta,binpref); hold on; plot(tx,fixedfsin(spp,tx,tempfreq,sr),'b','linewidth',2);
%     plot(ta,binprefl1,'r'); plot(tx,fixedfsin(sppl1,tx,tempfreq,sr),'r','linewidth',2);
%     title(['f1/f0 L0= ', num2str(f1f0(cell)) '   f1/f0 L1 = ' num2str(f1f0l1(cell))])
    
    
    figure
    subplot(2,2,1)
    boundedline(ta,binavg,binstd,'k');
%     hold on
%     boundedline(ta,binavg,binstd,lcol);
    mx = max(binavg);
    axis([-prestim,trialdur+poststim,-0.05,mx]);
    line([0,0],[0,mx],'color','k','linewidth',2);
    line([2000,2000],[0,mx],'color','k','linewidth',2);
    line([500,500],[0,mx],'color','b','linewidth',2)
    line([1500,1500],[0,mx],'color','b','linewidth',2);
    xlabel('time [ms]')
    ylabel('firingrate [Hz]')
    
    subplot(2,2,2)
    errorbar(oris,squeeze(condfr(:,fc)),squeeze(conderr(:,fc)),'ko-','markersize',8,'linewidth',2)
    xlabel('shown orientation')
    ylabel('Firing rate [Hz]')
    set(gca,'xtick',oris)
    title(['OSI: ' num2str(osi(cell))])
    
%     cresp = [controlfr(cell), squeeze(mean(condfr))];
%     xlevels = [0,clevels];
%     
%     subplot(2,2,3)
%     errorbar(xlevels,cresp,[controlerr, squeeze(mean(conderr))],'ko','markersize',8)
%     hold on
%     xlabel('shown contrast')
%     ylabel('Firing rate [Hz]')
%     set(gca,'xtick',clevels); %[0, clevels])
%     ax = axis;
%     axis([-.1,1.1,-.1,ax(4)])
%     title(['contrast response all orientations ']);
    
%     params = fit_crf_NR(xlevels,cresp);
%     plotx = [0:0.01:1];
%     plot(plotx,NakaRushton(params,plotx),'k','linewidth',2);
%     
%     rmax(cell) = params(1); 
%     c50(cell) = params(3); 
%     r0(cell) = params(4);
%     
    subplot(2,4,7)
    plot(spike,'linewidth',2)
    axis([0,40,-100,100])
    legend(['width: ' int2str(swidth(cell)) ' adiff: ' num2str(adiff(cell))])
    
%     subplot(2,4,8)
%     plot(ta,binctrl,'linewidth',2)
%     axis([0,2500,min([min(binctrl)]),max([max(binctrl),0.1])])
   
    % running figure
    runfr = mean(frs(oktrials));
    norunfr = mean(frs(stilltrials));
    runfrerr = std(frs(oktrials))./sqrt(length(oktrials));
    norunfrerr = std(frs(stilltrials))./sqrt(length(stilltrials));
    
    run1 = frs(oktrials);
    run0 = frs(stilltrials);
    rmi(cell) = (runfr-norunfr)/(runfr+norunfr);
    
    figure
    subplot(2,2,1)
    imagesc(speed);
    colorbar
    title(['oktrials: ' int2str(length(oktrials)) '/' int2str(size(speed,1))])
    xlabel('time [ms]')
    ylabel('trial number')
    
    subplot(2,2,2)
    errorbar(msta,mean(speed),std(speed)./sqrt(size(speed,1)),'b')
    xlabel('time [ms]')
    ylabel('average runspeed')
    
    subplot(2,2,3)
    plot(mean(speed(:,respwin),2),frs,'.')
    xlabel('average runspeed of trial')
    ylabel('average firing rate of trial')
    
    subplot(2,2,4)
    barweb([runfr,norunfr],...
        [runfrerr,norunfrerr],...
        [],[{'running only'};{'immobile only'}],[],...
        [],'firing rate [Hz]',[],[]);
    
    cell = cell+1;
    disp('');
    
end

figure
plot(swidth,adiff,'k.')
xlabel('spike width')
ylabel('amplitude diff')

pfs = find(swidth<15);
prs = find(swidth>=15);
pfsv = swidth<15;
prsv = swidth>=15;

% adjust depth according to penetration angle
depth(v1) = depth(v1).*cosd(22);
depth(v2) = depth(v2).*cosd(45);
% get STA and SFC spike rate equalized
if sfc
    for cell = 1:size(cellresp,1)

        respl0 = squeeze(cellresp(cell,find(result.light == 0),respwin));
        respl1 = squeeze(cellresp(cell,find(result.light == 1),respwin));
        lfprespl0 = squeeze(celllfpresp(cell,find(result.light == 0),respwin));
        lfprespl1 = squeeze(celllfpresp(cell,find(result.light == 1),respwin));
        s1 = find(respl0); s2 = find(respl1);
        if length(s1)>length(s2)
            rp = randperm(length(s1));
            respl0(s1(rp(1:(length(s1)-length(s2))))) = 0;
        else
            rp = randperm(length(s2));
            respl1(s2(rp(1:(length(s2)-length(s1))))) = 0;
        end   

%         for i = 1:size(respl0,1)
%             [stal0(i,:),avgsnipspecl0(i,:),staspecl0(i,:),...
%                 sfcoherl0(i,:),nsfcfax, nspikesl0(i), snippetsl0{i}] = getsfc(respl0(i,:),...
%                 lfprespl0(i,:),150,sr);
%             [stal1(i,:),avgsnipspecl1(i,:),staspecl1(i,:),...
%                 sfcoherl1(i,:),nsfcfax, nspikesl1(i), snippetsl1{i}] = getsfc(respl1(i,:),...
%                 lfprespl1(i,:),150,sr);
%             if ~isnan(nsfcfax), sfcfax = nsfcfax; end % in case last trial has no spikes
%         end

        figure
        plot(sfcfax,nanmean(sfcoherl0),'linewidth',2)
        hold on
        plot(sfcfax,nanmean(sfcoherl1),'r','linewidth',2)
        plot(sfcfax,nanmean(sfcoherl0)-(nanstd(sfcoherl0)./sqrt(size(sfcoherl0,1))))
        plot(sfcfax,nanmean(sfcoherl0)+(nanstd(sfcoherl0)./sqrt(size(sfcoherl0,1))))
        plot(sfcfax,nanmean(sfcoherl1)-(nanstd(sfcoherl1)./sqrt(size(sfcoherl1,1))),'r')
        plot(sfcfax,nanmean(sfcoherl1)+(nanstd(sfcoherl1)./sqrt(size(sfcoherl1,1))),'r')
        ax = axis;
        axis([0,110,ax(3),ax(4)]);
        title(['SFC: depth: ' int2str(depth(cell)) '  swidth: ' num2str(swidth(cell))])
        xlabel('frequency [Hz]')
        ylabel('Coherence')

        sfcoherencel0(cell,:) = nanmean(sfcoherl0);
        sfcoherencel1(cell,:) = nanmean(sfcoherl1);
        sfcfx(cell,:) = sfcfax;
    end
    
    figure
    g1 = find(sfcfx(1,:)>40,1);
    g2 = find(sfcfx(1,:)>70,1)-1;
    plot(nanmean(sfcoherencel1(:,g1:g2),2)-nanmean(sfcoherencel0(:,g1:g2),2),depth,'.')
    axis ij
    hold on
    plot(nanmean(sfcoherencel1(pfs,g1:g2),2)-nanmean(sfcoherencel0(pfs,g1:g2),2),depth(pfs),'o')
    line([0,0],[300,1000],'color','k')
end

% calculate lower layer rs cell correlations with and without light
infrs = intersect(find(depth>500),prs);
irsresp = cellresp(infrs,:,:);
ncells = size(irsresp,1);

if length(infrs) ~= 0
    %method 1
    for i = 1:420
        cr(i,:,:) = corr(squeeze(irsresp(:,i,:))');
    end

    % method 2
    pair = 1;
    for i = 1:ncells-1
        for j = i+1:ncells

            % equalize FR of pair
            resp1 = squeeze(irsresp(i,:,:));
            resp2 = squeeze(irsresp(j,:,:));
            s1 = find(resp1); s2 = find(resp2);
            if length(s1)>length(s2)
                rp = randperm(length(s1)); % delete random spikes from resp1
                resp1(s1(rp(1:(length(s1)-length(s2))))) = 0;
            else
                rp = randperm(length(s2)); % delete random spikes from resp2
                resp2(s2(rp(1:(length(s2)-length(s1))))) = 0;
            end
            % now resp1 and 2 should have equal no of spikes
            %bin
            ccbin = 25;
            for ii = 1:size(resp1,2)/ccbin;
                binresp1(:,ii) = sum(resp1(:,(ii-1)*ccbin+1:ii*ccbin),2);
                binresp2(:,ii) = sum(resp2(:,(ii-1)*ccbin+1:ii*ccbin),2);
            end

            for st = 1:size(irsresp,2)
                [r,p] = corrcoef(squeeze(binresp1(st,:)),squeeze(binresp2(st,:)));
                rval(pair,st) = r(1,2);
                pval(pair,st) = p(1,2);
            end
            pair = pair+1;
        end
    end

    rl0 = rval(:,find(~result.light));
    rl1 = rval(:,find(result.light));

    [s,p] = ttest(nanmean(rl0,2),nanmean(rl1,2));
    figure
    plot(nanmean(rl0,2),nanmean(rl1,2),'k.')
    axis([-.2,1,-.2,1])
    line([-.2,1],[-.2,1],'color','k')
    axis([-.2,1,-.2,1])
    axis square
    line([-.2,1],[0,0],'color','k')
    line([0,0],[-.2,1],'color','k')
    xlabel('mean corrcoef over light OFF conditions')
    ylabel('mean corrcoef over light ON conditions')
    title(['all pairwise correlations between lower layer RS cells p: ' num2str(p)]);
    if printyn
        figSize = [30 21];
        set(gcf,'paperunits','centimeters','papersize',figSize,'paperposition',[0 0 figSize])
        print(gcf,[popprintpath ,  '18_LLRScorrelations.pdf'], '-dpdf' );
    end
end

figure
plot(c50(prs),depth(prs),'b.','linewidth',2)
hold on
plot(c50(pfs),depth(pfs),'r.','linewidth',2)
axis ij
axis([-.1,1.1,150,1050])
title('C50 in depth')
xlabel('C50')
ylabel('cortical depth [mum]')

figure
plot((lc50(prs)-c50(prs))./(lc50(prs)+c50(prs)),depth(prs),'bo','linewidth',2)
axis ij
line([0,0],[150,1050],'color','k')
axis([-1.1,1.1,150,1050])
title('C50 changes in depth RS cells')
xlabel('(C50 light - C50 no light)/(C50 light + C50 no light)')
ylabel('cortical depth [mum]')
if printyn
    figSize = [30 21];
    set(gcf,'paperunits','centimeters','papersize',figSize,'paperposition',[0 0 figSize])
    print(gcf,[popprintpath ,  '07_RSC50depth.pdf'], '-dpdf' );
end

figure
plot((lc50(pfs)-c50(pfs))./(lc50(pfs)+c50(pfs)),depth(pfs),'ro','linewidth',2)
axis ij
line([0,0],[150,1050],'color','k')
axis([-1.1,1.1,150,1050])
title('C50 changes in depth FS cells')
xlabel('(C50 light - C50 no light)/(C50 light + C50 no light)')
ylabel('cortical depth [mum]')
if printyn
    figSize = [30 21];
    set(gcf,'paperunits','centimeters','papersize',figSize,'paperposition',[0 0 figSize])
    print(gcf,[popprintpath ,  '08_FSC50depth.pdf'], '-dpdf' );
end

figure
plot(c50(prs),lc50(prs),'bo','linewidth',2)
hold on
plot(c50(pfs),lc50(pfs),'ro','linewidth',2)
axis([-.1,1.1,-.1,1.1])
line([-.1,1.1],[-.1,1.1],'color','k')
axis square
[s,pvrs] = ttest(c50(prs),lc50(prs));
[s,pvfs] = ttest(c50(pfs),lc50(pfs));
xlabel('c50 no light')
ylabel('c50 light on')
title(['C50 changes with light. p RS: ' num2str(pvrs) '  p FS: ' num2str(pvfs)])
if printyn
    figSize = [30 21];
    set(gcf,'paperunits','centimeters','papersize',figSize,'paperposition',[0 0 figSize])
    print(gcf,[popprintpath ,  '09_C50scatter.pdf'], '-dpdf' );
end

r0(find(r0<0)) = 0;
lr0(find(lr0<0)) = 0;

figure
plot(r0(prs),depth(prs),'bo','linewidth',2)
hold on
plot(r0(pfs),depth(pfs),'ro','linewidth',2)
axis ij
axis([-.1,30,150,1050])
title('r0 in depth')
xlabel('r0')
ylabel('cortical depth [mum]')

figure
plot((lr0(prs)-r0(prs))./(lr0(prs)+r0(prs)),depth(prs),'bo','linewidth',2)
axis ij
% axis([-1.1,1.1,150,1050])
line([0,0],[150,1050],'color','k')
title('r0 changes in depth RS cells')
xlabel('(r0 light - r0 no light)/(r0 light + r0 no light)')
ylabel('cortical depth [mum]')
if printyn
    figSize = [30 21];
    set(gcf,'paperunits','centimeters','papersize',figSize,'paperposition',[0 0 figSize])
    print(gcf,[popprintpath ,  '10_RSR0depth.pdf'], '-dpdf' );
end

figure
plot((lr0(pfs)-r0(pfs))./(lr0(pfs)+r0(pfs)),depth(pfs),'ro','linewidth',2)
axis ij
% axis([-1.1,1.1,150,1050])
line([0,0],[150,1050],'color','k')
title('r0 changes in depth FS cells')
xlabel('(r0 light - r0 no light)/(r0 light + r0 no light)')
ylabel('cortical depth [mum]')
if printyn
    figSize = [30 21];
    set(gcf,'paperunits','centimeters','papersize',figSize,'paperposition',[0 0 figSize])
    print(gcf,[popprintpath ,  '11_FSR0depth.pdf'], '-dpdf' );
end

figure
plot(r0(prs),lr0(prs),'bo','linewidth',2)
hold on
plot(r0(pfs),lr0(pfs),'ro','linewidth',2)
line([0,30],[0,30],'color','k')
axis square
[s,pvrs] = ttest(r0(prs),lr0(prs));
[s,pvfs] = ttest(r0(pfs),lr0(pfs));
xlabel('r0 no light')
ylabel('r0 light on')
title(['r0 changes with light. p RS: ' num2str(pvrs) '  p FS: ' num2str(pvfs)])
if printyn
    figSize = [30 21];
    set(gcf,'paperunits','centimeters','papersize',figSize,'paperposition',[0 0 figSize])
    print(gcf,[popprintpath ,  '12_R0scatter.pdf'], '-dpdf' );
end

figure
plot(rmax,depth,'bo','linewidth',2)
hold on
plot(rmax(pfs),depth(pfs),'ro','linewidth',2)
axis ij
ax = axis;
axis([ax(1),ax(2),150,1050])
title('rmax in depth')
xlabel('rmax')
ylabel('cortical depth [mum]')

figure
plot((lrmax(prs)-rmax(prs))./(lrmax(prs)+rmax(prs)),depth(prs),'bo','linewidth',2)
axis ij
% axis([-1.1,1.1,150,1050])
line([0,0],[150,1050],'color','k')
title('rmax changes in depth RS cells')
xlabel('(rmax light - rmax no light)/(rmax light + rmax no light)')
ylabel('cortical depth [mum]')
if printyn
    figSize = [30 21];
    set(gcf,'paperunits','centimeters','papersize',figSize,'paperposition',[0 0 figSize])
    print(gcf,[popprintpath ,  '13_RSrmaxdepth.pdf'], '-dpdf' );
end

figure
plot((lrmax(pfs)-rmax(pfs))./(lrmax(pfs)+rmax(pfs)),depth(pfs),'ro','linewidth',2)
axis ij
% axis([-1.1,1.1,150,1050])
line([0,0],[150,1050],'color','k')
title('rmax changes in depth FS cells')
xlabel('(rmax light - rmax no light)/(rmax light + rmax no light)')
ylabel('cortical depth [mum]')
if printyn
    figSize = [30 21];
    set(gcf,'paperunits','centimeters','papersize',figSize,'paperposition',[0 0 figSize])
    print(gcf,[popprintpath ,  '14_FSrmaxdepth.pdf'], '-dpdf' );
end

figure
plot(rmax(prs),lrmax(prs),'bo','linewidth',2)
hold on
plot(rmax(pfs),lrmax(pfs),'ro','linewidth',2)
line([0,max([max(rmax),max(lrmax)])+2],[0,max([max(rmax),max(lrmax)])+2],'color','k')
axis square
% axis([0,80,0,80])
[s,pvrs] = ttest(rmax(prs),lrmax(prs));
[s,pvfs] = ttest(rmax(pfs),lrmax(pfs));
xlabel('rmax no light')
ylabel('rmax light on')
title(['rmax changes with light. p RS: ' num2str(pvrs) ' p FS: ' num2str(pvfs)])
if printyn
    figSize = [30 21];
    set(gcf,'paperunits','centimeters','papersize',figSize,'paperposition',[0 0 figSize])
    print(gcf,[popprintpath ,  '15_Rmaxscatter.pdf'], '-dpdf' );
end

figure
[s,pvrs] = ttest(osi(prs),losi(prs));
[s,pvfs] = ttest(osi(pfs),losi(pfs));
plot(osi(prs),losi(prs),'bo','linewidth',2)
hold on
plot(osi(pfs),losi(pfs),'ro','linewidth',2)
line([0,1],[0,1],'color','k')
axis square
xlabel('OSI control');
ylabel('OSI light');
title(['RS p: ' num2str(pvrs) ' FS p: ' num2str(pvfs)])
if printyn
    figSize = [30 21];
    set(gcf,'paperunits','centimeters','papersize',figSize,'paperposition',[0 0 figSize])
    print(gcf,[popprintpath ,  '16_OSIscatter.pdf'], '-dpdf' );
end

figure
[s,pvrs] = ttest(dsi(prs),ldsi(prs));
[s,pvfs] = ttest(dsi(pfs),ldsi(pfs));
plot(dsi(prs),ldsi(prs),'bo','linewidth',2)
hold on
plot(dsi(pfs),ldsi(pfs),'ro','linewidth',2)
line([0,1],[0,1],'color','k')
axis square
xlabel('DSI control');
ylabel('DSI light');
title(['RS p: ' num2str(pvrs) ' FS p: ' num2str(pvfs)])
if printyn
    figSize = [30 21];
    set(gcf,'paperunits','centimeters','papersize',figSize,'paperposition',[0 0 figSize])
    print(gcf,[popprintpath ,  '17_DSIscatter.pdf'], '-dpdf' );
end

figure
plot((preffr(prs,2)-preffr(prs,1))./(preffr(prs,2)+preffr(prs,1)),depth(prs),'bo','linewidth',2)
line([0,0],[0,1000],'color','k')
axis ij
axis([-1.1,1.1,0,1000])
ylabel('depth[mum]')
xlabel('(light-nolight)/(light+nolight)')
title('RS cells: preferred firing rate changes by layer 4 suppression')
if printyn
    figSize = [30 21];
    set(gcf,'paperunits','centimeters','papersize',figSize,'paperposition',[0 0 figSize])
    print(gcf,[popprintpath ,  '01_RSpreffrdepth.pdf'], '-dpdf' );
end

figure
plot((preffr(pfs,2)-preffr(pfs,1))./(preffr(pfs,2)+preffr(pfs,1)),depth(pfs),'ro','linewidth',2)
line([0,0],[0,1000],'color','k')
axis ij
axis([-1.1,1.1,0,1000])
ylabel('depth[mum]')
xlabel('(light-nolight)/(light+nolight)')
title('FS cells: preferred firing rate changes by layer 4 suppression')
if printyn
    figSize = [30 21];
    set(gcf,'paperunits','centimeters','papersize',figSize,'paperposition',[0 0 figSize])
    print(gcf,[popprintpath ,  '02_FSpreffrdepth.pdf'], '-dpdf' );
end

figure
[s,pvrs] = ttest(preffr(prs,1),preffr(prs,2));
[s,pvfs] = ttest(preffr(pfs,1),preffr(pfs,2));
plot(preffr(prs,1),preffr(prs,2),'bo','linewidth',2)
hold on
plot(preffr(pfs,1),preffr(pfs,2),'ro','linewidth',2)
line([0,30],[0,30],'color','k')
axis square
xlabel('preferred firing rate no light')
ylabel('preferred firing rate light on')
title(['pref firing rate changes RS p: ' num2str(pvrs) '  FS p: ' num2str(pvfs)]);
if printyn
    figSize = [30 21];
    set(gcf,'paperunits','centimeters','papersize',figSize,'paperposition',[0 0 figSize])
    print(gcf,[popprintpath ,  '03_preffrscatter.pdf'], '-dpdf' );
end

figure
plot((controlfr(prs,2)-controlfr(prs,1))./(controlfr(prs,2)+controlfr(prs,1)),depth(prs),'bo','linewidth',2)
line([0,0],[0,1000],'color','k')
axis ij
axis([-1.1,1.1,0,1000])
ylabel('depth[mum]')
xlabel('(light-nolight)/(light+nolight)')
title('RS cells: control firing rate changes by layer 4 suppression')
if printyn
    figSize = [30 21];
    set(gcf,'paperunits','centimeters','papersize',figSize,'paperposition',[0 0 figSize])
    print(gcf,[popprintpath ,  '04_RScontrolfrdepth.pdf'], '-dpdf' );
end

figure
plot((controlfr(pfs,2)-controlfr(pfs,1))./(controlfr(pfs,2)+controlfr(pfs,1)),depth(pfs),'ro','linewidth',2)
line([0,0],[0,1000],'color','k')
axis ij
axis([-1.1,1.1,0,1000])
ylabel('depth[mum]')
xlabel('(light-nolight)/(light+nolight)')
title('FS cells: control firing rate changes by layer 4 suppression')
if printyn
    figSize = [30 21];
    set(gcf,'paperunits','centimeters','papersize',figSize,'paperposition',[0 0 figSize])
    print(gcf,[popprintpath ,  '05_FScontrolfrdepth.pdf'], '-dpdf' );
end

figure
[s,pvrs] = ttest(controlfr(prs,1),controlfr(prs,2));
[s,pvfs] = ttest(controlfr(pfs,1),controlfr(pfs,2));
plot(controlfr(prs,1),controlfr(prs,2),'bo','linewidth',2)
hold on
plot(controlfr(pfs,1),controlfr(pfs,2),'ro','linewidth',2)
line([0,30],[0,30],'color','k')
axis square
xlabel('control firing rate no light')
ylabel('control firing rate light on')
title(['control firing rate changes RS p: ' num2str(pvrs) ' FS p: ' num2str(pvfs)]);
if printyn
    figSize = [30 21];
    set(gcf,'paperunits','centimeters','papersize',figSize,'paperposition',[0 0 figSize])
    print(gcf,[popprintpath ,  '06_controlfrscatter.pdf'], '-dpdf' );
end

disp('')
end

function params=fit_crf_NR(x,y)
% [rmax, n, c50, r0]

    nrfitmargin = .01;
    range = max(y)-min(y);
    if range ~= 0
        p0 = [range 2 .5 min(y)];
        lb = [(1-nrfitmargin)*range, .5, 0, min(y)-.1*range];
        ub = [(1+nrfitmargin)*range, 10, 1, min(y)+.1*range];
        warning off
        params = lsqcurvefit(@(p,x) NakaRushton(p,x),p0,x(:),y(:),lb,ub,optimset('Display','off'));
        warning on
    else
        params = [NaN,NaN,NaN,NaN];
    end
end

function val=NakaRushton(p,x)
    % parameters of Naka-Rushton function as in Disney et al., Neuron, 2007
    % [R_max, contrast Exponent n,  50%firing-Contrast, spontaneous rate sFR]

    val = p(4)+p(1)*((x.^p(2))./(x.^p(2)+p(3).^p(2)));
end

function p = fit_fixedfsin(x,y,f,sr)
    %[A, ph, offs]
    range = max(y)-min(y);
    if range~=0
        p0 = [range/2, pi, (max(y)+min(y))/2]; 
        lb = [.5*(range/2),0,min(y)];
        ub = [1.5*(range/2),2*pi,max(y)];
        p = lsqcurvefit(@(p,x) fixedfsin(p,x,f,sr), p0,x,y,lb,ub,optimset('Display','off'));
    else
        p = [NaN,NaN,NaN];
    end
end

function val = fixedfsin(p,x,f,sr)
    %[A f ph offs]
    val = p(1) * sin(x*((f*2*pi)/sr) + p(2)) + p(3);
end
