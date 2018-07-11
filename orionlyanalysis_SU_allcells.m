function orionlyanalysis_SU_allcells

animalid = '150916';
block = 3;
lcol = 'r';

onlymod = 0;
printyn = 0;
sfc = 1;


basepath = ['C:\Users\Julia\work\data\' animalid '\'];
supath = [basepath 'singleunits\'];
basename = [animalid '_block' int2str(block) '_tet'];

if printyn
%     if ~exist([basepath, 'pdfs/'],'dir'),mkdir([basepath, 'pdfs/']),end
%     if ~exist([basepath, 'pdfs/' 'LFPs/'],'dir'),mkdir([basepath, 'pdfs/' 'LFPs/']),end
%     if ~exist([basepath, 'pdfs/' 'spikes/'],'dir'),mkdir([basepath, 'pdfs/' 'spikes/']),end
%     if ~exist([basepath, 'pdfs/' 'running/'],'dir'),mkdir([basepath, 'pdfs/' 'running/']),end
end

lfpprintpath = [basepath 'pdfs\LFPs\'];
spikeprintpath = [basepath 'pdfs\spikes\'];
runprintpath = [basepath 'pdfs\running\'];
popprintpath = [basepath 'pdfs\'];

files = dir([supath, basename, '*.mat']);

prestim = 300;
poststim = 300;
respwin = 501:1500; % after stimulus onset
respwin = respwin+prestim;
freqbinwidth = 5;

cell = 1;
for fi = 1:length(files)
    
%     if strfind(files(fi).name, 'MU')
%         continue;
%     end
    
    load([supath, files(fi).name]);
    
    % calc spiking to see if includable
    msStimes = round(result.spikes);
    if ~isempty(msStimes) & msStimes(1) == 0, msStimes(1) = 1; end
    
    chan = zeros(1,length(result.lfp));
    chan(msStimes) = 1;
    
    isi = diff(msStimes);
    bursts = legendy_new3(isi,3,1000,3,0,15); %(ISI, fac, sr, min_length_of_burst, local_length, surprise_cutoff)
    surp = [bursts.surprise];
    bursts(surp == 100) = []; % delete probably wrong bursts
    burstbegs = [bursts.begin];
    burstchan = zeros(1,length(result.lfp));
    burstchan(msStimes(burstbegs)) = 1;
    
    wvchan = find(var(result.waveforms) == max(var(result.waveforms)));
    
    sr = 1000;
    lfp = result.lfp(:,wvchan)';
    nfft = 2^nextpow2(length(lfp));
    fax = sr/2*linspace(0,1,nfft/2+1);
    y = fft(lfp,nfft);
    lfpspectrum = abs(y(1:nfft/2+1));
    

%     plot(fax,lfpspectrum)
    gamma = eegfilt(lfp,sr,15,45);
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
%         msstamps(358) = []; % for 150603 block 4
%         result.msstamps = msstamps;
%         save([supath, files(fi).name],'result');
%         msstamps([83]) = []; % for 160125 block 7
%         msstamps([290]) = []; % for 150916 block 3
%         result.msstamps = msstamps;
%         save([supath, files(fi).name],'result');
        pause;
    end
%     trialnfft = 2^nextpow2(800);
%     trialfax = sr/2*linspace(0,1,trialnfft/2+1);

    clear lfpspect; clear lfpoffsspect; clear fax;
    for i = 1:length(msstamps)
        resp(i,:) = chan(msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim); 
        burstresp(i,:) = burstchan(msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim); 
        lfpresp(i,:) = result.lfp(msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim);
%         y = fft(lfpresp(i,1001:1800),trialnfft);
%         lfpspect(i,:) = abs(y(1:trialnfft/2+1));
        [lfpspect(i,:),trialfax] = pmtm(lfpresp(i,1001:1800),2,[],sr);
        [lfpoffspect(i,:), trialfax] = pmtm(lfpresp(i,1801:2600),3,[],sr);
        gammaresp(i,:) = gpow(msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim);
        gphaseresp(i,:) = gphas(msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim);
        bpresp(i,:) = gamma(msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim); % bandpassed
        isiresp{i} = diff(find(resp(i,respwin)));
        
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
    gp = mean(gammaresp(:,respwin),2);
    
    %determine if cell is visually modulated
    blfr = sum(resp(:,1:prestim),2);
    vrfr = sum(resp(:,prestim+40:2*prestim+40-1),2);
    vismod(fi) = ttest2(blfr,vrfr);
    visdriven(fi) = mean(vrfr)>=mean(blfr)+2; % average firing rate is increase at least 2Hz above baseline
    
    %determine if cell is modulated by light
    lightmod(fi) = ttest2(frs(find(result.light)),frs(find(result.light == 0)));
    
    if onlymod & ~vismod(fi) %~(lightmod(fi) & vismod(fi))
        continue;
    end
    
    cellname{cell} = files(fi).name;
    
    i = strfind(files(fi).name, 'tet');
    if strcmp(files(fi).name(i+4),'_')
        tetno = strread(files(fi).name(i+3)); % single character number
    else        
        tetno = strread(files(fi).name(i+3:i+4)); % number >10
    end
    tetnos(cell) = tetno;
    if tetno>8
        v1(cell) = logical(0); v2(cell) = logical(1); cellstr = 'V2';
    else
        v1(cell) = logical(1); v2(cell) = logical(0); cellstr = 'V1';
    end
    
    spike = result.waveforms(:,wvchan);
    interpspike = spline(1:32,spike,1:.1:32);
    [adiff(cell),swidth(cell),ptr(cell),endslope(cell),fwhh(cell)] = spikequant(interpspike);  %[adiff,swidth,ptr,endslope,fwhh]
    
    % phases
    tmp = zeros(size(gphaseresp));
    tmp(find(resp)) = gphaseresp(find(resp));
    l0phasemat = tmp(find(result.light == 0),:);
    l1phasemat = tmp(find(result.light == 1),:);
    l0phases{cell} = l0phasemat(find(l0phasemat));
    l1phases{cell} = l1phasemat(find(l1phasemat));
      
    rl0(cell) = circ_r(l0phases{cell});
    rl1(cell) = circ_r(l1phases{cell});
    cmeanl0(cell) = circ_mean(l0phases{cell});
    cmeanl1(cell) = circ_mean(l1phases{cell});
    
    if sfc
        tmpi = zeros(size(allphaseresp));
        for i = 1:size(allphaseresp,1)
            tmp = zeros(size(gphaseresp));
            tmp(find(resp)) = allphaseresp(i,find(resp));
            tmpi(i,:,:) = tmp;
        end
        alll0phasemat = tmpi(:,find(result.light == 0),:);
        alll1phasemat = tmpi(:,find(result.light == 1),:);
        for i = 1:size(allphaseresp,1)
            allphasesl0{i} = alll0phasemat(i,find(squeeze(alll0phasemat(i,:,:)))); 
            allphasesl1{i} = alll1phasemat(i,find(squeeze(alll1phasemat(i,:,:)))); 
            allrl0(cell,i) = circ_r(allphasesl0{i}');
            allrl1(cell,i) = circ_r(allphasesl1{i}');
            allcmeanl0(cell,i) = circ_mean(allphasesl0{i}');
            allcmeanl1(cell,i) = circ_mean(allphasesl1{i}');
        end
    end
    
    
    %%%%%%%%%%%%%%%%%%%%
%     % runspeed
%     if ~isfield(result,'runspeed')
%         x_t = cumsum(result.speed);
%         [smooth_win, FWHM] = MakeGaussWindow(round(1000),23.5/2, 1000); 
%         sw_len = length(smooth_win);
% 
%         x_t(end+1:end+sw_len) = x_t(end); %pad with last value for length of smoothing kernel
%         d_smooth_win = [0;diff(smooth_win)]/(1/1000);
%         dx_dt = conv(x_t,d_smooth_win,'same');
%         dx_dt(end-sw_len+1:end) = []; %remove values produced by convolving kernel with padded values
%         mouserad = 6; %cm
%         scalefact = 2*pi*mouserad/360;
%         result.runspeed = dx_dt.*scalefact;
%         save([supath, files(cell).name], 'result');
%     end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     result.msstamps = msstamps;
%     result.speed = speed;
%     result.speedta = speedta;
%     
%     maxspikechannel = find(var(result.waveforms) == max(var(result.waveforms)));
%     channo = (tetno-1)*4+maxspikechannel;
%     result.depth = filedepth - (channo-1)*contactspacing;
%     
%     chans = readTrodesFileChannels(filename,(tetno-1)*4+1:(tetno-1)*4+4);
%     lfp = eegfilt(chans.channelData',30000,0,200)';
%     result.lfp = resample(lfp,1,30);
%     
%     save([supath, files(cell).name], 'result');
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%iledepth - (channo-1)*contactspacing;
    
    depth(cell) = result.depth;
    cellresp(cell,:,:) = resp;
    celllfpresp(cell,:,:) = lfpresp;
    cellchan(cell,:) = chan;
    
    figure
    subplot(2,2,1)
    semilogy(trialfax,mean(lfpspect(find(~result.light),:)),'linewidth',2);
    hold on
    semilogy(trialfax,mean(lfpspect(find(result.light),:)),'r','linewidth',2);
    semilogy(trialfax,mean(lfpoffspect(find(result.light),:)),'c','linewidth',2);
    semilogy(trialfax,mean(lfpspect(find(~result.light),:))-(std(lfpspect(find(~result.light),:))./sqrt(length(find(~result.light)))));
    semilogy(trialfax,mean(lfpspect(find(~result.light),:))+(std(lfpspect(find(~result.light),:))./sqrt(length(find(~result.light)))));
    semilogy(trialfax,mean(lfpspect(find(result.light),:))-(std(lfpspect(find(result.light),:))./sqrt(length(find(result.light)))),'r');
    semilogy(trialfax,mean(lfpspect(find(result.light),:))+(std(lfpspect(find(result.light),:))./sqrt(length(find(result.light)))),'r');
    semilogy(trialfax,mean(lfpoffspect(find(result.light),:))-(std(lfpoffspect(find(result.light),:))./sqrt(length(find(result.light)))),'c');
    semilogy(trialfax,mean(lfpoffspect(find(result.light),:))+(std(lfpoffspect(find(result.light),:))./sqrt(length(find(result.light)))),'c');
    axis([0,120,...
        min([min(squeeze(mean(lfpspect(find(result.light == 0),1:125)))),min(squeeze(mean(lfpspect(find(result.light == 1),1:125))))]),...
        max([max(squeeze(mean(lfpspect(find(result.light == 0),1:125)))),max(squeeze(mean(lfpspect(find(result.light == 1),1:125))))])])
    legend({'light off', 'light on'});
    xlabel('frequency [Hz]')
    ylabel('spectral power')
    title('LFP spectrum during light on vs off')
    
    subplot(2,2,2)
    plot(mean(gammaresp(find(result.light == 0),:)))
    hold on
    plot(mean(gammaresp(find(result.light),:)),'r')
    title(['gamma power in time depth: ' int2str(depth(cell))])
    
    subplot(2,2,3)
    [to0,ro0] = rose(l0phases{cell}); [to1,ro1] = rose(l1phases{cell});
    if max(ro1)>max(ro0)
        polar(to1,ro1,'r')
        hold on
        polar(to0,ro0,'b')
    else
        polar(to0,ro0,'b')
        hold on
        polar(to1,ro1,'r')
    end
    title(['gamma phase locking of unit ' int2str(cell) ' spikewidth: ' int2str(swidth(cell))])
    
    subplot(2,2,4)
    plot(mean(lfpresp(find(result.light == 0),:)));
    hold on
    plot(mean(lfpresp(find(result.light),:)),'r')
    
%     if printyn
%         figSize = [30 21];
%         set(gcf,'paperunits','centimeters','papersize',figSize,'paperposition',[0 0 figSize])
%         if cell<10, printi = ['0', int2str(cell)]; else printi = int2str(cell); end
%         print(gcf,[lfpprintpath ,  printi '_' files(fi).name '.pdf'], '-dpdf' );
%     end
    
    msta = linspace(-prestim,trialdur+poststim,size(resp,2));
    
    lightresp = resp(find(result.light),:);
    nolightresp = resp(find(result.light == 0),:);
    lightlfpresp = lfpresp(find(result.light),:);
    nolightlfpresp = lfpresp(find(~result.light),:);
    
%     frs = sum(resp(:,respwin),2)./(length(respwin)/1000);
%     bl = sum(resp(:,1:prestim),2)./(prestim/1000);
%     
%     %determine if cell is visually modulated
%     blfr = sum(resp(:,1:prestim),2);
%     vrfr = sum(resp(:,prestim+40:2*prestim+40),2);
%     vismod(cell) = ttest2(blfr,vrfr);
%     
%     %determine if cell is modulated by light
%     lightmod(cell) = ttest2(frs(find(result.light)),frs(find(result.light == 0)));
    
    printname = files(fi).name;
    printname(find(printname=='_')) = ' ';
    
    nolbl = mean(bl(find(~result.light)));
    nolblerr = std(bl(find(~result.light)))./(sqrt(length(find(~result.light))));
    lbl = mean(bl(find(result.light)));
    lblerr = std(bl(find(result.light)))./(sqrt(length(find(result.light))));
    
    binwidth = 50;
    oris = unique(result.gratingInfo.Orientation); oris(find(oris == -1)) = [];
    clevels = unique(result.gratingInfo.Contrast); clevels(find(clevels == 0)) = []; %delete control condition
    for l = 1:length(unique(result.light))
        for ori = 1:length(oris)
            
            thisinds = find(result.gratingInfo.Orientation == oris(ori) &  result.light == l-1);
            condn(l,ori) = length(thisinds);
            condresp(l,ori,:) = nanmean(resp(thisinds,:),1);
            condresperr(l,ori,:) = nanstd(resp(thisinds,:),1,1)./sqrt(length(thisinds));
            if ~isnan(condresp(l,ori,:))
                [bincondresp(l,ori,:),bta] = binit(condresp(l,ori,:),binwidth);
            else
                bincondresp(l,ori,:) = binit(condresp(l,ori,:),binwidth);
            end
            binconderr(l,ori,:) = binit(condresperr(l,ori,:),binwidth);
            condfr(l,ori) = nanmean(frs(thisinds));
            conderr(l,ori) =nanstd(frs(thisinds))./sqrt(length(thisinds));
            
            condgp(l,ori) = nanmean(gp(thisinds));
            condlfpspect(l,ori,:) = nanmean(lfpspect(thisinds,:),1);
            condlfpspecterr(l,ori,:) = nanstd(lfpspect(thisinds,:))./sqrt(length(thisinds));
            
            fillspec(l,ori,:) = [squeeze(condlfpspect(l,ori,1:104)+condlfpspecterr(l,ori,1:104))',fliplr(squeeze(condlfpspect(l,ori,1:104)-condlfpspecterr(l,ori,1:104))')];
            
        end
    end
    bincondresp = bincondresp.*(1000/binwidth);
    bta = bta-prestim;
%     psthcondplot(bincondresp,binconderr,bta);

    fillx = [trialfax(1:104)',fliplr(trialfax(1:104)')];
    

    l0isi = []; l1isi = [];
    for i = find(result.light == 0)
        l0isi = [l0isi, isiresp{i}];
    end
    for i = find(result.light == 1)
        l1isi = [l1isi, isiresp{i}];
    end
    
    cfr(cell,:,:) = condfr;
    cellcondlfpspect(cell,:,:,:) = condlfpspect;
    
    if size(condfr,2)>1
        if size(condfr,2)>1
            [nlprefratio, nldirpreffratio, nlprefori, nlmeanori, nlosi(cell), nlmeandir(cell), nldsi(cell)] = getOSI(squeeze(condfr(1,:)),oris);
            if size(condfr,1)>1
                [lprefratio, ldirprefratio, lprefori, lmeanori, losi(cell), lmeandir(cell), ldsi(cell)] = getOSI(squeeze(condfr(2,:)),oris);
            else
                losi(cell) = NaN; dosi(cell) = NaN;
            end
        end
    end
    
    blscondfr = condfr-mean(bl);
    if size(condfr,2)>1
        oneorifr = mean(reshape(blscondfr,size(blscondfr,1),length(oris)/2,2),3);
        prefori = find(oneorifr(1,:) == max(oneorifr(1,:)),1);
        [preffr(cell,:), prefdir] = max(condfr(:,:),[],2); prefdir = prefdir(1);
        ortho = mod(prefori+2,length(oris)/2); if ortho == 0, ortho = length(oris)/2; end
    else
        prefdir = 1;
    end
    
    [binprefl0,bta] = binit(squeeze(condresp(1,prefdir,:)),binwidth); binprefl0 = binprefl0.*(1000/binwidth);
    if size(condresp,1)>1
        [binprefl1,bta] = binit(squeeze(condresp(2,prefdir,:)),binwidth); binprefl1 = binprefl1.*(1000/binwidth);
    else
        binprefl1 = nan(size(binprefl0));
    end
    [binavgl1,bta] = binit(mean(lightresp),binwidth); binavgl1 = binavgl1.*(1000/binwidth); % in Hz
    [binavgl0,bta] = binit(mean(nolightresp),binwidth); binavgl0 = binavgl0.*(1000/binwidth);
    binstdl1 = binit(std(lightresp)./sqrt(size(lightresp,1)),binwidth);
    binstdl0 = binit(std(nolightresp)./sqrt(size(nolightresp,1)),binwidth);
    ta = bta-prestim;

    
    % bin and fit a sine of correct temporal frequency
    binprefl0 = binprefl0-mean(bl);  % subtract baseline
    binprefl1 = binprefl1-mean(bl);
    stwin = find(ta>700&ta<1500);  % only take part of response after transient
    tx = ta(stwin);        % time axis for afer transient
    tempfreq = unique(result.gratingInfo.tFreq);
    if ~isnan(binprefl0)
        spsigl0 = binprefl0(stwin); 
        sppl0 = fit_fixedfsin(tx,spsigl0,tempfreq,sr); % the fit parameters ( Amplitude, Phase and Offset)
    else
        sppl0 = [NaN,NaN,NaN];
    end
    if ~isnan(binprefl1)
        spsigl1 = binprefl1(stwin); 
        sppl1 = fit_fixedfsin(tx,spsigl1,tempfreq,sr);
    else
        sppl1 = [NaN,NaN,NaN];
    end
    
    f1f0l0(cell) = (sppl0(1))/sppl0(3); % Amplitude = F1, Offset = F0
    f1f0l1(cell) = (sppl1(1))/sppl1(3);
    f1l0(cell) = sppl0(1); f1l1(cell) = sppl1(1); f0l0(cell) = sppl0(3); f0l1(cell) = sppl1(3);
    
%     figure
%     plot(ta,binprefl0); hold on; plot(tx,fixedfsin(sppl0,tx,tempfreq,sr),'b','linewidth',2);
%     plot(ta,binprefl1,'r'); plot(tx,fixedfsin(sppl1,tx,tempfreq,sr),'r','linewidth',2);
%     title(['f1/f0 L0= ', num2str(f1f0l0(cell)) '   f1/f0 L1 = ' num2str(f1f0l1(cell))])
   
    
    figure
    subplot(2,2,1)
    boundedline(ta,binavgl0,binstdl0,'k');
    hold on
    boundedline(ta,binavgl1,binstdl1,lcol);
    mx = max(binavgl0);
    axis([-prestim,trialdur+poststim,-0.05,mx]);
    line([0,0],[0,mx],'color','k','linewidth',2);
    line([2000,2000],[0,mx],'color','k','linewidth',2);
    line([500,500],[0,mx],'color','b','linewidth',2)
    line([1500,1500],[0,mx],'color','b','linewidth',2);
    legend({'Light OFF','Light ON'})
    xlabel('time [ms]')
    ylabel('firingrate [Hz]')
    title(['cell ' int2str(cell) ' depth: ' int2str(result.depth), 'cell ' printname ])
    
    subplot(2,2,2)
    if size(condfr,1)>1
    	errorbar(oris,squeeze(condfr(2,:)),squeeze(conderr(2,:)),'o-','color',lcol,'markersize',8,'linewidth',2)
    end
    hold on
    errorbar(oris,squeeze(condfr(1,:)),squeeze(conderr(1,:)),'ko-','markersize',8,'linewidth',2)
    xlabel('shown orientation')
    ylabel('Firing rate [Hz]')
    set(gca,'xtick',oris)
    legend({'Light ON','Light OFF'})
    if size(condfr,2)>1
        title(['OSI: ' num2str(nlosi(cell)) ' light OSI: ' num2str(losi(cell)) '   ' printname])
    end
    
    if size(condfr,2)>1
        oneorifr = nanmean(reshape(condfr,size(condfr,1),length(oris)/2,2),3);
        prefori = find(oneorifr(1,:) == max(oneorifr(1,:)),1);
        ortho = mod(prefori+2,length(oris)/2); if ortho == 0, ortho = length(oris)/2; end
    else
        prefori = 1;
    end
    
    
    subplot(2,4,7)
    plot(spike,'linewidth',2)
    axis([0,40,-100,100])
    legend(['width: ' int2str(swidth(cell)) ' adiff: ' num2str(adiff(cell))])
    
    if printyn
        figSize = [30 21];
        set(gcf,'paperunits','centimeters','papersize',figSize,'paperposition',[0 0 figSize])
        if cell<10, printi = ['0', int2str(cell)]; else printi = int2str(cell); end
        print(gcf,[spikeprintpath ,  printi '_' files(fi).name '.pdf'], '-dpdf' );
    end

%     o = 1; c = 3;
%     condl0 = (result.gratingInfo.Orientation == oris(o) | result.gratingInfo.Orientation == oris(o+length(oris)/2)) & (result.gratingInfo.Contrast == clevels(c) | result.gratingInfo.Contrast == clevels(c+1)) & result.light == 0;
%     condl1 = (result.gratingInfo.Orientation == oris(o) | result.gratingInfo.Orientation == oris(o+length(oris)/2)) & (result.gratingInfo.Contrast == clevels(c) | result.gratingInfo.Contrast == clevels(c+1)) & result.light == 1;
%     condl0 = result.gratingInfo.Orientation == oris(o) & result.gratingInfo.size == sizes(s) & result.light == 0;
%     condl1 = result.gratingInfo.Orientation == oris(o) & result.gratingInfo.size == sizes(s) & result.light == 1;
%     condl0 = result.light == 0; condl1 = result.light == 1;
%     figure, rasterplot(resp,condl0,condl1,msta);
    line([0,2000],[37,37],'color','k','linewidth',2)
    line([500,1500],[33,33],'color','r','linewidth',2)
    
    preffr(cell,:) = condfr(:,prefori);

    % running figure
    lfr = mean(frs(find(result.light)));
    nlfr = mean(frs(find(~result.light)));
    runlfr = mean(frs(intersect(find(result.light),oktrials)));
    runnlfr = mean(frs(intersect(find(~result.light),oktrials)));
    norunlfr = mean(frs(intersect(find(result.light),stilltrials)));
    norunnlfr = mean(frs(intersect(find(~result.light),stilltrials)));
    lfrerr = std(frs(find(result.light)))./sqrt(length(find(result.light)));
    nlfrerr = std(frs(find(~result.light)))./sqrt(length(find(~result.light)));
    runlfrerr = std(frs(intersect(find(result.light),oktrials)))./sqrt(length(intersect(find(result.light),oktrials)));
    runnlfrerr = std(frs(intersect(find(~result.light),oktrials)))./sqrt(length(intersect(find(~result.light),oktrials)));
    norunlfrerr = std(frs(intersect(find(result.light),stilltrials)))./sqrt(length(intersect(find(result.light),stilltrials)));
    norunnlfrerr = std(frs(intersect(find(~result.light),stilltrials)))./sqrt(length(intersect(find(~result.light),stilltrials)));
    
    l1r1 = frs(intersect(find(result.light),oktrials));
    l0r1 = frs(intersect(find(~result.light),oktrials));
    l1r0 = frs(intersect(find(result.light),stilltrials));
    l0r0 = frs(intersect(find(~result.light),stilltrials));
    anovavec = [l0r0;l0r1;l1r0;l1r1]; 
    g1 = [zeros(length(l0r0),1);zeros(length(l0r1),1);ones(length(l1r0),1);ones(length(l1r1),1)]; %light
    g2 = [zeros(length(l0r0),1);ones(length(l0r1),1);zeros(length(l1r0),1);ones(length(l1r1),1)]; %running
    [p,table,stats] = anovan(anovavec,{g1 g2},'model','full','display','off');
    lp(cell) = p(1); rp(cell) = p(2); rlip(cell) = p(3);
    
    r0omi(cell) = (norunlfr-norunnlfr)/(norunlfr+norunnlfr);
    r1omi(cell) = (runlfr-runnlfr)/(runlfr+runnlfr);
    l0rmi(cell) = (runnlfr-norunnlfr)/(runnlfr+norunnlfr);
    l1rmi(cell) = (runlfr-norunlfr)/(runlfr+norunlfr);
    
    figure
    subplot(2,2,1)
    imagesc(speed);
    colorbar
    title(['oktrials: ' int2str(length(oktrials)) '/' int2str(size(speed,1))])
    xlabel('time [ms]')
    ylabel('trial number')
    
    subplot(2,2,2)
    errorbar(msta,mean(speed(find(~result.light),:)),std(speed(find(~result.light),:))./sqrt(length(find(~result.light))),'b')
    hold on
    errorbar(msta,mean(speed(find(result.light),:)),std(speed(find(result.light),:))./sqrt(length(find(result.light))),'r')
    xlabel('time [ms]')
    ylabel('average runspeed')
    legend({'light off' 'light on'})
    
    subplot(2,2,3)
    plot(mean(speed(:,respwin),2),frs,'.')
    hold on
    plot(mean(speed(find(result.light),respwin),2),frs(find(result.light)),'r.')
    xlabel('average runspeed of trial')
    ylabel('average firing rate of trial')
    
%     subplot(2,2,4)
%     barweb([nlfr,lfr;runnlfr,runlfr;norunnlfr,norunlfr],...
%         [nlfrerr,lfrerr;runnlfrerr,runlfrerr;norunnlfrerr,norunlfrerr],...
%         [],[{'all'};{'running only'};{'immobile only'}],['ANOVA factor running p: ' num2str(p(2))],...
%         [],'firing rate [Hz]',[],[]);
    
%     if printyn
%         figSize = [30 21];
%         set(gcf,'paperunits','centimeters','papersize',figSize,'paperposition',[0 0 figSize])
%         if cell<10, printi = ['0', int2str(cell)]; else printi = int2str(cell); end
%         print(gcf,[runprintpath ,  printi '_' files(fi).name '.pdf'], '-dpdf' );
%     end
            
%     subplot(2,2,4)
%     errorbar([0,clevels],[lbl,squeeze(mean(condfr(2,[prefori,prefori+(length(oris)/2)],:),2))'],...
%         [lblerr, squeeze(mean(conderr(2,[prefori,prefori+(length(oris)/2)],:),2))'],'o-','color',lcol,'markersize',8,'linewidth',2)
%     hold on
%     errorbar([0,clevels],[nolbl,squeeze(mean(condfr(1,[prefori,prefori+(length(oris)/2)],:),2))'],...
%         [nolblerr, squeeze(mean(conderr(1,[prefori,prefori+(length(oris)/2)],:),2))'],'ko-','markersize',8,'linewidth',2)
%     errorbar([0,clevels],[lbl,squeeze(mean(condfr(2,[ortho,ortho+(length(oris)/2)],:),2))'],...
%         [lblerr, squeeze(mean(conderr(2,[ortho,ortho+(length(oris)/2)],:),2))'],'o--','color',lcol,'markersize',8,'linewidth',1)
%     hold on
%     errorbar([0,clevels],[nolbl,squeeze(mean(condfr(1,[ortho,ortho+(length(oris)/2)],:),2))'],...
%         [nolblerr, squeeze(mean(conderr(1,[ortho,ortho+(length(oris)/2)],:),2))'],'ko--','markersize',8,'linewidth',1)
%     xlabel('shown contrast')
%     ylabel('Firing rate [Hz]')
%     legend({'Light ON preferred','Light OFF preferred', 'Light ON orthogonal', 'Light OFF orthogonal'})    
%     set(gca,'xtick',clevels); %[0, clevels])
%     title(['contrast response preferred orientations ' printname]);
    
    cell = cell+1;
    disp('');
    
end

figure
plot(swidth,adiff,'k.')
xlabel('spike width')
ylabel('amplitude diff')

pfs = find(swidth<130);
prs = find(swidth>=130);

%coherence
cell1 = 12; cell2 = 21;
for i = 1:size(celllfpresp,2)
    [coh(i,:),cfx] = mscohere(squeeze(celllfpresp(cell1,i,respwin)),squeeze(celllfpresp(cell2,i,respwin)),[],[],512,1000);
end
for l = 1:2
    for ori = 1:length(oris)
        for co = 1:length(clevels)
            thisinds = find(result.gratingInfo.Orientation == oris(ori) &...0
                result.gratingInfo.Contrast ==clevels(co) & ...
                result.light == l-1);
            condcoher(l,ori,co,:) = nanmean(coh(thisinds,:),1);
        end
    end
end
            


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

        for i = 1:size(respl0,1)
            [stal0(i,:),avgsnipspecl0(i,:),staspecl0(i,:),...
                sfcoherl0(i,:),nsfcfax, nspikesl0(i), snippetsl0{i}] = getsfc(respl0(i,:),...
                lfprespl0(i,:),150,sr);
            [stal1(i,:),avgsnipspecl1(i,:),staspecl1(i,:),...
                sfcoherl1(i,:),nsfcfax, nspikesl1(i), snippetsl1{i}] = getsfc(respl1(i,:),...
                lfprespl1(i,:),150,sr);
            if ~isnan(nsfcfax), sfcfax = nsfcfax; end % in case last trial has no spikes
        end

        figure
        plot(sfcfax,nanmean(sfcoherl0),'linewidth',2)
        hold on
        plot(sfcfax,nanmean(sfcoherl1),'r','linewidth',2)
        plot(sfcfax,nanmean(sfcoherl0)-(nanstd(sfcoherl0)./sqrt(size(sfcoherl0,1))))
        plot(sfcfax,nanmean(sfcoherl0)+(nanstd(sfcoherl0)./sqrt(size(sfcoherl0,1))))
        plot(sfcfax,nanmean(sfcoherl1)-(nanstd(sfcoherl1)./sqrt(size(sfcoherl1,1))),'r')
        plot(sfcfax,nanmean(sfcoherl1)+(nanstd(sfcoherl1)./sqrt(size(sfcoherl1,1))),'r')
        ax = axis;
        axis([0,110,0,1]); %ax(3),ax(4)
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
    
    for i = 1:length(infrs)
        for j = 1:size(irsresp,2)
            binresp(i,j,:) = binit(squeeze(irsresp(i,j,:)),30);
        end
    end
    
    %method 1
    for i = 1:340
        cr(i,:,:) = corr(squeeze(irsresp(:,i,:))');
        bcr(i,:,:) = corr(squeeze(binresp(:,i,:))');
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
            ccbin = 10;
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
%     if printyn
%         figSize = [30 21];
%         set(gcf,'paperunits','centimeters','papersize',figSize,'paperposition',[0 0 figSize])
%         print(gcf,[popprintpath ,  '18_LLRScorrelations.pdf'], '-dpdf' );
%     end
    
    rl0c0 = nanmean(rval(:,result.light == 0 & result.gratingInfo.Contrast == 0),2);
    rl0c1 = nanmean(rval(:,result.light == 0 & result.gratingInfo.Contrast == 1),2);
    rl1c0 = nanmean(rval(:,result.light == 1 & result.gratingInfo.Contrast == 0),2);
    rl1c1 = nanmean(rval(:,result.light == 1 & result.gratingInfo.Contrast == 1),2);
    
    figure
    [s,p] = ttest(rl0c0,rl0c1);
    plot(rl0c0,rl0c1,'.')
    line([-.1,.7],[-.1,.7],'color','k')
    axis([-.1,.7,-.1,.7])
    axis square
    xlabel('average correlation zero contrast (light OFF)')
    ylabel('average correlation full contrast (light OFF)')
    title(['correlation changes with contrast light OFF: p: ' num2str(p)])
    
    figure
    [s,p] = ttest(rl1c0,rl1c1);
    plot(rl1c0,rl1c1,'.')
    line([-.1,.7],[-.1,.7],'color','k')
    axis([-.1,.7,-.1,.7])
    axis square
    xlabel('average correlation zero contrast (light ON)')
    ylabel('average correlation full contrast (light ON)')
    title(['correlation changes with contrast light ON: p: ' num2str(p)])
    
    figure
    [s,p] = ttest(rl0c0,rl1c0);
    plot(rl0c0,rl1c0,'.')
    line([-.1,.7],[-.1,.7],'color','k')
    axis([-.1,.7,-.1,.7])
    axis square
    xlabel('average correlation zero contrast (light OFF)')
    ylabel('average correlation zero contrast (light ON)')
    title(['correlation changes with light, zero contrast: p: ' num2str(p)])
    
    figure
    [s,p] = ttest(rl0c1,rl1c1);
    plot(rl0c1,rl1c1,'.')
    line([-.1,.7],[-.1,.7],'color','k')
    axis([-.1,.7,-.1,.7])
    axis square
    xlabel('average correlation full contrast (light OFF)')
    ylabel('average correlation full contrast (light ON)')
    title(['correlation changes with light, full contrast: p: ' num2str(p)])
    
    [s,p] = ttest(rl0c1-rl0c0,rl1c1-rl1c0)
    
end


figure
[s,pvrs] = ttest(nlosi(prs),losi(prs));
[s,pvfs] = ttest(nlosi(pfs),losi(pfs));
plot(nlosi(prs),losi(prs),'bo','linewidth',2)
hold on
plot(nlosi(pfs),losi(pfs),'ro','linewidth',2)
line([0,1],[0,1],'color','k')
axis square
xlabel('OSI control');
ylabel('OSI light');
title(['RS p: ' num2str(pvrs) ' FS p: ' num2str(pvfs)])

figure
[s,pvrs] = ttest(nldsi(prs),ldsi(prs));
[s,pvfs] = ttest(nldsi(pfs),ldsi(pfs));
plot(nldsi(prs),ldsi(prs),'bo','linewidth',2)
hold on
plot(nldsi(pfs),ldsi(pfs),'ro','linewidth',2)
line([0,1],[0,1],'color','k')
axis square
xlabel('DSI control');
ylabel('DSI light');
title(['RS p: ' num2str(pvrs) ' FS p: ' num2str(pvfs)])

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

disp('')
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
