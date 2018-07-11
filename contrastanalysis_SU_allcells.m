function contrastanalysis_SU_allcells

animalid = '180626';
block = 12;
lcol = 'r';

onlymod = 0;
printyn = 0;
sfc = 0;


basepath = ['C:\Users\Julia\work\data\' animalid '\'];
supath = [basepath 'singleunits\'];
basename = [animalid '_block' int2str(block) '_tet'];

if printyn
    if ~exist([basepath, 'pdfs/'],'dir'),mkdir([basepath, 'pdfs/']),end
    if ~exist([basepath, 'pdfs/' 'LFPs/'],'dir'),mkdir([basepath, 'pdfs/' 'LFPs/']),end
    if ~exist([basepath, 'pdfs/' 'spikes/'],'dir'),mkdir([basepath, 'pdfs/' 'spikes/']),end
    if ~exist([basepath, 'pdfs/' 'running/'],'dir'),mkdir([basepath, 'pdfs/' 'running/']),end
end

lfpprintpath = [basepath 'pdfs\LFPs\'];
spikeprintpath = [basepath 'pdfs\spikes\'];
runprintpath = [basepath 'pdfs\running\'];
popprintpath = [basepath 'pdfs\'];

files = dir([supath, basename, '*.mat']);

freqbinwidth = 5;

cell = 1;
for fi = 1:length(files)
    
    %     if strfind(files(fi).name, 'MU')
    %         continue;
    %     end
    
    load([supath, files(fi).name]);
    
    prestim = 300;
    poststim = 300;
    if result.stimduration == 2
        respwin = 501:1500; % after stimulus onset
    else
        respwin = 1:1000;
    end
    respwin = respwin+prestim;
    
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
        %         msstamps(358) = []; % for 150603 block 4
        %         result.msstamps = msstamps;
        %         save([supath, files(fi).name],'result');
        %         msstamps([83]) = []; % for 160125 block 7
        %         msstamps([6, 260]) = []; % for 160219 block 5
%         msstamps([89, 1010, 1142, 1177]) = []; % for 160902 block 4
%         result.msstamps = msstamps;
%         save([supath, files(fi).name],'result');
                pause;
    end
    
    % fix so it is usable with new multi-purpose grating stim script
    if isfield(result, 'contconds')
        allinds = sort(getSpecificIndices(result, 'contconds'));
        msstamps = result.msstamps(allinds);
        light = result.light(allinds);
        gratingInfo.Orientation = result.gratingInfo.Orientation(allinds);
        gratingInfo.Contrast = result.gratingInfo.Contrast(allinds);
        gratingInfo.tFreq = result.gratingInfo.tFreq(allinds);
        gratingInfo.size = result.gratingInfo.size(allinds);
    else
        msstamps = result.msstamps;
        light = result.light;
        gratingInfo = result.gratingInfo;
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
        [lfpspect(i,:),trialfax] = pmtm(lfpresp(i,respwin(1)+200:respwin(end)),3,[],sr);
        
        %         [lfpoffspect(i,:), trialfax] = pmtm(lfpresp(i,1801:2600),3,[],sr);
        gammaresp(i,:) = gpow(msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim);
        gphaseresp(i,:) = gphas(msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim);
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
    lightmod(fi) = signrank(frs(find(light)),frs(find(light == 0)));
    
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
    l0phasemat = tmp(find(light == 0),:);
    l1phasemat = tmp(find(light == 1),:);
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
        alll0phasemat = tmpi(:,find(light == 0),:);
        alll1phasemat = tmpi(:,find(light == 1),:);
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
    semilogy(trialfax,nanmean(lfpspect(find(~light),:)),'linewidth',2);
    hold on
    semilogy(trialfax,nanmean(lfpspect(find(light),:)),'r','linewidth',2);
%     semilogy(trialfax,nanmean(lfpoffspect(find(light),:)),'c','linewidth',2);
    semilogy(trialfax,nanmean(lfpspect(find(~light),:))-(std(lfpspect(find(~light),:))./sqrt(length(find(~light)))));
    semilogy(trialfax,nanmean(lfpspect(find(~light),:))+(std(lfpspect(find(~light),:))./sqrt(length(find(~light)))));
    semilogy(trialfax,nanmean(lfpspect(find(light),:))-(std(lfpspect(find(light),:))./sqrt(length(find(light)))),'r');
    semilogy(trialfax,nanmean(lfpspect(find(light),:))+(std(lfpspect(find(light),:))./sqrt(length(find(light)))),'r');
%     semilogy(trialfax,nanmean(lfpoffspect(find(light),:))-(std(lfpoffspect(find(light),:))./sqrt(length(find(light)))),'c');
%     semilogy(trialfax,nanmean(lfpoffspect(find(light),:))+(std(lfpoffspect(find(light),:))./sqrt(length(find(light)))),'c');
    axis([0,120,...
        min([min(squeeze(nanmean(lfpspect(find(light == 0),1:125)))),min(squeeze(nanmean(lfpspect(find(light == 1),1:125))))]),...
        max([max(squeeze(nanmean(lfpspect(find(light == 0),1:125)))),max(squeeze(nanmean(lfpspect(find(light == 1),1:125))))])])
    legend({'light off', 'light on'});
    xlabel('frequency [Hz]')
    ylabel('spectral power')
    title('LFP spectrum during light on vs off')
    
    subplot(2,2,2)
    plot(nanmean(gammaresp(find(light == 0),:)))
    hold on
    plot(nanmean(gammaresp(find(light),:)),'r')
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
    plot(nanmean(lfpresp(find(light == 0),:)));
    hold on
    plot(nanmean(lfpresp(find(light),:)),'r')
    
    if printyn
        figSize = [30 21];
        set(gcf,'paperunits','centimeters','papersize',figSize,'paperposition',[0 0 figSize])
        if cell<10, printi = ['0', int2str(cell)]; else printi = int2str(cell); end
        print(gcf,[lfpprintpath ,  printi '_' files(fi).name '.pdf'], '-dpdf' );
    end
    
    msta = linspace(-prestim,trialdur+poststim,size(resp,2));
    
    lightresp = resp(find(light),:);
    nolightresp = resp(find(light == 0),:);
    lightlfpresp = lfpresp(find(light),:);
    nolightlfpresp = lfpresp(find(~light),:);
    
    %     frs = sum(resp(:,respwin),2)./(length(respwin)/1000);
    %     bl = sum(resp(:,1:prestim),2)./(prestim/1000);
    %
    %     %determine if cell is visually modulated
    %     blfr = sum(resp(:,1:prestim),2);
    %     vrfr = sum(resp(:,prestim+40:2*prestim+40),2);
    %     vismod(cell) = ttest2(blfr,vrfr);
    %
    %     %determine if cell is modulated by light
    %     lightmod(cell) = ttest2(frs(find(light)),frs(find(light == 0)));
    
    printname = files(fi).name;
    printname(find(printname=='_')) = ' ';
    
    nolbl = nanmean(bl(find(~light)));
    nolblerr = std(bl(find(~light)))./(sqrt(length(find(~light))));
    lbl = nanmean(bl(find(light)));
    lblerr = std(bl(find(light)))./(sqrt(length(find(light))));
    
    binwidth = 50;
    oris = unique(gratingInfo.Orientation); oris(find(oris == -1)) = [];
    clevels = unique(gratingInfo.Contrast); clevels(find(clevels == 0)) = []; %delete control condition
    for l = 1:length(unique(light))
        for ori = 1:length(oris)
            for co = 1:length(clevels)
                thisinds = find(gratingInfo.Orientation == oris(ori) &...
                    gratingInfo.Contrast == clevels(co) & ...
                    light == l-1);
                condn(l,ori,co) = length(thisinds);
                condresp(l,ori,co,:) = nanmean(resp(thisinds,:),1);
                condresperr(l,ori,co,:) = nanstd(resp(thisinds,:),1,1)./sqrt(length(thisinds));
                if ~isnan(condresp(l,ori,co,:))
                    [bincondresp(l,ori,co,:),bta] = binit(condresp(l,ori,co,:),binwidth);
                else
                    bincondresp(l,ori,co,:) = binit(condresp(l,ori,co,:),binwidth);
                end
                binconderr(l,ori,co,:) = binit(condresperr(l,ori,co,:),binwidth);
                condfr(l,ori,co) = nanmean(frs(thisinds));
                conderr(l,ori,co) =nanstd(frs(thisinds))./sqrt(length(thisinds));
                
                condgp(l,ori,co) = nanmean(gp(thisinds));
                condlfpspect(l,ori,co,:) = nanmean(lfpspect(thisinds,:),1);
                condlfpspecterr(l,ori,co,:) = nanstd(lfpspect(thisinds,:))./sqrt(length(thisinds));
            end
        end
    end
    bincondresp = bincondresp.*(1000/binwidth);
    bta = bta-prestim;
    %     psthcondplot(bincondresp,binconderr,bta);
    
    l0isi = []; l1isi = [];
    for i = find(light == 0)
        l0isi = [l0isi, isiresp{i}];
    end
    for i = find(light == 1)
        l1isi = [l1isi, isiresp{i}];
    end
    
    contindsnl = find(gratingInfo.Contrast == 0 & light == 0);
    controlresp(1,:) = mean(resp(contindsnl,:),1);
    controlfr(cell,1) = mean(frs(contindsnl));
    controlerr(1) = std(frs(contindsnl))./sqrt(length(contindsnl));
    controlgp(cell,1) = mean(gp(contindsnl));
    controllfpspect(cell,1,:) = nanmean(lfpspect(contindsnl,:),1);
    
    contindsl = find(gratingInfo.Contrast == 0 & light == 1);
    controlresp(2,:) = mean(resp(contindsl,:),1);
    controlfr(cell,2) = mean(frs(contindsl));
    controlerr(2) = std(frs(contindsl))./sqrt(length(contindsl));
    controlgp(cell,2) = mean(gp(contindsl));
    controllfpspect(cell,2,:) = nanmean(lfpspect(contindsl,:),1);
    cfr(cell,:,:,:) = condfr;
    
    cellcondlfpspect(cell,:,:,:,:) = condlfpspect;
    
    [binctrl0,bta] = binit(controlresp(1,:),binwidth); binctrl0 = binctrl0.*(1000/binwidth);
    [binctrl1,bta] = binit(controlresp(2,:),binwidth); binctrl1 = binctrl1.*(1000/binwidth);
    
    fc = find(clevels == 1); if isempty(fc), fc = find(clevels == max(clevels)); end
    pc = find(squeeze(nanmedian(condfr(1,:,:),2)) == max(squeeze(nanmedian(condfr(1,:,:),2))),1,'last');
    if size(condfr,2)>2
        [nlprefratio, nldirpreffratio, nlprefori, nlmeanori, nlosi(cell), nlmeandir(cell), nldsi(cell)] = getOSI(squeeze(condfr(1,:,fc)),oris);
    else
        nlosi(cell) = NaN; nldsi(cell) = NaN;
    end
    if size(condfr,1)>1&size(condfr,2)>2
        [lprefratio, ldirprefratio, lprefori, lmeanori, losi(cell), lmeandir(cell), ldsi(cell)] = getOSI(squeeze(condfr(2,:,fc)),oris);
    else
        losi(cell) = NaN; dosi(cell) = NaN;
    end
    
    blscondfr = condfr-mean(bl);
    oneorifr = mean(reshape(blscondfr(:,:,pc),size(blscondfr,1),length(oris)/2,2),3);
    prefori = find(oneorifr(1,:) == max(oneorifr(1,:)),1);
    [preffr(cell,:), prefdir] = max(condfr(:,:,pc),[],2); prefdir = prefdir(1);
    ortho = mod(prefori+2,length(oris)/2); if ortho == 0, ortho = length(oris)/2; end
    
    [binprefl0,bta] = binit(squeeze(condresp(1,prefdir,pc,:)),binwidth); binprefl0 = binprefl0.*(1000/binwidth);
    [binprefl1,bta] = binit(squeeze(condresp(2,prefdir,pc,:)),binwidth); binprefl1 = binprefl1.*(1000/binwidth);
    
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
    tempfreq = unique(gratingInfo.tFreq);
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
    errorbar(oris,squeeze(condfr(2,:,fc)),squeeze(conderr(2,:,fc)),'o-','color',lcol,'markersize',8,'linewidth',2)
    hold on
    errorbar(oris,squeeze(condfr(1,:,fc)),squeeze(conderr(1,:,fc)),'ko-','markersize',8,'linewidth',2)
    xlabel('shown orientation')
    ylabel('Firing rate [Hz]')
    set(gca,'xtick',oris)
    legend({'Light ON','Light OFF'})
    title(['OSI: ' num2str(nlosi(cell)) ' light OSI: ' num2str(losi(cell)) '   ' printname])
    
    oneorifr = nanmean(reshape(condfr(:,:,pc),size(condfr,1),length(oris)/2,2),3);
    prefori = find(oneorifr(1,:) == max(oneorifr(1,:)),1);
    ortho = mod(prefori+2,length(oris)/2); if ortho == 0, ortho = length(oris)/2; end
    
    lcresp = [controlfr(cell,2), squeeze(nanmean(condfr(2,:,:),2))'];
    nolcresp = [controlfr(cell,1), squeeze(nanmean(condfr(1,:,:),2))'];
    xlevels = [0,clevels];
    
    subplot(2,2,3)
    errorbar(xlevels,lcresp,[controlerr(2), squeeze(mean(conderr(2,:,:),2))'],'o','color',lcol,'markersize',8)
    hold on
    errorbar(xlevels,nolcresp,[controlerr(1), squeeze(mean(conderr(1,:,:),2))'],'ko','markersize',8)
    xlabel('shown contrast')
    ylabel('Firing rate [Hz]')
    legend({'Light ON','Light OFF'})
    set(gca,'xtick',clevels); %[0, clevels])
    ax = axis;
    axis([-.1,1.1,-.1,ax(4)])
    title(['contrast response all orientations ' printname]);
    
    [nlparams, bu, nlrsq] = fit_crf_NR(xlevels,nolcresp);
    [lparams, bu, lrsq] = fit_crf_NR(xlevels,lcresp);
    plotx = [0:0.01:1];
    plot(plotx,NakaRushton(nlparams,plotx),'k','linewidth',2);
    plot(plotx,NakaRushton(lparams,plotx),'r','linewidth',2);
    
    nlrmax(cell) = nlparams(1); lrmax(cell) = lparams(1);
    nlc50(cell) = nlparams(3); lc50(cell) = lparams(3);
    nlr0(cell) = nlparams(4); lr0(cell) = lparams(4);
    
    subplot(2,4,7)
    plot(spike,'linewidth',2)
    axis([0,40,-100,100])
    legend(['width: ' int2str(swidth(cell)) ' adiff: ' num2str(adiff(cell))])
    
    subplot(2,4,8)
    plot(ta,binctrl0,'linewidth',2)
    hold on
    plot(ta,binctrl1,lcol,'linewidth',2)
    axis([0,2500,min([min(binctrl0),min(binctrl1),0]),max([max(binctrl0),max(binctrl1),0.1])])
    
    if printyn
        figSize = [30 21];
        set(gcf,'paperunits','centimeters','papersize',figSize,'paperposition',[0 0 figSize])
        if cell<10, printi = ['0', int2str(cell)]; else printi = int2str(cell); end
        print(gcf,[spikeprintpath ,  printi '_' files(fi).name '.pdf'], '-dpdf' );
    end
    
    %     o = 1; c = 3;
    %     condl0 = (gratingInfo.Orientation == oris(o) | gratingInfo.Orientation == oris(o+length(oris)/2)) & (gratingInfo.Contrast == clevels(c) | gratingInfo.Contrast == clevels(c+1)) & light == 0;
    %     condl1 = (gratingInfo.Orientation == oris(o) | gratingInfo.Orientation == oris(o+length(oris)/2)) & (gratingInfo.Contrast == clevels(c) | gratingInfo.Contrast == clevels(c+1)) & light == 1;
    %     condl0 = gratingInfo.Orientation == oris(o) & gratingInfo.size == sizes(s) & light == 0;
    %     condl1 = gratingInfo.Orientation == oris(o) & gratingInfo.size == sizes(s) & light == 1;
    %     condl0 = light == 0; condl1 = light == 1;
    %     figure, rasterplot(resp,condl0,condl1,msta);
    line([0,2000],[37,37],'color','k','linewidth',2)
    line([500,1500],[33,33],'color','r','linewidth',2)
    
    preffr(cell,:) = condfr(:,prefori);
    
    % running figure
    lfr = mean(frs(find(light)));
    nlfr = mean(frs(find(~light)));
    runlfr = mean(frs(intersect(find(light),oktrials)));
    runnlfr = mean(frs(intersect(find(~light),oktrials)));
    norunlfr = mean(frs(intersect(find(light),stilltrials)));
    norunnlfr = mean(frs(intersect(find(~light),stilltrials)));
    lfrerr = std(frs(find(light)))./sqrt(length(find(light)));
    nlfrerr = std(frs(find(~light)))./sqrt(length(find(~light)));
    runlfrerr = std(frs(intersect(find(light),oktrials)))./sqrt(length(intersect(find(light),oktrials)));
    runnlfrerr = std(frs(intersect(find(~light),oktrials)))./sqrt(length(intersect(find(~light),oktrials)));
    norunlfrerr = std(frs(intersect(find(light),stilltrials)))./sqrt(length(intersect(find(light),stilltrials)));
    norunnlfrerr = std(frs(intersect(find(~light),stilltrials)))./sqrt(length(intersect(find(~light),stilltrials)));
    
    l1r1 = frs(intersect(find(light),oktrials));
    l0r1 = frs(intersect(find(~light),oktrials));
    l1r0 = frs(intersect(find(light),stilltrials));
    l0r0 = frs(intersect(find(~light),stilltrials));
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
    errorbar(msta,mean(speed(find(~light),:)),std(speed(find(~light),:))./sqrt(length(find(~light))),'b')
    hold on
    errorbar(msta,mean(speed(find(light),:)),std(speed(find(light),:))./sqrt(length(find(light))),'r')
    xlabel('time [ms]')
    ylabel('average runspeed')
    legend({'light off' 'light on'})
    
    subplot(2,2,3)
    plot(mean(speed(:,respwin),2),frs,'.')
    hold on
    plot(mean(speed(find(light),respwin),2),frs(find(light)),'r.')
    xlabel('average runspeed of trial')
    ylabel('average firing rate of trial')
    
    %     subplot(2,2,4)
    %     barweb([nlfr,lfr;runnlfr,runlfr;norunnlfr,norunlfr],...
    %         [nlfrerr,lfrerr;runnlfrerr,runlfrerr;norunnlfrerr,norunlfrerr],...
    %         [],[{'all'};{'running only'};{'immobile only'}],['ANOVA factor running p: ' num2str(p(2))],...
    %         [],'firing rate [Hz]',[],[]);
    
    if printyn
        figSize = [30 21];
        set(gcf,'paperunits','centimeters','papersize',figSize,'paperposition',[0 0 figSize])
        if cell<10, printi = ['0', int2str(cell)]; else printi = int2str(cell); end
        print(gcf,[runprintpath ,  printi '_' files(fi).name '.pdf'], '-dpdf' );
    end
    
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
cell1 = 31; cell2 = 10;
for i = 1:size(celllfpresp,2)
    [coh(i,:),cfx] = mscohere(squeeze(celllfpresp(cell1,i,respwin)),squeeze(celllfpresp(cell2,i,respwin)),[],[],512,1000);
end
for l = 1:2
    for ori = 1:length(oris)
        for co = 1:length(clevels)
            thisinds = find(gratingInfo.Orientation == oris(ori) &...0
                gratingInfo.Contrast ==clevels(co) & ...
                light == l-1);
            condcoher(l,ori,co,:) = nanmean(coh(thisinds,:),1);
        end
    end
end
contcoher(1,:) = nanmean(coh(contindsnl,:),1);
contcoher(2,:) = nanmean(coh(contindsl,:),1);



% get STA and SFC spike rate equalized
if sfc
    for cell = 1:size(cellresp,1)
        
        respl0 = squeeze(cellresp(cell,find(light == 0),respwin));
        respl1 = squeeze(cellresp(cell,find(light == 1),respwin));
        lfprespl0 = squeeze(celllfpresp(cell,find(light == 0),respwin));
        lfprespl1 = squeeze(celllfpresp(cell,find(light == 1),respwin));
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
            binresp(i,j,:) = binit(squeeze(irsresp(i,j,:)),100);
        end
    end
    
    %method 1
    for i = 1:420
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
    
    rl0 = rval(:,find(~light));
    rl1 = rval(:,find(light));
    
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
    
    rl0c0 = nanmean(rval(:,light == 0 & gratingInfo.Contrast == 0),2);
    rl0c1 = nanmean(rval(:,light == 0 & gratingInfo.Contrast == 1),2);
    rl1c0 = nanmean(rval(:,light == 1 & gratingInfo.Contrast == 0),2);
    rl1c1 = nanmean(rval(:,light == 1 & gratingInfo.Contrast == 1),2);
    
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
plot(nlc50(prs),depth(prs),'b.','linewidth',2)
hold on
plot(nlc50(pfs),depth(pfs),'r.','linewidth',2)
axis ij
axis([-.1,1.1,150,1050])
title('C50 in depth')
xlabel('C50')
ylabel('cortical depth [mum]')

figure
plot((lc50(prs)-nlc50(prs))./(lc50(prs)+nlc50(prs)),depth(prs),'bo','linewidth',2)
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
plot((lc50(pfs)-nlc50(pfs))./(lc50(pfs)+nlc50(pfs)),depth(pfs),'ro','linewidth',2)
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
plot(nlc50(prs),lc50(prs),'bo','linewidth',2)
hold on
plot(nlc50(pfs),lc50(pfs),'ro','linewidth',2)
axis([-.1,1.1,-.1,1.1])
line([-.1,1.1],[-.1,1.1],'color','k')
axis square
[s,pvrs] = ttest(nlc50(prs),lc50(prs));
[s,pvfs] = ttest(nlc50(pfs),lc50(pfs));
xlabel('c50 no light')
ylabel('c50 light on')
title(['C50 changes with light. p RS: ' num2str(pvrs) '  p FS: ' num2str(pvfs)])
if printyn
    figSize = [30 21];
    set(gcf,'paperunits','centimeters','papersize',figSize,'paperposition',[0 0 figSize])
    print(gcf,[popprintpath ,  '09_C50scatter.pdf'], '-dpdf' );
end

nlr0(find(nlr0<0)) = 0;
lr0(find(lr0<0)) = 0;

figure
plot(nlr0(prs),depth(prs),'bo','linewidth',2)
hold on
plot(nlr0(pfs),depth(pfs),'ro','linewidth',2)
axis ij
axis([-.1,30,150,1050])
title('r0 in depth')
xlabel('r0')
ylabel('cortical depth [mum]')

figure
plot((lr0(prs)-nlr0(prs))./(lr0(prs)+nlr0(prs)),depth(prs),'bo','linewidth',2)
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
plot((lr0(pfs)-nlr0(pfs))./(lr0(pfs)+nlr0(pfs)),depth(pfs),'ro','linewidth',2)
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
plot(nlr0(prs),lr0(prs),'bo','linewidth',2)
hold on
plot(nlr0(pfs),lr0(pfs),'ro','linewidth',2)
line([0,30],[0,30],'color','k')
axis square
[s,pvrs] = ttest(nlr0(prs),lr0(prs));
[s,pvfs] = ttest(nlr0(pfs),lr0(pfs));
xlabel('r0 no light')
ylabel('r0 light on')
title(['r0 changes with light. p RS: ' num2str(pvrs) '  p FS: ' num2str(pvfs)])
if printyn
    figSize = [30 21];
    set(gcf,'paperunits','centimeters','papersize',figSize,'paperposition',[0 0 figSize])
    print(gcf,[popprintpath ,  '12_R0scatter.pdf'], '-dpdf' );
end

figure
plot(nlrmax,depth,'bo','linewidth',2)
hold on
plot(nlrmax(pfs),depth(pfs),'ro','linewidth',2)
axis ij
ax = axis;
axis([ax(1),ax(2),150,1050])
title('rmax in depth')
xlabel('rmax')
ylabel('cortical depth [mum]')

figure
plot((lrmax(prs)-nlrmax(prs))./(lrmax(prs)+nlrmax(prs)),depth(prs),'bo','linewidth',2)
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
plot((lrmax(pfs)-nlrmax(pfs))./(lrmax(pfs)+nlrmax(pfs)),depth(pfs),'ro','linewidth',2)
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
plot(nlrmax(prs),lrmax(prs),'bo','linewidth',2)
hold on
plot(nlrmax(pfs),lrmax(pfs),'ro','linewidth',2)
line([0,max([max(nlrmax),max(lrmax)])+2],[0,max([max(nlrmax),max(lrmax)])+2],'color','k')
axis square
% axis([0,80,0,80])
[s,pvrs] = ttest(nlrmax(prs),lrmax(prs));
[s,pvfs] = ttest(nlrmax(pfs),lrmax(pfs));
xlabel('rmax no light')
ylabel('rmax light on')
title(['rmax changes with light. p RS: ' num2str(pvrs) ' p FS: ' num2str(pvfs)])
if printyn
    figSize = [30 21];
    set(gcf,'paperunits','centimeters','papersize',figSize,'paperposition',[0 0 figSize])
    print(gcf,[popprintpath ,  '15_Rmaxscatter.pdf'], '-dpdf' );
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
if printyn
    figSize = [30 21];
    set(gcf,'paperunits','centimeters','papersize',figSize,'paperposition',[0 0 figSize])
    print(gcf,[popprintpath ,  '16_OSIscatter.pdf'], '-dpdf' );
end

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

function [params, resnorm, rsquared] = fit_crf_NR(x,y)
% [rmax, n, c50, r0]

nrfitmargin = .01; % .5
range = max(y)-min(y);
if range ~= 0
    p0 = [range 2 .5 min(y)];
    lb = [(1-nrfitmargin)*range, .5, 0, min(y)-.1*range];
    ub = [(1+nrfitmargin)*range, 10, 1, min(y)+.1*range];
    %         lb = [-max(y), .5, 0, -10*max(y)];
    %         ub = [10*max(y), 10, 1, max(y)];
    warning off
    [params,resnorm,residual,exitflag] = lsqcurvefit(@(p,x) NakaRushton(p,x),p0,x(:),y(:),lb,ub,optimset('Display','off'));
    warning on
    restot = sum((y-mean(y)).^2);
    rsquared = 1 - (resnorm/restot);
else
    params = [NaN,NaN,NaN,NaN];
    resnorm = NaN;
    residual = NaN;
    exitflag = NaN;
    rsquared = NaN;
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
