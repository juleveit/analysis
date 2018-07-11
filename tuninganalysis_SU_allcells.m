function tuninganalysis_SU_allcells

animalid = '151027';
block = 10;
lcol = 'r';

onlymod = 0;
printyn = 1;
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

prestim = 800;
poststim = 300;
respwin = 501:1500; % after stimulus onset
respwin = respwin+prestim;
freqbinwidth = 5;

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
    
    isi = diff(msStimes);
    bursts = legendy_new3(isi,3,1000,3,0,15); %(ISI, fac, sr, min_length_of_burst, local_length, surprise_cutoff)
    surp = [bursts.surprise];
    bursts(surp == 100) = []; % delete probably wrong bursts
    burstbegs = [bursts.begin];
    burstchan = zeros(1,length(result.lfp));
    burstchan(msStimes(burstbegs)) = 1;
    
    wvchan = find(var(result.waveforms) == max(var(result.waveforms)));
    
    % LFP GAMMA
%     lfp = lfp-mean(lfp(1,1:beg),2);
%     stdbl = std(lfp(:,1:beg)');
%     lfp = lfp./stdbl;
    
    sr = 1000;
    cm = [3,4,1,2]; % confusion matrix? take lfp from two electrodes away to not get too many spike related phase resets
    lfp = result.lfp(:,cm(wvchan))';
    nfft = 2^nextpow2(length(lfp));
    fax = sr/2*linspace(0,1,nfft/2+1);
    y = fft(lfp,nfft);
    lfpspectrum = abs(y(1:nfft/2+1));
    
    %find gamma peaks for this animal
    beta = [15,40];
    gamma = [50,70];
    nl = find(result.light == 0);
    for i = 1:length(nl)
        [pn(i,:),f] = pmtm(lfp(result.msstamps(nl(i)):result.msstamps(nl(i))+1000),3,[],1000);
    end
    b1 = find(f>beta(1),1); b2 = find(f>beta(2),1);
    g1 = find(f>gamma(1),1); g2 = find(f>gamma(2),1);
    bsig = nanmean(pn(:,b1:b2));
    gsig = nanmean(pn(:,g1:g2));
    if isempty(find(diff(bsig)>0)) % there is no clear beta peak
        bpi = round((b1+b2)/2);
    else
        peaks = find(diff(bsig)>0)+1;
        pvs = bsig(peaks);
        bpi = peaks(pvs == max(pvs));
        bpi = bpi+b1-1;
    end
    if isempty(find(diff(gsig)>0)) % there is no clear beta peak
        gpi = round((g1+g2)/2);
    else
        peaks = find(diff(gsig)>0)+1;
        pvs = gsig(peaks);
        gpi = peaks(pvs == max(pvs));
        gpi = gpi+g1-1;
    end
    
    gamma1 = eegfilt(lfp,sr,f(bpi)-2.5,f(bpi)+2.5);
    gamma2 = eegfilt(lfp,sr,f(gpi)-2.5,f(gpi)+2.5);
    h1 = hilbert(gamma1); gpow1 = abs(h1); gphas1 = angle(h1);
    h2 = hilbert(gamma2); gpow2 = abs(h2); gphas2 = angle(h2);

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
%         msstamps([304,321]) = []; % for 140703 block 5
%         result.msstamps = msstamps;
%         save([supath, files(fi).name],'result');
        pause;
    end
%     trialnfft = 2^nextpow2(800);
%     trialfax = sr/2*linspace(0,1,trialnfft/2+1);


    gamma1phases = zeros(length(msstamps),trialdur+poststim+prestim);
    gamma2phases = zeros(length(msstamps),trialdur+poststim+prestim);
    vst = zeros(12,length(msstamps),trialdur+poststim+prestim);
    clear lfpspect; clear lfpoffsspect; clear fax;
    for i = 1:length(msstamps)
        resp(i,:) = chan(msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim); 
        burstresp(i,:) = burstchan(msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim); 
        lfpresp(i,:) = result.lfp(msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim);
%         y = fft(lfpresp(i,1001:1800),trialnfft);
%         lfpspect(i,:) = abs(y(1:trialnfft/2+1));
        [lfpspect(i,:),trialfax] = pmtm(lfpresp(i,1001:1800),3,[],sr);
        [lfpoffspect(i,:), trialfax] = pmtm(lfpresp(i,1801:2600),3,[],sr);
        gammaresp(i,:) = gpow(msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim);
        gphaseresp(i,:) = gphas(msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim);
        
        gamma1resp(i,:) = gpow1(msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim);
        gamma2resp(i,:) = gpow2(msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim);
        gamma1phasresp(i,:) = gphas1(msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim);
        gamma2phasresp(i,:) = gphas2(msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim);
        gamma1phases(i,find(resp(i,:))) = gamma1phasresp(i,find(resp(i,:)));
        gamma2phases(i,find(resp(i,:))) = gamma2phasresp(i,find(resp(i,:)));
        for pb = 1:8
            vstfrg1(pb,i) = length(find(gamma1phases(i,respwin)> -pi+(pb-1)*pi/4 & gamma1phases(i,respwin)< -pi+pb*pi/4));
            vstfrg2(pb,i) = length(find(gamma2phases(i,respwin)> -pi+(pb-1)*pi/4 & gamma2phases(i,respwin)< -pi+pb*pi/4));
        end
       
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
    ppcl0(cell) = ppc(l0phases{cell});
    ppcl1(cell) = ppc(l1phases{cell});
    
    tmp1 = zeros(size(gamma1phasresp));
    tmp1(find(resp)) = gamma1phasresp(find(resp));
    l0g1phasemat = tmp1(find(result.light == 0),:);
    l1g1phasemat = tmp1(find(result.light == 1),:);
    l0g1phases{cell} = l0g1phasemat(find(l0g1phasemat));
    l1g1phases{cell} = l1g1phasemat(find(l1g1phasemat));
    
    rg1l0(cell) = circ_r(l0g1phases{cell});
    rg1l1(cell) = circ_r(l1g1phases{cell});
    cmeang1l0(cell) = circ_mean(l0g1phases{cell});
    cmeang1l1(cell) = circ_mean(l1g1phases{cell});
    ppcg1l0(cell) = ppc(l0g1phases{cell});
    ppcg1l1(cell) = ppc(l1g1phases{cell});
    
    tmp2 = zeros(size(gamma2phasresp));
    tmp2(find(resp)) = gamma2phasresp(find(resp));
    l0g2phasemat = tmp2(find(result.light == 0),:);
    l1g2phasemat = tmp2(find(result.light == 1),:);
    l0g2phases{cell} = l0g2phasemat(find(l0g2phasemat));
    l1g2phases{cell} = l1g2phasemat(find(l1g2phasemat));
    
    rg2l0(cell) = circ_r(l0g2phases{cell});
    rg2l1(cell) = circ_r(l1g2phases{cell});
    cmeang2l0(cell) = circ_mean(l0g2phases{cell});
    cmeang2l1(cell) = circ_mean(l1g2phases{cell});
    ppcg2l0(cell) = ppc(l0g2phases{cell});
    ppcg2l1(cell) = ppc(l1g2phases{cell});
    
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
    
    if printyn
        figSize = [30 21];
        set(gcf,'paperunits','centimeters','papersize',figSize,'paperposition',[0 0 figSize])
        if cell<10, printi = ['0', int2str(cell)]; else printi = int2str(cell); end
        print(gcf,[lfpprintpath ,  printi '_' files(fi).name '.pdf'], '-dpdf' );
    end
    
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
            thisinds = find(result.gratingInfo.Orientation == oris(ori) &...
                result.light == l-1);
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
            condlfpspect(l,ori,:) = nanmean(lfpspect(thisinds,:));
            
            condgp(l,ori) = nanmean(gp(thisinds));
            
            for pb = 1:8
                cvstfrg1(pb,l,ori) = mean(vstfrg1(pb,thisinds),2);
                cvstfrg2(pb,l,ori) = mean(vstfrg2(pb,thisinds),2);
            end
        end
    end
    bincondresp = bincondresp.*(1000/binwidth);
    bta = bta-prestim;
%     psthcondplot(bincondresp,binconderr,bta);

    l0isi = []; l1isi = [];
    for i = find(result.light == 0)
        l0isi = [l0isi, isiresp{i}];
    end
    for i = find(result.light == 1)
        l1isi = [l1isi, isiresp{i}];
    end
    
    contindsnl = find(result.gratingInfo.Contrast == 0 & result.light == 0);
    controlresp(1,:) = mean(resp(contindsnl,:),1);
    controlfr(cell,1) = mean(frs(contindsnl));
    controlerr(1) = std(frs(contindsnl))./sqrt(length(contindsnl));
    controlgp(cell,1) = mean(gp(contindsnl));
    contindsl = find(result.gratingInfo.Contrast == 0 & result.light == 1);
    controlresp(2,:) = mean(resp(contindsl,:),1);
    controlfr(cell,2) = mean(frs(contindsl));
    controlerr(2) = std(frs(contindsl))./sqrt(length(contindsl));
    controlgp(cell,2) = mean(gp(contindsl));
    cfr(cell,:,:,:) = condfr;
    
    cellcondlfpspect(cell,:,:,:,:) = condlfpspect;
    
    [binctrl0,bta] = binit(controlresp(1,:),binwidth); binctrl0 = binctrl0.*(1000/binwidth);
    [binctrl1,bta] = binit(controlresp(2,:),binwidth); binctrl1 = binctrl1.*(1000/binwidth);
    
    if size(condfr,2)>1
        [nlprefratio, nldirpreffratio, nlprefori, nlmeanori, nlosi(cell), nlmeandir(cell), nldsi(cell)] = getOSI(squeeze(condfr(1,:)),oris);
        if size(condfr,1)>1
            [lprefratio, ldirprefratio, lprefori, lmeanori, losi(cell), lmeandir(cell), ldsi(cell)] = getOSI(squeeze(condfr(2,:)),oris);
        else
            losi(cell) = NaN; dosi(cell) = NaN;
        end
    end
   
    if size(condfr,2)>1
        blscondfr = condfr-mean(bl);
        oneorifr = mean(reshape(blscondfr(1,:),8,2),2);
        prefori = find(oneorifr == max(oneorifr),1);
        [preffr(cell,:), prefdir] = max(condfr(1,:),[],2); prefdir = prefdir(1);
        ortho = mod(prefori+4,length(oris)/2); if ortho == 0, ortho = length(oris)/2; end
    end
    
    if size(condfr,2)>1
        [binprefl0,bta] = binit(squeeze(condresp(1,prefdir,:)),binwidth); binprefl0 = binprefl0.*(1000/binwidth);
        [binprefl1,bta] = binit(squeeze(condresp(2,prefdir,:)),binwidth); binprefl1 = binprefl1.*(1000/binwidth);
    else
        [binprefl0,bta] = binit(squeeze(condresp(1,1,:)),binwidth); binprefl0 = binprefl0.*(1000/binwidth);
        [binprefl1,bta] = binit(squeeze(condresp(2,1,:)),binwidth); binprefl1 = binprefl1.*(1000/binwidth);
    end
    
    [binavgl1,bta] = binit(mean(lightresp),binwidth); binavgl1 = binavgl1.*(1000/binwidth); % in Hz
    [binavgl0,bta] = binit(mean(nolightresp),binwidth); binavgl0 = binavgl0.*(1000/binwidth);
    binstdl1 = binit(std(lightresp)./sqrt(size(lightresp,1)),binwidth);
    binstdl0 = binit(std(nolightresp)./sqrt(size(nolightresp,1)),binwidth);
    ta = bta-prestim;

    
  
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
    errorbar(oris,squeeze(condfr(2,:)),squeeze(conderr(2,:)),'o-','color',lcol,'markersize',8,'linewidth',2)
    hold on
    errorbar(oris,squeeze(condfr(1,:)),squeeze(conderr(1,:)),'ko-','markersize',8,'linewidth',2)
    xlabel('shown orientation')
    ylabel('Firing rate [Hz]')
    set(gca,'xtick',oris)
    legend({'Light ON','Light OFF'})
    if size(condfr,2)>1
        title(['OSI: ' num2str(nlosi(cell)) ' light OSI: ' num2str(losi(cell)) '   ' printname])
    end
    
    lcresp = [controlfr(cell,2), squeeze(nanmean(condfr(2,:,:),2))'];
    nolcresp = [controlfr(cell,1), squeeze(nanmean(condfr(1,:,:),2))'];
    xlevels = [0,clevels];
    
    subplot(2,4,7)
    plot(spike,'linewidth',2)
    axis([0,40,-100,100])
    legend(['width: ' int2str(swidth(cell)) ' adiff: ' num2str(adiff(cell))])
       
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
    
    condvstfrg1(cell,:,:,:) = cvstfrg1;
    condvstfrg2(cell,:,:,:) = cvstfrg2;
    
%     figure
%     subplot(2,2,1)
%     imagesc(speed);
%     colorbar
%     title(['oktrials: ' int2str(length(oktrials)) '/' int2str(size(speed,1))])
%     xlabel('time [ms]')
%     ylabel('trial number')
%     
%     subplot(2,2,2)
%     errorbar(msta,mean(speed(find(~result.light),:)),std(speed(find(~result.light),:))./sqrt(length(find(~result.light))),'b')
%     hold on
%     errorbar(msta,mean(speed(find(result.light),:)),std(speed(find(result.light),:))./sqrt(length(find(result.light))),'r')
%     xlabel('time [ms]')
%     ylabel('average runspeed')
%     legend({'light off' 'light on'})
%     
%     subplot(2,2,3)
%     plot(mean(speed(:,respwin),2),frs,'.')
%     hold on
%     plot(mean(speed(find(result.light),respwin),2),frs(find(result.light)),'r.')
%     xlabel('average runspeed of trial')
%     ylabel('average firing rate of trial')
    
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

pfs = swidth<120;
prs = swidth>=120;


for i = 1:length(swidth)
    for pb = 1:8
        for l = 1:2
            [condoprefratio(i,pb,l),conddprefratio(i,pb,l),condprefori(i,pb,l),condmeanori(i,pb,l),condosig1(i,pb,l),condmeandir(i,pb,l),conddsig1(i,pb,l)] = getOSI(squeeze(condvstfrg1(i,pb,l,:))',oris);
            [condoprefratio(i,pb,l),conddprefratio(i,pb,l),condprefori(i,pb,l),condmeanori(i,pb,l),condosig2(i,pb,l),condmeandir(i,pb,l),conddsig2(i,pb,l)] = getOSI(squeeze(condvstfrg2(i,pb,l,:))',oris);
        end
    end
end

rpax = -pi+pi/8:pi/4:pi-pi/8;

disp('');
