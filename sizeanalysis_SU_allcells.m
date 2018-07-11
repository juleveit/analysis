function sizeanalysis_SU_allcells

animalid = '171205';
block = 7;
lcol = 'r'; %lasercolor

onlymod = 0;
printyn = 1;
sfc = 0;

addpath C:\Users\Julia\work\Matlab\others\fieldtrip-20160329
ft_defaults

supath = ['C:\Users\Julia\work\data\' animalid '\singleunits\'];
% supath = ['C:\Users\Julia\work\data\' animalid '\multiunits\'];
basename = [animalid '_block' int2str(block) '_tet'];

files = dir([supath, basename, '*.mat']);

freqbinwidth = 5;

% chronux parameters
params.tapers = [5,9]; params.Fs = 1000; params.err = [2, 0.05]; params.trialave = 1;

cll = 1;
for fi = 1:length(files)
%     
%     if strfind(files(fi).name, 'MU')
%         continue;
%     end
        
    load([supath, files(fi).name]);    
    
    prestim = 300;
    poststim = 700;
    if result.stimduration == 2
        respwin = 501:1500; % after stimulus onset
    else
        respwin = 1:1000;
    end
    respwin = respwin+prestim;
    
%     disp(['now analyzing file: ' files(cll).name]);
    cllname{cll} = files(fi).name;
    
    i = strfind(files(fi).name, 'tet');
    tetno(cll) = strread(files(fi).name(i+3));
    
    wvchan = find(var(result.waveforms) == max(var(result.waveforms)));
    
    sr = 1000;
    lfp = result.lfp(:,wvchan)';
    nfft = 2^nextpow2(length(lfp));
    fax = sr/2*linspace(0,1,nfft/2+1);
    y = fft(lfp,nfft);
    lfpspectrum = abs(y(1:nfft/2+1));

%     plot(fax,lfpspectrum)
    gamma = eegfilt(lfp,sr,15,35);  % adjust depending on peak
    h = hilbert(gamma); gpow = abs(h); gphas = angle(h);
    % %
    if sfc
        for i = 1:100/freqbinwidth
            filtmat(i,:) = eegfilt(lfp,sr,(i-1)*freqbinwidth+1,i*freqbinwidth);
            h = hilbert(filtmat(i,:));
            powmat(i,:) = abs(h); phasmat(i,:) = angle(h);
        end
    end
    filtlfp = eegfilt(lfp,sr,10,100);
    
    msStimes = round(result.spikes);
    if isempty(msStimes), msStimes(1) = 0; end
   if msStimes(1) == 0, msStimes(1) = 1; end  
    
    chan = zeros(1,length(result.lfp));
    chan(msStimes) = 1;
    
    trialdur = result.stimduration*1000;
    msstamps = result.msstamps;
    
     if length(msstamps)~=length(result.light)
         disp('');
%         msstamps([62,108,147]) = []; % for 140703 block 8
%         msstamps([161]) = []; % for 141204 block 3
%         msstamps([303]) = []; % for 150407 block 5
%         msstamps([169,336]) = []; % for 150523 block 11
%         msstamps([24]) = []; % for 150730 block 11
%         msstamps([207]) = []; % for 150730 block 11
%         msstamps([242]) = []; % for 151210 block 4
%         msstamps([412]) = []; % for 151210 block 3
%         msstamps([83]) = []; % for 160125 block 7
%         msstamps([302,340]) = []; % for 160328 block 2
%         msstamps([518]) = []; % for 160726 block 6
%         msstamps([390,631,875]) = []; % for 170426 block 2
%         result.msstamps = msstamps;
%         save([supath, files(fi).name],'result');
        pause;
     end    
    
    % fix so it is usable with new multi-purpose grating stim script
    if isfield(result, 'sizeconds')
        allinds = sort(getSpecificIndices(result, 'sizeconds'));
        msstamps = result.msstamps(allinds);
        light = result.light(allinds);
        gratingInfo.Orientation = result.gratingInfo.Orientation(allinds);
        gratingInfo.size = result.gratingInfo.size(allinds);
        gratingInfo.Contrast = result.gratingInfo.Contrast(allinds);
        gratingInfo.tFreq = result.gratingInfo.tFreq(allinds);
    else
        msstamps = result.msstamps;
        light = result.light;
        gratingInfo = result.gratingInfo;
    end
    iii = 1;
    for i = 1:length(msstamps)
        resp(i,:) = chan(msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim); 
        hh = find(resp(i,1001:1800))'; ptresp(i).times = hh./1000;
        lfpresp(i,:) = result.lfp(msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim, wvchan);
        filtlfpresp(i,:) = filtlfp(msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim);
        [lfpspect(i,:),trialfax] = pmtm(lfpresp(i,respwin(1)+200:respwin(end)),3,[],sr);
        [lfpsectrogram(i,:,:),ct,cf] = mtspecgramc(lfpresp(i,:)',[.25,.1],params);
        gammaresp(i,:) = gamma(msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim);
        gammapowresp(i,:) = gpow(msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim);
        gphaseresp(i,:) = gphas(msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim);
        bpresp(i,:) = gamma(msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim);
        
        s = find(resp(i,respwin));
        if ~isempty(s)
            for si = 1:length(s)
                lfps(si,:) = lfpresp(i,s(si)+respwin(1)-1-200:s(si)+respwin(1)-1+200);
            end
        else
            lfps = nan(1,401);
        end
        trialstalfp(i,:) = mean(lfps,1);
        nspkstrialsta(i) = size(lfps,1);
        clear lfps;
        [spkxcorr(i,:),lags] = xcorr(resp(i,respwin),resp(i,respwin));
        
        if sfc
            for j = 1:size(phasmat,1)
                allphaseresp(j,i,:) = phasmat(j, msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim);
                allpowresp(j,i,:) = powmat(j, msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim);
            end
        end
        
        speed(i,:) = result.runspeed(msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim);  
        
        a = find(resp(i,:));
        for ii = 1:length(a)
            speedsta(iii,:) = result.runspeed(msstamps(i)-prestim+1+a(ii)-500:msstamps(i)-prestim+1+a(ii)+500);
            iii = iii +1;
        end
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
       
%     a = find(gratingInfo.size == max(unique(gratingInfo.size)) & light == 0);
% %     b = find(gratingInfo.size == max(unique(gratingInfo.size)) & light == 1);
%     b = find((gratingInfo.size == 52 | gratingInfo.size == 36 )& light == 1);
%     for i = 1:40
%         phasdiff{i} = find(diff(gphaseresp(a(i),300:2300))<-5);
%         clear gamps;
%         for j = 1:length(phasdiff{i})
%             gamps(j) = bpresp(i,phasdiff{i}(j)+1);
%         end
%         gampl{i} = gamps;
%         intervals{i} = diff(phasdiff{i});
%     end
%     allints = [];
%     for i = 1:40
%         allints = [allints,intervals{i}];
%     end
%     meancycle = median(allints);
%     for i = 1:length(b)        
%         stimphas(i) = gphaseresp(b(i),800); % find out phase of pulse
%         % get all differences of one average phase interval starting 1/3 of an average cycle after stim pulse
%         phasdiff = stimphas(i)-gphaseresp(i,800+round(meancycle/3):800+round(meancycle/3)+meancycle);
%         [xx,ind] = min(abs(phasdiff)); % find minimum difference
%         nextphas = 800+round(meancycle/3)+ind-1;  % get index of next occurence of stim phase
%         stimcyclelen(i) = nextphas-800;
%     end

    
    clear stalfp; clear nspkssta; clear shufstalfp; clear nspksshufsta;
    sizes = unique(gratingInfo.size); sizes(sizes == 0) = [];
    for l = 1:length(unique(light))
        for sz = 1:length(sizes)
            trials = find(gratingInfo.size == sizes(sz) & light == l-1);
            
            %             [C(cll,l,sz,:),phi(cll,l,sz,:),S12(cll,l,sz,:),S1(cll,l,sz,:),...
            %                 S2(cll,l,sz,:),chfx,zerosp(cll,l,sz,:),confC(cll,l,sz),...
            %                 phistd(cll,l,sz,:),Cerr(cll,l,sz,:,:)] = coherencycpt(squeeze(lfpresp(trials,1001:1800))',...
            %                 ptresp(trials),params);
            
            sizelfpresp(cll,l,sz,:,:) = lfpresp(trials,:);
            [S,chf,Serr]=mtspectrumc(squeeze(lfpresp(trials,1001:1800))',params);
            condS(cll,l,sz,:) = S(1:150); condSerr(cll,l,sz,:,:) = Serr(:,1:150);
            
            thisruninds = intersect(trials,oktrials); thisstillinds = intersect(trials,stilltrials);
            if ~isempty(thisruninds)
                thisrunn(cll,l,sz) = length(thisruninds);
                [runS,chf,runSerr] = mtspectrumc(squeeze(lfpresp(thisruninds,1001:1800))',params);
                condrunS(cll,l,sz,:) = runS(1:150); condrunSerr(cll,l,sz,:,:) = runSerr(:,1:150);                
            else
                thisrunn(cll,l,sz) = 0;
                condrunS(cll,l,sz,:) = nan(1,150); condrunSerr(cll,l,sz,:,:) = nan(2,150);
            end
            if ~isempty(thisstillinds)
                thisstilln(cll,l,sz) = length(thisstillinds);
                [stillS,chf,stillSerr] = mtspectrumc(squeeze(lfpresp(thisstillinds,1001:1800))',params);
                condstillS(cll,l,sz,:) = stillS(1:150); condstillSerr(cll,l,sz,:,:) = stillSerr(:,1:150);
            else
                thisstilln(cll,l,sz) = 0;
                condstillS(cll,l,sz,:) = nan(1,150); condstillSerr(cll,l,sz,:,:) = nan(2,150);
            end
            
            rp = randperm(length(trials));
            for i = 1:length(trials)
                s = find(resp(trials(i),respwin)); % find spikes from the actual trial
                if ~isempty(s)
                    for si = 1:length(s) % take LFPs from random other smae size trial
                        lfps(si,:) = lfpresp(trials(i),s(si)+respwin(1)-1-200:s(si)+respwin(1)-1+200);
                    end
                    for si = 1:length(s) % take LFPs from random other smae size trial
                        shuflfps(si,:) = lfpresp(trials(rp(i)),s(si)+respwin(1)-1-200:s(si)+respwin(1)-1+200);
                    end
                else
                    lfps = nan(1,401);
                    shuflfps = nan(1,401);
                end
                
                stalfp(l,sz,i,:) = mean(lfps,1);
                nspkssta(l,sz,i) = size(lfps,1);
                shufstalfp(l,sz,i,:) = mean(shuflfps,1);
                nspksshufsta(l,sz,i) = size(shuflfps,1);
                clear lfps; clear shuflfps;
            end
            
            
            % field trip messing around
            for i = 1:length(trials)
                data.trial{i} = [lfpresp(trials(i),:);resp(trials(i),:)];
                time = -prestim+1:result.stimduration*1000+poststim;
                data.time{i} = time(1)./1000:.001:time(end)./1000;
%                 data.time{i} = -.299:.001:2.700;
            end
            data.fsample = 1000;
            data.trialinfo(:,1) = light(trials);
            data.trialinfo(:,2) = gratingInfo.size(trials);
            data.trialinfo(:,3) = gratingInfo.Orientation(trials);
            cfg = [];
            cfg.timwin = [-.25 .25];
            data.label{1,1} = 'lfp';
            data.label{2,1} = 'spikes';
            cfg.spikechannel = 'spikes';
            cfg.channel = 'lfp';
            cfg.latency = [0.5,1.5];
            
            cfg.method = 'mtmfft';
            cfg.foilim = [5,100];
            cfg.timwin = [-.15, .15];
            cfg.taper = 'hanning';
            cfg.spikechannel = 'spikes';
            cfg.channel = 'lfp';
            stsFFT           = ft_spiketriggeredspectrum(cfg, data);
            ang = angle(stsFFT.fourierspctrm{1});
            mag = abs(stsFFT.fourierspctrm{1});
            ftspect(cll,l,sz,:) = squeeze(nanmean(mag(:,1,:)));
            ftphases{cll,l,sz} = squeeze(ang);
            ftfax = stsFFT.freq;
            
            cfg               = [];
            cfg.method        = 'ral'; % compute the rayleigh test
            cfg.spikechannel  = stsFFT.label{1};
            cfg.channel       = stsFFT.lfplabel; % selected LFP channels
            cfg.avgoverchan   = 'unweighted'; % weight spike-LFP phases irrespective of LFP power
            cfg.timwin        = 'all'; % compute over all available spikes in the window
            cfg.latency       = [0.5 1.5]; % sustained visual stimulation period
            statSts           = ft_spiketriggeredspectrum_stat(cfg,stsFFT);
            ftralp(cll,l,sz,:) = statSts.ral;
            cfg.method = 'ppc0';
            statSts           = ft_spiketriggeredspectrum_stat(cfg,stsFFT);
            ftppc(cll,l,sz,:) = statSts.ppc0;
            cfg.method = 'plv';
            statSts           = ft_spiketriggeredspectrum_stat(cfg,stsFFT);
            ftplv(cll,l,sz,:) = statSts.plv;
        end
    end
    sizestalfp(cll,:,:,:) = nanmean(stalfp,3);
    shufsizestalfp(cll,:,:,:) = nanmean(shufstalfp,3);
    sfx = chf(1:150);
    
    % gamma in time figure
%     s1l0 = find(gratingInfo.size == sizes(1) & light == 0);
%     s1l1 = find(gratingInfo.size == sizes(1) & light == 1);
%     s3l0 = find(gratingInfo.size == sizes(3) & light == 0);
%     s3l1 = find(gratingInfo.size == sizes(3) & light == 1);
%     s5l0 = find(gratingInfo.size == sizes(5) & light == 0);
%     s5l1 = find(gratingInfo.size == sizes(5) & light == 1);
%     xfill = [-299:2300,fliplr(-299:2300)];
%     yfs1l0 = [mean(gammapowresp(s1l0,1:2600))+std(gammapowresp(s1l0,1:2600))./sqrt(40),fliplr(mean(gammapowresp(s1l0,1:2600))-std(gammapowresp(s1l0,1:2600))./sqrt(40))];
%     yfs1l1 = [mean(gammapowresp(s1l1,1:2600))+std(gammapowresp(s1l1,1:2600))./sqrt(40),fliplr(mean(gammapowresp(s1l1,1:2600))-std(gammapowresp(s1l1,1:2600))./sqrt(40))];
%     yfs3l0 = [mean(gammapowresp(s3l0,1:2600))+std(gammapowresp(s3l0,1:2600))./sqrt(40),fliplr(mean(gammapowresp(s3l0,1:2600))-std(gammapowresp(s3l0,1:2600))./sqrt(40))];
%     yfs3l1 = [mean(gammapowresp(s3l1,1:2600))+std(gammapowresp(s3l1,1:2600))./sqrt(40),fliplr(mean(gammapowresp(s3l1,1:2600))-std(gammapowresp(s3l1,1:2600))./sqrt(40))];
%     yfs5l0 = [mean(gammapowresp(s5l0,1:2600))+std(gammapowresp(s5l0,1:2600))./sqrt(40),fliplr(mean(gammapowresp(s5l0,1:2600))-std(gammapowresp(s5l0,1:2600))./sqrt(40))];
%     yfs5l1 = [mean(gammapowresp(s5l1,1:2600))+std(gammapowresp(s5l1,1:2600))./sqrt(40),fliplr(mean(gammapowresp(s5l1,1:2600))-std(gammapowresp(s5l1,1:2600))./sqrt(40))];
%     
%     figure
%     fill(xfill,yfs1l1,[1,.8,.8])
%     hold on
%     fill(xfill,yfs3l1,[1,.5,.5])
%     fill(xfill,yfs5l1,[1,.3,.3])
%     axis([-300,2300,0,180])
%     
%     figure
%     fill(xfill,yfs1l0,[.8,.8,1])
%     hold on
%     fill(xfill,yfs3l0,[.5,.5,1])
%     fill(xfill,yfs5l0,[.3,.3,1])
%     axis([-300,2500,0,180])

    depth(cll) = result.depth;
    
    spike = result.waveforms(:,wvchan);
    interpspike = spline(1:32,spike,1:.1:32);
    [adiff(cll),swidth(cll)] = spikequant(interpspike);
    
    cllrespl0(cll,:,:) = resp(find(light == 0),:);
    cllrespl1(cll,:,:) = resp(find(light == 1),:);
    cllresp(cll,:,:) = resp;
    clllfpresp(cll,:,:) = lfpresp;
    
    
    msta = linspace(-prestim,trialdur+poststim,size(resp,2));
    
    lightresp = resp(find(light),:);
    nolightresp = resp(find(light == 0),:);
    
    lightlfpresp = lfpresp(find(light),:);
    nolightlfpresp = lfpresp(find(~light),:);
    
    frs = sum(resp(:,respwin),2)./(length(respwin)/1000);
    bl = sum(resp(:,1:prestim),2)./(prestim/1000);
    sc = sum(resp(:,respwin),2);
        
     %determine if cll is visually modulated
    blfr = sum(resp(:,1:prestim),2);
    vrfr = sum(resp(:,prestim+40:2*prestim+40),2);
    vismod(cll) = ttest2(blfr,vrfr);
    
    %determine if cll is modulated by light
    lightmod(cll) = ttest2(frs(find(light)),frs(find(light == 0)));
    
    lfr(cll) = mean(frs(find(light)));
    nlfr(cll) = mean(frs(find(light == 0)));
    
    % phases
    tmp = zeros(size(gphaseresp));
    tmp(find(resp)) = gphaseresp(find(resp));
    l0phasemat = tmp(find(light == 0),:);
    l1phasemat = tmp(find(light == 1),:);
    l0phases{cll} = l0phasemat(find(l0phasemat));
    l1phases{cll} = l1phasemat(find(l1phasemat));
      
    rl0(cll) = circ_r(l0phases{cll});
    rl1(cll) = circ_r(l1phases{cll});
    cmeanl0(cll) = circ_mean(l0phases{cll});
    cmeanl1(cll) = circ_mean(l1phases{cll});
    
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
            allrl0(cll,i) = circ_r(allphasesl0{i}');
            allrl1(cll,i) = circ_r(allphasesl1{i}');
            allcmeanl0(cll,i) = circ_mean(allphasesl0{i}');
            allcmeanl1(cll,i) = circ_mean(allphasesl1{i}');
        end
    end
    
    
    binwidth = 30;
    [binnedlight,bta] = binit(mean(lightresp),binwidth);
    binnedlight = binnedlight.*(1000/binwidth);
    [binnednolight,bta] = binit(mean(nolightresp),binwidth);
    binnednolight = binnednolight.*(1000/binwidth);
    
    printname = files(fi).name;
    printname(find(printname=='_')) = ' ';
        
    sizes = unique(gratingInfo.size);  sizes(find(sizes == 0)) = []; %delete control condition
    oris = unique(gratingInfo.Orientation); oris(find(oris == -1)) = [];
    for l = 1:length(unique(light))
        for ori = 1:length(oris)
            for sz = 1:length(sizes)
                thisinds = find(gratingInfo.Orientation == oris(ori) &...0
                    gratingInfo.size == sizes(sz) & ...
                    light == l-1);
                condresp(l,ori,sz,:) = mean(resp(thisinds,:),1);
                condlfpspect(l,ori,sz,:) = nanmean(lfpspect(thisinds,:));
                condstspect(l,ori,sz,:) = pmtm(mean(resp(thisinds,1001:1800),1),3,[],sr);
                condfr(l,ori,sz) = mean(frs(thisinds));%-mean(bl);        
                conderr(l,ori,sz) =std(frs(thisinds))./sqrt(length(thisinds));
                
                condz(l,ori,sz) = {(sc(thisinds)-mean(sc(thisinds)))/std(sc(thisinds))}; %ecker 2010
                condsc(l,ori,sz) = {sc(thisinds)};
                ff(l,ori,sz) = var(sc(thisinds))/mean(sc(thisinds));
                
                condlfpresp(cll,l,ori,sz,:) = mean(lfpresp(thisinds,:),1);                
                
                condresperr(l,ori,sz,:) = nanstd(resp(thisinds,:),1,1)./sqrt(length(thisinds));
                if ~isnan(condresp(l,ori,sz,:))
                    [bincondresp(l,ori,sz,:),bta] = binit(condresp(l,ori,sz,:),binwidth);
                else
                    bincondresp(l,ori,sz,:) = binit(condresp(l,ori,sz,:),binwidth);
                end
                binconderr(l,ori,sz,:) = binit(condresperr(l,ori,sz,:),binwidth);
                
                mscc = []; bincc = [];
                for ii = 1:length(thisinds)-1
                    for jj = ii+1:length(thisinds)
                        help = corrcoef(resp(thisinds(ii),:),resp(thisinds(jj),:));
                        mscc = [mscc,help(1,2)];
                        help = corrcoef(binit(resp(thisinds(ii),:),binwidth),binit(resp(thisinds(jj),:),binwidth));
                        bincc = [bincc, help(1,2)];
                    end
                end
                msreliab(l,ori,sz) = nanmean(mscc);
                binreliab(l,ori,sz) = nanmean(bincc);
                eckerreliability(l,ori,sz) = var(frs(thisinds))/var(frs);

                thisruninds = intersect(thisinds,oktrials);
                if ~isempty(thisruninds)
                    runcondresp(l,ori,sz,:) = mean(resp(thisruninds,:),1);
                    runcondfr(l,ori,sz) = mean(frs(thisruninds));
                    runconderr(l,ori,sz) = std(frs(thisruninds))./sqrt(length(thisruninds));
                else
                    runcondresp(l,ori,sz,:) = nan(1,size(resp,2));
                    runcondfr(l,ori,sz) = NaN;
                    runconderr(l,ori,sz) = NaN;
                end  
                
                thisstillinds = intersect(thisinds,stilltrials);
                if ~isempty(thisstillinds)
                    stillcondresp(l,ori,sz,:) = mean(resp(thisstillinds,:),1);
                    stillcondfr(l,ori,sz) = mean(frs(thisstillinds));
                    stillconderr(l,ori,sz) = std(frs(thisstillinds))./sqrt(length(thisstillinds));
                else
                    stillcondresp(l,ori,sz,:) = nan(1,size(resp,2));
                    stillcondfr(l,ori,sz) = NaN;
                    stillconderr(l,ori,sz) = NaN;
                end  
                
            end
        end
    end
    
    bincondresp = bincondresp.*(1000/binwidth);
    bta = bta-prestim;
%     psthcondplot(bincondresp,binconderr,bta)
    cllbinresp(cll,:,:,:,:) = bincondresp;
    cllbinerr(cll,:,:,:,:) = binconderr;
    
    contindsnl = find(gratingInfo.size == 0 & light == 0);
    controlresp(1,:) = mean(resp(contindsnl,:),1);
    controlfr(cll,1) = mean(frs(contindsnl));
    controlerr(1) = std(frs(contindsnl))./sqrt(length(contindsnl));
    controllfpspect(cll,1,:) = nanmean(lfpspect(contindsnl,:));
    
    contindsl = find(gratingInfo.size == 0 & light == 1);
    controlresp(2,:) = mean(resp(contindsl,:),1);
    controlfr(cll,2) = mean(frs(contindsl));
    controlerr(2) = std(frs(contindsl))./sqrt(length(contindsl));
    controllfpspect(cll,2,:) = nanmean(lfpspect(contindsl,:));
    
    [binnedctrlight,bta] = binit(controlresp(2,:),binwidth);
    [binnedctrnolight,bta] = binit(controlresp(1,:),binwidth);

    binnedcllrespl0(cll,:) = binnednolight;
    binnedcllrespl1(cll,:) = binnedlight;
    
    nolmaxfr(cll) = max(max(condfr(1,:,:)));
    if size(condfr,1) == 2
        lmaxfr(cll) = max(max(condfr(2,:,:)));
    end
    
    cllz(cll,:,:,:) = condz;
    cllsc(cll,:,:,:) = condsc;
    cllff(cll,:,:,:) = ff;
    cllfr(cll,:,:,:) = condfr;
    cllerr(cll,:,:,:) = conderr;
    clleckerrely(cll,:,:,:) = eckerreliability;
    cllmsrely(cll,:,:,:) = msreliab;
    cllbinrely(cll,:,:,:) = binreliab;
    clllfpspect(cll,:,:,:,:) = condlfpspect;
    cllruncondfr(cll,:,:,:) = runcondfr;
    cllstillcondfr(cll,:,:,:) = stillcondfr;
%     figure
%     errorbar(sizes,squeeze(nanmean(condfr(2,:,:),2)),squeeze(nanmean(conderr(2,:,:),2)),'o-','color',lcol,'markersize',8,'linewidth',2)
%     hold on
%     errorbar(sizes,squeeze(nanmean(condfr(1,:,:),2)),squeeze(nanmean(conderr(1,:,:),2)),'ko-','markersize',8,'linewidth',2)
%     xlabel('shown patch size [vd]')
%     ylabel('Firing rate [Hz]')
%     legend({'Light ON','Light OFF'})    
%     set(gca,'xtick',sizes)  
%     title(['cll ' int2str(cll) ' average all orientations'])
    
    prefsize = find(mean(condfr(1,:,:),2) == max(mean(condfr(1,:,:),2)),1);
    
    [nloriprefratio(cll), nldirprefratio(cll), nlprefori, nlmeanori, nlosi(cll), nlmeandir, nldsi(cll)] = getOSI(squeeze(condfr(1,:,prefsize)),oris);
%     [loriprefratio(cll), ldirprefratio(cll), lprefori, lmeanori, losi(cll), lmeandir, ldsi(cll)] = getOSI(squeeze(condfr(2,:,prefsize)),oris);
    
    figure
    subplot(2,2,1)
    ta = bta-prestim;
    plot(ta,binnednolight,'k','linewidth',2);
    hold on
    plot(ta,binnedlight,lcol,'linewidth',2);
    mx = max([max(binnednolight),max(binnedlight),.01]);
    axis([-prestim,trialdur+poststim,0,mx]);
    line([0,0],[0,mx],'color','k','linewidth',2);
    line([2000,2000],[0,mx],'color','k','linewidth',2);
    line([500,500],[0,mx],'color','b','linewidth',2)
    line([1500,1500],[0,mx],'color','b','linewidth',2);
    legend({'Light OFF','Light ON'})
    xlabel('time [ms]')
    ylabel('firing rate [Hz]')
    title(['cll ' int2str(cll) ' depth: ' int2str(result.depth), 'cll ' printname ])
    
    subplot(2,2,2)
    if size(condfr,1) == 2
        errorbar(oris,squeeze(condfr(2,:,prefsize)),squeeze(conderr(2,:,prefsize)),'o-','color',lcol,'markersize',8,'linewidth',2)
    end
    hold on
    errorbar(oris,squeeze(condfr(1,:,prefsize)),squeeze(conderr(1,:,prefsize)),'ko-','markersize',8,'linewidth',2)
    xlabel('shown orientation')
    ylabel('Firing rate [Hz]')
    set(gca,'xtick',oris)
    legend({'Light ON','Light OFF'})
%     title([' OSI: ' num2str(nlosi(cll)) ' OSI Light: ' num2str(losi(cll))])
    
    if length(oris)<8
        oneorifr = squeeze(condfr(:,:,prefsize));
    else
        oneorifr = mean(reshape(condfr(:,:,prefsize),length(unique(light)),size(condfr,2)/2,2),3);
    end
    prefori = find(oneorifr(1,:) == max(oneorifr(1,:)),1);
    ortho = mod(prefori+2,length(oris)/2); if ortho == 0, ortho = length(oris)/2; end
    
    cllcondresppreforil0(cll,:,:) = squeeze(condresp(1,prefori,:,:));
    if size(condresp,1) == 2
        cllcondresppreforil1(cll,:,:) = squeeze(condresp(2,prefori,:,:));
    end
    
    preffr(cll,:) = oneorifr(:,prefori);
    
    subplot(2,2,3)
    if size(conderr,1) == 2
        if length(oris)>4
            errorbar(sizes,squeeze(nanmean(condfr(2,[prefori,prefori+(length(oris)/2)],:))),...
            squeeze(nanmean(conderr(2,[prefori,prefori+(length(oris)/2)],:))),'o-','color',lcol,'markersize',8,'linewidth',2);    
        else
            errorbar(sizes,squeeze(condfr(2,prefori,:)),...
            squeeze(conderr(2,prefori,:)),'ko-','markersize',8,'linewidth',2);
        end
    end
    hold on
    if length(oris)>4
        errorbar(sizes,squeeze(nanmean(condfr(1,[prefori,prefori+(length(oris)/2)],:))),...
            squeeze(nanmean(conderr(1,[prefori,prefori+(length(oris)/2)],:))),'ko-','markersize',8,'linewidth',2);
    else
        errorbar(sizes,squeeze(condfr(1,prefori,:)),...
            squeeze(conderr(1,prefori,:)),'ko-','markersize',8,'linewidth',2);
    end
    xlabel('shown patch size [vd]')
    ylabel('Firing rate [Hz]')
    legend({'Light ON','Light OFF'})    
    set(gca,'xtick',sizes)  
    title(['cll ' int2str(cll) ' preferred orientations' ' depth: ' int2str(result.depth) '  ' printname])
    if size(condfr,1) == 2
        if length(oris)>4
            sizetunelas = squeeze(nanmean(condfr(2,[prefori,prefori+(length(oris)/2)],:)));
        else
            sizetunelas = squeeze(nanmean(condfr(2,prefori,:)));
        end
    else
        sizetunelas = NaN;
    end
    if length(oris)>4
        sizetunenolas = squeeze(nanmean(condfr(1,[prefori,prefori+(length(oris)/2)],:)));
    else
        sizetunenolas = squeeze(condfr(1,prefori,:));
    end
    
%     sil(cll) = (sizetunelas(find(sizetunelas == max(sizetunelas),1))-sizetunelas(end))/sizetunelas(find(sizetunelas == max(sizetunelas),1));
    sinl(cll) = (sizetunenolas(find(sizetunenolas == max(sizetunenolas),1))-sizetunenolas(end))/sizetunenolas(find(sizetunenolas == max(sizetunenolas),1));
    
    subplot(2,4,7)
    plot(spike)
    axis([0,40,-100,100])
    legend(['width: ' int2str(swidth(cll)) ' adiff: ' num2str(adiff(cll))])
    
    subplot(2,4,8)
    plot(ta,binnedctrnolight)
    hold on
    if size(condfr,1) == 2
        plot(ta,binnedctrlight,lcol)
    end
    axis([0,2500,0,1])
    
    o = 2; s = 3;
%     condl0 = (gratingInfo.Orientation == oris(o) | gratingInfo.Orientation == oris(o+4)) & (gratingInfo.size == sizes(s) | gratingInfo.size == sizes(s+1)) & light == 0;
%     condl1 = (gratingInfo.Orientation == oris(o) | gratingInfo.Orientation == oris(o+4)) & (gratingInfo.size == sizes(s) | gratingInfo.size == sizes(s+1)) & light == 1;
%     condl0 = gratingInfo.Orientation == oris(o) & gratingInfo.size == sizes(s) & light == 0;
%     condl1 = gratingInfo.Orientation == oris(o) & gratingInfo.size == sizes(s) & light == 1;
    condl0 = light == 0; condl1 = light == 1;
%     figure, rasterplot(resp,condl0,condl1,msta);
%     line([0,2000],[37,37],'color','k','linewidth',2)
%     line([500,1500],[33,33],'color','r','linewidth',2)
    
    % running figure
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
    lp(cll) = p(1); rp(cll) = p(2); rlip(cll) = p(3);
    
    r0omi(cll) = (norunlfr-norunnlfr)/(norunlfr+norunnlfr);
    r1omi(cll) = (runlfr-runnlfr)/(runlfr+runnlfr);
    l0rmi(cll) = (runnlfr-norunnlfr)/(runnlfr+norunnlfr);
    l1rmi(cll) = (runlfr-norunlfr)/(runlfr+norunlfr);
    
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
    
    subplot(2,2,4)
    barweb([nlfr(cll),lfr(cll);runnlfr,runlfr;norunnlfr,norunlfr],...
        [nlfrerr,lfrerr;runnlfrerr,runlfrerr;norunnlfrerr,norunlfrerr],...
        [],[{'all'};{'running only'};{'immobile only'}],['ANOVA factor running p: ' num2str(p(2))],...
        [],'firing rate [Hz]',[],[]);
    
%     %running with LFP
%     for i = 1:20
%         a = squeeze(mean(allpowresp(i,:,respwin),3));
%         b = mean(speed(:,respwin),2);
%         [r,p] = corrcoef(a,b);
%         rr(i) = r(1,2);
%         pp(i) = p(1,2);
%     end
%     fax = 3:5:100
%     
%     figure
%     plot(fax,rr,'o','markerfacecolor','b')
%     hold on
%     ss = find(pp<.05)
%     plot(fax(ss),rr(ss),'ro','markerfacecolor','r')
    
    cll = cll + 1;
    disp('');
    
end


figure
plot(swidth,adiff,'k.')
xlabel('spike width')
ylabel('amplitude diff')

secpersamp = 1/30000;
interpf = secpersamp/10;
swidthms = swidth*interpf*1000;
pfs = find(swidth<125);

fsi = kmeans([swidth',adiff'],2); %kmeans([eslope',ptr',swidth',adiff'],2);
if mean(swidth(find(fsi==1)))<mean(swidth(find(fsi==2)))  %1 is FS
    pfs = find(fsi==1); prs = find(fsi==2);
else
    pfs = find(fsi==2); prs = find(fsi==1);
end

params.Fs = 1000; params.trialave = 1; params.err = [2 .05]; params.tapers = [5,9];

%coherence
cll1 = 30; cll2 = 9;
for i = 1:size(clllfpresp,2)
    [coh(i,:),cfx] = mscohere(squeeze(clllfpresp(cll1,i,respwin)),squeeze(clllfpresp(cll2,i,respwin)),[],[],512,1000);
end
for l = 1:2
    for ori = 1:length(oris)
        for sz = 1:length(sizes)
            thisinds = find(gratingInfo.Orientation == oris(ori) &...0
                gratingInfo.size == sizes(sz) & ...
                light == l-1);
            condcoher(l,ori,sz,:) = nanmean(coh(thisinds,:),1);
        end
    end
end
contcoher(1,:) = nanmean(coh(contindsnl,:),1);
contcoher(2,:) = nanmean(coh(contindsl,:),1);

for l = 1:2
    for sz = 1:length(sizes)        
        [a,b,c,d,e,ff,g,h,i] = coherencyc(squeeze(sizelfpresp(cll1,l,sz,:,respwin))',squeeze(sizelfpresp(cll2,l,sz,:,respwin))',params);
        Ccoh(l,sz,:) = a;
        phicc(l,sz,:) = b;
        S12cc(l,sz,:) = c;
        S1cc(l,sz,:) = d;
        S2cc(l,sz,:) = e;
        chfx = ff;
        confCcc(l,sz) = g;
        phistdcc(l,sz,:) = h;
        Ccoherr(l,sz,:,:) = i;
    end    
end

fillx = [sfx(1:124),fliplr(sfx(1:124))];
fillxcoh = [chfx(1:124),fliplr(chfx(1:124))];
for l = 1:length(unique(light))
    for i = 1:5
        filly(l,i,:) = [squeeze(condSerr(cll1,l,i,1,1:124))',fliplr(squeeze(condSerr(cll1,l,i,2,1:124))')];
%         fillcoh(l,i,:) = [squeeze(Ccoherr(l,i,1,1:124))',fliplr(squeeze(Ccoherr(l,i,2,1:124))')];
%         fillsfc(l,i,:) = [squeeze(Cerr(cll1,l,i,1,1:124))',fliplr(squeeze(Cerr(cll1,l,i,2,1:124))')];
    end
end

figure
fill(fillx,squeeze(filly(1,1,:)),[.8,.8,1]); %,'FaceAlpha',.5)
hold on
fill(fillx,squeeze(filly(1,3,:)),[.5,.5,1]); %,'FaceAlpha',.5)
fill(fillx,squeeze(filly(1,5,:)),[.3,.3,1]); %'FaceAlpha',.5)
% plot(sfx,squeeze(condS(cll1,1,1,:)),'color',[.7,.7,1])
% plot(sfx,squeeze(condS(cll1,1,3,:)),'color',[.4,.4,1])
% plot(sfx,squeeze(condS(cll1,1,5,:)),'color',[.2,.2,1])
set(gca,'yscale','log')
set(gca,'yscale','lin')

sizefillx = [sizes,fliplr(sizes)];
peakf = 27;
sizefilly = [squeeze(condSerr(cll1,1,:,1,27))',fliplr(squeeze(condSerr(cll1,1,:,2,27))')];
figure
fill(sizefillx,sizefilly,[.5,.5,1])
hold on
plot(sizes,squeeze(condS(cll1,1,:,27)),'linewidth',2)

for i = 1:2
    for j = 1:5
        blerrs(i,j,:,2) = squeeze(Ccoherr(i,j,2,:))-squeeze(Ccoh(i,j,:));
        blerrs(i,j,:,1) = squeeze(Ccoh(i,j,:))-squeeze(Ccoherr(i,j,1,:));
        blserrs(i,j,:,2) = squeeze(condSerr(cll1,i,j,2,:))-squeeze(condS(cll1,i,j,:));
        blserrs(i,j,:,1) = squeeze(condS(cll1,i,j,:))-squeeze(condSerr(cll1,i,j,1,:));
    end
end
% figure
% boundedline(chfx,squeeze(Ccoh(1,1,:)),squeeze(blerrs(1,1,:,:)),'k:','alpha')
% hold on
% boundedline(chfx,squeeze(Ccoh(1,5,:)),squeeze(blerrs(1,5,:,:)),'k','alpha')
% boundedline(chfx,squeeze(Ccoh(2,5,:)),squeeze(blerrs(2,5,:,:)),'r','alpha')
% boundedline(chfx,squeeze(Ccoh(2,1,:)),squeeze(blerrs(2,1,:,:)),'r:','alpha')
% 
% figure
% boundedline(sfx,squeeze(condS(cll1,1,1,:)),squeeze(blserrs(1,1,:,:)),'k:','alpha')
% hold on
% boundedline(sfx,squeeze(condS(cll1,1,5,:)),squeeze(blserrs(1,5,:,:)),'k','alpha')
% boundedline(sfx,squeeze(condS(cll1,2,5,:)),squeeze(blserrs(2,5,:,:)),'r','alpha')
% boundedline(sfx,squeeze(condS(cll1,2,1,:)),squeeze(blserrs(2,1,:,:)),'r:','alpha')

% cross correlations
cll1 = 9; cll2 = 30;
a = find(gratingInfo.size == 4 & light == 0);
b = find(gratingInfo.size == max(unique(gratingInfo.size)) & light == 0);

% spike shuffle bootstrap
for bs = 1:10000
    for i = 1:40
        s1 = find(squeeze(cllresp(cll1,a(i),801:1800)));
        s2 = find(squeeze(cllresp(cll2,a(i),801:1800)));
        l1 = find(squeeze(cllresp(cll1,b(i),801:1800)));
        l2 = find(squeeze(cllresp(cll2,b(i),801:1800)));
        alls = [s1',s2']; alll = [l1',l2'];
        
        rns = randperm(length(alls));
        new1 = zeros(1,1000); new2 = zeros(1,1000);
        new1(alls(rns(1:length(s1)))) = 1;
        new2(alls(rns(length(s1)+1:end))) = 1;
        [bssssxc(bs,i,:),lg] = xcorr(new1,new2);        
        
        rnl = randperm(length(alll));
        new1 = zeros(1,1000); new2 = zeros(1,1000);
        new1(alll(rnl(1:length(l1)))) = 1;
        new2(alll(rnl(length(l1)+1:end))) = 1;
        [bssslxc(bs,i,:),lg] = xcorr(new1,new2);
    end
    sssmxc(bs,:) = mean(bssssxc(bs,:,:),2); % spike shuffle small mean xcorr
    sslmxc(bs,:) = mean(bssslxc(bs,:,:),2); % spike shuffle large mean xcorr
end
sortedsss = sort(sssmxc,1);
sortedssl = sort(sslmxc,1);
       
% trial shuffle 
for bs = 1:10000
    c = randperm(40);
    d = randperm(40);
    for i = 1:40
        [bssxc(bs,i,:),lg] = xcorr(squeeze(cllresp(2,a(i),801:1800)),squeeze(cllresp(17,a(c(i)),801:1800)));
        [bslxc(bs,i,:),lg] = xcorr(squeeze(cllresp(2,b(i),801:1800)),squeeze(cllresp(17,b(d(i)),801:1800)));
    end
    tslmxc(bs,:) = mean(bslxc(bs,:,:),2); % trial shuffled large mean xcorr
    tssmxc(bs,:) = mean(bssxc(bs,:,:),2); % trial shuffled small mean xcorr
end
                
% correlations
infrs = intersect(find(depth>500),prs);
irsresp = cllresp(infrs,:,:);
nclls = size(irsresp,1);

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
    binsizes = [1,5,10,15,20,25,50,100];;
    for ccbin = 1:length(binsizes)
        pair = 1;
        for i = 1:nclls-1
            for j = i+1:nclls

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
%                 ccbin = 50;
                for ii = 1:size(resp1,2)/binsizes(ccbin);
                    binresp1(:,ii) = sum(resp1(:,(ii-1)*binsizes(ccbin)+1:ii*binsizes(ccbin)),2);
                    binresp2(:,ii) = sum(resp2(:,(ii-1)*binsizes(ccbin)+1:ii*binsizes(ccbin)),2);
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
        rsl0(ccbin,:) = nanmean(rl0,2);
        rsl1(ccbin,:) = nanmean(rl1,2);
    end

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
    title(['all pairwise correlations between lower layer RS clls p: ' num2str(p)]);
%     if printyn
%         figSize = [30 21];
%         set(gcf,'paperunits','centimeters','papersize',figSize,'paperposition',[0 0 figSize])
%         print(gcf,[popprintpath ,  '18_LLRScorrelations.pdf'], '-dpdf' );
%     end
    
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



%lifetime sparseness
N = 34; %1000ms binned 28:61 on bta (500:1500)
for cl = 1:size(cllbinresp,1)
    for l = 1:2
        for o = 1:8
            for s = 1:5
                slt(cl,l,o,s) = (1 - ((1/N)* ((sum(cllbinresp(cl,l,o,s,28:61),5).^2)/(sum(cllbinresp(cl,l,o,s,28:61).^2,5))) ))...
                 / (1-(1/N));
            end
        end
    end
end

%population sparseness
N = size(cllbinresp,1);
for bn = 1:34
    for l = 1:2
        for o = 1:8
            for s = 1:5
                spop(bn,l,o,s) = (1 - ((1/N)* ((sum(cllbinresp(:,l,o,s,bn+27),1).^2)/(sum(cllbinresp(:,l,o,s,bn+27).^2,1))) ))...
                 / (1-(1/N));
            end
        end
    end
end

ii = 1;
for i = 1:size(cllz,1)-1
    for j = i+1:size(cllz,1)
        for l = 1:2
            for o = 1:8
                for s = 1:5
                    cv = cov(cllz{i,l,o,s},cllz{j,l,o,s});
                    cvs(ii,l,o,s) = cv(1,2);
                    cvcount = cov(cllsc{i,l,o,s},cllsc{j,l,o,s});
                    if size(cvcount) == 1
                        cvcounts(ii,l,o,s) = NaN;
                    else
                        cvcounts(ii,l,o,s) = cvcount(1,2);
                    end
                    cc = corrcoef(cllz{i,l,o,s},cllz{j,l,o,s});
                    if isnan(cc)
                        cc(ii,l,o,s) = NaN;
                    else
                        ccs(ii,l,o,s) = cc(1,2);
                    end
                    pairinds(ii,:) = [i,j];
                end
            end
        end
        ii = ii+1;
    end
end

mscrl0 = nanmean(nanmean(cvs(:,1,:,:),3),4);
mscrl1 = nanmean(nanmean(cvs(:,2,:,:),3),4);
mscountrl0 = nanmean(nanmean(cvcounts(:,1,:,:),3),4);
mscountrl1 = nanmean(nanmean(cvcounts(:,2,:,:),3),4);
mccl0 = nanmean(nanmean(ccs(:,1,:,:),3),4);
mccl1 = nanmean(nanmean(ccs(:,2,:,:),3),4);
mffl0 = nanmean(nanmean(cllff(:,1,:,:),3),4);
mffl1 = nanmean(nanmean(cllff(:,2,:,:),3),4);



figure
[s,p] = ttest(nlosi,losi);
plot(nlosi,losi,'k.')
line([0,1],[0,1],'color','k')
axis square
xlabel('OSI control');
ylabel('OSI light');
title([' p: ' num2str(p)])

figure
[s,p] = ttest(nldsi,ldsi);
plot(nldsi,ldsi,'k.')
line([0,1],[0,1],'color','k')
axis square
xlabel('DSI control');
ylabel('DSI light');
title([' p: ' num2str(p)])

figure
plot(nolmaxfr,lmaxfr,'ko')
% axis([-5,25,-5,25])
axis square
line([-5,25],[-5,25],'color','k')
xlabel('max firing rate no light [Hz]')
ylabel('max firing rate light [Hz]')

figure
plot(nlfr,lfr,'k.')
hold on
plot(nlfr(pfs),lfr(pfs),'ko')
axis square
line([0,30],[0,30],'color','k')
xlabel('mean firing rate no light [Hz]')
ylabel('mean firing rate light [Hz]')

figure
plot((lfr-nlfr)./(lfr+nlfr),depth,'k.')
hold on
plot((lfr(pfs)-nlfr(pfs))./(lfr(pfs)+nlfr(pfs)),depth(pfs),'ko')
line([0,0],[0,1000],'color','k')
axis ij
axis([-1.1,1.1,0,1000])
ylabel('depth[mum]')
xlabel('(light-nolight)/(light+nolight)')
title('mean firing rate changes by layer 4 suppression')

figure
plot((preffr(:,2)-preffr(:,1))./(preffr(:,2)+preffr(:,1)),depth,'k.')
hold on
plot((preffr(pfs,2)-preffr(pfs,1))./(preffr(pfs,2)+preffr(pfs,1)),depth(pfs),'ko')
line([0,0],[0,1000],'color','k')
axis ij
axis([-1.1,1.1,0,1000])
ylabel('depth[mum]')
xlabel('(light-nolight)/(light+nolight)')
title('preferred firing rate changes by layer 4 suppression')

figure
plot((controlfr(:,2)-controlfr(:,1))./(controlfr(:,2)+controlfr(:,1)),depth,'k.')
hold on
plot((controlfr(pfs,2)-controlfr(pfs,1))./(controlfr(pfs,2)+controlfr(pfs,1)),depth(pfs),'ko')
line([0,0],[0,1000],'color','k')
axis ij
axis([-1.1,1.1,0,1000])
ylabel('depth[mum]')
xlabel('(light-nolight)/(light+nolight)')
title('control firing rate changes by layer 4 suppression')

figure
plot(sinl,sil,'k.')
hold on
plot(sinl(pfs),sil(pfs),'ko')
axis([0,1,0,1])
axis square
line([0,1],[0,1],'color','k')
xlabel('SI ((max-largest)/max) no light')
ylabel('SI ((max-largest)/max) light')

figure
plot(sil-sinl,depth,'k.')
hold on
plot(sil(pfs)-sinl(pfs),depth(pfs),'ko')
axis ij
line([0,0],[250,750],'color','k')
ylabel('depth[mum]')
xlabel('SI light - SI no light')

disp('');