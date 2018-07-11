function figuregroundanalysis_SU_allcells

animalid = '180605';
block = 8;
rfblock = 9;
lcol = 'r'; %lasercolor

onlymod = 0;
printyn = 0;
sfc = 0;

addpath C:\Users\Julia\work\Matlab\others\fieldtrip-20160329
ft_defaults

supath = ['C:\Users\Julia\work\data\' animalid '\singleunits\'];
% supath = ['C:\Users\Julia\work\data\' animalid '\multiunits\'];
basename = [animalid '_block' int2str(block) '_tet'];
snbasename = [animalid '_block' int2str(rfblock) '_tet'];

files = dir([supath, basename, '*.mat']);
rffiles = dir([supath, snbasename, '*.mat']);

freqbinwidth = 5;

% chronux parameters
params.tapers = [5,9]; params.Fs = 1000; params.err = [2, 0.05]; params.trialave = 1;

%Scott Kernel
%%%%%%%%%%%%kernal propeties%%%%%%%%%%%%%%%
kernel_width_eval_s = 0.15;
sdf_freq_hz = 3700;
exp_growth_ms = 2;
exp_decay_ms = 12;
exp_growth_s = exp_growth_ms/1000;
exp_decay_s = exp_decay_ms/1000;
eval_kernel_x_s = 0:1/sdf_freq_hz:kernel_width_eval_s;
exp_kernel = eval_kernel_x_s;
exp_kernel = (1-(exp(-(exp_kernel./exp_growth_s)))).*(exp(-(exp_kernel./exp_decay_s)));
excit_kernel = exp_kernel/sum(exp_kernel); 
%%%%%%%%%%kernal properties%%%%%%%%%%%%%%%

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
    if strcmp(files(fi).name(i+4),'_')
        tetno = strread(files(fi).name(i+3)); % single character number
    else
        tetno = strread(files(fi).name(i+3:i+4)); % number >10
    end
    
    wvchan = find(var(result.waveforms) == max(var(result.waveforms)));
    
    sr = 1000;
    lfp = result.lfp(:,wvchan)';
    nfft = 2^nextpow2(length(lfp));
    fax = sr/2*linspace(0,1,nfft/2+1);
    y = fft(lfp,nfft);
    lfpspectrum = abs(y(1:nfft/2+1));
    
    if strcmp (animalid, '180605')
        bipderlfp = result.lfp(:,1)-result.lfp(:,3);
    else
        if ~(wvchan == 4)
            bipderlfp = result.lfp(:,wvchan)-result.lfp(:,wvchan+1);
        else
            bipderlfp = result.lfp(:,3)-result.lfp(:,4);
        end
    end

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
%          msstamps([62,108,147]) = []; % for 140703 block 8
%          newmsstamps = msstamps(1:319); % for 170712 block 3
%          newmsstamps(320) = newmsstamps(319)+3683;
%          newmsstamps(321) = newmsstamps(320)+3683;
%          newmsstamps(322) = newmsstamps(321)+3683;
%          newmsstamps = [newmsstamps; msstamps(320:end)];
%             %terrible hack 20171201
%             mscopy(1:384) = msstamps(1:384);
%             mscopy(385) = msstamps(384)+2355;
%             mscopy(386:418) = msstamps(385:end);
%             mscopy2(1:396) = mscopy(1:396);
%             mscopy2(397) = mscopy(396)+2355;
%             mscopy2(398:419) = mscopy(397:418);
%             mscopy3(1:416) = mscopy2(1:416);
%             mscopy3(417) = mscopy2(416)+2355;
%             mscopy3(418:420) = mscopy2(417:419);
%             newmsstamps = mscopy3;
    
%          msstamps = newmsstamps;
%          result.msstamps = msstamps;
%          save([supath, files(fi).name],'result');
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
    
    for i = 1:length(msstamps)
        resp(i,:) = chan(msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim); 
        hh = find(resp(i,1001:1800))'; ptresp(i).times = hh./1000;
        lfpresp(i,:) = result.lfp(msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim, wvchan);
        bdlfpresp(i,:) = bipderlfp(msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim);
        filtlfpresp(i,:) = filtlfp(msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim);
        [lfpspect(i,:),trialfax] = pmtm(lfpresp(i,respwin(1)+200:respwin(end)),3,[],sr);
        [lfpsectrogram(i,:,:),ct,cf] = mtspecgramc(lfpresp(i,:)',[.25,.1],params);
        gammaresp(i,:) = gamma(msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim);
        gammapowresp(i,:) = gpow(msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim);
        gphaseresp(i,:) = gphas(msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim);
        bpresp(i,:) = gamma(msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim);        
        
        speed(i,:) = result.runspeed(msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim);  
        
        clear lfps; clear speeds;
        s = find(resp(i,respwin));
        if ~isempty(s)
            for si = 1:length(s)
                lfps(si,:) = lfpresp(i,s(si)+respwin(1)-1-200:s(si)+respwin(1)-1+200);
                speeds(si,:) = speed(i,s(si)+respwin(1)-1-200:s(si)+respwin(1)-1+200);
            end
        else
            lfps = nan(1,401);
            speeds = nan(1,401);
        end
        trialstalfp(i,:) = mean(lfps,1);
        trialstaspeed(i,:) = mean(speeds,1);
        nspkstrialsta(i) = size(lfps,1);
        clear lfps;
        [spkxcorr(i,:),lags] = xcorr(resp(i,respwin),resp(i,respwin));
        
        if sfc
            for j = 1:size(phasmat,1)
                allphaseresp(j,i,:) = phasmat(j, msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim);
                allpowresp(j,i,:) = powmat(j, msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim);
            end
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
       

    depth(cll) = result.depth;
    
    spike = result.waveforms(:,wvchan);
    interpspike = spline(1:32,spike,1:.1:32);
    [adiff(cll),swidth(cll)] = spikequant(interpspike);
    
    cllresp(cll,:,:) = resp;
    clllfpresp(cll,:,:) = lfpresp;
    cllbdlfpresp(cll,:,:) = bdlfpresp;
    
    msta = linspace(-prestim,trialdur+poststim,size(resp,2));
    
    frs = sum(resp(:,respwin),2)./(length(respwin)/1000);
    bl = sum(resp(:,1:prestim),2)./(prestim/1000);
    sc = sum(resp(:,respwin),2);
        
     %determine if cll is visually modulated
    blfr = sum(resp(:,1:prestim),2);
    vrfr = sum(resp(:,prestim+40:2*prestim+40),2);
    vismod(cll) = ttest2(blfr,vrfr);
    
    %determine if cll is modulated by light
%     lightmod(cll) = ttest2(frs(find(light)),frs(find(light == 0)));

    %get RF of cell, too
    [on,off,gaussfiton,gaussfitoff,rson(cll),rsoff(cll),xax,yax] = get_rf([supath, rffiles(fi).name]);
    

    binwidth = 40;
    clear stalfp; clear nspkssta; clear shufstalfp; clear nspksshufsta;
    sizes = unique(gratingInfo.sizecenter);
    orisc = unique(gratingInfo.Orientation_center);
    oriss = unique(gratingInfo.Orientation_surround);
    for lg = 1:length(unique(light))
        for sz = 1:length(sizes)
            for oc = 1:length(orisc)
                for os = 1:length(oriss)
                    trials = find(gratingInfo.sizecenter == sizes(sz) &...
                        gratingInfo.Orientation_center == orisc(oc) & ...
                        gratingInfo.Orientation_surround == oriss(os) & light == lg-1);

                    
                    condresp(cll,lg,sz,oc,os,:) = mean(resp(trials,:),1);
                    condlfpspect(cll,lg,sz,oc,os,:) = nanmean(lfpspect(trials,:));
                    condstspect(cll,lg,sz,oc,os,:) = pmtm(mean(resp(trials,1001:1800),1),3,[],sr);
                    condfr(cll,lg,sz,oc,os) = mean(frs(trials));%-mean(bl);
                    conderr(cll,lg,sz,oc,os) =std(frs(trials))./sqrt(length(trials));
                    
                    condfiltresp(cll,lg,sz,oc,os,:) = filter(excit_kernel,1,mean(resp(trials,:),1));
                    
                    condz(cll,lg,sz,oc,os) = {(sc(trials)-mean(sc(trials)))/std(sc(trials))}; %ecker 2010
                    condsc(cll,lg,sz,oc,os) = {sc(trials)};
                    ff(cll,lg,sz,oc,os) = var(sc(trials))/mean(sc(trials));
                    
                    condlfpresp(cll,lg,sz,oc,os,:) = mean(lfpresp(trials,:),1);
                    
                    condresperr(cll,lg,sz,oc,os,:) = nanstd(resp(trials,:),1,1)./sqrt(length(trials));
                    if ~isnan(condresp(cll,lg,sz,oc,os,:))
                        [bincondresp(cll,lg,sz,oc,os,:),bta] = binit(condresp(cll,lg,sz,oc,os,:),binwidth);
                    else
                        bincondresp(cll,lg,sz,oc,os,:) = binit(condresp(cll,lg,sz,oc,os,:),binwidth);
                    end
                    binconderr(cll,lg,sz,oc,os,:) = binit(condresperr(cll,lg,sz,oc,os,:),binwidth);
                    
                    
                    if ~isempty(trials)
                        [S,chf,Serr]=mtspectrumc(squeeze(lfpresp(trials,1001:1800))',params);
                        condS(cll,lg,sz,oc,os,:) = S(1:150); condSerr(cll,lg,sz,oc,os,:,:) = Serr(:,1:150);
                    else
                        condS(cll,lg,sz,oc,os,:) = nan(1,150);
                    end

                    thisruninds = intersect(trials,oktrials); thisstillinds = intersect(trials,stilltrials);
                    if ~isempty(thisruninds)
                        runcondresp(cll,lg,sz,oc,os,:) = mean(resp(thisruninds,:),1);
                        runcondfr(cll,lg,sz,oc,os) = mean(frs(thisruninds));
                        runconderr(cll,lg,sz,oc,os) = std(frs(thisruninds))./sqrt(length(thisruninds));
                        
                        thisrunn(cll,lg,sz,oc,os) = length(thisruninds);
                        [runS,chf,runSerr] = mtspectrumc(squeeze(lfpresp(thisruninds,1001:1800))',params);
                        condrunS(cll,lg,sz,oc,os,:) = runS(1:150); condrunSerr(cll,lg,sz,oc,os,:,:) = runSerr(:,1:150);                
                    else
                        runcondresp(cll,lg,sz,oc,os,:) = nan(1,size(resp,2));
                        runcondfr(cll,lg,sz,oc,os) = NaN;
                        runconderr(cll,lg,sz,oc,os) = NaN;
                        
                        thisrunn(cll,lg,sz,oc,os) = 0;
                        condrunS(cll,lg,sz,oc,os,:) = nan(1,150); condrunSerr(cll,lg,sz,oc,os,:,:) = nan(2,150);
                    end
                    if ~isempty(thisstillinds)
                        stillcondresp(cll,lg,sz,oc,os,:) = mean(resp(thisstillinds,:),1);
                        stillcondfr(cll,lg,sz,oc,os) = mean(frs(thisstillinds));
                        stillconderr(cll,lg,sz,oc,os) = std(frs(thisstillinds))./sqrt(length(thisstillinds));
                        
                        thisstilln(cll,lg,oc,os,sz) = length(thisstillinds);
                        [stillS,chf,stillSerr] = mtspectrumc(squeeze(lfpresp(thisstillinds,1001:1800))',params);
                        condstillS(cll,lg,oc,os,sz,:) = stillS(1:150); condstillSerr(cll,lg,oc,os,sz,:,:) = stillSerr(:,1:150);
                    else
                        stillcondresp(cll,lg,sz,oc,os,:) = nan(1,size(resp,2));
                        stillcondfr(cll,lg,sz,oc,os) = NaN;
                        stillconderr(cll,lg,sz,oc,os) = NaN;
                        
                        thisstilln(cll,lg,sz,oc,os) = 0;
                        condstillS(cll,lg,sz,oc,os,:) = nan(1,150); condstillSerr(cll,lg,sz,oc,os,:,:) = nan(2,150);
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

                        stalfp(cll,lg,sz,oc,os,i,:) = mean(lfps,1);
                        nspkssta(cll,lg,sz,oc,os,i) = size(lfps,1);
                        shufstalfp(cll,lg,sz,oc,os,i,:) = mean(shuflfps,1);
                        nspksshufsta(cll,lg,sz,oc,os,i) = size(shuflfps,1);
                        clear lfps; clear shuflfps;
                    end


                    % field trip messing around
                    if ~isempty(trials)
                        for i = 1:length(trials)
                            data.trial{i} = [lfpresp(trials(i),:);resp(trials(i),:)];
                            time = -prestim+1:result.stimduration*1000+poststim;
                            data.time{i} = time(1)./1000:.001:time(end)./1000;
            %                 data.time{i} = -.299:.001:2.700;
                        end
                        data.fsample = 1000;
                        data.trialinfo(:,1) = light(trials);
                        data.trialinfo(:,2) = gratingInfo.sizecenter(trials);
                        data.trialinfo(:,3) = gratingInfo.Orientation_center(trials);
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
                        ftspect(cll,lg,sz,:) = squeeze(nanmean(mag(:,1,:)));
                        ftphases{cll,lg,sz} = squeeze(ang);
                        ftfax = stsFFT.freq;

                        cfg               = [];
                        cfg.method        = 'ral'; % compute the rayleigh test
                        cfg.spikechannel  = stsFFT.label{1};
                        cfg.channel       = stsFFT.lfplabel; % selected LFP channels
                        cfg.avgoverchan   = 'unweighted'; % weight spike-LFP phases irrespective of LFP power
                        cfg.timwin        = 'all'; % compute over all available spikes in the window
                        cfg.latency       = [0.5 1.5]; % sustained visual stimulation period
                        statSts           = ft_spiketriggeredspectrum_stat(cfg,stsFFT);
                        ftralp(cll,lg,sz,oc,os,:) = statSts.ral;
                        cfg.method = 'ppc0';
                        statSts           = ft_spiketriggeredspectrum_stat(cfg,stsFFT);
                        ftppc(cll,lg,sz,oc,os,:) = statSts.ppc0;
                        cfg.method = 'plv';
                        statSts           = ft_spiketriggeredspectrum_stat(cfg,stsFFT);
                        ftplv(cll,lg,sz,oc,os,:) = statSts.plv;
                    else
                        ftralp(cll,lg,sz,oc,os,:) = nan(1,30);
                        ftppc(cll,lg,sz,oc,os,:) = nan(1,30);
                        ftplv(cll,lg,sz,oc,os,:) = nan(1,30);
                    end
                end
            end
        end
    end
    sfx = chf(1:150);
    
    bta = bta-prestim;
    
    printname = files(fi).name;
    printname(find(printname=='_')) = ' ';
    
    % size tuning (cll,l,sz,orisc,oriss)
    figure    
    errorbar(sizes,squeeze(nanmean(condfr(cll,1,:,2:3,1),4)),squeeze(nanmean(conderr(cll,1,:,2:3,1),4)),'ko-','markersize',8,'linewidth',2)
    hold on
    if size(condfr,2) == 2
        errorbar(sizes,squeeze(nanmean(condfr(cll,2,:,2:3,1),4)),squeeze(nanmean(conderr(cll,2,:,2:3,1),4)),'ko-','color',lcol,'markersize',8,'linewidth',2)
    end
    xlabel('shown patch size [vd]')
    ylabel('Firing rate [Hz]')
    legend({'Light OFF','Light ON'})    
    set(gca,'xtick',sizes)  
    title(['cll ' int2str(cll) ' average all orientations'])
    
   
    % BL, sm md lg bgonly smbg mdbg lgbg smap mdap lgap
    for lg = 1:size(condfr,2)
        vec(cll,lg,1) = condfr(cll,lg,1,1,1); vecerr(cll,lg,1) = conderr(cll,lg,1,1,1);
        vec(cll,lg,2) = nanmean(condfr(cll,lg,1,2:3,1),4); vecerr(cll,lg,2) = nanmean(conderr(cll,lg,1,2:3,1),4);
        vec(cll,lg,3) = nanmean(condfr(cll,lg,2,2:3,1),4); vecerr(cll,lg,3) = nanmean(conderr(cll,lg,2,2:3,1),4);
        vec(cll,lg,4) = nanmean(condfr(cll,lg,3,2:3,1),4); vecerr(cll,lg,4) = nanmean(conderr(cll,lg,3,2:3,1),4);
        vec(cll,lg,5) = nanmean([condfr(cll,lg,1,2,2),condfr(cll,lg,1,3,3)]); vecerr(cll,lg,5) = nanmean([conderr(cll,lg,1,2,2),conderr(cll,lg,1,3,3)]);
        vec(cll,lg,6) = nanmean([condfr(cll,lg,1,2,3),condfr(cll,lg,1,3,2)]); vecerr(cll,lg,6) = nanmean([conderr(cll,lg,1,2,3),conderr(cll,lg,1,3,2)]);
        vec(cll,lg,7) = nanmean([condfr(cll,lg,2,2,3),condfr(cll,lg,2,3,2)]); vecerr(cll,lg,7) = nanmean([conderr(cll,lg,2,2,3),conderr(cll,lg,2,3,2)]);
        vec(cll,lg,8) = nanmean([condfr(cll,lg,3,2,3),condfr(cll,lg,3,3,2)]); vecerr(cll,lg,8) = nanmean([conderr(cll,lg,3,2,3),conderr(cll,lg,3,3,2)]);
        vec(cll,lg,9) = nanmean(condfr(cll,lg,1,1,2:3),5); vecerr(cll,lg,9) = nanmean(conderr(cll,lg,1,1,2:3),5);
        vec(cll,lg,10) = nanmean(condfr(cll,lg,2,1,2:3),5); vecerr(cll,lg,10) = nanmean(conderr(cll,lg,2,1,2:3),5);
        vec(cll,lg,11) = nanmean(condfr(cll,lg,3,1,2:3),5); vecerr(cll,lg,11) = nanmean(conderr(cll,lg,3,1,2:3),5);
        
        rvec(cll,lg,1) = runcondfr(cll,lg,1,1,1); rvecerr(cll,lg,1) = runconderr(cll,lg,1,1,1);
        rvec(cll,lg,2) = nanmean(runcondfr(cll,lg,1,2:3,1),4); rvecerr(cll,lg,2) = nanmean(runconderr(cll,lg,1,2:3,1),4);
        rvec(cll,lg,3) = nanmean(runcondfr(cll,lg,2,2:3,1),4); rvecerr(cll,lg,3) = nanmean(runconderr(cll,lg,2,2:3,1),4);
        rvec(cll,lg,4) = nanmean(runcondfr(cll,lg,3,2:3,1),4); rvecerr(cll,lg,4) = nanmean(runconderr(cll,lg,3,2:3,1),4);
        rvec(cll,lg,5) = nanmean([runcondfr(cll,lg,1,2,2),runcondfr(cll,lg,1,3,3)]); rvecerr(cll,lg,5) = nanmean([runconderr(cll,lg,1,2,2),runconderr(cll,lg,1,3,3)]);
        rvec(cll,lg,6) = nanmean([runcondfr(cll,lg,1,2,3),runcondfr(cll,lg,1,3,2)]); rvecerr(cll,lg,6) = nanmean([runconderr(cll,lg,1,2,3),runconderr(cll,lg,1,3,2)]);
        rvec(cll,lg,7) = nanmean([runcondfr(cll,lg,2,2,3),runcondfr(cll,lg,2,3,2)]); rvecerr(cll,lg,7) = nanmean([runconderr(cll,lg,2,2,3),runconderr(cll,lg,2,3,2)]);
        rvec(cll,lg,8) = nanmean([runcondfr(cll,lg,3,2,3),runcondfr(cll,lg,3,3,2)]); rvecerr(cll,lg,8) = nanmean([runconderr(cll,lg,3,2,3),runconderr(cll,lg,3,3,2)]);
        rvec(cll,lg,9) = nanmean(runcondfr(cll,lg,1,1,2:3),5); rvecerr(cll,lg,9) = nanmean(runconderr(cll,lg,1,1,2:3),5);
        rvec(cll,lg,10) = nanmean(runcondfr(cll,lg,2,1,2:3),5); rvecerr(cll,lg,10) = nanmean(runconderr(cll,lg,2,1,2:3),5);
        rvec(cll,lg,11) = nanmean(runcondfr(cll,lg,3,1,2:3),5); rvecerr(cll,lg,11) = nanmean(runconderr(cll,lg,3,1,2:3),5);
  
        svec(cll,lg,1) = stillcondfr(cll,lg,1,1,1); svecerr(cll,lg,1) = stillconderr(cll,lg,1,1,1);
        svec(cll,lg,2) = nanmean(stillcondfr(cll,lg,1,2:3,1),4); svecerr(cll,lg,2) = nanmean(stillconderr(cll,lg,1,2:3,1),4);
        svec(cll,lg,3) = nanmean(stillcondfr(cll,lg,2,2:3,1),4); svecerr(cll,lg,3) = nanmean(stillconderr(cll,lg,2,2:3,1),4);
        svec(cll,lg,4) = nanmean(stillcondfr(cll,lg,3,2:3,1),4); svecerr(cll,lg,4) = nanmean(stillconderr(cll,lg,3,2:3,1),4);
        svec(cll,lg,5) = nanmean([stillcondfr(cll,lg,1,2,2),stillcondfr(cll,lg,1,3,3)]); svecerr(cll,lg,5) = nanmean([stillconderr(cll,lg,1,2,2),stillconderr(cll,lg,1,3,3)]);
        svec(cll,lg,6) = nanmean([stillcondfr(cll,lg,1,2,3),stillcondfr(cll,lg,1,3,2)]); svecerr(cll,lg,6) = nanmean([stillconderr(cll,lg,1,2,3),stillconderr(cll,lg,1,3,2)]);
        svec(cll,lg,7) = nanmean([stillcondfr(cll,lg,2,2,3),stillcondfr(cll,lg,2,3,2)]); svecerr(cll,lg,7) = nanmean([stillconderr(cll,lg,2,2,3),stillconderr(cll,lg,2,3,2)]);
        svec(cll,lg,8) = nanmean([stillcondfr(cll,lg,3,2,3),stillcondfr(cll,lg,3,3,2)]); svecerr(cll,lg,8) = nanmean([stillconderr(cll,lg,3,2,3),stillconderr(cll,lg,3,3,2)]);
        svec(cll,lg,9) = nanmean(stillcondfr(cll,lg,1,1,2:3),5); svecerr(cll,lg,9) = nanmean(stillconderr(cll,lg,1,1,2:3),5);
        svec(cll,lg,10) = nanmean(stillcondfr(cll,lg,2,1,2:3),5); svecerr(cll,lg,10) = nanmean(stillconderr(cll,lg,2,1,2:3),5);
        svec(cll,lg,11) = nanmean(stillcondfr(cll,lg,3,1,2:3),5); svecerr(cll,lg,11) = nanmean(stillconderr(cll,lg,3,1,2:3),5);
  
        ppcvec(cll,lg,1,:) =  squeeze(ftppc(cll,lg,1,1,1,:)); 
        ppcvec(cll,lg,2,:) = nanmean(ftppc(cll,lg,1,2:3,1,:),4);
        ppcvec(cll,lg,3,:) = nanmean(ftppc(cll,lg,2,2:3,1,:),4);
        ppcvec(cll,lg,4,:) = nanmean(ftppc(cll,lg,3,2:3,1,:),4);
        ppcvec(cll,lg,5,:) = nanmean([ftppc(cll,lg,1,2,2,:),ftppc(cll,lg,1,3,3,:)]);
        ppcvec(cll,lg,6,:) = nanmean([ftppc(cll,lg,1,2,3,:),ftppc(cll,lg,1,3,2,:)]); 
        ppcvec(cll,lg,7,:) = nanmean([ftppc(cll,lg,2,2,3,:),ftppc(cll,lg,2,3,2,:)]);
        ppcvec(cll,lg,8,:) = nanmean([ftppc(cll,lg,3,2,3,:),ftppc(cll,lg,3,3,2,:)]); 
        ppcvec(cll,lg,9,:) = nanmean(ftppc(cll,lg,1,1,2:3,:),5);
        ppcvec(cll,lg,10,:) = nanmean(ftppc(cll,lg,2,1,2:3,:),5);
        ppcvec(cll,lg,11,:) = nanmean(ftppc(cll,lg,3,1,2:3,:),5);
    end
    
%     figure
%     errorbar(squeeze(vec(cll,1,:)),squeeze(vecerr(cll,1,:)),'k.')
%     set(gca,'xtick',1:11)
%     set(gca,'xticklabel',{'BL','sm','md','lg','bgonly','smbg','mdbg','lgbg','smap','mdap','lgap'})
        
    figure
    errorbar(squeeze(vec(cll,1,:)),squeeze(vecerr(cll,1,:)),'k.')
    hold on    
    errorbar(squeeze(vec(cll,2,:)),squeeze(vecerr(cll,2,:)),'r.')
    set(gca,'xtick',1:11)
    set(gca,'xticklabel',{'BL','sm','md','lg','bgonly','smbg','mdbg','lgbg','smap','mdap','lgap'})
    legend('light off','light on')
    
    figure
    errorbar(1:3,squeeze(vec(cll,1,2:4)),squeeze(vecerr(cll,1,2:4)),'k.-','linewidth',2);
    hold on
    errorbar(1:3,squeeze(vec(cll,1,6:8)),squeeze(vecerr(cll,1,6:8)),'r.-','linewidth',2);
    errorbar(1:3,squeeze(vec(cll,1,9:11)),squeeze(vecerr(cll,1,9:11)),'color',[.75,.75,.75]);
    errorbar(0,squeeze(vec(cll,1,1)),squeeze(vecerr(cll,1,1)),'k.','linewidth',2);
    errorbar(0,squeeze(vec(cll,1,5)),squeeze(vecerr(cll,1,5)),'r.','linewidth',2);
    ax = axis;
    axis([-.5,3.5,ax(3),ax(4)]);
    xlabel('sizes')
    set(gca,'xtick',0:1:3);
    set(gca,'xticklabel',{'BG','8','20','40'})
    legend('gray background','drifting background','aperture')
    
%     figure
%     plot(bta, nanmean([squeeze(bincondresp(cll,1,1,2,2,:)),squeeze(bincondresp(cll,1,1,3,3,:))],2),'k','linewidth',2);
%     hold on
%     plot(bta, nanmean([squeeze(bincondresp(cll,1,1,2,3,:)),squeeze(bincondresp(cll,1,1,3,2,:))],2),'b','linewidth',2);
%     plot(bta, nanmean([squeeze(bincondresp(cll,1,2,2,3,:)),squeeze(bincondresp(cll,1,2,3,2,:))],2),'c','linewidth',2);
%     plot(bta, nanmean([squeeze(bincondresp(cll,1,3,2,3,:)),squeeze(bincondresp(cll,1,3,3,2,:))],2),'g','linewidth',2);
%     legend('BG','SM','MD','LG')
%     title('with background')
%     
%     figure
%     plot(bta, squeeze(bincondresp(cll,1,1,1,1,:)),'k','linewidth',2);
%     hold on
%     plot(bta, squeeze(nanmean(bincondresp(cll,1,1,2:3,1,:),4)),'b','linewidth',2);
%     plot(bta, squeeze(nanmean(bincondresp(cll,1,2,2:3,1,:),4)),'c','linewidth',2);
%     plot(bta, squeeze(nanmean(bincondresp(cll,1,3,2:3,1,:),4)),'g','linewidth',2);
%     legend('BL','SM','MD','LG')
%     title('no background')
    
%     figure
%     subplot(3,1,1)
%     plot(bta, squeeze(nanmean(bincondresp(cll,1,1,2:3,1,:),4)),'b','linewidth',2);
%     hold on
%     plot(bta, nanmean([squeeze(bincondresp(cll,1,1,2,3,:)),squeeze(bincondresp(cll,1,1,3,2,:))],2),'r','linewidth',2);
%     legend('gray BG','with BG')
%     ylabel('small')
%     
%     subplot(3,1,2)
%     plot(bta, squeeze(nanmean(bincondresp(cll,1,2,2:3,1,:),4)),'b','linewidth',2);
%     hold on
%     plot(bta, nanmean([squeeze(bincondresp(cll,1,2,2,3,:)),squeeze(bincondresp(cll,1,2,3,2,:))],2),'r','linewidth',2);
%     ylabel('medium')
%     
%     subplot(3,1,3)    
%     plot(bta, squeeze(nanmean(bincondresp(cll,1,3,2:3,1,:),4)),'b','linewidth',2);
%     hold on
%     plot(bta, nanmean([squeeze(bincondresp(cll,1,3,2,3,:)),squeeze(bincondresp(cll,1,3,3,2,:))],2),'r','linewidth',2);
%     ylabel('large')
    
    figure
    subplot(2,2,1)
    imagesc(xax,yax,on)
    hold on
    plot_orrf_absdeg(gaussfiton,1,'w',2)
    axis square
    plot_circle(result.position(1),result.position(2),min(sizes)/2,'k',2)
    plot_circle(result.position(1),result.position(2),sizes(2)/2,'k',2)
    plot_circle(result.position(1),result.position(2),sizes(3)/2,'k',2)
    
    subplot(2,2,2)
    imagesc(xax,yax,off)
    hold on
    plot_orrf_absdeg(gaussfitoff,1,'w',2)
    axis square
    plot_circle(result.position(1),result.position(2),min(sizes)/2,'k',2)
    plot_circle(result.position(1),result.position(2),sizes(2)/2,'k',2)
    plot_circle(result.position(1),result.position(2),sizes(3)/2,'k',2)
    
    subplot(2,2,3)
    imagesc(xax,yax,on)
    hold on
    plot_orrf_absdeg(gaussfiton,1,'w',2)
    axis square
    plot_circle(result.position(1),result.position(2),min(sizes)/2,'k',2)
    plot_circle(result.position(1),result.position(2),sizes(2)/2,'k',2)
    plot_circle(result.position(1),result.position(2),sizes(3)/2,'k',2)
    axis([-50,50,-30,30])
    axis equal
    
    subplot(2,2,4)
    imagesc(xax,yax,off)
    hold on
    plot_orrf_absdeg(gaussfitoff,1,'w',2)
    axis square
    plot_circle(result.position(1),result.position(2),min(sizes)/2,'k',2)
    plot_circle(result.position(1),result.position(2),sizes(2)/2,'k',2)
    plot_circle(result.position(1),result.position(2),sizes(3)/2,'k',2)
    axis([-50,50,-30,30])
    axis equal
    
%     cllvec(cll,:,:) = vec;
%     cllvecerr(cll,:,:) = vecerr;
    
    cll = cll + 1;
    disp('');
    
end
bincondresp = bincondresp.*(1000/binwidth);
binconderr = binconderr.*(1000/binwidth);

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


wo = 60/(1000/2); bw = wo/35;
[b,a] = iirnotch(wo,bw);

params.Fs = 1000; params.trialave = 1; params.err = [2 .05]; params.tapers = [5,9];
cell1 = 19; cell2 = 6;
for i = 1:size(clllfpresp,2)
    flbd1(i,:) = filtfilt(b,a,squeeze(cllbdlfpresp(cell1,i,:)));
    flbd2(i,:) = filtfilt(b,a,squeeze(cllbdlfpresp(cell2,i,:)));
    fl1(i,:) = filtfilt(b,a,squeeze(clllfpresp(cell1,i,:)));
    fl2(i,:) = filtfilt(b,a,squeeze(clllfpresp(cell2,i,:)));
    
%     fl1(i,:) = squeeze(cllbdlfpresp(cell1,i,:));
%     fl2(i,:) = squeeze(cllbdlfpresp(cell2,i,:));
%     [coh(i,:),cfxx] = mscohere(fl1(i,respwin),fl2(i,respwin),[],[],512,1000);
end
for lg = 1:length(unique(light))
    for sz = 1:length(sizes)
        for oc = 1:length(orisc)
            for os = 1:length(oriss)
                trials = find(gratingInfo.sizecenter == sizes(sz) &...
                    gratingInfo.Orientation_center == orisc(oc) & ...
                    gratingInfo.Orientation_surround == oriss(os) & light == lg-1);
                if ~isempty(trials)
                    [p1(lg,sz,oc,os,:),chf,p1err(lg,sz,oc,os,:,:)] = mtspectrumc(fl1(trials,respwin)',params);
                    [p2(lg,sz,oc,os,:),chf,p2err(lg,sz,oc,os,:,:)] = mtspectrumc(fl2(trials,respwin)',params);
                    [C(lg,sz,oc,os,:),phi(lg,sz,oc,os,:),S12(lg,sz,oc,os,:),S1(lg,sz,oc,os,:),S2(lg,sz,oc,os,:),...
                        cfx,confC(lg,sz,oc,os),phistd(lg,sz,oc,os,:),Cerr(lg,sz,oc,os,:,:)] = coherencyc(fl1(trials,respwin)',fl2(trials,respwin)',params);    
                    [Cbd(lg,sz,oc,os,:),phibd(lg,sz,oc,os,:),S12bd(lg,sz,oc,os,:),S1bd(lg,sz,oc,os,:),S2bd(lg,sz,oc,os,:),...
                        cfx,confCbd(lg,sz,oc,os),phistdbd(lg,sz,oc,os,:),Cerrbd(lg,sz,oc,os,:,:)] = coherencyc(flbd1(trials,respwin)',flbd2(trials,respwin)',params);
                   
                    % field trip causality
                    for i = 1:length(trials)
                        data.trial{i} = [flbd1(trials(i),respwin);flbd2(trials(i),respwin)];
                        time = 1:1000;
%                         time = -prestim+1:result.stimduration*1000+poststim;
                        data.time{i} = time(1)./1000:.001:time(end)./1000;
                        %                 data.time{i} = -.299:.001:2.700;
                    end
                    data.fsample = 1000;
                    data.trialinfo(:,1) = light(trials);
                    data.trialinfo(:,2) = gratingInfo.sizecenter(trials);
                    data.trialinfo(:,3) = gratingInfo.Orientation_center(trials);
                    cfg = [];
                    data.label{1,1} = 'v1';
                    data.label{2,1} = 'v2';
                    cfg.latency = [0.5,1.5];
                              
                    % non-parametric
                    cfg           = [];
                    cfg.method    = 'mtmfft';
                    cfg.taper     = 'dpss';
                    cfg.output    = 'fourier';
                    cfg.tapsmofrq = 2;
                    freq          = ft_freqanalysis(cfg, data);
                    % coherence for both parametric and non-parametirc
                    cfg           = [];
                    cfg.method    = 'coh';
                    coh           = ft_connectivityanalysis(cfg, freq);
                    ftcoh(lg,sz,oc,os,:) = squeeze(coh.cohspctrm(1,2,:));
                    
                    cfg           = [];
                    cfg.method    = 'granger';
                    granger       = ft_connectivityanalysis(cfg, freq);
                    ftgcv2v1(lg,sz,oc,os,:) = squeeze(granger.grangerspctrm(2,1,:));
                    ftgcv1v2(lg,sz,oc,os,:) = squeeze(granger.grangerspctrm(1,2,:));
                                  
                else
                    p1(lg,sz,oc,os,:) = nan(1,513);
                    p1err(lg,sz,oc,os,:,:) = nan(2,513);
                    C(lg,sz,oc,os,:) = nan(1,513);
                    phi(lg,sz,oc,os,:) = nan(1,513);
                    S12(lg,sz,oc,os,:) = nan(1,513);
                    S1(lg,sz,oc,os,:) = nan(1,513);
                    S2(lg,sz,oc,os,:) = nan(1,513);
                    confC(lg,sz,oc,os) = NaN;
                    phistd(lg,sz,oc,os,:) = nan(1,513);
                    Cerr(lg,sz,oc,os,:,:) = nan(2,513);
                    Cbd(lg,sz,oc,os,:) = nan(1,513);
                    phibd(lg,sz,oc,os,:) = nan(1,513);
                    S12bd(lg,sz,oc,os,:) = nan(1,513);
                    S1bd(lg,sz,oc,os,:) = nan(1,513);
                    S2bd(lg,sz,oc,os,:) = nan(1,513);
                    confCbd(lg,sz,oc,os) = NaN;
                    phistdbd(lg,sz,oc,os,:) = nan(1,513);
                    Cerrbd(lg,sz,oc,os,:,:) = nan(2,513);
                    ftcoh(lg,sz,oc,os,:) = nan(1,501);
                    ftgcv1v2(lg,sz,oc,os,:) = nan(1,501);
                    ftgcv2v1(lg,sz,oc,os,:) = nan(1,501);
                end
            end
        end
    end
end
for l = 1:2
    bd1grayspect(l,:) = p1(l,1,1,1,:);
    bd1lgonlyspect(l,:) = nanmean(p1(l,3,2:3,1,:),3);
    bd1mdonlyspect(l,:) = nanmean(p1(l,2,2:3,1,:),3);
    bd1smonlyspect(l,:) = nanmean(p1(l,1,2:3,1,:),3);
    bd1lgonbgspect(l,:) = nanmean([squeeze(p1(l,3,2,3,:)),squeeze(p1(l,3,3,2,:))],2);
    bd1mdonbgspect(l,:) = nanmean([squeeze(p1(l,2,2,3,:)),squeeze(p1(l,2,3,2,:))],2);
    bd1smonbgspect(l,:) = nanmean([squeeze(p1(l,1,2,3,:)),squeeze(p1(l,1,3,2,:))],2);
    bd1bgonlyspect(l,:) = nanmean([squeeze(p1(l,1,2,2,:)),squeeze(p1(l,1,3,3,:))],2);
    bd2grayspect(l,:) = p2(l,1,1,1,:);
    bd2lgonlyspect(l,:) = nanmean(p2(l,3,2:3,1,:),3);
    bd2mdonlyspect(l,:) = nanmean(p2(l,2,2:3,1,:),3);
    bd2smonlyspect(l,:) = nanmean(p2(l,1,2:3,1,:),3);
    bd2lgonbgspect(l,:) = nanmean([squeeze(p2(l,3,2,3,:)),squeeze(p2(l,3,3,2,:))],2);
    bd2mdonbgspect(l,:) = nanmean([squeeze(p2(l,2,2,3,:)),squeeze(p2(l,2,3,2,:))],2);
    bd2smonbgspect(l,:) = nanmean([squeeze(p2(l,1,2,3,:)),squeeze(p2(l,1,3,2,:))],2);
    bd2bgonlyspect(l,:) = nanmean([squeeze(p2(l,1,2,2,:)),squeeze(p2(l,1,3,3,:))],2);
     
    spontcoher(l,:) = squeeze(Cbd(l,1,1,1,:));
    lgonlycoher(l,:) = nanmean(Cbd(l,3,2:3,1,:),3);
    mdonlycoher(l,:) = nanmean(Cbd(l,2,2:3,1,:),3);
    smonlycoher(l,:) = nanmean(Cbd(l,1,2:3,1,:),3);
    lgonbgcoher(l,:) = nanmean([squeeze(Cbd(l,3,2,3,:)),squeeze(Cbd(l,3,3,2,:))],2);
    mdonbgcoher(l,:) = nanmean([squeeze(Cbd(l,2,2,3,:)),squeeze(Cbd(l,2,3,2,:))],2);
    smonbgcoher(l,:) = nanmean([squeeze(Cbd(l,1,2,3,:)),squeeze(Cbd(l,1,3,2,:))],2);
    bgonlycoher(l,:) = nanmean([squeeze(Cbd(l,1,2,2,:)),squeeze(Cbd(l,1,3,3,:))],2);
    
    spontcoherft(l,:) = squeeze(ftcoh(l,1,1,1,:));
    lgonlycoherft(l,:) = nanmean(ftcoh(l,3,2:3,1,:),3);
    mdonlycoherft(l,:) = nanmean(ftcoh(l,2,2:3,1,:),3);
    smonlycoherft(l,:) = nanmean(ftcoh(l,1,2:3,1,:),3);
    lgonbgcoherft(l,:) = nanmean([squeeze(ftcoh(l,3,2,3,:)),squeeze(ftcoh(l,3,3,2,:))],2);
    mdonbgcoherft(l,:) = nanmean([squeeze(ftcoh(l,2,2,3,:)),squeeze(ftcoh(l,2,3,2,:))],2);
    smonbgcoherft(l,:) = nanmean([squeeze(ftcoh(l,1,2,3,:)),squeeze(ftcoh(l,1,3,2,:))],2);
    bgonlycoherft(l,:) = nanmean([squeeze(ftcoh(l,1,2,2,:)),squeeze(ftcoh(l,1,3,3,:))],2);
    
    spontgcff(l,:) = squeeze(ftgcv1v2(l,1,1,1,:));
    lgonlygcff(l,:) = nanmean(ftgcv1v2(l,3,2:3,1,:),3);
    mdonlygcff(l,:) = nanmean(ftgcv1v2(l,2,2:3,1,:),3);
    smonlygcff(l,:) = nanmean(ftgcv1v2(l,1,2:3,1,:),3);
    lgonbggcff(l,:) = nanmean([squeeze(ftgcv1v2(l,3,2,3,:)),squeeze(ftgcv1v2(l,3,3,2,:))],2);
    mdonbggcff(l,:) = nanmean([squeeze(ftgcv1v2(l,2,2,3,:)),squeeze(ftgcv1v2(l,2,3,2,:))],2);
    smonbggcff(l,:) = nanmean([squeeze(ftgcv1v2(l,1,2,3,:)),squeeze(ftgcv1v2(l,1,3,2,:))],2);
    bgonlygcff(l,:) = nanmean([squeeze(ftgcv1v2(l,1,2,2,:)),squeeze(ftgcv1v2(l,1,3,3,:))],2);
    
    spontgcfb(l,:) = squeeze(ftgcv2v1(l,1,1,1,:));
    lgonlygcfb(l,:) = nanmean(ftgcv2v1(l,3,2:3,1,:),3);
    mdonlygcfb(l,:) = nanmean(ftgcv2v1(l,2,2:3,1,:),3);
    smonlygcfb(l,:) = nanmean(ftgcv2v1(l,1,2:3,1,:),3);
    lgonbggcfb(l,:) = nanmean([squeeze(ftgcv2v1(l,3,2,3,:)),squeeze(ftgcv2v1(l,3,3,2,:))],2);
    mdonbggcfb(l,:) = nanmean([squeeze(ftgcv2v1(l,2,2,3,:)),squeeze(ftgcv2v1(l,2,3,2,:))],2);
    smonbggcfb(l,:) = nanmean([squeeze(ftgcv2v1(l,1,2,3,:)),squeeze(ftgcv2v1(l,1,3,2,:))],2);
    bgonlygcfb(l,:) = nanmean([squeeze(ftgcv2v1(l,1,2,2,:)),squeeze(ftgcv2v1(l,1,3,3,:))],2);
    
    
    lgonlycoherr(l,:,:) = nanmean(Cerr(l,3,2:3,1,:,:),3);
end

% try granger causality - see mvgc_demo

ntrials   = 210;     % number of trials
nobs      = 1000;   % number of observations per trial
regmode   = 'OLS';  % VAR model estimation regression mode ('OLS', 'LWR' or empty for default)
icregmode = 'LWR';  % information criteria regression mode ('OLS', 'LWR' or empty for default)
morder    = 'AIC';  % model order to use ('actual', 'AIC', 'BIC' or supplied numerical value)
momax     = 20;     % maximum model order for model order estimation
acmaxlags = 1000;   % maximum autocovariance lags (empty for automatic calculation)
tstat     = '';     % statistical test for MVGC:  'F' for Granger's F-test (default) or 'chi2' for Geweke's chi2 test
alpha     = 0.05;   % significance level for significance test
mhtc      = 'FDR';  % multiple hypothesis test correction (see routine 'significance')
fs        = 1000;    % sample rate (Hz)
fres      = [];     % frequency resolution (empty for automatic calculation)

X(1,:,:) = fl1(trials,:);
X(2,:,:) = fl2(trials,:);
[AIC,BIC,moAIC,moBIC] = tsdata_to_infocrit(X,momax,icregmode);

figure(1); clf;
plot_tsdata([AIC BIC]',{'AIC','BIC'},1/fs);
title('Model order estimation');

fprintf('\nbest model order (AIC) = %d\n',moAIC);
fprintf('best model order (BIC) = %d\n',moBIC);
morder = 20; % TODO
[A,SIG] = tsdata_to_var(X,morder,regmode);
[G,info] = var_to_autocov(A,SIG,acmaxlags);
var_info(info,true);
F = autocov_to_pwcgc(G);
pval = mvgc_pval(F,morder,nobs,ntrials,1,1,nvars-2,tstat); % take careful note of arguments!
sig  = significance(pval,alpha,mhtc);

figure(2); clf;
subplot(1,3,1);
plot_pw(F);
title('Pairwise-conditional GC');
subplot(1,3,2);
plot_pw(pval);
title('p-values');
subplot(1,3,3);
plot_pw(sig);
title(['Significant at p = ' num2str(alpha)])

f = autocov_to_spwcgc(G,fres);

figure(3); clf;
plot_spw(f,fs);
freqs = linspace(0,fs/2,size(f,3))';

for i = 1:length(depth)
    for l = 1:2
        grayspect(i,l,:) = condlfpspect(i,l,1,1,1,:);
        lgonlyspect(i,l,:) = nanmean(condlfpspect(i,l,3,2:3,1,:),4);
        mdonlyspect(i,l,:) = nanmean(condlfpspect(i,l,2,2:3,1,:),4);
        smonlyspect(i,l,:) = nanmean(condlfpspect(i,l,1,2:3,1,:),4);
        lgonbgspect(i,l,:) = nanmean([squeeze(condlfpspect(i,l,3,2,3,:)),squeeze(condlfpspect(i,l,3,3,2,:))],2);
        mdonbgspect(i,l,:) = nanmean([squeeze(condlfpspect(i,l,2,2,3,:)),squeeze(condlfpspect(i,l,2,3,2,:))],2);
        smonbgspect(i,l,:) = nanmean([squeeze(condlfpspect(i,l,1,2,3,:)),squeeze(condlfpspect(i,l,1,3,2,:))],2);
        bgonlyspect(i,l,:) = nanmean([squeeze(condlfpspect(i,l,1,2,2,:)),squeeze(condlfpspect(i,l,1,3,3,:))],2);
        
        
        lgonlybinresp(i,l,:) = nanmean(bincondresp(i,l,3,2:3,1,:),4);
        mdonlybinresp(i,l,:) = nanmean(bincondresp(i,l,2,2:3,1,:),4);
        smonlybinresp(i,l,:) = nanmean(bincondresp(i,l,1,2:3,1,:),4);
        bgonlybinresp(i,l,:) = nanmean([squeeze(bincondresp(i,l,1,2,2,:)),squeeze(bincondresp(i,l,1,3,3,:))],2);
        lgfigbinresp(i,l,:) = nanmean([squeeze(bincondresp(i,l,3,2,3,:)),squeeze(bincondresp(i,l,3,3,2,:))],2);
        mdfigbinresp(i,l,:) = nanmean([squeeze(bincondresp(i,l,2,2,3,:)),squeeze(bincondresp(i,l,2,3,2,:))],2);
        smfigbinresp(i,l,:) = nanmean([squeeze(bincondresp(i,l,1,2,3,:)),squeeze(bincondresp(i,l,1,3,2,:))],2);
        
        bgonlybinerr(i,l,:) = nanmean([squeeze(binconderr(i,l,1,2,2,:)),squeeze(binconderr(i,l,1,3,3,:))],2);
        lgfigbinerr(i,l,:) = nanmean([squeeze(binconderr(i,l,3,2,3,:)),squeeze(binconderr(i,l,3,3,2,:))],2);
        
        bgonlyresp(i,l,:) = nanmean([squeeze(condresp(i,l,1,2,2,:)),squeeze(condresp(i,l,1,3,3,:))],2);
        lgfigresp(i,l,:) = nanmean([squeeze(condresp(i,l,3,2,3,:)),squeeze(condresp(i,l,3,3,2,:))],2);
        mdfigresp(i,l,:) = nanmean([squeeze(condresp(i,l,2,2,3,:)),squeeze(condresp(i,l,2,3,2,:))],2);
    end
end

plot(bta,mean(bgonlybinresp))
hold on
plot(bta,mean(lgfigbinresp),'r')

for i = 1:length(depth)
    figure
    plot(bta,bgonlybinresp(i,:),'b')
    hold on
    plot(bta,lgfigbinresp(i,:),'r')
end

for i = 1:length(depth)
    figure
    plot(bta,bgonlybinresp(i,:),'b')
    hold on
    plot(bta,mdfigbinresp(i,:),'r')
end

disp('');