function positionalysis_SU_allcells

animalid = '170414';
block = 2;
lcol = 'b'; %lasercolor

onlymod = 0;
printyn = 1;
sfc = 0;

addpath C:\Users\Julia\work\Matlab\others\fieldtrip-20160329
ft_defaults

supath = ['C:\Users\Julia\work\data\' animalid '\singleunits\'];
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
    
    prestim = 0;
    poststim = 0;    
    respwin = 1500:2250;
    respwin = respwin-prestim;
    blwin = 1:1000;
    
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
    
    trialdur = result.sweeplength;
    msstamps = result.msstamps;
    
%     if strcmp(animalid, '170414')
%         randconds = result.randconds;
%         randconds1 = randconds(:,1:85);
%         randconds2 = randconds(:,86:end);
%         newrandconds = [randconds1,[0;0],randconds2];
%         randconds1 = newrandconds(:,1:171);
%         randconds2 = newrandconds(:,172:end);
%         newrandconds2 = [randconds1,[0;0],randconds2];
%         randconds1 = newrandconds2(:,1:174);
%         randconds2 = newrandconds2(:,175:end);
%         newrandconds3 = [randconds1,[0;0],randconds2];
%         randconds1 = newrandconds3(:,1:191);
%         randconds2 = newrandconds3(:,192:end);
%         newrandconds4 = [randconds1,[0;0],randconds2];
%         randconds1 = newrandconds4(:,1:225);
%         randconds2 = newrandconds4(:,226:end);        
%         newrandconds5 = [randconds1,[0;0],randconds2];
%         randconds1 = newrandconds5(:,1:233);
%         randconds2 = newrandconds5(:,234:end);
%         newrandconds6 = [randconds1,[0;0],randconds2];
%         randconds1 = newrandconds6(:,1:428);
%         randconds2 = newrandconds6(:,429:end);
%         newrandconds7 = [randconds1,[0;0],randconds2];
%         newrandconds8 = newrandconds7(:,1:450);
%         if fi ~= 1
%             result.randconds = newrandconds8;
%             save([supath, files(fi).name],'result');
%         end
%     end
    if strcmp(animalid, '170328')
        result.randconds = result.randconds(:,1:396);
        msstamps = msstamps(1:396);
    end
    if length(msstamps)~=size(result.randconds,2)
%         disp('');
% %                 msstamps([132]) = []; % for 170516 block 5
% %                 msstamps([49]) = []; % for 170601 block 3
%                 msstamps([38]) = []; % for 170518 block 7
%                 result.msstamps = msstamps;
%                 save([supath, files(fi).name],'result');
        pause;
    end
    
    
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
    end
    
    
    % figure out sufficiently high and nonvariable runspeed trials
    meanspeed = mean(speed(:,respwin),2);
    stdspeed = std(speed(:,respwin),1,2);
    notstill = find(meanspeed>1);
    okspeed = find(meanspeed>( mean(meanspeed(notstill))-(1.5*std(meanspeed(notstill))) ) & meanspeed>1);
    okvar = find(stdspeed<( mean(stdspeed(notstill))+(1.5*std(stdspeed(notstill)))) & stdspeed>.5);
    oktrials = intersect(okspeed,okvar);
    nonoktrials = 1:size(speed,1); nonoktrials(oktrials) = [];
    stilltrials = 1:size(speed,1); stilltrials(notstill) = [];
    okcond = zeros(1,size(resp,1)); okcond(oktrials) = 1;
    
    % important stuff
    depth(cll) = result.depth;    
    light = result.randconds(2,:); position = result.randconds(1,:);
    ta = 1:3000;
    
    spike = result.waveforms(:,wvchan);
    interpspike = spline(1:32,spike,1:.1:32);
    [adiff(cll),swidth(cll)] = spikequant(interpspike);
       
    msta = linspace(-blwin(end),trialdur-blwin(end),size(resp,2));
    
    frs = sum(resp(:,respwin),2)./(length(respwin)/1000);
    bl = sum(resp(:,blwin),2)./(length(blwin)/1000);
    sc = sum(resp(:,respwin),2); % spike count
        
     %determine if cll is touch modulated
%     touchfr = sum(resp(okcond&light == 0,blwin(end)+40:blwin(end)+40+200),2);
%     shortbl = sum(resp(okcond&light == 0,blwin(end)-40-200:blwin(end)-40),2);
    touchfr = sum(resp(okcond&light == 0,blwin(end)+151:blwin(end)+150+200),2);
    shortbl = sum(resp(okcond&light == 0,blwin(end)-149-200:blwin(end)-150),2);
    [p, touchmod(cll)] = signrank(shortbl,touchfr);
    
    %determine if cll is modulated by light
    lightmod(cll) = ttest2(frs(find(light)),frs(find(light == 0)));
    
    lfr(cll) = mean(frs(find(light)));
    nlfr(cll) = mean(frs(find(light == 0)));
    
    binwidth = 20;
    [binnedlight,bta] = binit(mean(resp(find(light),:)),binwidth);
    binnedlight = binnedlight.*(1000/binwidth);
    [binnednolight,bta] = binit(mean(resp(find(light==0),:)),binwidth);
    binnednolight = binnednolight.*(1000/binwidth);
    
    printname = files(fi).name;
    printname(find(printname=='_')) = ' ';
       
    clear stalfp; clear nspkssta; clear shufstalfp; clear nspksshufsta;
    for l = 1:length(unique(result.light))
        for pos = 1:length(result.positions)
            trials = find(position == result.positions(pos) & light == l-1);
            
            condresp(cll,l,pos,:) = mean(resp(trials,:),1);
            condresperr(cll,l,pos,:) = nanstd(resp(trials,:),1,1)./sqrt(length(trials));
            condlfpresp(cll,l,pos,:) = mean(lfpresp(trials,:),1);
            
            condlfpspect(cll,l,pos,:) = nanmean(lfpspect(trials,:));            
            [S,chf,Serr]=mtspectrumc(squeeze(lfpresp(trials,1750:2250))',params);
            condS(cll,l,pos,:) = S(1:60); condSerr(cll,l,pos,:,:) = Serr(:,1:60);
            
            condfr(cll,l,pos) = mean(frs(trials));%-mean(bl);
            conderr(cll,l,pos) =std(frs(trials))./sqrt(length(trials));
            
%             condz(cll,l,pos,:) = (sc(trials)-mean(sc(trials)))/std(sc(trials)); %ecker 2010
%             condsc(cll,l,pos,:) = sc(trials);
            ff(cll,l,pos) = var(sc(trials))/mean(sc(trials));
            
            hlp = ones(size(condresp,4),1);
            [xx,bta] = binit(hlp,binwidth);
            bincondresp(cll,l,pos,:) = binit(condresp(cll,l,pos,:),binwidth).*(1000/binwidth);
            binconderr(cll,l,pos,:) = binit(condresperr(cll,l,pos,:),binwidth).*(1000/binwidth);
            
            mscc = []; bincc = [];
            for ii = 1:length(trials)-1
                for jj = ii+1:length(trials)
                    help = corrcoef(resp(trials(ii),:),resp(trials(jj),:));
                    mscc = [mscc,help(1,2)];
                    help = corrcoef(binit(resp(trials(ii),:),binwidth),binit(resp(trials(jj),:),binwidth));
                    bincc = [bincc, help(1,2)];
                end
            end
            msreliab(cll,l,pos) = nanmean(mscc);
            binreliab(cll,l,pos) = nanmean(bincc);
            eckerreliability(cll,l,pos) = var(frs(trials))/var(frs);
            
            thisruninds = intersect(trials,oktrials);
            if ~isempty(thisruninds)
                thisrunn(cll,l,pos) = length(thisruninds);
                runcondresp(cll,l,pos,:) = mean(resp(thisruninds,:),1);
                runcondresperr(cll,l,pos,:) = nanstd(resp(thisruninds,:),1,1)./sqrt(length(thisruninds));
                runcondfr(cll,l,pos) = mean(frs(thisruninds));
                runconderr(cll,l,pos) = std(frs(thisruninds))./sqrt(length(thisruninds));
                runbincondresp(cll,l,pos,:) = binit(runcondresp(cll,l,pos,:),binwidth).*(1000/binwidth);
                runbinconderr(cll,l,pos,:) = binit(runcondresperr(cll,l,pos,:),binwidth).*(1000/binwidth);
                [runS,chf,runSerr] = mtspectrumc(squeeze(lfpresp(thisruninds,1750:2250))',params);
                condrunS(cll,l,pos,:) = runS(1:60); condrunSerr(cll,l,pos,:,:) = runSerr(:,1:60);  
            else
                thisrunn(cll,l,pos) = 0;
                runcondresp(cll,l,pos,:) = nan(1,size(resp,2));
                runcondresperr(cll,l,pos,:) = nan(1,size(resp,2));
                runcondfr(cll,l,pos) = NaN;
                runconderr(cll,l,pos) = NaN;
                runbincondresp(cll,l,pos,:) = nan(1,length(bta));
                runbinconderr(cll,l,pos,:) = nan(1,length(bta));
                condrunS(cll,l,pos,:) = nan(1,60); condrunSerr(cll,l,pos,:,:) = nan(2,60);
            end
            
            thisstillinds = intersect(trials,stilltrials);
            if ~isempty(thisstillinds)
                thisstilln(cll,l,pos) = length(thisstillinds);
                stillcondresp(cll,l,pos,:) = mean(resp(thisstillinds,:),1);
                stillcondresperr(cll,l,pos,:) = nanstd(resp(thisstillinds,:),1,1)./sqrt(length(thisstillinds));
                stillcondfr(cll,l,pos) = mean(frs(thisstillinds));
                stillconderr(cll,l,pos) = std(frs(thisstillinds))./sqrt(length(thisstillinds));
                stillbincondresp(cll,l,pos,:) = binit(stillcondresp(cll,l,pos,:),binwidth).*(1000/binwidth);
                stillbinconderr(cll,l,pos,:) = binit(stillcondresperr(cll,l,pos,:),binwidth).*(1000/binwidth);
                [stillS,chf,stillSerr] = mtspectrumc(squeeze(lfpresp(thisstillinds,1750:2250))',params);
                condstillS(cll,l,pos,:) = stillS(1:60); condstillSerr(cll,l,pos,:,:) = stillSerr(:,1:60);
            else
                thisstilln(cll,l,pos) = 0;
                stillcondresp(cll,l,pos,:) = nan(1,size(resp,2));
                stillcondresperr(cll,l,pos,:) = nan(1,size(resp,2));
                stillcondfr(cll,l,pos) = NaN;
                stillconderr(cll,l,pos) = NaN;
                stillbincondresp(cll,l,pos,:) = nan(1,length(bta));
                stillbinconderr(cll,l,pos,:) = nan(1,length(bta));
                condstillS(cll,l,pos,:) = nan(1,60); condstillSerr(cll,l,pos,:,:) = nan(2,60);
            end
            
            % field trip shitz fo SFC spectra
            clear data;
            for i = 1:length(trials)
                data.trial{i} = [lfpresp(trials(i),:);resp(trials(i),:)];
                time = 1:result.sweeplength;
                data.time{i} = time(1)./1000:.001:time(end)./1000;
%                 data.time{i} = -.299:.001:2.700;
            end
            data.fsample = 1000;
            data.trialinfo(:,1) = light(trials);
            data.trialinfo(:,2) = position(trials);
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
            ftspect(cll,l,pos,:) = squeeze(nanmean(mag(:,1,:)));
            ftphases{cll,l,pos} = squeeze(ang);
            ftfax = stsFFT.freq;
            
            cfg               = [];
            cfg.method        = 'ral'; % compute the rayleigh test
            cfg.spikechannel  = stsFFT.label{1};
            cfg.channel       = stsFFT.lfplabel; % selected LFP channels
            cfg.avgoverchan   = 'unweighted'; % weight spike-LFP phases irrespective of LFP power
            cfg.timwin        = 'all'; % compute over all available spikes in the window
            cfg.latency       = [0.5 1.5]; % sustained visual stimulation period
            statSts           = ft_spiketriggeredspectrum_stat(cfg,stsFFT);
            ftralp(cll,l,pos,:) = statSts.ral;
            cfg.method = 'ppc0';
            statSts           = ft_spiketriggeredspectrum_stat(cfg,stsFFT);
            ftppc(cll,l,pos,:) = statSts.ppc0;
            cfg.method = 'plv';
            statSts           = ft_spiketriggeredspectrum_stat(cfg,stsFFT);
            ftplv(cll,l,pos,:) = statSts.plv;
        end
    end
    sfx = chf(1:60);
    
       
    ta = bta-blwin(end);
    
    figure
    subplot(2,2,1)
    plot(ta,squeeze(mean(runbincondresp(cll,1,:,:),3)),'k','linewidth',2);
    hold on
    plot(ta,squeeze(mean(runbincondresp(cll,2,:,:),3)),lcol,'linewidth',2);
    mx = max([max(squeeze(mean(runbincondresp(cll,1,:,:),3))),max(squeeze(mean(runbincondresp(cll,2,:,:),3))),.01]);
    axis([-blwin(end),trialdur-blwin(end),0,mx]);
    line([0,0],[0,mx],'color','k','linewidth',2);
    line([1500,1500],[0,mx],'color','k','linewidth',2);
    line([500,500],[0,mx],'color','b','linewidth',2)
    line([1250,1250],[0,mx],'color','b','linewidth',2);
    legend({'Light OFF','Light ON'})
    xlabel('time [ms]')
    ylabel('firing rate [Hz]')
    title(['cll ' int2str(cll) ' depth: ' int2str(result.depth), 'cll ' printname ])
    
    subplot(2,2,2)
    errorbar(result.positions,squeeze(runcondfr(cll,1,:)),squeeze(runconderr(cll,1,:)),'ko-','markerfacecolor','k','linewidth',2)
    hold on
    errorbar(result.positions,squeeze(runcondfr(cll,2,:)),squeeze(runconderr(cll,2,:)),'ro-','markerfacecolor','r','linewidth',2)
    ax = axis;
    axis([-1,9,ax(3),ax(4)]);
    legend('light off','light on')
    xlabel('bar position')
    ylabel('firing rate [Hz]')
        
    subplot(2,2,3)
    plot(ta,squeeze(mean(stillbincondresp(cll,1,:,:),3)),'k','linewidth',2);
    hold on
    plot(ta,squeeze(mean(stillbincondresp(cll,2,:,:),3)),lcol,'linewidth',2);
    mx = max([max(squeeze(mean(stillbincondresp(cll,1,:,:),3))),max(squeeze(mean(stillbincondresp(cll,2,:,:),3))),.01]);
    axis([-blwin(end),trialdur-blwin(end),0,mx]);
    line([0,0],[0,mx],'color','k','linewidth',2);
    line([1500,1500],[0,mx],'color','k','linewidth',2);
    line([500,500],[0,mx],'color','b','linewidth',2)
    line([1250,1250],[0,mx],'color','b','linewidth',2);
    legend({'Light OFF','Light ON'})
    xlabel('time [ms]')
    ylabel('firing rate [Hz]')
    title('non-running')
    
    subplot(2,4,7)
    plot(spike)
    axis([0,40,-100,100])
    legend(['width: ' int2str(swidth(cll)) ' adiff: ' num2str(adiff(cll))])
    
    
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
    barweb([squeeze(mean(condfr(cll,1,:),3)),squeeze(mean(condfr(cll,2,:)));...
        squeeze(mean(runcondfr(cll,1,:),3)),squeeze(mean(runcondfr(cll,2,:)));...        
        squeeze(mean(stillcondfr(cll,1,:),3)),squeeze(mean(stillcondfr(cll,2,:)))],...
        [squeeze(mean(conderr(cll,1,:),3)),squeeze(mean(conderr(cll,2,:)));...
        squeeze(mean(runconderr(cll,1,:),3)),squeeze(mean(runconderr(cll,2,:)));...        
        squeeze(mean(stillconderr(cll,1,:),3)),squeeze(mean(stillconderr(cll,2,:)))],...
        [],[{'all'};{'running only'};{'immobile only'}],[],...
        [],'firing rate [Hz]',[],[]);
    
    % LFP figure
    fillx = [sfx,fliplr(sfx)];
    fillyrun0 = [squeeze(nanmean(condrunSerr(cll,1,:,1,:),3))',fliplr(squeeze(nanmean(condrunSerr(cll,1,:,2,:)))')];
    fillyrun1 = [squeeze(nanmean(condrunSerr(cll,2,:,1,:),3))',fliplr(squeeze(nanmean(condrunSerr(cll,2,:,2,:)))')];
    fillystill0 = [squeeze(nanmean(condstillSerr(cll,1,:,1,:),3))',fliplr(squeeze(nanmean(condstillSerr(cll,1,:,2,:)))')];
    
    figure
    
    subplot(2,2,1)
    fill(fillx,fillyrun0,'b')
    hold on
    fill(fillx,fillyrun1,'r')
    set(gca,'yscale','log')
    xlabel('frequency [Hz]')
    ylabel('power')
    legend('L0','L1')
    
    subplot(2,2,2)
    fill(fillx,fillyrun0,'c')
    hold on
    fill(fillx,fillystill0,'b')
    set(gca,'yscale','log')
    xlabel('frequency [Hz]')
    ylabel('power')
    legend('running','still')
    
    subplot(2,2,3)
    plot(ftfax,squeeze(nanmean(ftppc(cll,1,:,:),3)),'b')
    hold on    
    plot(ftfax,squeeze(nanmean(ftppc(cll,2,:,:),3)),'r')
    xlabel('frequency')
    ylabel('PPC')
    legend('L0','L1')   
    
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

prsv = swidthms>.36; pfsv = swidthms<=.36;
prs = find(prsv); pfs = find(pfsv);

pangle = 30;
depth = depth.*cosd(45-pangle);

for i = 1:length(depth)
    mfr(i) = mean(runcondfr(i,1,2:9),3);
    omi(i) = (nanmean(runcondfr(i,2,2:9),3)-nanmean(runcondfr(i,1,2:9),3))./(nanmean(runcondfr(i,2,2:9),3)+nanmean(runcondfr(i,1,2:9),3));
    ominr(i) = (nanmean(stillcondfr(i,2,2:9),3)-nanmean(stillcondfr(i,1,2:9),3))./(nanmean(stillcondfr(i,2,2:9),3)+nanmean(stillcondfr(i,1,2:9),3));
    deltaspikes(i) = nanmean(runcondfr(i,2,:),3)-nanmean(runcondfr(i,1,:),3);
    for l = 1:2
        ssi(i,l) = get_ssi(squeeze(runcondfr(i,l,2:8))); % spatial selectivity
    end
end



disp('');


function ssi = get_ssi(curve)
ssi  = 1 - (((norm(curve)/max(curve)) - 1)./((sqrt(length(curve)))-1));