function speedanalysis_SU_allcells

animalid = '170404';
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
    respwin = 2650:3350;
    respwin = respwin-prestim;
    strespwin = 2850:3050;
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
  
    if length(msstamps)~=size(result.randconds,2)
        disp('');
        %         msstamps([518]) = []; % for 160726 block 6
        %         result.msstamps = msstamps;
        %         save([supath, files(fi).name],'result');
        pause;
    end
    
    
    for i = 1:length(msstamps)
        resp(i,:) = chan(msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim); 
        lfpresp(i,:) = result.lfp(msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim, wvchan);
        filtlfpresp(i,:) = filtlfp(msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim);
        [lfpspect(i,:),trialfax] = pmtm(lfpresp(i,respwin(1)+100:respwin(end)),3,[],sr);
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
        
        runspeed(i,:) = result.runspeed(msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim);     
    end
    
    
    % figure out sufficiently high and nonvariable runspeed trials
    meanrunspeed = mean(runspeed(:,respwin),2);
    stdrunspeed = std(runspeed(:,respwin),1,2);
    notstill = find(meanrunspeed>1);
    okrunspeed = find(meanrunspeed>( mean(meanrunspeed(notstill))-(1.5*std(meanrunspeed(notstill))) ) & meanrunspeed>1);
    okvar = find(stdrunspeed<( mean(stdrunspeed(notstill))+(1.5*std(stdrunspeed(notstill)))) & stdrunspeed>.5);
    oktrials = intersect(okrunspeed,okvar);
    nonoktrials = 1:size(runspeed,1); nonoktrials(oktrials) = [];
    stilltrials = 1:size(runspeed,1); stilltrials(notstill) = [];
    
    % important stuff
    depth(cll) = result.depth;    
    light = result.randconds(2,:); barspeed = result.randconds(1,:);
    
    spike = result.waveforms(:,wvchan);
    interpspike = spline(1:32,spike,1:.1:32);
    [adiff(cll),swidth(cll)] = spikequant(interpspike);
       
    msta = linspace(-blwin(end),trialdur-blwin(end),size(resp,2));
    
    frs = sum(resp(:,respwin),2)./(length(respwin)/1000);
    stfrs = sum(resp(:,strespwin),2)./(length(strespwin)/1000);
    bl = sum(resp(:,blwin),2)./(length(blwin)/1000);
    sc = sum(resp(:,respwin),2); % spike count
        
     %determine if cll is touch modulated
    touchfr = sum(resp(:,strespwin),2);
    shortbl = sum(resp(:,blwin(end)-length(strespwin):blwin(end)-1),2);
    touchmod(cll) = ttest2(shortbl,touchfr);
    
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
        for sp = 1:length(result.speeds)
            trials = find(barspeed == result.speeds(sp) & light == l-1);
            
            condresp(cll,l,sp,:) = mean(resp(trials,:),1);
            condresperr(cll,l,sp,:) = nanstd(resp(trials,:),1,1)./sqrt(length(trials));
            condlfpresp(cll,l,sp,:) = mean(lfpresp(trials,:),1);
            
            condlfpspect(cll,l,sp,:) = nanmean(lfpspect(trials,:));            
            [S,chf,Serr]=mtspectrumc(squeeze(lfpresp(trials,respwin))',params);
            condS(cll,l,sp,:) = S(1:103); condSerr(cll,l,sp,:,:) = Serr(:,1:103);
            
            condfr(cll,l,sp) = mean(frs(trials));%-mean(bl);
            conderr(cll,l,sp) =std(frs(trials))./sqrt(length(trials));
            condstfr(cll,l,sp) = mean(stfrs(trials));
            condsterr(cll,l,sp) = std(frs(trials))./sqrt(length(trials));
            
            condz(cll,l,sp,:) = (sc(trials)-mean(sc(trials)))/std(sc(trials)); %ecker 2010
            condsc(cll,l,sp,:) = sc(trials);
            ff(cll,l,sp) = var(sc(trials))/mean(sc(trials));
            
            hlp = ones(size(condresp,4),1);
            [xx,bta] = binit(hlp,binwidth);
            bincondresp(cll,l,sp,:) = binit(condresp(cll,l,sp,:),binwidth).*(1000/binwidth);
            binconderr(cll,l,sp,:) = binit(condresperr(cll,l,sp,:),binwidth).*(1000/binwidth);
            
            mscc = []; bincc = [];
            for ii = 1:length(trials)-1
                for jj = ii+1:length(trials)
                    help = corrcoef(resp(trials(ii),:),resp(trials(jj),:));
                    mscc = [mscc,help(1,2)];
                    help = corrcoef(binit(resp(trials(ii),:),binwidth),binit(resp(trials(jj),:),binwidth));
                    bincc = [bincc, help(1,2)];
                end
            end
            msreliab(cll,l,sp) = nanmean(mscc);
            binreliab(cll,l,sp) = nanmean(bincc);
            eckerreliability(cll,l,sp) = var(frs(trials))/var(frs);
            
            thisruninds = intersect(trials,oktrials);
            if ~isempty(thisruninds)
                thisrunn(cll,l,sp) = length(thisruninds);
                runcondresp(cll,l,sp,:) = mean(resp(thisruninds,:),1);
                runcondresperr(cll,l,sp,:) = nanstd(resp(thisruninds,:),1,1)./sqrt(length(thisruninds));
                runcondfr(cll,l,sp) = mean(frs(thisruninds));
                runconderr(cll,l,sp) = std(frs(thisruninds))./sqrt(length(thisruninds));
                runcondstfr(cll,l,sp) = mean(stfrs(thisruninds));
                runcondsterr(cll,l,sp) = std(stfrs(thisruninds))./sqrt(length(thisruninds));
                runbincondresp(cll,l,sp,:) = binit(runcondresp(cll,l,sp,:),binwidth).*(1000/binwidth);
                runbinconderr(cll,l,sp,:) = binit(runcondresperr(cll,l,sp,:),binwidth).*(1000/binwidth);
                [runS,chf,runSerr] = mtspectrumc(squeeze(lfpresp(thisruninds,respwin))',params);
                condrunS(cll,l,sp,:) = runS(1:103); condrunSerr(cll,l,sp,:,:) = runSerr(:,1:103);  
            else
                thisrunn(cll,l,sp) = 0;
                runcondresp(cll,l,sp,:) = nan(1,size(resp,2));
                runcondresperr(cll,l,sp,:) = nan(1,size(resp,2));
                runcondfr(cll,l,sp) = NaN;
                runconderr(cll,l,sp) = NaN;
                runcondstfr(cll,l,sp) = NaN;
                runcondsterr(cll,l,sp) = NaN;
                runbincondresp(cll,l,sp,:) = nan(1,length(bta));
                runbinconderr(cll,l,sp,:) = nan(1,length(bta));
                condrunS(cll,l,sp,:) = nan(1,103); condrunSerr(cll,l,sp,:,:) = nan(2,103);
            end
            
            thisstillinds = intersect(trials,stilltrials);
            if ~isempty(thisstillinds)
                thisstilln(cll,l,sp) = length(thisstillinds);
                stillcondresp(cll,l,sp,:) = mean(resp(thisstillinds,:),1);
                stillcondresperr(cll,l,sp,:) = nanstd(resp(thisstillinds,:),1,1)./sqrt(length(thisstillinds));
                stillcondfr(cll,l,sp) = mean(frs(thisstillinds));
                stillconderr(cll,l,sp) = std(frs(thisstillinds))./sqrt(length(thisstillinds));
                stillcondstfr(cll,l,sp) = mean(stfrs(thisstillinds));
                stillcondsterr(cll,l,sp) = std(stfrs(thisstillinds))./sqrt(length(thisstillinds));
                stillbincondresp(cll,l,sp,:) = binit(stillcondresp(cll,l,sp,:),binwidth).*(1000/binwidth);
                stillbinconderr(cll,l,sp,:) = binit(stillcondresperr(cll,l,sp,:),binwidth).*(1000/binwidth);
                [stillS,chf,stillSerr] = mtspectrumc(squeeze(lfpresp(thisstillinds,respwin))',params);
                condstillS(cll,l,sp,:) = stillS(1:103); condstillSerr(cll,l,sp,:,:) = stillSerr(:,1:103);
            else
                thisstilln(cll,l,sp) = 0;
                stillcondresp(cll,l,sp,:) = nan(1,size(resp,2));
                stillcondresperr(cll,l,sp,:) = nan(1,size(resp,2));
                stillcondfr(cll,l,sp) = NaN;
                stillconderr(cll,l,sp) = NaN;
                stillcondstfr(cll,l,sp) = NaN;
                stillcondsterr(cll,l,sp) = NaN;
                stillbincondresp(cll,l,sp,:) = nan(1,length(bta));
                stillbinconderr(cll,l,sp,:) = nan(1,length(bta));
                condstillS(cll,l,sp,:) = nan(1,60); condstillSerr(cll,l,sp,:,:) = nan(2,103);
            end
            
            % field trip shitz fo SFC spectra
            for i = 1:length(trials)
                data.trial{i} = [lfpresp(trials(i),:);resp(trials(i),:)];
                time = 1:result.sweeplength;
                data.time{i} = time(1)./1000:.001:time(end)./1000;
%                 data.time{i} = -.299:.001:2.700;
            end
            data.fsample = 1000;
            data.trialinfo(:,1) = light(trials);
            data.trialinfo(:,2) = barspeed(trials);
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
            ftspect(cll,l,sp,:) = squeeze(nanmean(mag(:,1,:)));
            ftphases{cll,l,sp} = squeeze(ang);
            ftfax = stsFFT.freq;
            
            cfg               = [];
            cfg.method        = 'ral'; % compute the rayleigh test
            cfg.spikechannel  = stsFFT.label{1};
            cfg.channel       = stsFFT.lfplabel; % selected LFP channels
            cfg.avgoverchan   = 'unweighted'; % weight spike-LFP phases irrespective of LFP power
            cfg.timwin        = 'all'; % compute over all available spikes in the window
            cfg.latency       = [0.5 1.5]; % sustained visual stimulation period
            statSts           = ft_spiketriggeredspectrum_stat(cfg,stsFFT);
            ftralp(cll,l,sp,:) = statSts.ral;
            cfg.method = 'ppc0';
            statSts           = ft_spiketriggeredspectrum_stat(cfg,stsFFT);
            ftppc(cll,l,sp,:) = statSts.ppc0;
            cfg.method = 'plv';
            statSts           = ft_spiketriggeredspectrum_stat(cfg,stsFFT);
            ftplv(cll,l,sp,:) = statSts.plv;
        end
    end
    sfx = chf(1:103);
    
       
    ta = bta-blwin(end);
    
    figure
    subplot(2,2,1)
    plot(ta,squeeze(mean(runbincondresp(cll,1,:,:),3)),'k','linewidth',2);
    hold on
    plot(ta,squeeze(mean(runbincondresp(cll,2,:,:),3)),lcol,'linewidth',2);
    mx = max([max(squeeze(mean(runbincondresp(cll,1,:,:),3))),max(squeeze(mean(runbincondresp(cll,2,:,:),3))),.01]);
    axis([-blwin(end),trialdur-blwin(end),0,mx]);
    line([2650,2650],[0,mx],'color','b','linewidth',2)
    line([3350,3350],[0,mx],'color','b','linewidth',2);
    legend({'Light OFF','Light ON'})
    xlabel('time [ms]')
    ylabel('firing rate [Hz]')
    title(['cll ' int2str(cll) ' depth: ' int2str(result.depth), 'cll ' printname ])
    
    subplot(2,2,2)
    errorbar(result.speeds,squeeze(runcondstfr(cll,1,:)),squeeze(runcondsterr(cll,1,:)),'ko-','markerfacecolor','k','linewidth',2)
    hold on
    errorbar(result.speeds,squeeze(runcondstfr(cll,2,:)),squeeze(runcondsterr(cll,2,:)),'ro-','markerfacecolor','r','linewidth',2)
    ax = axis;
    axis([-.1,1.1,ax(3),ax(4)]);
    legend('light off','light on')
    xlabel('speed')
    ylabel('firing rate [Hz]')
        
    subplot(2,2,3)
    plot(ta,squeeze(mean(stillbincondresp(cll,1,:,:),3)),'k','linewidth',2);
    hold on
    plot(ta,squeeze(mean(stillbincondresp(cll,2,:,:),3)),lcol,'linewidth',2);
    mx = max([max(squeeze(mean(stillbincondresp(cll,1,:,:),3))),max(squeeze(mean(stillbincondresp(cll,2,:,:),3))),.01]);
    axis([-blwin(end),trialdur-blwin(end),0,mx]);
    line([2650,2650],[0,mx],'color','b','linewidth',2)
    line([3350,3350],[0,mx],'color','b','linewidth',2);
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
    imagesc(runspeed);
    colorbar
    title(['oktrials: ' int2str(length(oktrials)) '/' int2str(size(runspeed,1))])
    xlabel('time [ms]')
    ylabel('trial number')
    
    subplot(2,2,2)
    errorbar(msta,mean(runspeed(find(~light),:)),std(runspeed(find(~light),:))./sqrt(length(find(~light))),'b')
    hold on
    errorbar(msta,mean(runspeed(find(light),:)),std(runspeed(find(light),:))./sqrt(length(find(light))),'r')
    xlabel('time [ms]')
    ylabel('average runspeed')
    legend({'light off' 'light on'})
    
    subplot(2,2,3)
    plot(mean(runspeed(:,respwin),2),frs,'.')
    hold on
    plot(mean(runspeed(find(light),respwin),2),frs(find(light)),'r.')
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

for i = 1:length(runcondstfr)
    [sfr,ind] = sort(runcondstfr(i,1,:));
    rankfr(i,1,:) = runcondstfr(i,1,ind);
    rankfr(i,2,:) = runcondstfr(i,2,ind);
    normrankfr(i,1,:) = squeeze(rankfr(i,1,:))./max(squeeze(rankfr(i,1,:)));
    normrankfr(i,2,:) = squeeze(rankfr(i,2,:))./max(squeeze(rankfr(i,1,:)));
end

cond = prs;
figure
errorbar(nanmean(normrankfr(cond,1,:),1),nanmean(normrankfr(cond,2,:),1),nanstd(normrankfr(cond,2,:))./sqrt(length(cond)),'.')
hold on
herrorbar(nanmean(normrankfr(cond,1,:),1),nanmean(normrankfr(cond,2,:),1),nanstd(normrankfr(cond,1,:))./sqrt(length(cond)),'.')
refline(1,0)