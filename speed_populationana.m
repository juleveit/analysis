function speed_populationana



% X94 ChR2
animalids = {'170328','170330'};
blocks    = [ 2,       3];
animal    = [ 1,       2];
electrodes =[[1,32];  [1,32]];
penangle =  [ 30,      30];
printpath = 'C:\Users\Julia\work\data\populations\X94_ChR2_S1\speed\units\';
runprintpath = 'C:\Users\Julia\work\data\populations\X94_ChR2_S1\speed\running\';
lfpprintpath = 'C:\Users\Julia\work\data\populations\X94_ChR2_S1\speed\LFP\';
popfile = 'C:\Users\Julia\work\data\populations\X94_ChR2_S1\speed\speed_population.mat';

tic
lcol = 'r'; %lasercolor

recalculate = 0;
printyn = 1
sfc = 1;

%Scott Kernel
%%%%%%%%%%%%kernel propeties%%%%%%%%%%%%%%%
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
%%%%%%%%%%kernel properties%%%%%%%%%%%%%%%
%%%%% Fieldtrip
addpath C:\Users\Julia\work\Matlab\others\fieldtrip-20160329
ft_defaults
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% chronux parameters
params.tapers = [5,9]; params.Fs = 1000; params.err = [2, 0.05]; params.trialave = 1;

freqbinwidth = 5;
sr = 1000;
prestim = 0;
poststim = 0;
respwin = 2850:3050;
respwin = respwin+prestim;
lightwin = 2626:3375;
blwin = 1:1000;

if ~exist(popfile) || recalculate

    cll = 1;
    for blck = 1:length(blocks)

        supath = ['C:\Users\Julia\work\data\' animalids{blck} '\singleunits\'];
        basename = [animalids{blck} '_block' int2str(blocks(blck)) '_tet'];

        files = dir([supath, basename, '*.mat']);
        
        for fi = 1:length(files)

            if strfind(files(fi).name, 'MU')
                continue;
            end
            
            load([supath, files(fi).name]);            
            
            i = strfind(files(fi).name, 'tet');
            if strcmp(files(fi).name(i+4),'_')
                tetno = strread(files(fi).name(i+3)); % single character number
            else
                tetno = strread(files(fi).name(i+3:i+4)); % number >10
            end
            if tetno*4<electrodes(blck,1) || tetno*4>electrodes(blck,2) % assure we're only getting V1 in the population
                continue;
            end
            
            cllname{cll} = files(fi).name;
            printname = files(fi).name;
            printname(find(printname=='_')) = ' ';
                        
            % important stuff
            depth(cll) = result.depth;
            pangle(cll) = penangle(blck);
            recording(cll) = blck;
            animalno(cll) = animal(blck);
            
            light = result.randconds(2,:); barspeed = result.randconds(1,:);
            
            wvchan = find(var(result.waveforms) == max(var(result.waveforms)));
            lfp = result.lfp(:,wvchan)';
            spike = result.waveforms(:,wvchan);
            interpspike = spline(1:32,spike,1:.1:32);
            [adiff(cll),swidth(cll)] = spikequant(interpspike);
            waveform(cll,:) = spike;
            clustqual(cll) = result.clusterquality;

            % get spiketimes
            msStimes = round(result.spikes);
            if ~isempty(msStimes) && msStimes(1) == 0, msStimes(1) = 1; end

            chan = zeros(1,length(result.lfp));
            chan(msStimes) = 1;                
            
            gamma = eegfilt(lfp,sr,15,35);  % adjust depending on peak
            h = hilbert(gamma); gpow = abs(h); gphas = angle(h);
            % %
            if sfc
                clear filtmat; clear powmat; clear phasmat;
                for i = 1:100/freqbinwidth
                    filtmat(i,:) = eegfilt(lfp,sr,(i-1)*freqbinwidth+1,i*freqbinwidth);
                    h = hilbert(filtmat(i,:));
                    powmat(i,:) = abs(h); phasmat(i,:) = angle(h);
                end
            end
            
            trialdur = result.sweeplength;
            msstamps = result.msstamps;
            
            if length(msstamps)~=size(result.randconds,2)
                disp('');
                %         msstamps([518]) = []; 
                %         result.msstamps = msstamps;
                %         save([supath, files(fi).name],'result');
                pause;
            end
            
            for i = 1:length(msstamps)
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
            
            nrun(cll) = length(oktrials);
            nstill(cll) = length(stilltrials);
            if ~isempty(oktrials)
                cellmeanspeed(cll) = mean(meanrunspeed(oktrials));
                cellminspeed(cll) = min(meanrunspeed(oktrials));
            else
                cellmeanspeed(cll) = NaN;
                cellminspeed(cll) = NaN;
            end
            
            for i = 1:length(msstamps)
                resp(i,:) = chan(msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim);
                lfpresp(i,:) = lfp(msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim);               
                msl(i) = mean(find(resp(i,respwin)));
                [lfpspect(i,:),trialfax] = pmtm(lfpresp(i,lightwin),3,[],sr);
                [lfpsectrogram(i,:,:),ct,cf] = mtspecgramc(lfpresp(i,:)',[.25,.1],params);
                gammaresp(i,:) = gamma(msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim);
                gammapowresp(i,:) = gpow(msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim);
                gphaseresp(i,:) = gphas(msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim);
                
                % Spike triggered LFP
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
                
                % for all frequency bands
                if sfc
                    for j = 1:size(phasmat,1)
                        allphaseresp(j,i,:) = phasmat(j, msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim);
                        allpowresp(j,i,:) = powmat(j, msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim);
                    end
                end
            end            
                        
            frs = sum(resp(:,respwin),2)./(length(respwin)/1000);
            bl = sum(resp(:,blwin),2)./(length(blwin)/1000);
            sc = sum(resp(:,respwin),2); % spike count            
            
            %determine if cll is touch modulated
            touchfr = sum(resp(:,respwin),2);
            shortbl = sum(resp(:,blwin(end)-length(respwin):blwin(end)-1),2);
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
            
            clear stalfp; clear nspkssta; clear shufstalfp; clear nspksshufsta;
            for l = 1:length(unique(result.light))
                for sp = 1:length(result.speeds)
                    trials = find(barspeed == result.speeds(sp) & light == l-1);
                    
                    condresp(cll,l,sp,:) = mean(resp(trials,:),1);
                    condresperr(cll,l,sp,:) = nanstd(resp(trials,:),1,1)./sqrt(length(trials));
                    condlfpresp(cll,l,sp,:) = mean(lfpresp(trials,:),1);
                    
                    condlfpspect(cll,l,sp,:) = nanmean(lfpspect(trials,:));
                    [S,chf,Serr]=mtspectrumc(squeeze(lfpresp(trials,lightwin))',params);
                    condS(cll,l,sp,:) = S(1:103); condSerr(cll,l,sp,:,:) = Serr(:,1:103);
                    
                    condfr(cll,l,sp) = mean(frs(trials));%-mean(bl);
                    conderr(cll,l,sp) =std(frs(trials))./sqrt(length(trials));
                    
                    condz{cll,l,sp,:} = (sc(trials)-mean(sc(trials)))/std(sc(trials)); %ecker 2010
                    condsc{cll,l,sp,:} = sc(trials);
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
                        runbincondresp(cll,l,sp,:) = binit(runcondresp(cll,l,sp,:),binwidth).*(1000/binwidth);
                        runbinconderr(cll,l,sp,:) = binit(runcondresperr(cll,l,sp,:),binwidth).*(1000/binwidth);
                        [runS,chf,runSerr] = mtspectrumc(squeeze(lfpresp(thisruninds,lightwin))',params);
                        condrunS(cll,l,sp,:) = runS(1:103); condrunSerr(cll,l,sp,:,:) = runSerr(:,1:103);
                    else
                        thisrunn(cll,l,sp) = 0;
                        runcondresp(cll,l,sp,:) = nan(1,size(resp,2));
                        runcondresperr(cll,l,sp,:) = nan(1,size(resp,2));
                        runcondfr(cll,l,sp) = NaN;
                        runconderr(cll,l,sp) = NaN;
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
                        stillbincondresp(cll,l,sp,:) = binit(stillcondresp(cll,l,sp,:),binwidth).*(1000/binwidth);
                        stillbinconderr(cll,l,sp,:) = binit(stillcondresperr(cll,l,sp,:),binwidth).*(1000/binwidth);
                        [stillS,chf,stillSerr] = mtspectrumc(squeeze(lfpresp(thisstillinds,lightwin))',params);
                        condstillS(cll,l,sp,:) = stillS(1:103); condstillSerr(cll,l,sp,:,:) = stillSerr(:,1:103);
                    else
                        thisstilln(cll,l,sp) = 0;
                        stillcondresp(cll,l,sp,:) = nan(1,size(resp,2));
                        stillcondresperr(cll,l,sp,:) = nan(1,size(resp,2));
                        stillcondfr(cll,l,sp) = NaN;
                        stillconderr(cll,l,sp) = NaN;
                        stillbincondresp(cll,l,sp,:) = nan(1,length(bta));
                        stillbinconderr(cll,l,sp,:) = nan(1,length(bta));
                        condstillS(cll,l,sp,:) = nan(1,103); condstillSerr(cll,l,sp,:,:) = nan(2,103);
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
            
            %  figure
            clf;
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
            errorbar(result.speeds,squeeze(runcondfr(cll,1,:)),squeeze(runconderr(cll,1,:)),'ko-','markerfacecolor','k','linewidth',2)
            hold on
            errorbar(result.speeds,squeeze(runcondfr(cll,2,:)),squeeze(runconderr(cll,2,:)),'ro-','markerfacecolor','r','linewidth',2)
            ax = axis;
            axis([-1,9,ax(3),ax(4)]);
            legend('light off','light on')
            xlabel('bar speed')
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
                
            if printyn
                figSize = [30 21];
                set(gcf,'paperunits','centimeters','papersize',figSize,'paperposition',[0 0 figSize])
                if cll<10, printi = ['0', int2str(cll)]; else printi = int2str(cll); end
                print([printpath ,  printi '__' files(fi).name '.pdf'],'-dpdf')
            end
            
            % running figure
            runlfr(cll) = mean(frs(intersect(find(result.light),oktrials)));
            runnlfr(cll) = mean(frs(intersect(find(~result.light),oktrials)));
            norunlfr(cll) = mean(frs(intersect(find(result.light),stilltrials)));
            norunnlfr(cll) = mean(frs(intersect(find(~result.light),stilltrials)));
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
            lp(cll) = p(1); rp(cll) = p(2); rlip(cll) = p(3);
            
            clf;
            subplot(2,2,1)
            imagesc(runspeed);
            colorbar
            title(['oktrials: ' int2str(length(oktrials)) '/' int2str(size(runspeed,1))])
            xlabel('time [ms]')
            ylabel('trial number')
            
            subplot(2,2,2)
            errorbar(mean(runspeed(find(~light),:)),std(runspeed(find(~light),:))./sqrt(length(find(~light))),'b')
            hold on
            errorbar(mean(runspeed(find(light),:)),std(runspeed(find(light),:))./sqrt(length(find(light))),'r')
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
            
            if printyn
                figSize = [30 21];
                set(gcf,'paperunits','centimeters','papersize',figSize,'paperposition',[0 0 figSize])
                if cll<10, printi = ['0', int2str(cll)]; else printi = int2str(cll); end
                print([runprintpath ,  printi '__' files(fi).name '.pdf'],'-dpdf')
            end
            
            clf;% LFP figure
            fillx = [sfx,fliplr(sfx)];
            fillyrun0 = [squeeze(nanmean(condrunSerr(cll,1,:,1,:),3))',fliplr(squeeze(nanmean(condrunSerr(cll,1,:,2,:)))')];
            fillyrun1 = [squeeze(nanmean(condrunSerr(cll,2,:,1,:),3))',fliplr(squeeze(nanmean(condrunSerr(cll,2,:,2,:)))')];
            fillystill0 = [squeeze(nanmean(condstillSerr(cll,1,:,1,:),3))',fliplr(squeeze(nanmean(condstillSerr(cll,1,:,2,:)))')];
                        
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
            
            if printyn
                figSize = [30 21];
                set(gcf,'paperunits','centimeters','papersize',figSize,'paperposition',[0 0 figSize])
                if cll<10, printi = ['0', int2str(cll)]; else printi = int2str(cll); end
                print([lfpprintpath ,  printi '__' files(fi).name '.pdf'],'-dpdf')
            end
            
            disp([files(fi).name '   done'])
            cll = cll + 1;
            
        end
    end
    save(popfile, '-v7.3');
else
    load(popfile);   
end

toc

clear cell;
%spike classification
% kmeansind = kmeans([eslope',ptr',swidth',adiff'],4);
kmeansind = kmeans([swidth',adiff'],3);
secpersamp = 1/30000;
interpf = secpersamp/10;
swidthms = swidth*interpf*1000;

m1 = mean(swidth(kmeansind==1)); m2 = mean(swidth(kmeansind==2)); m3 = mean(swidth(kmeansind==3));
fm = min([m1,m2,m3]); fmax = max([m1,m2,m3]);
if fm == m1
    pfs = find(kmeansind==1); pfsv = kmeansind==1;
elseif fm == m2
    pfs = find(kmeansind==2); pfsv = kmeansind == 2;
else
    pfs = find(kmeansind==3); pfsv = kmeansind == 3;
end
if fmax == m1
    prs = find(kmeansind==1); prsv = kmeansind==1;
elseif fmax == m2
    prs = find(kmeansind==2); prsv = kmeansind == 2;
else
    prs = find(kmeansind==3); prsv = kmeansind == 3;
end

% m1 = mean(swidth(kmeansind==1)); m2 = mean(swidth(kmeansind==2)); m3 = mean(swidth(kmeansind==3));
% if find(kmeansind == 4), m4 = mean(swidth(kmeansind == 4)); end
% fm = min([m1,m2,m3,m4]);
% if fm == m1
%     pfs = find(kmeansind==1); prs = find(kmeansind==2|kmeansind==3|kmeansind==4); 
%     pfsv = kmeansind==1; prsv = kmeansind==2|kmeansind==3|kmeansind==4;
% elseif fm == m2
%     pfs = find(kmeansind==2); prs = find(kmeansind==1|kmeansind==3|kmeansind==4);
%     pfsv = kmeansind==2; prsv = kmeansind==1|kmeansind==3|kmeansind==4;
% elseif fm == m3
%     pfs = find(kmeansind==3); prs = find(kmeansind==1|kmeansind==2|kmeansind==4);
%     pfsv = kmeansind==3; prsv = kmeansind==1|kmeansind==2|kmeansind==4;
% else
%     pfs = find(kmeansind==4); prs = find(kmeansind==1|kmeansind==2|kmeansind==3);
%     pfsv = kmeansind==4; prsv = kmeansind==1|kmeansind==2|kmeansind==3;
% end

% if mean(swidth(find(kmeansind==1)))<mean(swidth(find(kmeansind==2)))  %1 is FS
%     pfs = find(kmeansind==1); prs = find(kmeansind==2); pfsv = kmeansind==1;
% else
%     pfs = find(kmeansind==2); prs = find(kmeansind==1); pfsv = kmeansind==2;
% end

figure
plot(swidthms(kmeansind==1),adiff(kmeansind==1),'b.')
xlabel('spike width')
ylabel('amplitude diff')
hold on
plot(swidthms(kmeansind==2),adiff(kmeansind==2),'r.')
if ~isempty(find(kmeansind==3))
    plot(swidthms(kmeansind==3),adiff(kmeansind==3),'g.')
    plot(swidthms(pfsv),adiff(pfsv),'ro');
    plot(swidthms(prsv),adiff(prsv),'o');
end

% figure
% plot(eslope(kmeansind==1),adiff(kmeansind==1),'b.')
% xlabel('end slope')
% ylabel('amplitude diff')
% hold on
% plot(eslope(kmeansind==2),adiff(kmeansind==2),'r.')
% if ~isempty(find(kmeansind==3))
%     plot(eslope(kmeansind==3),adiff(kmeansind==3),'g.')
%     plot(eslope(pfsv),adiff(pfsv),'ro');
%     plot(eslope(prsv),adiff(prsv),'o');
% end
%     % axis([5.5,20.5,-.9,.7])

prsv = swidthms>=.38; pfsv = swidthms<=.36;
prs = find(prsv); pfs = find(pfsv);

swamp = max(waveform,[],2)-min(waveform,[],2);
okwv = swamp>42;

for i = 1:length(okwv)
    isodist(i) = clustqual(i).IsolationDistance;
    lratio(i) = clustqual(i).L_Ratio.Lratio;
end

% ok = vismod&okwv'&nlfr>1;
% ok = okwv'&isodist>8;
ok = okwv';

% adjust depth according to penetration angle
% pangle = 10; %TODO temp fix
% depth = depth.*cosd(22).*cosd(pangle);
depth = depth.*cosd(pangle);

phe = zeros(1,length(depth));

% %putative halo expressing for X94 ChR2
% phe([]) = 1; 

phe = logical(phe);

[sd,si] = sort(depth);
[rsd,rsi] = sort(depth(prs));
[fsd,fsi] = sort(depth(pfs));

l23 = depth<375;
l4 = depth>=375&depth<=550;
l5 = depth>550&depth<=800;
l5a = depth>550&depth<=650;
l5b = depth>650&depth<=800;
l6 = depth>800;
l23rs = l23&prsv&~phe&ok;
l23fs = l23&pfsv&~phe&ok;
l4rs = l4&prsv&~phe&ok;
l4fs = l4&pfsv&~phe&ok;
l5rs = l5&prsv&~phe&ok;
l5fs = l5&pfsv&~phe&ok;
l6rs = l6&prsv&~phe&ok;
l6fs = l6&pfsv&~phe&ok;
okrs = prsv&~phe&ok; 

for i = 1:length(depth)
    omi(i) = (nanmean(runcondfr(i,2,:),3)-nanmean(runcondfr(i,1,:),3))./(nanmean(runcondfr(i,2,:),3)+nanmean(runcondfr(i,1,:),3));
    ominr(i) = (nanmean(stillcondfr(i,2,:),3)-nanmean(stillcondfr(i,1,:),3))./(nanmean(stillcondfr(i,2,:),3)+nanmean(stillcondfr(i,1,:),3));
    for l = 1:2
        ssi(i,l) = get_ssi(squeeze(runcondfr(i,l,2:6)));
    end
end

% running OMI in depth
figure
plot(omi(prs),depth(prs),'ko','markerfacecolor','k')
hold on
plot(omi(pfs),depth(pfs),'go','markerfacecolor','g')
axis ij
line([0,0],[100,900],'color','k')
line([-1,1],[550,550],'color','k','linestyle',':')
line([-1,1],[375,375],'color','k','linestyle',':')
legend('RS','FS')
xlabel('OMI')
ylabel('cortical depth')
title('OMI running')

% non-running OMI in depth
figure
plot(ominr(prs),depth(prs),'ko','markerfacecolor','k')
hold on
plot(ominr(pfs),depth(pfs),'go','markerfacecolor','g')
axis ij
line([0,0],[100,900],'color','k')
line([-1,1],[550,550],'color','k','linestyle',':')
line([-1,1],[375,375],'color','k','linestyle',':')
legend('RS','FS')
xlabel('OMI')
ylabel('cortical depth')
title('OMI non-running')

% spatial selectivity
figure
plot(ssi(prs,1),depth(prs),'ko','markerfacecolor','k')
hold on
plot(ssi(pfs,1),depth(pfs),'go','markerfacecolor','g')
axis ij
line([0,1],[550,550],'color','k','linestyle',':')
line([0,1],[375,375],'color','k','linestyle',':')
legend('RS','FS')
xlabel('SSI')
ylabel('cortical depth')
title('SSI running')

% spatial selectivity change
figure
plot(ssi(prs,2)-ssi(prs,1),depth(prs),'ko','markerfacecolor','k')
hold on
plot(ssi(pfs,2)-ssi(pfs,1),depth(pfs),'go','markerfacecolor','g')
axis ij
line([0,0],[100,900],'color','k')
line([-.6,.6],[550,550],'color','k','linestyle',':')
line([-.6,.6],[375,375],'color','k','linestyle',':')
axis([-.6,.6,100,900])
legend('RS','FS')
xlabel('SSI')
ylabel('cortical depth')
title('change in SSI with light')



function ssi = get_ssi(curve)
ssi  = 1 - (((norm(curve)/max(curve)) - 1)./((sqrt(length(curve)))-1));