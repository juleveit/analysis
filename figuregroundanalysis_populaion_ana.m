function figureground_populationana

% % any animal
% animalids = {'170712','170804','170807','170808','170810','170814','170817','171201','180125','180330','180404','180410','180412','180413'};
% blocks    = [3,        5,       6,       3,       5,       7,       6,       6,       5,       12,      6,       8,       8,       7];
% rfblocks  = [2,        4,       5,       2,       4,       6,       5,       7,       8,       9,       8,       7,       7,       6];
% animal    = [1,        2,       3,       4,       5,       6,       7,       8,       9,       10,      11,      12,      13,      14];
% electrodes =[[1,16];   [1,16];  [1,16];  [1,16];  [1,16];  [1,16];  [1,16]; [1,16];  [1,16];  [1,16];  [1,16];  [1,16];  [1,16];  [1,16]];
% penangle =  [30,       30,      30,      30,      30,      30,      30,      30,      30,      30,      30,      30,      30,      30];
% printpath = 'C:\Users\Julia\work\data\populations\control\figureground\units\';
% lfpprintpath = 'C:\Users\Julia\work\data\populations\control\figureground\lfp\';
% rfprintpath = 'C:\Users\Julia\work\data\populations\control\figureground\rfs\';
% popfile = 'C:\Users\Julia\work\data\populations\control\figureground\figureground_population.mat';

% VIP Halo
animalids = {'171201','180125','180404','180410','180412','180413'};
blocks    = [ 6,       5,       6,       8,       8,       7];
rfblocks  = [ 7,       8,       8,       7,       7,       6];
animal    = [ 8,       9,       11,      12,      13,      14];
electrodes =[[1,16];  [1,16];  [1,16];  [1,16];  [1,16];  [1,16]];
penangle =  [ 30,      30,      30,      30,      30,      30];
printpath = 'C:\Users\Julia\work\data\populations\VIP_Halo\figureground\units\';
lfpprintpath = 'C:\Users\Julia\work\data\populations\VIP_Halo\figureground\lfp\';
rfprintpath = 'C:\Users\Julia\work\data\populations\VIP_Halo\figureground\rfs\';
popfile = 'C:\Users\Julia\work\data\populations\VIP_Halo\figureground\figureground_population.mat';

% % SOM Halo
% animalids = {'180330'};
% blocks    = [ 12];
% rfblocks  = [ 9];
% animal    = [ 10];
% electrodes =[[1,16]];
% penangle =  [ 30];
% printpath = 'C:\Users\Julia\work\data\populations\SOM_Halo_later\figureground\units\';
% lfpprintpath = 'C:\Users\Julia\work\data\populations\SOM_Halo_later\figureground\lfp\';
% rfprintpath = 'C:\Users\Julia\work\data\populations\SOM_Halo_later\figureground\rfs\';
% popfile = 'C:\Users\Julia\work\data\populations\SOM_Halo_later\figureground\figureground_population.mat';

tic
lcol = 'r'; %lasercolor

recalculate = 0;
printyn = 1;
sfc = 0;

addpath C:\Users\Julia\work\Matlab\others\fieldtrip-20160329
ft_defaults

freqbinwidth = 5;

% chronux parameters
params.tapers = [2,5]; params.Fs = 1000; params.err = [2, 0.05]; params.trialave = 1;

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

sr = 1000;

if ~exist(popfile) || recalculate

    cll = 1;
    for blck = 1:length(blocks)
        
        supath = ['C:\Users\Julia\work\data\' animalids{blck} '\singleunits\'];
        basename = [animalids{blck} '_block' int2str(blocks(blck)) '_tet'];
        snbasename = [animalids{blck} '_block' int2str(rfblocks(blck)) '_tet'];

        files = dir([supath, basename, '*.mat']);
        rffiles = dir([supath, snbasename, '*.mat']);
        
        clear filtmat; clear powmat; clear phasmat;

        for fi = 1:length(files)

            if strfind(files(fi).name, 'MU')
                continue;
            end
            
            load([supath, files(fi).name]);   

            prestim = 300;
            poststim = 700;
            if result.stimduration == 2
                respwin = 501:1500; % after stimulus onset
            else
                respwin = 1:1000;
            end
            respwin = respwin+prestim;
            
            cllname{cll} = files(fi).name;
            printname = files(fi).name;
            printname(find(printname=='_')) = ' ';
            animalno(cll) = animal(blck);
            recording(cll) = blck;
            depth(cll) = result.depth;
            pangle(cll) = penangle(blck);
            
            i = strfind(files(fi).name, 'tet');
            tetno(cll) = strread(files(fi).name(i+3));
                        
            wvchan = find(var(result.waveforms) == max(var(result.waveforms)));
            spike = result.waveforms(:,wvchan);
            interpspike = spline(1:32,spike,1:.1:32);
            [adiff(cll),swidth(cll),ptr(cll),eslope(cll)] = spikequant(interpspike);
            
            waveform(cll,:) = spike;
            clustqual(cll) = result.clusterquality;
            
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
                %          msstamps([62,108,147]) = []; % for 140703 block 8
                %          newmsstamps = msstamps(1:319); % for 170712 block 3
                %          newmsstamps(320) = newmsstamps(319)+3683;
                %          newmsstamps(321) = newmsstamps(320)+3683;
                %          newmsstamps(322) = newmsstamps(321)+3683;
                %          newmsstamps = [newmsstamps; msstamps(320:end)];
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
                filtlfpresp(i,:) = filtlfp(msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim);
                [lfpspect(i,:),trialfax] = pmtm(lfpresp(i,respwin(1)+200:respwin(end)),3,[],sr);
                [lfpsectrogram(i,:,:),ct,cf] = mtspecgramc(lfpresp(i,:)',[.25,.1],params);
                gammaresp(i,:) = gamma(msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim);
                gammapowresp(i,:) = gpow(msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim);
                gphaseresp(i,:) = gphas(msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim);
                bpresp(i,:) = gamma(msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim);
                
                speed(i,:) = result.runspeed(msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim);
                
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
            
            cllresp(cll,:,:) = resp;
            clllfpresp(cll,:,:) = lfpresp;
            
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
                            condresperr(cll,lg,sz,oc,os,:) = nanstd(resp(trials,:),1,1)./sqrt(length(trials));
                            condlfpspect(cll,lg,sz,oc,os,:) = nanmean(lfpspect(trials,:));
                            condstspect(cll,lg,sz,oc,os,:) = pmtm(mean(resp(trials,1001:1800),1),3,[],sr);
                            condfr(cll,lg,sz,oc,os) = mean(frs(trials));%-mean(bl);
                            conderr(cll,lg,sz,oc,os) =std(frs(trials))./sqrt(length(trials));
                            
                            condfiltresp(cll,lg,sz,oc,os,:) = filter(excit_kernel,1,mean(resp(trials,:),1));
                            
                            condz(cll,lg,sz,oc,os) = {(sc(trials)-mean(sc(trials)))/std(sc(trials))}; %ecker 2010
                            condsc(cll,lg,sz,oc,os) = {sc(trials)};
                            ff(cll,lg,sz,oc,os) = var(sc(trials))/mean(sc(trials));
                            
                            condlfpresp(cll,lg,sz,oc,os,:) = mean(lfpresp(trials,:),1);
                            
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
                            
%                             
%                             % field trip messing around
%                             if ~isempty(trials)
%                                 for i = 1:length(trials)
%                                     data.trial{i} = [lfpresp(trials(i),:);resp(trials(i),:)];
%                                     time = -prestim+1:result.stimduration*1000+poststim;
%                                     data.time{i} = time(1)./1000:.001:time(end)./1000;
%                                     %                 data.time{i} = -.299:.001:2.700;
%                                 end
%                                 data.fsample = 1000;
%                                 data.trialinfo(:,1) = light(trials);
%                                 data.trialinfo(:,2) = gratingInfo.sizecenter(trials);
%                                 data.trialinfo(:,3) = gratingInfo.Orientation_center(trials);
%                                 cfg = [];
%                                 cfg.timwin = [-.25 .25];
%                                 data.label{1,1} = 'lfp';
%                                 data.label{2,1} = 'spikes';
%                                 cfg.spikechannel = 'spikes';
%                                 cfg.channel = 'lfp';
%                                 cfg.latency = [0.5,1.5];
%                                 
%                                 cfg.method = 'mtmfft';
%                                 cfg.foilim = [5,100];
%                                 cfg.timwin = [-.15, .15];
%                                 cfg.taper = 'hanning';
%                                 cfg.spikechannel = 'spikes';
%                                 cfg.channel = 'lfp';
%                                 stsFFT           = ft_spiketriggeredspectrum(cfg, data);
%                                 ang = angle(stsFFT.fourierspctrm{1});
%                                 mag = abs(stsFFT.fourierspctrm{1});
%                                 ftspect(cll,lg,sz,:) = squeeze(nanmean(mag(:,1,:)));
%                                 ftphases{cll,lg,sz} = squeeze(ang);
%                                 ftfax = stsFFT.freq;
%                                 
%                                 cfg               = [];
%                                 cfg.method        = 'ral'; % compute the rayleigh test
%                                 cfg.spikechannel  = stsFFT.label{1};
%                                 cfg.channel       = stsFFT.lfplabel; % selected LFP channels
%                                 cfg.avgoverchan   = 'unweighted'; % weight spike-LFP phases irrespective of LFP power
%                                 cfg.timwin        = 'all'; % compute over all available spikes in the window
%                                 cfg.latency       = [0.5 1.5]; % sustained visual stimulation period
%                                 statSts           = ft_spiketriggeredspectrum_stat(cfg,stsFFT);
%                                 ftralp(cll,lg,sz,oc,os,:) = statSts.ral;
%                                 cfg.method = 'ppc0';
%                                 statSts           = ft_spiketriggeredspectrum_stat(cfg,stsFFT);
%                                 ftppc(cll,lg,sz,oc,os,:) = statSts.ppc0;
%                                 cfg.method = 'plv';
%                                 statSts           = ft_spiketriggeredspectrum_stat(cfg,stsFFT);
%                                 ftplv(cll,lg,sz,oc,os,:) = statSts.plv;
%                             else
%                                 ftralp(cll,lg,sz,oc,os,:) = nan(1,30);
%                                 ftppc(cll,lg,sz,oc,os,:) = nan(1,30);
%                                 ftplv(cll,lg,sz,oc,os,:) = nan(1,30);
%                             end
                        end
                    end
                end
            end
            sfx = chf(1:150);
            
            bta = bta-prestim;
            
            printname = files(fi).name;
            printname(find(printname=='_')) = ' ';
            
            
            % BL, sm md lg bgonly smbg mdbg lgbg smap mdap lgap
            for lg = 1:size(condfr,2)
                vecfr(cll,lg,1) = condfr(cll,lg,1,1,1); vecerr(cll,lg,1) = conderr(cll,lg,1,1,1);
                vecfr(cll,lg,2) = nanmean(condfr(cll,lg,1,2:3,1),4); vecerr(cll,lg,2) = nanmean(conderr(cll,lg,1,2:3,1),4);
                vecfr(cll,lg,3) = nanmean(condfr(cll,lg,2,2:3,1),4); vecerr(cll,lg,3) = nanmean(conderr(cll,lg,2,2:3,1),4);
                vecfr(cll,lg,4) = nanmean(condfr(cll,lg,3,2:3,1),4); vecerr(cll,lg,4) = nanmean(conderr(cll,lg,3,2:3,1),4);
                vecfr(cll,lg,5) = nanmean([condfr(cll,lg,1,2,2),condfr(cll,lg,1,3,3)]); vecerr(cll,lg,5) = nanmean([conderr(cll,lg,1,2,2),conderr(cll,lg,1,3,3)]);
                vecfr(cll,lg,6) = nanmean([condfr(cll,lg,1,2,3),condfr(cll,lg,1,3,2)]); vecerr(cll,lg,6) = nanmean([conderr(cll,lg,1,2,3),conderr(cll,lg,1,3,2)]);
                vecfr(cll,lg,7) = nanmean([condfr(cll,lg,2,2,3),condfr(cll,lg,2,3,2)]); vecerr(cll,lg,7) = nanmean([conderr(cll,lg,2,2,3),conderr(cll,lg,2,3,2)]);
                vecfr(cll,lg,8) = nanmean([condfr(cll,lg,3,2,3),condfr(cll,lg,3,3,2)]); vecerr(cll,lg,8) = nanmean([conderr(cll,lg,3,2,3),conderr(cll,lg,3,3,2)]);
                vecfr(cll,lg,9) = nanmean(condfr(cll,lg,1,1,2:3),5); vecerr(cll,lg,9) = nanmean(conderr(cll,lg,1,1,2:3),5);
                vecfr(cll,lg,10) = nanmean(condfr(cll,lg,2,1,2:3),5); vecerr(cll,lg,10) = nanmean(conderr(cll,lg,2,1,2:3),5);
                vecfr(cll,lg,11) = nanmean(condfr(cll,lg,3,1,2:3),5); vecerr(cll,lg,11) = nanmean(conderr(cll,lg,3,1,2:3),5);
            end
            
            clf;
            subplot(2,2,1)
            errorbar(vecfr(cll,1,:),vecerr(cll,1,:),'k.')
            set(gca,'xtick',1:11)
            set(gca,'xticklabel',{'BL','sm','md','lg','bgonly','smbg','mdbg','lgbg','smap','mdap','lgap'})  
            
            subplot(2,2,2)
            errorbar(1:3,vecfr(cll,1,2:4),vecerr(cll,1,2:4),'k.-','linewidth',2);
            hold on
            errorbar(1:3,vecfr(cll,1,6:8),vecerr(cll,1,6:8),'r.-','linewidth',2);
            errorbar(1:3,vecfr(cll,1,9:11),vecerr(cll,1,9:11),'color',[.75,.75,.75]);
            errorbar(0,vecfr(cll,1,1),vecerr(cll,1,1),'k.','linewidth',2);
            errorbar(0,vecfr(cll,1,5),vecerr(cll,1,5),'r.','linewidth',2);
            ax = axis;
            axis([-.5,3.5,ax(3),ax(4)]);
            xlabel('sizes')
            set(gca,'xtick',0:1:3);
            set(gca,'xticklabel',{'BG','8','20','40'})
            legend('gray background','drifting background','aperture') 
            
            subplot(2,2,3)
            plot(bta, nanmean([squeeze(bincondresp(cll,1,1,2,2,:)),squeeze(bincondresp(cll,1,1,3,3,:))],2),'k','linewidth',2);
            hold on
            plot(bta, nanmean([squeeze(bincondresp(cll,1,1,2,3,:)),squeeze(bincondresp(cll,1,1,3,2,:))],2),'b','linewidth',2);
            plot(bta, nanmean([squeeze(bincondresp(cll,1,2,2,3,:)),squeeze(bincondresp(cll,1,2,3,2,:))],2),'c','linewidth',2);
            plot(bta, nanmean([squeeze(bincondresp(cll,1,3,2,3,:)),squeeze(bincondresp(cll,1,3,3,2,:))],2),'g','linewidth',2);
            legend('BG','SM','MD','LG')
            title('with background')
            
            subplot(2,2,4)
            plot(bta, nanmean([squeeze(bincondresp(cll,1,1,2,2,:)),squeeze(bincondresp(cll,1,1,3,3,:))],2),'k','linewidth',2);
            hold on
            plot(bta, nanmean([squeeze(bincondresp(cll,1,3,2,3,:)),squeeze(bincondresp(cll,1,3,3,2,:))],2),'g','linewidth',2);
            legend('BL','LG')
            
            if printyn
                figSize = [30 21];
                set(gcf,'paperunits','centimeters','papersize',figSize,'paperposition',[0 0 figSize])
                if cll<10, printi = ['0', int2str(cll)]; else printi = int2str(cll); end
                print([printpath ,  printi '__' files(fi).name '.pdf'],'-dpdf')
            end  
            
            clf;
            semilogy(sfx, squeeze(nanmean(condS(cll,1,3,2:3,1,:),4)),'b','linewidth',2);
            hold on
            semilogy(sfx,nanmean([squeeze(condS(cll,1,3,2,3,:)),squeeze(condS(cll,1,3,3,2,:))],2),'g','linewidth',2);
            semilogy(sfx,nanmean([squeeze(condS(cll,1,1,2,2,:)),squeeze(condS(cll,1,1,3,3,:))],2),'k','linewidth',2);
            legend('LG on gray','LG on BG','BG')
            if printyn
                figSize = [30 21];
                set(gcf,'paperunits','centimeters','papersize',figSize,'paperposition',[0 0 figSize])
                if cll<10, printi = ['0', int2str(cll)]; else printi = int2str(cll); end
                print([lfpprintpath ,  printi '__' files(fi).name '.pdf'],'-dpdf')
            end  
            
            clf;
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
            if printyn
                figSize = [30 21];
                set(gcf,'paperunits','centimeters','papersize',figSize,'paperposition',[0 0 figSize])
                if cll<10, printi = ['0', int2str(cll)]; else printi = int2str(cll); end
                print([rfprintpath ,  printi '__' files(fi).name '.pdf'],'-dpdf')
            end

            disp(['done with ' cllname{cll}]);
            
            cll = cll + 1;
        end
    end    
    bincondresp = bincondresp.*(1000/binwidth);
    binconderr = binconderr.*(1000/binwidth);
    
    save(popfile, '-v7.3'); 
else
    load(popfile);   
end


%spike classification
kmeansind = kmeans([swidth',adiff'],3);
% kmeansind = kmeans([eslope',adiff'],2);
secpersamp = 1/30000;
interpf = secpersamp/10;
swidthms = swidth*interpf*1000;

% m1 = mean(swidth(kmeansind==1)); m2 = mean(swidth(kmeansind==2)); m3 = mean(swidth(kmeansind==3));
% fm = min([m1,m2,m3]); fmax = max([m1,m2,m3]);
% if fm == m1
%     pfs = find(kmeansind==1); pfsv = kmeansind==1;
% elseif fm == m2
%     pfs = find(kmeansind==2); pfsv = kmeansind == 2;
% else
%     pfs = find(kmeansind==3); pfsv = kmeansind == 3;
% end
% if fmax == m1
%     prs = find(kmeansind==1); prsv = kmeansind==1;
% elseif fmax == m2
%     prs = find(kmeansind==2); prsv = kmeansind == 2;
% else
%     prs = find(kmeansind==3); prsv = kmeansind == 3;
% end

m1 = mean(swidth(kmeansind==1)); m2 = mean(swidth(kmeansind==2)); m3 = mean(swidth(kmeansind==3));
fm = min([m1,m2,m3]);
if fm == m1
    pfs = find(kmeansind==1); prs = find(kmeansind==2|kmeansind==3); 
    pfsv = kmeansind==1; prsv = kmeansind==2|kmeansind==3;
elseif fm == m2
    pfs = find(kmeansind==2); prs = find(kmeansind==1|kmeansind==3);
    pfsv = kmeansind==2; prsv = kmeansind==1|kmeansind==3;
else
    pfs = find(kmeansind==3); prs = find(kmeansind==1|kmeansind==2);
    pfsv = kmeansind==3; prsv = kmeansind==1|kmeansind==2;
end

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

prsv = swidthms>=.38; pfsv = swidthms<=.36;
prs = find(prsv); pfs = find(pfsv);

swamp = max(waveform,[],2)-min(waveform,[],2);
okwv = swamp>42;
vismod(isnan(vismod)) = 0;

for i = 1:length(okwv)
    isodist(i) = clustqual(i).IsolationDistance;
    lratio(i) = clustqual(i).L_Ratio.Lratio;
end

% ok = vismod&okwv'&nlfr>1;
% ok = okwv'&isodist>8;
ok = okwv';

% adjust depth according to penetration angle
depth = depth.*cosd(22).*cosd(pangle);

% halo expressing cells
phe = zeros(1,length(depth));

% % putative halo expressing for SOM later pop: 
% phe([6,35,64,75,172]) = 1;

% putative halo expressing for PV Halo pop: 
% phe([111,185,213,224]) = 1;

phe = logical(phe)';

[sd,si] = sort(depth);
[rsd,rsi] = sort(depth(prs));
[fsd,fsi] = sort(depth(pfs));

% fsl4 = find(depth(pfs)>=375 & depth(pfs)<=500);
% rsl4 = find(depth(prs)>=375 & depth(prs)<=500);
% fsl23 = find(depth(pfs)<375);
% rsl23 = find(depth(prs)<375);
% fsl5 = find(depth(pfs)>500 & depth(pfs)<=800);
% rsl5 = find(depth(prs)>500 & depth(prs)<=800);
l23 = depth<375;
l4 = depth>=375&depth<=550;
l5 = depth>550&depth<=800;
l5a = depth>550&depth<=650;
l5b = depth>650&depth<=800;
l6 = depth>800;
l23rs = l23&prsv&~phe'&ok;
l23fs = l23&pfsv&~phe'&ok;
l4rs = l4&prsv&~phe'&ok;
l4fs = l4&pfsv&~phe'&ok;
l5rs = l5&prsv&~phe'&ok;
l5fs = l5&pfsv&~phe'&ok;
l6rs = l6&prsv&~phe'&ok;
l6fs = l6&pfsv&~phe'&ok;
okrs = prsv&~phe'&ok;

targetdepth = 300;
anmls = unique(animalno);
for i = 1:length(anmls)
    ais = find(animalno == anmls(i));
    [xx,j] = min(abs(depth(ais)-300));
    lfpinds(i) = ais(j);
end


for i = 1:length(depth)
    for j = 1:11
        omi(i,j) = (vecfr(i,2,j)-vecfr(i,1,j))./(vecfr(i,2,j)+vecfr(i,1,j));
        logratio(i,j) = log10(vecfr(i,2,j)./vecfr(i,1,j));
        frdiff(i,j) = vecfr(i,2,j)-vecfr(i,1,j);
    end
    spontfr(i) = condfr(i,1,1,1,1);
    bgfr(i) = nanmean([condfr(i,1,1,2,2),condfr(i,1,1,3,3)]); 
    bgerr(i) = nanmean([conderr(i,1,1,2,2),conderr(i,1,1,3,3)]);
    lgfr(i) = nanmean([condfr(i,1,3,2,3),condfr(i,1,3,3,2)]); 
    lgerr(i) = nanmean([conderr(i,1,3,2,3),conderr(i,1,3,3,2)]);
    mdfr(i) = nanmean([condfr(i,1,2,2,3),condfr(i,1,2,3,2)]); 
    mderr(i) = nanmean([conderr(i,1,2,2,3),conderr(i,1,2,3,2)]);
    
    lgonlyspect(i,:) = nanmean(condlfpspect(i,1,3,2:3,1,:),4);
    mdonlyspect(i,:) = nanmean(condlfpspect(i,1,2,2:3,1,:),4);
    smonlyspect(i,:) = nanmean(condlfpspect(i,1,1,2:3,1,:),4);
    lgonbgspect(i,:) = nanmean([squeeze(condlfpspect(i,1,3,2,3,:)),squeeze(condlfpspect(i,1,3,3,2,:))],2);
    mdonbgspect(i,:) = nanmean([squeeze(condlfpspect(i,1,2,2,3,:)),squeeze(condlfpspect(i,1,2,3,2,:))],2);
    smonbgspect(i,:) = nanmean([squeeze(condlfpspect(i,1,1,2,3,:)),squeeze(condlfpspect(i,1,1,3,2,:))],2);
    bgonlyspect(i,:) = nanmean([squeeze(condlfpspect(i,1,1,2,2,:)),squeeze(condlfpspect(i,1,1,3,3,:))],2);
    
    spontbinresp(i,:) = squeeze(bincondresp(i,1,1,1,1,:));
    bgonlybinresp(i,:) = nanmean([squeeze(bincondresp(i,1,1,2,2,:)),squeeze(bincondresp(i,1,1,3,3,:))],2);
    lgfigbinresp(i,:) = nanmean([squeeze(bincondresp(i,1,3,2,3,:)),squeeze(bincondresp(i,1,3,3,2,:))],2);
    mdfigbinresp(i,:) = nanmean([squeeze(bincondresp(i,1,2,2,3,:)),squeeze(bincondresp(i,1,2,3,2,:))],2);
    bgonlybinerr(i,:) = nanmean([squeeze(binconderr(i,1,1,2,2,:)),squeeze(binconderr(i,1,1,3,3,:))],2);
    lgfigbinerr(i,:) = nanmean([squeeze(binconderr(i,1,3,2,3,:)),squeeze(binconderr(i,1,3,3,2,:))],2);
    mdfigbinerr(i,:) = nanmean([squeeze(binconderr(i,1,2,2,3,:)),squeeze(binconderr(i,1,2,3,2,:))],2);
    
    spontfiltresp(i,:) = squeeze(condfiltresp(i,1,1,1,1,:));
    bgonlyfiltresp(i,:) = nanmean([squeeze(condfiltresp(i,1,1,2,2,:)),squeeze(condfiltresp(i,1,1,3,3,:))],2);
    lgfigfiltresp(i,:) = nanmean([squeeze(condfiltresp(i,1,3,2,3,:)),squeeze(condfiltresp(i,1,3,3,2,:))],2);
    mdfigfiltresp(i,:) = nanmean([squeeze(condfiltresp(i,1,2,2,3,:)),squeeze(condfiltresp(i,1,2,3,2,:))],2);
    
    spontresp(i,:) = squeeze(condresp(i,1,1,1,1));
    bgonlyresp(i,:) = nanmean([squeeze(condresp(i,1,1,2,2,:)),squeeze(condresp(i,1,1,3,3,:))],2);
    lgfigresp(i,:) = nanmean([squeeze(condresp(i,1,3,2,3,:)),squeeze(condresp(i,1,3,3,2,:))],2);
    mdfigresp(i,:) = nanmean([squeeze(condresp(i,1,2,2,3,:)),squeeze(condresp(i,1,2,3,2,:))],2);
end

i = 170
figure
plot(bta,bgonlybinresp(i,:),'b','linewidth',2)
hold on
plot(bta,mdfigbinresp(i,:),'r','linewidth',2)
plot(bta,spontbinresp(i,:),'k','linewidth',2)