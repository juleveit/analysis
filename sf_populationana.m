function contrast_populationana


% SOM population 
animalids = {'151023', '151023', '151204', '151204', '151209'};
blocks    = [16,        16,       9,        9,        10];
animal =    [1          1         2         2         3];
electrodes =[[1,16];   [17,32];  [1,16];   [17,32];  [17,32]];
penangle =  [25,        25,       25,       25,       25];
printpath =     'C:\Users\Julia\work\data\populations\SOM_Halo\SF\units\';
runprintpath =  'C:\Users\Julia\work\data\populations\SOM_Halo\SF\running\';
condprintpath = 'C:\Users\Julia\work\data\populations\SOM_Halo\SF\conditions\';
popfile =       'C:\Users\Julia\work\data\populations\SOM_Halo\SF\SF_population.mat';


% % SOM and control - look at no light population 
% animalids = {'151023', '151023', '151204', '151204', '151209', '160125', '160125', '160128', '160129', '160129'};
% blocks    = [16,        16,       9,        9,        10,        8,        8,        11,      6,        6];
% animal =    [1          1         2         2         3          5         5,        6,       7,        7];
% electrodes =[[1,16];   [17,32];  [1,16];   [17,32];  [17,32];   [1,16];   [17,32];  [1,16];  [1,16];   [17,32]];
% penangle =  [25,        25,       25,       25,       25,        25,       25,       25,      25,       25];
% printpath =     'C:\Users\Julia\work\data\populations\control\SF\units\';
% runprintpath =  'C:\Users\Julia\work\data\populations\control\SF\running\';
% condprintpath = 'C:\Users\Julia\work\data\populations\control\SF\conditions\';
% popfile =       'C:\Users\Julia\work\data\populations\control\SF\SF_population.mat';
% % popfile = 'C:\Users\Julia\work\data\populations\SOM_Halo\contrast\contrast_population_onlymale.mat';


tic
lcol = 'r'; %lasercolor

recalculate = 0;
printyn = 1;

screendeg = 62;
screencm = 28;
degpercm = screendeg/screencm;

% chronux parameters
params.tapers = [2,5]; params.Fs = 1000; params.err = [2, 0.05]; params.trialave = 1;

onlymod = 0; % run only for visually and laser modulated
sr = 1000;

if ~exist(popfile) || recalculate
    
    cll = 1;
    for blck = 1: length(blocks)

        supath = ['C:\Users\Julia\work\data\' animalids{blck} '\singleunits\'];
        basename = [animalids{blck} '_block' int2str(blocks(blck)) '_tet'];

        files = dir([supath, basename, '*.mat']);

        prestim = 300;
        poststim = 700;
        respwin = 501:1500; % after stimulus onset
        respwin = respwin+prestim;
        freqbinwidth = 5;
        
        clear filtmat; clear powmat; clear phasmat; 

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
            
            % get spiketimes
            msStimes = round(result.spikes);
            if ~isempty(msStimes) & msStimes(1) == 0, msStimes(1) = 1; end

            chan = zeros(1,length(result.lfp));
            chan(msStimes) = 1;
            
            wvchan = find(var(result.waveforms) == max(var(result.waveforms)));
            cm = [3,4,1,2]; % confusion matrix? take lfp from two electrodes away to not get too many spike related phase resets
            lfp = result.lfp(:,cm(wvchan))';      
            
            trialdur = result.stimduration*1000;
            msstamps = result.msstamps;
            
            if length(msstamps)~=length(result.light)
%             disp('');
%             msstamps(16) = []; % for 150414 block 10
%             result.msstamps = msstamps;
%             save([supath, SUfiles(cl).name],'result');            
                pause;
            end
            
            for i = 1:length(msstamps)
                speed(i,:) = result.runspeed(msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim);
            end
            
            meanspeed = mean(speed(:,respwin),2);
            stdspeed = std(speed(:,respwin),1,2);
            notstill = find(meanspeed>1);
            okspeed = find(meanspeed>( mean(meanspeed(notstill))-(1.5*std(meanspeed(notstill))) ) & meanspeed>1);
            okvar = find(stdspeed<( mean(stdspeed(notstill))+(1.5*std(stdspeed(notstill)))) & stdspeed>.5);
            oktrials = intersect(okspeed,okvar);
            nonoktrials = 1:size(speed,1); nonoktrials(oktrials) = [];
            stilltrials = 1:size(speed,1); stilltrials(notstill) = [];
            if ~isempty(oktrials)
                cllmeanspeed(cll) = mean(meanspeed(oktrials));
                cllminspeed(cll) = min(meanspeed(oktrials));
            else
                cllmeanspeed(cll) = NaN;
                cllminspeed(cll) = NaN;
            end
            
            %find gamma peaks for this animal
            beta = [20,40];
            gamma = [50,70];
            high = find(result.gratingInfo.spFreq == 0.04 & result.light == 0);
            lr1 = intersect(high,oktrials);
            for i = 1:length(high)
                [pl(i,:),f] = pmtm(lfp(result.msstamps(high(i)):result.msstamps(high(i))+1000),3,[],1000);
            end      
            plr1 = nan(1,size(pl,2));
            if length(lr1>=5)
                for i = 1:length(lr1)
                    [plr1(i,:),f] = mtspectrumc(lfp(result.msstamps(lr1(i))+700:result.msstamps(lr1(i))+1500),params);
                end
            end
            b1 = find(f>beta(1),1); b2 = find(f>beta(2),1);
            
            if ~isnan(plr1(1))
                bsig = nanmean(plr1(:,b1:b2));
            else
                bsig = nanmean(pl(:,b1:b2));
            end
            if isempty(find(diff(bsig)>0)) % there is no clear beta peak
                bpi = round((b1+b2)/2);
            else
                peaks = find(diff(bsig)>0)+1;
                pvs = bsig(peaks);
                bpi = peaks(pvs == max(pvs));
                bpi = bpi+b1-1;
            end
            
            gamma1 = eegfilt(lfp,sr,f(bpi)-2.5,f(bpi)+2.5);
            h1 = hilbert(gamma1); gpow1 = abs(h1); gphas1 = angle(h1);
            
            for i = 1:100/freqbinwidth
                filtmat(i,:) = eegfilt(lfp,sr,(i-1)*freqbinwidth+1,i*freqbinwidth);
                h = hilbert(filtmat(i,:));
                powmat(i,:) = abs(h); phasmat(i,:) = angle(h);
            end           
            
            clear lfpspect; clear fax;            
            for i = 1:length(msstamps)
                resp(i,:) = chan(msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim);
                lfpresp(i,:) = result.lfp(msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim,wvchan);
                gamma1resp(i,:) = gpow1(msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim);                
                gamma1phasresp(i,:) = gphas1(msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim); 
                speed(i,:) = result.runspeed(msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim);
                [lfpspect(i,:),fax] = pmtm(lfpresp(i,respwin),2,[],1000);
                allS(i,:)=mtspectrumc(squeeze(lfpresp(i,1001:1800))',params);                
                for j = 1:size(phasmat,1)
                    allphaseresp(j,i,:) = phasmat(j, msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim);
                    allpowresp(j,i,:) = powmat(j, msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim);
                end                
                isiresp{i} = diff(find(resp(i,respwin)));
                msl(i) = mean(find(resp(i,respwin))); %mean spike latency
            end
            
            lfpspect = lfpspect(:,1:150); fax = fax(1:150);
            
            msta = linspace(-prestim,trialdur+poststim,size(resp,2));

            lightresp = resp(find(result.light),:);
            nolightresp = resp(find(result.light == 0),:);
            lightlfpresp = lfpresp(find(result.light),:);
            nolightlfpresp = lfpresp(find(~result.light),:);

            frs = sum(resp(:,respwin),2)./(length(respwin)/1000);
            bl = sum(resp(:,1:prestim),2)./(prestim/1000);
            sc = sum(resp(:,respwin),2);
            gpower1 = mean(gamma1resp(:,respwin),2);
            
%             %determine if cll is visually modulated
%             blfr = sum(resp(:,1:prestim),2);
%             vrfr = sum(resp(:,prestim+40:2*prestim+40-1),2);
%             vm(fi) = ttest2(blfr,vrfr);
%             if isnan(vm(fi)), vm(fi) = 0; end
            
            tmpi = zeros(size(allphaseresp));
            for i = 1:size(allphaseresp,1)
                tmp = zeros(size(resp));
                tmp(find(resp)) = allphaseresp(i,find(resp));
                tmpi(i,:,:) = tmp;
            end
            g1 = zeros(size(gamma1phasresp));
            g1(find(resp)) = gamma1phasresp(find(resp));
            
            
            binwidth = 20;
            spFreqs = unique(result.gratingInfo.spFreq);   %delete control condition
            tFreqs = unique(result.gratingInfo.tFreq);
            spFreqscm = 1./(spFreqs.*degpercm); % in cm/cycle (* .5 = stripe thickness)
            for i = 1:length(spFreqs)
                tFreqscms(:,i) = spFreqscm(i).*tFreqs;
            end
            for l = 1:2
                for tf = 1:length(tFreqs)
                    for sf = 1:length(spFreqs)
                        thisinds = find(result.gratingInfo.tFreq == tFreqs(tf) &...
                            result.gratingInfo.spFreq == spFreqs(sf) & ...
                            result.light == l-1);
                        
                        condresp(l,tf,sf,:) = mean(resp(thisinds,:),1);
                        cgamma1resp(l,tf,sf,:) = mean(gamma1resp(thisinds,:),1);
                        condlfpresp(l,tf,sf,:) = mean(lfpresp(thisinds,:));
                        clfpspect(l,tf,sf,:) = nanmean(lfpspect(thisinds,:));
                        clfpspecterr(l,tf,sf,:) = nanstd(lfpspect(thisinds,:))./sqrt(length(thisinds));
                        cgpow1(l,tf,sf) = mean(gpower1(thisinds));
                        cfr(l,tf,sf) = mean(frs(thisinds));
                        allfr(l,tf,sf,:) = frs(thisinds);
                        cerr(l,tf,sf) =std(frs(thisinds))./sqrt(length(thisinds));
                        [bincondresp(l,tf,sf,:),bta] = binit(squeeze(condresp(l,tf,sf,:)),binwidth);
                        cmsl(l,tf,sf) = mean(msl(thisinds)); % mean spike latency
                        
                        for b = 1:size(phasmat,1)
                            condpowresp(l,tf,sf,b,:) = squeeze(mean(allpowresp(b,thisinds,:),2));
                            condpow(l,tf,sf,b) = squeeze(mean(mean(allpowresp(b,thisinds,respwin),2),3));
                        end
                        
                        %phase locking to all bands and the two peak bands
                        phasemat = tmpi(:,thisinds,:);
                        for b = 1:size(allphaseresp,1)
                            allphases = squeeze(phasemat(b,find(squeeze(phasemat(b,:,:)))));
                            condr(cll,b,l,tf,sf) = circ_r(allphases');
                            condppc(cll,b,l,tf,sf) = ppc(allphases);
                            condcmean(cll,b,l,tf,sf) = circ_mean(allphases');
                        end
                        g1phasmat = g1(thisinds,:);
                        allg1phases = squeeze(g1phasmat(find(g1phasmat)));
                        if isrow(allg1phases) allg1phases = allg1phases'; end
                        condg1r(cll,l,tf,sf) = circ_r(allg1phases);
                        condg1ppc(cll,l,tf,sf) = ppc(allg1phases);
                        condg1cmean(cll,l,tf,sf) = circ_mean(allg1phases);
                        
                        %running
                        thisruninds = intersect(thisinds,oktrials);
                        if ~isempty(thisruninds)
                            runcondresp(l,tf,sf,:) = mean(resp(thisruninds,:),1);
                            runcondfr(l,tf,sf) = mean(frs(thisruninds));
                            runconderr(l,tf,sf) = std(frs(thisruninds))./sqrt(length(thisruninds));
                            runclfpspect(l,tf,sf,:) = nanmean(lfpspect(thisruninds,:),1);
                            for b = 1:size(phasmat,1)
                                runcondpow(l,tf,sf,b) = squeeze(mean(mean(allpowresp(b,thisruninds,respwin),2),3));
                            end
                            runcondg1pow(l,tf,sf) = mean(gpower1(thisruninds));
                            g1runphasmat = g1(thisruninds,:);
                            allrung1phases = squeeze(g1runphasmat(find(g1runphasmat)));
                            if isrow(allrung1phases) allrung1phases = allrung1phases'; end
                            runcondg1ppc(l,tf,sf) = ppc(allrung1phases);
                            runcondg1cmean(l,tf,sf) = circ_mean(allrung1phases);
                        else
                            runcondresp(l,tf,sf,:) = nan(1,size(resp,2));
                            runcondfr(l,tf,sf) = NaN;
                            runconderr(l,tf,sf) = NaN;
                            runclfpspect(l,tf,sf,:) = nan(1,size(lfpspect,2));
                            for b = 1:size(phasmat,1)
                                runcondpow(l,tf,sf,b) = NaN;
                            end
                            runcondg1pow(l,tf,sf) = NaN;
                            runcondg1ppc(l,tf,sf) = NaN;
                            runcondg1cmean(l,tf,sf) = NaN;
                        end
                        
                        thisstillinds = intersect(thisinds,stilltrials);
                        if ~isempty(thisstillinds)
                            stillcondresp(l,tf,sf,:) = mean(resp(thisstillinds,:),1);
                            stillcondfr(l,tf,sf) = mean(frs(thisstillinds));
                            stillconderr(l,tf,sf) = std(frs(thisstillinds))./sqrt(length(thisstillinds));
                            stillclfpspect(l,tf,sf,:) = nanmean(lfpspect(thisstillinds,:),1);
                            for b = 1:size(phasmat,1)
                                stillcondpow(l,tf,sf,b) = squeeze(mean(mean(allpowresp(b,thisstillinds,respwin),2),3));
                            end
                            stillcondg1pow(l,tf,sf) = mean(gpower1(thisstillinds));
                            g1stillphasmat = g1(thisstillinds,:);
                            allstillg1phases = g1stillphasmat(find(g1stillphasmat));
                            if isrow(allstillg1phases) allstillg1phases = allstillg1phases'; end
                            stillcondg1ppc(l,tf,sf) = ppc(allstillg1phases);
                            stillcondg1cmean(l,tf,sf) = circ_mean(allstillg1phases);
                        else
                            stillcondresp(l,tf,sf,:) = nan(1,size(resp,2));
                            stillcondfr(l,tf,sf) = NaN;
                            stillconderr(l,tf,sf) = NaN;
                            stillclfpspect(l,tf,sf,:) = nan(1,size(lfpspect,2));
                            for b = 1:size(phasmat,1)
                                stillcondpow(l,tf,sf,b) = NaN;
                            end
                            stillcondg1pow(l,tf,sf) = NaN;
                            stillcondg1ppc(l,tf,sf) = NaN;
                            stillcondg1cmean(l,tf,sf) = NaN;
                        end
                        
                            [S,chf,Serr]=mtspectrumc(squeeze(lfpresp(thisinds,1001:1800))',params);
                            condS(cll,l,tf,sf,:) = S(1:150); condSerr(cll,l,tf,sf,:,:) = Serr(:,1:150);
                    
                            % running
                            thisruninds = intersect(thisinds,oktrials); thisstillinds = intersect(thisinds,stilltrials);
                            if ~isempty(thisruninds)
                                thisrunn(cll,l,tf,sf) = length(thisruninds);
                                [runS,chf,runSerr] = mtspectrumc(squeeze(lfpresp(thisruninds,1001:1800))',params);
                                condrunS(cll,l,tf,sf,:) = runS(1:150); condrunSerr(cll,l,tf,sf,:,:) = runSerr(:,1:150);
                                r1condallS(cll,l,tf,sf,:) = nanmean(allS(thisruninds,:),1);
                                r1condallSerr(cll,l,tf,sf,:) = nanstd(allS(thisruninds,:),1,1)./sqrt(length(thisruninds));

                            else
                                thisrunn(cll,l,tf,sf) = 0;
                                condrunS(cll,l,tf,sf,:) = nan(1,150); condrunSerr(cll,l,tf,sf,:,:) = nan(2,150);
                                r1condallS(cll,l,tf,sf,:) = nan(1,size(allS,2));
                                r1condallSerr(cll,l,tf,sf,:) = nan(1,size(allS,2));
                            end
                            if ~isempty(thisstillinds)
                                thisstilln(cll,l,tf,sf) = length(thisstillinds);
                                [stillS,chf,stillSerr] = mtspectrumc(squeeze(lfpresp(thisstillinds,1001:1800))',params);
                                condstillS(cll,l,tf,sf,:) = stillS(1:150); condstillSerr(cll,l,tf,sf,:,:) = stillSerr(:,1:150);
                                r0condallS(cll,l,tf,sf,:) = nanmean(allS(thisstillinds,:),1);
                                r0condallSerr(cll,l,tf,sf,:) = nanstd(allS(thisstillinds,:),1,1)./sqrt(length(thisstillinds));
                            else
                                thisstilln(cll,l,tf,sf) = 0;
                                condstillS(cll,l,tf,sf,:) = nan(1,150); condstillSerr(cll,l,tf,sf,:,:) = nan(2,150);
                                r0condallS(cll,l,tf,sf,:) = nan(1,size(allS,2));
                                r0condallSerr(cll,l,tf,sf,:) = nan(1,size(allS,2));
                            end


                            % field trip values
                            % field trip messing around
                            [ftspect(cll,l,tf,sf,:),ftphases{cll,l,tf,sf},ftfax,ftralp(cll,l,tf,sf,:),...
                                ftppc(cll,l,tf,sf,:),ftplv(cll,l,tf,sf,:)] = get_ft_spectstats(lfpresp,resp,thisinds);

                            % running only
                            [ftr1spect(cll,l,tf,sf,:),ftr1phases{cll,l,tf,sf},ftr1fax,ftr1ralp(cll,l,tf,sf,:),...
                                ftr1ppc(cll,l,tf,sf,:),ftr1plv(cll,l,tf,sf,:)] = get_ft_spectstats(lfpresp,resp,thisruninds);

                            % still only
                            [ftr0spect(cll,l,tf,sf,:),ftr0phases{cll,l,tf,sf},ftr0fax,ftr0ralp(cll,l,tf,sf,:),...
                                ftr0ppc(cll,l,tf,sf,:),ftr0plv(cll,l,tf,sf,:)] = get_ft_spectstats(lfpresp,resp,thisstillinds);
                            
                    end
                end
            end
            
            bincondresp = bincondresp.*(1000/binwidth);
            maxdriven = max(max(squeeze(cfr(1,:,:))))-mean(bl);             
            
            prefspfreq = find(mean(cfr(1,:,:),2) == max(mean(cfr(1,:,:),2)),1);
            preftfreq = find(mean(cfr(1,:,:),3) == max(mean(cfr(1,:,:),3)),1);            
            psf(cll) = spFreqs(prefspfreq);
            ptf(cll) = tFreqs(preftfreq);
            
            tfkwp(cll) = kruskalwallis(squeeze(allfr(1,:,prefspfreq,:))',[],'off');
            sffkwp(cll) = kruskalwallis(squeeze(allfr(1,preftfreq,:,:))',[],'off');
            condallfr(cll,:,:,:,:) = allfr;
            
            bincondcllresp(cll,:,:,:,:) = bincondresp;
            
            nok(cll) = length(oktrials);
            nstill(cll) = length(stilltrials);
            
            g1centerfreq(cll) = f(bpi);

            condfr(cll,:,:,:) = cfr;
            conderr(cll,:,:,:) = cerr;
            condmsl(cll,:,:,:) = cmsl;
            condgpow1(cll,:,:,:) = cgpow1;
            condlfpspect(cll,:,:,:,:) = clfpspect;
            condlfpspecterr(cll,:,:,:,:) = clfpspecterr; 
            
            condbandpow(cll,:,:,:,:) = condpow;
            condgamma1resp(cll,:,:,:,:) = cgamma1resp;
            
            r0condfr(cll,:,:,:) = stillcondfr;
            r0conderr(cll,:,:,:) = stillconderr;
            r1condfr(cll,:,:,:) = runcondfr;
            r1conderr(cll,:,:,:) = runconderr;  
            ntrialsr0(cll) = length(stilltrials);
            ntrialsr1(cll) = length(oktrials);  
            r0condlfpspect(cll,:,:,:,:) = stillclfpspect;
            r1condlfpspect(cll,:,:,:,:) = runclfpspect;
            runcondbandpow(cll,:,:,:,:) = runcondpow;
            stillcondbandpow(cll,:,:,:,:) = stillcondpow;
            r1g1pow(cll,:,:,:) = runcondg1pow;
            r0g1pow(cll,:,:,:) = stillcondg1pow;
            r1g1pcc(cll,:,:,:) = runcondg1ppc;
            r0g1pcc(cll,:,:,:) = stillcondg1ppc;
            r1g1cmean(cll,:,:,:) = runcondg1cmean;
            r0g1cmean(cll,:,:,:) = stillcondg1cmean;
            
            l0power(cll,:) = squeeze(mean(mean(allpowresp(:,find(result.light == 0),respwin),2),3));
            l1power(cll,:) = squeeze(mean(mean(allpowresp(:,find(result.light == 1),respwin),2),3));
            r1power(cll,:) = squeeze(mean(mean(allpowresp(:,oktrials,respwin),2),3));
            r0power(cll,:) = squeeze(mean(mean(allpowresp(:,stilltrials,respwin),2),3));
            
            l0lfp(cll,:) = squeeze(mean(lfpresp(find(result.light == 0),:)));
            l1lfp(cll,:) = squeeze(mean(lfpresp(find(result.light == 1),:)));
            r1lfp(cll,:) = squeeze(mean(lfpresp(oktrials,:)));
            r0lfp(cll,:) = squeeze(mean(lfpresp(stilltrials,:)));
            
            cllname{cll} = files(fi).name;
            animalno(cll) = animal(blck);
            depth(cll) = result.depth;
            recording(cll) = blck;
            pangle(cll) = penangle(blck);
         
            %determine if cll is modulated by light
            lightmod(cll) = ttest2(frs(find(result.light)),frs(find(result.light == 0)));
            vismod(cll) = ttest(frs(find(result.light == 0 & result.gratingInfo.Contrast ~= 0)),bl(find(result.light == 0 & result.gratingInfo.Contrast ~= 0)));
            if ~isnan(vismod(cll)) & vismod(cll) & mean(frs(find(result.light == 0 &...
                    result.gratingInfo.Contrast ~= 0)))<mean(bl(find(result.light == 0 & result.gratingInfo.Contrast ~= 0)))
                vismod(cll) = -1;
            end

            wvchan = find(var(result.waveforms) == max(var(result.waveforms)));
            spike = result.waveforms(:,wvchan);
            interpspike = spline(1:32,spike,1:.1:32);
            [adiff(cll),swidth(cll),ptr(cll),eslope(cll),fwhh(cll)] = spikequant(interpspike);
            
            waveform(cll,:) = spike;
            clustqual(cll) = result.clusterquality;
            
            bfr(cll) = mean(bl);
            lfr(cll) = mean(frs(find(result.light)));
            nlfr(cll) = mean(frs(find(result.light == 0)));
            
            alll0phasemat = tmpi(:,find(result.light == 0),:);
            alll1phasemat = tmpi(:,find(result.light == 1),:);
            for i = 1:size(allphaseresp,1)
                allphasesl0{i} = alll0phasemat(i,find(squeeze(alll0phasemat(i,:,:))));
                allphasesl1{i} = alll1phasemat(i,find(squeeze(alll1phasemat(i,:,:))));
                allrl0(cll,i) = circ_r(allphasesl0{i}');
                allrl1(cll,i) = circ_r(allphasesl1{i}');
                allcmeanl0(cll,i) = circ_mean(allphasesl0{i}');
                allcmeanl1(cll,i) = circ_mean(allphasesl1{i}');
            end
            
            for i = 1:20
                [ccoef,p] = corrcoef(squeeze(mean(allpowresp(i,:,respwin),3)),mean(speed(:,respwin),2));
                speedpowcc(cll,i) = ccoef(1,2);
                speedpowp(cll,i) = p(1,2);                
                
                [l0ccoef,l0p] = corrcoef(squeeze(mean(allpowresp(i,find(result.light ==0),respwin),3)),mean(speed(find(result.light == 0),respwin),2));
                l0speedpowcc(cll,i) = l0ccoef(1,2);
                l0speedpowp(cll,i) = l0p(1,2);
                
                [l1ccoef,l1p] = corrcoef(squeeze(mean(allpowresp(i,find(result.light ==1),respwin),3)),mean(speed(find(result.light == 1),respwin),2));
                l1speedpowcc(cll,i) = l1ccoef(1,2);
                l1speedpowp(cll,i) = l1p(1,2);
            end

            [binnedlight,bta] = binit(mean(lightresp),binwidth); binnedlight = binnedlight.*(1000/binwidth);
            [binnednolight,bta] = binit(mean(nolightresp),binwidth); binnednolight =  binnednolight.*(1000/binwidth);
            
            printname = files(fi).name;
            printname(find(printname=='_')) = ' ';

            nolbl = mean(bl(find(~result.light)));
            nolblerr = std(bl(find(~result.light)))./(sqrt(length(find(~result.light))));
            lbl = mean(bl(find(result.light)));
            lblerr = std(bl(find(result.light)))./(sqrt(length(find(result.light))));

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%

            %baselinesubtracted condfr
            blscondfr = squeeze(condfr(cll,:,:,:))-bfr(cll);

            %how much does the firing rate change for each condition
            condchange = squeeze(blscondfr(2,:,:))-squeeze(blscondfr(1,:,:));
            
%             figure
            subplot(2,2,1)
            ta = bta;
            plot(ta,binnednolight,'k','linewidth',2);
            hold on
            plot(ta,binnedlight,lcol,'linewidth',2);
            mx = max([max(binnednolight),max(binnedlight),0.1]);
            axis([-prestim,trialdur+poststim,0,mx]);
            line([0,0],[0,mx],'color','k','linewidth',2);
            line([2000,2000],[0,mx],'color','k','linewidth',2);
            %     line([500,500],[0,mx],'color','b','linewidth',2)
            %     line([1500,1500],[0,mx],'color','b','linewidth',2);
            %     legend({'Light OFF','Light ON'})
            xlabel('time [ms]')
            ylabel('firing rate [Hz]')
            title(['cell ' int2str(cll) ' depth: ' int2str(result.depth), 'cell ' printname ])
            
            subplot(2,2,2)
            errorbar(tFreqs,squeeze(condfr(cll,2,:,prefspfreq)),squeeze(conderr(cll,2,:,prefspfreq)),'o-','color',lcol,'markersize',8,'linewidth',2)
            hold on
            errorbar(tFreqs,squeeze(condfr(cll,1,:,prefspfreq)),squeeze(conderr(cll,1,:,prefspfreq)),'ko-','markersize',8,'linewidth',2)
            xlabel('shown temporal frequency')
            ylabel('Firing rate [Hz]')
            set(gca,'xtick',tFreqs)
            %     legend({'Light ON','Light OFF'})
            title([' TF Tuning at preferred SF '])
            
            subplot(2,2,3)
            errorbar(spFreqs,squeeze(condfr(cll,2,preftfreq,:)),...
                squeeze(conderr(cll,2,preftfreq,:)),'o-','color',lcol,'markersize',8,'linewidth',2);
            hold on
            errorbar(spFreqs,squeeze(condfr(cll,1,preftfreq,:)),...
                squeeze(conderr(cll,1,preftfreq,:)),'ko-','markersize',8,'linewidth',2);
            xlabel('shown spfreq')
            ylabel('Firing rate [Hz]')
            %     legend({'Light ON','Light OFF'})
            set(gca,'xtick',spFreqs)
            title(['SF Tuning at preferred TF'])
            
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
                        
            cll = cll + 1;
            disp([files(fi).name '   done'])

        end
    end
    
    save(popfile, '-v7.3');
else
    load(popfile);   
end

toc

%spike classification
kmeansind = kmeans([eslope',ptr',swidth',adiff'],3);
% kmeansind = kmeans([eslope',adiff'],2);
secpersamp = 1/30000;
interpf = secpersamp/10;
swidthms = swidth*interpf*1000;

m1 = mean(swidth(kmeansind==1)); m2 = mean(swidth(kmeansind==2)); m3 = mean(swidth(kmeansind==3));
fm = min([m1,m2,m3]); fmax = max([m1,m2,m3]);

%three clusters, use only extremes
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

% three clusters, two wide ones RS
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

% %just two clusters
% if mean(swidth(find(kmeansind==1)))<mean(swidth(find(kmeansind==2)))  %1 is FS
%     pfs = find(kmeansind==1); prs = find(kmeansind==2); pfsv = kmeansind==1; prsv = kmeansind==2;
% else
%     pfs = find(kmeansind==2); prs = find(kmeansind==1); pfsv = kmeansind==2; prsv = kmeansind==1;
% end

figure
plot(swidthms(kmeansind==1),adiff(kmeansind==1),'b.')
xlabel('spike width')
ylabel('amplitude diff')
hold on
plot(swidthms(kmeansind==2),adiff(kmeansind==2),'r.')
plot(swidthms(kmeansind==3),adiff(kmeansind==3),'g.')
plot(swidthms(pfsv),adiff(pfsv),'ro');
plot(swidthms(prsv),adiff(prsv),'o');
% axis([5.5,20.5,-.9,.7])

swamp = max(waveform,[],2)-min(waveform,[],2);
okwv = swamp>42;

vismod(isnan(vismod)) = 0;
% ok = okwv'&vismod;
% ok = okwv'&nlfr>1;
ok = okwv';%&nlfr>.5;
% ok = vismod&okwv'&nlfr>1;
% ok(460) = 0; ok(288) = 0; ok(284) = 0;

% adjust depth according to penetration angle 22 in coronal plane and penangle in saggital
depth = depth.*cosd(22).*cosd(pangle);
phe = zeros(1,length(depth));

prsv = zeros(1,length(depth));
% prsv(swidth>15) = 1;
prsv(prs) = 1; 
prsv = logical(prsv'); phe = logical(phe)';

% % for SOM later pop, around 430, none, 320, 320, none and 340
% lfpinds = [31, 82, 94, 114];

% for SOM and control pop - shoot for 336 - animal 4 no cells
lfpinds = [7, 34, 47, 53, 64, 76];


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
l6 = depth>800;
l23rs = l23&prsv'&~phe'&ok;
l23fs = l23&pfsv'&~phe'&ok;
l4rs = l4&prsv'&~phe'&ok;
l4fs = l4&pfsv'&~phe'&ok;
l5rs = l5&prsv'&~phe'&ok;
l5fs = l5&pfsv'&~phe'&ok;
l6rs = l6&prsv'&~phe'&ok;

end


function [ftspect,ftphases,ftfax,ftralp,ftppc,ftplv] = get_ft_spectstats(lfpresp,resp,thisinds)
 
     if isempty(thisinds) || isempty(find(resp(thisinds,801:1800)))
         ftspect = nan(1,30); ftphases = NaN; ftfax = nan(1,30);
         ftralp = nan(1,30); ftppc = nan(1,30); ftplv = nan(1,30);
     else
         clear data;
         for i = 1:length(thisinds)
             data.trial{i} = [lfpresp(thisinds(i),:);resp(thisinds(i),:)];
             data.time{i} = -.299:.001:2.700;
         end
         data.fsample = 1000;
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
         ftspect = squeeze(nanmean(mag(:,1,:)));
         ftphases = squeeze(ang);
         ftfax = stsFFT.freq;

         cfg               = [];
         cfg.method        = 'ral'; % compute the rayleigh test
         cfg.spikechannel  = stsFFT.label{1};
         cfg.channel       = stsFFT.lfplabel; % selected LFP channels
         cfg.avgoverchan   = 'unweighted'; % weight spike-LFP phases irrespective of LFP power
         cfg.timwin        = 'all'; % compute over all available spikes in the window
         cfg.latency       = [0.5 1.5]; % sustained visual stimulation period
         statSts           = ft_spiketriggeredspectrum_stat(cfg,stsFFT);
         ftralp = statSts.ral;
         cfg.method = 'ppc0';
         statSts           = ft_spiketriggeredspectrum_stat(cfg,stsFFT);
         ftppc = statSts.ppc0;
         cfg.method = 'plv';
         statSts           = ft_spiketriggeredspectrum_stat(cfg,stsFFT);
         ftplv = statSts.plv;
     end
end
