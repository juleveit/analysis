function orionly_populationana


% SOM later population
animalids = {'150529','150602','150603','160623','150625','150825','150831','150902','150909','150915','150916','150916','151023','151109','151109','151110','151209','151221','160422'};
blocks    = [6,        8,       5,       8,       11,      6,       5,       7,       5,       5,       3,       3,       12,      13,      13,      14,      9,       11,      5];
animal    = [1,        2,       3,       4,       5,       6,       7,       8,       9,       10,      11,      11,      12,      13,      13,      14,      15,      16,      17];
electrodes =[[1,32];  [1,32];  [1,16];  [1,16];  [1,16];  [17,32]; [1,16];  [1,16];  [1,16];  [17,32]; [1,16];  [17,32]; [17,32]; [1,16];  [17,32]; [1,16];  [17,32]; [1,16];  [1,16]];
penangle =  [25,       25,      10,      10,      10,      25,      25,      25,      25,      25,      25,      25,      25,      25,      25,      25,      25,      25,      25];
printpath = 'C:\Users\Julia\work\data\populations\SOM_Halo_later\orientation\units\';
runprintpath = 'C:\Users\Julia\work\data\populations\SOM_Halo_later\orientation\running\';
condprintpath = 'C:\Users\Julia\work\data\populations\SOM_Halo_later\orientation\conditions\';
popfile = 'C:\Users\Julia\work\data\populations\SOM_Halo_later\orientation\orionly_population.mat';


tic
lcol = 'r'; %lasercolor

recalculate = 1;
printyn = 1;

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
        offsetwin = 1501:2500;
        respwin = respwin+prestim;
        offsetwin = offsetwin+prestim;
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
            
            isi = diff(msStimes);
            bursts = legendy_new3(isi,4,1000,3,0,15); %(ISI, fac, sr, min_length_of_burst, local_length, surprise_cutoff)
            surp = [bursts.surprise];
            bursts(surp == 100) = []; % delete probably wrong bursts
            burstbegs = [bursts.begin];
            burstchan = zeros(1,length(result.lfp));
            burstchan(msStimes(burstbegs)) = 1;
            
            wvchan = find(var(result.waveforms) == max(var(result.waveforms)));
            cm = [3,4,1,2]; % confusion matrix? take lfp from two electrodes away to not get too many spike related phase resets
            lfp = result.lfp(:,cm(wvchan))';      
            
            trialdur = result.stimduration*1000;
            msstamps = result.msstamps;
            
            if length(msstamps)~=length(result.light)
%             disp('');
%             msstamps(16) = []; % for 150414 block 10
%             result.msstamps = msstamps;
%             save([supath, files(fi).name],'result');            
                pause;
            end
            
            for i = 1:length(msstamps)
                speed(i,:) = result.runspeed(msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim);
            end
            
            meanspeed = mean(speed(:,respwin),2);
            stdspeed = nanstd(speed(:,respwin),1,2);
            notstill = find(meanspeed>1);
            okspeed = find(meanspeed>( mean(meanspeed(notstill))-(1.5*nanstd(meanspeed(notstill))) ) & meanspeed>1);
            okvar = find(stdspeed<( mean(stdspeed(notstill))+(1.5*nanstd(stdspeed(notstill)))) & stdspeed>.5);
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
            high = find(result.gratingInfo.size == max(unique(result.gratingInfo.size)) & result.light == 0);
            low = find(result.gratingInfo.size == min(unique(result.gratingInfo.size)) & result.light == 0);
            lr1 = intersect(high,oktrials);
            sr1 = intersect(low,oktrials);
            for i = 1:length(low)
                [pl(i,:),f] = pmtm(lfp(result.msstamps(high(i)):result.msstamps(high(i))+1000),3,[],1000);
                [ps(i,:),f] = pmtm(lfp(result.msstamps(low(i)):result.msstamps(low(i))+1000),3,[],1000);
            end      
            plr1 = nan(1,size(pl,2)); psr1 = nan(1,size(pl,2));
            if length(lr1>=5)
                for i = 1:length(lr1)
                    [plr1(i,:),f] = mtspectrumc(lfp(result.msstamps(lr1(i))+700:result.msstamps(lr1(i))+1500),params);
%                     [plr1(i,:),f] = pmtm(lfp(result.msstamps(lr1(i))+300:result.msstamps(lr1(i))+1300),3,[],1000);
                end
            end
            if length(sr1>=5)
                for i = 1:length(sr1)
                    [psr1(i,:),f] = mtspectrumc(lfp(result.msstamps(sr1(i))+700:result.msstamps(sr1(i))+1500),params); 
                    y(i,:) = fft(lfp(result.msstamps(sr1(i)):result.msstamps(sr1(i))+2000),2048);
%                     [psr1(i,:),f] = pmtm(lfp(result.msstamps(sr1(i))+300:result.msstamps(sr1(i))+1300),3,[],1000);
                end
            end          
            b1 = find(f>beta(1),1); b2 = find(f>beta(2),1);
            g1 = find(f>gamma(1),1); g2 = find(f>gamma(2),1);
            
            if ~isnan(plr1(1))
                bsig = nanmean(plr1(:,b1:b2));
            else
                bsig = nanmean(pl(:,b1:b2));
            end
            if ~isnan(psr1(1))
                gsig = nanmean(psr1(:,g1:g2));
            else
                gsig = nanmean(ps(:,g1:g2));
            end
            
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
            gamma3 = eegfilt(lfp,sr,f(gpi)+2.5,100);
            h1 = hilbert(gamma1); gpow1 = abs(h1); gphas1 = angle(h1);
            h2 = hilbert(gamma2); gpow2 = abs(h2); gphas2 = angle(h2);
            h3 = hilbert(gamma3); gpow3 = abs(h3); gphas3 = angle(h3);
            
            for i = 1:100/freqbinwidth
                filtmat(i,:) = eegfilt(lfp,sr,(i-1)*freqbinwidth+1,i*freqbinwidth);
                h = hilbert(filtmat(i,:));
                powmat(i,:) = abs(h); phasmat(i,:) = angle(h);
            end
           
            
            clear lfpspect; clear lfpoffsspect; clear fax;            
            for i = 1:length(msstamps)
                resp(i,:) = chan(msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim);
                burstresp(i,:) = burstchan(msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim); 
                lfpresp(i,:) = result.lfp(msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim,wvchan);
                gamma1resp(i,:) = gpow1(msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim); 
                gamma2resp(i,:) = gpow2(msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim); 
                gamma3resp(i,:) = gpow3(msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim); 
                gamma1phasresp(i,:) = gphas1(msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim); 
                gamma2phasresp(i,:) = gphas2(msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim); 
                gamma3phasresp(i,:) = gphas3(msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim); 
                speed(i,:) = result.runspeed(msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim);
                [lfpspect(i,:),fax] = pmtm(lfpresp(i,respwin),2,[],1000);
                lfpoffsspect(i,:) = pmtm(lfpresp(i,offsetwin),3,[],1000);
                
                allS(i,:)=mtspectrumc(squeeze(lfpresp(i,1001:1800))',params);
                
                for j = 1:size(phasmat,1)
                    allphaseresp(j,i,:) = phasmat(j, msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim);
                    allpowresp(j,i,:) = powmat(j, msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim);
                end
                
                isiresp{i} = diff(find(resp(i,respwin)));
                msl(i) = nanmean(find(resp(i,respwin))); %mean spike latency

            end
            
            lfpspect = lfpspect(:,1:150); fax = fax(1:150);
            lfpoffsspect = lfpoffsspect(:,1:150);
            
            msta = linspace(-prestim,trialdur+poststim,size(resp,2));

            lightresp = resp(find(result.light),:);
            nolightresp = resp(find(result.light == 0),:);
            lightlfpresp = lfpresp(find(result.light),:);
            nolightlfpresp = lfpresp(find(~result.light),:);

            frs = sum(resp(:,respwin),2)./(length(respwin)/1000);
            brs = sum(burstresp(:,respwin),2);
            bl = sum(resp(:,1:prestim),2)./(prestim/1000);
            sc = sum(resp(:,respwin),2);
            gpower1 = nanmean(gamma1resp(:,respwin),2);
            gpower2 = nanmean(gamma2resp(:,respwin),2);
            gpower3 = nanmean(gamma3resp(:,respwin),2);

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
            g2 = zeros(size(gamma2phasresp));
            g2(find(resp)) = gamma2phasresp(find(resp));
            g3 = zeros(size(gamma3phasresp));
            g3(find(resp)) = gamma3phasresp(find(resp));
            
            
            binwidth = 20;
            condsize(cll) = max(unique(result.gratingInfo.size));
            oris = unique(result.gratingInfo.Orientation); oris(find(oris == -1)) = [];
            for l = 1:length(unique(result.light))
                for ori = 1:length(oris)
                    thisinds = find(result.gratingInfo.Orientation == oris(ori) & result.light == l-1);
                    
                    
                    [S,chf,Serr]=mtspectrumc(squeeze(lfpresp(thisinds,1001:1800))',params);
                    condS(cll,l,ori,:) = S(1:150); condSerr(cll,l,ori,:,:) = Serr(:,1:150);
                    
                    % running
                    thisruninds = intersect(thisinds,oktrials); thisstillinds = intersect(thisinds,stilltrials);
                    if ~isempty(thisruninds)
                        thisrunn(cll,l,ori) = length(thisruninds);
                        [runS,chf,runSerr] = mtspectrumc(squeeze(lfpresp(thisruninds,1001:1800))',params);
                        condrunS(cll,l,ori,:) = runS(1:150); condrunSerr(cll,l,ori,:,:) = runSerr(:,1:150);
                        r1condallS(cll,l,ori,:) = nanmean(allS(thisruninds,:),1);
                        r1condallSerr(cll,l,ori,:) = nanstd(allS(thisruninds,:),1,1)./sqrt(length(thisruninds));
                        
                    else
                        thisrunn(cll,l,ori) = 0;
                        condrunS(cll,l,ori,:) = nan(1,150); condrunSerr(cll,l,ori,:,:) = nan(2,150);
                        r1condallS(cll,l,ori,:) = nan(1,size(allS,2));
                        r1condallSerr(cll,l,ori,:) = nan(1,size(allS,2));
                    end
                    if ~isempty(thisstillinds)
                        thisstilln(cll,l,ori) = length(thisstillinds);
                        [stillS,chf,stillSerr] = mtspectrumc(squeeze(lfpresp(thisstillinds,1001:1800))',params);
                        condstillS(cll,l,ori,:) = stillS(1:150); condstillSerr(cll,l,ori,:,:) = stillSerr(:,1:150);
                        r0condallS(cll,l,ori,:) = nanmean(allS(thisstillinds,:),1);
                        r0condallSerr(cll,l,ori,:) = nanstd(allS(thisstillinds,:),1,1)./sqrt(length(thisstillinds));
                    else
                        thisstilln(cll,l,ori) = 0;
                        condstillS(cll,l,ori,:) = nan(1,150); condstillSerr(cll,l,ori,:,:) = nan(2,150);
                        r0condallS(cll,l,ori,:) = nan(1,size(allS,2));
                        r0condallSerr(cll,l,ori,:) = nan(1,size(allS,2));
                    end
                    
                    % field trip values
                    % field trip messing around
                    [ftspect(cll,l,ori,:),ftphases{cll,l,ori},ftfax,ftralp(cll,l,ori,:),...
                        ftppc(cll,l,ori,:),ftplv(cll,l,ori,:)] = get_ft_spectstats(lfpresp,resp,thisinds);
                    
                    % running only
                    [ftr1spect(cll,l,ori,:),ftr1phases{cll,l,ori},ftr1fax,ftr1ralp(cll,l,ori,:),...
                        ftr1ppc(cll,l,ori,:),ftr1plv(cll,l,ori,:)] = get_ft_spectstats(lfpresp,resp,thisruninds);
                    
                    % still only
                    [ftr0spect(cll,l,ori,:),ftr0phases{cll,l,ori},ftr0fax,ftr0ralp(cll,l,ori,:),...
                        ftr0ppc(cll,l,ori,:),ftr0plv(cll,l,ori,:)] = get_ft_spectstats(lfpresp,resp,thisstillinds);
                    
                    
                    
                    condn(cll,l,ori) = length(thisinds);
                    condresp(cll,l,ori,:) = nanmean(resp(thisinds,:),1);
                    condgamma1resp(cll,l,ori,:) = nanmean(gamma1resp(thisinds,:),1);
                    condgamma2resp(cll,l,ori,:) = nanmean(gamma2resp(thisinds,:),1);
                    condgamma3resp(cll,l,ori,:) = nanmean(gamma3resp(thisinds,:),1);
                    condlfpresp(cll,l,ori,:) = nanmean(lfpresp(thisinds,:));
                    condlfpspect(cll,l,ori,:) = nanmean(lfpspect(thisinds,:));
                    condlfpspecterr(cll,l,ori,:) = nanstd(lfpspect(thisinds,:))./sqrt(length(thisinds));
                    condlfpoffsspect(cll,l,ori,:) = nanmean(lfpoffsspect(thisinds,:));
                    condgpow1(cll,l,ori) = nanmean(gpower1(thisinds));
                    condgpow2(cll,l,ori) = nanmean(gpower2(thisinds));
                    condgpow3(cll,l,ori) = nanmean(gpower3(thisinds));
                    condfr(cll,l,ori) = nanmean(frs(thisinds));
                    condbr(cll,l,ori) = nanmean(brs(thisinds));
                    allfr{cll,l,ori} = frs(thisinds);
                    conderr(cll,l,ori) =nanstd(frs(thisinds))./sqrt(length(thisinds));
                    [bincondresp(cll,l,ori,:),bta] = binit(squeeze(condresp(cll,l,ori,:)),binwidth);
                    condmsl(cll,l,ori) = nanmean(msl(thisinds)); % mean spike latency
                    
                    lfptbtvar(cll,l,ori) = mean(var(lfpresp(thisinds,respwin)),2);
                    frtbtvar(cll,l,ori) = var(frs(thisinds));
                    
                    condisi{cll,l,ori} = [];
                    for i = 1:length(thisinds)
                        condisi{cll,l,ori} = [condisi{cll,l,ori},isiresp{thisinds(i)}];
                    end
                    % test each cell for each condition ranksum
                    if l == 1
                        l1inds = find(result.gratingInfo.Orientation == oris(ori) &...
                            result.light == 1);
                        if ~isempty(l1inds)
                            cdiffp(ori) = ranksum(frs(thisinds),frs(l1inds));
                        else
                            cdiffp(ori) = NaN;
                        end
                    end
                    
                    condz{cll,l,ori} = {(sc(thisinds)-mean(sc(thisinds)))/nanstd(sc(thisinds))}; %ecker 2010
                    condsc(cll,l,ori) = {sc(thisinds)};
                    condff(cll,l,ori) = var(sc(thisinds))/mean(sc(thisinds));                   
                    
                    mscc = []; bincc = [];
                    for ii = 1:length(thisinds)-1
                        for jj = ii+1:length(thisinds)
                            help = corrcoef(resp(thisinds(ii),:),resp(thisinds(jj),:));
                            mscc = [mscc,help(1,2)];
                            help = corrcoef(binit(resp(thisinds(ii),:),binwidth),binit(resp(thisinds(jj),:),binwidth));
                            bincc = [bincc, help(1,2)];
                        end
                    end
                    msrely(cll,l,ori) = nanmean(mscc);
                    binrely(cll,l,ori) = nanmean(bincc);
                    eckerrely(cll,l,ori) = var(frs(thisinds))/var(frs);
            
                    for b = 1:size(phasmat,1)
                        condpowresp(cll,l,ori,b,:) = squeeze(nanmean(allpowresp(b,thisinds,:),2));
                        condbandpow(cll,l,ori,b) = squeeze(nanmean(nanmean(allpowresp(b,thisinds,respwin),2),3));
                    end
                    
                    %phase locking to all bands and the two peak bands
                    phasemat = tmpi(:,thisinds,:);
                    for b = 1:size(allphaseresp,1)
                        allphases = squeeze(phasemat(b,find(squeeze(phasemat(b,:,:)))));
                        condr(cll,b,l,ori) = circ_r(allphases');
                        condppc(cll,b,l,ori) = ppc(allphases);
                        condcmean(cll,b,l,ori) = circ_mean(allphases');
                    end
                    g1phasmat = g1(thisinds,:);
                    g2phasmat = g2(thisinds,:);
                    g3phasmat = g3(thisinds,:);
                    allg1phases = squeeze(g1phasmat(find(g1phasmat)));
                    allg2phases = squeeze(g2phasmat(find(g2phasmat)));
                    allg3phases = squeeze(g3phasmat(find(g2phasmat)));
                    if isrow(allg1phases) allg1phases = allg1phases'; end
                    if isrow(allg2phases) allg2phases = allg2phases'; end
                    if isrow(allg3phases) allg3phases = allg3phases'; end
                    condg1r(cll,l,ori) = circ_r(allg1phases);
                    condg1ppc(cll,l,ori) = ppc(allg1phases);
                    condg1cmean(cll,l,ori) = circ_mean(allg1phases);
                    condg2r(cll,l,ori) = circ_r(allg2phases);
                    condg2ppc(cll,l,ori) = ppc(allg2phases);
                    condg2cmean(cll,l,ori) = circ_mean(allg2phases);
                    condg3r(cll,l,ori) = circ_r(allg3phases);
                    condg3ppc(cll,l,ori) = ppc(allg3phases);
                    condg3cmean(cll,l,ori) = circ_mean(allg3phases);
                    
                    
                    %running
                    thisruninds = intersect(thisinds,oktrials);
                    if ~isempty(thisruninds)
                        r1condresp(cll,l,ori,:) = nanmean(resp(thisruninds,:),1);
                        r1condfr(cll,l,ori) = nanmean(frs(thisruninds),1);
                        r1conderr(cll,l,ori) = nanstd(frs(thisruninds))./sqrt(length(thisruninds));
                        r1clfpspect(cll,l,ori,:) = nanmean(lfpspect(thisruninds,:),1);
                        for b = 1:size(phasmat,1)
                            r1condbandpow(cll,l,ori,b) = squeeze(nanmean(nanmean(allpowresp(b,thisruninds,respwin),2),3));
                        end
                        r1g1pow(cll,l,ori) = nanmean(gpower1(thisruninds));
                        r1g2pow(cll,l,ori) = nanmean(gpower2(thisruninds));
                        r1g3pow(cll,l,ori) = nanmean(gpower3(thisruninds));
                        g1runphasmat = g1(thisruninds,:);
                        g2runphasmat = g2(thisruninds,:);
                        g3runphasmat = g3(thisruninds,:);
                        allrung1phases = squeeze(g1runphasmat(find(g1runphasmat)));
                        allrung2phases = squeeze(g2runphasmat(find(g2runphasmat)));
                        allrung3phases = squeeze(g3runphasmat(find(g3runphasmat)));
                        if isrow(allrung1phases) allrung1phases = allrung1phases'; end
                        if isrow(allrung2phases) allrung2phases = allrung2phases'; end
                        if isrow(allrung3phases) allrung3phases = allrung3phases'; end
                        runcondg1ppc(l,ori) = ppc(allrung1phases);
                        runcondg1cmean(l,ori) = circ_mean(allrung1phases);
                        runcondg2ppc(l,ori) = ppc(allrung2phases);
                        runcondg2cmean(l,ori) = circ_mean(allrung2phases);
                        runcondg3ppc(l,ori) = ppc(allrung3phases);
                        runcondg3cmean(l,ori) = circ_mean(allrung3phases);
                    else
                        r1condresp(cll,l,ori,:) = nan(1,size(resp,2));
                        r1condfr(cll,l,ori) = NaN;
                        r1conderr(cll,l,ori) = NaN;
                        r1clfpspect(cll,l,ori,:) = nan(1,size(lfpspect,2));
                        for b = 1:size(phasmat,1)
                            r1condbandpow(cll,l,ori,b) = NaN;
                        end
                        r1g1pow(cll,l,ori) = NaN;
                        r1g2pow(cll,l,ori) = NaN;
                        r1g3pow(cll,l,ori) = NaN;
                        r1g1ppc(cll,l,ori) = NaN;
                        r1g1cmean(cll,l,ori) = NaN;
                        r1g2ppc(cll,l,ori) = NaN;
                        r1g2cmean(cll,l,ori) = NaN;
                        r1g3ppc(cll,l,ori) = NaN;
                        r1g3cmean(cll,l,ori) = NaN;
                    end
                    
                    thisstillinds = intersect(thisinds,stilltrials);
                    if ~isempty(thisstillinds)
                        r0condresp(cll,l,ori,:) = nanmean(resp(thisstillinds,:),1);
                        r0condfr(cll,l,ori) = nanmean(frs(thisstillinds),1);
                        r0conderr(cll,l,ori) = nanstd(frs(thisstillinds),1)./sqrt(length(thisstillinds));
                        r0clfpspect(cll,l,ori,:) = nanmean(lfpspect(thisstillinds,:),1);
                        for b = 1:size(phasmat,1)
                            r0condbandpow(cll,l,ori,b) = squeeze(nanmean(nanmean(allpowresp(b,thisstillinds,respwin),2),3));
                        end
                        r0g1pow(cll,l,ori) = nanmean(gpower1(thisstillinds));
                        r0g2pow(cll,l,ori) = nanmean(gpower2(thisstillinds));
                        r0g3pow(cll,l,ori) = nanmean(gpower3(thisstillinds));
                        g1stillphasmat = g1(thisstillinds,:);
                        g2stillphasmat = g2(thisstillinds,:);
                        g3stillphasmat = g3(thisstillinds,:);
                        allstillg1phases = g1stillphasmat(find(g1stillphasmat));
                        allstillg2phases = g2stillphasmat(find(g2stillphasmat));
                        allstillg3phases = g3stillphasmat(find(g3stillphasmat));
                        if isrow(allstillg1phases) allstillg1phases = allstillg1phases'; end
                        if isrow(allstillg2phases) allstillg2phases = allstillg2phases'; end
                        if isrow(allstillg3phases) allstillg3phases = allstillg3phases'; end
                        r0g1ppc(cll,l,ori) = ppc(allstillg1phases);
                        r0g1cmean(cll,l,ori) = circ_mean(allstillg1phases);
                        r0g2ppc(cll,l,ori) = ppc(allstillg2phases);
                        r0g2cmean(cll,l,ori) = circ_mean(allstillg2phases);
                        r0g3ppc(cll,l,ori) = ppc(allstillg3phases);
                        r0g3cmean(cll,l,ori) = circ_mean(allstillg3phases);
                    else
                        r0condresp(cll,l,ori,:) = nan(1,size(resp,2));
                        r0condfr(cll,l,ori) = NaN;
                        r0conderr(cll,l,ori) = NaN;
                        r0clfpspect(cll,l,ori,:) = nan(1,size(lfpspect,2));
                        for b = 1:size(phasmat,1)
                            r0condbandpow(cll,l,ori,b) = NaN;
                        end
                        r0g1pow(cll,l,ori) = NaN;
                        r0g2pow(cll,l,ori) = NaN;
                        r0g3pow(cll,l,ori) = NaN;
                        r0g1ppc(cll,l,ori) = NaN;
                        r0g1cmean(cll,l,ori) = NaN;
                        r0g2ppc(cll,l,ori) = NaN;
                        r0g2cmean(cll,l,ori) = NaN;
                        r0g3ppc(cll,l,ori) = NaN;
                        r0g3cmean(cll,l,ori) = NaN;
                    end                    
                end                
                [oriprefratio(cll,l), dirprefratio(cll,l), po, meanori(cll,l), osi(cll,l), meandir(cll,l), dsi(cll,l)] = getOSI(squeeze(condfr(cll,l,:))',oris);
                [g1oriprefratio(cll,l), g1dirprefratio(cll,l), po, g1meanori(cll,l), g1osi(cll,l), g1meandir(cll,l), g1dsi(cll,l)] = getOSI(squeeze(condgpow1(cll,l,:))',oris);
                [g2oriprefratio(cll,l), g2dirprefratio(cll,l), po, g2meanori(cll,l), g2osi(cll,l), g2meandir(cll,l), g2dsi(cll,l)] = getOSI(squeeze(condgpow2(cll,l,:))',oris);            
            end           
            bincondresp = bincondresp.*(1000/binwidth);
            
            
            % now bootstrap OSI significance
            for jj = 1:10000
                a = randperm(size(resp,1));
                shuffr = frs(a);
                shufg1 = gpower1(a);
                shufg2 = gpower2(a);
                for l = 1:length(unique(result.light))
                    for ori = 1:length(oris)
                        thisinds = find(result.gratingInfo.Orientation == oris(ori) &...
                            result.light == l-1);
                        cshuffr(l,ori) = nanmean(shuffr(thisinds));
                        cshufg1(l,ori) = nanmean(shufg1(thisinds));
                        cshufg2(l,ori) = nanmean(shufg2(thisinds));
                    end
                    [x, x, x, x, shufosi(jj,l), x, x] = getOSI(squeeze(cshuffr(l,:)),oris);
                    [x, x, po, x, shufg1osi(jj,l), x, x] = getOSI(squeeze(cshufg1(l,:)),oris);
                    [x, x, po, x, shufg2osi(jj,l), x, x] = getOSI(squeeze(cshufg2(l,:)),oris);
                end
            end
            for l = 1:length(unique(result.light))
                hhh = sort(shufosi(:,l));
                osiconf(cll,l) = hhh(9500);
                osisig(cll,l) = hhh(9500)<osi(cll,l);
                hhh = sort(shufg1osi(:,l));
                g1osiconf(cll,l) = hhh(9500);
                g1osisig(cll,l) = hhh(9500)<g1osi(cll,l);
                hhh = sort(shufg2osi(:,l));
                g2osiconf(cll,l) = hhh(9500);
                g2osisig(cll,l) = hhh(9500)<g2osi(cll,l);
            end
            
            l0isis = []; l1isis = [];
            for i = find(result.light == 0)
                l0isis = [l0isis, isiresp{i}];
            end
            for i = find(result.light == 1)
                l1isis = [l1isis, isiresp{i}];
            end           
            l0isi{cll} = l0isis; l1isi{cll} = l1isis;
            
            ntrialsr0(cll) = length(stilltrials);
            ntrialsr1(cll) = length(oktrials); 
            
            g1centerfreq(cll) = f(bpi);
            g2centerfreq(cll) = f(gpi);
            
            nbursts(cll) = length(bursts);
            nspikes{cll} = [bursts.num_spikes];
            
            l0power(cll,:) = squeeze(nanmean(nanmean(allpowresp(:,find(result.light == 0),respwin),2),3));
            l1power(cll,:) = squeeze(nanmean(nanmean(allpowresp(:,find(result.light == 1),respwin),2),3));
            r1power(cll,:) = squeeze(nanmean(nanmean(allpowresp(:,oktrials,respwin),2),3));
            r0power(cll,:) = squeeze(nanmean(nanmean(allpowresp(:,stilltrials,respwin),2),3));
            
            l0lfp(cll,:) = squeeze(nanmean(lfpresp(find(result.light == 0),:)));
            l1lfp(cll,:) = squeeze(nanmean(lfpresp(find(result.light == 1),:)));
            r1lfp(cll,:) = squeeze(nanmean(lfpresp(oktrials,:)));
            r0lfp(cll,:) = squeeze(nanmean(lfpresp(stilltrials,:)));
            
            cllname{cll} = files(fi).name;
            animalno(cll) = animal(blck);
            depth(cll) = result.depth;
            recording(cll) = blck;
            pangle(cll) = penangle(blck);
            
            %determine if cll is modulated by light
            [p,lightmod(cll)] = signrank(frs(find(result.light)),frs(find(result.light == 0)));
            [p,vismod(cll)] = signrank(frs(find(result.light == 0 & result.gratingInfo.Orientation ~= -1)),bl(find(result.light == 0 & result.gratingInfo.Orientation ~= -1)));
            if ~isnan(vismod(cll)) & vismod(cll) & nanmean(frs(find(result.light == 0 &...
                    result.gratingInfo.Orientation ~= -1)))<nanmean(bl(find(result.light == 0 & result.gratingInfo.Orientation ~= -1)))
                vismod(cll) = -1;
            end

            wvchan = find(var(result.waveforms) == max(var(result.waveforms)));
            spike = result.waveforms(:,wvchan);
            interpspike = spline(1:32,spike,1:.1:32);
            [adiff(cll),swidth(cll),ptr(cll),eslope(cll),fwhh(cll)] = spikequant(interpspike);
            
            waveform(cll,:) = spike;
            clustqual(cll) = result.clusterquality;
            
            bfr(cll) = nanmean(bl);
            lfr(cll) = nanmean(frs(find(result.light)));
            nlfr(cll) = nanmean(frs(find(result.light == 0)));
            
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
                [ccoef,p] = corrcoef(squeeze(nanmean(allpowresp(i,:,respwin),3)),nanmean(speed(:,respwin),2));
                speedpowcc(cll,i) = ccoef(1,2);
                speedpowp(cll,i) = p(1,2);                
                
                [l0ccoef,l0p] = corrcoef(squeeze(nanmean(allpowresp(i,find(result.light ==0),respwin),3)),nanmean(speed(find(result.light == 0),respwin),2));
                l0speedpowcc(cll,i) = l0ccoef(1,2);
                l0speedpowp(cll,i) = l0p(1,2);
                
                [l1ccoef,l1p] = corrcoef(squeeze(nanmean(allpowresp(i,find(result.light ==1),respwin),3)),nanmean(speed(find(result.light == 1),respwin),2));
                l1speedpowcc(cll,i) = l1ccoef(1,2);
                l1speedpowp(cll,i) = l1p(1,2);
            end

            [binnedlight,bta] = binit(nanmean(lightresp),binwidth); binnedlight = binnedlight.*(1000/binwidth);
            [binnednolight,bta] = binit(nanmean(nolightresp),binwidth); binnednolight =  binnednolight.*(1000/binwidth);
            
            printname = files(fi).name;
            printname(find(printname=='_')) = ' ';

            nolbl = nanmean(bl(find(~result.light)));
            nolblerr = nanstd(bl(find(~result.light)))./(sqrt(length(find(~result.light))));
            lbl = nanmean(bl(find(result.light)));
            lblerr = nanstd(bl(find(result.light)))./(sqrt(length(find(result.light))));

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%

            %baselinesubtracted condfr
            blscondfr(cll,:,:) = squeeze(condfr(cll,:,:))-bfr(cll);

            %how much does the BLS firing rate change for each condition
            condchange = squeeze(blscondfr(cll,2,:))-squeeze(blscondfr(cll,1,:));

            contindsnl = find(result.gratingInfo.Orientation == -1 & result.light == 0);
            controlresp(1,:) = nanmean(resp(contindsnl,:),1);
            controlfr(cll,1) = nanmean(frs(contindsnl));
            blscontrolfr(cll,1) = controlfr(cll,1)-bfr(cll);
            controlerr(cll,1) = nanstd(frs(contindsnl))./sqrt(length(contindsnl));
            controllfpspect(cll,1,:) = nanmean(lfpspect(contindsnl,:));
            controllfpspecterr(cll,1,:) = nanstd(lfpspect(contindsnl,:))./sqrt(length(contindsnl));
            controlgamma1pow(cll,1) = nanmean(gpower1(contindsnl,:));
            controlgamma2pow(cll,1) = nanmean(gpower2(contindsnl,:));
            controlgamma3pow(cll,1) = nanmean(gpower3(contindsnl,:));
            controlgamma1powerr(cll,1) = nanstd(gpower1(contindsnl,:))./sqrt(length(contindsnl));
            controlgamma2powerr(cll,1) = nanstd(gpower2(contindsnl,:))./sqrt(length(contindsnl));
            controlgamma3powerr(cll,1) = nanstd(gpower3(contindsnl,:))./sqrt(length(contindsnl));
            
            contindsl = find(result.gratingInfo.Orientation == -1 & result.light == 1);
            controlresp(2,:) = nanmean(resp(contindsl,:),1);
            controlfr(cll,2) = nanmean(frs(contindsl));
            blscontrolfr(cll,2) = controlfr(cll,2)-bfr(cll);
            controlerr(cll,2) = nanstd(frs(contindsl))./sqrt(length(contindsl));
            controllfpspect(cll,2,:) = nanmean(lfpspect(contindsl,:));
            controllfpspecterr(cll,2,:) = nanstd(lfpspect(contindsl,:))./sqrt(length(contindsl)); 
            controlgamma1pow(cll,2) = nanmean(gpower1(contindsl,:));
            controlgamma2pow(cll,2) = nanmean(gpower2(contindsl,:));
            controlgamma3pow(cll,2) = nanmean(gpower3(contindsl,:));
            controlgamma1powerr(cll,2) = nanstd(gpower1(contindsl,:))./sqrt(length(contindsl));
            controlgamma2powerr(cll,2) = nanstd(gpower2(contindsl,:))./sqrt(length(contindsl));  
            controlgamma3powerr(cll,2) = nanstd(gpower3(contindsl,:))./sqrt(length(contindsl));  
            
            r0l0continds = intersect(contindsnl,stilltrials);
            r0l1continds = intersect(contindsl,stilltrials);
            r1l0continds = intersect(contindsnl,oktrials);
            r1l1continds = intersect(contindsl,oktrials);
            
            if ~isempty(r0l0continds)
                [contS,chf,contSerr] = mtspectrumc(squeeze(lfpresp(r0l0continds,1001:1800))',params);
                r0l0contS(cll,:) = contS(1:150); r0l0contSerr(cll,:,:) = contSerr(:,1:150);
                r0l0contallS(cll,:) = nanmean(allS(r0l0continds,:),1); 
                r0l0contallSerr(cll,:) = nanstd(allS(r0l0continds,:),1,1)./sqrt(length(r0l0continds));
            else
                r0l0contS(cll,:) = nan(1,150); r0l0contSerr(cll,:,:) = nan(2,150);
                r0l0contallS(cll,:) = nan(1,513); r0l0contallSerr(cll,:,:) = nan(1,513);
            end
            if ~isempty(r0l1continds)
                [contS,chf,contSerr] = mtspectrumc(squeeze(lfpresp(r0l1continds,1001:1800))',params);
                r0l1contS(cll,:) = contS(1:150); r0l1contSerr(cll,:,:) = contSerr(:,1:150);
                r0l1contallS(cll,:) = nanmean(allS(r0l1continds,:),1); 
                r0l1contallSerr(cll,:) = nanstd(allS(r0l1continds,:),1,1)./sqrt(length(r0l1continds)); 
            else
                r0l1contS(cll,:) = nan(1,150); r0l1contSerr(cll,:,:) = nan(2,150);
                r0l1contallS(cll,:) = nan(1,513); r0l1contallSerr(cll,:,:) = nan(1,513);
            end
            if ~isempty(r1l0continds)
                [contS,chf,contSerr] = mtspectrumc(squeeze(lfpresp(r1l0continds,1001:1800))',params);
                r1l0contS(cll,:) = contS(1:150); r1l0contSerr(cll,:,:) = contSerr(:,1:150);
                r1l0contallS(cll,:) = nanmean(allS(r1l0continds,:),1); 
                r1l0contallSerr(cll,:) = nanstd(allS(r1l0continds,:),1,1)./sqrt(length(r1l0continds));
            else
                r1l0contS(cll,:) = nan(1,150); r1l0contSerr(cll,:,:) = nan(2,150);
                r1l0contallS(cll,:) = nan(1,513); r1l0contallSerr(cll,:,:) = nan(1,513);
            end
            if ~isempty(r1l1continds)
                [contS,chf,contSerr] = mtspectrumc(squeeze(lfpresp(r1l1continds,1001:1800))',params);
                r1l1contS(cll,:) = contS(1:150); r1l1contSerr(cll,:,:) = contSerr(:,1:150);   
                r1l1contallS(cll,:) = nanmean(allS(r1l1continds,:),1); 
                r1l1contallSerr(cll,:) = nanstd(allS(r1l1continds,:),1,1)./sqrt(length(r1l1continds)); 
            else
                r1l1contS(cll,:) = nan(1,150); r1l1contSerr(cll,:,:) = nan(2,150);
                r1l1contallS(cll,:) = nan(1,513); r1l1contallSerr(cll,:,:) = nan(1,513);
            end
            
            r0controllfpspect(cll,1,:) = nan(1,size(lfpspect,2));
            r0controllfpspect(cll,1,:) = nanmean(lfpspect(r0l0continds,:));
            r0controllfpspect(cll,2,:) = nan(1,size(lfpspect,2));
            r0controllfpspect(cll,2,:) = nanmean(lfpspect(r0l1continds,:));
            r1controllfpspect(cll,1,:) = nan(1,size(lfpspect,2));
            r1controllfpspect(cll,1,:) = nanmean(lfpspect(r1l0continds,:));
            r1controllfpspect(cll,2,:) = nan(1,size(lfpspect,2));
            r1controllfpspect(cll,2,:) = nanmean(lfpspect(r1l1continds,:));
            
            r1controlg1pow(cll,1) = nanmean(gpower1(r1l0continds));
            r0controlg1pow(cll,1) = nanmean(gpower1(r0l0continds));
            r1controlg1powerr(cll,1) = nanstd(gpower1(r1l0continds))./sqrt(length(r1l0continds));
            r0controlg1powerr(cll,1) = nanstd(gpower1(r0l0continds))./sqrt(length(r0l0continds));
            r1controlg2pow(cll,1) = nanmean(gpower2(r1l0continds));
            r0controlg2pow(cll,1) = nanmean(gpower2(r0l0continds));
            r1controlg2powerr(cll,1) = nanmean(gpower2(r1l0continds))./sqrt(length(r1l0continds));
            r0controlg2powerr(cll,1) = nanmean(gpower2(r0l0continds))./sqrt(length(r0l0continds));
            r1controlg3pow(cll,1) = nanmean(gpower3(r1l0continds));
            r0controlg3pow(cll,1) = nanmean(gpower3(r0l0continds));
            r1controlg3powerr(cll,1) = nanmean(gpower3(r1l0continds))./sqrt(length(r1l0continds));
            r0controlg3powerr(cll,1) = nanmean(gpower3(r0l0continds))./sqrt(length(r0l0continds));
            r1controlg1pow(cll,2) = nanmean(gpower1(r1l1continds));
            r0controlg1pow(cll,2) = nanmean(gpower1(r0l1continds));
            r1controlg1powerr(cll,2) = nanstd(gpower1(r1l1continds))./sqrt(length(r1l1continds));
            r0controlg1powerr(cll,2) = nanstd(gpower1(r0l1continds))./sqrt(length(r0l1continds));
            r1controlg2pow(cll,2) = nanmean(gpower2(r1l1continds));
            r0controlg2pow(cll,2) = nanmean(gpower2(r0l1continds));
            r1controlg2powerr(cll,2) = nanmean(gpower2(r1l1continds))./sqrt(length(r1l1continds));
            r0controlg2powerr(cll,2) = nanmean(gpower2(r0l1continds))./sqrt(length(r0l1continds));
            r1controlg3pow(cll,2) = nanmean(gpower3(r1l1continds));
            r0controlg3pow(cll,2) = nanmean(gpower3(r0l1continds));
            r1controlg3powerr(cll,2) = nanmean(gpower3(r1l1continds))./sqrt(length(r1l1continds));
            r0controlg3powerr(cll,2) = nanmean(gpower3(r0l1continds))./sqrt(length(r0l1continds));
            
            spontdiffp(cll) = ranksum(frs(contindsnl),frs(contindsl));
            conddiffp(cll,:,:) = cdiffp;
            
            [nloriprefratio(cll), nldirprefratio(cll), nlprefori, meanoril0(cll), nlosi(cll), meandirl0(cll), nldsi(cll)] = getOSI(squeeze(condfr(cll,1,:))',oris);
            [loriprefratio(cll), ldirprefratio(cll), lprefori, meanoril1(cll), losi(cll), meandirl1(cll), ldsi(cll)] = getOSI(squeeze(condfr(cll,2,:))',oris);

            [preffr(cll,:), prefdir] = max(condfr(cll,1,:),[],3); prefdir = prefdir(1);
            [binprefl0,bta] = binit(squeeze(condresp(cll,1,prefdir,:)),binwidth); binprefl0 = binprefl0.*(1000/binwidth);
            [binprefl1,bta] = binit(squeeze(condresp(cll,2,prefdir,:)),binwidth); binprefl1 = binprefl1.*(1000/binwidth);
            ta = bta-prestim;
            
            l0prefresp(cll,:) = binprefl0;
            l1prefresp(cll,:) = binprefl1;
            l0meanresp(cll,:) = binnednolight;
            l1meanresp(cll,:) = binnedlight;
            respta = ta;
                       
            
            % bin and fit a sine of correct temporal frequency
            binprefl0 = binprefl0-nanmean(bl);  % subtract baseline
            binprefl1 = binprefl1-nanmean(bl);
            stwin = find(ta>700&ta<1500);  % only take part of response after transient
            spsigl0 = binprefl0(stwin); spsigl1 = binprefl1(stwin);
            tx = ta(stwin);        % time axis for afer transient
            tempfreq = unique(result.gratingInfo.tFreq);
            sppl0 = fit_fixedfsin(tx,spsigl0,tempfreq,sr); % the fit parameters ( Amplitude, Phase and Offset)
            sppl1 = fit_fixedfsin(tx,spsigl1,tempfreq,sr);
            f1f0l0(cll) = (sppl0(1))/sppl0(3); % Amplitude = F1, Offset = F0
            f1f0l1(cll) = (sppl1(1))/sppl1(3);
            f1l0(cll) = sppl0(1); f1l1(cll) = sppl1(1); f0l0(cll) = sppl0(3); f0l1(cll) = sppl1(3);

%             figure
            clf;
            subplot(2,2,1)
            ta = bta-prestim;
            plot(ta,binnednolight,'k','linewidth',2);
            hold on
            plot(ta,binnedlight,'color',lcol,'linewidth',2);
            mx = max(binnedlight);
            axis([-prestim,trialdur+poststim,-0.05,mx]);
            line([0,0],[0,mx],'color','k','linewidth',2);
            line([2000,2000],[0,mx],'color','k','linewidth',2);
            line([500,500],[0,mx],'color','b','linewidth',2)
            line([1500,1500],[0,mx],'color','b','linewidth',2);
            legend({'Light OFF','Light ON'})
            xlabel('time [ms]')
            ylabel('firingrate [Hz]')
            title(['cll ' int2str(cll) ' depth: ' int2str(result.depth), 'cll ' printname ])

            
            subplot(2,2,2)
%             if length(oris) == 16
                errorbar(oris,squeeze(blscondfr(cll,2,:)),squeeze(conderr(cll,2,:)),'o-','color',lcol,'markersize',8,'linewidth',2)
                hold on
                errorbar(oris,squeeze(blscondfr(cll,1,:)),squeeze(conderr(cll,1,:)),'ko-','markersize',8,'linewidth',2)
%             else
%                 errorbar(oris,squeeze(blscondfr(cll,2,1:4)),squeeze(conderr(cll,2,1:4)),'o-','color',lcol,'markersize',8,'linewidth',2)
%                 hold on
%                 errorbar(oris,squeeze(blscondfr(cll,1,1:4)),squeeze(conderr(cll,1,1:4)),'ko-','markersize',8,'linewidth',2)
%             end
            xlabel('shown orientation')
            ylabel('Firing rate [Hz]')
            set(gca,'xtick',oris)
            line([-5,320],[0,0],'color','k')
            legend({'Light ON','Light OFF'})
%             title(['OSI: ' num2str(nlosi(cll)) ' light OSI: ' num2str(losi(cll)) ' DSI: ' num2str(nldsi(cll)) ' light DSI: ' num2str(ldsi(cll)) '    ' prefcstring])

            
            subplot(2,2,3)
            errorbar(oris,squeeze(condfr(cll,2,:)),squeeze(conderr(cll,2,:)),'o-','color',lcol,'markersize',8,'linewidth',2)
            hold on
            errorbar(oris,squeeze(condfr(cll,1,:)),squeeze(conderr(cll,1,:)),'ko-','markersize',8,'linewidth',2)
            xlabel('shown orientation')
            ylabel('Firing rate [Hz]')
            set(gca,'xtick',oris)
            line([-5,320],[0,0],'color','k')
            legend({'Light ON','Light OFF'})
            
            subplot(2,2,4)
            plot(interpspike)
            title(['swidth: ' int2str(swidth(cll))]);
            
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
            lfrerr = nanstd(frs(find(result.light)))./sqrt(length(find(result.light)));
            nlfrerr = nanstd(frs(find(~result.light)))./sqrt(length(find(~result.light)));
            runlfrerr = nanstd(frs(intersect(find(result.light),oktrials)))./sqrt(length(intersect(find(result.light),oktrials)));
            runnlfrerr = nanstd(frs(intersect(find(~result.light),oktrials)))./sqrt(length(intersect(find(~result.light),oktrials)));
            norunlfrerr = nanstd(frs(intersect(find(result.light),stilltrials)))./sqrt(length(intersect(find(result.light),stilltrials)));
            norunnlfrerr = nanstd(frs(intersect(find(~result.light),stilltrials)))./sqrt(length(intersect(find(~result.light),stilltrials)));
            
            l1r1 = frs(intersect(find(result.light),oktrials));
            l0r1 = frs(intersect(find(~result.light),oktrials));
            l1r0 = frs(intersect(find(result.light),stilltrials));
            l0r0 = frs(intersect(find(~result.light),stilltrials));
            anovavec = [l0r0;l0r1;l1r0;l1r1];
            g1 = [zeros(length(l0r0),1);zeros(length(l0r1),1);ones(length(l1r0),1);ones(length(l1r1),1)]; %light
            g2 = [zeros(length(l0r0),1);ones(length(l0r1),1);zeros(length(l1r0),1);ones(length(l1r1),1)]; %running
            [p,table,stats] = anovan(anovavec,{g1 g2},'model','full','display','off');
            lp(cll) = p(1); rp(cll) = p(2); rlip(cll) = p(3);
            
            r0omi(cll) = (norunlfr(cll)-norunnlfr(cll))/(norunlfr(cll)+norunnlfr(cll));
            r1omi(cll) = (runlfr(cll)-runnlfr(cll))/(runlfr(cll)+runnlfr(cll));
            l0rmi(cll) = (runnlfr(cll)-norunnlfr(cll))/(runnlfr(cll)+norunnlfr(cll));
            l1rmi(cll) = (runlfr(cll)-norunlfr(cll))/(runlfr(cll)+norunlfr(cll));
            
%             figure
            clf
            subplot(2,2,1)
            imagesc(speed);
            colorbar
            title(['oktrials: ' int2str(length(oktrials)) '/' int2str(size(speed,1))])
            xlabel('time [ms]')
            ylabel('trial number')
            
            subplot(2,2,2)
            errorbar(msta,mean(speed(find(~result.light),:)),nanstd(speed(find(~result.light),:))./sqrt(length(find(~result.light))),'b')
            hold on
            errorbar(msta,mean(speed(find(result.light),:)),nanstd(speed(find(result.light),:))./sqrt(length(find(result.light))),'r')
            xlabel('time [ms]')
            ylabel('average runspeed')
            legend({'light off' 'light on'})
            
            subplot(2,2,3)
            plot(mean(speed(:,respwin),2),frs,'.')
            hold on
            plot(mean(speed(find(result.light),respwin),2),frs(find(result.light)),'r.')
            xlabel('average runspeed of trial')
            ylabel('average firing rate of trial')
            
            subplot(2,2,4)
            barweb([nlfr(cll),lfr(cll);runnlfr(cll),runlfr(cll);norunnlfr(cll),norunlfr(cll)],...
                [nlfrerr,lfrerr;runnlfrerr,runlfrerr;norunnlfrerr,norunlfrerr],...
                [],[{'all'};{'running only'};{'immobile only'}],'firing rates with running',...
                [],'firing rate [Hz]',[],[]);
            
            if printyn
                figSize = [30 21];
                set(gcf,'paperunits','centimeters','papersize',figSize,'paperposition',[0 0 figSize])
                if cll<10, printi = ['0', int2str(cll)]; else printi = int2str(cll); end
                print([runprintpath ,  printi '__' files(fi).name '.pdf'],'-dpdf')
            end
            
            cll = cll + 1;
            disp([files(fi).name '   done'])

        end
    end
    
    save(popfile, '-v7.3');
%     ,'cllname','depth','adiff','swidth','ptr','eslope','bfr','lfr','nlfr','clevels',...
%         'controlfr','controlerr','blscontrolfr', 'nlosi','losi','nldsi','ldsi','nloriprefratio','loriprefratio',...
%         'meanoril0','meanoril1','meandirl0','meandirl1',...
%         'nldirprefratio','ldirprefratio','lcresp','nolcresp',...
%         'nlrmax', 'lrmax', 'nlc50', 'lc50', 'nlr0', 'lr0', 'nrparams','condfr', 'conderr','preffr',...
%         'nlresnorm', 'lresnorm','nlrsq','lrsq','lightmod','animalno','allrl0','allrl1',...
%         'allcmeanl0','allcmeanl1','f1f0l0','f1f0l1','f1l0','f1l1','f0l0','f0l1','xlevels','oris',...
%         'rp','lp','rlip','r0omi','r1omi','l0rmi','l1rmi','nok','nstill','runnlfr','runlfr','norunnlfr','norunlfr',...
%         'r0condfr','r0conderr','r1condfr','r1conderr','tvalues','contrastsig',...
%         'speedpowcc','speedpowp','l0speedpowcc','l0speedpowp','l1speedpowcc','l1speedpowp',...
%         'l0prefresp','l1prefresp','l0meanresp','l1meanresp','respta',...
%         'l0power','l1power','r0power','r1power','l1lfp','l0lfp','r0lfp','r1lfp','bincondcllresp',...
%         'condbandpow','stillcondbandpow','runcondbandpow','condr','condcmean','waveform','clustqual',...
%         'trialfrl0','trialfrl1','trialbl','vismod','condvismod','cllz','cllsc','cllff','recording',...
%         'clleckerrely','cllbinrely','cllmsrely','condbr','nbursts','nspikes','l0isi','l1isi','conddiffp','spontdiffp',...
%         'ta','ntrialsr0','ntrialsr1','lightp','contrastp','clip',...
%         'condmsl','condgpow','condlfpspect','condlfpspecterr','condlfpoffsspect','controllfpspect','controllfpspecterr','fax',...
%         'r0controllfpspect','r1controllfpspect','r0condlfpspect','r1condlfpspect','condisi');  %,'condbandpowresp'
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

prsv = swidthms>=.38; pfsv = swidthms<=.36;
prs = find(prsv); pfs = find(pfsv);

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

% %putative halo expressing for Scnn pop responsive only: 
% phe([17,54,81,83,143,174,175]) = 1;

%putative halo expressing for Scnn pop Scott paper: 
% phe([16,24,71,85,106,111,113,216,260,261,321,368,382,386]) = 1;

% % SCNN with new animals
% phe([16,24,71,85,106,111,113,216,260,261,321,368,382,383]) = 1;

% SCNN eArch
% phe([]) = 1;

% %putative halo expressing for SOM P0 pop: 
% phe([]) = 1;

%putative halo expressing for SOM later pop: 
phe([2,23,36,147,194 ]) = 1;

% putative halo for PV Halo
% phe([12,62,72]) = 1;

phe = logical(phe)';

% % % for SOM later pop, around
% [400,none,300,315,none,300,none,320,336,336,300,336,336,336
% lfpinds = [31, 82, 93, 115, 125, 134, 141, 154, 164, 175, 184];

% for SOM P0 pop - shoot for 320, animal 4,6,7 433 uppermost 5,8,10,11 too deep
% lfpinds = [];

% % for PV Halo pop, around 319 or 336
% lfpinds = [6,13,30,46,59,77];


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
l23rs = l23&prsv&~phe'&ok;
l23fs = l23&pfsv&~phe'&ok;
l4rs = l4&prsv&~phe'&ok;
l4fs = l4&pfsv&~phe'&ok;
l5rs = l5&prsv&~phe'&ok;
l5fs = l5&pfsv&~phe'&ok;
l6rs = l6&prsv&~phe'&ok;


for i = 1:length(depth)
    for l = 1:2
        [r1prefratio(i), r1dirprefratio(i), prefori, r1meanori(i), r1osi(i), r1meandir(i), r1dsi(i)] = getOSI(squeeze(r1condfr(i,1,:))',oris);
        [r0prefratio(i), r0dirprefratio(i), prefori, r0meanori(i), r0osi(i), r0meandir(i), r0dsi(i)] = getOSI(squeeze(r0condfr(i,1,:))',oris);
    end
end


cond = l23rs;
figure
r0mean = squeeze(nanmean(r0osi(cond,1,:),1)); r1mean = squeeze(nanmean(r1osi(cond,1,:),1));
for c = 1:5
    r0stde(c) = squeeze(nanstd(r0osi(cond,1,c),1))./sqrt(length(find(~isnan(r0osi(cond,1,c)))));
    r1stde(c) = squeeze(nanstd(r1osi(cond,1,c),1))./sqrt(length(find(~isnan(r1osi(cond,1,c)))));
end
errorbar(clevels,r0mean,r0stde,'bo-','markerfacecolor','b','linewidth',2)
hold on
errorbar(clevels,r1mean,r1stde,'ro-','markerfacecolor','r','linewidth',2)
legend('non-running','running')
set(gca,'xtick',clevels);
xlabel('contrast levels')
ylabel('OSI')
axis([0,1.1,.25,.5])

figure
errorbar(clevels,squeeze(nanmean(nanmean(r0condfr(cond,1,:,:),3),1)),squeeze(nanstd(nanmean(r0condfr(cond,1,:,:),3),1))./sqrt(sum(cond)),'b')
hold on
errorbar(clevels,squeeze(nanmean(nanmean(r1condfr(cond,1,:,:),3),1)),squeeze(nanstd(nanmean(r1condfr(cond,1,:,:),3),1))./sqrt(sum(cond)),'r')
set(gca,'xtick',clevels);
xlabel('contrast levels')
ylabel('firing rate')
legend('non-running','running')

% individual gamma peaks
for i = 1:length(lfpinds)
    beta = [20,40];
    b1 = find(f>beta(1),1); b2 = find(f>beta(2),1);
    g1 = find(f>gamma(1),1); g2 = find(f>gamma(2),1);
    for s = 1:5
        for l = 1:2
            smoothspec(i,s,l,:) = smooth(squeeze(condS(lfpinds(i),l,s,:)));

            bsig = squeeze(smoothspec(i,s,l,b1:b2));        
            gsig = squeeze(smoothspec(i,s,l,g1:g2));
        
            if isempty(find(diff(bsig)>0)) % there is no clear beta peak
                slgi(i,s,l) = NaN;
            else
                peaks = find(diff(bsig)>0)+1;
                pvs = bsig(peaks);
                bpi = peaks(pvs == max(pvs));
                bpi = bpi+b1-1;
                slgi(i,s,l) = bpi;
            end
        
            if isempty(find(diff(gsig)>0)) % there is no clear beta peak
                shgi(i,s,l) = NaN;
            else
                peaks = find(diff(gsig)>0)+1;
                pvs = gsig(peaks);
                gpi = peaks(pvs == max(pvs));
                gpi = gpi+g1-1;
                shgi(i,s,l) = gpi;
            end
        
            if ~isnan(slgi(i,s,l))
                indpeakgammapow(i,s,l) = smoothspec(i,s,l,slgi(i,s,l));
                dipsig = squeeze(smoothspec(i,s,l,slgi(i,s,l)-17:slgi(i,s,l)-3));
                [xx,ind] = min(abs(diff(dipsig(diff(dipsig)>0))));
                hlp = find(diff(dipsig)>0);
                if ~isempty(ind)
                    dipi(i,s,l) = slgi(i,s,l)-17+hlp(ind)-1;
                    inddipgammapow(i,s,l) = smoothspec(i,s,l,dipi(i,s,l));
                else
                    dipi(i,s,l) = NaN;
                    inddipgammapow(i,s,l) = NaN;
                end
            else
                indpeakgammapow(i,s,l) = NaN;
                inddipgammapow(i,s,l) = NaN;
                dipi(i,s,l) = NaN;
            end
        
            if ~isnan(shgi(i,s,l))
                indhighgammapow(i,s,l) = smoothspec(i,s,l,shgi(i,s,l));
                dipsig = squeeze(smoothspec(i,s,l,shgi(i,s,l)-17:shgi(i,s,l)-3));
                [xx,ind] = min(abs(diff(dipsig(diff(dipsig)>0))));
                hlp = find(diff(dipsig)>0);
                if ~isempty(ind)
                    dipihg(i,s,l) = shgi(i,s,l)-17+hlp(ind)-1;
                    indhighgammadippow(i,s,l) = smoothspec(i,s,l,dipihg(i,s,l));
                else
                    dipihg(i,s,l) = NaN;
                    indhighgammadippow(i,s,l) = NaN;
                end
            else
                indhighgammapow(i,s,l) = NaN;
                indhighgammadippow(i,s,l) = NaN;
                dipihg(i,s,l) = NaN;
            end
        end
    end
    
    for l = 1:2
        if l == 1
            smoothc0spec(i,l,:) = smooth(squeeze(r1l0contS(lfpinds(i),:))); 
        else
            smoothc0spec(i,l,:) = smooth(squeeze(r1l1contS(lfpinds(i),:)));
        end
        bsig = squeeze(smoothc0spec(i,l,b1:b2));
        gsig = squeeze(smoothc0spec(i,l,g1:g2));
        if isempty(find(diff(bsig)>0)) % there is no clear beta peak
            slgic0(i,l) = NaN;
        else
            peaks = find(diff(bsig)>0)+1;
            pvs = bsig(peaks);
            bpi = peaks(pvs == max(pvs));
            bpi = bpi+b1-1;
            slgic0(i,l) = bpi;
        end
    
        if isempty(find(diff(gsig)>0)) % there is no clear beta peak
            shgic0(i,l) = NaN;
        else
            peaks = find(diff(gsig)>0)+1;
            pvs = gsig(peaks);
            gpi = peaks(pvs == max(pvs));
            gpi = gpi+g1-1;
            shgic0(i,l) = gpi;
        end
    
        if ~isnan(slgic0(i,l))
            indc0peakgammapow(i,l) = smoothc0spec(i,l,slgic0(i,l));
        else
            indc0peakgammapow(i,l) = NaN;
        end
        if ~isnan(shgic0(i,l))
            indc0highgammapow(i,l) = smoothc0spec(i,l,shgic0(i,l));
        else
            indc0highgammapow(i,l) = NaN;
        end
        % gamma power at individual peaks
        clear hlp;
        hlp = [indc0peakgammapow(i,l), squeeze(indpeakgammapow(i,:,l))];
        indnc0gp(i,l,:) = hlp./max(hlp);
        indnnc0gp(i,l,:) = hlp;
        clear hlp;
        hlp = [indc0highgammapow(i,l), squeeze(indhighgammapow(i,:,l))];
        indnc0hgp(i,l,:) = hlp./max(hlp);
    end
    
    figure
    semilogy(squeeze(smoothspec(i,1,1,:)),'b')
    hold on
    semilogy(squeeze(smoothspec(i,2,1,:)),'c')
    semilogy(squeeze(smoothspec(i,3,1,:)),'g')
    semilogy(squeeze(smoothspec(i,4,1,:)),'m')
    semilogy(squeeze(smoothspec(i,5,1,:)),'r')
    if ~isnan(slgi(i,1,1))
        plot(slgi(i,1,1),smoothspec(i,1,1,slgi(i,1,1)),'b*')
    end
    if ~isnan(shgi(i,1,1))
        plot(shgi(i,1,1),smoothspec(i,1,1,shgi(i,1,1)),'b*')
    end
    if ~isnan(dipi(i,1,1))
        plot(dipi(i,1,1),smoothspec(i,1,1,dipi(i,1,1)),'bo')
    end
    if ~isnan(dipihg(i,1,1))
        plot(dipihg(i,1,1),smoothspec(i,1,1,dipihg(i,1,1)),'bo')
    end
    if ~isnan(slgi(i,2,1))
        plot(slgi(i,2,1),smoothspec(i,2,1,slgi(i,2,1)),'c*')
    end
    if ~isnan(shgi(i,2,1))
        plot(shgi(i,2,1),smoothspec(i,2,1,shgi(i,2,1)),'c*')
    end
    if ~isnan(dipi(i,2,1))
        plot(dipi(i,2,1),smoothspec(i,2,1,dipi(i,2,1)),'co')
    end
    if ~isnan(dipihg(i,2,1))
        plot(dipihg(i,2,1),smoothspec(i,2,1,dipihg(i,2,1)),'co')
    end
    if ~isnan(slgi(i,3,1))
        plot(slgi(i,3,1),smoothspec(i,3,1,slgi(i,3,1)),'g*')
    end
    if ~isnan(shgi(i,3,1))
        plot(shgi(i,3,1),smoothspec(i,3,1,shgi(i,3,1)),'g*')
    end
    if ~isnan(dipi(i,3,1))
        plot(dipi(i,3,1),smoothspec(i,3,1,dipi(i,3,1)),'go')
    end
    if ~isnan(dipihg(i,3,1))
        plot(dipihg(i,3,1),smoothspec(i,3,1,dipihg(i,3,1)),'go')
    end
    if ~isnan(slgi(i,4,1))
        plot(slgi(i,4,1),smoothspec(i,4,1,slgi(i,4,1)),'m*')
    end
    if ~isnan(shgi(i,4,1))
        plot(shgi(i,4,1),smoothspec(i,4,1,shgi(i,4,1)),'m*')
    end
    if ~isnan(dipi(i,4,1))
        plot(dipi(i,4,1),smoothspec(i,4,1,dipi(i,4,1)),'mo')
    end
    if ~isnan(dipihg(i,4,1))
        plot(dipihg(i,4,1),smoothspec(i,4,1,dipihg(i,4,1)),'mo')
    end
    if ~isnan(slgi(i,5,1))
        plot(slgi(i,5,1),smoothspec(i,5,1,slgi(i,5,1)),'r*')
    end
    if ~isnan(shgi(i,5,1))
        plot(shgi(i,5,1),smoothspec(i,5,1,shgi(i,5,1)),'r*')
    end
    if ~isnan(dipi(i,5,1))
        plot(dipi(i,5,1),smoothspec(i,5,1,dipi(i,5,1)),'ro')
    end
    if ~isnan(dipihg(i,5,1))
        plot(dipihg(i,5,1),smoothspec(i,5,1,dipihg(i,5,1)),'ro')
    end
    
%     P = findpeaksG(chf(1:104),squeeze(r1condallS(lfpinds(i),1,1,1:104)),0,3,3,3,1)
end


for i = 1:size(shgi,1)
    if ~isnan(shgi(i,1,1))
        gcurve(i,:) = [smoothc0spec(i,1,shgi(i,1,1)), smoothspec(i,:,1,shgi(i,1,1))];
    else
        gcurve(i,:) = NaN(1,6);
    end
    normcurve(i,:) = gcurve(i,:)./max(gcurve(i,:));
end
% SOM P0 population in paper
[p,t,s] = kruskalwallis(normcurve)
n = length(find(~isnan(normcurve(:,1))));
figure
errorbar([0,clevels],nanmean(normcurve),nanstd(normcurve)./sqrt(n),'ko-','linewidth',2,'markerfacecolor','k')
set(gca,'xtick',[0,clevels])
xlabel('contrast level')
ylabel('normalized high gamma power')
title(['contrast dependence of high gamma power n = ' int2str(n) '  p: ' num2str(p)]) % for SOM P0 pop


% normalized high gamma
for i = 1:length(lfpinds)
    nihgp(i,:) = indhighgammapow(i,:,1)./max(indhighgammapow(i,:,1));
end

indnc0hgp(6,1,3) = NaN; % only value not NaN for this animal so take it out of the equation
for i = 1:size(indnc0hgp,3)
    indnc0hgperr(i) = nanstd(indnc0hgp(:,1,i),1,1)./sqrt(length(find(~isnan(indnc0hgp(:,1,i)))));
end
clevels = [0, .1, .18, .32, .56, 1];
figure
errorbar(clevels,squeeze(nanmean(indnc0hgp(:,1,:))),indnc0hgperr,'ko-','linewidth',2,'markerfacecolor','k')
set(gca,'xtick',clevels)
xlabel('contrast level')
ylabel('normalized high gamma power')
title('contrast dependence of high gamma power n = 8') % for SOM P0 pop

% orientation tuning with size
for i = 1:length(depth)
    for l = 1:2
        for c = 1:5
            [condoprefratio(i,l,c),conddprefratio(i,l,c),condprefori(i,l,c),condmeanori(i,l,c),condosi(i,l,c),condmeandir(i,l,c),conddsi(i,l,c)] = getOSI(squeeze(condfr(i,l,:,c))',oris);
        end
    end
end
cond = prsv;
anovavec = condosi(cond,:,:); anovavec = anovavec(:);
help = repmat(clevels',1,2*sum(cond))'; 
gs = help(:); %contrast
help = [zeros(1,sum(cond)),ones(1,sum(cond))]'; %light
gl = repmat(help,5,1);
[p,table,stats] = anovan(anovavec,{gs,gl},'model','full');
multcompare(stats)
multcompare(stats,'dimension',2)

%phase locking reloaded
gllol0 = squeeze(nanmean(condr(:,12,1,:,1),4)); % 12 is 56:60 - around the center of the gamma peak
gllol1 = squeeze(nanmean(condr(:,12,2,:,1),4));
glhil0 = squeeze(nanmean(condr(:,12,1,:,5),4));
glhil1 = squeeze(nanmean(condr(:,12,2,:,5),4));

%across cll types
cond = gllol0;
bars = [nanmean(cond(l23rs)),nanmean(cond(l23fs));...
    nanmean(cond(l4rs)),nanmean(cond(l4fs));...
    nanmean(cond(l5rs)),nanmean(cond(l5fs))];
errorbars = [nanstd(cond(l23rs))./sqrt(sum(l23rs)),...
    nanstd(cond(l23fs))./sqrt(sum(l23fs));...
    nanstd(cond(l4rs))./sqrt(sum(l4rs)),...
    nanstd(cond(l4fs))./sqrt(sum(l4fs));...
    nanstd(cond(l5rs))./sqrt(sum(l5rs)),...
    nanstd(cond(l5fs))./sqrt(sum(l5fs))];
figure
barweb(bars, errorbars, [], [{'L23'},{'L4'},{'L5'}], 'Gamma phase locking', [], 'resultant vector length', [], [], [{'RS'},{'FS'}],[],'axis');

%across contrast
condlo = gllol1;
condhi = glhil1;
bars = [nanmean(condlo(l23rs)),nanmean(condhi(l23rs));...
    nanmean(condlo(l4rs)),nanmean(condhi(l4rs));...
    nanmean(condlo(l5rs)),nanmean(condhi(l5rs))];
errorbars = [nanstd(condlo(l23rs))./sqrt(sum(l23rs)),...
    nanstd(condhi(l23rs))./sqrt(sum(l23rs));...
    nanstd(condlo(l4rs))./sqrt(sum(l4rs)),...
    nanstd(condhi(l4rs))./sqrt(sum(l4rs));...
    nanstd(condlo(l5rs))./sqrt(sum(l5rs)),...
    nanstd(condhi(l5rs))./sqrt(sum(l5rs))];
figure
barweb(bars, errorbars, [], [{'L23'},{'L4'},{'L5'}], 'Gamma phase locking', [], 'resultant vector length', [], [], [{'low'},{'high'}],[],'axis');

%across light
cond0 = glhil0;
cond1 = glhil1;
bars = [nanmean(cond0(l23rs)),nanmean(cond1(l23rs));...
    nanmean(cond0(l4rs)),nanmean(cond1(l4rs));...
    nanmean(cond0(l5rs)),nanmean(cond1(l5rs))];
errorbars = [nanstd(cond0(l23rs))./sqrt(sum(l23rs)),...
    nanstd(cond1(l23rs))./sqrt(sum(l23rs));...
    nanstd(cond0(l4rs))./sqrt(sum(l4rs)),...
    nanstd(cond1(l4rs))./sqrt(sum(l4rs));...
    nanstd(cond0(l5rs))./sqrt(sum(l5rs)),...
    nanstd(cond1(l5rs))./sqrt(sum(l5rs))];
figure
barweb(bars, errorbars, [], [{'L23'},{'L4'},{'L5'}], 'Gamma phase locking', [], 'resultant vector length', [], [], [{'L0'},{'L1'}],[],'axis');


%bursts
for i = 1:length(l0isi)
    [n1(i,:),c1(i,:)] = hist(l1isi{i}(l1isi{i}<=50),1:50);
    [n0(i,:),c0(i,:)] = hist(l0isi{i}(l0isi{i}<=50),1:50);
    nn0(i,:) = n0(i,:)./max(n0(i,:));
    nn1(i,:) = n1(i,:)./max(n1(i,:));
    cv0(i) = std(l0isi{i})/mean(l0isi{i});
    cv1(i) = std(l1isi{i})/mean(l1isi{i});
    percl0(i) = length(find(l0isi{i}<10))./length(l0isi{i});
    percl1(i) = length(find(l1isi{i}<10))./length(l1isi{i});
    cbr0(i) = squeeze(mean(mean(condbr(i,1,:,:),3),4));
    cbr1(i) = squeeze(mean(mean(condbr(i,2,:,:),3),4));
end

cond = l5rs;
figure
subplot(2,2,1)
[s,p] = ttest(cbr0(cond),cbr1(cond));
plot(cbr0(cond),cbr1(cond),'o','markersize',4,'markerfacecolor','b')
refline(1,0)
axis square
title(['poisson surprise: p = ' num2str(p)]);
xlabel('poisson surprise light off')
ylabel('poisson surprise light on')

subplot(2,2,2)
[s,p] = ttest(cv0(cond),cv1(cond));
plot(cv0(cond),cv1(cond),'o','markersize',4,'markerfacecolor','b')
a = axis;
axis([min(a),max(a),min(a),max(a)]);
refline(1,0)
axis square
xlabel('CV light off')
ylabel('CV light on')
title(['coefficient of variation: p = ' num2str(p)]);

subplot(2,2,3)
[s,p] = ttest(percl0(cond),percl1(cond));
plot(percl0(cond),percl1(cond),'o','markersize',4,'markerfacecolor','b');
refline(1,0)
axis square
title(['%ISI<10ms: p = ' num2str(p)]);
xlabel('%ISI<10ms light off')
ylabel('%ISI<10ms light on')



disp('')
% % gammapath = 'C:\Users\Julia\work\data\populations\contrast\gamma\';
% % gammapath = 'C:\Users\Julia\work\data\populations\SOM_Halo\contrast\gamma\';
% gammapath = 'C:\Users\Julia\work\data\populations\SOM_Halo\contrast\gamma\running\';
% bandx = 3:5:98;
% for i = 1:size(condbandpow,1)
%     clf;
% %     plot(bandx,squeeze(nanmean(condbandpow(i,1,:,1,:),3)),'b','linewidth',2);
% %     hold on
% %     plot(bandx,squeeze(nanmean(condbandpow(i,1,:,2,:),3)),'c','linewidth',2,'linewidth',2);
% %     plot(bandx,squeeze(nanmean(condbandpow(i,1,:,3,:),3)),'k','linewidth',2);
% %     plot(bandx,squeeze(nanmean(condbandpow(i,1,:,4,:),3)),'m','linewidth',2);
% %     plot(bandx,squeeze(nanmean(condbandpow(i,1,:,5,:),3)),'r','linewidth',2);
%     semilogy(bandx,squeeze(nanmean(runcondbandpow(i,1,:,1,:),3)),'m','linewidth',2)
%     hold on
%     semilogy(bandx,squeeze(nanmean(runcondbandpow(i,1,:,5,:),3)),'r','linewidth',2)
%     semilogy(bandx,squeeze(nanmean(stillcondbandpow(i,1,:,1,:),3)),'c','linewidth',2)
%     semilogy(bandx,squeeze(nanmean(stillcondbandpow(i,1,:,5,:),3)),'b','linewidth',2)
%     xlabel('frequency band')
%     ylabel('power')
% %     legend([{'lowest'},{''},{''},{''},{'highest'}])
%     legend([{'running low'},{'running high'},{'still low'},{'still high'}])
%     title([' i:  ' int2str(i) '  depth: ' num2str(depth(i)) '  width: ' int2str(swidth(i))])
%     figSize = [30 21];
%     set(gcf,'paperunits','centimeters','papersize',figSize,'paperposition',[0 0 figSize])
%     if i<10, printi = ['0', int2str(i)]; else printi = int2str(i); end
%     print([gammapath ,  printi '__' cllname{i} '.pdf'],'-dpdf')
% %     pause
% end

for rec = 1: max(recording)
    recclls = length(find(recording == rec));
    recinds = find(recording == rec);
    
    
    %choristers, solists - correlation with what the rest of pop is doing
    fan = condfr(recinds,1,:,:);
    rfan = reshape(fan,recclls,40);
    tl0 = trialfrl0(recinds,:); tl1 = trialfrl1(recinds,:); tbl = trialbl(recinds,:);
    for i = 1:recclls
        nrfan(i,:) = rfan(i,:)-min(rfan(i,:));
        nrfan(i,:) = nrfan(i,:)./max(nrfan(i,:));
        ntl0(i,:) = tl0(i,:)-min(tl0(i,:));
        ntl0(i,:) = ntl0(i,:)./max(ntl0(i,:));
        ntl1(i,:) = tl1(i,:)-min(tl1(i,:));
        ntl1(i,:) = ntl1(i,:)./max(ntl1(i,:));
        ntbl(i,:) = tbl(i,:)-min(tbl(i,:));
        ntbl(i,:) = ntbl(i,:)./max(ntbl(i,:));
    end
    for i = 1:recclls
        tmpnrfan = nrfan; tmpnrfan(i,:) = nan(1,40); %get rid of this cll for caclulation of pop avg
        [pc(i),pp(i)] = nancorr(rfan(i,:),nanmean(tmpnrfan));
        tmpl0 = ntl0; tmpl0(i,:) = nan(1,210);
        [l0pc(i),l0pp(i)] = nancorr(tl0(i,:),nanmean(tmpl0));
        tmpl1 = ntl1; tmpl1(i,:) = nan(1,210);
        [l1pc(i),l1pp(i)] = nancorr(tl1(i,:),nanmean(tmpl1));
        tmpbl = ntbl; tmpbl(i,:) = nan(1,420);
        [blpc(i),blpp(i)] = nancorr(tbl(i,:),nanmean(tmpbl));
    end
    popcorr(recinds) = pc;
    popp(recinds) = pp;
    popl0corr(recinds) = l0pc;
    popl1corr(recinds) = l1pc;
    popblcorr(recinds) = blpc;
    clear pc; clear pp; clear l0pc; clear l0pp; clear l1pc; clear l1pp; clear blpp; clear blpc;
    
    % lifetime and population kurtosis - large kurtosis means sparser code
    l0fr = condfr(recinds,1,:,:);
    l1fr = condfr(recinds,2,:,:);
    l0fr = reshape(l0fr,recclls,40);
    l1fr = reshape(l1fr,recclls,40);
    ltkurtosisl0(recinds) = kurtosis(l0fr,0,2);
    ltkurtosisl1(recinds) = kurtosis(l1fr,0,2);
    popkurtosisl0(rec,:) = kurtosis(l0fr,0,1);
    popkurtosisl1(rec,:) = kurtosis(l1fr,0,1);
    N = 40;
    sltl0(recinds) = (1 - ((1/N)* ((sum(l0fr,2).^2)./(sum(l0fr.^2,2))) )) / (1-(1/N));
    sltl1(recinds) = (1 - ((1/N)* ((sum(l1fr,2).^2)./(sum(l1fr.^2,2))) )) / (1-(1/N));
    N = recclls;
    spopl0(rec,:) = (1 - ((1/N)* ((sum(l0fr,1).^2)./(sum(l0fr.^2,1))) )) / (1-(1/N));
    spopl1(rec,:) = (1 - ((1/N)* ((sum(l1fr,1).^2)./(sum(l1fr.^2,1))) )) / (1-(1/N));
    
    ii = 1;
    for i = 1:recclls-1
        for j = i+1:recclls
            for l = 1:2
                for o = 1:8
                    for s = 1:5
                        cv = cov(cllz{recinds(i),l,o,s},cllz{recinds(j),l,o,s});
                        cvs(ii,l,o,s) = cv(1,2);
                        cvcount = cov(cllsc{recinds(i),l,o,s},cllsc{recinds(j),l,o,s});
                        if size(cvcount) == 1
                            cvcounts(ii,l,o,s) = NaN;
                        else
                            cvcounts(ii,l,o,s) = cvcount(1,2);
                        end
                        cc = corrcoef(cllz{recinds(i),l,o,s},cllz{recinds(j),l,o,s});
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
    mscrl0{rec} = nanmean(nanmean(cvs(:,1,:,:),3),4); % mean spike count covariance on z scores
    mscrl1{rec} = nanmean(nanmean(cvs(:,2,:,:),3),4);
    mscountrl0{rec} = nanmean(nanmean(cvcounts(:,1,:,:),3),4); %on raw spike counts
    mscountrl1{rec} = nanmean(nanmean(cvcounts(:,2,:,:),3),4);
    mccl0{rec} = nanmean(nanmean(ccs(:,1,:,:),3),4); % correlation coefficients
    mccl1{rec} = nanmean(nanmean(ccs(:,2,:,:),3),4);
    pairs{rec} = recinds(pairinds);
    recs{rec} = recinds;
    clear pairinds ccs cvcounts cvs;
%     %lifetime sparseness
%     N = 34; %1000ms binned 28:61 on bta (500:1500)
%     for cl = 1:recclls
%         for l = 1:2
%             for o = 1:8
%                 for s = 1:5
%                     slt(cl,l,o,s) = (1 - ((1/N)* ((sum(bincondcllresp(recinds(cl),l,o,s,28:61),5).^2)/(sum(bincondcllresp(recinds(cl),l,o,s,28:61).^2,5))) ))...
%                         / (1-(1/N));
%                 end
%             end
%         end
%     end
%     
%     %population sparseness
%     N = size(bincondcllresp,1);
%     for bn = 1:34
%         for l = 1:2
%             for o = 1:8
%                 for s = 1:5
%                     spop(bn,l,o,s) = (1 - ((1/N)* ((sum(bincondcllresp(recinds,l,o,s,bn+27),1).^2)/(sum(bincondcllresp(recinds,l,o,s,bn+27).^2,1))) ))...
%                         / (1-(1/N));
%                 end
%             end
%         end
%     end
%     lifetimesparseness{rec} = slt;
%     popsparseness{rec} = spop;
%     clear slt; clear spop;
    
end


bars = [mean(popl0corr(l23rs)),mean(popl0corr(l23fs));...
    mean(popl0corr(l4rs)),mean(popl0corr(l4fs));...
    mean(popl0corr(l5rs)),mean(popl0corr(l5fs))];
errorbars = [std(popl0corr(l23rs))./sqrt(sum(l23rs)),...
    std(popl0corr(l23fs))./sqrt(sum(l23fs));...
    std(popl0corr(l4rs))./sqrt(sum(l4rs)),...
    std(popl0corr(l4fs))./sqrt(sum(l4fs));...
    std(popl0corr(l5rs))./sqrt(sum(l5rs)),...
    std(popl0corr(l5fs))./sqrt(sum(l5fs))];
figure
barweb(bars, errorbars, [], [{'L23'},{'L4'},{'L5'}], 'correlation with population light OFF', [], 'corr coef', [], [], [{'RS'},{'FS'}],[],'axis');
    
bars = [mean(popl1corr(l23rs)),mean(popl1corr(l23fs));...
    mean(popl1corr(l4rs)),mean(popl1corr(l4fs));...
    mean(popl1corr(l5rs)),mean(popl1corr(l5fs))];
errorbars = [std(popl1corr(l23rs))./sqrt(sum(l23rs)),...
    std(popl1corr(l23fs))./sqrt(sum(l23fs));...
    std(popl1corr(l4rs))./sqrt(sum(l4rs)),...
    std(popl1corr(l4fs))./sqrt(sum(l4fs));...
    std(popl1corr(l5rs))./sqrt(sum(l5rs)),...
    std(popl1corr(l5fs))./sqrt(sum(l5fs))];
figure
barweb(bars, errorbars, [], [{'L23'},{'L4'},{'L5'}], 'correlation with population light ON', [], 'corr coef', [], [], [{'RS'},{'FS'}],[],'axis');
  
% lifetime sparseness
bars = [mean(sltl0(l23rs)),mean(sltl0(l23fs));...
    mean(sltl0(l4rs)),mean(sltl0(l4fs));...
    mean(sltl0(l5rs)),mean(sltl0(l5fs))];
errorbars = [std(sltl0(l23rs))./sqrt(sum(l23rs)),...
    std(sltl0(l23fs))./sqrt(sum(l23fs));...
    std(sltl0(l4rs))./sqrt(sum(l4rs)),...
    std(sltl0(l4fs))./sqrt(sum(l4fs));...
    std(sltl0(l5rs))./sqrt(sum(l5rs)),...
    std(sltl0(l5fs))./sqrt(sum(l5fs))];
figure
barweb(bars, errorbars, [], [{'L23'},{'L4'},{'L5'}], 'Lifetime sparseness light OFF', [], 'corr coef', [], [], [{'RS'},{'FS'}],[],'axis');

% lifetime sparseness
bars = [mean(sltl0(l23rs)),mean(sltl1(l23rs));...
    mean(sltl0(l4rs)),mean(sltl1(l4rs));...
    mean(sltl0(l5rs)),mean(sltl1(l5rs))];
errorbars = [std(sltl0(l23rs))./sqrt(sum(l23rs)),...
    std(sltl1(l23rs))./sqrt(sum(l23rs));...
    std(sltl0(l4rs))./sqrt(sum(l4rs)),...
    std(sltl1(l4rs))./sqrt(sum(l4rs));...
    std(sltl0(l5rs))./sqrt(sum(l5rs)),...
    std(sltl1(l5rs))./sqrt(sum(l5rs))];
figure
barweb(bars, errorbars, [], [{'L23'},{'L4'},{'L5'}], 'lifetime sparseness RS clls', [], 'corr coef', [], [], [{'L0'},{'L1'}],[],'axis');

% lifetime sparseness
bars = [mean(sltl0(l23fs)),mean(sltl1(l23fs));...
    mean(sltl0(l4fs)),mean(sltl1(l4fs));...
    mean(sltl0(l5fs)),mean(sltl1(l5fs))];
errorbars = [std(sltl0(l23fs))./sqrt(sum(l23fs)),...
    std(sltl1(l23fs))./sqrt(sum(l23fs));...
    std(sltl0(l4fs))./sqrt(sum(l4fs)),...
    std(sltl1(l4fs))./sqrt(sum(l4fs));...
    std(sltl0(l5fs))./sqrt(sum(l5fs)),...
    std(sltl1(l5fs))./sqrt(sum(l5fs))];
figure
barweb(bars, errorbars, [], [{'L23'},{'L4'},{'L5'}], 'lifetime sparseness FS clls', [], 'corr coef', [], [], [{'L0'},{'L1'}],[],'axis');



figure
subplot(2,2,1)
semilogy(squeeze(nanmean(nanmean(condbandpow(l4rs,1,:,1,:),1),3)),'linewidth',2)
hold on
semilogy(squeeze(nanmean(nanmean(condbandpow(l4rs,1,:,5,:),1),3)),'c','linewidth',2)
semilogy(squeeze(nanmean(nanmean(condbandpow(l4rs,2,:,1,:),1),3)),'r','linewidth',2)
semilogy(squeeze(nanmean(nanmean(condbandpow(l4rs,2,:,5,:),1),3)),'m','linewidth',2)
legend([{'l0 low'},{'l0 high'},{'l1 low'},{'l1 high'}]);
title('RS L4')

subplot(2,2,2)
semilogy(squeeze(nanmean(nanmean(condbandpow(l23rs,1,:,1,:),1),3)),'linewidth',2)
hold on
semilogy(squeeze(nanmean(nanmean(condbandpow(l23rs,1,:,5,:),1),3)),'c','linewidth',2)
semilogy(squeeze(nanmean(nanmean(condbandpow(l23rs,2,:,1,:),1),3)),'r','linewidth',2)
semilogy(squeeze(nanmean(nanmean(condbandpow(l23rs,2,:,5,:),1),3)),'m','linewidth',2)
legend([{'l0 low'},{'l0 high'},{'l1 low'},{'l1 high'}]);
title('RS L2/3')

subplot(2,2,3)
semilogy(squeeze(nanmean(nanmean(condbandpow(l5rs,1,:,1,:),1),3)),'linewidth',2)
hold on
semilogy(squeeze(nanmean(nanmean(condbandpow(l5rs,1,:,5,:),1),3)),'c','linewidth',2)
semilogy(squeeze(nanmean(nanmean(condbandpow(l5rs,2,:,1,:),1),3)),'r','linewidth',2)
semilogy(squeeze(nanmean(nanmean(condbandpow(l5rs,2,:,5,:),1),3)),'m','linewidth',2)
legend([{'l0 low'},{'l0 high'},{'l1 low'},{'l1 high'}]);
title('RS L5')

figure
subplot(2,2,1)
semilogy(squeeze(nanmean(nanmean(condbandpow(l4fs,1,:,1,:),1),3)),'linewidth',2)
hold on
semilogy(squeeze(nanmean(nanmean(condbandpow(l4fs,1,:,5,:),1),3)),'c','linewidth',2)
semilogy(squeeze(nanmean(nanmean(condbandpow(l4fs,2,:,1,:),1),3)),'r','linewidth',2)
semilogy(squeeze(nanmean(nanmean(condbandpow(l4fs,2,:,5,:),1),3)),'m','linewidth',2)
legend([{'l0 low'},{'l0 high'},{'l1 low'},{'l1 high'}]);
title('FS L4')

subplot(2,2,2)
semilogy(squeeze(nanmean(nanmean(condbandpow(l23fs,1,:,1,:),1),3)),'linewidth',2)
hold on
semilogy(squeeze(nanmean(nanmean(condbandpow(l23fs,1,:,5,:),1),3)),'c','linewidth',2)
semilogy(squeeze(nanmean(nanmean(condbandpow(l23fs,2,:,1,:),1),3)),'r','linewidth',2)
semilogy(squeeze(nanmean(nanmean(condbandpow(l23fs,2,:,5,:),1),3)),'m','linewidth',2)
legend([{'l0 low'},{'l0 high'},{'l1 low'},{'l1 high'}]);
title('FS L2/3')

subplot(2,2,3)
semilogy(squeeze(nanmean(nanmean(condbandpow(l5fs,1,:,1,:),1),3)),'linewidth',2)
hold on
semilogy(squeeze(nanmean(nanmean(condbandpow(l5fs,1,:,5,:),1),3)),'c','linewidth',2)
semilogy(squeeze(nanmean(nanmean(condbandpow(l5fs,2,:,1,:),1),3)),'r','linewidth',2)
semilogy(squeeze(nanmean(nanmean(condbandpow(l5fs,2,:,5,:),1),3)),'m','linewidth',2)
legend([{'l0 low'},{'l0 high'},{'l1 low'},{'l1 high'}]);
title('FS L5')

% % badfits = union(find((nlresnorm./(nlfr.^2))>5), find((lresnorm./(lfr.^2))>5)); % normalize error by firing rate otherwise throws out many relatively good high fr fits
% badfits = find(nlresnorm>1 | lresnorm>1);
% nlc50(badfits) = NaN; lc50(badfits) = NaN;
% nlrmax(badfits) = NaN; lrmax(badfits) = NaN;
% nlr0(badfits) = NaN; lr0(badfits) = NaN;
for an = 1: max(animalno)
    anclls = length(find(animalno == an));
    aninds = find(animalno == an);
    
    fan = condfr(aninds,1,:,:);
    rfan = reshape(fan,anclls,40);
    tl0 = trialfrl0(aninds,:); tl1 = trialfrl1(aninds,:); tbl = trialbl(aninds,:);
    for i = 1:anclls
        nrfan(i,:) = rfan(i,:)-min(rfan(i,:));
        nrfan(i,:) = nrfan(i,:)./max(nrfan(i,:));
        ntl0(i,:) = tl0(i,:)-min(tl0(i,:));
        ntl0(i,:) = ntl0(i,:)./max(ntl0(i,:));
        ntl1(i,:) = tl1(i,:)-min(tl1(i,:));
        ntl1(i,:) = ntl1(i,:)./max(ntl1(i,:));
        ntbl(i,:) = tbl(i,:)-min(tbl(i,:));
        ntbl(i,:) = ntbl(i,:)./max(ntbl(i,:));
    end
    for i = 1:anclls
        tmpnrfan = nrfan; tmpnrfan(i,:) = nan(1,40); %get rid of this cll for caclulation of pop avg
        [pc(i),pp(i)] = nancorr(rfan(i,:),nanmean(tmpnrfan));
        tmpl0 = ntl0; tmpl0(i,:) = nan(1,210);
        [l0pc(i),l0pp(i)] = nancorr(tl0(i,:),nanmean(tmpl0));
        tmpl1 = ntl1; tmpl1(i,:) = nan(1,210);
        [l1pc(i),l1pp(i)] = nancorr(tl1(i,:),nanmean(tmpl1));
        tmpbl = ntbl; tmpbl(i,:) = nan(1,420);
        [blpc(i),blpp(i)] = nancorr(tbl(i,:),nanmean(tmpbl));
    end
    popcorr(aninds) = pc;
    popp(aninds) = pp;
    popl0corr(aninds) = l0pc;
    popl1corr(aninds) = l1pc;
    popblcorr(aninds) = blpc;
    clear pc; clear pp; clear l0pc; clear l0pp; clear l1pc; clear l1pp; clear blpp; clear blpc;
end


normwin = find(respta>500&respta<1500);
for i = 1:size(condfr,1)
    pc(i) = find(nanmean(condfr(i,1,:,:),3) == max(nanmean(condfr(i,1,:,:))),1,'last'); %pref contrast at mean ori
    prefori(i) = find(condfr(i,1,:,pc(i)) == max(condfr(i,1,:,pc(i))),1);       %pref ori at pref contrast
    nonprefori(i) = find(condfr(i,1,:,pc(i)) == min(condfr(i,1,:,pc(i))),1);
    preffrl0(i) = condfr(i,1,prefori(i),pc(i)); %absolute pref fr
    preffrl1(i) = condfr(i,2,prefori(i),pc(i)); 
    preforirangel0(i) = (max(condfr(i,1,:,pc(i)))-min(condfr(i,1,:,pc(i)))); % fr range over oris at pref contrast
    preforirangel1(i) = (max(condfr(i,2,:,pc(i)))-min(condfr(i,2,:,pc(i))));
    prefcontrangel0(i) = (max(condfr(i,1,prefori(i),:))-min(condfr(i,1,prefori(i),:))); %fr range over contrasts at pref ori
    prefcontrangel1(i) = (max(condfr(i,2,prefori(i),:))-min(condfr(i,2,prefori(i),:)));
    meanorirangel0(i) = (max(nanmean(condfr(i,1,:,:),4))-min(nanmean(condfr(i,1,:,:),4))); % fr range over oris at contrast mean
    meanorirangel1(i) = (max(nanmean(condfr(i,2,:,:),4))-min(nanmean(condfr(i,2,:,:),4)));
    meancontrangel0(i) = (max(nanmean(condfr(i,1,:,:),3))-min(nanmean(condfr(i,1,:,:),3))); % fr range over contrast at mean ori
    meancontrangel1(i) = (max(nanmean(condfr(i,2,:,:),3))-min(nanmean(condfr(i,2,:,:),3)));
    normpreforirangel0(i) = (max(condfr(i,1,:,pc(i)))-min(condfr(i,1,:,pc(i))))./max(condfr(i,1,:,pc(i))); % normalized fr range over oris
    normpreforirangel1(i) = (max(condfr(i,2,:,pc(i)))-min(condfr(i,2,:,pc(i))))./max(condfr(i,2,:,pc(i)));
    normprefcontrangel0(i) = (max(condfr(i,1,prefori(i),:))-min(condfr(i,1,prefori(i),:)))./max(condfr(i,1,prefori(i),:)); %normalized fr range over contrasts
    normprefcontrangel1(i) = (max(condfr(i,2,prefori(i),:))-min(condfr(i,2,prefori(i),:)))./max(condfr(i,2,prefori(i),:));
    normmeanorirangel0(i) = (max(nanmean(condfr(i,1,:,:),4))-min(nanmean(condfr(i,1,:,:),4)))./max(nanmean(condfr(i,1,:,:),4));
    normmeanorirangel1(i) = (max(nanmean(condfr(i,2,:,:),4))-min(nanmean(condfr(i,2,:,:),4)))./max(nanmean(condfr(i,2,:,:),4));
    normmeancontrangel0(i) = (max(nanmean(condfr(i,1,:,:),3))-min(nanmean(condfr(i,1,:,:),3)))./max(nanmean(condfr(i,1,:,:),3));
    normmeancontrangel1(i) = (max(nanmean(condfr(i,2,:,:),3))-min(nanmean(condfr(i,2,:,:),3)))./max(nanmean(condfr(i,2,:,:),3));
    condomi(i,:,:) = (condfr(i,2,:,:)-condfr(i,1,:,:))./(condfr(i,2,:,:)+condfr(i,1,:,:)); % omi for each condition 
    controlomi(i) = (controlfr(i,2)-controlfr(i,1))./(controlfr(i,2)+controlfr(i,1));       % omi for control condition
    
    condspikediff(i,:,:) = condfr(i,2,:,:)-condfr(i,1,:,:); % difference in fr for each condition
    help = condfr(i,1,:,:); condspikesl0(i,:) = help(:); condspikesl0wc(i,:) = [condspikesl0(i,:), controlfr(i,1)]; % all 40 or 41 conditions
    help = condfr(i,2,:,:); condspikesl1(i,:) = help(:); condspikesl1wc(i,:) = [condspikesl1(i,:), controlfr(i,2)];
    normcondspikesl0(i,:) = condspikesl0(i,:)./max(condspikesl0(i,:)); % normalized frs over all conditions
    normcondspikesl1(i,:) = condspikesl1(i,:)./max(condspikesl1(i,:));
    
    cafrl0(i,:) = squeeze(nanmean(condfr(i,1,:,:),3)); % fr for contrast conditions - average of oris
    cafrl1(i,:) = squeeze(nanmean(condfr(i,2,:,:),3));
%     carespl0(i,:,:) = squeeze(nanmean(bincondcllresp(i,1,:,:,:),3)); % response to each contrast conditions
%     carespl1(i,:,:) = squeeze(nanmean(bincondcllresp(i,2,:,:,:),3));
    [cas,cai] = sort(cafrl0(i,:));
    scafrl0(i,:) = cafrl0(i,cai); scafrl1(i,:) = cafrl1(i,cai); % each condition sorted by fr
%     scarespl0(i,:,:) = carespl0(i,cai,:); scarespl1(i,:,:) = carespl1(i,cai,:); % each response sorted
%     for j = 1:5 %and normalized
%         nscarespl0(i,j,:) = scarespl0(i,j,:)./nanmean(scarespl0(i,j,normwin),3);
%         nscarespl1(i,j,:) = scarespl1(i,j,:)./nanmean(scarespl0(i,j,normwin),3);
%         ntomscarespl0(i,j,:) = scarespl0(i,j,:)./nanmean(scarespl0(i,5,normwin),3);
%         ntomscarespl1(i,j,:) = scarespl1(i,j,:)./nanmean(scarespl0(i,5,normwin),3);
%     end
    scaomi(i,:) = (scafrl1(i,:)-scafrl0(i,:))./(scafrl1(i,:)+scafrl0(i,:)); % omi per contrast condition sorted
    sncafrl0(i,:) = scafrl0(i,:)./max(scafrl0(i,:)); % normalized sorted contrast fr
    sncafrl1(i,:) = scafrl1(i,:)./max(scafrl0(i,:));
    cafrangel0(i) = (max(cafrl0(i,:))-min(cafrl0(i,:)))/max(cafrl0(i,:)); % normalized range of fr over contrasts
    cafrangel1(i) = (max(cafrl1(i,:))-min(cafrl1(i,:)))/max(cafrl1(i,:));
    cafrangennl0(i) = (max(cafrl0(i,:))-min(cafrl0(i,:)));                  % non normalized range over contrast response
    cafrangennl1(i) = (max(cafrl1(i,:))-min(cafrl1(i,:)));
    caparams(i,:) = polyfit(cafrl0(i,:),cafrl1(i,:),1); % fit to slope of scatter l0 vs l1
    ncafrl0(i,:) = cafrl0(i,:)/max(cafrl0(i,:));
    ncafrl1(i,:) = cafrl1(i,:)/max(cafrl0(i,:));
    siminl0(i) = max(cafrl0(i,:))/nanmean(cafrl0(i,:),2);  % simons index
    siminl1(i) = max(cafrl1(i,:))/nanmean(cafrl1(i,:),2);

    long = [condspikesl0(i,:),condspikesl1(i,:)]; long = long(randperm(length(long)));
    shufl0(i,:) = long(1:40); shufl1(i,:) = long(41:80);
    [shufs,shufi] = sort(shufl0(i,:));
    nshufl0(i,:) = shufl0(i,shufi)./max(shufl0(i,:));  %randomly shuffled normalized sorted fr per cll
    nshufl1(i,:) = shufl1(i,shufi)./max(shufl0(i,:));
    
    range(i) = max(condspikesl0(i,:))-min(condspikesl0(i,:));
    rangel1(i) = max(condspikesl1(i,:))-min(condspikesl1(i,:));
    normrange(i) = (max(condspikesl0wc(i,:))-min(condspikesl0wc(i,:)))./max(condspikesl0wc(i,:));
    normrangel1(i) = (max(condspikesl1wc(i,:))-min(condspikesl1wc(i,:)))./max(condspikesl1wc(i,:));
    mincondl0(i) = min(condspikesl0(i,:)); mincondl1(i) = min(condspikesl1(i,:));
    maxcondl0(i) = max(condspikesl0(i,:)); maxcondl1(i) = max(condspikesl1(i,:));
    [s,in] = sort(condspikesl0(i,:));
    sortedl0(i,:) = condspikesl0(i,in)./maxcondl0(i);  %sorted 40 conditions per cll
    sortedl1(i,:) = condspikesl1(i,in)./maxcondl0(i);
    ncontrolfr(i,:) = controlfr(i,:)./maxcondl0(i);
    sortedparamsl0(i,:) = polyfit(1:40,sortedl0(i,:),1);
    sortedparamsl1(i,:) = polyfit(1:40,sortedl1(i,:),1);
    
%     l0ff(i,:) = reshape(cllff(i,1,:,:),1,40);
%     l1ff(i,:) = reshape(cllff(i,2,:,:),1,40);
%     sl0ff(i,:) = l0ff(i,in);
%     sl1ff(i,:) = l1ff(i,in);
%     prefffl0(i) = cllff(i,1,prefori(i),pc(i));
%     prefffl1(i) = cllff(i,2,prefori(i),pc(i));

    [ss,iin] = sort(condspikesl1(i,:));
    sortedlightl0(i,:) = condspikesl0(i,iin)./maxcondl1(i); %sorted to the light condition
    sortedlightl1(i,:) = condspikesl1(i,iin)./maxcondl1(i);
    [r(i),p(i)] = nancorr(condspikesl0(i,:),condspikesl1(i,:)-condspikesl0(i,:));    %spikes added related to l0 fr?
    fitparams(i,:) = polyfit(condspikesl0(i,:),condspikesl1(i,:)-condspikesl0(i,:),1);
    normparams(i,:) = polyfit(normcondspikesl0(i,:),normcondspikesl1(i,:)-normcondspikesl0(i,:),1);
    zerocross(i) = -fitparams(i,2)./fitparams(i,1);
    if zerocross(i)>min(condspikesl0(i,:)) & zerocross(i)<max(condspikesl0(i,:)) %how many switch from fascilitation to suppression along their preferredness axes
        switches(i) = 1; else switches(i) = 0; end
    
    avgcrfs(i,:,:) = [controlfr(i,:)', squeeze(nanmean(condfr(i,:,:,:),3))];
    navgcrfs(i,:,:) = avgcrfs(i,:,:)./max(avgcrfs(i,1,:),[],3);
    
    meancrf(i,1,:) = [controlfr(i,1)',cafrl0(i,:)];
    meancrf(i,2,:) = [controlfr(i,2)',cafrl1(i,:)];
    meancrferr(i,1,:) = [controlerr(i,1)',squeeze(nanmean(conderr(i,1,:,:),3))'];
    meancrferr(i,2,:) = [controlerr(i,2)',squeeze(nanmean(conderr(i,2,:,:),3))'];
    nmeancrf(i,:,:) = meancrf(i,:,:)./max(meancrf(i,1,:),[],3);
    
    prefcrfs(i,:,:) = [controlfr(i,:)',squeeze(condfr(i,:,prefori(i),:))]; % preferred ori crf with baseline
    preftunecurve(i,:,:) = [controlfr(i,:)',squeeze(condfr(i,:,:,pc(i)))]; % preferred contrast tuning curve
    nprefcrfs(i,:,:) = prefcrfs(i,:,:)./max(prefcrfs(i,1,:),[],3);
    npreftunecurve(i,:,:) = preftunecurve(i,:,:)./max(preftunecurve(i,1,:),[],3);
    [ss,iin] = sort(prefcrfs(i,1,:));
    rankcrfl0(i,:) = prefcrfs(i,1,iin); nrankcrfl0(i,:) = rankcrfl0(i,:)./rankcrfl0(i,6);
    rankcrfl1(i,:) = prefcrfs(i,2,iin); nrankcrfl1(i,:) = rankcrfl1(i,:)./rankcrfl0(i,6);
    
    nonprefcrfs(i,:,:) = [controlfr(i,:)',squeeze(condfr(i,:,nonprefori(i),:))]; % preferred ori crf with baseline
    nnonprefcrfs(i,:,:) = nonprefcrfs(i,:,:)./max(nonprefcrfs(i,1,:),[],3);
    [ss,iin] = sort(nonprefcrfs(i,1,:));
    ranknonprefcrfl0(i,:) = nonprefcrfs(i,1,iin); nranknonprefcrfl0(i,:) = ranknonprefcrfl0(i,:)./ranknonprefcrfl0(i,6);
    ranknonprefcrfl1(i,:) = nonprefcrfs(i,2,iin); nranknonprefcrfl1(i,:) = ranknonprefcrfl1(i,:)./ranknonprefcrfl0(i,6);
    
    scatterparams(i,:) = polyfit(condspikesl0(i,:),condspikesl1(i,:),1);
    unitycross(i) = scatterparams(i,2)/(1-scatterparams(i,1));
    if unitycross(i)>min(condspikesl0(i,:)) & unitycross(i)<max(condspikesl0(i,:))
        unityswitch(i) = 1; else unityswitch(i) = 0; end   %how many cross unity when fitting their l0 vs l1 for all conditions
    
    %normalize average PSTHs to average of no light response period
    l1meanrespn(i,:) = l1meanresp(i,:)./nanmean(l0meanresp(i,normwin),2);
    l0meanrespn(i,:) = l0meanresp(i,:)./nanmean(l0meanresp(i,normwin),2);
    l1prefrespn(i,:) = l1prefresp(i,:)./nanmean(l0prefresp(i,normwin),2);
    l0prefrespn(i,:) = l0prefresp(i,:)./nanmean(l0prefresp(i,normwin),2);
    
    
    % running  
    crespl0r0(i,:) = squeeze(nanmean(r0condfr(i,1,:,:),3));
    crespl1r0(i,:) = squeeze(nanmean(r0condfr(i,2,:,:),3));
    crespl0r1(i,:) = squeeze(nanmean(r1condfr(i,1,:,:),3));
    crespl1r1(i,:) = squeeze(nanmean(r1condfr(i,2,:,:),3));
    ncrespl0r0(i,:)= crespl0r0(i,:)./max(crespl0r1(i,:));
    ncrespl0r1(i,:)= crespl0r1(i,:)./max(crespl0r1(i,:));
    ncrespl1r0(i,:)= crespl1r0(i,:)./max(crespl0r1(i,:));
    ncrespl1r1(i,:)= crespl1r1(i,:)./max(crespl0r1(i,:));
    [ss,si] = sort(ncrespl0r1(i,:));
    sncrespl0r1(i,:) = ncrespl0r1(i,si);
    sncrespl0r0(i,:) = ncrespl0r0(i,si);
    sncrespl1r1(i,:) = ncrespl1r1(i,si);
    sncrespl1r0(i,:) = ncrespl1r0(i,si);    
    
end
meancrfs = cat(3,controlfr,squeeze(mean(condfr,3)));
prefomi = ((preffrl1-preffrl0)./(preffrl1+preffrl0));

%determine which go up
for i = 1:length(nlfr)
    if controlfr(i,1)>nlfr(i)
        up(i) = -1;
    else
        up(i) = 1;
    end
end
cond = l23rs&vismod&up==1;

% average across cond
figure
errorbar([0,clevels],squeeze(nanmean(nmeancrf(cond,1,:))),nanstd(nmeancrf(cond,1,:))./sqrt(sum(cond)),'ko-','linewidth',2,'markerfacecolor','k')

% every single crf
a = find(cond);
for i = 1:length(a)
figure
errorbar([0,clevels],squeeze(meancrf(a(i),1,:)),squeeze(meancrferr(a(i),1,:)))
end


%running vs no rrunning crfs
cond = l23rs;
figure
errorbar(xlevels(2:6),nanmean(ncrespl0r0(cond,:)),nanstd(ncrespl0r0(cond,:))./sqrt(sum(cond)),'b');
hold on
errorbar(xlevels(2:6),nanmean(ncrespl1r0(cond,:)),nanstd(ncrespl1r0(cond,:))./sqrt(sum(cond)),'r');
errorbar(xlevels(2:6),nanmean(ncrespl1r1(cond,:)),nanstd(ncrespl1r1(cond,:))./sqrt(sum(cond)),'m');
errorbar(xlevels(2:6),nanmean(ncrespl0r1(cond,:)),nanstd(ncrespl0r1(cond,:))./sqrt(sum(cond)),'c');

cond = l5rs;
figure
plot(nanmean(ncrespl0r0(cond,:)),nanmean(ncrespl1r0(cond,:)),'bo');
hold on
plot(nanmean(ncrespl0r1(cond,:)),nanmean(ncrespl1r1(cond,:)),'ro');

%fr vs fr plots
cond = l5rs;
figure
plot(nanmean(sortedl0(cond,:)),nanmean(sortedl1(cond,:)),'go','markersize',5,'markerfacecolor','g');
axis([0,1.3,0,1.3])
refline(1,0)
lsline
hold on
plot(nanmean(ncontrolfr(cond,1)),nanmean(ncontrolfr(cond,2)),'co','markersize',5,'markerfacecolor','c')
errorbar(nanmean(sortedl0(cond,:)),nanmean(sortedl1(cond,:)),nanstd(sortedl1(cond,:))./sqrt(sum(cond)),'g.')
herrorbar(nanmean(sortedl0(cond,:)),nanmean(sortedl1(cond,:)),nanstd(sortedl0(cond,:))./sqrt(sum(cond)),'g.')
errorbar(nanmean(ncontrolfr(cond,1)),nanmean(ncontrolfr(cond,2)),nanstd(ncontrolfr(cond,2))./sqrt(sum(cond)),'c.')
herrorbar(nanmean(ncontrolfr(cond,1)),nanmean(ncontrolfr(cond,2)),nanstd(ncontrolfr(cond,1))./sqrt(sum(cond)),'c.')
axis([0,1.3,0,1.3])
axis square
xlabel('normalized firing rate light OFF')
ylabel('normalized firing rate light ON')

cond = l5rs;
figure
plot(nanmean(nrankcrfl0(cond,:)),nanmean(nrankcrfl1(cond,:)),'ko','markersize',5,'markerfacecolor','k');
axis([0,1.3,0,1.3])
refline(1,0)
lsline
hold on
errorbar(nanmean(nrankcrfl0(cond,:)),nanmean(nrankcrfl1(cond,:)),nanstd(nrankcrfl1(cond,:))./sqrt(sum(cond)),'k.')
herrorbar(nanmean(nrankcrfl0(cond,:)),nanmean(nrankcrfl1(cond,:)),nanstd(nrankcrfl0(cond,:))./sqrt(sum(cond)),'k.')
axis([0,1.3,0,1.3])
axis square
xlabel('normalized firing rate light OFF')
ylabel('normalized firing rate light ON')

cond = l5rs;
figure
plot(nanmean(nranknonprefcrfl0(cond,:)),nanmean(nranknonprefcrfl1(cond,:)),'ko','markersize',5,'markerfacecolor','k');
axis([0,1.3,0,1.3])
refline(1,0)
lsline
hold on
errorbar(nanmean(nranknonprefcrfl0(cond,:)),nanmean(nranknonprefcrfl1(cond,:)),nanstd(nranknonprefcrfl1(cond,:))./sqrt(sum(cond)),'k.')
herrorbar(nanmean(nranknonprefcrfl0(cond,:)),nanmean(nranknonprefcrfl1(cond,:)),nanstd(nranknonprefcrfl0(cond,:))./sqrt(sum(cond)),'k.')
axis([0,1.3,0,1.3])
axis square
xlabel('normalized firing rate light OFF')
ylabel('normalized firing rate light ON')

cond = l5fs;
errorbar(mean(rankcrfl1(cond,:)-rankcrfl0(cond,:)),std(rankcrfl1(cond,:)-rankcrfl0(cond,:))./sqrt(sum(cond)),'.')
hold on
errorbar(mean(ranknonprefcrfl1(cond,:)-ranknonprefcrfl0(cond,:)),std(ranknonprefcrfl1(cond,:)-ranknonprefcrfl0(cond,:))./sqrt(sum(cond)),'.')
line([0,7],[0,0])
xlabel('rank according to contrast level')
ylabel('average number of spikes added')
legend([{'best orientation'},{'worst orientation'}],'location','ne')
title(['Layer 5 FS clls n = ' int2str(sum(cond))]);

for i = 1:size(condfr,1)
    l0oricurve = squeeze(condfr(i,1,:,pc(i))); %-controlfr(i,1);
    l1oricurve = squeeze(condfr(i,2,:,pc(i))); %-controlfr(i,2);
    prefol0(i) = find(l0oricurve == max(l0oricurve),1);
    prefol1(i) = find(l1oricurve == max(l1oricurve),1);
    
    shiftcurvel0(i,:) = circshift(l0oricurve,-(prefol0(i)-1));
    shiftcurvel1(i,:) = circshift(l1oricurve,-(prefol0(i)-1));
    shiftcurvectrll0(i) = controlfr(i,1)./shiftcurvel0(i,1);
    shiftcurvectrll1(i) = controlfr(i,2)./shiftcurvel0(i,1);
    shiftcurvel1(i,:) = shiftcurvel1(i,:)./shiftcurvel0(i,1);
    shiftcurvel0(i,:) = shiftcurvel0(i,:)./shiftcurvel0(i,1);
    scl0center(i,:) = circshift(shiftcurvel0(i,:),3,2);
    scl1center(i,:) = circshift(shiftcurvel1(i,:),3,2);
    
    for c = 1:5
        tmpcrvl0 = squeeze(condfr(i,1,:,c)); %-controlfr(i,1);
        tmpcrvl1 = squeeze(condfr(i,2,:,c)); %-controlfr(i,2);
        prefoclsl0(i,c) = find(tmpcrvl0 == max(tmpcrvl0),1);
        prefoclsl1(i,c) = find(tmpcrvl1 == max(tmpcrvl1),1);
        shiftcurveclsl0(i,c,:) = circshift(tmpcrvl0,-(prefoclsl0(i,c)-1));
        shiftcurveclsl1(i,c,:) = circshift(tmpcrvl1,-(prefoclsl0(i,c)-1));
        shiftcurveclsl1(i,c,:) = shiftcurveclsl1(i,c,:)./shiftcurveclsl0(i,c,1);
        shiftcurveclsl0(i,c,:) = shiftcurveclsl0(i,c,:)./shiftcurveclsl0(i,c,1);
    end   
    
    normmeancrfs(i,:,:) = meancrfs(i,:,:);%-repmat(meancrfs(i,1,1),1,2,6);
    normmeancrfsctrl(i,:) = controlfr(i,:)./normmeancrfs(i,1,pc(i)+1);
    normmeancrfs(i,2,:) = normmeancrfs(i,2,:)./normmeancrfs(i,1,pc(i)+1);
    normmeancrfs(i,1,:) = normmeancrfs(i,1,:)./normmeancrfs(i,1,pc(i)+1);
    normprefcrfs(i,:,:) = prefcrfs(i,:,:);%-repmat(prefcrfs(i,1,1),1,2,6);
    normprefcrfsctrl(i,:) = controlfr(i,:)./normprefcrfs(i,1,pc(i)+1);
    normprefcrfs(i,2,:) = normprefcrfs(i,2,:)./normprefcrfs(i,1,pc(i)+1);
    normprefcrfs(i,1,:) = normprefcrfs(i,1,:)./normprefcrfs(i,1,pc(i)+1);
    
    %both normalized
    bnmeancrfsctrl(i,1,:) = controlfr(i,1)./meancrfs(i,1,pc(i)+1);
    bnmeancrfsctrl(i,2,:) = controlfr(i,2)./meancrfs(i,2,pc(i)+1);
    bnmeancrfs(i,1,:) = meancrfs(i,1,:)./meancrfs(i,1,pc(i)+1);
    bnmeancrfs(i,2,:) = meancrfs(i,2,:)./meancrfs(i,2,pc(i)+1);
    bnprefcrfsctrl(i,1,:) = controlfr(i,1)./prefcrfs(i,1,pc(i)+1);
    bnprefcrfsctrl(i,2,:) = controlfr(i,2)./prefcrfs(i,2,pc(i)+1);
    bnprefcrfs(i,1,:) = prefcrfs(i,1,:)./prefcrfs(i,1,pc(i)+1);
    bnprefcrfs(i,2,:) = prefcrfs(i,2,:)./prefcrfs(i,2,pc(i)+1);
    
end


% oris = [0,45,90,135,180,225,270,315];
% controlerr = ones(size(conderr,1),2);
% for cll = find(l5&prsv')
%     figure
%     subplot(2,2,1)
%     errorbar(oris,squeeze(condfr(cll,2,:,pc(cll))),squeeze(conderr(cll,2,:,pc(cll))),'ro-','markersize',8,'linewidth',2)
%     hold on
%     errorbar(oris,squeeze(condfr(cll,1,:,pc(cll))),squeeze(conderr(cll,1,:,pc(cll))),'ko-','markersize',8,'linewidth',2)
%     xlabel('shown orientation')
%     ylabel('Firing rate [Hz]')
%     set(gca,'xtick',oris)
%     line([-5,320],[0,0],'color','k')
%     legend({'Light ON','Light OFF'})
%     ax = axis;
%     axis([-10,320,ax(3),ax(4)])
%     title(['preferred contrast orientation tuning cll: ' int2str(cll)])
% 
%     subplot(2,2,2)
%     errorbar(oris,squeeze(mean(condfr(cll,2,:,:),4)),squeeze(mean(conderr(cll,2,:,:),4)),'ro-','markersize',8,'linewidth',2)
%     hold on
%     errorbar(oris,squeeze(mean(condfr(cll,1,:,:),4)),squeeze(mean(conderr(cll,1,:,:),4)),'ko-','markersize',8,'linewidth',2)
%     xlabel('shown orientation')
%     ylabel('Firing rate [Hz]')
%     set(gca,'xtick',oris)
%     line([-5,320],[0,0],'color','k')
%     legend({'Light ON','Light OFF'})
%     ax = axis;
%     axis([-10,320,ax(3),ax(4)])
%     title(['mean orientation tuning depth: ' int2str(depth(cll))])
%     
%     subplot(2,2,3)
%     errorbar(xlevels,[controlfr(cll,2);squeeze(condfr(cll,2,prefori(cll),:))],[controlerr(cll,2);squeeze(conderr(cll,2,prefori(cll),:))],'ro-','markersize',8,'linewidth',2);
%     hold on
%     errorbar(xlevels,[controlfr(cll,1);squeeze(condfr(cll,1,prefori(cll),:))],[controlerr(cll,1);squeeze(conderr(cll,1,prefori(cll),:))],'ko-','markersize',8,'linewidth',2);
%     xlabel('shown contrast')
%     ylabel('Firing rate [Hz]')
%     legend({'Light ON','Light OFF'})
%     set(gca,'xtick',xlevels); %[0, clevels])
%     ax = axis;
%     axis([-.1,1.1,ax(3),ax(4)])
%     line([-0.1,1.1],[0,0],'color','k')
%     title(['preferred orientation contrast response ' cllname{cll}]);
%     
%     subplot(2,2,4)
%     errorbar(xlevels,lcresp(cll,:),[controlerr(cll,2), squeeze(mean(conderr(cll,2,:,:),3))'],'o-','color',lcol,'markersize',8,'linewidth',2)
%     hold on
%     errorbar(xlevels,nolcresp(cll,:),[controlerr(cll,1), squeeze(mean(conderr(cll,1,:,:),3))'],'ko-','markersize',8,'linewidth',2)
%     xlabel('shown contrast')
%     ylabel('Firing rate [Hz]')
%     legend({'Light ON','Light OFF'})
%     set(gca,'xtick',xlevels); %[0, clevels])
%     ax = axis;
%     axis([-.1,1.1,ax(3),ax(4)])
%     line([-0.1,1.1],[0,0],'color','k')
%     title(['contrast response all orientations ']);
%     
% %     figSize = [30 21];
% %     set(gcf,'paperunits','centimeters','papersize',figSize,'paperposition',[0 0 figSize])
% %     if cll<10, printi = ['0', int2str(cll)]; else printi = int2str(cll); end
% %     print([l5RSprintpath ,  printi '__' cllname{cll} '.pdf'],'-dpdf')
% end

% cond = l5&prsv'&~phe';
% anovavec = [scafrl0(cond,1);scafrl0(cond,2);scafrl0(cond,3);scafrl0(cond,4);scafrl0(cond,5);...
%     scafrl1(cond,1);scafrl1(cond,2);scafrl1(cond,3);scafrl1(cond,4);scafrl1(cond,5)];
% gl = [ones(sum(cond)*5,1);ones(sum(cond)*5,1).*2];
% h = [ones(sum(cond),1)];
% gp = [h.*1;h.*2;h.*3;h.*4;h.*5;h.*1;h.*2;h.*3;h.*4;h.*5];
% [p,table,stats] = anovan(anovavec,{gl,gp});



cond = l4&prsv'&~phe'&ok;

figure
boundedline(respta,mean(l0meanrespn(cond,:)),std(l0meanrespn(cond,:))./sqrt(sum(cond)),'b')
hold on
boundedline(respta,mean(l1meanrespn(cond,:)),std(l1meanrespn(cond,:))./sqrt(sum(cond)),'r')
line([0,0],[.5,2],'color','k')
line([2000,2000],[.5,2],'color','k')
line([500,500],[.5,2],'color','k','linestyle',':')
line([1500,1500],[.5,2],'color','k','linestyle',':')

%ranked psths
cond = l4rs;
rank = 1;
figure
boundedline(respta,squeeze(nanmean(nscarespl0(cond,rank,:))),squeeze(nanstd(nscarespl0(cond,rank,:)))./sqrt(sum(cond)),'k')
hold on
boundedline(respta,squeeze(nanmean(nscarespl1(cond,rank,:))),squeeze(nanstd(nscarespl1(cond,rank,:)))./sqrt(sum(cond)),'r')
line([0,0],[.5,2],'color','k')
line([2000,2000],[.5,2],'color','k')
line([500,500],[.5,2],'color','k','linestyle',':')
line([1500,1500],[.5,2],'color','k','linestyle',':')

%ranked psths normalized to mean of best rank

for i = 1:5
    subplot(3,5,i)
    cond = l23fs;
    rank = i;
    boundedline(respta,squeeze(nanmean(ntomscarespl0(cond,rank,:))),squeeze(nanstd(ntomscarespl0(cond,rank,:)))./sqrt(sum(cond)),'k')
    hold on
    boundedline(respta,squeeze(nanmean(ntomscarespl1(cond,rank,:))),squeeze(nanstd(ntomscarespl1(cond,rank,:)))./sqrt(sum(cond)),'r')
    line([0,0],[0,2],'color','k')
    line([2000,2000],[0,2],'color','k')
    line([500,500],[0,2],'color','k','linestyle',':')
    line([1500,1500],[0,2],'color','k','linestyle',':')
    axis([-300,2300,0,2])
    if i == 1, ylabel('normalized FR [Hz]'); end
    
    subplot(3,5,5+i)
    cond = l4fs;
    rank = i;
    boundedline(respta,squeeze(nanmean(ntomscarespl0(cond,rank,:))),squeeze(nanstd(ntomscarespl0(cond,rank,:)))./sqrt(sum(cond)),'k')
    hold on
    boundedline(respta,squeeze(nanmean(ntomscarespl1(cond,rank,:))),squeeze(nanstd(ntomscarespl1(cond,rank,:)))./sqrt(sum(cond)),'r')
    line([0,0],[0,2],'color','k')
    line([2000,2000],[0,2],'color','k')
    line([500,500],[0,2],'color','k','linestyle',':')
    line([1500,1500],[0,2],'color','k','linestyle',':')
    axis([-300,2300,0,2])
    if i == 1, ylabel('normalized FR [Hz]'); end
    
    subplot(3,5,2*5+i)
    cond = l5fs;
    rank = i;
    boundedline(respta,squeeze(nanmean(ntomscarespl0(cond,rank,:))),squeeze(nanstd(ntomscarespl0(cond,rank,:)))./sqrt(sum(cond)),'k')
    hold on
    boundedline(respta,squeeze(nanmean(ntomscarespl1(cond,rank,:))),squeeze(nanstd(ntomscarespl1(cond,rank,:)))./sqrt(sum(cond)),'r')
    line([0,0],[0,2],'color','k')
    line([2000,2000],[0,2],'color','k')
    line([500,500],[0,2],'color','k','linestyle',':')
    line([1500,1500],[0,2],'color','k','linestyle',':')
    axis([-300,2300,0,2])
    xlabel('time [ms]')
    if i == 1, ylabel('normalized FR [Hz]'); end
end

% cond = l23&prsv'&~phe'&ok;
% figure
% errorbar(mean(sortedl0(cond,:)),std(sortedl0(cond,:))./sqrt(sum(cond)),'.')
% hold on
% errorbar(mean(sortedl1(cond,:)),std(sortedl1(cond,:))./sqrt(sum(cond)),'r.')
% xlabel('rank')
% ylabel('normalized fr')
% title('ranked firing rates')
% legend([{'L0'},{'L1'}],'location','nw')
% 
% figure
% errorbar(mean(sortedlightl0(cond,:)),std(sortedlightl0(cond,:))./sqrt(sum(cond)),'.')
% hold on
% errorbar(mean(sortedlightl1(cond,:)),std(sortedlightl1(cond,:))./sqrt(sum(cond)),'r.')
% xlabel('rank')
% ylabel('normalized fr')
% title('fr ranked by light condition')
% legend([{'L0'},{'L1'}],'location','nw')
% 
% figure
% errorbar(mean(nshufl0(cond,:)),std(nshufl0(cond,:))./sqrt(sum(cond)),'.')
% hold on
% errorbar(mean(nshufl1(cond,:)),std(nshufl1(cond,:))./sqrt(sum(cond)),'r.')
% xlabel('rank')
% ylabel('normalized fr')
% title('l0/l1 shuffled frs ranked')
% legend([{'L0'},{'L1'}],'location','nw')

%% CRFs

% normalized at preferred ori
cond = l5rs&(pc==5);
figure
errorbar(xlevels,squeeze(nanmean(nprefcrfs(cond,1,:))),squeeze(nanstd(nprefcrfs(cond,1,:)))./sqrt(sum(cond)),'k.-','linewidth',2)
hold on
errorbar(xlevels,squeeze(nanmean(nprefcrfs(cond,2,:))),squeeze(nanstd(nprefcrfs(cond,2,:)))./sqrt(sum(cond)),'r.-','linewidth',2)
% errorbar(xlevels,squeeze(nanmean(nnonprefcrfs(cond,1,:))),squeeze(nanstd(nnonprefcrfs(cond,1,:)))./sqrt(sum(cond)),'k.-')
% errorbar(xlevels,squeeze(nanmean(nnonprefcrfs(cond,2,:))),squeeze(nanstd(nnonprefcrfs(cond,2,:)))./sqrt(sum(cond)),'r.-')
xlabel('contrast level')
ylabel('normalized FR')
axis([-.2,1.2,0,1.4])

% normalized at preferred ori - RS vs FS
cond1 = l4rs; cond2 = l4fs;
figure
errorbar(xlevels,squeeze(mean(nprefcrfs(cond1,1,:))),squeeze(std(nprefcrfs(cond1,1,:)))./sqrt(sum(cond1)),'.-')
hold on
errorbar(xlevels,squeeze(mean(nprefcrfs(cond2,1,:))),squeeze(std(nprefcrfs(cond2,1,:)))./sqrt(sum(cond2)),'r.-')
xlabel('contrast level')
ylabel('normalized FR')
axis([-.2,1.2,0,1])
legend([{'RS'},{'FS'}],'location','nw')

% ranked averaged over oris
cond = l23fs;
figure
errorbar(mean(sncafrl0(cond,:)),std(sncafrl0(cond,:))./sqrt(sum(cond)),'.-')
hold on
errorbar(mean(sncafrl1(cond,:)),std(sncafrl1(cond,:))./sqrt(sum(cond)),'r.-')
xlabel('rank')
ylabel('normalized FR')
axis([.5,5.5,0,1.1])

% ranked averaged over oris - RS vs FS
figure
errorbar(mean(sncafrl0(cond,:)),mean(sncafrl1(cond,:)),std(sncafrl1(cond,:))./sqrt(sum(cond)),'r.-')
hold on
herrorbar(mean(sncafrl0(cond,:)),mean(sncafrl1(cond,:)),std(sncafrl0(cond,:))./sqrt(sum(cond)),'r.-')
line([0,1.2],[0,1.2])
axis([0,1.2,0,1.2])
xlabel('normalized ranked firing rate light OFF')
ylabel('normalized ranked firing rate light ON')
title('contrast level ranked')
axis square

% normalized mean over oris
cond = l23fs;
figure
errorbar(xlevels,squeeze(mean(normmeancrfs(cond,1,:))),squeeze(std(normmeancrfs(cond,1,:)))./sqrt(sum(cond)),'k.-','linewidth',2)
hold on
errorbar(xlevels,squeeze(mean(normmeancrfs(cond,2,:))),squeeze(std(normmeancrfs(cond,2,:)))./sqrt(sum(cond)),'r.-','linewidth',2)
xlabel('contrast level')
ylabel('normalized FR')
legend([{'light off'},{'light on'}],'location','nw')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% this one
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% normalized mean over oris - RS vs FS
cond1 = l4rs; cond2 = l23fs;
l = 0;
figure
errorbar(xlevels,squeeze(mean(normmeancrfs(cond1,l+1,:))),squeeze(std(normmeancrfs(cond1,l+1,:)))./sqrt(sum(cond1)),'.-','linewidth',2)
hold on
errorbar(xlevels,squeeze(mean(normmeancrfs(cond2,l+1,:))),squeeze(std(normmeancrfs(cond2,l+1,:)))./sqrt(sum(cond2)),'r.-','linewidth',2)
xlabel('contrast level')
ylabel('normalized FR')
axis([-.2,1.2,0,1])
legend([{'RS'},{'FS'}],'location','nw')

% both normalized mean over oris
cond = l23rs&pc==5;
figure
errorbar(xlevels,squeeze(mean(bnmeancrfs(cond,1,:))),squeeze(std(bnmeancrfs(cond,1,:)))./sqrt(sum(cond)),'.-','linewidth',2)
hold on
errorbar(xlevels,squeeze(mean(bnmeancrfs(cond,2,:))),squeeze(std(bnmeancrfs(cond,2,:)))./sqrt(sum(cond)),'r.-','linewidth',2)
xlabel('contrast level')
ylabel('normalized FR')
legend([{'light off'},{'light on'}],'location','nw')

% both normalized preferred oris
cond = l23&pc==5;
figure
errorbar(xlevels,squeeze(mean(bnprefcrfs(cond,1,:))),squeeze(std(bnprefcrfs(cond,1,:)))./sqrt(sum(cond)),'.-','linewidth',2)
hold on
errorbar(xlevels,squeeze(mean(bnprefcrfs(cond,2,:))),squeeze(std(bnprefcrfs(cond,2,:)))./sqrt(sum(cond)),'r.-','linewidth',2)
xlabel('contrast level')
ylabel('normalized FR')
legend([{'light off'},{'light on'}],'location','nw')

% normalized pref ori - RS vs FS
cond1 = l23rs&(pc==4|pc==5); cond2 = l23fs&(pc==4|pc==5);
figure
errorbar(xlevels,squeeze(mean(bnprefcrfs(cond1,1,:))),squeeze(std(bnprefcrfs(cond1,1,:)))./sqrt(sum(cond1)),'.-','linewidth',2)
hold on
errorbar(xlevels,squeeze(mean(bnprefcrfs(cond2,1,:))),squeeze(std(bnprefcrfs(cond2,1,:)))./sqrt(sum(cond2)),'r.-','linewidth',2)
xlabel('contrast level')
ylabel('normalized FR')
axis([-.2,1.2,0,1])
legend([{'RS'},{'FS'}],'location','nw')

% normalized mean over oris - RS vs FS
cond1 = l23'&kmeansind==1;%&(pc==4|pc==5); 
cond2 = l23'&kmeansind==2; %&(pc==4|pc==5);
cond3 = l23'&kmeansind==3; %&(pc==4|pc==5);
figure
errorbar(xlevels,squeeze(mean(bnmeancrfs(cond1,1,:))),squeeze(std(bnmeancrfs(cond1,1,:)))./sqrt(sum(cond1)),'.-','linewidth',2)
hold on
errorbar(xlevels,squeeze(mean(bnmeancrfs(cond2,1,:))),squeeze(std(bnmeancrfs(cond2,1,:)))./sqrt(sum(cond2)),'r.-','linewidth',2)
errorbar(xlevels,squeeze(mean(bnmeancrfs(cond3,1,:))),squeeze(std(bnmeancrfs(cond3,1,:)))./sqrt(sum(cond2)),'g.-','linewidth',2)
xlabel('contrast level')
ylabel('normalized FR')
axis([-.2,1.2,0,1])
legend([{'RS'},{'FS'},{'middle'}],'location','nw')

% non-normalized average over orientations
cond = l23rs&(pc==5);
figure
errorbar(xlevels(2:6),mean(ncafrl0(cond,:)),std(ncafrl0(cond,:))./sqrt(sum(cond)),'k','linewidth',2)
hold on
errorbar(xlevels(2:6),mean(ncafrl1(cond,:)),std(ncafrl1(cond,:))./sqrt(sum(cond)),'r','linewidth',2)
xlabel('contrast')
ylabel('normalized firing rate')

%% orienation tuning

oriscenter = circshift(oris,3,2);
oriscenter(3) = -45; oriscenter(2) = -90; oriscenter(1) = -135;

cond = l5fs;
figure
errorbar(oriscenter,mean(scl0center(cond,:)),std(scl0center(cond,:))./sqrt(length(find(cond))),'k','linewidth',2)
hold on
errorbar(oriscenter,mean(scl1center(cond,:)),std(scl1center(cond,:))./sqrt(length(find(cond))),'r','linewidth',2)
errorbar(225, mean(shiftcurvectrll0(cond)),std(shiftcurvectrll0(cond))./sqrt(sum(cond)),'k','linewidth',2)
errorbar(225, mean(shiftcurvectrll1(cond)),std(shiftcurvectrll1(cond))./sqrt(sum(cond)),'r','linewidth',2)
xlabel('degrees from optimal')
ylabel('normalized FR')
% title('L5 RS clls')
set(gca,'xtick',[-135,-90,-45,0,45,90,135,180])

conds = l23fs;
clvl = 5;
figure
errorbar(squeeze(nanmean(shiftcurveclsl0(conds,clvl,:))),...
    squeeze(nanstd(shiftcurveclsl0(conds,clvl,:)))./sqrt(length(find(conds))));
hold on
errorbar(squeeze(nanmean(shiftcurveclsl1(conds,clvl,:))),...
    squeeze(nanstd(shiftcurveclsl1(conds,clvl,:)))./sqrt(length(find(conds))),'r');

cond = l23rs
figure
errorbar(oris,mean(shiftcurvel0(cond,:)),std(shiftcurvel0(cond,:))./sqrt(length(find(cond))))
hold on
errorbar(oris,mean(shiftcurvel1(cond,:)),std(shiftcurvel1(cond,:))./sqrt(length(find(cond))),'r')
xlabel('degrees from optimal')
ylabel('normalized FR')
title('L5 RS clls')


%% firing rate changes

cond = l5rs
figure
errorbar(nanmean(scaomi(cond,:)),nanstd((scaomi(cond,:)))./(sqrt(sum(cond))),'.')
hold on
errorbar(0,nanmean(controlomi(cond)),nanstd(controlomi(cond))./sqrt(sum(cond)),'.');
line([-1,6],[0,0],'color','k')
xlabel('rank')
ylabel('OMI')

percica = scafrl1(cond,:)./scafrl0(cond,:);
percica(isinf(percica)) = NaN;
figure
errorbar(nanmean(percica),nanstd(percica)./sqrt(sum(cond)),'.-')
line([-1,6],[1,1],'color','k')
xlabel('rank')
ylabel('percent change')

spikesadded = scafrl1(cond,:)-scafrl0(cond,:);
figure
errorbar(nanmean(spikesadded),nanstd(spikesadded)./sqrt(sum(cond)),'.-')
line([-1,6],[0,0],'color','k')
xlabel('rank')
ylabel('spikes added')


cond = l5rs;
a = find(cond);
figure
hold on
for i = a
    plot([min(condspikesl0(i)),max(condspikesl1(i))],[min(condspikesl0(i))*params(i,1)+params(i,2),max(condspikesl1(i))*params(i,1)+params(i,2)])
end
xlabel('firing rate range no light')
ylabel('spikes added/subtracted by light')
title('slope of spikes added per spike')

figure
hold on
for i = a
    plot([min(normcondspikesl0(i)),max(normcondspikesl1(i))],[min(normcondspikesl0(i))*normparams(i,1)+normparams(i,2),max(normcondspikesl1(i))*normparams(i,1)+normparams(i,2)])
end
xlabel('firing rate range no light')
ylabel('spikes added/subtracted by light')
title('slope of spikes added per spike - normalized data')







% %non normalized
% figure
% errorbar(xlevels,nanmean(nolcresp(prsv'&sigcontdec'&l5&~phe',:)),...
%     nanstd(nolcresp(prsv'&sigcontdec'&l5&~phe',:))./sqrt(sum(prsv'&sigcontdec'&l5&~phe')),'.-','linewidth',2);
% hold on
% errorbar(xlevels,nanmean(lcresp(prsv'&sigcontdec'&l5&~phe',:)),...
%     nanstd(lcresp(prsv'&sigcontdec'&l5&~phe',:))./sqrt(sum(prsv'&sigcontdec'&l5&~phe')),'r.-','linewidth',2);
% xlabel('contrast level');
% ylabel('normalized firing rate')
% legend([{'light off'},{'light on'}],'location','nw')
% 
% figure
% errorbar(xlevels,nanmean(nolcresp(prsv'&sigcontinc'&depth<=500&depth>=375,:)),...
%     nanstd(nolcresp(prsv'&sigcontinc'&depth<=500&depth>=375,:))./sqrt(sum(prsv'&sigcontinc'&depth<=500&depth>=375)),'.-','linewidth',2);
% hold on
% errorbar(xlevels,nanmean(lcresp(prsv'&sigcontinc'&depth<=500&depth>=375,:)),...
%     nanstd(lcresp(prsv'&sigcontinc'&depth<=500&depth>=375,:))./sqrt(sum(prsv'&sigcontinc'&depth<=500&depth>=375)),'r.-','linewidth',2);
% xlabel('contrast level');
% ylabel('normalized firing rate')
% legend([{'light off'},{'light on'}],'location','nw')
% 
% figure
% errorbar(xlevels,nanmean(nolcresp(prsv'&sigcontinc'&depth>500&depth<=800,:)),...
%     nanstd(nolcresp(prsv'&sigcontinc'&depth>500&depth<=800,:))./sqrt(sum(prsv'&sigcontinc'&depth>500&depth<=800)),'.-','linewidth',2);
% hold on
% errorbar(xlevels,nanmean(lcresp(prsv'&sigcontinc'&depth>500&depth<=800,:)),...
%     nanstd(lcresp(prsv'&sigcontinc'&depth>500&depth<=800,:))./sqrt(sum(prsv'&sigcontinc'&depth>500&depth<=800)),'r.-','linewidth',2);
% xlabel('contrast level');
% ylabel('normalized firing rate')
% legend([{'light off'},{'light on'}],'location','nw')


% figure
% errorbar(xlevels,nanmean(nolcresp(pfsv&sigcontinc,:)),nanstd(nolcresp(pfsv&sigcontinc,:))./sqrt(sum(pfsv&sigcontinc)),'.-','linewidth',2);
% hold on
% errorbar(xlevels,nanmean(lcresp(pfsv&sigcontinc,:)),nanstd(lcresp(pfsv&sigcontinc,:))./sqrt(sum(pfsv&sigcontinc)),'r.-','linewidth',2);
% xlabel('contrast level');
% ylabel('normalized firing rate')
% legend([{'light off'},{'light on'}],'location','nw')
% 
% figure
% errorbar(xlevels,nanmean(nolcresp(prsv&sigcontdec,:)),nanstd(nolcresp(prsv&sigcontdec,:))./sqrt(sum(prsv&sigcontdec)),'.-','linewidth',2);
% hold on
% errorbar(xlevels,nanmean(lcresp(prsv&sigcontdec,:)),nanstd(lcresp(prsv&sigcontdec,:))./sqrt(sum(prsv&sigcontdec)),'r.-','linewidth',2);
% xlabel('contrast level');
% ylabel('normalized firing rate')
% legend([{'light off'},{'light on'}],'location','nw')
% 
% figure
% errorbar(xlevels,nanmean(nolcresp(pfsv&sigcontdec,:)),nanstd(nolcresp(pfsv&sigcontdec,:))./sqrt(sum(pfsv&sigcontdec)),'.-','linewidth',2);
% hold on
% errorbar(xlevels,nanmean(lcresp(pfsv&sigcontdec,:)),nanstd(lcresp(pfsv&sigcontdec,:))./sqrt(sum(pfsv&sigcontdec)),'r.-','linewidth',2);
% xlabel('contrast level');
% ylabel('normalized firing rate')
% legend([{'light off'},{'light on'}],'location','nw')
% 
% figure
% errorbar(xlevels,nanmean(nolcresp(prsv&~sigcontdec&~sigcontinc,:)),nanstd(nolcresp(prsv&~sigcontdec&~sigcontinc,:))./sqrt(sum(prsv&sigcontdec)),'.-','linewidth',2);
% hold on
% errorbar(xlevels,nanmean(lcresp(prsv&~sigcontdec&~sigcontinc,:)),nanstd(lcresp(prsv&~sigcontdec&~sigcontinc,:))./sqrt(sum(prsv&sigcontdec)),'r.-','linewidth',2);
% xlabel('contrast level');
% ylabel('normalized firing rate')
% legend([{'light off'},{'light on'}],'location','nw')
% 
% figure
% errorbar(xlevels,nanmean(nolcresp(pfsv&~sigcontdec&~sigcontinc,:)),nanstd(nolcresp(pfsv&~sigcontdec&~sigcontinc,:))./sqrt(sum(pfsv&sigcontdec)),'.-','linewidth',2);
% hold on
% errorbar(xlevels,nanmean(lcresp(pfsv&~sigcontdec&~sigcontinc,:)),nanstd(lcresp(pfsv&~sigcontdec&~sigcontinc,:))./sqrt(sum(pfsv&sigcontdec)),'r.-','linewidth',2);
% xlabel('contrast level');
% ylabel('normalized firing rate')
% legend([{'light off'},{'light on'}],'location','nw')



% population PSTHs
figure
boundedline(respta,mean(l0meanrespn(l23&prsv',:)),std(l0meanrespn(l23&prsv',:))./sqrt(sum(l23&prsv')))
hold on
boundedline(respta,mean(l1meanrespn(l23&prsv',:)),std(l1meanrespn(l23&prsv',:))./sqrt(sum(l23&prsv')),'r')

figure
boundedline(respta,mean(l0meanrespn(l23&pfsv',:)),std(l0meanrespn(l23&pfsv',:))./sqrt(sum(l23&pfsv')),'b')
hold on
boundedline(respta,mean(l1meanrespn(l23&pfsv',:)),std(l1meanrespn(l23&pfsv',:))./sqrt(sum(l23&pfsv')),'r')

figure
boundedline(respta,mean(l0prefrespn(l23&prsv',:)),std(l0prefrespn(l23&prsv',:))./sqrt(sum(l23&prsv')))
hold on
boundedline(respta,mean(l1prefrespn(l23&prsv',:)),std(l1prefrespn(l23&prsv',:))./sqrt(sum(l23&prsv')),'r')

figure
boundedline(respta,mean(l0meanrespn(l4&pfsv',:)),std(l0meanrespn(l4&pfsv',:))./sqrt(sum(l4&pfsv')),'b')
hold on
boundedline(respta,mean(l1meanrespn(l4&pfsv',:)),std(l1meanrespn(l4&pfsv',:))./sqrt(sum(l4&pfsv')),'r')

figure
boundedline(respta,mean(l0meanrespn(l4&prsv',:)),std(l0meanrespn(l4&prsv',:))./sqrt(sum(l4&prsv')),'b')
hold on
boundedline(respta,mean(l1meanrespn(l4&prsv',:)),std(l1meanrespn(l4&prsv',:))./sqrt(sum(l4&prsv')),'r')

figure
boundedline(respta,mean(l0meanrespn(l5&prsv'&~phe',:)),std(l0meanrespn(l5&prsv'&~phe',:))./sqrt(sum(l5&prsv'&~phe')),'b')
hold on
boundedline(respta,mean(l1meanrespn(l5&prsv'&~phe',:)),std(l1meanrespn(l5&prsv'&~phe',:))./sqrt(sum(l5&prsv'&~phe')),'r')

figure
boundedline(respta,mean(l0meanrespn(l5&pfsv'&~phe',:)),std(l0meanrespn(l5&pfsv'&~phe',:))./sqrt(sum(l5&pfsv'&~phe')),'b')
hold on
boundedline(respta,mean(l1meanrespn(l5&pfsv'&~phe',:)),std(l1meanrespn(l5&pfsv'&~phe',:))./sqrt(sum(l5&pfsv'&~phe')),'r')

figure
boundedline(respta,mean(l0prefrespn(l5&prsv'&~phe',:)),std(l0prefrespn(l5&prsv'&~phe',:))./sqrt(sum(l5&prsv'&~phe')),'b')
hold on
boundedline(respta,mean(l1prefrespn(l5&prsv'&~phe',:)),std(l1prefrespn(l5&prsv'&~phe',:))./sqrt(sum(l5&prsv'&~phe')),'r')

figure
boundedline(respta,mean(l0prefrespn(l5&pfsv'&~phe',:)),std(l0prefrespn(l5&pfsv'&~phe',:))./sqrt(sum(l5&pfsv'&~phe')),'b')
hold on
boundedline(respta,mean(l1prefrespn(l5&pfsv'&~phe',:)),std(l1prefrespn(l5&pfsv'&~phe',:))./sqrt(sum(l5&pfsv'&~phe')),'r')

% fsl0cresps = nolcresp(pfs,:);
% fsl1cresps = lcresp(pfs,:);
% rsl0cresps = nolcresp(prs,:);
% rsl1cresps = lcresp(prs,:);
% for i = 1:size(fsl0cresps,1)
%     fsl0cresps(i,:) = fsl0cresps(i,:)-fsl0cresps(i,1);
%     fsl0cresps(i,:) = fsl0cresps(i,:)./max(fsl0cresps(i,:));
%     if find(isnan(fsl0cresps(i,:))), fsl0cresps(i,find(isinf(fsl0cresps(i,:)))) = NaN; end
%     fsl1cresps(i,:) = fsl1cresps(i,:)-fsl1cresps(i,1);
%     fsl1cresps(i,:) = fsl1cresps(i,:)./max(fsl1cresps(i,:));
%     if find(isnan(fsl1cresps(i,:))), fsl1cresps(i,find(isinf(fsl1cresps(i,:)))) = NaN; end
%     if find(fsl1cresps(i,:)<-10), fsl1cresps(i,find(fsl1cresps(i,:)<-10)) = NaN; end
% end
% for i = 1:size(rsl0cresps,1)
%     rsl0cresps(i,:) = rsl0cresps(i,:)-rsl0cresps(i,1);
%     rsl0cresps(i,:) = rsl0cresps(i,:)./max(rsl0cresps(i,:));
%     if find(isnan(rsl0cresps(i,:))), rsl0cresps(i,find(isinf(rsl0cresps(i,:)))) = NaN; end
%     rsl1cresps(i,:) = rsl1cresps(i,:)-rsl1cresps(i,1);
%     rsl1cresps(i,:) = rsl1cresps(i,:)./max(rsl1cresps(i,:));
%     if find(isnan(rsl1cresps(i,:))), rsl1cresps(i,find(isinf(rsl1cresps(i,:)))) = NaN; end
%     if find(rsl1cresps(i,:)<-10), rsl1cresps(i,find(rsl1cresps(i,:)<-10)) = NaN; end
% end
% 
% figure
% errorbar(xlevels,nanmean(rsl0cresps,1),nanstd(rsl0cresps)./(sqrt(size(rsl0cresps,1)-length(find(isnan(rsl0cresps(:,1)))))),'.-')
% hold on
% errorbar(xlevels,nanmean(fsl0cresps,1),nanstd(fsl0cresps)./(sqrt(size(fsl0cresps,1)-length(find(isnan(fsl0cresps(:,1)))))),'r.-')
% title('no light average normalized contrast response curves')
% xlabel('contrast level')
% ylabel('normalized firing rate')
% legend([{'RS'},{'FS'}],'location','nw')
% 
% figure
% errorbar(xlevels,nanmean(rsl1cresps,1),nanstd(rsl1cresps)./(sqrt(size(rsl1cresps,1)-length(find(isnan(rsl1cresps(:,1)))))),'b.-')
% hold on
% errorbar(xlevels,nanmean(fsl1cresps,1),nanstd(fsl1cresps)./(sqrt(size(fsl1cresps,1)-length(find(isnan(fsl1cresps(:,1)))))),'r.-')
% title('light on average normalized contrast response curves')
% xlabel('contrast level')
% ylabel('normalized firing rate')
% legend([{'RS'},{'FS'}],'location','nw')

% Scotts paper tuning
vl5rs = l5rs&vismod;
sbcc = vl5rs'&controlfr(:,1)>nlfr(:,1); % suppressed by contrast clls
gl5rs = vl5rs'&controlfr(:,1)<nlfr(:,1); %good layer 5 rs


% Scotts paper
comi = (controlfr(:,2)-controlfr(:,1))./(controlfr(:,2)+controlfr(:,1));
omi = (lfr-nlfr)./(lfr+nlfr);
rss = prsv' & ~phe' & ok;
fss = pfsv' &ok;
lightmod(isnan(lightmod)) = 0;

figure
plot(comi(rss),depth(rss),'o','markersize',4,'markerfacecolor','b')
hold on
plot(comi(rss&lightmod),depth(rss&lightmod),'ro','markersize',4,'markerfacecolor','r')
axis ij
line([-1,1],[375,375],'color','k','linestyle',':');
line([-1,1],[550,550],'color','k','linestyle',':');
line([-1,1],[800,800],'color','k','linestyle',':');
line([0,0],[0,1000],'color','k','linewidth',2);
title('spontaneous OMI RS clls')
xlabel('OMI')
ylabel('depth')

figure
plot(comi(fss),depth(fss),'o','markersize',4,'markerfacecolor','b')
hold on
plot(comi(fss&lightmod),depth(fss&lightmod),'ro','markersize',4,'markerfacecolor','r')
axis ij
line([-1,1],[375,375],'color','k','linestyle',':');
line([-1,1],[550,550],'color','k','linestyle',':');
line([-1,1],[800,800],'color','k','linestyle',':');
line([0,0],[0,1000],'color','k','linewidth',2);
title('spontaneous OMI FS clls')
xlabel('OMI')
ylabel('depth')

figure
plot(omi(rss),depth(rss),'o','markersize',4,'markerfacecolor','b')
hold on
plot(omi(rss&lightmod),depth(rss&lightmod),'ro','markersize',4,'markerfacecolor','r')
axis ij
line([-1,1],[375,375],'color','k','linestyle',':');
line([-1,1],[550,550],'color','k','linestyle',':');
line([-1,1],[800,800],'color','k','linestyle',':');
line([0,0],[0,1000],'color','k','linewidth',2);
title('average visual OMI RS clls')
xlabel('OMI')
ylabel('depth')

figure
plot(omi(fss),depth(fss),'o','markersize',4,'markerfacecolor','b')
hold on
plot(omi(fss&lightmod),depth(fss&lightmod),'ro','markersize',4,'markerfacecolor','r')
axis ij
line([-1,1],[375,375],'color','k','linestyle',':');
line([-1,1],[550,550],'color','k','linestyle',':');
line([-1,1],[800,800],'color','k','linestyle',':');
line([0,0],[0,1000],'color','k','linewidth',2);
title('average visual OMI FS clls')
xlabel('OMI')
ylabel('depth')

% one bar per layer plot
bardepthsnew = [0,375;375,550;550,800;800,1200];
rsdepth = depth(prsv'&~phe'); fsdepth = depth(pfsv'&~phe');
for i = 1:size(bardepthsnew,1)
    rsbars(i) = nanmean(rsomi(rsdepth>=bardepthsnew(i,1) & rsdepth<bardepthsnew(i,2)));
    rsbarerr(i) = nanstd(rsomi(rsdepth>=bardepthsnew(i,1) & rsdepth<bardepthsnew(i,2)))./sqrt(length(find(rsdepth>=bardepthsnew(i,1) & rsdepth<bardepthsnew(i,2))));
    rsn(i) = length(find(rsdepth>=bardepthsnew(i,1) & rsdepth<bardepthsnew(i,2)));
    fsbars(i) = mean(rsomi(fsdepth>=bardepthsnew(i,1) & fsdepth<bardepthsnew(i,2)));
    fsbarerr(i) = std(rsomi(fsdepth>=bardepthsnew(i,1) & fsdepth<bardepthsnew(i,2)))./sqrt(length(find(fsdepth>=bardepthsnew(i,1) & fsdepth<bardepthsnew(i,2))));
end

figure
barh(1:4,rsbars,'r','EdgeColor','r')
hold on
herrorbar(rsbars,1:4,rsbarerr,'k.')
% line([-.3,.3],[500,500],'color','k','linestyle',':');
% line([-.3,.3],[375,375],'color','k','linestyle',':');

omi = (lfr-nlfr)./(lfr+nlfr);
% bardepths = 250:50:800;
binwidth = 50;
% bardepths = 175:binwidth:925;
bardepths = 250:binwidth:850;
for i = 1:length(bardepths)
    rsbars(i) = nanmean(omi(depth>=bardepths(i) & depth<bardepths(i)+binwidth & prsv' & ~phe' & ok));
    rsbarerr(i) = nanstd(omi(depth>=bardepths(i) & depth<bardepths(i)+binwidth & prsv' & ~phe' & ok))./sqrt(length(find(depth>=bardepths(i) & depth<bardepths(i)+50 & prsv' & ~phe' & ok)));
    rsn(i) = sum(depth>=bardepths(i) & depth<bardepths(i)+binwidth & prsv' & ~phe' & ok);
    fsbars(i) = nanmean(omi(depth>=bardepths(i) & depth<bardepths(i)+binwidth & pfsv' & ok));
    fsbarerr(i) = nanstd(omi(depth>=bardepths(i) & depth<bardepths(i)+binwidth & pfsv' & ok))./sqrt(length(find(depth>=bardepths(i) & depth<bardepths(i)+50 & pfsv' & ok)));
end
bardepths = bardepths + 25; % to center for plot

figure
barh(bardepths,rsbars,'r','EdgeColor','r')
hold on
herrorbar(rsbars,bardepths,rsbarerr,'k.')
line([-.6,.6],[550,550],'color','k','linestyle',':');
line([-.6,.6],[375,375],'color','k','linestyle',':');
line([-.6,.6],[800,800],'color','k','linestyle',':');
axis ij
axis([-.6,.6,250,900])
xlabel('average OMI')
ylabel('cortical depth [mum]')
set(gca,'box','off')
title('RS clls')


% Scott's Paper
rsomi = (lfr(prsv&~phe)-nlfr(prsv&~phe))./(lfr(prsv&~phe)+nlfr(prsv&~phe));
% rsomi = (lfr(prsv&~phe&ok')-nlfr(prsv&~phe&ok'))./(lfr(prsv&~phe&ok')+nlfr(prsv&~phe&ok'));
% rsvdomi = (lfr(prsv&~phe&sigcontinc)-nlfr(prsv&~phe&sigcontinc))./(lfr(prsv&~phe&sigcontinc)+nlfr(prsv&~phe&sigcontinc));
% plot(rsomi,depth(prsv&~phe&ok'),'o','markersize',4,'markerfacecolor',[.3,.3,.3],'color',[.3,.3,.3])
plot(rsomi,depth(prsv&~phe),'o','markersize',4,'markerfacecolor',[.3,.3,.3],'color',[.3,.3,.3])
line([0,0],[0,1000],'color','k')
line([-1,1,],[375,375],'color','k','linestyle',':')
line([-1,1,],[500,500],'color','k','linestyle',':')
line([-1,1,],[800,800],'color','k','linestyle',':')
hold on
% [x,y,xerr] = runningMedian(depth(prsv&~phe&ok'),rsomi,0,12);
[x,y,xerr] = runningMedian(depth(prsv&~phe),rsomi,0,12);
plot(x,y,'linewidth',2,'color',[.3,.3,.3]);
plot(x+xerr,y,'color',[.3,.3,.3])
plot(x-xerr,y,'color',[.3,.3,.3])
% plot((lfr(prsv&phe&ok')-nlfr(prsv&phe&ok'))./(lfr(prsv&phe&ok')+nlfr(prsv&phe&ok')),depth(prsv&phe&ok'),'mo','markersize',4,'markerfacecolor','m')
plot((lfr(prsv&phe)-nlfr(prsv&phe))./(lfr(prsv&phe)+nlfr(prsv&phe)),depth(prsv&phe),'mo','markersize',4,'markerfacecolor','m')
axis ij
axis([-1.1,1.1,150,950])
ylabel('depth[mum]')
xlabel('OMI of average visual response rate')
title('RS clls: average firing rate changes by layer 4 suppression')
set(gca,'box','off')

figure
% fsomi = (lfr(pfsv)-nlfr(pfsv))./(lfr(pfsv)+nlfr(pfsv));
fsomi = (lfr(pfsv&ok')-nlfr(pfsv&ok'))./(lfr(pfsv&ok')+nlfr(pfsv&ok'));
% fsvdomi = (lfr(pfsv&~phe&sigcontinc)-nlfr(pfsv&~phe&sigcontinc))./(lfr(pfsv&~phe&sigcontinc)+nlfr(pfsv&~phe&sigcontinc));
plot(fsomi,depth(pfsv&ok'),'o','markersize',4,'markerfacecolor','r','color','r')
line([0,0],[0,1000],'color','k')
line([-1,1,],[375,375],'color','k','linestyle',':')
line([-1,1,],[500,500],'color','k','linestyle',':')
hold on
[x,y,xerr] = runningMedian(depth(pfsv&ok'),fsomi,0,12);
plot(x,y,'linewidth',2,'color','r');
plot(x+xerr,y,'color','r')
plot(x-xerr,y,'color','r')
axis ij
axis([-1.1,1.1,150,950])
ylabel('depth[mum]')
xlabel('OMI of average visual response rate')
title('FS clls: average firing rate changes by layer 4 suppression')
set(gca,'box','off')

rsdepth = depth(prsv&~phe&ok'); fsdepth = depth(pfsv&ok');
bardepths = 250:50:800;
for i = 1:length(bardepths)
    rsbars(i) = mean(rsomi(rsdepth>=bardepths(i) & rsdepth<bardepths(i)+50));
    rsbarerr(i) = std(rsomi(rsdepth>=bardepths(i) & rsdepth<bardepths(i)+50))./sqrt(length(find(rsdepth>=bardepths(i) & rsdepth<bardepths(i)+50)));
    fsbars(i) = mean(fsomi(fsdepth>=bardepths(i) & fsdepth<bardepths(i)+50));
    fsbarerr(i) = std(fsomi(fsdepth>=bardepths(i) & fsdepth<bardepths(i)+50))./sqrt(length(find(fsdepth>=bardepths(i) & fsdepth<bardepths(i)+50)));
end
    
figure
barh(bardepths,rsbars,'r','EdgeColor','r')
hold on
herrorbar(rsbars,bardepths,rsbarerr,'k.')
line([-.3,.3],[500,500],'color','k','linestyle',':');
line([-.3,.3],[375,375],'color','k','linestyle',':');
axis ij
axis([-.3,.3,200,800])
xlabel('average OMI')
ylabel('cortical depth [mum]')
set(gca,'box','off')
title('RS clls')

figure
barh(bardepths,fsbars,'r')
hold on
herrorbar(fsbars,bardepths,fsbarerr,'k.')
line([-.5,.5],[500,500],'color','k','linestyle',':');
line([-.5,.5],[375,375],'color','k','linestyle',':');
axis ij
xlabel('average OMI')
ylabel('cortical depth [mum]')
set(gca,'box','off')
title('FS clls')

clear rsbars; clear fsbars; clear rsbarerr; clear fsbarerr;
bardepthsnew = [0,375;375,500;500,800];
for i = 1:size(bardepthsnew,1)
    rsbars(i) = mean(rsomi(rsdepth>=bardepthsnew(i,1) & rsdepth<bardepthsnew(i,2)));
    rsbarerr(i) = std(rsomi(rsdepth>=bardepthsnew(i,1) & rsdepth<bardepthsnew(i,2)))./sqrt(length(find(rsdepth>=bardepthsnew(i,1) & rsdepth<bardepthsnew(i,2))));
    fsbars(i) = mean(rsomi(fsdepth>=bardepthsnew(i,1) & fsdepth<bardepthsnew(i,2)));
    fsbarerr(i) = std(rsomi(fsdepth>=bardepthsnew(i,1) & fsdepth<bardepthsnew(i,2)))./sqrt(length(find(fsdepth>=bardepthsnew(i,1) & fsdepth<bardepthsnew(i,2))));
end

figure
barh([300,437,700],rsbars,'r','EdgeColor','r')
hold on
herrorbar(rsbars,[300,437,700],rsbarerr,'k.')
% line([-.3,.3],[500,500],'color','k','linestyle',':');
% line([-.3,.3],[375,375],'color','k','linestyle',':');
axis ij
% axis([-.3,.3,200,800])
xlabel('average OMI')
ylabel('cortical depth [mum]')
set(gca,'box','off')
title('RS clls')

    
figure
barh(1:3,fsbars,'r')
hold on
herrorbar(fsbars,1:3,fsbarerr,'k.')
% line([-.5,.5],[500,500],'color','k','linestyle',':');
% line([-.5,.5],[375,375],'color','k','linestyle',':');
axis ij
xlabel('average OMI')
ylabel('cortical depth [mum]')
set(gca,'box','off')
title('FS clls')

 % preferred stimulus
figure
rsomi = (preffrl1(prsv&~phe)-preffrl0(prsv&~phe))./(preffrl1(prsv&~phe)+preffrl0(prsv&~phe));
plot(rsomi,depth(prsv&~phe),'o','markersize',4,'markerfacecolor',[.3,.3,.3],'color',[.3,.3,.3])
line([0,0],[0,1000],'color','k')
line([-1,1,],[375,375],'color','k','linestyle',':')
line([-1,1,],[500,500],'color','k','linestyle',':')
hold on
[x,y,xerr] = runningMedian(depth(prsv&~phe),rsomi,0,12);
plot(x,y,'linewidth',2,'color',[.3,.3,.3]);
plot(x+xerr,y,'color',[.3,.3,.3])
plot(x-xerr,y,'color',[.3,.3,.3])
plot((preffrl1(prsv&phe)-preffrl0(prsv&phe))./(preffrl1(prsv&phe)+preffrl0(prsv&phe)),depth(prsv&phe),'mo','markersize',4,'markerfacecolor','m')
axis ij
axis([-1.1,1.1,150,950])
ylabel('depth[mum]')
xlabel('OMI of preferred condition visual response rate')
title('RS clls: preferred firing rate changes by layer 4 suppression')
set(gca,'box','off')

figure
fsomi = (preffrl1(pfs)-preffrl0(pfs))./(preffrl1(pfs)+preffrl0(pfs));
plot(fsomi,depth(pfs),'o','markersize',4,'markerfacecolor','r','color','r')
line([0,0],[0,1000],'color','k')
line([-1,1,],[375,375],'color','k','linestyle',':')
line([-1,1,],[500,500],'color','k','linestyle',':')
hold on
[x,y,xerr] = runningMedian(depth(pfs),fsomi,0,12);
plot(x,y,'linewidth',2,'color','r');
plot(x+xerr,y,'color','r')
plot(x-xerr,y,'color','r')
axis ij
axis([-1.1,1.1,150,950])
ylabel('depth[mum]')
xlabel('OMI of preferred condition visual response rate')
title('FS clls: preferred firing rate changes by layer 4 suppression')
set(gca,'box','off')


% running

figure
plot(r0omi(prsv),depth(prsv),'o','markersize',4,'markerfacecolor','b','color','b')
line([0,0],[0,1000],'color','k')
hold on
[x,y,xerr] = runningMedian(depth(prsv&~phe),r0omi(prsv&~phe),0,12);
plot(x,y,'linewidth',2,'color','b');
plot(x+xerr,y,'color','b')
plot(x-xerr,y,'color','b')
line([-1,1,],[375,375],'color','k','linestyle',':')
line([-1,1,],[500,500],'color','k','linestyle',':')
plot(r0omi(prsv&phe),depth(prsv&phe),'mo','markersize',4,'markerfacecolor','m')
axis ij
axis([-1.1,1.1,150,950])
ylabel('depth[mum]')
xlabel('OMI of average visual response rate')
title('Quiet: RS clls: average firing rate changes by layer 4 suppression')
set(gca,'box','off')

figure
plot(r1omi(prsv),depth(prsv),'o','markersize',4,'markerfacecolor','b','color','b')
line([0,0],[0,1000],'color','k')
hold on
[x,y,xerr] = runningMedian(depth(prsv&~phe),r1omi(prsv&~phe),0,12);
plot(x,y,'linewidth',2,'color','b');
plot(x+xerr,y,'color','b')
plot(x-xerr,y,'color','b')
line([-1,1,],[375,375],'color','k','linestyle',':')
line([-1,1,],[500,500],'color','k','linestyle',':')
plot(r1omi(prsv&phe),depth(prsv&phe),'mo','markersize',4,'markerfacecolor','m')
axis ij
axis([-1.1,1.1,150,950])
ylabel('depth[mum]')
xlabel('OMI of average visual response rate')
title('Running: RS clls: average firing rate changes by layer 4 suppression')
set(gca,'box','off')

figure
plot(r0omi(pfs),depth(pfs),'o','markersize',4,'markerfacecolor','r','color','r')
line([0,0],[0,1000],'color','k')
hold on
[x,y,xerr] = runningMedian(depth(pfs),r0omi(pfs),0,12);
plot(x,y,'linewidth',2,'color','r');
plot(x+xerr,y,'color','r')
plot(x-xerr,y,'color','r')
line([-1,1,],[375,375],'color','k','linestyle',':')
line([-1,1,],[500,500],'color','k','linestyle',':')
axis ij
axis([-1.1,1.1,150,950])
ylabel('depth[mum]')
xlabel('OMI of average visual response rate')
title('Quiet: FS clls: average firing rate changes by layer 4 suppression')
set(gca,'box','off')

figure
plot(r1omi(pfs),depth(pfs),'o','markersize',4,'markerfacecolor','r','color','r')
line([0,0],[0,1000],'color','k')
hold on
[x,y,xerr] = runningMedian(depth(pfs),r1omi(pfs),0,12);
plot(x,y,'linewidth',2,'color','r');
plot(x+xerr,y,'color','r')
plot(x-xerr,y,'color','r')
line([-1,1,],[375,375],'color','k','linestyle',':')
line([-1,1,],[500,500],'color','k','linestyle',':')
axis ij
axis([-1.1,1.1,150,950])
ylabel('depth[mum]')
xlabel('OMI of average visual response rate')
title('Running: FS clls: average firing rate changes by layer 4 suppression')
set(gca,'box','off')

% RMI

figure
plot(l0rmi(prsv),depth(prsv),'o','markersize',4,'markerfacecolor','b','color','b')
line([0,0],[0,1000],'color','k')
hold on
[x,y,xerr] = runningMedian(depth(prsv&~phe),l0rmi(prsv&~phe),0,12);
plot(x,y,'linewidth',2,'color','b');
plot(x+xerr,y,'color','b')
plot(x-xerr,y,'color','b')
line([-1,1,],[375,375],'color','k','linestyle',':')
line([-1,1,],[500,500],'color','k','linestyle',':')
plot(l0rmi(prsv&phe),depth(prsv&phe),'mo','markersize',4,'markerfacecolor','m')
axis ij
axis([-1.1,1.1,150,950])
ylabel('depth[mum]')
xlabel('RMI of average visual response rate')
title('NO light: RS clls: average firing rate changes by layer 4 suppression')
set(gca,'box','off')

figure
plot(l1rmi(prsv),depth(prsv),'o','markersize',4,'markerfacecolor','b','color','b')
line([0,0],[0,1000],'color','k')
hold on
[x,y,xerr] = runningMedian(depth(prsv&~phe),l1rmi(prsv&~phe),0,12);
plot(x,y,'linewidth',2,'color','b');
plot(x+xerr,y,'color','b')
plot(x-xerr,y,'color','b')
line([-1,1,],[375,375],'color','k','linestyle',':')
line([-1,1,],[500,500],'color','k','linestyle',':')
plot(l1rmi(prsv&phe),depth(prsv&phe),'mo','markersize',4,'markerfacecolor','m')
axis ij
axis([-1.1,1.1,150,950])
ylabel('depth[mum]')
xlabel('RMI of average visual response rate')
title('Light ON: RS clls: average firing rate changes by layer 4 suppression')
set(gca,'box','off')

figure
plot(l0rmi(pfs),depth(pfs),'o','markersize',4,'markerfacecolor','r','color','r')
line([0,0],[0,1000],'color','k')
hold on
[x,y,xerr] = runningMedian(depth(pfs),l0rmi(pfs),0,12);
plot(x,y,'linewidth',2,'color','r');
plot(x+xerr,y,'color','r')
plot(x-xerr,y,'color','r')
line([-1,1,],[375,375],'color','k','linestyle',':')
line([-1,1,],[500,500],'color','k','linestyle',':')
axis ij
axis([-1.1,1.1,150,950])
ylabel('depth[mum]')
xlabel('RMI of average visual response rate')
title('NO Light: FS clls: average firing rate changes by layer 4 suppression')
set(gca,'box','off')

figure
plot(l1rmi(pfs),depth(pfs),'o','markersize',4,'markerfacecolor','r','color','r')
line([0,0],[0,1000],'color','k')
hold on
[x,y,xerr] = runningMedian(depth(pfs),l1rmi(pfs),0,12);
plot(x,y,'linewidth',2,'color','r');
plot(x+xerr,y,'color','r')
plot(x-xerr,y,'color','r')
line([-1,1,],[375,375],'color','k','linestyle',':')
line([-1,1,],[500,500],'color','k','linestyle',':')
axis ij
axis([-1.1,1.1,150,950])
ylabel('depth[mum]')
xlabel('RMI of average visual response rate')
title('Light ON: FS clls: average firing rate changes by layer 4 suppression')
set(gca,'box','off')


% phase locking

% LFP phas stuff
% RS clls
figure
subplot(2,2,1)
imagesc(3:5:98,rsd,allrl0(prs(rsi),:));
caxis([0,.5])
colorbar
title('RS clls: r of spike phases - light OFF')
xlabel('frequency bands')
ylabel('cll number - sorted by depth')

subplot(2,2,2)
imagesc(3:5:98,rsd,allrl1(prs(rsi),:));
caxis([0,.5])
colorbar
title('RS clls: r of spike phases - light ON')
xlabel('frequency bands')
ylabel('cll number - sorted by depth')

subplot(2,2,3)
imagesc(3:5:98,fsd,allrl0(pfs(fsi),:));
caxis([0,.5])
colorbar
title('FS clls: r of spike phases - light OFF')
xlabel('frequency bands')
ylabel('cll number - sorted by depth')

subplot(2,2,4)
imagesc(3:5:98,fsd,allrl1(pfs(fsi),:));
caxis([0,.5])
colorbar
title('FS clls: r of spike phases - light ON')
xlabel('frequency bands')
ylabel('cll number - sorted by depth')


figure
subplot(2,2,1)
errorbar(3:5:98,nanmean(allrl0(pfs(l4fs),:)),nanstd(allrl0(pfs(l4fs),:))./sqrt(length(l4fs)),'r')
hold on
errorbar(3:5:98,nanmean(allrl0(prs(l4rs),:)),nanstd(allrl0(prs(l4rs),:))./sqrt(length(l4rs)),'b')
legend([{'FS'},{'RS'}]);
xlabel('frequency [Hz]')
ylabel('average phase of spike');
title('Light OFF: Layer 4');
axis([0,100,0,.3])

subplot(2,2,2)
errorbar(3:5:98,mean(allrl0(pfs(l23fs),:)),std(allrl0(pfs(l23fs),:))./sqrt(length(l23fs)),'r')
hold on
errorbar(3:5:98,mean(allrl0(prs(l23rs),:)),std(allrl0(prs(l23rs),:))./sqrt(length(l23rs)),'b')
legend([{'FS'},{'RS'}]);
xlabel('frequency [Hz]')
ylabel('average phase of spike');
title('Light OFF: Layer 2/3');
axis([0,100,0,.3])

subplot(2,2,3)
errorbar(3:5:98,mean(allrl0(pfs(l5fs),:)),std(allrl0(pfs(l5fs),:))./sqrt(length(l5fs)),'r')
hold on
errorbar(3:5:98,mean(allrl0(prs(l5rs),:)),std(allrl0(prs(l5rs),:))./sqrt(length(l5rs)),'b')
legend([{'FS'},{'RS'}]);
xlabel('frequency [Hz]')
ylabel('average phase of spike');
title('Light OFF: Layer 5');
axis([0,100,0,.3])


figure
subplot(2,2,1)
errorbar(3:5:98,mean(allrl1(pfs(l4fs),:)),std(allrl1(pfs(l4fs),:))./sqrt(length(l4fs)),'r')
hold on
errorbar(3:5:98,mean(allrl1(prs(l4rs),:)),std(allrl1(prs(l4rs),:))./sqrt(length(l4rs)),'b')
legend([{'FS'},{'RS'}]);
xlabel('frequency [Hz]')
ylabel('average phase of spike');
title('Light ON: Layer 4');
axis([0,100,0,.3])

subplot(2,2,2)
errorbar(3:5:98,mean(allrl1(pfs(l23fs),:)),std(allrl1(pfs(l23fs),:))./sqrt(length(l23fs)),'r')
hold on
errorbar(3:5:98,mean(allrl1(prs(l23rs),:)),std(allrl1(prs(l23rs),:))./sqrt(length(l23rs)),'b')
legend([{'FS'},{'RS'}]);
xlabel('frequency [Hz]')
ylabel('average phase of spike');
title('Light ON: Layer 2/3');
axis([0,100,0,.3])

subplot(2,2,3)
errorbar(3:5:98,mean(allrl1(pfs(l5fs),:)),std(allrl1(pfs(l5fs),:))./sqrt(length(l5fs)),'r')
hold on
errorbar(3:5:98,mean(allrl1(prs(l5rs),:)),std(allrl1(prs(l5rs),:))./sqrt(length(l5rs)),'b')
legend([{'FS'},{'RS'}]);
xlabel('frequency [Hz]')
ylabel('average phase of spike');
title('Light ON: Layer 5');
axis([0,100,0,.3])


figure
subplot(2,2,1)
errorbar(3:5:98,nanmean(allrl0(prs(l4rs),:)),nanstd(allrl0(prs(l4fs),:))./sqrt(length(l4rs)),'b')
hold on
errorbar(3:5:98,nanmean(allrl1(prs(l4rs),:)),nanstd(allrl1(prs(l4rs),:))./sqrt(length(l4rs)),'r')
legend([{'L0'},{'L1'}]);
xlabel('frequency [Hz]')
ylabel('average phase of spike');
title('RS: Layer 4');
axis([0,100,0,.3])

subplot(2,2,2)
errorbar(3:5:98,mean(allrl0(prs(l23rs),:)),std(allrl0(prs(l23rs),:))./sqrt(length(l23rs)),'b')
hold on
errorbar(3:5:98,mean(allrl1(prs(l23rs),:)),std(allrl1(prs(l23rs),:))./sqrt(length(l23rs)),'r')
legend([{'L0'},{'L1'}]);
xlabel('frequency [Hz]')
ylabel('average phase of spike');
title('RS: Layer 2/3');
axis([0,100,0,.3])

subplot(2,2,3)
errorbar(3:5:98,mean(allrl0(prs(l5rs),:)),std(allrl0(prs(l5rs),:))./sqrt(length(l5rs)),'b')
hold on
errorbar(3:5:98,mean(allrl1(prs(l5rs),:)),std(allrl1(prs(l5rs),:))./sqrt(length(l5rs)),'r')
legend([{'L0'},{'L1'}]);
xlabel('frequency [Hz]')
ylabel('average phase of spike');
title('RS: Layer 5');
axis([0,100,0,.3])

figure
subplot(2,2,1)
errorbar(3:5:98,nanmean(allrl0(pfs(l4fs),:)),nanstd(allrl0(pfs(l4fs),:))./sqrt(length(l4fs)),'b')
hold on
errorbar(3:5:98,nanmean(allrl1(pfs(l4fs),:)),nanstd(allrl1(pfs(l4fs),:))./sqrt(length(l4fs)),'r')
legend([{'L0'},{'L1'}]);
xlabel('frequency [Hz]')
ylabel('average phase of spike');
title('FS: Layer 4');
axis([0,100,0,.3])

subplot(2,2,2)
errorbar(3:5:98,mean(allrl0(pfs(l23fs),:)),std(allrl0(pfs(l23fs),:))./sqrt(length(l23fs)),'b')
hold on
errorbar(3:5:98,mean(allrl1(pfs(l23fs),:)),std(allrl1(pfs(l23fs),:))./sqrt(length(l23fs)),'r')
legend([{'L0'},{'L1'}]);
xlabel('frequency [Hz]')
ylabel('average phase of spike');
title('FS: Layer 2/3');
axis([0,100,0,.3])

subplot(2,2,3)
errorbar(3:5:98,mean(allrl0(pfs(l5fs),:)),std(allrl0(pfs(l5fs),:))./sqrt(length(l5fs)),'b')
hold on
errorbar(3:5:98,mean(allrl1(pfs(l5fs),:)),std(allrl1(pfs(l5fs),:))./sqrt(length(l5fs)),'r')
legend([{'L0'},{'L1'}]);
xlabel('frequency [Hz]')
ylabel('average phase of spike');
title('FS: Layer 5');
axis([0,100,0,.3])

figure
subplot(2,2,1)
imagesc(3:5:98,rsd,allrl1(prs(rsi),:)-allrl0(prs(rsi),:));
caxis([-.2,.2])
colorbar
title('RS clls: delta r of spike phases - light ON-OFF')
xlabel('frequency bands')
ylabel('cll number - sorted by depth')

subplot(2,2,2)
imagesc(3:5:98,fsd,allrl1(pfs(fsi),:)-allrl0(pfs(fsi),:));
caxis([-.2,.2])
colorbar
title('FS clls: delta r of spike phases - light ON-OFF')
xlabel('frequency bands')
ylabel('cll number - sorted by depth')

subplot(2,2,3)
errorbar(3:5:98,nanmean(allrl1(prs,:)-allrl0(prs,:)),nanstd(allrl1(prs,:)-allrl0(prs,:))./sqrt(size(allrl0(prs,:),1)),'.')
hold on
errorbar(3:5:98,nanmean(allrl1(pfs,:)-allrl0(pfs,:)),nanstd(allrl1(pfs,:)-allrl0(pfs,:))./sqrt(size(allrl0(pfs,:),1)),'r.')
plot(3:5:98,nanmean(allrl1(prs,:)-allrl0(prs,:)),'linewidth',2)
plot(3:5:98,nanmean(allrl1(pfs,:)-allrl0(pfs,:)),'r','linewidth',2)
axis([-1,100,-.04,.04])
line([-1,100],[0,0],'color','k')


figure
subplot(2,2,1)
imagesc(3:5:98,sd,uncircle(allcmeanl0(prs(rsi),:)));
caxis([0,2*pi])
colormap hsv
colorbar
title('RS: average spike phases - light OFF')
xlabel('frequency bands')
ylabel('cll number - sorted by depth')

subplot(2,2,2)
imagesc(3:5:98,sd,uncircle(allcmeanl1(prs(rsi),:)));
caxis([0,2*pi])
colormap hsv
colorbar
title('RS: average spike phases - light ON')
xlabel('frequency bands')
ylabel('cll number - sorted by depth')

subplot(2,2,3)
imagesc(3:5:98,sd,uncircle(allcmeanl0(pfs(fsi),:)));
caxis([0,2*pi])
colormap hsv
colorbar
title('FS: average spike phases - light OFF')
xlabel('frequency bands')
ylabel('cll number - sorted by depth')

subplot(2,2,4)
imagesc(3:5:98,sd,uncircle(allcmeanl1(pfs(fsi),:)));
caxis([0,2*pi])
colormap hsv
colorbar
title('FS: average spike phases - light ON')
xlabel('frequency bands')
ylabel('cll number - sorted by depth')

%phase advance

figure
subplot(2,2,1)
errorbar(3:5:98,mean(uncircle(allcmeanl0(pfs(l4fs),:))),std(uncircle(allcmeanl0(pfs(l4fs),:)))./sqrt(length(l4fs)),'r')
hold on
errorbar(3:5:98,mean(uncircle(allcmeanl0(prs(l4rs),:))),std(uncircle(allcmeanl0(prs(l4rs),:)))./sqrt(length(l4rs)),'b')
line([0,100],[pi,pi],'color','k')
legend([{'FS'},{'RS'}]);
xlabel('frequency [Hz]')
ylabel('average phase of spike');
title('Light OFF: Layer 4');

subplot(2,2,2)
errorbar(3:5:98,mean(uncircle(allcmeanl0(pfs(l23fs),:))),std(uncircle(allcmeanl0(pfs(l23fs),:)))./sqrt(length(l23fs)),'r')
hold on
errorbar(3:5:98,mean(uncircle(allcmeanl0(prs(l23rs),:))),std(uncircle(allcmeanl0(prs(l23rs),:)))./sqrt(length(l23rs)),'b')
line([0,100],[pi,pi],'color','k')
legend([{'FS'},{'RS'}]);
xlabel('frequency [Hz]')
ylabel('average phase of spike');
title('Light OFF: Layer 2/3');

subplot(2,2,3)
errorbar(3:5:98,mean(uncircle(allcmeanl0(pfs(l5fs),:))),std(uncircle(allcmeanl0(pfs(l5fs),:)))./sqrt(length(l5fs)),'r')
hold on
errorbar(3:5:98,mean(uncircle(allcmeanl0(prs(l5rs),:))),std(uncircle(allcmeanl0(prs(l5rs),:)))./sqrt(length(l5rs)),'b')
line([0,100],[pi,pi],'color','k')
legend([{'FS'},{'RS'}]);
xlabel('frequency [Hz]')
ylabel('average phase of spike');
title('Light OFF: Layer 5');

figure

subplot(2,2,1)
errorbar(3:5:98,mean(uncircle(allcmeanl1(pfs(l4fs),:))),std(uncircle(allcmeanl1(pfs(l4fs),:)))./sqrt(length(l4fs)),'r')
hold on
errorbar(3:5:98,mean(uncircle(allcmeanl1(prs(l4rs),:))),std(uncircle(allcmeanl1(prs(l4rs),:)))./sqrt(length(l4rs)),'b')
line([0,100],[pi,pi],'color','k')
legend([{'FS'},{'RS'}]);
xlabel('frequency [Hz]')
ylabel('average phase of spike');
title('Light ON: Layer 4');

subplot(2,2,2)
errorbar(3:5:98,mean(uncircle(allcmeanl1(pfs(l23fs),:))),std(uncircle(allcmeanl1(pfs(l23fs),:)))./sqrt(length(l23fs)),'r')
hold on
errorbar(3:5:98,mean(uncircle(allcmeanl1(prs(l23rs),:))),std(uncircle(allcmeanl1(prs(l23rs),:)))./sqrt(length(l23rs)),'b')
line([0,100],[pi,pi],'color','k')
legend([{'FS'},{'RS'}]);
xlabel('frequency [Hz]')
ylabel('average phase of spike');
title('Light ON: Layer 2/3');

subplot(2,2,3)
errorbar(3:5:98,mean(uncircle(allcmeanl1(pfs(l5fs),:))),std(uncircle(allcmeanl1(pfs(l5fs),:)))./sqrt(length(l5fs)),'r')
hold on
errorbar(3:5:98,mean(uncircle(allcmeanl1(prs(l5rs),:))),std(uncircle(allcmeanl1(prs(l5rs),:)))./sqrt(length(l5rs)),'b')
line([0,100],[pi,pi],'color','k')
legend([{'FS'},{'RS'}]);
xlabel('frequency [Hz]')
ylabel('average phase of spike');
title('Light ON: Layer 5');

figure
subplot(2,2,1)
errorbar(3:5:98,mean(uncircle(allcmeanl0(prs(l4rs),:))),std(uncircle(allcmeanl0(prs(l4rs),:)))./sqrt(length(l4rs)),'b')
hold on
errorbar(3:5:98,mean(uncircle(allcmeanl1(prs(l4rs),:))),std(uncircle(allcmeanl1(prs(l4rs),:)))./sqrt(length(l4rs)),'r')
line([0,100],[pi,pi],'color','k')
legend([{'L0'},{'L1'}]);
xlabel('frequency [Hz]')
ylabel('average phase of spike');
title('RS: Layer 4');

subplot(2,2,2)
errorbar(3:5:98,mean(uncircle(allcmeanl0(prs(l23rs),:))),std(uncircle(allcmeanl0(prs(l23rs),:)))./sqrt(length(l23rs)),'b')
hold on
errorbar(3:5:98,mean(uncircle(allcmeanl1(prs(l23rs),:))),std(uncircle(allcmeanl1(prs(l23rs),:)))./sqrt(length(l23rs)),'r')
line([0,100],[pi,pi],'color','k')
legend([{'L0'},{'L1'}]);
xlabel('frequency [Hz]')
ylabel('average phase of spike');
title('RS: Layer 2/3');

subplot(2,2,3)
errorbar(3:5:98,mean(uncircle(allcmeanl0(prs(l5rs),:))),std(uncircle(allcmeanl0(prs(l5rs),:)))./sqrt(length(l5rs)),'b')
hold on
errorbar(3:5:98,mean(uncircle(allcmeanl1(prs(l5rs),:))),std(uncircle(allcmeanl1(prs(l5rs),:)))./sqrt(length(l5rs)),'r')
line([0,100],[pi,pi],'color','k')
legend([{'L0'},{'L1'}]);
xlabel('frequency [Hz]')
ylabel('average phase of spike');
title('RS: Layer 5');

figure
subplot(2,2,1)
errorbar(3:5:98,mean(uncircle(allcmeanl0(pfs(l4fs),:))),std(uncircle(allcmeanl0(pfs(l4fs),:)))./sqrt(length(l4fs)),'b')
hold on
errorbar(3:5:98,mean(uncircle(allcmeanl1(pfs(l4fs),:))),std(uncircle(allcmeanl1(pfs(l4fs),:)))./sqrt(length(l4fs)),'r')
line([0,100],[pi,pi],'color','k')
legend([{'L0'},{'L1'}]);
xlabel('frequency [Hz]')
ylabel('average phase of spike');
title('FS: Layer 4');

subplot(2,2,2)
errorbar(3:5:98,mean(uncircle(allcmeanl0(pfs(l23fs),:))),std(uncircle(allcmeanl0(pfs(l23fs),:)))./sqrt(length(l23fs)),'b')
hold on
errorbar(3:5:98,mean(uncircle(allcmeanl1(pfs(l23fs),:))),std(uncircle(allcmeanl1(pfs(l23fs),:)))./sqrt(length(l23fs)),'r')
line([0,100],[pi,pi],'color','k')
legend([{'L0'},{'L1'}]);
xlabel('frequency [Hz]')
ylabel('average phase of spike');
title('FS: Layer 2/3');

subplot(2,2,3)
errorbar(3:5:98,mean(uncircle(allcmeanl0(pfs(l5fs),:))),std(uncircle(allcmeanl0(pfs(l5fs),:)))./sqrt(length(l5fs)),'b')
hold on
errorbar(3:5:98,mean(uncircle(allcmeanl1(pfs(l5fs),:))),std(uncircle(allcmeanl1(pfs(l5fs),:)))./sqrt(length(l5fs)),'r')
line([0,100],[pi,pi],'color','k')
legend([{'L0'},{'L1'}]);
xlabel('frequency [Hz]')
ylabel('average phase of spike');
title('FS: Layer 5');

figure
subplot(2,2,1)
imagesc(3:5:98,sd,circ_dist(allcmeanl1(pfs(fsi),:),allcmeanl0(pfs(fsi),:)));
caxis([-pi,pi])
colorbar
title('FS clls: delta spike phases - light ON-OFF')
xlabel('frequency bands')
ylabel('cll number - sorted by depth')

subplot(2,2,2)
imagesc(3:5:98,sd,circ_dist(allcmeanl1(prs(rsi),:),allcmeanl0(prs(rsi),:)));
caxis([-pi,pi])
colorbar
title('RS clls: delta spike phases - light ON-OFF')
xlabel('frequency bands')
ylabel('cll number - sorted by depth')

subplot(2,2,3)
errorbar(3:5:98,circ_mean(circ_dist(allcmeanl1(prs,:),allcmeanl0(prs,:))),circ_std(circ_dist(allcmeanl1(prs,:),allcmeanl0(prs,:)))./sqrt(size(allcmeanl0(prs,:),1)),'.')
hold on
errorbar(3:5:98,circ_mean(circ_dist(allcmeanl1(pfs,:),allcmeanl0(pfs,:))),circ_std(circ_dist(allcmeanl1(pfs,:),allcmeanl0(pfs,:)))./sqrt(size(allcmeanl0(pfs,:),1)),'r.')
plot(3:5:98,circ_mean(circ_dist(allcmeanl1(prs,:),allcmeanl0(prs,:))),'linewidth',2)
plot(3:5:98,circ_mean(circ_dist(allcmeanl1(pfs,:),allcmeanl0(pfs,:))),'r','linewidth',2)
axis([-1,100,-.2,.2])
line([-1,100],[0,0],'color','k')


figure
% plot(nlc50,depth,'k.')
% hold on
plot(nlc50(pfs),depth(pfs),'r.')
axis ij
axis([-.1,1.1,150,1050])
title('Fast spiking C50 in depth')
xlabel('C50')
ylabel('cortical depth [mum]')

figure
plot(nlc50(prs),depth(prs),'b.')
axis ij
axis([-.1,1.1,150,1050])
title('Regular spiking C50 in depth')
xlabel('C50')
ylabel('cortical depth [mum]')

figure
% plot((lc50-nlc50)./(lc50+nlc50),depth,'k.')
% hold on
plot((lc50(pfs)-nlc50(pfs))./(lc50(pfs)+nlc50(pfs)),depth(pfs),'ro','markersize',4,'markerfacecolor','r')
hold on
[x,y,xerr] = runningMedian(depth(pfs),(lc50(pfs)-nlc50(pfs))./(lc50(pfs)+nlc50(pfs)),0,12);
plot(x,y,'r','linewidth',2);
plot(x+xerr,y,'r')
plot(x-xerr,y,'r')
line([-1,1,],[375,375],'color','k','linestyle',':')
line([-1,1,],[500,500],'color','k','linestyle',':')
axis ij
line([0,0],[150,1050],'color','k')
axis([-1.1,1.1,150,1050])
title('Fast spiking C50 changes in depth')
xlabel('(C50 light - C50 no light)/(C50 light + C50 no light)')
ylabel('cortical depth [mum]')

figure
plot((lc50(prs)-nlc50(prs))./(lc50(prs)+nlc50(prs)),depth(prs),'bo','markersize',4,'markerfacecolor','b')
hold on
plot((lc50(prsv&phe)-nlc50(prsv&phe))./(lc50(prsv&phe)+nlc50(prsv&phe)),depth(prsv&phe),'mo','markersize',4,'markerfacecolor','m')
[x,y,xerr] = runningMedian(depth(prsv&~phe),(lc50(prsv&~phe)-nlc50(prsv&~phe))./(lc50(prsv&~phe)+nlc50(prsv&~phe)),0,12);
plot(x,y,'b','linewidth',2);
plot(x+xerr,y,'b')
plot(x-xerr,y,'b')
line([-1,1,],[375,375],'color','k','linestyle',':')
line([-1,1,],[500,500],'color','k','linestyle',':')
axis ij
line([0,0],[150,1050],'color','k')
axis([-1.1,1.1,150,1050])
title('Regular spiking C50 changes in depth')
xlabel('(C50 light - C50 no light)/(C50 light + C50 no light)')
ylabel('cortical depth [mum]')

figure
% plot(nlc50,lc50,'k.')
% hold on
plot(nlc50(pfs),lc50(pfs),'r.')
axis([-.1,1.1,-.1,1.1])
line([-.1,1.1],[-.1,1.1],'color','k')
axis square
[s,p] = ttest(nlc50,lc50);
xlabel('c50 no light')
ylabel('c50 light on')
title(['fast spiking C50 changes with light. p: ' num2str(p)])

figure
plot(nlc50(prs),lc50(prs),'b.')
axis([-.1,1.1,-.1,1.1])
line([-.1,1.1],[-.1,1.1],'color','k')
axis square
[s,p] = ttest(nlc50,lc50);
xlabel('c50 no light')
ylabel('c50 light on')
title(['Regular spiking C50 changes with light. p: ' num2str(p)])

% nlr0(find(nlr0<0)) = 0;
% lr0(find(lr0<0)) = 0;

nlr0(find(nlr0<0)) = 0; 
lr0(find(lr0<0)) = 0;

figure
% plot(nlr0,depth,'k.')
% hold on
plot(nlr0(pfs),depth(pfs),'r.')
axis ij
% axis([-15,15,150,1050])
% line([0,0],[150,1050],'color','k')
title('FS: r0 in depth')
xlabel('r0')
ylabel('cortical depth [mum]')

figure
plot(nlr0(prs),depth(prs),'b.')
axis ij
% axis([-15,15,150,1050])
% line([0,0],[150,1050],'color','k')
title('RS: r0 in depth')
xlabel('r0')
ylabel('cortical depth [mum]')

figure
% plot((lr0-nlr0)./(lr0+nlr0),depth,'k.')
% hold on
plot((lr0(pfs)-nlr0(pfs))./(lr0(pfs)+nlr0(pfs)),depth(pfs),'ro','markersize',4,'markerfacecolor','r')
hold on
[x,y,xerr] = runningMedian(depth(pfs),(lr0(pfs)-nlr0(pfs))./(lr0(pfs)+nlr0(pfs)),0,12);
plot(x,y,'r','linewidth',2);
plot(x+xerr,y,'r')
plot(x-xerr,y,'r')
line([-1,1,],[375,375],'color','k','linestyle',':')
line([-1,1,],[500,500],'color','k','linestyle',':')
axis ij
axis([-1.1,1.1,150,1050])
line([0,0],[150,1050],'color','k')
title('FS: r0 changes in depth')
xlabel('(r0 light - r0 no light)/(r0 light + r0 no light)')
ylabel('cortical depth [mum]')

figure
plot((lr0(prs)-nlr0(prs))./(lr0(prs)+nlr0(prs)),depth(prs),'bo','markersize',4,'markerfacecolor','b')
hold on
plot((lr0(prsv&phe)-nlr0(prsv&phe))./(lr0(prsv&phe)+nlr0(prsv&phe)),depth(prsv&phe),'mo','markersize',4,'markerfacecolor','m')
[x,y,xerr] = runningMedian(depth(prsv&~phe),(lr0(prsv&~phe)-nlr0(prsv&~phe))./(lr0(prsv&~phe)+nlr0(prsv&~phe)),0,12);
plot(x,y,'b','linewidth',2);
plot(x+xerr,y,'b')
plot(x-xerr,y,'b')
line([-1,1,],[375,375],'color','k','linestyle',':')
line([-1,1,],[500,500],'color','k','linestyle',':')
axis ij
axis([-1.1,1.1,150,1050])
line([0,0],[150,1050],'color','k')
title('RS: r0 changes in depth')
xlabel('(r0 light - r0 no light)/(r0 light + r0 no light)')
ylabel('cortical depth [mum]')

figure
% plot(nlr0,lr0,'k.')
% hold on
plot(nlr0(pfs),lr0(pfs),'ro','markersize',4,'markerfacecolor','r')
line([0,30],[0,30],'color','k')
axis square
[s,p] = ttest(nlr0(pfs),lr0(pfs));
xlabel('r0 no light')
ylabel('r0 light on')
title(['FS: r0 changes with light. p: ' num2str(p)])

figure
plot(nlr0(prs),lr0(prs),'b.')
line([0,30],[0,30],'color','k')
axis square
[s,p] = ttest(nlr0(prs),lr0(prs));
xlabel('r0 no light')
ylabel('r0 light on')
title(['RS: r0 changes with light. p: ' num2str(p)])

figure
% plot(nlrmax,depth,'k.')
% hold on
plot(nlrmax(pfs),depth(pfs),'r.')
axis ij
axis([-.1,30,150,1050])
title('FS: rmax in depth')
xlabel('rmax')
ylabel('cortical depth [mum]')

figure
plot(nlrmax(prs),depth(prs),'b.')
axis ij
axis([-.1,30,150,1050])
title('RS: rmax in depth')
xlabel('rmax')
ylabel('cortical depth [mum]')

figure
% plot((lrmax-nlrmax)./(lrmax+nlrmax),depth,'k.')
% hold on
plot((lrmax(pfs)-nlrmax(pfs))./(lrmax(pfs)+nlrmax(pfs)),depth(pfs),'ro','markersize',4,'markerfacecolor','r')
hold on
[x,y,xerr] = runningMedian(depth(pfs),(lrmax(pfs)-nlrmax(pfs))./(lrmax(pfs)+nlrmax(pfs)),0,12);
plot(x,y,'r','linewidth',2);
plot(x+xerr,y,'r')
plot(x-xerr,y,'r')
line([-1,1,],[375,375],'color','k','linestyle',':')
line([-1,1,],[500,500],'color','k','linestyle',':')
axis ij
line([0,0],[150,1050],'color','k')
axis([-1.1,1.1,150,1050])
title('FS: rmax changes in depth')
xlabel('(rmax light - rmax no light)/(rmax light + rmax no light)')
ylabel('cortical depth [mum]')

figure
plot((lrmax(prs)-nlrmax(prs))./(lrmax(prs)+nlrmax(prs)),depth(prs),'bo','markersize',4,'markerfacecolor','b')
hold on
plot((lrmax(prsv&phe)-nlrmax(prsv&phe))./(lrmax(prsv&phe)+nlrmax(prsv&phe)),depth(prsv&phe),'mo','markersize',4,'markerfacecolor','m')
[x,y,xerr] = runningMedian(depth(prsv&~phe),(lrmax(prsv&~phe)-nlrmax(prsv&~phe))./(lrmax(prsv&~phe)+nlrmax(prsv&~phe)),0,12);
plot(x,y,'b','linewidth',2);
plot(x+xerr,y,'b')
plot(x-xerr,y,'b')
line([-1,1,],[375,375],'color','k','linestyle',':')
line([-1,1,],[500,500],'color','k','linestyle',':')
axis ij
line([0,0],[150,1050],'color','k')
axis([-1.1,1.1,150,1050])
title('RS: rmax changes in depth')
xlabel('(rmax light - rmax no light)/(rmax light + rmax no light)')
ylabel('cortical depth [mum]')

figure
% plot(nlrmax,lrmax,'k.')
% hold on
plot(nlrmax(pfs),lrmax(pfs),'r.')
line([0,80],[0,80],'color','k')
axis square
axis([0,80,0,80])
[s,p] = ttest(nlrmax(pfs),lrmax(pfs));
xlabel('rmax no light')
ylabel('rmax light on')
title(['rmax changes with light. p: ' num2str(p)])

figure
plot(nlrmax(prs),lrmax(prs),'b.')
line([0,80],[0,80],'color','k')
axis square
axis([0,80,0,80])
[s,p] = ttest(nlrmax(prs),lrmax(prs));
xlabel('rmax no light')
ylabel('rmax light on')
title(['rmax changes with light. p: ' num2str(p)])

figure
plot(nloriprefratio(prs),depth(prs),'bo','markersize',4,'markerfacecolor','b')
hold on
[x,y,xerr] = runningMedian(depth(prs),nloriprefratio(prs),0,12);
plot(x,y,'b','linewidth',2);
plot(x+xerr,y,'b')
plot(x-xerr,y,'b')
line([0,1],[375,375],'color','k','linestyle',':')
line([0,1],[500,500],'color','k','linestyle',':')
axis ij
axis([0,1.1,150,1050])
title('RS: OSI in depth')
xlabel('OSI')
ylabel('cortical depth [mum]')

figure
plot(nloriprefratio(pfs),depth(pfs),'ro','markersize',4,'markerfacecolor','r')
hold on
[x,y,xerr] = runningMedian(depth(pfs),nloriprefratio(pfs),0,12);
plot(x,y,'r','linewidth',2);
plot(x+xerr,y,'r')
plot(x-xerr,y,'r')
line([0,1],[375,375],'color','k','linestyle',':')
line([0,1],[500,500],'color','k','linestyle',':')
axis ij
axis([-.1,1.1,150,1050])
title('FS: OSI in depth')
xlabel('OSI')
ylabel('cortical depth [mum]')

figure
plot(nldirprefratio(prs),depth(prs),'bo','markersize',4,'markerfacecolor','b')
hold on
[x,y,xerr] = runningMedian(depth(prs),nldirprefratio(prs),0,12);
plot(x,y,'b','linewidth',2);
plot(x+xerr,y,'b')
plot(x-xerr,y,'b')
line([0,1],[375,375],'color','k','linestyle',':')
line([0,1],[500,500],'color','k','linestyle',':')
axis ij
axis([0,1.1,150,1050])
title('RS: DSI in depth')
xlabel('DSI')
ylabel('cortical depth [mum]')

figure
plot(nldirprefratio(pfs),depth(pfs),'ro','markersize',4,'markerfacecolor','r')
hold on
[x,y,xerr] = runningMedian(depth(pfs),nldirprefratio(pfs),0,12);
plot(x,y,'r','linewidth',2);
plot(x+xerr,y,'r')
plot(x-xerr,y,'r')
line([0,1],[375,375],'color','k','linestyle',':')
line([0,1],[500,500],'color','k','linestyle',':')
axis ij
axis([-.1,1.1,150,1050])
title('FS: DSI in depth')
xlabel('DSI')
ylabel('cortical depth [mum]')

figure
plot((loriprefratio(prs)-nloriprefratio(prs))./(loriprefratio(prs)+nloriprefratio(prs)),depth(prs),'bo','markersize',4,'markerfacecolor','b')
hold on
[x,y,xerr] = runningMedian(depth(prs),(loriprefratio(prs)-nloriprefratio(prs))./(loriprefratio(prs)+nloriprefratio(prs)),0,12);
plot(x,y,'b','linewidth',2);
plot(x+xerr,y,'b')
plot(x-xerr,y,'b')
line([-1,1,],[375,375],'color','k','linestyle',':')
line([-1,1,],[500,500],'color','k','linestyle',':')
axis ij
axis([-1.1,1.1,150,1050])
line([0,0],[150,1050],'color','k')
title('RS: OSI changes in depth')
xlabel('(OSI light - OSI no light)/(OSI light + OSI no light)')
ylabel('cortical depth [mum]')

figure
plot((loriprefratio(pfs)-nloriprefratio(pfs))./(loriprefratio(pfs)+nloriprefratio(pfs)),depth(pfs),'ro','markersize',4,'markerfacecolor','r')
hold on
[x,y,xerr] = runningMedian(depth(pfs),(loriprefratio(pfs)-nloriprefratio(pfs))./(loriprefratio(pfs)+nloriprefratio(pfs)),0,12);
plot(x,y,'r','linewidth',2);
plot(x+xerr,y,'r')
plot(x-xerr,y,'r')
line([-1,1,],[375,375],'color','k','linestyle',':')
line([-1,1,],[500,500],'color','k','linestyle',':')
axis ij
axis([-1.1,1.1,150,1050])
line([0,0],[150,1050],'color','k')
title('FS: OSI changes in depth')
xlabel('(OSI light - OSI no light)/(OSI light + OSI no light)')
ylabel('cortical depth [mum]')

figure
plot((ldirprefratio(prs)-nldirprefratio(prs))./(ldirprefratio(prs)+nldirprefratio(prs)),depth(prs),'b.')
hold on
[x,y,xerr] = runningMedian(depth(prs),(ldirprefratio(prs)-nldirprefratio(prs))./(ldirprefratio(prs)+nldirprefratio(prs)),0,12);
plot(x,y,'b','linewidth',2);
plot(x+xerr,y,'b')
plot(x-xerr,y,'b')
line([-1,1,],[375,375],'color','k','linestyle',':')
line([-1,1,],[500,500],'color','k','linestyle',':')
axis ij
axis([-1.1,1.1,150,1050])
line([0,0],[150,1050],'color','k')
title('RS: dsi changes in depth')
xlabel('(dsi light - dsi no light)/(dsi light + dsi no light)')
ylabel('cortical depth [mum]')

figure
plot((ldirprefratio(pfs)-nldirprefratio(pfs))./(ldirprefratio(pfs)+nldirprefratio(pfs)),depth(pfs),'r.')
hold on
[x,y,xerr] = runningMedian(depth(pfs),(ldirprefratio(pfs)-nldirprefratio(pfs))./(ldirprefratio(pfs)+nldirprefratio(pfs)),0,12);
plot(x,y,'r','linewidth',2);
plot(x+xerr,y,'r')
plot(x-xerr,y,'r')
line([-1,1,],[375,375],'color','k','linestyle',':')
line([-1,1,],[500,500],'color','k','linestyle',':')
axis ij
axis([-1.1,1.1,150,1050])
line([0,0],[150,1050],'color','k')
title('FS: dsi changes in depth')
xlabel('(dsi light - dsi no light)/(dsi light + DSI no light)')
ylabel('cortical depth [mum]')

figure
[s,p] = ttest(nloriprefratio(prs),loriprefratio(prs));
plot(nloriprefratio(prs),loriprefratio(prs),'b.')
line([0,1],[0,1],'color','k')
axis square
xlabel('OSI control');
ylabel('OSI light');
title(['RS clls: p: ' num2str(p)])

figure
[s,p] = ttest(nloriprefratio(pfs),loriprefratio(pfs));
plot(nloriprefratio(pfs),loriprefratio(pfs),'r.')
line([0,1],[0,1],'color','k')
axis square
xlabel('OSI control');
ylabel('OSI light');
title(['FS clls p: ' num2str(p)])

figure
[s,p] = ttest(nldirprefratio(prs),ldirprefratio(prs));
plot(nldirprefratio(prs),ldirprefratio(prs),'b.')
line([0,1],[0,1],'color','k')
axis square
xlabel('DSI control');
ylabel('DSI light');
title(['RS clls p: ' num2str(p)])

figure
[s,p] = ttest(nldirprefratio(pfs),ldirprefratio(pfs));
plot(nldirprefratio(pfs),ldirprefratio(pfs),'r.')
line([0,1],[0,1],'color','k')
axis square
xlabel('DSI control');
ylabel('DSI light');
title(['FS clls p: ' num2str(p)])

%%%%%%%%%%%%%%%%%%%%%%%%%% baseline subtracted preffr %%%%%%%%%%%%%%%%%%%%%%%%%% 
blspreffr(:,1) = preffrl0-bfr'; blspreffr(:,2) = preffrl1-bfr'; 
blspreffr(find(blspreffr(:,2)<0),2) = 0;

figure
plot(blspreffr(prs,1),depth(prs),'b.')
axis ij
% axis([-.1,1.1,150,1050])
title('RS: BLS Preferred FR in depth')
xlabel('Firing Rate [Hz]')
ylabel('cortical depth [mum]')

figure
plot(blspreffr(pfs,2),depth(pfs),'r.')
axis ij
% axis([-.1,1.1,150,1050])
title('FS: BLS Preferred FR in depth')
xlabel('Firing Rate [Hz]')
ylabel('cortical depth [mum]')

figure
plot((blspreffr(prs,2)-blspreffr(prs,1))./(blspreffr(prs,2)+blspreffr(prs,1)),depth(prs),'bo','markersize',4,'markerfacecolor','b')
hold on
plot((blspreffr(prsv&phe,2)-blspreffr(prsv&phe,1))./(blspreffr(prsv&phe,2)+blspreffr(prsv&phe,1)),depth(prsv&phe),'mo','markersize',4,'markerfacecolor','m')
[x,y,xerr] = runningMedian(depth(prsv&~phe),(blspreffr(prsv&~phe,2)-blspreffr(prsv&~phe,1))./(blspreffr(prsv&~phe,2)+blspreffr(prsv&~phe,1)),0,12);
plot(x,y,'b','linewidth',2);
plot(x+xerr,y,'b')
plot(x-xerr,y,'b')
line([-1,1,],[375,375],'color','k','linestyle',':')
line([-1,1,],[500,500],'color','k','linestyle',':')
line([0,0],[150,950],'color','k')
axis ij
axis([-1.1,1.1,150,950])
ylabel('depth[mum]')
xlabel('(light-nolight)/(light+nolight)')
title('RS clls: BLS preferred firing rate changes by layer 4 suppression')

figure
plot((blspreffr(pfs,2)-blspreffr(pfs,1))./(blspreffr(pfs,2)+blspreffr(pfs,1)),depth(pfs),'ro','markersize',4,'markerfacecolor','r')
hold on
[x,y,xerr] = runningMedian(depth(pfs),(blspreffr(pfs,2)-blspreffr(pfs,1))./(blspreffr(pfs,2)+blspreffr(pfs,1)),0,12);
plot(x,y,'r','linewidth',2);
plot(x+xerr,y,'r')
plot(x-xerr,y,'r')
line([-1,1,],[375,375],'color','k','linestyle',':')
line([-1,1,],[500,500],'color','k','linestyle',':')
line([0,0],[150,950],'color','k')
axis ij
axis([-1.1,1.1,150,950])
ylabel('depth[mum]')
xlabel('(light-nolight)/(light+nolight)')
title('FS clls: BLS preferred firing rate changes by layer 4 suppression')

figure
plot(blspreffr(prs,1),blspreffr(prs,2),'b.')
line([0,50],[0,50],'color','k')
axis square
xlabel('preferred firing rate no light')
ylabel('preferred firing rate light on')
title('RS clls: BLS pref firing rate changes');

figure
plot(blspreffr(pfs,1),blspreffr(pfs,2),'r.')
line([0,80],[0,80],'color','k')
axis square
xlabel('preferred firing rate no light')
ylabel('preferred firing rate light on')
title('FS clls: BLS pref firing rate changes');

%%%%%%%%%%%%%%%%%%%%%%%% preferred firing rates %%%%%%%%%%%%%%%%%%%%%%
figure
plot(preffr(prs,1),depth(prs),'b.')
axis ij
% axis([-.1,1.1,150,1050])
title('RS: Preferred FR in depth')
xlabel('Firing Rate [Hz]')
ylabel('cortical depth [mum]')

figure
plot(preffr(pfs,2),depth(pfs),'r.')
axis ij
% axis([-.1,1.1,150,1050])
title('FS: Preferred FR in depth')
xlabel('Firing Rate [Hz]')
ylabel('cortical depth [mum]')

figure
plot((preffr(prs,2)-preffr(prs,1))./(preffr(prs,2)+preffr(prs,1)),depth(prs),'bo','markersize',4,'markerfacecolor','b')
hold on
plot((preffr(prsv&phe,2)-preffr(prsv&phe,1))./(preffr(prsv&phe,2)+preffr(prsv&phe,1)),depth(prsv&phe),'mo','markersize',4,'markerfacecolor','m')
[x,y,xerr] = runningMedian(depth(prsv&~phe),(preffr(prsv&~phe,2)-preffr(prsv&~phe,1))./(preffr(prsv&~phe,2)+preffr(prsv&~phe,1)),0,12);
plot(x,y,'b','linewidth',2);
plot(x+xerr,y,'b')
plot(x-xerr,y,'b')
line([-1,1,],[375,375],'color','k','linestyle',':')
line([-1,1,],[500,500],'color','k','linestyle',':')
line([0,0],[150,950],'color','k')
axis ij
axis([-1.1,1.1,0,1000])
ylabel('depth[mum]')
xlabel('(light-nolight)/(light+nolight)')
title('RS clls: preferred firing rate changes by layer 4 suppression')

figure
plot((preffr(pfs,2)-preffr(pfs,1))./(preffr(pfs,2)+preffr(pfs,1)),depth(pfs),'ro','markersize',4,'markerfacecolor','r')
hold on
[x,y,xerr] = runningMedian(depth(pfs),(preffr(pfs,2)-preffr(pfs,1))./(preffr(pfs,2)+preffr(pfs,1)),0,12);
plot(x,y,'r','linewidth',2);
plot(x+xerr,y,'r')
plot(x-xerr,y,'r')
line([-1,1,],[375,375],'color','k','linestyle',':')
line([-1,1,],[500,500],'color','k','linestyle',':')
line([0,0],[150,950],'color','k')
axis ij
axis([-1.1,1.1,0,1000])
ylabel('depth[mum]')
xlabel('(light-nolight)/(light+nolight)')
title('FS clls: preferred firing rate changes by layer 4 suppression')

figure
plot(preffr(prs,1),preffr(prs,2),'b.')
line([0,50],[0,50],'color','k')
axis square
xlabel('preferred firing rate no light')
ylabel('preferred firing rate light on')
title('RS clls: pref firing rate changes');

figure
plot(preffr(pfs,1),preffr(pfs,2),'r.')
line([0,80],[0,80],'color','k')
axis square
xlabel('preferred firing rate no light')
ylabel('preferred firing rate light on')
title('FS clls: pref firing rate changes');

%%%%%%%%%%%%%%%%%%%%%%%% average firing rates %%%%%%%%%%%%%%%%%%%%%%
figure
plot(nlfr(prs),depth(prs),'b.')
axis ij
% axis([-.1,1.1,150,1050])
title('RS: Average FR in depth')
xlabel('Firing Rate [Hz]')
ylabel('cortical depth [mum]')

figure
plot(nlfr(pfs),depth(pfs),'r.')
axis ij
% axis([-.1,1.1,150,1050])
title('FS: Average FR in depth')
xlabel('Firing Rate [Hz]')
ylabel('cortical depth [mum]')

figure
plot((lfr(prs)-nlfr(prs))./(lfr(prs)+nlfr(prs)),depth(prs),'bo','markersize',4,'markerfacecolor','b')
hold on
plot((lfr(prsv&phe)-nlfr(prsv&phe))./(lfr(prsv&phe)+nlfr(prsv&phe)),depth(prsv&phe),'mo','markersize',4,'markerfacecolor','m')
[x,y,xerr] = runningMedian(depth(prsv&~phe),(lfr(prsv&~phe)-nlfr(prsv&~phe))./(lfr(prsv&~phe)+nlfr(prsv&~phe)),0,12);
plot(x,y,'b','linewidth',2);
plot(x+xerr,y,'b')
plot(x-xerr,y,'b')
line([-1,1,],[375,375],'color','k','linestyle',':')
line([-1,1,],[500,500],'color','k','linestyle',':')
line([0,0],[150,950],'color','k')
line([0,0],[0,1000],'color','k')
axis ij
axis([-1.1,1.1,150,950])
ylabel('depth[mum]')
xlabel('(light-nolight)/(light+nolight)')
title('RS clls: average firing rate changes by layer 4 suppression')

figure
plot((lfr(pfs)-nlfr(pfs))./(lfr(pfs)+nlfr(pfs)),depth(pfs),'ro','markersize',4,'markerfacecolor','r')
hold on
[x,y,xerr] = runningMedian(depth(pfsv),(lfr(pfsv)-nlfr(pfsv))./(lfr(pfsv)+nlfr(pfsv)),0,12);
plot(x,y,'r','linewidth',2);
plot(x+xerr,y,'r')
plot(x-xerr,y,'r')
line([-1,1,],[375,375],'color','k','linestyle',':')
line([-1,1,],[500,500],'color','k','linestyle',':')
line([0,0],[150,950],'color','k')
line([0,0],[0,1000],'color','k')
axis ij
axis([-1.1,1.1,150,950])
ylabel('depth[mum]')
xlabel('(light-nolight)/(light+nolight)')
title('FS clls: average firing rate changes by layer 4 suppression')

figure
plot(nlfr(prs),lfr(prs),'b.')
line([0,50],[0,50],'color','k')
axis square
xlabel('preferred firing rate no light')
ylabel('preferred firing rate light on')
title('RS clls:average firing rate changes');

figure
plot(nlfr(pfs),lfr(pfs),'r.')
line([0,80],[0,80],'color','k')
axis square
xlabel('preferred firing rate no light')
ylabel('preferred firing rate light on')
title('FS clls: average firing rate changes');

%%%%%%%%%%%%%%%% Control FR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
plot(controlfr(prs,1),depth(prs),'b.')
axis ij
% axis([-.1,1.1,150,1050])
title('RS: Control FR in depth')
xlabel('Firing Rate [Hz]')
ylabel('cortical depth [mum]')

figure
plot(controlfr(pfs,2),depth(pfs),'r.')
axis ij
% axis([-.1,1.1,150,1050])
title('FS: Control FR in depth')
xlabel('Firing Rate [Hz]')
ylabel('cortical depth [mum]')

figure
plot((controlfr(prs,2)-controlfr(prs,1))./(controlfr(prs,2)+controlfr(prs,1)),depth(prs),'bo','markersize',4,'markerfacecolor','b');
hold on
plot((controlfr(prsv&phe,2)-controlfr(prsv&phe,1))./(controlfr(prsv&phe,2)+controlfr(prsv&phe,1)),depth(prsv&phe),'mo','markersize',4,'markerfacecolor','m')
[x,y,xerr] = runningMedian(depth(prsv&~phe),(controlfr(prsv&~phe,2)-controlfr(prsv&~phe,1))./(controlfr(prsv&~phe,2)+controlfr(prsv&~phe,1)),0,12);
plot(x,y,'b','linewidth',2);
plot(x+xerr,y,'b')
plot(x-xerr,y,'b')
line([-1,1,],[375,375],'color','k','linestyle',':')
line([-1,1,],[500,500],'color','k','linestyle',':')
line([0,0],[150,950],'color','k')
axis ij
axis([-1.1,1.1,0,1000])
ylabel('depth[mum]')
xlabel('(light-nolight)/(light+nolight)')
title('RS clls: control firing rate changes by layer 4 suppression')

figure
plot((controlfr(pfs,2)-controlfr(pfs,1))./(controlfr(pfs,2)+controlfr(pfs,1)),depth(pfs),'ro','markersize',4,'markerfacecolor','r');
hold on
[x,y,xerr] = runningMedian(depth(pfs),(controlfr(pfs,2)-controlfr(pfs,1))./(controlfr(pfs,2)+controlfr(pfs,1)),0,12);
plot(x,y,'r','linewidth',2);
plot(x+xerr,y,'r')
plot(x-xerr,y,'r')
line([-1,1,],[375,375],'color','k','linestyle',':')
line([-1,1,],[500,500],'color','k','linestyle',':')
line([0,0],[150,950],'color','k')
line([0,0],[0,1000],'color','k')
axis ij
axis([-1.1,1.1,0,1000])
ylabel('depth[mum]')
xlabel('(light-nolight)/(light+nolight)')
title('FS clls: control firing rate changes by layer 4 suppression')

figure
plot(controlfr(prs,1),controlfr(prs,2),'b.')
line([0,30],[0,30],'color','k')
axis square
xlabel('control firing rate no light')
ylabel('control firing rate light on')
title('RS clls: control firing rate changes');

figure
plot(controlfr(pfs,1),controlfr(pfs,2),'r.')
line([0,30],[0,30],'color','k')
axis square
xlabel('control firing rate no light')
ylabel('control firing rate light on')
title('FS clls: control firing rate changes');

% animals:
c = colormap;
skip = floor(size(c,1)/max(animalno));

figure
hold on
for i = 1:max(animalno)
    aprs = intersect(prs,find(animalno==i)); apfs = intersect(pfs,find(animalno==i));
    plot((blspreffr(aprs,2)-blspreffr(aprs,1))./(blspreffr(aprs,2)+blspreffr(aprs,1)),depth(aprs),'o','linewidth',2,'color',c((i-1)*skip+1,:))
end
line([0,0],[0,1100],'color','k')
axis ij
axis([-1.1,1.1,150,1100])
ylabel('depth[mum]')
xlabel('(light-nolight)/(light+nolight)')
title('RS clls: BLS preferred firing rate changes by layer 4 suppression')
a = 1:max(animalno);
legend(int2str(a'));

figure
hold on
for i = 1:max(animalno)
    aprs = intersect(prs,find(animalno==i)); apfs = intersect(pfs,find(animalno==i));
    plot(nloriprefratio(aprs),loriprefratio(aprs),'o','linewidth',2,'color',c((i-1)*skip+1,:))
end
line([0,1],[0,1],'color','k')
axis square
ylabel('OSI light on')
xlabel('OSI light off')
title('RS clls: OSI changes by lower layer SOM suppression')
a = 1:max(animalno);
legend(int2str(a'));

figure
hold on
for i = 1:max(animalno)
    aprs = intersect(prs,find(animalno==i)); apfs = intersect(pfs,find(animalno==i));
    plot(nldirprefratio(aprs),ldirprefratio(aprs),'o','linewidth',2,'color',c((i-1)*skip+1,:))
end
line([0,1],[0,1],'color','k')
axis square
ylabel('DSI light on')
xlabel('DSI light off')
title('RS clls: OSI changes by lower layer SOM suppression')
a = 1:max(animalno);
legend(int2str(a'));

disp('')
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

function [params, resnorm, rsquared] = fit_crf_NR(x,y)
% [rmax, n, c50, r0]

    nrfitmargin = .5;
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

function depthplot(value,depth,formatstr,formatstr2)
plot(value,depth,formatstr)
axis ij
line([0,0],[0,1200],'color','k')
hold on
[x,y,xerr] = runningAverage(depth,value,0,7);
plot(x,y,formatstr2,'linewidth',2)
plot(x-xerr,y,formatstr2)
plot(x+xerr,y,formatstr2)
end

function out = uncircle(in)
    in(find(in<0)) = in(find(in<0))+2*pi;
    out = in;
end

function [avg,dispersion,amplitude,offset,rmsd] = fitGauss(yvals,xvals)

    % [mu, sigma, amplitude, spontrate]
    p0 = [mean(xvals) 4 max(yvals)-min(yvals) min(yvals)];
    lb = [0,  1,  0    ,0];
    ub = [max(xvals),8,max(yvals),max(yvals)];
    
    if lb(3) == ub(3), ub(3) = ub(3)+1; end
    if lb(4) == ub(4), ub(4) = ub(4)+1; end
    warning off
    options=optimset('TolFun',1e-8,'Display','off');
    [params,resnorm,residual,exitflag] = lsqcurvefit(@(p,x) gauss(p,x),p0,...
        xvals,yvals,lb,ub,options);
    %     [0 0 mean(yvals(:)) 0],...
    %     [180 90 2*max(yvals(:)) mean(yvals(:))],options);
    %warning on
    avg = params(1);
    dispersion = params(2);
    amplitude = params(3);
    offset = params(4);
    rmsd = sqrt(mean(residual.^2));
end

function val = gauss(params,x)
   val = exp(-((x-params(1)).^2)./(2*params(2)^2));
   val=params(4) + params(3)*val(:);
end