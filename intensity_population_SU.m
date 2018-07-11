function intensity_population_SU

% % SOM later population
% animalids = {'150331', '150401','150602','150603','150825','150902', '151110'};
% blocks    = [5,         9,       9,       8,       16,      15,       16];
% animal    = [1,         2,       3,       4,       5,       6,        7];
% electrodes =[[1,32];    [1,32];  [1,32];  [1,16]; [17,32]; [1,16];   [1,16]];
% penangle =  [25,        25,      25,      10,      25,      25,       25];
% spacing =   [25,        25,      25,      25,      25,      25,       25];
% color =     [3,         3,       3,       3,       3,       3,        3];
% popfile = 'C:\Users\Julia\work\data\populations\SOM_Halo_later\intensity\SOM_intensity_SU.mat';

% % PV Halo up to 50% population
% animalids = {'150804', '150804'};
% blocks    = [15,        15];
% animal    = [1,         1];
% electrodes =[[1,16];    [17,32]];
% penangle =  [25,        25];
% spacing =   [25,        25];
% color =     [3,         3];
% popfile = 'C:\Users\Julia\work\data\populations\PV_Halo\intensity\PV50_intensity_SU.mat';

% % PV Halo up to 100% population
% animalids = {'150820', '150820'};
% blocks    = [16,         16];
% animal    = [1,          1];
% electrodes =[[1,16];    [17,32]];
% penangle =  [25,        25];
% spacing =   [25,        25];
% color =     [3,         3];
% popfile = 'C:\Users\Julia\work\data\populations\PV_Halo\intensity\PV100_intensity_SU.mat';
% % 
% % PV eArch population
% animalids = {'151211', '151211'};
% blocks    = [5,         6];
% animal    = [1,         1];
% electrodes =[[1,16];    [1,16]];
% penangle =  [25,        25];
% spacing =   [25,        25];
% color =     [3,         4];
% popfile = 'C:\Users\Julia\work\data\populations\PV_eArch\intensity\PV_intensity_SU.mat';

% % PV Arch reporter population - GERMLINE
% animalids = {'151216', '151216', '151216'};
% blocks    = [5,         6,        7];
% animal    = [1,         1,        1];
% electrodes =[[1,16];    [1,16];  [1,16]];
% penangle =  [25,        25,       25];
% spacing =   [25,        25,       25];
% color =     [4,         3,        4];
% popfile = 'C:\Users\Julia\work\data\populations\PV_Archtransg\intensity\PV_intensity_SU.mat';
% 
% % control population - red
% animalids = {'151214', '151215', '151217'};
% blocks    = [4,         8,        5];
% animal    = [1,         2,        3];
% electrodes =[[1,16];    [1,16];   [1,16]];
% penangle =  [25,        25,       25];
% spacing =   [25,        25,       25];
% color =     [3,         3,        3];
% popfile = 'C:\Users\Julia\work\data\populations\control\intensity\intensityred_SU.mat';

% % control population - green
% animalids = {'151214', '151215', '151215', '151217', '151217'};
% blocks    = [5,        9,        10,        7,        8];
% animal    = [1,        2,        2,         3,        3];
% electrodes =[[1,16];  [1,16];   [1,16];    [1,16];   [1,16]];
% penangle =  [25,       25,       25,        25,       25];
% spacing =   [25,       25,       25,        25,       25];
% color =     [4,        4,        4,         4,        4];
% popfile = 'C:\Users\Julia\work\data\populations\control\intensity\intensitygreen_SU.mat';


% 14, 20 and 40 interpolated for red
intensitiesred =   [5,    9,    14,   16,   20,   28,   40,    50, 10,   18,   32,  56,   100];
muwattsred       = [1.44, 2.69, 4.35, 4.84, 6.05, 8.56, 11.66, 15, 3.05, 5.44, 9.8, 16.7, 27.8];
intensitiesgreen = [2,    4,    6,    8,    10,   18,    32,   56,   100];
muwattsgreen     = [0.33, 1.48, 2.79, 4.18, 5.43, 11.47, 24.6, 50.7, 91.2];


recalculate_pop = 0;
recalculate_muafile = 0;

% chronux parameters
params.tapers = [5,9]; params.Fs = 1000; params.err = [2, 0.05]; params.trialave = 1;


if ~exist(popfile) || recalculate_pop
    
    cll = 1;
    for blck = 1:length(blocks)
        
        
        supath = ['C:\Users\Julia\work\data\' animalids{blck} '\singleunits\'];
        basename = [animalids{blck} '_block' int2str(blocks(blck)) '_tet'];
        
        files = dir([supath, basename, '*.mat']);
        
        prestim = 300;
        poststim = 700;
        respwin = 501:1500; % after stimulus onset
        offsetwin = 1501:2500;
        respwin = respwin+prestim;
        
        
        for fi = 1:length(files)
            
            if strfind(files(fi).name, 'MU')
                continue;
            end
            
            load([supath, files(fi).name]);
            
            % figure out light intensity
            lightlevels = unique(result.light);
            for i = 1:length(lightlevels)
                if lightlevels(i) == 0
                    watts(blck,i) = 0;
                elseif color(blck) == 3
                    a = find(intensitiesred == lightlevels(i));
                    watts(blck,i) = muwattsred(a);
                elseif color(blck) == 4
                    a = find(intensitiesgreen == lightlevels(i));
                    watts(blck,i) = muwattsgreen(a);
                end
            end
            
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
            if ~isempty(msStimes) && msStimes(1) == 0, msStimes(1) = 1; end
            
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
            
            disp(['Block ' int2str(blck) '/' int2str(length(blocks)) '   file ' int2str(fi) '/' int2str(length(electrodes(blck,1):electrodes(blck,2)))]);
            
            trialdur = result.stimduration*1000;
            msstamps = result.msstamps;
            if length(msstamps)~=length(result.light)
                %                 msstamps([193,291]) = []; % for 151023 block 5
                %                 result.msstamps = msstamps;
                %                 save(file,'result');
                pause;
            end
            
            for i = 1:length(msstamps)
                speed(i,:) = result.runspeed(msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim);
            end
            
            % figure out sufficiently high and nonvariable runspeed trials
            meanspeed = mean(speed(:,respwin),2);
            stdspeed = std(speed(:,respwin),1,2);
            notstill = find(meanspeed>1);
            okspeed = find(meanspeed>( mean(meanspeed(notstill))-(1.5*std(meanspeed(notstill))) ) );
            okvar = find(stdspeed<( mean(stdspeed(notstill))+(1.5*std(stdspeed(notstill)))) & stdspeed>.5);
            oktrials = intersect(okspeed,okvar);
            nonoktrials = 1:size(speed,1); nonoktrials(oktrials) = [];
            stilltrials = 1:size(speed,1); stilltrials(notstill) = [];
            
            % find low and high gamma peaks in spectra
            beta = [15,40];
            gamma = [50,70];
            sizes = unique(result.gratingInfo.size);
            large = find(result.gratingInfo.size == max(sizes) & result.light == 0);
            small = find(result.gratingInfo.size == min(sizes) & result.light == 0);
            lr1 = intersect(large,oktrials);
            sr1 = intersect(small,oktrials);
            for i = 1:length(large)
                [pl(i,:),f] = mtspectrumc(lfp(result.msstamps(large(i))+700:result.msstamps(large(i))+1500),params);
                [ps(i,:),f] = mtspectrumc(lfp(result.msstamps(small(i))+700:result.msstamps(small(i))+1500),params);
            end
            plr1 = nan(1,size(pl,2)); psr1 = nan(1,size(pl,2));
            if length(lr1>=5)
                for i = 1:length(lr1)
                    [plr1(i,:),f] = mtspectrumc(lfp(result.msstamps(lr1(i))+700:result.msstamps(lr1(i))+1500),params);
                    %                     [plr1(i,:),f] = pmtm(centresult.lfp(ch,centresult.msstamps(lr1(i))+300:centresult.msstamps(lr1(i))+1300),3,[],1000);
                end
            end
            if length(sr1>=5)
                for i = 1:length(sr1)
                    [psr1(i,:),f] = mtspectrumc(lfp(result.msstamps(sr1(i))+700:result.msstamps(sr1(i))+1500),params);
                    %                     [psr1(i,:),f] = pmtm(centresult.lfp(ch,centresult.msstamps(sr1(i))+300:centresult.msstamps(sr1(i))+1300),3,[],1000);
                end
            end
            b1 = find(f>beta(1),1); b2 = find(f>beta(2),1);
            g1 = find(f>gamma(1),1); g2 = find(f>gamma(2),1);
            if ~isnan(plr1(1)) % take the peak from the running spectra if there is more than a couple of tunning trials
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
                g1peak(blck) = 0;
            else
                peaks = find(diff(bsig)>0)+1;
                pvs = bsig(peaks);
                bpi = peaks(pvs == max(pvs));
                bpi = bpi+b1-1;
                g1peak(blck) = 1;
            end
            if isempty(find(diff(gsig)>0)) % there is no clear beta peak
                gpi = round((g1+g2)/2);
                g2peak(blck) = 0;
            else
                peaks = find(diff(gsig)>0)+1;
                pvs = gsig(peaks);
                gpi = peaks(pvs == max(pvs));
                gpi = gpi+g1-1;
                g2peak(blck) = 1;
            end
            lgi(cll) = bpi; hgi(cll) = gpi;
            
            sr = 1000;
            gbandwidth = 20;
            gamma1 = eegfilt(lfp,sr,f(bpi)-gbandwidth/2,f(bpi)+gbandwidth/2);
            gamma2 = eegfilt(lfp,sr,f(gpi)-gbandwidth/2,f(gpi)+gbandwidth/2);
            h1 = hilbert(gamma1); gpow1 = abs(h1); gphas1 = angle(h1);
            h2 = hilbert(gamma2); gpow2 = abs(h2); gphas2 = angle(h2);
            
            freqbinwidth = 10; clear filtmat; clear phasmat; clear powmat;
            for i = 1:100/freqbinwidth
                filtmat(i,:) = eegfilt(lfp,sr,(i-1)*freqbinwidth+1,i*freqbinwidth);
                h = hilbert(filtmat(i,:));
                powmat(i,:) = abs(h); phasmat(i,:) = angle(h);
            end
                       
            clear gamma1phases; clear gamma2phases; clear gamma1powers; clear gamma2powers;
            clear bandphases; clear bandpowers;
            gamma1phases = zeros(length(msstamps),trialdur+poststim+prestim);
            gamma2phases = zeros(length(msstamps),trialdur+poststim+prestim);
            bandphases = zeros(size(powmat,1),length(msstamps),trialdur+poststim+prestim);
            bandpowers = zeros(size(powmat,1),length(msstamps),trialdur+poststim+prestim);
            for i = 1:length(msstamps)
                resp(i,:) = chan(msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim);
                hh = find(resp(i,1001:1795))'; ptresp(i).times = hh./1000;
                lfpresp(i,:) = lfp(msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim);
                lfpspecttbt(i,:) = mtspectrumc(squeeze(lfpresp(i,1001:1800))',params);
                
                gamma1powresp(i,:) = gpow1(msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim);
                gamma2powresp(i,:) = gpow2(msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim);
                gamma1resp(i,:) = gamma1(msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim);
                gamma2resp(i,:) = gamma2(msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim);
                gamma1phasresp(i,:) = gphas1(msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim);
                gamma2phasresp(i,:) = gphas2(msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim);
                gamma1phases(i,find(resp(i,:))) = gamma1phasresp(i,find(resp(i,:)));
                gamma2phases(i,find(resp(i,:))) = gamma2phasresp(i,find(resp(i,:)));
                gamma1powers(i,find(resp(i,:))) = gamma1powresp(i,find(resp(i,:)));
                gamma2powers(i,find(resp(i,:))) = gamma2powresp(i,find(resp(i,:)));
                
                for j = 1:size(phasmat,1)
                    bandphaseresp(j,i,:) = phasmat(j, msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim);
                    bandpowresp(j,i,:) = powmat(j, msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim);
                    bandphases(j,i,find(resp(i,:))) = bandphaseresp(j,i,find(resp(i,:)));
                    bandpowers(j,i,find(resp(i,:))) = bandpowresp(j,i,find(resp(i,:)));
                end
            end
            
            frs = sum(resp(:,respwin),2)./(length(respwin)/1000);
            bl = sum(resp(:,1:prestim),2)./(prestim/1000);
            g1pow = mean(gamma1powresp(:,respwin),2);
            g2pow = mean(gamma2powresp(:,respwin),2);
            
            lg = []; sg = []; % light group and stim group
            for l = 1:length(lightlevels)
                for d = 1:length(sizes)
                    
                    thisinds = find(result.gratingInfo.size == sizes(d) & result.light == lightlevels(l));
                    thisruninds = intersect(thisinds,oktrials); thisstillinds = intersect(thisinds,stilltrials);
                    
                    condresp(cll,l,d,:) = nanmean(resp(thisinds,:),1);
                    condfr(cll,l,d) = nanmean(frs(thisinds));
                    condg1pow(cll,l,d) = nanmean(g1pow(thisinds));
                    condg2pow(cll,l,d) = nanmean(g2pow(thisinds));
                    
                    condlfpresp(cll,l,d,:) = nanmean(lfpresp(thisinds,:));
                    [lfpspect(cll,l,d,:),chf,lfpspecterr(cll,l,d,:,:)] = mtspectrumc(squeeze(lfpresp(thisinds,1001:1800))',params);
                    tbtlfpspect(cll,l,d,:) = nanmean(lfpspecttbt(thisinds,:),1);
                    tbtlfpspecterr(cll,l,d,:) = nanstd(lfpspecttbt(thisinds,:),1,1)./sqrt(length(thisinds));
                    %                 [sfcoher(cll,l,d,:),a,b,c,de,sfcfx,e,f,g,sfcoherr_c(cll,l,d,:,:)] = coherencycpt(lfpresp(thisinds,1001:1800)',ptresp(thisinds),params);
                    
                    condg1powresp(cll,l,d,:) = mean(gamma1powresp(thisinds,:));
                    condg2powresp(cll,l,d,:) = mean(gamma2powresp(thisinds,:));
                    
                    clear g1phasmatl; clear allg1phases;
                    g1phasmat = gamma1phases(thisinds,respwin);
                    allg1phases = squeeze(g1phasmat(find(g1phasmat)));
                    if isrow(allg1phases), allg1phases = allg1phases'; end
                    condg1phases{cll,l,d} = allg1phases;
                    prefg1phase(cll,l,d) = circ_mean(allg1phases);
                    g1lockpval(cll,l,d) = circ_rtest(allg1phases);
                    g1plv(cll,l,d) = circ_r(allg1phases);
                    g1ppc(cll,l,d) = ppc(allg1phases);
                    [g1angdev(cll,l,d),g1cstd(cll,l,d)] = circ_std(allg1phases);
                    
                    clear g2phasmatl; clear allg2phases;
                    g2phasmat = gamma2phases(thisinds,respwin);
                    allg2phases = squeeze(g2phasmat(find(g2phasmat)));
                    if isrow(allg2phases), allg2phases = allg2phases'; end
                    condg2phases{cll,l,d} = allg2phases;
                    prefg2phase(cll,l,d) = circ_mean(allg2phases);
                    g2lockpval(cll,l,d) = circ_rtest(allg2phases);
                    g2plv(cll,l,d) = circ_r(allg2phases);
                    g2ppc(cll,l,d) = ppc(allg2phases);
                    [g2angdev(cll,l,d),g2cstd(cll,l,d)] = circ_std(allg2phases);
                    
                    clear phasmat; clear allphases;
                    for b = 1:size(bandphases,1)
                        phasmat = squeeze(bandphases(b,thisinds,respwin));
                        allphases = squeeze(phasmat(find(phasmat)));
                        if isrow(allphases), allphases = allphases'; end
                        condbandphases{cll,b,l,d} = allphases;
                        bandprefphase(cll,b,l,d) = circ_mean(allphases);
                        bandlockpval(cll,b,l,d) = circ_rtest(allphases);
                        bandplv(cll,b,l,d) = circ_r(allphases);
                        bandppc(cll,b,l,d) = ppc(allphases);
                        [bandangdev(cll,b,l,d),bandcstd(cll,b,l,d)] = circ_std(allphases);
                        
                        powmat = squeeze(bandpowers(b,thisinds,respwin));
                        allpowers = squeeze(powmat(find(powmat)));
                        if isrow(allpowers), allpowers = allpowers'; end
                        condbandpowers{cll,b,l,d} = allpowers;
                    end
                    
                    gnspikes(cll,l,d) = length(allg1phases);
                    
                    clear r1phasmat; clear allr1phases; clear phasmat; clear allphases;
                    if ~isempty(thisruninds)
                        r1fr(cll,l,d) = nanmean(frs(thisruninds));
                        [r1_lfpspect(cll,l,d,:),chf,r1_lfpspecterr(cll,l,d,:,:)] = mtspectrumc(squeeze(lfpresp(thisruninds,1001:1800))',params);
                        %                     [r1_sfcoher(cll,l,d,:),a,b,c,de,sfcfx,e,f,g,r1_sfcoherr_c(cll,l,d,:,:)] = coherencycpt(lfpresp(thisruninds,1001:1800)',ptresp(thisruninds),params);
                        r1_tbtlfpspect(cll,l,d,:) = nanmean(lfpspecttbt(thisruninds,:),1);
                        r1_tbtlfpspecterr(cll,l,d,:) = nanstd(lfpspecttbt(thisruninds,:),1,1)./sqrt(length(thisruninds));
                        lg = [lg; ones(length(thisruninds),1).*l];
                        sg = [sg; ones(length(thisruninds),1).*d];
                        r1_ntrials(cll,l,d) = length(thisruninds);
                        r1phasmat = gamma1phases(thisruninds,respwin);
                        allr1phases = squeeze(r1phasmat(find(r1phasmat)));
                        if isrow(allr1phases), allr1phases = allr1phases'; end
                        condr1phases{cll,l,d} = allr1phases;
                        prefr1phase(cll,l,d) = circ_mean(allr1phases);
                        r1lockpval(cll,l,d) = circ_rtest(allr1phases);
                        r1plv(cll,l,d) = circ_r(allr1phases);
                        r1ppc(cll,l,d) = ppc(allr1phases);
                        [r1angdev(cll,l,d),r1cstd(cll,l,d)] = circ_std(allr1phases);
                        for b = 1:size(bandphases,1)
                            phasmat = squeeze(bandphases(b,thisruninds,respwin));
                            allphases = squeeze(phasmat(find(phasmat)));
                            if isrow(allphases), allphases = allphases'; end
                            condr1bandphases{cll,b,l,d} = allphases;
                            prefr1bandphase(cll,b,l,d) = circ_mean(allphases);
                            r1bandlockpval(cll,b,l,d) = circ_rtest(allphases);
                            r1bandplv(cll,b,l,d) = circ_r(allphases);
                            r1bandppc(cll,b,l,d) = ppc(allphases);
                            [r1bandangdev(cll,b,l,d),r1bandcstd(cll,b,l,d)] = circ_std(allphases);
                            
                            powmat = squeeze(bandpowers(b,thisruninds,respwin));
                            allpowers = squeeze(powmat(find(powmat)));
                            if isrow(allpowers), allpowers = allpowers'; end
                            condr1bandpowers{cll,b,l,d} = allpowers;
                        end
                        r1nspikes(cll,l,d) = length(allr1phases);
                    else
                        r1fr(cll,l,d) = NaN;
                        r1_lfpspect(cll,l,d,:) = nan(1,length(lfpspect(cll,l,d,:)));
                        r1_tbtlfpspect(cll,l,d,:) = nan(1,length(lfpspect(cll,l,d,:)));
                        r1_tbtlfpspecterr(cll,l,d,:) = nan(1,length(lfpspect(cll,l,d,:)));
                        r1_lfpspecterr(cll,l,d,:,:) = nan(2,length(lfpspect(cll,l,d,:)));
                        %                     r1_sfcoher(cll,l,d,:) = nan(1,length(sfcoher(cll,l,d,:)));
                        %                     r1_sfcohererr(cll,l,d,:,:) = nan(2,length(sfcoher(cll,l,d,:)));
                        r1_ntrials(cll,l,d) = 0;
                        condr1phases{cll,l,d} = NaN;
                        r1_phases{cll,l,d} = NaN;
                        prefr1phase(cll,l,d) = NaN;
                        r1lockpval(cll,l,d) = NaN;
                        r1plv(cll,l,d) = NaN;
                        r1ppc(cll,l,d) = NaN;
                        r1angdev(cll,l,d) = NaN;
                        r1cstd(cll,l,d) = NaN;
                        for b = 1:size(bandphases,1)
                            prefr1bandphase(cll,b,l,d) = NaN;
                            r1bandlockpval(cll,b,l,d) = NaN;
                            r1bandplv(cll,b,l,d) = NaN;
                            r1bandppc(cll,b,l,d) = NaN;
                            r1bandangdev(cll,b,l,d) = NaN;
                            r1bandcstd(cll,b,l,d) = NaN;
                            condr1bandphases{cll,b,l,d} = NaN;
                            condr1bandpowers{cll,b,l,d} = NaN;
                        end
                        r1nspikes(cll,l,d) = NaN;
                    end
                    
                    clear r0phasmat; clear allr0phases; clear phasmat; clear allphases;
                    if ~isempty(thisstillinds)
                        r0fr(cll,l,d) = nanmean(frs(thisstillinds));
                        [r0_lfpspect(cll,l,d,:),chf,r0_lfpspecterr(cll,l,d,:,:)] = mtspectrumc(squeeze(lfpresp(thisstillinds,1001:1800))',params);
                        %                     [r0_sfcoher(cll,l,d,:),a,b,c,de,sfcfx,e,f,g,r0_sfcoherr_c(cll,l,d,:,:)] = coherencycpt(lfpresp(thisstillinds,1001:1800)',ptresp(thisstillinds),params);
                        r0_tbtlfpspect(cll,l,d,:) = nanmean(lfpspecttbt(thisstillinds,:),1);
                        r0_tbtlfpspecterr(cll,l,d,:) = nanstd(lfpspecttbt(thisstillinds,:),1,1)./sqrt(length(thisstillinds));
                        r0_ntrials(cll,l,d) = length(thisstillinds);
                        r0phasmat = gamma1phases(thisstillinds,respwin);
                        allr0phases = squeeze(r0phasmat(find(r0phasmat)));
                        if isrow(allr0phases), allr0phases = allr0phases'; end
                        condr0phases{cll,l,d} = allr0phases;
                        prefr0phase(cll,l,d) = circ_mean(allr0phases);
                        r0lockpval(cll,l,d) = circ_rtest(allr0phases);
                        r0plv(cll,l,d) = circ_r(allr0phases);
                        r0ppc(cll,l,d) = ppc(allr0phases);
                        [r0angdev(cll,l,d),r0cstd(cll,l,d)] = circ_std(allr0phases);
                        for b = 1:size(bandphases,1)
                            phasmat = squeeze(bandphases(b,thisstillinds,respwin));
                            allphases = squeeze(phasmat(find(phasmat)));
                            if isrow(allphases), allphases = allphases'; end
                            condr0bandphases{cll,b,l,d} = allphases;
                            prefr0bandphase(cll,b,l,d) = circ_mean(allphases);
                            r0bandlockpval(cll,b,l,d) = circ_rtest(allphases);
                            r0bandplv(cll,b,l,d) = circ_r(allphases);
                            r0bandppc(cll,b,l,d) = ppc(allphases);
                            [r0bandangdev(cll,b,l,d),r0bandcstd(cll,b,l,d)] = circ_std(allphases);
                            
                            powmat = squeeze(bandpowers(b,thisstillinds,respwin));
                            allpowers = squeeze(powmat(find(powmat)));
                            if isrow(allpowers), allpowers = allpowers'; end
                            condr0bandpowers{cll,b,l,d} = allpowers;
                        end
                        r0nspikes(cll,l,d) = length(allr0phases);
                    else
                        r0fr(cll,l,d) = NaN;
                        r0_lfpspect(cll,l,d,:) = nan(1,length(lfpspect(cll,l,d,:)));
                        r0_lfpspecterr(cll,l,d,:,:) = nan(2,length(lfpspect(cll,l,d,:)));
                        r0_tbtlfpspect(cll,l,d,:) = nan(1,length(lfpspect(cll,l,d,:)));
                        r0_tbtlfpspecterr(cll,l,d,:) = nan(1,length(lfpspect(cll,l,d,:)));
                        %                     r0_sfcoher(cll,l,d,:) = nan(1,length(sfcoher(cll,l,d,:)));
                        %                     r0_sfcohererr(cll,l,d,:,:) = nan(2,length(sfcoher(cll,l,d,:)));
                        r0_ntrials(cll,l,d) = 0;
                        condr0phases{cll,l,d} = NaN;
                        r0_phases{cll,l,d} = NaN;
                        prefr0phase(cll,l,d) = NaN;
                        r0lockpval(cll,l,d) = NaN;
                        r0plv(cll,l,d) = NaN;
                        r0ppc(cll,l,d) = NaN;
                        r0angdev(cll,l,d) = NaN;
                        r0cstd(cll,l,d) = NaN;
                        for b = 1:size(bandphases,1)
                            prefr0bandphase(cll,b,l,d) = NaN;
                            r0bandlockpval(cll,b,l,d) = NaN;
                            r0bandplv(cll,b,l,d) = NaN;
                            r0bandppc(cll,b,l,d) = NaN;
                            r0bandangdev(cll,b,l,d) = NaN;
                            r0bandcstd(cll,b,l,d) = NaN;
                            condr0bandphases{cll,b,l,d} = NaN;
                            condr0bandpowers{cll,b,l,d} = NaN;
                        end
                        r0nspikes(cll,l,d) = NaN;
                    end
                end
            end
            depth(cll) = result.depth;
            cll = cll+1;
        end
    end
    save(popfile, '-v7.3');
else
    load(popfile);
end

for i = 1:size(r1fr,1)
    lgpow(i,:,:) = squeeze(lfpspect(i,:,:,lgi(i))); % at second contact only
    r1lgpow(i,:,:) = squeeze(r1_tbtlfpspect(i,:,:,lgi(i)));
    r1lgpowerr(i,:,:) = squeeze(r1_tbtlfpspecterr(i,:,:,lgi(i)));
    normlgpow(i,:) = r1lgpow(i,:,2)./r1lgpow(i,1,2);
    normr1fr(i,:) = r1fr(i,:,2)./r1fr(i,1,2);
    for l = 1:6
        for d = 1:2
            r1tbtfillspecy(i,l,d,:) = [(squeeze(r1_tbtlfpspect(i,l,d,1:104))+squeeze(r1_tbtlfpspecterr(i,l,d,1:104)))',fliplr((squeeze(r1_tbtlfpspect(i,l,d,1:104))-squeeze(r1_tbtlfpspecterr(i,l,d,1:104)))')];
        end
    end    
end
fillx = [chf(1:104),fliplr(chf(1:104))];

% spectrum small vs large
for i =1:size(r1tbtfillspecy,1)
    figure
    fill(fillx,squeeze(r1tbtfillspecy(i,1,2,:)),'b')
    hold on
    fill(fillx,squeeze(r1tbtfillspecy(i,2,2,:)),'c')
    fill(fillx,squeeze(r1tbtfillspecy(i,3,2,:)),'g')
    fill(fillx,squeeze(r1tbtfillspecy(i,4,2,:)),'k')
    fill(fillx,squeeze(r1tbtfillspecy(i,5,2,:)),'m')
    fill(fillx,squeeze(r1tbtfillspecy(i,6,2,:)),'r')
    set(gca,'yscale','log')
end

figure
hold on
for i = 1:size(r1lgpow,1)
    errorbar(lightlevels,r1lgpow(i,:,2),r1lgpowerr(i,:,2))
end
title('non-normalized gamma power')
ylabel('gamma power')
xlabel('light intensity')



figure
hold on
for i = 1:size(watts,1)
    if color(i) == 3, pc = 'r'; else pc = 'g'; end
    plot(watts(i,:),normlgpow(i,:),'g')
end
title('normalized gamma power')
ylabel('gamma power')
xlabel('light intensity')



