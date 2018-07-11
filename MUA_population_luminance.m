function MUA_population_luminance


% temp population 
animalids = {'160319', '160320', '160321', '160323', '160324', '160328'};
blocks    = [ 9,        7,        6,        5,        5,        6];
animal    = [ 1,        2,        3,        4,        5,        6];
electrodes =[[17,32];  [1,16];   [1,16];   [1,16];   [1,16];   [1,16]];
bldepth  =  [ 500,      550,      550,      500,      500,      500];
penangle =  [ 25,       25,       25,       25,       25,       25];
spacing =   [ 25,       25,       25,       25,       25,       25];
popfile = 'C:\Users\Julia\work\data\populations\templum.mat';

% % temp population 
% animalids = {'160323'};
% blocks    = [ 5];
% animal    = [ 1];
% electrodes =[[1,16]];
% bldepth  =  [ 500];
% penangle =  [ 25];
% spacing =   [ 25];
% popfile = 'C:\Users\Julia\work\data\populations\templum.mat';


depth = bldepth.*(cosd(penangle)*cosd(22));
spacing = spacing.*(cosd(penangle)*cosd(22));
evaldepth = 300;
for i = 1:length(depth)
    for j = 1:length(electrodes(i,1):electrodes(i,2))
        dm(i,j) = depth(i)-((j-1)*spacing(i)); % depth matrix
    end
    [c,di(i)] = min(abs(dm(i,:)-evaldepth)); % get depth index, least distance to evaldepth
end

recalculate_pop = 0;
recalculate_muafile = 0;

% chronux parameters
params.tapers = [2,5]; params.Fs = 1000; params.err = [2, 0.05]; params.trialave = 1;

if ~exist(popfile) || recalculate_pop

    cll = 1;
    for blck = 1:length(blocks)
        
        basepath = strcat('C:\Users\Julia\work\data\', animalids{blck}, '\');
        file = strcat(basepath, 'muaresult_', int2str(blocks(blck)), '_', int2str(electrodes(blck,1)), '-', int2str(electrodes(blck,2)), '.mat');
        oldfile = strcat(basepath, 'muaresult_', int2str(blocks(blck)), '_', int2str(electrodes(blck,1)), ':', int2str(electrodes(blck,2)), '.mat');
        
        clear result;
        if ~exist(file) || recalculate_muafile
            if ~exist(oldfile)
                result = MUAdataprepare(basepath,animalids{blck},blocks(blck),electrodes(blck,1):electrodes(blck,2));
                save(file,'result')
            else
                load(oldfile)
                save(file,'result');
            end
        else
            clear centresult; clear surrresult;
            load(file);
            if exist('centresult')
                result = centresult;
            elseif exist('surrresult')
                result = surrresult;
            end            
        end      
       
        
        prestim = 300;
        poststim = 700;
        respwin = 501:1500; % after stimulus onset
        offsetwin = 1501:2500;
        respwin = respwin+prestim;        
        
        
        ch = di(blck);
        
        disp(['Block ' int2str(blck) '/' int2str(length(blocks)) '   channel ' int2str(ch) '/' int2str(length(electrodes(blck,1):electrodes(blck,2)))]);
              
        trialdur = result.stimduration*1000;
        msstamps = result.msstamps;
        if length(msstamps)~=length(result.light)
%             msstamps([193,291]) = []; % for 151023 block 5
%             msstamps([52, 98]) = []; % for 160210 block 13
%             msstamps([87]) = []; % for 160321 block 6
%             msstamps([37]) = []; % for 160319 block 9
%             msstamps([40,152]) = []; % for 160320 block 7
%             msstamps([37]) = []; % for 160319 block 9
%             msstamps([297]) = []; % for 160319 block 9
%             result.msstamps = msstamps;
%             save(file,'result');
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
        large = find(result.gratingInfo.size == max(sizes) & result.light == 0 & result.gratingInfo.Contrast ~= 0);
        small = find(result.gratingInfo.Contrast == 0 & result.light == 0);
        lr1 = intersect(large,oktrials);
        sr1 = intersect(small,oktrials);
        for i = 1:min(length(large),length(small))
            [pl(i,:),f] = mtspectrumc(result.lfp(ch,result.msstamps(large(i))+700:result.msstamps(large(i))+1500),params);
            [ps(i,:),f] = mtspectrumc(result.lfp(ch,result.msstamps(small(i))+700:result.msstamps(small(i))+1500),params);
        end
        plr1 = nan(1,size(pl,2)); psr1 = nan(1,size(pl,2));
        if length(lr1>=5)
            for i = 1:length(lr1)
                [plr1(i,:),f] = mtspectrumc(result.lfp(ch,result.msstamps(lr1(i))+700:result.msstamps(lr1(i))+1500),params);
                %                     [plr1(i,:),f] = pmtm(centresult.lfp(ch,centresult.msstamps(lr1(i))+300:centresult.msstamps(lr1(i))+1300),3,[],1000);
            end
        end
        if length(sr1>=5)
            for i = 1:length(sr1)
                [psr1(i,:),f] = mtspectrumc(result.lfp(ch,result.msstamps(sr1(i))+700:result.msstamps(sr1(i))+1500),params);
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
        lgi(blck) = bpi; hgi(blck) = gpi;
        
        lfp = result.lfp(ch,:); sr = 1000;
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
        
        msStimes = round(result.msStimes{ch});
        if ~isempty(msStimes) & msStimes(1) == 0, msStimes(1) = 1; end
        
        chan = zeros(1,size(result.lfp,2));
        chan(msStimes) = 1;
        
        clear gamma1phases; clear gamma2phases; clear gamma1powers; clear gamma2powers;
        clear bandphases; clear bandpowers;
        for i = 1:length(msstamps)
            resp(i,:) = chan(msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim);
            hh = find(resp(i,1001:1795))'; ptresp(i).times = hh./1000;
            lfpresp(i,:) = result.lfp(ch,msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim);
            lfpspecttbt(i,:) = mtspectrumc(squeeze(lfpresp(i,1001:1800))',params);
            muacresp(i,:) =  result.muac(ch,msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim).^2;
                        
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
        muacmean = mean(muacresp(:,respwin),2);
        
        r1anovavec = []; 
        lg = []; sg = []; % light group and stim group
        lightlevels = unique(result.light);
        luminances = unique(result.gratingInfo.luminance);
        contrasts = unique(result.gratingInfo.Contrast);
%         for l = 1:length(lightlevels)
%             for d = 1:length(sizes)
%                 for g = 1:length(luminances)
%                     for c = 1:length(contrasts)
        for l = 1
            for d = 1:length(sizes)
                for g = 1:length(luminances)
                    for c = 1
                
                        thisinds = find(result.gratingInfo.size == sizes(d) & result.light == lightlevels(l) & result.gratingInfo.luminance == luminances(g) & result.gratingInfo.Contrast == contrasts(c));
                        thisruninds = intersect(thisinds,oktrials); thisstillinds = intersect(thisinds,stilltrials);

                        condresp(blck,l,d,g,c,:) = nanmean(resp(thisinds,:),1);
                        condfr(blck,l,d,g,c) = nanmean(frs(thisinds));
                        condmuacresp(blck,l,d,g,c,:) = nanmean(muacresp(thisinds,:));
                        condmua(blck,l,d,g,c) = nanmean(muacmean(thisinds));
                        condg1pow(blck,l,d,g,c) = nanmean(g1pow(thisinds));
                        condg2pow(blck,l,d,g,c) = nanmean(g2pow(thisinds));

                        condlfpresp(blck,l,d,g,c,:) = nanmean(lfpresp(thisinds,:));                
                        [lfpspect(blck,l,d,g,c,:),chf,lfpspecterr(blck,l,d,g,c,:,:)] = mtspectrumc(squeeze(lfpresp(thisinds,1001:1800))',params);
                        [muacspect(blck,l,d,g,c,:),chf,muacspecterr_c(blck,l,d,g,c,:,:)] = mtspectrumc(squeeze(muacresp(thisinds,1001:1800))',params);
                        tbtlfpspect(blck,l,d,g,c,:) = nanmean(lfpspecttbt(thisinds,:),1);
                        tbtlfpspecterr(blck,l,d,g,c,:) = nanstd(lfpspecttbt(thisinds,:),1,1)./sqrt(length(thisinds));
        %                 [sfcoher(blck,l,d,:),a,b,c,de,sfcfx,e,f,g,sfcoherr_c(blck,l,d,:,:)] = coherencycpt(lfpresp(thisinds,1001:1800)',ptresp(thisinds),params);

                        condg1powresp(blck,l,d,g,c,:) = mean(gamma1powresp(thisinds,:));
                        condg2powresp(blck,l,d,g,c,:) = mean(gamma2powresp(thisinds,:));

                        clear g1phasmatl; clear allg1phases;
                        g1phasmat = gamma1phases(thisinds,respwin);
                        allg1phases = squeeze(g1phasmat(find(g1phasmat)));
                        if isrow(allg1phases), allg1phases = allg1phases'; end
                        condg1phases{blck,l,d,g,c} = allg1phases;
                        prefg1phase(blck,l,d,g,c) = circ_mean(allg1phases);
                        g1lockpval(blck,l,d,g,c) = circ_rtest(allg1phases);
                        g1plv(blck,l,d,g,c) = circ_r(allg1phases);
                        g1ppc(blck,l,d,g,c) = ppc(allg1phases);
                        [g1angdev(blck,l,d,g,c),g1cstd(blck,l,d,g,c)] = circ_std(allg1phases);

                        clear g2phasmatl; clear allg2phases;
                        g2phasmat = gamma2phases(thisinds,respwin);
                        allg2phases = squeeze(g2phasmat(find(g2phasmat)));
                        if isrow(allg2phases), allg2phases = allg2phases'; end
                        condg2phases{blck,l,d,g,c} = allg2phases;
                        prefg2phase(blck,l,d,g,c) = circ_mean(allg2phases);
                        g2lockpval(blck,l,d,g,c) = circ_rtest(allg2phases);
                        g2plv(blck,l,d,g,c) = circ_r(allg2phases);
                        g2ppc(blck,l,d,g,c) = ppc(allg2phases);
                        [g2angdev(blck,l,d,g,c),g2cstd(blck,l,d,g,c)] = circ_std(allg2phases);

                        clear phasmat; clear allphases;
                        for b = 1:size(bandphases,1)
                            phasmat = squeeze(bandphases(b,thisinds,respwin));
                            allphases = squeeze(phasmat(find(phasmat)));
                            if isrow(allphases), allphases = allphases'; end
                            condbandphases{blck,b,l,d,g} = allphases;
                            bandprefphase(blck,b,l,d,g) = circ_mean(allphases);
                            bandlockpval(blck,b,l,d,g) = circ_rtest(allphases);
                            bandplv(blck,b,l,d,g) = circ_r(allphases);
                            bandppc(blck,b,l,d,g) = ppc(allphases);
                            [bandangdev(blck,b,l,d,g),bandcstd(blck,b,l,d,g)] = circ_std(allphases);                    

                            powmat = squeeze(bandpowers(b,thisinds,respwin));
                            allpowers = squeeze(powmat(find(powmat)));
                            if isrow(allpowers), allpowers = allpowers'; end
                            condbandpowers{blck,b,l,d,g} = allpowers;
                        end

                        gnspikes(blck,l,d) = length(allg1phases);

                        clear r1phasmat; clear allr1phases; clear phasmat; clear allphases;
                        if ~isempty(thisruninds)
                            r1condfr(blck,l,d,g,c) = nanmean(frs(thisruninds));
                            [r1_lfpspect(blck,l,d,g,c,:),chf,r1_lfpspecterr(blck,l,d,g,c,:,:)] = mtspectrumc(squeeze(lfpresp(thisruninds,1001:1800))',params);
        %                     [r1_sfcoher(blck,l,d,:),a,b,c,de,sfcfx,e,f,g,r1_sfcoherr_c(blck,l,d,:,:)] = coherencycpt(lfpresp(thisruninds,1001:1800)',ptresp(thisruninds),params);
                            r1_tbtlfpspect(blck,l,d,g,c,:) = nanmean(lfpspecttbt(thisruninds,:),1);                    
                            r1_tbtlfpspecterr(blck,l,d,g,c,:) = nanstd(lfpspecttbt(thisruninds,:),1,1)./sqrt(length(thisruninds));
                            r1anovavec = [r1anovavec; lfpspecttbt(thisruninds,lgi(blck))];
                            lg = [lg; ones(length(thisruninds),1).*l];
                            sg = [sg; ones(length(thisruninds),1).*d];
                            r1_ntrials(blck,l,d,g,c) = length(thisruninds);
                            r1phasmat = gamma1phases(thisruninds,respwin);
                            allr1phases = squeeze(r1phasmat(find(r1phasmat)));
                            if isrow(allr1phases), allr1phases = allr1phases'; end
                            condr1phases{blck,l,d,g,c} = allr1phases;
                            prefr1phase(blck,l,d,g,c) = circ_mean(allr1phases);
                            r1lockpval(blck,l,d,g,c) = circ_rtest(allr1phases);
                            r1plv(blck,l,d,g,c) = circ_r(allr1phases);
                            r1ppc(blck,l,d,g,c) = ppc(allr1phases);
                            [r1angdev(blck,l,d,g,c),r1cstd(blck,l,d,g,c)] = circ_std(allr1phases);
                            for b = 1:size(bandphases,1)
                                phasmat = squeeze(bandphases(b,thisruninds,respwin));
                                allphases = squeeze(phasmat(find(phasmat)));
                                if isrow(allphases), allphases = allphases'; end
                                condr1bandphases{blck,b,l,d,g} = allphases;
                                prefr1bandphase(blck,b,l,d,g) = circ_mean(allphases);
                                r1bandlockpval(blck,b,l,d,g) = circ_rtest(allphases);
                                r1bandplv(blck,b,l,d,g) = circ_r(allphases);
                                r1bandppc(blck,b,l,d,g) = ppc(allphases);
                                [r1bandangdev(blck,b,l,d,g),r1bandcstd(blck,b,l,d,g)] = circ_std(allphases);

                                powmat = squeeze(bandpowers(b,thisruninds,respwin));
                                allpowers = squeeze(powmat(find(powmat)));
                                if isrow(allpowers), allpowers = allpowers'; end
                                condr1bandpowers{blck,b,l,d,g} = allpowers;
                            end
                            r1nspikes(blck,l,d,g,c) = length(allr1phases);
                        else
                            r1condfr(blck,l,d,g,c) = NaN;
                            r1_lfpspect(blck,l,d,g,c,:) = nan(1,length(lfpspect(blck,l,d,g,c,:)));
                            r1_tbtlfpspect(blck,l,d,g,c,:) = nan(1,length(lfpspect(blck,l,d,g,c,:)));                    
                            r1_tbtlfpspecterr(blck,l,d,g,c,:) = nan(1,length(lfpspect(blck,l,d,g,c,:))); 
                            r1_lfpspecterr(blck,l,d,g,c,:,:) = nan(2,length(lfpspect(blck,l,d,g,c,:)));
        %                     r1_sfcoher(blck,l,d,:) = nan(1,length(sfcoher(blck,l,d,:)));
        %                     r1_sfcohererr(blck,l,d,:,:) = nan(2,length(sfcoher(blck,l,d,:)));
                            r1_ntrials(blck,l,d,g,c) = 0;
                            condr1phases{blck,l,d,g,c} = NaN;
                            r1_phases{blck,l,d,g,c} = NaN;
                            prefr1phase(blck,l,d,g,c) = NaN;
                            r1lockpval(blck,l,d,g,c) = NaN;
                            r1plv(blck,l,d,g,c) = NaN;
                            r1ppc(blck,l,d,g,c) = NaN;
                            r1angdev(blck,l,d,g,c) = NaN;
                            r1cstd(blck,l,d,g,c) = NaN;
                            for b = 1:size(bandphases,1)
                                prefr1bandphase(blck,b,l,d,g) = NaN;
                                r1bandlockpval(blck,b,l,d,g) = NaN;
                                r1bandplv(blck,b,l,d,g) = NaN;
                                r1bandppc(blck,b,l,d,g) = NaN;
                                r1bandangdev(blck,b,l,d,g) = NaN;
                                r1bandcstd(blck,b,l,d,g) = NaN;
                                condr1bandphases{blck,b,l,d,g} = NaN;
                                condr1bandpowers{blck,b,l,d,g} = NaN;
                            end
                            r1nspikes(blck,l,d,g,c) = NaN;
                        end

                        clear r0phasmat; clear allr0phases; clear phasmat; clear allphases;
                        if ~isempty(thisstillinds)
                            r0condfr(blck,l,d,g,c) = nanmean(frs(thisstillinds));
                            [r0_lfpspect(blck,l,d,g,c,:),chf,r0_lfpspecterr(blck,l,d,g,c,:,:)] = mtspectrumc(squeeze(lfpresp(thisstillinds,1001:1800))',params);
        %                     [r0_sfcoher(blck,l,d,:),a,b,c,de,sfcfx,e,f,g,r0_sfcoherr_c(blck,l,d,:,:)] = coherencycpt(lfpresp(thisstillinds,1001:1800)',ptresp(thisstillinds),params);
                            r0_tbtlfpspect(blck,l,d,g,c,:) = nanmean(lfpspecttbt(thisstillinds,:),1);                    
                            r0_tbtlfpspecterr(blck,l,d,g,c,:) = nanstd(lfpspecttbt(thisstillinds,:),1,1)./sqrt(length(thisstillinds));
                            r0_ntrials(blck,l,d,g,c) = length(thisstillinds);
                            r0phasmat = gamma1phases(thisstillinds,respwin);
                            allr0phases = squeeze(r0phasmat(find(r0phasmat)));
                            if isrow(allr0phases), allr0phases = allr0phases'; end
                            condr0phases{blck,l,d,g,c} = allr0phases;
                            prefr0phase(blck,l,d,g,c) = circ_mean(allr0phases);
                            r0lockpval(blck,l,d,g,c) = circ_rtest(allr0phases);
                            r0plv(blck,l,d,g,c) = circ_r(allr0phases);
                            r0ppc(blck,l,d,g,c) = ppc(allr0phases);
                            [r0angdev(blck,l,d,g,c),r0cstd(blck,l,d,g,c)] = circ_std(allr0phases);
                            for b = 1:size(bandphases,1)
                                phasmat = squeeze(bandphases(b,thisstillinds,respwin));
                                allphases = squeeze(phasmat(find(phasmat)));
                                if isrow(allphases), allphases = allphases'; end
                                condr0bandphases{blck,b,l,d,g} = allphases;
                                prefr0bandphase(blck,b,l,d,g) = circ_mean(allphases);
                                r0bandlockpval(blck,b,l,d,g) = circ_rtest(allphases);
                                r0bandplv(blck,b,l,d,g) = circ_r(allphases);
                                r0bandppc(blck,b,l,d,g) = ppc(allphases);
                                [r0bandangdev(blck,b,l,d,g),r0bandcstd(blck,b,l,d,g)] = circ_std(allphases);

                                powmat = squeeze(bandpowers(b,thisstillinds,respwin));
                                allpowers = squeeze(powmat(find(powmat)));
                                if isrow(allpowers), allpowers = allpowers'; end
                                condr0bandpowers{blck,b,l,d,g} = allpowers;
                            end
                            r0nspikes(blck,l,d) = length(allr0phases);
                        else
                            r0condfr(blck,l,d,g,c) = NaN;
                            r0_lfpspect(blck,l,d,g,c,:) = nan(1,length(lfpspect(blck,l,d,g,c,:)));
                            r0_lfpspecterr(blck,l,d,g,c,:,:) = nan(2,length(lfpspect(blck,l,d,g,c,:)));
                            r0_tbtlfpspect(blck,l,d,g,c,:) = nan(1,length(lfpspect(blck,l,d,g,c,:)));                    
                            r0_tbtlfpspecterr(blck,l,d,g,c,:) = nan(1,length(lfpspect(blck,l,d,g,c,:))); 
        %                     r0_sfcoher(blck,l,d,:) = nan(1,length(sfcoher(blck,l,d,:)));
        %                     r0_sfcohererr(blck,l,d,:,:) = nan(2,length(sfcoher(blck,l,d,:)));
                            r0_ntrials(blck,l,d,g,c) = 0;
                            condr0phases{blck,l,d,g,c} = NaN;
                            r0_phases{blck,l,d,g,c} = NaN;
                            prefr0phase(blck,l,d,g,c) = NaN;
                            r0lockpval(blck,l,d,g,c) = NaN;
                            r0plv(blck,l,d,g,c) = NaN;
                            r0ppc(blck,l,d,g,c) = NaN;
                            r0angdev(blck,l,d,g,c) = NaN;
                            r0cstd(blck,l,d,g,c) = NaN;
                            for b = 1:size(bandphases,1)
                                prefr0bandphase(blck,b,l,d,g) = NaN;
                                r0bandlockpval(blck,b,l,d,g) = NaN;
                                r0bandplv(blck,b,l,d,g) = NaN;
                                r0bandppc(blck,b,l,d,g) = NaN;
                                r0bandangdev(blck,b,l,d,g) = NaN;
                                r0bandcstd(blck,b,l,d,g) = NaN;
                                condr0bandphases{blck,b,l,d,g} = NaN;
                                condr0bandpowers{blck,b,l,d,g} = NaN;
                            end
                            r0nspikes(blck,l,d,g,c) = NaN;
                        end
                    end
                end
            end
        end
        
        
        
%         figure
%         title('change with contrast and size at highest luminance')
%         semilogy(f,squeeze(mean(tbtlfpspect(1,1,:,5,1,:),3)),'k','linewidth',2)
%         hold on
%         semilogy(f,squeeze(tbtlfpspect(1,1,1,5,2,:)),'c','linewidth',2)
%         semilogy(f,squeeze(tbtlfpspect(1,1,1,5,3,:)),'b','linewidth',2)
%         semilogy(f,squeeze(tbtlfpspect(1,1,2,5,2,:)),'g','linewidth',2)
%         semilogy(f,squeeze(tbtlfpspect(1,1,2,5,3,:)),'r','linewidth',2)
%         legend('no contrast','small low','small high','large low','large high');
%         axis([0,100,0,500])
%         xlabel('frequency')
%         ylabel('power')
        
%         if ~isempty(r1anovavec)
%             [anovap(blck,:),table,stats] = anovan(r1anovavec,{lg,sg}); % 1 = light, 2 = stimulus
%         end
    end
    save(popfile, '-v7.3');
else
    load(popfile);
end

validinds = [1,3,4,5];

for i = validinds
    figure
    title('plain luminance')
    semilogy(f,smooth(squeeze(mean(tbtlfpspect(i,1,:,1,1,:),3))),'k','linewidth',3)
    hold on
    semilogy(f,smooth(squeeze(mean(tbtlfpspect(i,1,:,2,1,:),3))),'color',[.3,.3,.3],'linewidth',3)
    semilogy(f,smooth(squeeze(mean(tbtlfpspect(i,1,:,3,1,:),3))),'color',[.5,.5,.5],'linewidth',3)
    semilogy(f,smooth(squeeze(mean(tbtlfpspect(i,1,:,4,1,:),3))),'color',[.7,.7,.7],'linewidth',3)
    semilogy(f,smooth(squeeze(mean(tbtlfpspect(i,1,:,5,1,:),3))),'color',[.8,.8,.8],'linewidth',3)
    legend('0','0.25','0.5','0.75','1')
    axis([0,100,2,800])
    xlabel('frequency')
    ylabel('power')
end

for i = 1:6
    hgpow(i,:) = squeeze(nanmean(tbtlfpspect(i,1,:,:,(hgi(i))),3));
end
figure
[p,s] = signrank(hgpow(validinds,1),hgpow(validinds,5));
plot(1,hgpow(validinds,1),'ko','markerfacecolor','k');
hold on
plot(2,hgpow(validinds,5),'ko','markerfacecolor','k');
for i = validinds
    plot([1,2],[hgpow(i,1),hgpow(i,5)],'k')
end
set(gca,'xticklabel',{'black','white'})
set(gca,'xtick',[1,2]);
axis([0,3,0,50])
ylabel('high gamma power')
set(gcf,'OuterPosition',[573   504   242   513])
title(['n: 4 p: ' num2str(p)]);

    

for i = 1:length(validinds)
    
    b1 = find(f>beta(1),1); b2 = find(f>beta(2),1);
    g1 = find(f>gamma(1),1); g2 = find(f>gamma(2),1);
    for s = 1:5 % luminances here
        smoothspec(i,s,:) = smooth(squeeze(mean(tbtlfpspect(validinds(i),1,:,s,1,:),3)));
        
        bsig = squeeze(smoothspec(i,s,b1:b2));        
        gsig = squeeze(smoothspec(i,s,g1:g2));
        
        if isempty(find(diff(bsig)>0)) % there is no clear beta peak
            slgi(i,s) = NaN;
        else
            peaks = find(diff(bsig)>0)+1;
            pvs = bsig(peaks);
            bpi = peaks(pvs == max(pvs));
            bpi = bpi+b1-1;
            slgi(i,s) = bpi;
        end
        
        if isempty(find(diff(gsig)>0)) % there is no clear beta peak
            shgi(i,s) = NaN;
        else
            peaks = find(diff(gsig)>0)+1;
            pvs = gsig(peaks);
            gpi = peaks(pvs == max(pvs));
            gpi = gpi+g1-1;
            shgi(i,s) = gpi;
        end
    end
    
end

for i = 1:size(shgi,1)
    figure
    title('plain luminance')
    semilogy(f,squeeze(smoothspec(i,1,:)),'k','linewidth',2)
    hold on
    semilogy(f,squeeze(smoothspec(i,2,:)),'color',[.3,.3,.3],'linewidth',2)
    semilogy(f,squeeze(smoothspec(i,3,:)),'color',[.5,.5,.5],'linewidth',2)
    semilogy(f,squeeze(smoothspec(i,4,:)),'color',[.7,.7,.7],'linewidth',2)
    semilogy(f,squeeze(smoothspec(i,5,:)),'color',[.8,.8,.8],'linewidth',2)
    if ~isnan(shgi(i,1))
        plot(f(shgi(i,1)),smoothspec(i,1,shgi(i,1)),'ko','markerfacecolor','k')
    end
    if ~isnan(shgi(i,2))
        plot(f(shgi(i,2)),smoothspec(i,2,shgi(i,2)),'o','color',[.3,.3,.3],'markerfacecolor',[.3,.3,.3])
    end
    if ~isnan(shgi(i,3))
        plot(f(shgi(i,3)),smoothspec(i,3,shgi(i,3)),'o','color',[.5,.5,.5],'markerfacecolor',[.5,.5,.5])
    end
    plot(f(shgi(i,4)),smoothspec(i,4,shgi(i,4)),'o','color',[.7,.7,.7],'markerfacecolor',[.7,.7,.7])
    plot(f(shgi(i,5)),smoothspec(i,5,shgi(i,5)),'o','color',[.8,.8,.8],'markerfacecolor',[.8,.8,.8])
    legend('0','0.25','0.5','0.75','1')
    axis([0,100,0.5,500])
    xlabel('frequency')
    ylabel('power')
end

for i = 1:size(shgi,1)
    gcurve(i,:) = smoothspec(i,:,shgi(i,5));
    normcurve(i,:) = smoothspec(i,:,shgi(i,5))./smoothspec(i,5,shgi(i,5));
end

figure
plot(luminances,normcurve')

[p,t,s] = kruskalwallis(normcurve);
figure
errorbar(luminances,mean(normcurve),std(normcurve)./sqrt(size(normcurve,1)),'ko-','linewidth',2,'markerfacecolor','k')
xlabel('luminance')
ylabel('norm high gamma amplitude')
set(gca,'xtick',[0,.25,.5,.75,1])
axis([-.2,1.2,.2,1.1])
title(['high gamma with luminance n = 4, p: ' num2str(p)])

disp('');
