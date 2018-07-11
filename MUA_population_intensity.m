function MUA_population_intensity

% % SOM later population
% animalids = {'150331', '150401','150602','150603','150825','150902', '151110', '151209', '151221'};
% blocks    = [5,         9,       9,       8,       16,      15,       16,       14,       2];
% animal    = [1,         2,       3,       4,       5,       6,        7,        8,        9];
% electrodes =[[1,32];    [1,32];  [1,32];  [1,16]; [17,32]; [1,16];   [1,16];   [17,32];  [1,16]];
% bldepth  =  [1050,      1000,    1000,    550,     400,     400,      500,      500,      500];
% penangle =  [25,        25,      25,      10,      25,      25,       25,       25,       25];
% spacing =   [25,        25,      25,      25,      25,      25,       25,       25,       25];
% color =     [3,         3,       3,       3,       3,       3,        3,        3,        3];
% popfile = 'C:\Users\Julia\work\data\populations\SOM_Halo_later\intensity\SOM_intensity_MUA.mat';

% % PV Halo up to 50% population
% animalids = {'150804', '150804'};
% blocks    = [15,        15];
% animal    = [1,         1];
% electrodes =[[1,16];    [17,32]];
% bldepth  =  [400,       400];
% penangle =  [25,        25];
% spacing =   [25,        25];
% color =     [3,         3];
% popfile = 'C:\Users\Julia\work\data\populations\PV_Halo\intensity\PV50_intensity_MUA.mat';
% % 
% % PV Halo up to 100% population
% animalids = {'150820', '150820', '150823', '160320', '160324', '160328'};
% blocks    = [16,         16,      14,       6,        6,        7];
% animal    = [1,          1,       2,        3,        4,        5];
% electrodes =[[1,16];    [17,32]; [1,16];   [1,16];   [1,16];   [1,16]];
% bldepth  =  [400,       400,      400,      550,      500,      500];
% penangle =  [25,        25,       25,       25,       25,       25];
% spacing =   [25,        25,       25,       25,       25,       25];
% color =     [3,         3,        3,        3,        3,        3];
% popfile = 'C:\Users\Julia\work\data\populations\PV_Halo\intensity\PV100_intensity_MUA.mat';
% % 
% % PV eArch population - red up to 40%
% animalids = {'151211', '160114', '160115', '160205', '160210'};
% blocks    = [5,         7,        6,        14,       11];
% animal    = [1,         2,        3,        4,        5];
% electrodes =[[1,16];   [1,16];   [1,16];   [17,32];  [17,32]];
% bldepth  =  [500,       500,      500,      500,      500];
% penangle =  [25,        25,       25,       25,       25];
% spacing =   [25,        25,       25,       25,       25];
% color =     [3,         3,        3,        3,        3];
% popfile = 'C:\Users\Julia\work\data\populations\PV_eArch\intensity\PV_intensityred_MUA.mat';

% % PV eArch population - green low 
% animalids = {'151211', '160114', '160204', '160205', '160210'};
% blocks    = [6,         12,       11,       13,       10];
% animal    = [1,         2,        3,        4,        5];
% electrodes =[[1,16];   [1,16];   [1,16];   [17,32];  [17,32]];
% bldepth  =  [500,       500,      500,      500,      500];
% penangle =  [25,        25,       25,       25,       25];
% spacing =   [25,        25,       25,       25,       25];
% color =     [4,         4,        4,        4,        4];
% popfile = 'C:\Users\Julia\work\data\populations\PV_eArch\intensity\PV_intensitygreen_MUA.mat';

% % PV eArch population - green up to 40% 
% animalids = {'160114', '160115', '160204', '160205', '160210'};
% blocks    = [8,         8,        16,       15,       13];
% animal    = [1,         2,        3,        4,        5];
% electrodes =[[1,16];   [1,16];   [1,16];   [17,32];  [17,32]];
% bldepth  =  [500,       500,      500,      500,      500];
% penangle =  [25,        25,       25,       25,       25];
% spacing =   [25,        25,       25,       25,       25];
% color =     [4,         4,        4,        4,        4];
% popfile = 'C:\Users\Julia\work\data\populations\PV_eArch\intensity\PV_intensitygreen40_MUA.mat';

% % PV Arch reporter population - GERMLINE
% animalids = {'151216', '151216', '151216'};
% blocks    = [5,         6,        7];
% animal    = [1,         1,        1];
% electrodes =[[1,16];    [1,16];  [1,16]];
% bldepth  =  [500,       500,      500];
% penangle =  [25,        25,       25];
% spacing =   [25,        25,       25];
% color =     [4,         3,        4];
% popfile = 'C:\Users\Julia\work\data\populations\PV_Archtransg\intensity\PV_intensity_MUA.mat';
% 
% % control population - red
% animalids = {'151214', '151215', '151217', '151222', '151223', '160113'};  %160112 too noisy
% blocks    = [4,         8,        5,        9,        9,        7];
% animal    = [1,         2,        3,        4,        5,        6];
% electrodes =[[1,16];    [1,16];   [1,16];  [1,16];   [1,16];   [1,16]];
% bldepth  =  [500,       500,      500,      500,      500,      500];
% penangle =  [25,        25,       25,       25,       25,       25];
% spacing =   [25,        25,       25,       25,       25,       25];
% color =     [3,         3,        3,        3,        3,        3];
% popfile = 'C:\Users\Julia\work\data\populations\control\intensity\intensityred_MUA.mat';

% % control population - green
% animalids = {'151214', '151215', '151215', '151217', '151217'};
% blocks    = [5,        9,        10,        7,        8];
% animal    = [1,        2,        2,         3,        3];
% electrodes =[[1,16];  [1,16];   [1,16];    [1,16];   [1,16]];
% bldepth  =  [500,      500,      500,       500,      500];
% penangle =  [25,       25,       25,        25,       25];
% spacing =   [25,       25,       25,        25,       25];
% color =     [4,        4,        4,         4,        4];
% popfile = 'C:\Users\Julia\work\data\populations\control\intensity\intensitygreen_MUA.mat';

% % temp population 
% animalids = {'160210', '160210', '160210', '160210', '160210'};
% blocks    = [ 10,       11,       13,       14,       15];
% animal    = [ 1,        1,        1,        1,        1];
% electrodes =[[1,16];   [1,16];   [1,16];   [1,16];   [1,16]];
% bldepth  =  [ 500,      500,      500,      500,      500];
% penangle =  [ 25,       25,       25,       25,       25];
% spacing =   [ 25,       25,       25,       25,       25];
% color =     [ 4,        3,        4,        3,        4];
% popfile = 'C:\Users\Julia\work\data\populations\temppop.mat';

% temp population 
animalids = {'160114'};
blocks    = [ 17];
animal    = [ 1];
electrodes =[[1,16]];
bldepth  =  [ 500];
penangle =  [ 25];
spacing =   [ 25];
color =     [ 3];
popfile = 'C:\Users\Julia\work\data\populations\tempone.mat';


% 14, 20 and 40 interpolated for red
intensitiesred =   [5,    9,    14,   16,   20,   28,   40,    50, 10,   18,   32,  56,   100];
muwattsred       = [1.44, 2.69, 4.35, 4.84, 6.05, 8.56, 11.66, 15, 3.05, 5.44, 9.8, 16.7, 27.8];
% 14, 20,   28, 40 interpolated for green
intensitiesgreen = [2,    4,    6,    8,    10,   18,    32,   56,   100, 14, 20,   28, 40];
muwattsgreen     = [0.33, 1.48, 2.79, 4.18, 5.43, 11.47, 24.6, 50.7, 91.2, 9, 15.5, 23, 34.5];

depth = bldepth.*(cosd(penangle)*cosd(22));
spacing = spacing.*(cosd(penangle)*cosd(22));
evaldepth = 300;
for i = 1:length(depth)
    for j = 1:length(electrodes(i,1):electrodes(i,2))
        dm(i,j) = depth(i)-((j-1)*spacing(i)); % depth matrix
    end
    [c,di(i)] = min(abs(dm(i,:)-evaldepth)); % get depth index, least distance to evaldepth
end

recalculate_pop = 1;
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
        
        ch = di(blck);
        
        disp(['Block ' int2str(blck) '/' int2str(length(blocks)) '   channel ' int2str(ch) '/' int2str(length(electrodes(blck,1):electrodes(blck,2)))]);
              
        trialdur = result.stimduration*1000;
        msstamps = result.msstamps;
        if length(msstamps)~=length(result.light)
%             msstamps([193,291]) = []; % for 151023 block 5
%             msstamps([52, 98]) = []; % for 160210 block 13
%             msstamps([93]) = []; % for 160217 block 13
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
        large = find(result.gratingInfo.size == max(sizes) & result.light == 0);
        small = find(result.gratingInfo.size == min(sizes) & result.light == 0);
        lr1 = intersect(large,oktrials);
        sr1 = intersect(small,oktrials);
        for i = 1:length(large)
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
        for l = 1:length(lightlevels)
            for d = 1:length(sizes)
                
                thisinds = find(result.gratingInfo.size == sizes(d) & result.light == lightlevels(l));
                thisruninds = intersect(thisinds,oktrials); thisstillinds = intersect(thisinds,stilltrials);
                
                condresp(blck,l,d,:) = nanmean(resp(thisinds,:),1);
                condfr(blck,l,d) = nanmean(frs(thisinds));
                condmuacresp(blck,l,d,:) = nanmean(muacresp(thisinds,:));
                condmua(blck,l,d) = nanmean(muacmean(thisinds));
                condg1pow(blck,l,d) = nanmean(g1pow(thisinds));
                condg2pow(blck,l,d) = nanmean(g2pow(thisinds));
                
                condlfpresp(blck,l,d,:) = nanmean(lfpresp(thisinds,:));                
                [lfpspect(blck,l,d,:),chf,lfpspecterr(blck,l,d,:,:)] = mtspectrumc(squeeze(lfpresp(thisinds,1001:1800))',params);
                [muacspect(blck,l,d,:),chf,muacspecterr_c(blck,l,d,:,:)] = mtspectrumc(squeeze(muacresp(thisinds,1001:1800))',params);
                tbtlfpspect(blck,l,d,:) = nanmean(lfpspecttbt(thisinds,:),1);
                tbtlfpspecterr(blck,l,d,:) = nanstd(lfpspecttbt(thisinds,:),1,1)./sqrt(length(thisinds));
%                 [sfcoher(blck,l,d,:),a,b,c,de,sfcfx,e,f,g,sfcoherr_c(blck,l,d,:,:)] = coherencycpt(lfpresp(thisinds,1001:1800)',ptresp(thisinds),params);
                
                condg1powresp(blck,l,d,:) = mean(gamma1powresp(thisinds,:));
                condg2powresp(blck,l,d,:) = mean(gamma2powresp(thisinds,:));
                
                clear g1phasmatl; clear allg1phases;
                g1phasmat = gamma1phases(thisinds,respwin);
                allg1phases = squeeze(g1phasmat(find(g1phasmat)));
                if isrow(allg1phases), allg1phases = allg1phases'; end
                condg1phases{blck,l,d} = allg1phases;
                prefg1phase(blck,l,d) = circ_mean(allg1phases);
                g1lockpval(blck,l,d) = circ_rtest(allg1phases);
                g1plv(blck,l,d) = circ_r(allg1phases);
                g1ppc(blck,l,d) = ppc(allg1phases);
                [g1angdev(blck,l,d),g1cstd(blck,l,d)] = circ_std(allg1phases);
                
                clear g2phasmatl; clear allg2phases;
                g2phasmat = gamma2phases(thisinds,respwin);
                allg2phases = squeeze(g2phasmat(find(g2phasmat)));
                if isrow(allg2phases), allg2phases = allg2phases'; end
                condg2phases{blck,l,d} = allg2phases;
                prefg2phase(blck,l,d) = circ_mean(allg2phases);
                g2lockpval(blck,l,d) = circ_rtest(allg2phases);
                g2plv(blck,l,d) = circ_r(allg2phases);
                g2ppc(blck,l,d) = ppc(allg2phases);
                [g2angdev(blck,l,d),g2cstd(blck,l,d)] = circ_std(allg2phases);
                
                clear phasmat; clear allphases;
                for b = 1:size(bandphases,1)
                    phasmat = squeeze(bandphases(b,thisinds,respwin));
                    allphases = squeeze(phasmat(find(phasmat)));
                    if isrow(allphases), allphases = allphases'; end
                    condbandphases{blck,b,l,d} = allphases;
                    bandprefphase(blck,b,l,d) = circ_mean(allphases);
                    bandlockpval(blck,b,l,d) = circ_rtest(allphases);
                    bandplv(blck,b,l,d) = circ_r(allphases);
                    bandppc(blck,b,l,d) = ppc(allphases);
                    [bandangdev(blck,b,l,d),bandcstd(blck,b,l,d)] = circ_std(allphases);                    
                    
                    powmat = squeeze(bandpowers(b,thisinds,respwin));
                    allpowers = squeeze(powmat(find(powmat)));
                    if isrow(allpowers), allpowers = allpowers'; end
                    condbandpowers{blck,b,l,d} = allpowers;
                end
                
                gnspikes(blck,l,d) = length(allg1phases);
                
                clear r1phasmat; clear allr1phases; clear phasmat; clear allphases;
                if ~isempty(thisruninds)
                    r1condfr(blck,l,d) = nanmean(frs(thisruninds));
                    [r1_lfpspect(blck,l,d,:),chf,r1_lfpspecterr(blck,l,d,:,:)] = mtspectrumc(squeeze(lfpresp(thisruninds,1001:1800))',params);
%                     [r1_sfcoher(blck,l,d,:),a,b,c,de,sfcfx,e,f,g,r1_sfcoherr_c(blck,l,d,:,:)] = coherencycpt(lfpresp(thisruninds,1001:1800)',ptresp(thisruninds),params);
                    r1_tbtlfpspect(blck,l,d,:) = nanmean(lfpspecttbt(thisruninds,:),1);                    
                    r1_tbtlfpspecterr(blck,l,d,:) = nanstd(lfpspecttbt(thisruninds,:),1,1)./sqrt(length(thisruninds));
                    r1anovavec = [r1anovavec; lfpspecttbt(thisruninds,lgi(blck))];
                    lg = [lg; ones(length(thisruninds),1).*l];
                    sg = [sg; ones(length(thisruninds),1).*d];
                    r1_ntrials(blck,l,d) = length(thisruninds);
                    r1phasmat = gamma1phases(thisruninds,respwin);
                    allr1phases = squeeze(r1phasmat(find(r1phasmat)));
                    if isrow(allr1phases), allr1phases = allr1phases'; end
                    condr1phases{blck,l,d} = allr1phases;
                    prefr1phase(blck,l,d) = circ_mean(allr1phases);
                    r1lockpval(blck,l,d) = circ_rtest(allr1phases);
                    r1plv(blck,l,d) = circ_r(allr1phases);
                    r1ppc(blck,l,d) = ppc(allr1phases);
                    [r1angdev(blck,l,d),r1cstd(blck,l,d)] = circ_std(allr1phases);
                    for b = 1:size(bandphases,1)
                        phasmat = squeeze(bandphases(b,thisruninds,respwin));
                        allphases = squeeze(phasmat(find(phasmat)));
                        if isrow(allphases), allphases = allphases'; end
                        condr1bandphases{blck,b,l,d} = allphases;
                        prefr1bandphase(blck,b,l,d) = circ_mean(allphases);
                        r1bandlockpval(blck,b,l,d) = circ_rtest(allphases);
                        r1bandplv(blck,b,l,d) = circ_r(allphases);
                        r1bandppc(blck,b,l,d) = ppc(allphases);
                        [r1bandangdev(blck,b,l,d),r1bandcstd(blck,b,l,d)] = circ_std(allphases);
                        
                        powmat = squeeze(bandpowers(b,thisruninds,respwin));
                        allpowers = squeeze(powmat(find(powmat)));
                        if isrow(allpowers), allpowers = allpowers'; end
                        condr1bandpowers{blck,b,l,d} = allpowers;
                    end
                    r1nspikes(blck,l,d) = length(allr1phases);
                else
                    r1condfr(blck,l,d) = NaN;
                    r1_lfpspect(blck,l,d,:) = nan(1,length(lfpspect(blck,l,d,:)));
                    r1_tbtlfpspect(blck,l,d,:) = nan(1,length(lfpspect(blck,l,d,:)));                    
                    r1_tbtlfpspecterr(blck,l,d,:) = nan(1,length(lfpspect(blck,l,d,:))); 
                    r1_lfpspecterr(blck,l,d,:,:) = nan(2,length(lfpspect(blck,l,d,:)));
%                     r1_sfcoher(blck,l,d,:) = nan(1,length(sfcoher(blck,l,d,:)));
%                     r1_sfcohererr(blck,l,d,:,:) = nan(2,length(sfcoher(blck,l,d,:)));
                    r1_ntrials(blck,l,d) = 0;
                    condr1phases{blck,l,d} = NaN;
                    r1_phases{blck,l,d} = NaN;
                    prefr1phase(blck,l,d) = NaN;
                    r1lockpval(blck,l,d) = NaN;
                    r1plv(blck,l,d) = NaN;
                    r1ppc(blck,l,d) = NaN;
                    r1angdev(blck,l,d) = NaN;
                    r1cstd(blck,l,d) = NaN;
                    for b = 1:size(bandphases,1)
                        prefr1bandphase(blck,b,l,d) = NaN;
                        r1bandlockpval(blck,b,l,d) = NaN;
                        r1bandplv(blck,b,l,d) = NaN;
                        r1bandppc(blck,b,l,d) = NaN;
                        r1bandangdev(blck,b,l,d) = NaN;
                        r1bandcstd(blck,b,l,d) = NaN;
                        condr1bandphases{blck,b,l,d} = NaN;
                        condr1bandpowers{blck,b,l,d} = NaN;
                    end
                    r1nspikes(blck,l,d) = NaN;
                end
                
                clear r0phasmat; clear allr0phases; clear phasmat; clear allphases;
                if ~isempty(thisstillinds)
                    r0condfr(blck,l,d) = nanmean(frs(thisstillinds));
                    [r0_lfpspect(blck,l,d,:),chf,r0_lfpspecterr(blck,l,d,:,:)] = mtspectrumc(squeeze(lfpresp(thisstillinds,1001:1800))',params);
%                     [r0_sfcoher(blck,l,d,:),a,b,c,de,sfcfx,e,f,g,r0_sfcoherr_c(blck,l,d,:,:)] = coherencycpt(lfpresp(thisstillinds,1001:1800)',ptresp(thisstillinds),params);
                    r0_tbtlfpspect(blck,l,d,:) = nanmean(lfpspecttbt(thisstillinds,:),1);                    
                    r0_tbtlfpspecterr(blck,l,d,:) = nanstd(lfpspecttbt(thisstillinds,:),1,1)./sqrt(length(thisstillinds));
                    r0_ntrials(blck,l,d) = length(thisstillinds);
                    r0phasmat = gamma1phases(thisstillinds,respwin);
                    allr0phases = squeeze(r0phasmat(find(r0phasmat)));
                    if isrow(allr0phases), allr0phases = allr0phases'; end
                    condr0phases{blck,l,d} = allr0phases;
                    prefr0phase(blck,l,d) = circ_mean(allr0phases);
                    r0lockpval(blck,l,d) = circ_rtest(allr0phases);
                    r0plv(blck,l,d) = circ_r(allr0phases);
                    r0ppc(blck,l,d) = ppc(allr0phases);
                    [r0angdev(blck,l,d),r0cstd(blck,l,d)] = circ_std(allr0phases);
                    for b = 1:size(bandphases,1)
                        phasmat = squeeze(bandphases(b,thisstillinds,respwin));
                        allphases = squeeze(phasmat(find(phasmat)));
                        if isrow(allphases), allphases = allphases'; end
                        condr0bandphases{blck,b,l,d} = allphases;
                        prefr0bandphase(blck,b,l,d) = circ_mean(allphases);
                        r0bandlockpval(blck,b,l,d) = circ_rtest(allphases);
                        r0bandplv(blck,b,l,d) = circ_r(allphases);
                        r0bandppc(blck,b,l,d) = ppc(allphases);
                        [r0bandangdev(blck,b,l,d),r0bandcstd(blck,b,l,d)] = circ_std(allphases);
                        
                        powmat = squeeze(bandpowers(b,thisstillinds,respwin));
                        allpowers = squeeze(powmat(find(powmat)));
                        if isrow(allpowers), allpowers = allpowers'; end
                        condr0bandpowers{blck,b,l,d} = allpowers;
                    end
                    r0nspikes(blck,l,d) = length(allr0phases);
                else
                    r0condfr(blck,l,d) = NaN;
                    r0_lfpspect(blck,l,d,:) = nan(1,length(lfpspect(blck,l,d,:)));
                    r0_lfpspecterr(blck,l,d,:,:) = nan(2,length(lfpspect(blck,l,d,:)));
                    r0_tbtlfpspect(blck,l,d,:) = nan(1,length(lfpspect(blck,l,d,:)));                    
                    r0_tbtlfpspecterr(blck,l,d,:) = nan(1,length(lfpspect(blck,l,d,:))); 
%                     r0_sfcoher(blck,l,d,:) = nan(1,length(sfcoher(blck,l,d,:)));
%                     r0_sfcohererr(blck,l,d,:,:) = nan(2,length(sfcoher(blck,l,d,:)));
                    r0_ntrials(blck,l,d) = 0;
                    condr0phases{blck,l,d} = NaN;
                    r0_phases{blck,l,d} = NaN;
                    prefr0phase(blck,l,d) = NaN;
                    r0lockpval(blck,l,d) = NaN;
                    r0plv(blck,l,d) = NaN;
                    r0ppc(blck,l,d) = NaN;
                    r0angdev(blck,l,d) = NaN;
                    r0cstd(blck,l,d) = NaN;
                    for b = 1:size(bandphases,1)
                        prefr0bandphase(blck,b,l,d) = NaN;
                        r0bandlockpval(blck,b,l,d) = NaN;
                        r0bandplv(blck,b,l,d) = NaN;
                        r0bandppc(blck,b,l,d) = NaN;
                        r0bandangdev(blck,b,l,d) = NaN;
                        r0bandcstd(blck,b,l,d) = NaN;
                        condr0bandphases{blck,b,l,d} = NaN;
                        condr0bandpowers{blck,b,l,d} = NaN;
                    end
                    r0nspikes(blck,l,d) = NaN;
                end
            end
        end
%         if ~isempty(r1anovavec)
%             [anovap(blck,:),table,stats] = anovan(r1anovavec,{lg,sg}); % 1 = light, 2 = stimulus
%         end
    disp('')
    end
    save(popfile, '-v7.3');
else
    load(popfile);
end

for i = 1:length(blocks)
    lgpow(i,:,:) = squeeze(lfpspect(i,:,:,lgi(i))); % at second contact only
    r1lgpow(i,:,:) = squeeze(r1_tbtlfpspect(i,:,:,lgi(i)));
    r1lgpowerr(i,:,:) = squeeze(r1_tbtlfpspecterr(i,:,:,lgi(i)));
    rellgpow(i,:) = squeeze(r1_tbtlfpspect(i,:,2,lgi(i)))./squeeze(nanmean(r1_tbtlfpspect(i,:,2,1:104),4));
    normlgpow(i,:) = r1lgpow(i,:,2)./r1lgpow(i,1,2);
    normrellgpow(i,:) = rellgpow(i,:)./rellgpow(i,1);
    r1fr(i,:) = r1condfr(i,:,2);
    normfr(i,:) = r1fr(i,:)./r1fr(i,1);
    for l = 1:6
        for d = 1:2
            r1tbtfillspecy(i,l,d,:) = [(squeeze(r1_tbtlfpspect(i,l,d,1:104))+squeeze(r1_tbtlfpspecterr(i,l,d,1:104)))',fliplr((squeeze(r1_tbtlfpspect(i,l,d,1:104))-squeeze(r1_tbtlfpspecterr(i,l,d,1:104)))')];
            tbtfillspecy(i,l,d,:) = [(squeeze(tbtlfpspect(i,l,d,1:104))+squeeze(tbtlfpspecterr(i,l,d,1:104)))',fliplr((squeeze(tbtlfpspect(i,l,d,1:104))-squeeze(tbtlfpspecterr(i,l,d,1:104)))')];
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

figure
hold on
for i = 1:size(watts,1)
    if color(i) == 3, pc = 'r'; else pc = 'g'; end
    plot(watts(i,:),rellgpow(i,:),'g')
end
title('normalized gamma power')
ylabel('gamma power')
xlabel('light intensity')

errorbar(watts(1,1:5),nanmean(normlgpow(:,1:5),1),nanstd(normlgpow(:,1:5),1,1),'go-','markersize',3,'markerfacecolor','g','linewidth',2)
errorbar(watts(1,1:5),nanmean(normlgpow(:,1:5),1),nanstd(normlgpow(:,1:5),1,1),'bo-','markersize',3,'markerfacecolor','b','linewidth',2)
errorbar(watts(1,1:5),nanmean(normlgpow(:,1:5),1),nanstd(normlgpow(:,1:5),1,1),'ro-','markersize',3,'markerfacecolor','r','linewidth',2)
errorbar(watts(1,1:5),nanmean(normlgpow(:,1:5),1),nanstd(normlgpow(:,1:5),1,1),'ko-','markersize',3,'markerfacecolor','k','linewidth',2)

errorbar(watts(1,1:5),nanmean(normrellgpow(:,1:5),1),nanstd(normrellgpow(:,1:5),1,1),'ro-','markersize',3,'markerfacecolor','r','linewidth',2)
errorbar(watts(1,1:5),nanmean(normrellgpow(:,1:5),1),nanstd(normrellgpow(:,1:5),1,1),'bo-','markersize',3,'markerfacecolor','b','linewidth',2)
errorbar(watts(1,1:5),nanmean(normrellgpow(:,1:5),1),nanstd(normrellgpow(:,1:5),1,1),'go-','markersize',3,'markerfacecolor','g','linewidth',2)
errorbar(watts(1,1:5),nanmean(normrellgpow(:,1:5),1),nanstd(normrellgpow(:,1:5),1,1),'ko-','markersize',3,'markerfacecolor','k','linewidth',2)

legend('PV eArch red n = 3','PV Halo n = 2', 'SOM Halo n = 7','control n = 6','location','nw')
legend('SOM Halo n = 9', 'PV Halo n = 6', 'PV eArch red n = 5','control n = 6','location','nw')
legend('control n = 6','PV eArch red n = 3','PV Halo high n = 2','PV Halo low n = 2', 'SOM Halo n = 7','location','sw')
xlabel('light intensity [mW]')
ylabel('normalized gamma power')