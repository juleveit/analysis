function MUA_population_twogratingcohere


% % PV Halo population
% animalids =    { '150818', '150820', '150823', '150824'};
% blocks    =     [7,        12,       11.        9];
% animal    =     [1,         2,       3,         4];
% centelectrodes =[[1,16];   [1,16];   [1,16];    [17,32]];
% centdepth    =  [400,      400,      400,       400];
% surrelectrodes =[[17,32];  [17,32];  [17,32];   [1,16]];
% surrdepth    =  [400,      400,      400,       400];
% penangle =      [25,       25,       25,        25];
% spacing =       [25,       25,       25,        25];
% popfile = 'C:\Users\Julia\work\data\populations\PV_Halo\twogratingcoher\MUA_population.mat';

% % SOM Halo population
% animalids =    {'150825', '150831', '150902', '150907', '150909', '150915', '151022', '151023', '151027', '151109', '151110', '151127'};
% blocks    =     [10,        9,       4,        4,        3,        4,        9,        8,        9,        4,        7,        2];
% animal    =     [1,         2,       3,        4,        5,        6,        7,        8,        9,        10,       11,       12];
% centelectrodes =[[17,32];   [1,16];  [1,16];   [1,16];   [1,16];   [17,32];  [1,16];   [17,32];  [17,32];  [17,32];  [1,16];   [1,16]];
% centdepth    =  [400,       400,     400,      400,      400,      400,      400,      400,      500,      500,      500,      500];
% surrelectrodes =[[1,16];    [17,32]; [17,32];  [17,32];  [17,32];  [1,16];   [17,32];  [1,16];   [1,16];   [1,16];   [17,32];  [17,32]];
% surrdepth    =  [400,       400,     400,      400,      400,      400,      400,      400,      500,      500,      500,      500];
% penangle =      [25,        25,      25,       25,       25,       25,       25,       25,       25,       25,       25,       25];
% spacing =       [25,        25,      25,       25,       25,       25,       25,       25,       25,       25,       25,       25];
% popfile = 'C:\Users\Julia\work\data\populations\SOM_Halo_later\twogratingcoher\MUA_population.mat';

% % SOM PV population
% animalids =    {'150825', '150831', '150902', '150907', '150909', '150915', '151022', '151023', '151027', '151109', '151110', '150818', '150820', '150823', '150824'};
% blocks    =     [10,        9,       4,        4,        3,        4,        9,        8,        9,        4,        7,        7,        12,       11.       9];
% animal    =     [1,         2,       3,        4,        5,        6,        7,        8,        9,        10,       11,       12,       13,       14,       15];
% centelectrodes =[[17,32];   [1,16];  [1,16];   [1,16];   [1,16];   [17,32];  [1,16];   [17,32];  [17,32];  [17,32];  [1,16];  [1,16];   [1,16];   [1,16];   [17,32]];
% centdepth    =  [400,       400,     400,      400,      400,      400,      400,      400,      500,      500,      500,      400,      400,      400,      400];
% surrelectrodes =[[1,16];    [17,32]; [17,32];  [17,32];  [17,32];  [1,16];   [17,32];  [1,16];   [1,16];   [1,16];   [17,32]; [17,32];  [17,32];  [17,32];  [1,16]];
% surrdepth    =  [400,       400,     400,      400,      400,      400,      400,      400,      500,      500,      500,      400,      400,      400,      400];
% penangle =      [25,        25,      25,       25,       25,       25,       25,       25,       25,       25,       25,       25,       25,       25,       25];
% spacing =       [25,        25,      25,       25,       25,       25,       25,       25,       25,       25,       25,       25,       25,       25,       25];
% popfile = 'C:\Users\Julia\work\data\populations\SOMPVcombined\twogratingcoher\MUA_population.mat';

% animalids =    { '151204', '151204', '151204'};
% blocks    =     [5,         6,        7];
% animal    =     [1,         2,        3];
% centelectrodes =[[1,16];   [1,16];   [1,16]];
% centdepth    =  [500,       500,      500];
% surrelectrodes =[[17,32];  [17,32];   [17,32]];
% surrdepth    =  [500,       500,      500];
% penangle =      [25,        25,       25];
% spacing =       [25,        25,       25];
% popfile = 'C:\Users\Julia\work\data\populations\testpop.mat';
% % 
% animalids =    { '151215', '151215'};
% blocks    =     [3,         4];
% animal    =     [1,         2];
% centelectrodes =[[1,16];   [1,16]];
% centdepth    =  [500,       500];
% surrelectrodes =[[17,32];  [17,32]];
% surrdepth    =  [500,       500];
% penangle =      [25,        25];
% spacing =       [25,        25];
% popfile = 'C:\Users\Julia\work\data\populations\testpop2.mat';

% animalids =    { '151222', '151222', '151222', '151222'};
% blocks    =     [4,         5,        6,        7];
% animal    =     [1,         1,        1,        1];
% centelectrodes =[[1,16];   [1,16];   [1,16];   [1,16]];
% centdepth    =  [600,       600,      600,      600];
% surrelectrodes =[[17,32];  [17,32];  [17,32];  [17,32]];
% surrdepth    =  [500,       500,      500,      500];
% penangle =      [25,        25,       25,       25];
% spacing =       [25,        25,       25,       25];
% popfile = 'C:\Users\Julia\work\data\populations\testpop3.mat';


% animalids =    { '151223', '151223', '151223', '151223'};
% blocks    =     [5,         6,        7         8];
% animal    =     [1,         1,        1,        1];
% centelectrodes =[[1,16];   [1,16];   [1,16];   [1,16]];
% centdepth    =  [500,       500,      500,      500];
% surrelectrodes =[[17,32];  [17,32];  [17,32];  [17,32]];
% surrdepth    =  [500,       500,      500,      500];
% penangle =      [25,        25,       25,       25];
% spacing =       [25,        25,       25,       25];
% popfile = 'C:\Users\Julia\work\data\populations\testpop4.mat';


% animalids =    { '160107', '160107'};
% blocks    =     [4,         5];
% animal    =     [1,         2];
% centelectrodes =[[1,16];   [1,16]];
% centdepth    =  [400,       400];
% surrelectrodes =[[17,32];  [17,32]];
% surrdepth    =  [500,       500];
% penangle =      [25,        25];
% spacing =       [25,        25];
% popfile = 'C:\Users\Julia\work\data\populations\testpop5.mat';


% animalids =    { '160113', '160113'};
% blocks    =     [3,         4];
% animal    =     [1,         2];
% centelectrodes =[[1,16];   [1,16]];
% centdepth    =  [500,       500];
% surrelectrodes =[[17,32];  [17,32]];
% surrdepth    =  [500,       500];
% penangle =      [25,        25];
% spacing =       [25,        25];
% popfile = 'C:\Users\Julia\work\data\populations\testpop6.mat';


% animalids =    {'160114',  '160114'};
% blocks    =     [6,         5];
% animal    =     [1,         2];
% centelectrodes =[[1,16];   [1,16]];
% centdepth    =  [500,       500];
% surrelectrodes =[[17,32];  [17,32]];
% surrdepth    =  [500,       500];
% penangle =      [25,        25];
% spacing =       [25,        25];
% popfile = 'C:\Users\Julia\work\data\populations\testpop7.mat';
% 
% % 
% animalids =    { '160115', '160115'};
% blocks    =     [5,         4];
% animal    =     [1,         2];
% centelectrodes =[[1,16];   [1,16]];
% centdepth    =  [500,       500];
% surrelectrodes =[[17,32];  [17,32]];
% surrdepth    =  [500,       500];
% penangle =      [25,        25];
% spacing =       [25,        25];
% popfile = 'C:\Users\Julia\work\data\populations\testpop8.mat';

% animalids =    { '160204', '160204'};
% blocks    =     [4,         5];
% animal    =     [1,         2];
% centelectrodes =[[1,16];   [1,16]];
% centdepth    =  [500,       500];
% surrelectrodes =[[17,32];  [17,32]];
% surrdepth    =  [500,       500];
% penangle =      [25,        25];
% spacing =       [25,        25];
% popfile = 'C:\Users\Julia\work\data\populations\testpop9.mat';

% animalids =    { '151127', '151127'};    % 3 is vertical, 2 is horizontal as usual
% blocks    =     [2,         3];
% animal    =     [1,         2];
% centelectrodes =[[1,16];   [1,16]];
% centdepth    =  [500,       500];
% surrelectrodes =[[17,32];  [17,32]];
% surrdepth    =  [500,       500];
% penangle =      [25,        25];
% spacing =       [25,        25];
% popfile = 'C:\Users\Julia\work\data\populations\testpop10.mat';

% animalids =    { '160310', '160310'};    
% blocks    =     [11,         10];
% animal    =     [1,         2];
% centelectrodes =[[1,16];   [1,16]];
% centdepth    =  [500,       500];
% surrelectrodes =[[17,32];  [17,32]];
% surrdepth    =  [500,       500];
% penangle =      [25,        25];
% spacing =       [25,        25];
% popfile = 'C:\Users\Julia\work\data\populations\testpop11.mat';

animalids =    { '160311'};    
blocks    =     [6];
animal    =     [1];
centelectrodes =[[1,16]];
centdepth    =  [500];
surrelectrodes =[[17,32]];
surrdepth    =  [500];
penangle =      [25];
spacing =       [25];
popfile = 'C:\Users\Julia\work\data\populations\testpop12.mat';


recalculate_pop = 0;
recalculate_muafile = 0;

% chronux parameters
params.tapers = [5,9]; params.Fs = 1000; params.err = [2, 0.05]; params.trialave = 1;


if ~exist(popfile) || recalculate_pop

    cll = 1;
    for blck = 1:length(blocks)
        
        basepath = strcat('C:\Users\Julia\work\data\', animalids{blck}, '\');
        centfile = strcat(basepath, 'muaresult_', int2str(blocks(blck)), '_', int2str(centelectrodes(blck,1)), '-', int2str(centelectrodes(blck,2)), '.mat');
        surrfile = strcat(basepath, 'muaresult_', int2str(blocks(blck)), '_', int2str(surrelectrodes(blck,1)), '-', int2str(surrelectrodes(blck,2)), '.mat');
        oldcentfile = strcat(basepath, 'muaresult_', int2str(blocks(blck)), '_', int2str(centelectrodes(blck,1)), ':', int2str(centelectrodes(blck,2)), '.mat');
        oldsurrfile = strcat(basepath, 'muaresult_', int2str(blocks(blck)), '_', int2str(surrelectrodes(blck,1)), ':', int2str(surrelectrodes(blck,2)), '.mat');
        
        if ~exist(centfile) || recalculate_muafile
            if ~exist(oldcentfile)
                centresult = MUAdataprepare(basepath,animalids{blck},blocks(blck),centelectrodes(blck,1):centelectrodes(blck,2));
                save(centfile,'centresult')
            else
                load(oldcentfile);
                save(centfile, 'centresult');
            end
        else
            load(centfile);
        end        
        if ~exist(surrfile) || recalculate_muafile
            if ~exist(oldsurrfile)
                surrresult = MUAdataprepare(basepath,animalids{blck},blocks(blck),surrelectrodes(blck,1):surrelectrodes(blck,2));
                save(surrfile,'surrresult')
            else
                load(oldsurrfile);
                save(surrfile, 'surrresult');
            end
        else
            load(surrfile);
        end
        
        prestim = 300;
        poststim = 700;
        respwin = 501:1495; % after stimulus onset
        offsetwin = 1501:2500;
        respwin = respwin+prestim;
        
        for ch = 1: length(centelectrodes(blck,1):centelectrodes(blck,2))
            
            disp(['Block ' int2str(blck) '/' int2str(length(blocks)) '   channel ' int2str(ch) '/' int2str(length(centelectrodes(blck,1):centelectrodes(blck,2)))]);
            
            depthc(ch) = centdepth(blck)-(ch-1)*spacing(blck);
            depths(ch) = surrdepth(blck)-(ch-1)*spacing(blck);    
            
            trialdur = centresult.stimduration*1000;
            msstamps = centresult.msstamps;
            if length(msstamps)~=length(centresult.light)
% %                 msstamps([193,291]) = []; % for 151023 block 5
%                     msstamps(213) = []; % 160114 block 5
%                 centresult.msstamps = msstamps;
%                 surrresult.msstamps = msstamps;
%                 save(surrfile,'surrresult');
%                 save(centfile,'centresult');
                pause;
            end
            
            for i = 1:length(msstamps)
                speed(i,:) = centresult.runspeed(msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim);
            end
            
            % figure out sufficiently high and nonvariable runspeed trials
            meanspeed = mean(speed(:,respwin),2);
            stdspeed = std(speed(:,respwin),1,2);
            notstill = find(meanspeed>1);
            okspeed = find(meanspeed>( mean(meanspeed(notstill))-(1.5*std(meanspeed(notstill))) ) & meanspeed>1 );
            okvar = find(stdspeed<( mean(stdspeed(notstill))+(1.5*std(stdspeed(notstill)))) & stdspeed>.5);
            oktrials = intersect(okspeed,okvar);
            nonoktrials = 1:size(speed,1); nonoktrials(oktrials) = [];
            stilltrials = 1:size(speed,1); stilltrials(notstill) = [];
            
            % find low and high gamma peaks in spectra            
            oris1 = unique(centresult.gratingInfo.Orientation1);
            oris2 = unique(centresult.gratingInfo.Orientation2);
            clevels1 = unique(centresult.gratingInfo.Contrast1);
            clevels2 = unique(centresult.gratingInfo.Contrast2);
            cmax = max(clevels1);
            beta = [15,40];
            gamma = [50,70];
            large = find(centresult.gratingInfo.Contrast1 == cmax & centresult.gratingInfo.Contrast2 == cmax & centresult.light == 0);
            small = find(centresult.gratingInfo.Contrast1 == 0 & centresult.gratingInfo.Contrast2 == 0 & centresult.light == 0);
            largel1 = find(centresult.gratingInfo.Contrast1 == cmax & centresult.gratingInfo.Contrast2 == cmax & centresult.light == 1);
            lr1 = intersect(large,oktrials);
            sr1 = intersect(small,oktrials);
            lr1l1 = intersect(largel1,oktrials);
            for i = 1:length(large)
                [pl(i,:),f] = pmtm(centresult.lfp(ch,centresult.msstamps(large(i))+300:centresult.msstamps(large(i))+1300),3,[],1000);
                [ps(i,:),f] = pmtm(centresult.lfp(ch,centresult.msstamps(small(i))+300:centresult.msstamps(small(i))+1300),3,[],1000);
            end 
            plr1 = nan(1,size(pl,2)); psr1 = nan(1,size(pl,2));
            if length(lr1>=5)
                for i = 1:length(lr1)
                    [plr1(i,:),f] = mtspectrumc(centresult.lfp(ch,centresult.msstamps(lr1(i))+700:centresult.msstamps(lr1(i))+1500),params);
%                     [plr1(i,:),f] = pmtm(centresult.lfp(ch,centresult.msstamps(lr1(i))+300:centresult.msstamps(lr1(i))+1300),3,[],1000);
                end
            end
            if length(sr1>=5)
                for i = 1:length(sr1)
                    [psr1(i,:),f] = mtspectrumc(centresult.lfp(ch,centresult.msstamps(sr1(i))+700:centresult.msstamps(sr1(i))+1500),params);
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
                g1peak(blck,ch) = 0;
            else
                peaks = find(diff(bsig)>0)+1;
                pvs = bsig(peaks);
                bpi = peaks(pvs == max(pvs));
                bpi = bpi+b1-1;
                g1peak(blck,ch) = 1;
            end
            if isempty(find(diff(gsig)>0)) % there is no clear beta peak
                gpi = round((g1+g2)/2);
                g2peak(blck,ch) = 0;
            else
                peaks = find(diff(gsig)>0)+1;
                pvs = gsig(peaks);
                gpi = peaks(pvs == max(pvs));
                gpi = gpi+g1-1;
                g2peak(blck,ch) = 1;
            end
            lgi(blck,ch) = bpi; hgi(blck,ch) = gpi;

            msStimes_c = round(centresult.msStimes{ch});
            msStimes_s = round(surrresult.msStimes{ch});
            if ~isempty(msStimes_c) & msStimes_c(1) == 0, msStimes_c(1) = 1; end
            if ~isempty(msStimes_s) & msStimes_s(1) == 0, msStimes_s(1) = 1; end

            chan_c = zeros(1,size(centresult.lfp,2));
            chan_s = zeros(1,size(surrresult.lfp,2));
            chan_c(msStimes_c) = 1;
            chan_s(msStimes_s) = 1;
    
            for i = 1:length(msstamps)
                resp_c(i,:) = chan_c(msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim); 
                resp_s(i,:) = chan_s(msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim);  
                [xc(i,:),lags] = xcorr(resp_c(i,respwin),resp_s(i,respwin),200,'coeff');  
                hc = xcorr(resp_c(i,respwin),resp_s(i,respwin),10);
                corstrength(i) = sum(hc)/sqrt((sum(resp_c(i,respwin))^2+sum(resp_s(i,respwin))^2)/2);
                hh = find(resp_c(i,1001:1800))'; ptresp_c(i).times = hh./1000;                
                hh = find(resp_s(i,1001:1800))'; ptresp_s(i).times = hh./1000;
                lfpresp_c(i,:) = centresult.lfp(ch,msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim);
                lfpresp_s(i,:) = surrresult.lfp(ch,msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim);
                lfpspecttbt_c(i,:) = mtspectrumc(squeeze(lfpresp_c(i,1001:1800))',params);
                lfpspecttbt_s(i,:) = mtspectrumc(squeeze(lfpresp_s(i,1001:1800))',params);
                tbtlfpcoher(i,:) = coherencyc(lfpresp_c(i,respwin)',lfpresp_s(i,respwin)',params);
                muacresp_c(i,:) =  centresult.muac(ch,msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim).^2;
                muacresp_s(i,:) =  surrresult.muac(ch,msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim).^2;
                [mcxc(i,:),lags] = xcorr(muacresp_c(i,respwin),muacresp_s(i,respwin),200,'coeff');   
                muaac_c(i,:) = xcorr(resp_c(i,respwin),resp_c(i,respwin),200,'coeff');
                muaac_s(i,:) = xcorr(resp_s(i,respwin),resp_s(i,respwin),200,'coeff');
                spksc = find(resp_c(i,respwin)); spksc = spksc+respwin(1)-1;
                spkss = find(resp_s(i,respwin)); spkss = spkss+respwin(1)-1;
                clear snips_s; clear snipspect_s; 
                clear snips_c; clear snipspect_c;
                if ~isempty(spksc)
                    for s = 1:length(spksc)
                        snips_cs(s,:) = lfpresp_s(i,spksc(s)-250:spksc(s)+250);
                        snipspect_cs(s,:) = mtspectrumc(snips_cs(s,:)',params);
                        
                        snips_c(s,:) = lfpresp_c(i,spksc(s)-250:spksc(s)+250);
                        snipspect_c(s,:) = mtspectrumc(snips_c(s,:)',params);
                    end                    
                    nspikes_c(i) = size(snips_cs,1);
                    sta_cs(i,:) = nanmean(snips_cs)./nspikes_c(i);
                    [staspect_cs(i,:),staf] = mtspectrumc(sta_cs(i,:)',params);
                    allsnipspect_cs(i,:) = nanmean(snipspect_cs,1);
                    sfc_cs(i,:) = staspect_cs(i,:)./allsnipspect_cs(i,:);
                    
                    sta_c(i,:) = nanmean(snips_c)./nspikes_c(i);
                    [staspect_c(i,:),staf] = mtspectrumc(sta_c(i,:)',params);
                    allsnipspect_c(i,:) = nanmean(snipspect_c,1);
                    sfc_c(i,:) = staspect_c(i,:)./allsnipspect_c(i,:);
                else
                    snips_cs = nan(1,501);
                    snipspect_cs = nan(1,257);
                    nspikes_c(i) = 0;
                    sta_cs(i,:) = nan(1,501);
                    staspect_cs(i,:) = nan(1,257);
                    allsnipspect_cs(i,:) = nan(1,257);
                    sfc_cs(i,:) = nan(1,257);                    
                    
                    snips_c = nan(1,501);
                    snipspect_c = nan(1,257);
                    sta_c(i,:) = nan(1,501);
                    staspect_c(i,:) = nan(1,257);
                    allsnipspect_c(i,:) = nan(1,257);
                    sfc_c(i,:) = nan(1,257);
                end
                
                if ~isempty(spkss)
                    for s = 1:length(spkss)
                        snips_sc(s,:) = lfpresp_c(i,spkss(s)-250:spkss(s)+250);
                        snipspect_sc(s,:) = mtspectrumc(snips_sc(s,:)',params);
                        
                        snips_s(s,:) = lfpresp_c(i,spkss(s)-250:spkss(s)+250);
                        snipspect_s(s,:) = mtspectrumc(snips_s(s,:)',params);
                    end                    
                    nspikes_s(i) = size(snips_sc,1);
                    sta_sc(i,:) = nanmean(snips_sc)./nspikes_s(i);
                    [staspect_sc(i,:),staf] = mtspectrumc(sta_sc(i,:)',params);
                    allsnipspect_sc(i,:) = nanmean(snipspect_sc,1);
                    sfc_sc(i,:) = staspect_sc(i,:)./allsnipspect_sc(i,:);
                    
                    sta_s(i,:) = nanmean(snips_s)./nspikes_s(i);
                    [staspect_s(i,:),staf] = mtspectrumc(sta_s(i,:)',params);
                    allsnipspect_s(i,:) = nanmean(snipspect_s,1);
                    sfc_s(i,:) = staspect_s(i,:)./allsnipspect_s(i,:);
                else
                    snips_sc = nan(1,501);
                    snipspect_sc = nan(1,257);
                    nspikes_s(i) = 0;
                    sta_sc(i,:) = nan(1,501);
                    staspect_sc(i,:) = nan(1,257);
                    allsnipspect_sc(i,:) = nan(1,257);
                    sfc_sc(i,:) = nan(1,257);                    
                    
                    snips_s = nan(1,501);
                    snipspect_s = nan(1,257);
                    sta_s(i,:) = nan(1,501);
                    staspect_s(i,:) = nan(1,257);
                    allsnipspect_s(i,:) = nan(1,257);
                    sfc_s(i,:) = nan(1,257);
                end
            end
            
            [leffect(blck,ch),pleffect(blck,ch)] = ttest2(lfpspecttbt_c(lr1,bpi),lfpspecttbt_c(lr1l1,bpi));
            
            frs_c = sum(resp_c(:,respwin),2)./(length(respwin)/1000);
            frs_s = sum(resp_s(:,respwin),2)./(length(respwin)/1000);
            bl_c = sum(resp_c(:,1:prestim),2)./(prestim/1000);
            bl_s = sum(resp_s(:,1:prestim),2)./(prestim/1000);
            muacmean_c = mean(muacresp_c(:,respwin),2);
            muacmean_s = mean(muacresp_s(:,respwin),2);
            
            ordf = [0,90,-1,-2];  % 1: same, 2: different, 3: cent only, 4: surround only
            for l = 1:length(unique(centresult.light))
                for d = 1:length(ordf)
                    if ordf(d) == -1
                        thisinds = find( centresult.gratingInfo.Contrast1 == cmax &...
                            centresult.gratingInfo.Contrast2 == 0 &...
                            centresult.light == l-1);                        
                    elseif ordf(d) == -2
                        thisinds = find( centresult.gratingInfo.Contrast1 == 0 &...
                            centresult.gratingInfo.Contrast2 == cmax &...
                            centresult.light == l-1);
                    elseif ordf(d) == 0
                        thisinds = find(centresult.gratingInfo.Orientation1 == ...
                            centresult.gratingInfo.Orientation2 &...
                            centresult.gratingInfo.Contrast1 == cmax &...
                            centresult.gratingInfo.Contrast2 == cmax &...
                            centresult.light == l-1);
                    else
                        thisinds = find(centresult.gratingInfo.Orientation1 ~= ...
                            centresult.gratingInfo.Orientation2 &...
                            centresult.gratingInfo.Contrast1 == cmax &...
                            centresult.gratingInfo.Contrast2 == cmax &...
                            centresult.light == l-1);
                    end             
                    
                    thisruninds = intersect(thisinds,oktrials); thisstillinds = intersect(thisinds,stilltrials);
                  
                    condresp_c(blck,ch,l,d,:) = nanmean(resp_c(thisinds,:),1);
                    condresp_s(blck,ch,l,d,:) = nanmean(resp_s(thisinds,:),1);
                    condfr_c(blck,ch,l,d) = nanmean(frs_c(thisinds),1);
                    condfr_s(blck,ch,l,d) = nanmean(frs_s(thisinds),1);
                    condmuacresp_c(blck,ch,l,d,:) = nanmean(muacresp_c(thisinds,:),1);
                    condmuacresp_s(blck,ch,l,d,:) = nanmean(muacresp_s(thisinds,:),1);
                    condmua_c(blck,ch,l,d) = nanmean(muacmean_c(thisinds),1);
                    condmua_s(blck,ch,l,d) = nanmean(muacmean_s(thisinds),1);
                    condlfpresp_c(blck,ch,l,d,:) = nanmean(lfpresp_c(thisinds,:),1);
                    condlfpresp_s(blck,ch,l,d,:) = nanmean(lfpresp_s(thisinds,:),1);
                    
                    [lfpspect_c(blck,ch,l,d,:),chf,lfpspecterr_c(blck,ch,l,d,:,:)] = mtspectrumc(squeeze(lfpresp_c(thisinds,1001:1800))',params);
                    [lfpspect_s(blck,ch,l,d,:),chf,lfpspecterr_s(blck,ch,l,d,:,:)] = mtspectrumc(squeeze(lfpresp_s(thisinds,1001:1800))',params);
                    tbtlfpspect_c(blck,ch,l,d,:) = nanmean(lfpspecttbt_c(thisinds,:),1);
                    tbtlfpspect_s(blck,ch,l,d,:) = nanmean(lfpspecttbt_s(thisinds,:),1);
                    tbtlfpspecterr_c(blck,ch,l,d,:) = nanstd(lfpspecttbt_c(thisinds,:),1,1)./sqrt(length(thisinds));
                    tbtlfpspecterr_s(blck,ch,l,d,:) = nanstd(lfpspecttbt_s(thisinds,:),1,1)./sqrt(length(thisinds));
                    [muacspect_c(blck,ch,l,d,:),chf,muacspecterr_c(blck,ch,l,d,:,:)] = mtspectrumc(squeeze(muacresp_c(thisinds,1001:1800))',params);
                    [muacspect_s(blck,ch,l,d,:),chf,muacspecterr_s(blck,ch,l,d,:,:)] = mtspectrumc(squeeze(muacresp_s(thisinds,1001:1800))',params);
                    
                    [lfpcoherence(blck,ch,l,d,:),a,b,c,de,cfx,e,f,lfpcoherr(blck,ch,l,d,:,:)] = coherencyc(lfpresp_c(thisinds,respwin)',lfpresp_s(thisinds,respwin)',params);      
                    condtbtlfpcoherence(blck,ch,l,d,:) = nanmean(tbtlfpcoher(thisinds,:),1);
                    condtbtlfpcoherr(blck,ch,l,d,:) = nanstd(tbtlfpcoher(thisinds,:),1,1)./sqrt(length(thisinds));
                    condxc(blck,ch,l,d,:) = squeeze(mean(xc(thisinds,:),1));
                    condmcxc(blck,ch,l,d,:) = squeeze(mean(mcxc(thisinds,:),1));
                    condmuaac_c(blck,ch,l,d,:) = squeeze(nanmean(muaac_c(thisinds,:),1));
                    condmuaac_s(blck,ch,l,d,:) = squeeze(nanmean(muaac_s(thisinds,:),1));
%                     [sfcoher_c(blck,ch,l,d,:),a,b,c,de,sfcfx,e,f,g,sfcoherr_c(blck,ch,l,d,:,:)] = coherencycpt(lfpresp_c(thisinds,1001:1800)',ptresp_c(thisinds),params);
%                     [sfcoher_s(blck,ch,l,d,:),a,b,c,de,sfcfx,e,f,g,sfcoherr_s(blck,ch,l,d,:,:)] = coherencycpt(lfpresp_s(thisinds,1001:1800)',ptresp_s(thisinds),params);
                    condcorstrength(blck,ch,l,d,:) = mean(corstrength(thisinds));
                    
                    if ~isempty(thisruninds)
                        [r1_lfpspect_c(blck,ch,l,d,:),chf,r1_lfpspecterr_c(blck,ch,l,d,:,:)] = mtspectrumc(squeeze(lfpresp_c(thisruninds,1001:1800))',params);
                        [r1_lfpspect_s(blck,ch,l,d,:),chf,r1_lfpspecterr_s(blck,ch,l,d,:,:)] = mtspectrumc(squeeze(lfpresp_s(thisruninds,1001:1800))',params);
                        r1_tbtlfpspect_c(blck,ch,l,d,:) = nanmean(lfpspecttbt_c(thisruninds,:),1);
                        r1_tbtlfpspect_s(blck,ch,l,d,:) = nanmean(lfpspecttbt_s(thisruninds,:),1);
                        r1_tbtlfpspecterr_c(blck,ch,l,d,:) = nanstd(lfpspecttbt_c(thisruninds,:),1,1)./sqrt(length(thisruninds));
                        r1_tbtlfpspecterr_s(blck,ch,l,d,:) = nanstd(lfpspecttbt_s(thisruninds,:),1,1)./sqrt(length(thisruninds));
                        [r1_lfpcoherence(blck,ch,l,d,:),a,b,c,de,cfx,e,f,r1_lfpcoherr(blck,ch,l,d,:,:)] = coherencyc(lfpresp_c(thisruninds,respwin)',lfpresp_s(thisruninds,respwin)',params);
                        r1_tbtlfpcoherence(blck,ch,l,d,:) = nanmean(tbtlfpcoher(thisruninds,:),1);
                        r1_tbtlfpcoherr(blck,ch,l,d,:) = nanstd(tbtlfpcoher(thisruninds,:),1,1)./sqrt(length(thisruninds));
                        r1condxc(blck,ch,l,d,:) = squeeze(mean(xc(thisruninds,:),1));
%                         [r1_sfcoher_c(blck,ch,l,d,:),a,b,c,de,sfcfx,e,f,g,r1_sfcoherr_c(blck,ch,l,d,:,:)] = coherencycpt(lfpresp_c(thisruninds,1001:1800)',ptresp_c(thisruninds),params);
%                         [r1_sfcoher_s(blck,ch,l,d,:),a,b,c,de,sfcfx,e,f,g,r1_sfcoherr_s(blck,ch,l,d,:,:)] = coherencycpt(lfpresp_s(thisruninds,1001:1800)',ptresp_s(thisruninds),params);
                        r1_condcorstrength(blck,ch,l,d) = mean(corstrength(thisruninds));
                        r1_condmuaac_c(blck,ch,l,d,:) = squeeze(nanmean(muaac_c(thisruninds,:),1));
                        r1_condmuaac_s(blck,ch,l,d,:) = squeeze(nanmean(muaac_s(thisruninds,:),1));
                        r1_ntrials(blck,ch,l,d) = length(thisruninds);
                    else
                        r1_lfpspect_c(blck,ch,l,d,:) = nan(1,length(lfpspect_c(blck,ch,l,d,:)));
                        r1_lfpspect_s(blck,ch,l,d,:) = nan(1,length(lfpspect_s(blck,ch,l,d,:)));
                        r1_lfpspecterr_c(blck,ch,l,d,:,:) = nan(2,length(lfpspect_c(blck,ch,l,d,:)));
                        r1_lfpspecterr_s(blck,ch,l,d,:,:) = nan(2,length(lfpspect_s(blck,ch,l,d,:)));
                        r1_tbtlfpspect_c(blck,ch,l,d,:) = nan(1,length(lfpspect_c(blck,ch,l,d,:)));
                        r1_tbtlfpspect_s(blck,ch,l,d,:) = nan(1,length(lfpspect_c(blck,ch,l,d,:)));
                        r1_tbtlfpspecterr_c(blck,ch,l,d,:) = nan(1,length(lfpspect_c(blck,ch,l,d,:)));
                        r1_tbtlfpspecterr_s(blck,ch,l,d,:) = nan(1,length(lfpspect_c(blck,ch,l,d,:)));
                        r1_lfpcoherence(blck,ch,l,d,:) = nan(1,length(lfpcoherence(blck,ch,l,d,:)));
                        r1_lfpcoherr(blck,ch,l,d,:,:) = nan(2,length(lfpcoherence(blck,ch,l,d,:)));
                        r1_tbtlfpcoherence(blck,ch,l,d,:) = nan(1,length(lfpcoherence(blck,ch,l,d,:)));
                        r1_tbtlfpcoherr(blck,ch,l,d,:) = nan(1,length(lfpcoherence(blck,ch,l,d,:)));
                        r1condxc(blck,ch,l,d,:) = nan(1,size(xc,2));
                        r1_condmuaac_c(blck,ch,l,d,:) = nan(1,size(xc,2));
                        r1_condmuaac_s(blck,ch,l,d,:) = nan(1,size(xc,2));
%                         r1_sfcoher_c(blck,ch,l,d,:) = nan(1,length(sfcoher_c(blck,ch,l,d,:)));
%                         r1_sfcoher_s(blck,ch,l,d,:) = nan(1,length(sfcoher_s(blck,ch,l,d,:)));
%                         r1_sfcohererr_c(blck,ch,l,d,:,:) = nan(2,length(sfcoher_c(blck,ch,l,d,:)));
%                         r1_sfcohererr_s(blck,ch,l,d,:,:) = nan(2,length(sfcoher_s(blck,ch,l,d,:)));
                        r1_condcorstrength(blck,ch,l,d) = nan;
                        r1_ntrials(blck,ch,l,d) = 0;
                    end
                    
                    if ~isempty(thisstillinds)                    
                        [r0_lfpspect_c(blck,ch,l,d,:),chf,r0_lfpspecterr_c(blck,ch,l,d,:,:)] = mtspectrumc(squeeze(lfpresp_c(thisstillinds,1001:1800))',params);
                        [r0_lfpspect_s(blck,ch,l,d,:),chf,r0_lfpspecterr_s(blck,ch,l,d,:,:)] = mtspectrumc(squeeze(lfpresp_s(thisstillinds,1001:1800))',params);
                        r0_tbtlfpspect_c(blck,ch,l,d,:) = nanmean(lfpspecttbt_c(thisstillinds,:),1);
                        r0_tbtlfpspect_s(blck,ch,l,d,:) = nanmean(lfpspecttbt_s(thisstillinds,:),1);
                        r0_tbtlfpspecterr_c(blck,ch,l,d,:) = nanstd(lfpspecttbt_c(thisstillinds,:),1,1)./sqrt(length(thisstillinds));
                        r0_tbtlfpspecterr_s(blck,ch,l,d,:) = nanstd(lfpspecttbt_s(thisstillinds,:),1,1)./sqrt(length(thisstillinds));
                        [r0_lfpcoherence(blck,ch,l,d,:),a,b,c,de,cfx,e,f,r0_lfpcoherr(blck,ch,l,d,:,:)] = coherencyc(lfpresp_c(thisstillinds,respwin)',lfpresp_s(thisstillinds,respwin)',params);
                        r0_tbtlfpcoherence(blck,ch,l,d,:) = nanmean(tbtlfpcoher(thisstillinds,:),1);
                        r0_tbtlfpcoherr(blck,ch,l,d,:) = nanstd(tbtlfpcoher(thisstillinds,:),1,1)./sqrt(length(thisstillinds));
                        r0condxc(blck,ch,l,d,:) = squeeze(mean(xc(thisstillinds,:),1));
%                         [r0_sfcoher_c(blck,ch,l,d,:),a,b,c,de,sfcfx,e,f,g,r0_sfcoherr_c(blck,ch,l,d,:,:)] = coherencycpt(lfpresp_c(thisstillinds,1001:1800)',ptresp_c(thisstillinds),params);
%                         [r0_sfcoher_s(blck,ch,l,d,:),a,b,c,de,sfcfx,e,f,g,r0_sfcoherr_s(blck,ch,l,d,:,:)] = coherencycpt(lfpresp_s(thisstillinds,1001:1800)',ptresp_s(thisstillinds),params);
                        r0_condcorstrength(blck,ch,l,d,:) = mean(corstrength(thisstillinds));
                        r0_condmuaac_c(blck,ch,l,d,:) = squeeze(nanmean(muaac_c(thisstillinds,:),1));
                        r0_condmuaac_s(blck,ch,l,d,:) = squeeze(nanmean(muaac_s(thisstillinds,:),1));
                        r0_ntrials(blck,ch,l,d) = length(thisstillinds);
                    else
                        r0_lfpspect_c(blck,ch,l,d,:) = nan(1,length(lfpspect_c(blck,ch,l,d,:)));
                        r0_lfpspect_s(blck,ch,l,d,:) = nan(1,length(lfpspect_s(blck,ch,l,d,:)));
                        r0_lfpspecterr_c(blck,ch,l,d,:,:) = nan(2,length(lfpspect_c(blck,ch,l,d,:)));
                        r0_lfpspecterr_s(blck,ch,l,d,:,:) = nan(2,length(lfpspect_s(blck,ch,l,d,:)));
                        r0_lfpcoherence(blck,ch,l,d,:) = nan(1,length(lfpcoherence(blck,ch,l,d,:)));
                        r0_lfpcoherr(blck,ch,l,d,:,:) = nan(2,length(lfpcoherence(blck,ch,l,d,:)));
                        r0_tbtlfpcoherence(blck,ch,l,d,:) = nan(1,length(lfpcoherence(blck,ch,l,d,:)));
                        r0_tbtlfpcoherr(blck,ch,l,d,:) = nan(1,length(lfpcoherence(blck,ch,l,d,:)));
                        r0condxc(blck,ch,l,d,:) = nan(1,size(xc,2));
                        r0_condmuaac_c(blck,ch,l,d,:) = nan(1,size(xc,2));
                        r0_condmuaac_s(blck,ch,l,d,:) = nan(1,size(xc,2));
%                         r0_sfcoher_c(blck,ch,l,d,:) = nan(1,length(sfcoher_c(blck,ch,l,d,:)));
%                         r0_sfcoher_s(blck,ch,l,d,:) = nan(1,length(sfcoher_s(blck,ch,l,d,:)));
%                         r0_sfcohererr_c(blck,ch,l,d,:,:) = nan(2,length(sfcoher_c(blck,ch,l,d,:)));
%                         r0_sfcohererr_s(blck,ch,l,d,:,:) = nan(2,length(sfcoher_s(blck,ch,l,d,:)));
                        r0_condcorstrength(blck,ch,l,d) = nan;
                        r0_ntrials(blck,ch,l,d) = 0;
                    end
                end
            end
        end
    end
    save(popfile, '-v7.3');
else
    load(popfile);
end

centdepth = centdepth.*(cosd(penangle)*cosd(22));
surrdepth = surrdepth.*(cosd(penangle)*cosd(22));
spacing = spacing.*(cosd(penangle)*cosd(22));
evaldepth = 275;
for i = 1:length(centdepth)
    for j = 1:length(centelectrodes(i,1):centelectrodes(i,2))
        cdm(i,j) = centdepth(i)-((j-1)*spacing(i)); % center depth matrix
        sdm(i,j) = surrdepth(i)-((j-1)*spacing(i)); % surround depth matrix
    end
    [~,di(i)] = min(abs(cdm(i,:)-evaldepth)); % get depth index, least distance to evaldepth
    valid(i) = leffect(i,di(i));              % valid if light has a significant effect on large size stimuli
end
% valid(isnan(valid)) = 0;
valid = ones(size(di));
valid = logical(valid);
vi = find(valid);

for i = 1:length(vi)
%     lgpow(i,:,:) = squeeze(lfpspect_c(vi(i),di(vi(i)),:,:,lgi(vi(i),di(vi(i))))); % at second contact only
    r1lgpow(i,:,:) = squeeze(r1_lfpspect_c(vi(i),di(vi(i)),:,:,lgi(vi(i),di(vi(i)))));
%     lglfpcoher(i,:,:) = squeeze(lfpcoherence(vi(i),di(vi(i)),:,:,lgi(vi(i),di(vi(i)))));
%     r1lglfpcoher(i,:,:) = squeeze(r1_lfpcoherence(vi(i),di(vi(i)),:,:,lgi(vi(i),di(vi(i)))));
    r1tbtlglfpcoher(i,:,:) = squeeze(r1_tbtlfpcoherence(vi(i),di(vi(i)),:,:,lgi(vi(i),di(vi(i)))));
%     lgsfcoher(i,:,:) = squeeze(sfcoher_c(vi(i),di(vi(i)),:,:,lgi(vi(i),di(vi(i)))));
%     r1lgsfcoher(i,:,:) = squeeze(r1_sfcoher_c(vi(i),di(vi(i)),:,:,lgi(vi(i),di(vi(i)))));
    
    for l = 1:2
        for d = 1:4
            r1fillspecy(i,l,d,:) = [squeeze(r1_lfpspecterr_c(vi(i),di(vi(i)),l,d,1,1:104))',fliplr(squeeze(r1_lfpspecterr_c(vi(i),2,l,d,2,1:104))')];            
            r1filltbtspecy(i,l,d,:) = [(squeeze(r1_tbtlfpspect_c(i,di(i),l,d,1:104))+squeeze(r1_tbtlfpspecterr_c(i,di(i),l,d,1:104)))',fliplr((squeeze(r1_tbtlfpspect_c(i,di(i),l,d,1:104))-squeeze(r1_tbtlfpspecterr_c(i,di(i),l,d,1:104)))')]';
%             r1fillcohery(i,l,d,:) = [squeeze(r1_lfpcoherr(vi(i),di(vi(i)),l,d,1,1:104))',fliplr(squeeze(r1_lfpcoherr(vi(i),2,l,d,2,1:104))')]';
            r1filltbtcohery(i,l,d,:) = [(squeeze(r1_tbtlfpcoherence(vi(i),di(vi(i)),l,d,1:104))+squeeze(r1_tbtlfpcoherr(vi(i),di(vi(i)),l,d,1:104)))',fliplr((squeeze(r1_tbtlfpcoherence(vi(i),di(vi(i)),l,d,1:104))-squeeze(r1_tbtlfpcoherr(vi(i),di(vi(i)),l,d,1:104)))')]';
        end
    end    
end
fillx = [chf(1:104),fliplr(chf(1:104))];

% spectrum iso vs cross
for i =1:size(r1filltbtspecy,1)
    figure
    fill(fillx,squeeze(r1filltbtspecy(i,1,1,:)),[.3,.3,1])
    hold on
    fill(fillx,squeeze(r1filltbtspecy(i,1,2,:)),[1,.3,.3])
    set(gca,'yscale','log')
    legend('iso','cross')
    title(animalids{i})
end

% power iso vs cross scatter
[p,s] = signrank(r1lgpow(:,1,1),r1lgpow(:,1,2));
figure
plot(r1lgpow(:,1,1),r1lgpow(:,1,2),'ko','markerfacecolor','k')
line([0,350],[0,350],'color','k')
axis square
xlabel('gamma power same orientation')
ylabel('gamma power 90 off orientation')
title(['signrank p = ' num2str(p)]);
n = size(r1lgpow,1);
disp(['same ori power: ' num2str(nanmean(r1lgpow(:,1,1))) '+-' num2str(std(r1lgpow(:,1,1))./sqrt(n))]);
disp(['diff ori power: ' num2str(nanmean(r1lgpow(:,1,2))) '+-' num2str(std(r1lgpow(:,1,2))./sqrt(n))]);
pcc = (1-r1lgpow(:,1,2)./r1lgpow(:,1,1)).*100;
disp(['percent change cross/iso = ' num2str(mean(pcc)) '+-' num2str(std(pcc)./sqrt(n))])

% power iso vs cross connected columns
figure
plot(1,r1lgpow(:,1,1),'ko','markerfacecolor','k')
hold on
plot(2,r1lgpow(:,1,2),'o','color',[.5,.5,.5],'markerfacecolor',[.5,.5,.5])
axis([.5,2.5,0,300])
for i = 1:size(r1lgpow,1)    
    plot([1,2],[r1lgpow(i,1,1),r1lgpow(i,1,2)],'k')
end
set(gca,'xtick',[1,2]);
set(gca,'xticklabel',{'iso','cross'})
title('change in gamma power')
ylabel('gamma power')
set(gcf,'OuterPosition',[186   300   240   513])


% spectrum L0 vs L1
for i =1:size(r1filltbtspecy,1)
    figure
    fill(fillx,squeeze(r1filltbtspecy(i,1,1,:)),[.3,.3,1])
    hold on
    fill(fillx,squeeze(r1filltbtspecy(i,2,1,:)),[1,.3,.3])
    set(gca,'yscale','log')
    legend('L0','L1')
    title(animalids{i})
end

% power L0 vs L1 scatter
[p,s] = signrank(r1lgpow(:,1,1),r1lgpow(:,2,1));
figure
plot(r1lgpow(:,1,1),r1lgpow(:,2,1),'ko','markerfacecolor','k')
line([0,350],[0,350],'color','k')
axis square
xlabel('gamma power L0')
ylabel('gamma power L1')
title(['signrank p = ' num2str(p)]);
n = size(r1lgpow,1);
disp(['same ori power: ' num2str(nanmean(r1lgpow(:,1,1))) '+-' num2str(std(r1lgpow(:,1,1))./sqrt(n))]);
disp(['diff ori power: ' num2str(nanmean(r1lgpow(:,2,1))) '+-' num2str(std(r1lgpow(:,2,1))./sqrt(n))]);
pcc = (1-r1lgpow(:,2,1)./r1lgpow(:,1,1)).*100;
disp(['percent change L2/L1 = ' num2str(mean(pcc)) '+-' num2str(std(pcc)./sqrt(n))])

% power L0 vs L1 connected columns
figure
plot(1,r1lgpow(:,1,1),'bo','markerfacecolor','b')
hold on
plot(2,r1lgpow(:,2,1),'ro','markerfacecolor','r')
axis([.5,2.5,0,300])
for i = 1:size(r1lgpow,1)    
    plot([1,2],[r1lgpow(i,1,1),r1lgpow(i,2,1)],'k')
end
set(gca,'xtick',[1,2]);
set(gca,'xticklabel',{'L0','L1'})
title('change in gamma power')
ylabel('gamma power')
set(gcf,'OuterPosition',[186   300   240   513])


% coherence iso vs cross, trial by trial standard error
for i =1:size(r1filltbtcohery,1)
    figure
    fill(fillx,squeeze(r1filltbtcohery(i,1,1,:)),[.3,.3,1])
    hold on
    fill(fillx,squeeze(r1filltbtcohery(i,1,2,:)),[1,.3,.3])
    legend('iso','cross')
end

% coherence iso vs cross scatter
[p,s] = signrank(r1tbtlglfpcoher(:,1,1),r1tbtlglfpcoher(:,1,2));
n = size(r1tbtlglfpcoher,1); 
figure
plot(r1tbtlglfpcoher(:,1,1),r1tbtlglfpcoher(:,1,2),'ko','markerfacecolor','k')
line([0,1],[0,1],'color','k')
axis square
xlabel('gamma coherence same orientation')
ylabel('gamma coherence 90 off orientation')
title(['signrank p = ' num2str(p)]);
disp(['same ori coherency: ' num2str(nanmean(r1tbtlglfpcoher(:,1,1))) '+-' num2str(std(r1tbtlglfpcoher(:,1,1))./sqrt(n))]);
disp(['diff ori coherency: ' num2str(nanmean(r1tbtlglfpcoher(:,1,2))) '+-' num2str(std(r1tbtlglfpcoher(:,1,2))./sqrt(n))]);

%coherence iso vs. cross connected columns
figure
plot(1,r1tbtlglfpcoher(:,1,1),'ko','markerfacecolor','k')
hold on
plot(2,r1tbtlglfpcoher(:,1,2),'o','color',[.5,.5,.5],'markerfacecolor',[.5,.5,.5])
axis([.5,2.5,.25,.9])
for i = 1:size(r1tbtlglfpcoher,1)    
    plot([1,2],[r1tbtlglfpcoher(i,1,1),r1tbtlglfpcoher(i,1,2)],'k')
end
set(gca,'xtick',[1,2]);
set(gca,'xticklabel',{'iso','cross'})
title('coherence between sites')
ylabel('coherence')
set(gcf,'OuterPosition',[186   438   240   513])


% coherence light off vs on, trial by trial standard error
for i =1:size(r1filltbtcohery,1)
    figure
    fill(fillx,squeeze(r1filltbtcohery(i,1,1,:)),[.3,.3,1])
    hold on
    fill(fillx,squeeze(r1filltbtcohery(i,2,1,:)),[1,.3,.3])
    legend('L0','L1')
end

% coherence L0 vs L1 scatter
[p,s] = signrank(r1tbtlglfpcoher(:,1,1),r1tbtlglfpcoher(:,2,1));
n = size(r1tbtlglfpcoher,1); 
figure
plot(r1tbtlglfpcoher(:,1,1),r1tbtlglfpcoher(:,2,1),'ko','markerfacecolor','k')
line([0,1],[0,1],'color','k')
axis square
xlabel('gamma power same orientation')
ylabel('gamma power 90 off orientation')
title(['signrank p = ' num2str(p)]);
disp(['same ori L0 coherency: ' num2str(nanmean(r1tbtlglfpcoher(:,1,1))) '+-' num2str(std(r1tbtlglfpcoher(:,1,1))./sqrt(n))]);
disp(['same ori L1 coherency: ' num2str(nanmean(r1tbtlglfpcoher(:,2,1))) '+-' num2str(std(r1tbtlglfpcoher(:,2,1))./sqrt(n))]);

% coherence L0 vs L1 connected columns
figure
plot(1,r1tbtlglfpcoher(:,1,1),'ko','markerfacecolor','k')
hold on
plot(2,r1tbtlglfpcoher(:,2,1),'ro','markerfacecolor','r')
axis([.5,2.5,.25,.9])
for i = 1:size(r1tbtlglfpcoher,1)    
    plot([1,2],[r1tbtlglfpcoher(i,1,1),r1tbtlglfpcoher(i,2,1)],'k')
end
set(gca,'xtick',[1,2]);
set(gca,'xticklabel',{'L0','L1'})
title('coherence between sites')
ylabel('coherence')
set(gcf,'OuterPosition',[186   438   240   513])


n = size(r1tbtlglfpcoher,1);
anovavec = [squeeze(r1tbtlglfpcoher(:,1,1));squeeze(r1tbtlglfpcoher(:,1,2));squeeze(r1tbtlglfpcoher(:,2,1));squeeze(r1tbtlglfpcoher(:,2,2))];
gc = [zeros(n,1);ones(n,1);zeros(n,1);ones(n,1)]; %size
gl = [zeros(n*2,1);ones(n*2,1)];
[p,table,stats] = anovan(anovavec,{gc,gl},'model','full');
multcompare(stats)
multcompare(stats,'dimension',2)
[s,p] = ttest(r1tbtlglfpcoher(:,1,1),r1tbtlglfpcoher(:,1,2))
[s,p] = ttest(r1tbtlglfpcoher(:,1,1),r1tbtlglfpcoher(:,2,1))


% same for running only
for i = 1:length(blocks)
    figure('name',animalids{i});
    subplot(2,2,1)
    semilogy(chf,squeeze(r1_lfpspect_c(i,3,1,1,:)),'g')
    hold on    
    semilogy(chf,squeeze(r1_lfpspect_c(i,3,1,2,:)),'r')
    semilogy(chf,squeeze(r1_lfpspect_c(i,3,1,3,:)),'b')
    semilogy(chf,squeeze(r1_lfpspect_c(i,3,1,4,:)),'c')
    ax = axis;
    axis([0,120,1,ax(4)])
    legend('same','90deg','on field only','off field only')
    xlabel('frequency')
    ylabel('psd')
    
    subplot(2,2,2)
    plot(chf,squeeze(r1_lfpcoherence(i,3,1,1,:)),'g');
    hold on
    plot(chf,squeeze(r1_lfpcoherence(i,3,1,2,:)),'r');
    plot(chf,squeeze(r1_lfpcoherence(i,3,1,3,:)),'b');
    plot(chf,squeeze(r1_lfpcoherence(i,3,1,4,:)),'c');
    axis([0,120,0,1]);
    xlabel('frequency')
    ylabel('site coherence')    
    
%     subplot(2,3,3)
%     plot(chf,squeeze(r1_sfcoher_c(i,3,1,1,:)),'g');
%     hold on
%     plot(chf,squeeze(r1_sfcoher_c(i,3,1,2,:)),'r');
%     plot(chf,squeeze(r1_sfcoher_c(i,3,1,3,:)),'b');
%     plot(chf,squeeze(r1_sfcoher_c(i,3,1,4,:)),'c');
%     axis([0,120,0,1]);
%     xlabel('frequency')
%     ylabel('spike field coherence')
    
    subplot(2,2,3)
    semilogy(chf,squeeze(r1_lfpspect_c(i,3,1,3,:)),'b')
    hold on    
    semilogy(chf,squeeze(r1_lfpspect_c(i,3,1,1,:)),'c')
    semilogy(chf,squeeze(r1_lfpspect_c(i,3,2,3,:)),'r')
    semilogy(chf,squeeze(r1_lfpspect_c(i,3,2,1,:)),'m')
    ax = axis;
    axis([0,120,1,ax(4)])
    legend('small l0','same l0','small l1','same l1')
    xlabel('frequency')
    ylabel('psd')
    
    subplot(2,2,4)
    plot(chf,squeeze(r1_lfpcoherence(i,3,1,3,:)),'b');
    hold on
    plot(chf,squeeze(r1_lfpcoherence(i,3,1,1,:)),'c');
    plot(chf,squeeze(r1_lfpcoherence(i,3,2,3,:)),'r');
    plot(chf,squeeze(r1_lfpcoherence(i,3,2,1,:)),'m');
    axis([0,120,0,1]);
    xlabel('frequency')
    ylabel('site coherence')
    
%     subplot(2,3,6)
%     plot(chf,squeeze(r1_sfcoher_c(i,3,1,3,:)),'b');
%     hold on
%     plot(chf,squeeze(r1_sfcoher_c(i,3,1,1,:)),'c');
%     plot(chf,squeeze(r1_sfcoher_c(i,3,2,3,:)),'r');
%     plot(chf,squeeze(r1_sfcoher_c(i,3,2,1,:)),'m');
%     axis([0,120,0,1]);
%     xlabel('frequency')
%     ylabel('spike field coherence')
    
    set(gcf,'OuterPosition',[50,100,1200,700])
end
    
for i = 1:length(blocks)
    
    depthax = centdepth(i):-spacing(i):centdepth(i)-(length(centelectrodes(i,1):centelectrodes(i,2))-1)*spacing(i);
    
    figure
    subplot(3,3,1)
    imagesc(chf,depthax,log(squeeze(r1_lfpspect_c(i,:,1,3,:))))
    axis([1,120,depthax(end)-10,depthax(1)+10])
    caxis([-2,9])
    title('power center')
    
    subplot(3,3,2)
    imagesc(chf,depthax,log(squeeze(r1_lfpspect_c(i,:,1,1,:))))
    axis([1,120,depthax(end)-10,depthax(1)+10])
    caxis([-2,9])
    title('power iso')
    
    subplot(3,3,3)
    imagesc(chf,depthax,log(squeeze(r1_lfpspect_c(i,:,1,2,:))))
    axis([1,120,depthax(end)-10,depthax(1)+10])
    caxis([-2,9])
    title('power cross')
    
    subplot(3,3,4)
    imagesc(chf,depthax,log(squeeze(r1_lfpspect_c(i,:,1,1,:)))-log(squeeze(r1_lfpspect_c(i,:,1,3,:))));
    axis([1,120,depthax(end)-10,depthax(1)+10])
    caxis([-2,2])
    title('large - small')
    
    subplot(3,3,5)
    imagesc(chf,depthax,log(squeeze(r1_lfpspect_c(i,:,1,2,:)))-log(squeeze(r1_lfpspect_c(i,:,1,1,:))));
    axis([1,120,depthax(end)-10,depthax(1)+10])
    caxis([-2,2])
    title('cross - iso')    
    
    subplot(3,3,7)
    imagesc(chf,depthax,log(squeeze(r1_lfpspect_c(i,:,2,3,:)))-log(squeeze(r1_lfpspect_c(i,:,1,3,:))));
    axis([1,120,depthax(end)-10,depthax(1)+10])
    caxis([-2,2])
    title('small L1-L0')
    
    subplot(3,3,8)
    imagesc(chf,depthax,log(squeeze(r1_lfpspect_c(i,:,2,1,:)))-log(squeeze(r1_lfpspect_c(i,:,1,1,:))));
    axis([1,120,depthax(end)-10,depthax(1)+10])
    caxis([-2,2])
    title('iso L1-L0')
    
    subplot(3,3,9)
    imagesc(chf,depthax,log(squeeze(r1_lfpspect_c(i,:,2,2,:)))-log(squeeze(r1_lfpspect_c(i,:,1,2,:))));
    axis([1,120,depthax(end)-10,depthax(1)+10])
    caxis([-2,2])
    title('cross L1-L0')
    
    
    figure
    subplot(3,3,1)
    imagesc(chf,depthax,squeeze(r1_tbtlfpcoherence(i,:,1,3,:)))
    axis([1,120,depthax(end)-10,depthax(1)+10])
    caxis([0,1])
    title('coherence RF only')
    
    subplot(3,3,2)
    imagesc(chf,depthax,squeeze(r1_tbtlfpcoherence(i,:,1,1,:)))
    axis([1,120,depthax(end)-10,depthax(1)+10])
    caxis([0,1])
    title('match surround')    
    
    subplot(3,3,3)
    imagesc(chf,depthax,squeeze(r1_tbtlfpcoherence(i,:,1,2,:)))
    axis([1,120,depthax(end)-10,depthax(1)+10])
    caxis([0,1])
    title('rotated surround')    
    
    subplot(3,3,4)
    imagesc(chf,depthax,squeeze(r1_tbtlfpcoherence(i,:,1,1,:))-squeeze(r1_tbtlfpcoherence(i,:,1,3,:)))
    axis([1,120,depthax(end)-10,depthax(1)+10])
    caxis([-.5,.5])
    title('iso large - small')    
    
    subplot(3,3,5)
    imagesc(chf,depthax,squeeze(r1_tbtlfpcoherence(i,:,1,2,:))-squeeze(r1_tbtlfpcoherence(i,:,1,1,:)))
    axis([1,120,depthax(end)-10,depthax(1)+10])
    caxis([-.3,.3])
    title('cross surround - iso surround')
    
    subplot(3,3,7)
    imagesc(chf,depthax,squeeze(r1_tbtlfpcoherence(i,:,2,3,:))-squeeze(r1_tbtlfpcoherence(i,:,1,3,:)))
    axis([1,120,depthax(end)-10,depthax(1)+10])
    caxis([-.3,.3])
    title('small L1-L0')       
    
    subplot(3,3,8)
    imagesc(chf,depthax,squeeze(r1_tbtlfpcoherence(i,:,2,1,:))-squeeze(r1_tbtlfpcoherence(i,:,1,1,:)))
    axis([1,120,depthax(end)-10,depthax(1)+10])
    caxis([-.3,.3])
    title('iso L1-L0')    
    
    subplot(3,3,9)
    imagesc(chf,depthax,squeeze(r1_tbtlfpcoherence(i,:,2,2,:))-squeeze(r1_tbtlfpcoherence(i,:,1,2,:)))
    axis([1,120,depthax(end)-10,depthax(1)+10])
    caxis([-.3,.3])
    title('cross L1-L0')    
    
%     figure
%     subplot(3,3,1)
%     imagesc(chf,depthax,squeeze(sfcoher_c(i,:,1,3,:)))
%     axis([1,120,depthax(end)-10,depthax(1)+10])
%     caxis([0,.6])
%     title('SFC small')
%     
%     subplot(3,3,2)
%     imagesc(chf,depthax,squeeze(sfcoher_c(i,:,1,1,:)))
%     axis([1,120,depthax(end)-10,depthax(1)+10])
%     caxis([0,.6])
%     title('SFC iso')
%     
%     subplot(3,3,3)
%     imagesc(chf,depthax,squeeze(sfcoher_c(i,:,1,2,:)))
%     axis([1,120,depthax(end)-10,depthax(1)+10])
%     caxis([0,.6])
%     title('SFC cross')    
%     
%     subplot(3,3,4)
%     imagesc(chf,depthax,squeeze(sfcoher_c(i,:,1,1,:))-squeeze(sfcoher_c(i,:,1,3,:)))
%     axis([1,120,depthax(end)-10,depthax(1)+10])
%     caxis([-.2,.2])
%     title('SFC large - small')
%     
%     subplot(3,3,5)
%     imagesc(chf,depthax,squeeze(sfcoher_c(i,:,1,2,:))-squeeze(sfcoher_c(i,:,1,1,:)))
%     axis([1,120,depthax(end)-10,depthax(1)+10])
%     caxis([-.2,.2])
%     title('SFC cross - iso')
%     
%     subplot(3,3,7)    
%     imagesc(chf,depthax,squeeze(sfcoher_c(i,:,2,3,:))-squeeze(sfcoher_c(i,:,1,3,:)))
%     axis([1,120,depthax(end)-10,depthax(1)+10])
%     caxis([-.2,.2])
%     title('small L1-L0')
%     
%     subplot(3,3,8)    
%     imagesc(chf,depthax,squeeze(sfcoher_c(i,:,2,1,:))-squeeze(sfcoher_c(i,:,1,1,:)))
%     axis([1,120,depthax(end)-10,depthax(1)+10])
%     caxis([-.2,.2])
%     title('iso L1-L0')
%     
%     subplot(3,3,9)    
%     imagesc(chf,depthax,squeeze(sfcoher_c(i,:,2,2,:))-squeeze(sfcoher_c(i,:,1,2,:)))
%     axis([1,120,depthax(end)-10,depthax(1)+10])
%     caxis([-.2,.2])
%     title('cross L1-L0')

end


 disp('');