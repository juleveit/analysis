function MUA_population_phassurroundcohere

% % 
% % PV Halo population
% animalids =    {'150818', '150820', '150824','151104'};
% blocks    =     [5,        11,       7,       8];
% animal    =     [1,         2,       3,       4];
% centelectrodes =[[1,16];   [1,16];   [17,32]; [17,32]];
% centdepth    =  [400,      400,      400,      400];
% surrelectrodes =[[17,32];  [17,32];  [1,16];  [1,16]];
% surrdepth    =  [400,      400,      400,      400];
% penangle =      [25,       25,       25,       25];
% spacing =       [25,       25,       25,       25];
% popfile = 'C:\Users\Julia\work\data\populations\PV_Halo\phassurcoher\MUA_population.mat';
% % 
% % SOM Halo population
% animalids =    {'150825', '150902', '150909', '151022', '151023', '151027', '151109', '151110'};
% blocks    =     [9,        11,       7,        5,        6,        8,        7,        10];
% animal    =     [1,        2,        3,        4,        5,        6,        7,        8];
% centelectrodes =[[17,32];  [1,16];   [1,16];   [1,16];   [17,32];  [17,32];  [17,32];  [1,16]];
% centdepth    =  [400,      400,      400,      400,      400,      500,      500,      500];
% surrelectrodes =[[1,16];   [17,32];  [17,32];  [17,32];  [1,16];   [1,16];   [1,16];   [17,32]];
% surrdepth    =  [400,      400,      400,      400,      400,      500,      500,      500];
% penangle =      [25,       25,       25,       25,       25,       25,       25,       25];
% spacing =       [25,       25,       25,       25,       25,       25,       25,       25];
% popfile = 'C:\Users\Julia\work\data\populations\SOM_Halo_later\phassurcoher\MUA_population.mat';

% SOM PV Mix population
animalids =    {'150825', '150902', '150909', '151022', '151023', '151027', '151109', '151110', '150818', '150820', '150824','151104'};
blocks    =     [9,        11,       7,        5,        6,        8,        7,        10,       5,        11,       7,       8];
animal    =     [1,        2,        3,        4,        5,        6,        7,        8,        9,        10,       11,      12];
centelectrodes =[[17,32];  [1,16];   [1,16];   [1,16];   [17,32];  [17,32];  [17,32];  [1,16];  [1,16];   [1,16];   [17,32]; [17,32]];
centdepth    =  [400,      400,      400,      400,      400,      500,      500,      500,      400,      400,      400,      400];
surrelectrodes =[[1,16];   [17,32];  [17,32];  [17,32];  [1,16];   [1,16];   [1,16];   [17,32]; [17,32];  [17,32];  [1,16];  [1,16]];
surrdepth    =  [400,      400,      400,      400,      400,      500,      500,      500,      400,      400,      400,      400];
penangle =      [25,       25,       25,       25,       25,       25,       25,       25,       25,       25,       25,       25];
spacing =       [25,       25,       25,       25,       25,       25,       25,       25,       25,       25,       25,       25];
popfile = 'C:\Users\Julia\work\data\populations\SOMPVcombined\phassurcoher\MUA_population.mat';

% % single thing
% animalids =    {'150915'};
% blocks    =     [11];
% animal    =     [1];
% centelectrodes =[[17,32]];
% centdepth    =  [400];
% surrelectrodes =[[1,16]];
% surrdepth    =  [400];
% penangle =      [25];
% spacing =       [25];
% popfile = 'C:\Users\Julia\work\data\populations\temp1.mat';

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
        respwin = 501:1500; % after stimulus onset
        offsetwin = 1501:2500;
        respwin = respwin+prestim;
        
        clear bipderlfp_c; clear bipderlfp_s;
        for ch = 1: length(centelectrodes(blck,1):centelectrodes(blck,2))
            
            disp(['Block ' int2str(blck) '/' int2str(length(blocks)) '   channel ' int2str(ch) '/' int2str(length(centelectrodes(blck,1):centelectrodes(blck,2)))]);
            
            depthc(ch) = centdepth(blck)-(ch-1)*spacing(blck);
            depths(ch) = surrdepth(blck)-(ch-1)*spacing(blck);
            
            trialdur = centresult.stimduration*1000;
            msstamps = centresult.msstamps;
            if length(msstamps)~=length(centresult.light)
% %                 msstamps([193,291]) = []; % for 151023 block 5
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
            okspeed = find(meanspeed>( mean(meanspeed(notstill))-(1.5*std(meanspeed(notstill))) ) & meanspeed>1);
            okvar = find(stdspeed<( mean(stdspeed(notstill))+(1.5*std(stdspeed(notstill)))) & stdspeed>.5);
            oktrials = intersect(okspeed,okvar);
            nonoktrials = 1:size(speed,1); nonoktrials(oktrials) = [];
            stilltrials = 1:size(speed,1); stilltrials(notstill) = [];
            
            % bipolar derivation like fries
            if ~(ch == length(centelectrodes(blck,1):centelectrodes(blck,2)))
                bipderlfp_c(ch,:) = centresult.lfp(ch,:)-centresult.lfp(ch+1,:);
                bipderlfp_s(ch,:) = surrresult.lfp(ch,:)-surrresult.lfp(ch+1,:);
            else
                bipderlfp_c(ch,:) = centresult.lfp(ch,:);
                bipderlfp_s(ch,:) = surrresult.lfp(ch,:);
            end            
            
            % find low and high gamma peaks in spectra
            beta = [15,40];
            gamma = [50,70];
            large = find(centresult.gratingInfo.Phase_surround == 0 & centresult.gratingInfo.Orientation ~= -1 & centresult.light == 0);
            small = find(centresult.gratingInfo.Phase_surround == -1 & centresult.gratingInfo.Orientation ~= -1 & centresult.light == 0);
            lr1 = intersect(large,oktrials);
            sr1 = intersect(small,oktrials);
            for i = 1:length(large)
                [pl(i,:),f] = pmtm(centresult.lfp(ch,centresult.msstamps(large(i))+300:centresult.msstamps(large(i))+1300),3,[],1000);
                [ps(i,:),f] = pmtm(centresult.lfp(ch,centresult.msstamps(small(i))+300:centresult.msstamps(small(i))+1300),3,[],1000);
                [plbd(i,:),f] = pmtm(bipderlfp_c(ch,centresult.msstamps(large(i))+300:centresult.msstamps(large(i))+1300),3,[],1000);
                [psbd(i,:),f] = pmtm(bipderlfp_c(ch,centresult.msstamps(small(i))+300:centresult.msstamps(small(i))+1300),3,[],1000);
            end 
            plr1 = nan(1,size(pl,2)); psr1 = nan(1,size(pl,2));
            if length(lr1>=5)
                for i = 1:length(lr1)
                    [plr1(i,:),f] = mtspectrumc(centresult.lfp(ch,centresult.msstamps(lr1(i))+700:centresult.msstamps(lr1(i))+1500),params);
                    [plbdr1(i,:),f] = mtspectrumc(bipderlfp_c(ch,centresult.msstamps(lr1(i))+700:centresult.msstamps(lr1(i))+1500),params);
%                     [plr1(i,:),f] = pmtm(centresult.lfp(ch,centresult.msstamps(lr1(i))+300:centresult.msstamps(lr1(i))+1300),3,[],1000);
                end
            end
            if length(sr1>=5)
                for i = 1:length(sr1)
                    [psr1(i,:),f] = mtspectrumc(centresult.lfp(ch,centresult.msstamps(sr1(i))+700:centresult.msstamps(sr1(i))+1500),params);
                    [psbdr1(i,:),f] = mtspectrumc(bipderlfp_c(ch,centresult.msstamps(sr1(i))+700:centresult.msstamps(sr1(i))+1500),params);
%                     [psr1(i,:),f] = pmtm(centresult.lfp(ch,centresult.msstamps(sr1(i))+300:centresult.msstamps(sr1(i))+1300),3,[],1000);
                end
            end                       
            b1 = find(f>beta(1),1); b2 = find(f>beta(2),1);
            g1 = find(f>gamma(1),1); g2 = find(f>gamma(2),1);
            if ~isnan(plr1(1)) % take the peak from the running spectra if there is more than a couple of tunning trials
                bsig = nanmean(plr1(:,b1:b2));
                bdbsig = nanmean(plbdr1(:,b1:b2));
            else
                bsig = nanmean(pl(:,b1:b2));
                bdbsig = nanmean(plbd(:,b1:b2));
            end
            if ~isnan(psr1(1))
                gsig = nanmean(psr1(:,g1:g2));
                bdgsig = nanmean(psbdr1(:,g1:g2));
            else
                gsig = nanmean(ps(:,g1:g2));
                bdgsig = nanmean(psbd(:,g1:g2));
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
            if isempty(find(diff(gsig)>0)) % there is no clear gamma peak
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
             
            if isempty(find(diff(bdbsig)>0)) % there is no clear beta peak
                bdbpi = round((b1+b2)/2);
                g1bdpeak(blck,ch) = 0;
            else
                peaks = find(diff(bdbsig)>0)+1;
                pvs = bdbsig(peaks);
                bdbpi = peaks(pvs == max(pvs));
                bdbpi = bdbpi+b1-1;
                g1bdpeak(blck,ch) = 1;
            end
            if isempty(find(diff(bdgsig)>0)) % there is no clear beta peak
                bdgpi = round((g1+g2)/2);
                g2bdpeak(blck,ch) = 0;
            else
                peaks = find(diff(bdgsig)>0)+1;
                pvs = bdgsig(peaks);
                bdgpi = peaks(pvs == max(pvs));
                bdgpi = bdgpi+g1-1;
                g2bdpeak(blck,ch) = 1;
            end
             bdlgi(blck,ch) = bdbpi; bdhgi(blck,ch) = bdgpi;


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
                
                lfprespbd_c(i,:) = bipderlfp_c(ch,msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim);
                lfprespbd_s(i,:) = bipderlfp_s(ch,msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim);
                lfpspecttbtbd_c(i,:) = mtspectrumc(squeeze(lfprespbd_c(i,1001:1800))',params);
                lfpspecttbtbd_s(i,:) = mtspectrumc(squeeze(lfprespbd_s(i,1001:1800))',params);
                tbtlfpcoherbd(i,:) = coherencyc(lfprespbd_c(i,respwin)',lfprespbd_s(i,respwin)',params);
                
                muacresp_c(i,:) =  centresult.muac(ch,msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim).^2;
                muacresp_s(i,:) =  surrresult.muac(ch,msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim).^2;
                [mcxc(i,:),lags] = xcorr(muacresp_c(i,respwin),muacresp_s(i,respwin),200,'coeff');   
                muaac_c(i,:) = xcorr(resp_c(i,respwin),resp_c(i,respwin),200,'coeff');
                muaac_s(i,:) = xcorr(resp_s(i,respwin),resp_s(i,respwin),200,'coeff'); 
                
            end

            frs_c = sum(resp_c(:,respwin),2)./(length(respwin)/1000);
            frs_s = sum(resp_s(:,respwin),2)./(length(respwin)/1000);
            bl_c = sum(resp_c(:,1:prestim),2)./(prestim/1000);
            bl_s = sum(resp_s(:,1:prestim),2)./(prestim/1000);
            muacmean_c = mean(muacresp_c(:,respwin),2);
            muacmean_s = mean(muacresp_s(:,respwin),2);

            surroundphases = unique(centresult.gratingInfo.Phase_surround);
            for l = 1:length(unique(centresult.light))
                
                for d = 1:length(surroundphases)+1
                    if d == 4
                        % surround only
                        thisinds = find(centresult.gratingInfo.Orientation == -1 & centresult.light == l-1);
                    else
                        thisinds = find(centresult.gratingInfo.Orientation ~= -1 & ...
                            centresult.gratingInfo.Phase_surround == surroundphases(d) & ...
                            centresult.light == l-1);
                    end
                    
                    thisruninds = intersect(thisinds,oktrials); thisstillinds = intersect(thisinds,stilltrials);
                                    
                    condresp_c(blck,ch,l,d,:) = nanmean(resp_c(thisinds,:),1);
                    condresp_s(blck,ch,l,d,:) = nanmean(resp_s(thisinds,:),1);
                    condfr_c(blck,ch,l,d) = nanmean(frs_c(thisinds));
                    condfr_s(blck,ch,l,d) = nanmean(frs_s(thisinds));
                    condmuacresp_c(blck,ch,l,d,:) = nanmean(muacresp_c(thisinds,:),1);
                    condmuacresp_s(blck,ch,l,d,:) = nanmean(muacresp_s(thisinds,:),1);
                    condmua_c(blck,ch,l,d) = nanmean(muacmean_c(thisinds),1);
                    condmua_s(blck,ch,l,d) = nanmean(muacmean_s(thisinds),1);
                    condlfpresp_c(blck,ch,l,d,:) = nanmean(lfpresp_c(thisinds,:),1);
                    condlfpresp_s(blck,ch,l,d,:) = nanmean(lfpresp_s(thisinds,:),1);
                    condlfprespbd_c(blck,ch,l,d,:) = nanmean(lfprespbd_c(thisinds,:),1);
                    condlfprespbd_s(blck,ch,l,d,:) = nanmean(lfprespbd_s(thisinds,:),1);
                    
                    [lfpspect_c(blck,ch,l,d,:),chf,lfpspecterr_c(blck,ch,l,d,:,:)] = mtspectrumc(squeeze(lfpresp_c(thisinds,1001:1800))',params);
                    [lfpspect_s(blck,ch,l,d,:),chf,lfpspecterr_s(blck,ch,l,d,:,:)] = mtspectrumc(squeeze(lfpresp_s(thisinds,1001:1800))',params);
                    tbtlfpspect_c(blck,ch,l,d,:) = nanmean(lfpspecttbt_c(thisinds,:),1);
                    tbtlfpspect_s(blck,ch,l,d,:) = nanmean(lfpspecttbt_s(thisinds,:),1);
                    tbtlfpspecterr_c(blck,ch,l,d,:) = nanstd(lfpspecttbt_c(thisinds,:),1,1)./sqrt(length(thisinds));
                    tbtlfpspecterr_s(blck,ch,l,d,:) = nanstd(lfpspecttbt_s(thisinds,:),1,1)./sqrt(length(thisinds));
                    [muacspect_c(blck,ch,l,d,:),chf,muacspecterr_c(blck,ch,l,d,:,:)] = mtspectrumc(squeeze(muacresp_c(thisinds,1001:1800))',params);
                    [muacspect_s(blck,ch,l,d,:),chf,muacspecterr_s(blck,ch,l,d,:,:)] = mtspectrumc(squeeze(muacresp_s(thisinds,1001:1800))',params);
                                        
                    [lfpspectbd_c(blck,ch,l,d,:),chf,lfpspecterrbd_c(blck,ch,l,d,:,:)] = mtspectrumc(squeeze(lfprespbd_c(thisinds,1001:1800))',params);
                    [lfpspectbd_s(blck,ch,l,d,:),chf,lfpspecterrbd_s(blck,ch,l,d,:,:)] = mtspectrumc(squeeze(lfprespbd_s(thisinds,1001:1800))',params);
                    tbtlfpspectbd_c(blck,ch,l,d,:) = nanmean(lfpspecttbtbd_c(thisinds,:),1);
                    tbtlfpspectbd_s(blck,ch,l,d,:) = nanmean(lfpspecttbtbd_s(thisinds,:),1);
                    tbtlfpspecterrbd_c(blck,ch,l,d,:) = nanstd(lfpspecttbtbd_c(thisinds,:),1,1)./sqrt(length(thisinds));
                    tbtlfpspecterrbd_s(blck,ch,l,d,:) = nanstd(lfpspecttbtbd_s(thisinds,:),1,1)./sqrt(length(thisinds));                    
                    
                    [lfpcoherence(blck,ch,l,d,:),a,b,c,de,cfx,e,f,lfpcoherr(blck,ch,l,d,:,:)] = coherencyc(lfpresp_c(thisinds,respwin)',lfpresp_s(thisinds,respwin)',params);        
                    condtbtlfpcoherence(blck,ch,l,d,:) = nanmean(tbtlfpcoher(thisinds,:),1);
                    condtbtlfpcoherr(blck,ch,l,d,:) = nanstd(tbtlfpcoher(thisinds,:),1,1)./sqrt(length(thisinds));
                    
                    [lfpcoherencebd(blck,ch,l,d,:),a,b,c,de,cfx,e,f,lfpcoherrbd(blck,ch,l,d,:,:)] = coherencyc(lfprespbd_c(thisinds,respwin)',lfprespbd_s(thisinds,respwin)',params);  
                    condtbtlfpcoherencebd(blck,ch,l,d,:) = nanmean(tbtlfpcoherbd(thisinds,:),1);
                    condtbtlfpcoherrbd(blck,ch,l,d,:) = nanstd(tbtlfpcoherbd(thisinds,:),1,1)./sqrt(length(thisinds));
                                        
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
                        
                        [r1_lfpspectbd_c(blck,ch,l,d,:),chf,r1_lfpspecterrbd_c(blck,ch,l,d,:,:)] = mtspectrumc(squeeze(lfprespbd_c(thisruninds,1001:1800))',params);
                        [r1_lfpspectbd_s(blck,ch,l,d,:),chf,r1_lfpspecterrbd_s(blck,ch,l,d,:,:)] = mtspectrumc(squeeze(lfprespbd_s(thisruninds,1001:1800))',params);
                        r1_tbtlfpspectbd_c(blck,ch,l,d,:) = nanmean(lfpspecttbtbd_c(thisruninds,:),1);
                        r1_tbtlfpspectbd_s(blck,ch,l,d,:) = nanmean(lfpspecttbtbd_s(thisruninds,:),1);
                        r1_tbtlfpspecterrbd_c(blck,ch,l,d,:) = nanstd(lfpspecttbtbd_c(thisruninds,:),1,1)./sqrt(length(thisruninds));
                        r1_tbtlfpspecterrbd_s(blck,ch,l,d,:) = nanstd(lfpspecttbtbd_s(thisruninds,:),1,1)./sqrt(length(thisruninds));
                        [r1_lfpcoherencebd(blck,ch,l,d,:),a,b,c,de,cfx,e,f,r1_lfpcoherrbd(blck,ch,l,d,:,:)] = coherencyc(lfprespbd_c(thisruninds,respwin)',lfpresp_s(thisruninds,respwin)',params);
                        r1_tbtlfpcoherencebd(blck,ch,l,d,:) = nanmean(tbtlfpcoherbd(thisruninds,:),1);
                        r1_tbtlfpcoherrbd(blck,ch,l,d,:) = nanstd(tbtlfpcoherbd(thisruninds,:),1,1)./sqrt(length(thisruninds));
                        
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
                        
                        r1_lfpspectbd_c(blck,ch,l,d,:) = nan(1,length(lfpspect_c(blck,ch,l,d,:)));
                        r1_lfpspectbd_s(blck,ch,l,d,:) = nan(1,length(lfpspect_s(blck,ch,l,d,:)));
                        r1_lfpspecterrbd_c(blck,ch,l,d,:,:) = nan(2,length(lfpspect_c(blck,ch,l,d,:)));
                        r1_lfpspecterrbd_s(blck,ch,l,d,:,:) = nan(2,length(lfpspect_s(blck,ch,l,d,:)));
                        r1_tbtlfpspectbd_c(blck,ch,l,d,:) = nan(1,length(lfpspect_c(blck,ch,l,d,:)));
                        r1_tbtlfpspectbd_s(blck,ch,l,d,:) = nan(1,length(lfpspect_c(blck,ch,l,d,:)));
                        r1_tbtlfpspecterrbd_c(blck,ch,l,d,:) = nan(1,length(lfpspect_c(blck,ch,l,d,:)));
                        r1_tbtlfpspecterrbd_s(blck,ch,l,d,:) = nan(1,length(lfpspect_c(blck,ch,l,d,:)));
                        r1_lfpcoherencebd(blck,ch,l,d,:) = nan(1,length(lfpcoherence(blck,ch,l,d,:)));
                        r1_lfpcoherrbd(blck,ch,l,d,:,:) = nan(2,length(lfpcoherence(blck,ch,l,d,:)));
                        r1_tbtlfpcoherencebd(blck,ch,l,d,:) = nan(1,length(lfpcoherence(blck,ch,l,d,:)));
                        r1_tbtlfpcoherrbd(blck,ch,l,d,:) = nan(1,length(lfpcoherence(blck,ch,l,d,:)));
                        
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
                        r0_tbtlfpcoherr(blck,ch,l,d,:) = nanstd(tbtlfpcoher(thisstillinds,:),1)./sqrt(length(thisstillinds));
                        
                        [r0_lfpspectbd_c(blck,ch,l,d,:),chf,r0_lfpspecterrbd_c(blck,ch,l,d,:,:)] = mtspectrumc(squeeze(lfprespbd_c(thisstillinds,1001:1800))',params);
                        [r0_lfpspectbd_s(blck,ch,l,d,:),chf,r0_lfpspecterrbd_s(blck,ch,l,d,:,:)] = mtspectrumc(squeeze(lfprespbd_s(thisstillinds,1001:1800))',params);
                        r0_tbtlfpspectbd_c(blck,ch,l,d,:) = nanmean(lfpspecttbtbd_c(thisstillinds,:),1);
                        r0_tbtlfpspectbd_s(blck,ch,l,d,:) = nanmean(lfpspecttbtbd_s(thisstillinds,:),1);
                        r0_tbtlfpspecterrbd_c(blck,ch,l,d,:) = nanstd(lfpspecttbtbd_c(thisstillinds,:),1,1)./sqrt(length(thisstillinds));
                        r0_tbtlfpspecterrbd_s(blck,ch,l,d,:) = nanstd(lfpspecttbtbd_s(thisstillinds,:),1,1)./sqrt(length(thisstillinds));
                        [r0_lfpcoherencebd(blck,ch,l,d,:),a,b,c,de,cfx,e,f,r0_lfpcoherrbd(blck,ch,l,d,:,:)] = coherencyc(lfprespbd_c(thisstillinds,respwin)',lfprespbd_s(thisstillinds,respwin)',params);
                        r0_tbtlfpcoherencebd(blck,ch,l,d,:) = nanmean(tbtlfpcoherbd(thisstillinds,:),1);
                        r0_tbtlfpcoherrbd(blck,ch,l,d,:) = nanstd(tbtlfpcoherbd(thisstillinds,:),1,1)./sqrt(length(thisstillinds));
                        
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
                        r0_tbtlfpspect_c(blck,ch,l,d,:) = nan(1,length(lfpspect_c(blck,ch,l,d,:)));
                        r0_tbtlfpspect_s(blck,ch,l,d,:) = nan(1,length(lfpspect_c(blck,ch,l,d,:)));
                        r0_tbtlfpspecterr_c(blck,ch,l,d,:) = nan(1,length(lfpspect_c(blck,ch,l,d,:)));
                        r0_tbtlfpspecterr_s(blck,ch,l,d,:) = nan(1,length(lfpspect_c(blck,ch,l,d,:)));
                        r0_lfpcoherence(blck,ch,l,d,:) = nan(1,length(lfpcoherence(blck,ch,l,d,:)));
                        r0_lfpcoherr(blck,ch,l,d,:,:) = nan(2,length(lfpcoherence(blck,ch,l,d,:)));
                        r0_tbtlfpcoherence(blck,ch,l,d,:) = nan(1,length(lfpcoherence(blck,ch,l,d,:)));
                        r0_tbtlfpcoherr(blck,ch,l,d,:) = nan(1,length(lfpcoherence(blck,ch,l,d,:)));
                        
                        r0_lfpspectbd_c(blck,ch,l,d,:) = nan(1,length(lfpspect_c(blck,ch,l,d,:)));
                        r0_lfpspectbd_s(blck,ch,l,d,:) = nan(1,length(lfpspect_s(blck,ch,l,d,:)));
                        r0_lfpspecterrbd_c(blck,ch,l,d,:,:) = nan(2,length(lfpspect_c(blck,ch,l,d,:)));
                        r0_lfpspecterrbd_s(blck,ch,l,d,:,:) = nan(2,length(lfpspect_s(blck,ch,l,d,:)));
                        r0_tbtlfpspectbd_c(blck,ch,l,d,:) = nan(1,length(lfpspect_c(blck,ch,l,d,:)));
                        r0_tbtlfpspectbd_s(blck,ch,l,d,:) = nan(1,length(lfpspect_c(blck,ch,l,d,:)));
                        r0_tbtlfpspecterrbd_c(blck,ch,l,d,:) = nan(1,length(lfpspect_c(blck,ch,l,d,:)));
                        r0_tbtlfpspecterrbd_s(blck,ch,l,d,:) = nan(1,length(lfpspect_c(blck,ch,l,d,:)));
                        r0_lfpcoherencebd(blck,ch,l,d,:) = nan(1,length(lfpcoherence(blck,ch,l,d,:)));
                        r0_lfpcoherrbd(blck,ch,l,d,:,:) = nan(2,length(lfpcoherence(blck,ch,l,d,:)));
                        r0_tbtlfpcoherencebd(blck,ch,l,d,:) = nan(1,length(lfpcoherence(blck,ch,l,d,:)));
                        r0_tbtlfpcoherrbd(blck,ch,l,d,:) = nan(1,length(lfpcoherence(blck,ch,l,d,:)));
                        
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
    [c,di(i)] = min(abs(cdm(i,:)-evaldepth)); % get depth index, least distance to evaldepth
    
    lgind(i) = lgi(i,di(i));
end
valid = ones(size(di));
valid(8) = 0; % huge noise
valid = logical(valid);
vi = find(valid);
% vi = 1:length(centdepth);

for i = 1:length(vi)
%     lgpow(i,:,:) = squeeze(lfpspect_c(i,di(i),:,:,lgi(i,di(i)))); % at second contact only
    r1lgpow(i,:,:) = squeeze(r1_lfpspect_c(vi(i),di(vi(i)),:,:,lgi(vi(i),di(vi(i)))));
%     lglfpcoher(i,:,:) = squeeze(lfpcoherence(i,di(i),:,:,lgi(i,di(i))));
%     r1lglfpcoher(i,:,:) = squeeze(r1_lfpcoherence(i,di(i),:,:,lgi(i,di(i))));
    r1tbtlglfpcoher(i,:,:) = squeeze(r1_tbtlfpcoherence(vi(i),di(vi(i)),:,:,lgi(vi(i),di(vi(i)))));
%     lgsfcoher(i,:,:) = squeeze(sfcoher_c(i,2,:,:,lgi(i,2)));
%     r1lgsfcoher(i,:,:) = squeeze(r1_sfcoher_c(i,2,:,:,lgi(i,2)));

    % bipolar derivation coherence
    r1lgpowbd(i,:,:) = squeeze(r1_lfpspectbd_c(i,di(i),:,:,lgi(i,di(i))));
    r0lgpowbd(i,:,:) = squeeze(r0_lfpspectbd_c(i,di(i),:,:,lgi(i,di(i))));
    r1tbtlglfpcoherbd(i,:,:) = squeeze(r1_tbtlfpcoherencebd(i,di(i),:,:,lgi(i,di(i))));
    r0tbtlglfpcoherbd(i,:,:) = squeeze(r0_tbtlfpcoherencebd(i,di(i),:,:,lgi(i,di(i))));
    
    % with jack knifed error bars
    r1lglfpcoherbd(i,:,:) = squeeze(r1_lfpcoherencebd(i,di(i),:,:,lgi(i,di(i))));
    r0lglfpcoherbd(i,:,:) = squeeze(r0_lfpcoherencebd(i,di(i),:,:,lgi(i,di(i))));
    
    for l = 1:2
        for d = 1:4
            r1fillspecy(i,l,d,:) = [squeeze(r1_lfpspecterr_c(vi(i),di(vi(i)),l,d,1,12:103))',fliplr(squeeze(r1_lfpspecterr_c(vi(i),di(vi(i)),l,d,2,12:103))')];
            r1filltbtspecy(i,l,d,:) = [(squeeze(r1_tbtlfpspect_c(vi(i),di(vi(i)),l,d,12:103))+squeeze(r1_tbtlfpspecterr_c(vi(i),di(vi(i)),l,d,12:103)))',fliplr((squeeze(r1_tbtlfpspect_c(vi(i),di(vi(i)),l,d,12:103))-squeeze(r1_tbtlfpspecterr_c(vi(i),di(vi(i)),l,d,12:103)))')]';
            r1fillcohery(i,l,d,:) = [squeeze(r1_lfpcoherr(i,di(i),l,d,1,12:103))',fliplr(squeeze(r1_lfpcoherr(i,di(i),l,d,2,12:103))')]';
            r1filltbtcohery(i,l,d,:) = [(squeeze(r1_tbtlfpcoherence(vi(i),di(vi(i)),l,d,12:103))+squeeze(r1_tbtlfpcoherr(vi(i),di(vi(i)),l,d,12:103)))',fliplr((squeeze(r1_tbtlfpcoherence(vi(i),di(vi(i)),l,d,12:103))-squeeze(r1_tbtlfpcoherr(vi(i),di(vi(i)),l,d,12:103)))')]';
            
            r1fillspecybd(i,l,d,:) = [squeeze(r1_lfpspecterrbd_c(i,di(i),l,d,1,12:103))',fliplr(squeeze(r1_lfpspecterrbd_c(i,di(i),l,d,2,12:103))')];
            r1filltbtspecybd(i,l,d,:) = [(squeeze(r1_tbtlfpspectbd_c(i,di(i),l,d,12:103))+squeeze(r1_tbtlfpspecterrbd_c(i,di(i),l,d,12:103)))',fliplr((squeeze(r1_tbtlfpspectbd_c(i,di(i),l,d,12:103))-squeeze(r1_tbtlfpspecterrbd_c(i,di(i),l,d,12:103)))')]';
            r1fillcoherybd(i,l,d,:) = [squeeze(r1_lfpcoherrbd(i,di(i),l,d,1,12:103))',fliplr(squeeze(r1_lfpcoherrbd(i,di(i),l,d,2,12:103))')]';
            r1filltbtcoherybd(i,l,d,:) = [(squeeze(r1_tbtlfpcoherencebd(i,di(i),l,d,12:103))+squeeze(r1_tbtlfpcoherrbd(i,di(i),l,d,12:103)))',fliplr((squeeze(r1_tbtlfpcoherencebd(i,di(i),l,d,12:103))-squeeze(r1_tbtlfpcoherrbd(i,di(i),l,d,12:103)))')]';        
        end
    end    
end
fillx = [chf(12:103),fliplr(chf(12:103))];


% spectrum iso vs cross
for i =1:size(r1filltbtspecy,1)
    figure
    fill(fillx,squeeze(r1filltbtspecy(i,1,2,:)),[0,0,0])
    hold on
    fill(fillx,squeeze(r1filltbtspecy(i,1,3,:)),[.5,.5,.5])
    plot(chf(12:103), squeeze(r1_tbtlfpspect_c(vi(i),di(vi(i)),1,2,12:103)),'w')
    plot(chf(12:103), squeeze(r1_tbtlfpspect_c(vi(i),di(vi(i)),1,3,12:103)),'k')
    set(gca,'yscale','log')
    legend('same','offset')
    title(animalids{i})
end

% power iso vs cross scatter
[p,s,stats] = signrank(r1lgpow(:,1,2),r1lgpow(:,1,3));
figure
plot(r1lgpow(:,1,2),r1lgpow(:,1,3),'ko','markerfacecolor','k')
line([0,350],[0,350],'color','k')
axis square
xlabel('gamma power same phase')
ylabel('gamma power offset phase')
title(['signrank p = ' num2str(p)]);
n = size(r1lgpow,1);
disp(['same phase power: ' num2str(nanmean(r1lgpow(:,1,2))) '+-' num2str(std(r1lgpow(:,1,2))./sqrt(n))]);
disp(['diff phase power: ' num2str(nanmean(r1lgpow(:,1,3))) '+-' num2str(std(r1lgpow(:,1,3))./sqrt(n))]);
pcc = (1-r1lgpow(:,1,3)./r1lgpow(:,1,2)).*100;
disp(['percent change cross/iso = ' num2str(mean(pcc)) '+-' num2str(std(pcc)./sqrt(n))])

% power iso vs cross connected columns
figure
plot(1,r1lgpow(:,1,2),'ko','markerfacecolor','k')
hold on
plot(2,r1lgpow(:,1,3),'o','color',[.5,.5,.5],'markerfacecolor',[.5,.5,.5])
axis([0,3,0,300])
for i = 1:size(r1lgpow,1)    
    plot([1,2],[r1lgpow(:,1,2),r1lgpow(:,1,3)],'k')
end
set(gca,'xtick',[1,2]);
set(gca,'xticklabel',{'same','offset'})
title('change in gamma power')
ylabel('gamma power')
set(gcf,'OuterPosition',[186   300   240   513])


% spectrum L0 vs L1
for i =1:size(r1filltbtspecy,1)
    figure
    fill(fillx,squeeze(r1filltbtspecy(i,1,2,:)),[.3,.3,1])
    hold on
    fill(fillx,squeeze(r1filltbtspecy(i,2,2,:)),[1,.3,.3])
    set(gca,'yscale','log')
    legend('L0','L1')
    title(animalids{i})
end

% power L0 vs L1 scatter
[p,s] = signrank(r1lgpow(:,1,2),r1lgpow(:,2,2));
figure
plot(r1lgpow(:,1,2),r1lgpow(:,2,2),'ko','markerfacecolor','k')
line([0,350],[0,350],'color','k')
axis square
xlabel('gamma power L0')
ylabel('gamma power L1')
title(['signrank p = ' num2str(p)]);
n = size(r1lgpow,1);
disp(['same ori power: ' num2str(nanmean(r1lgpow(:,1,1))) '+-' num2str(std(r1lgpow(:,1,1))./sqrt(n))]);
disp(['diff ori power: ' num2str(nanmean(r1lgpow(:,2,1))) '+-' num2str(std(r1lgpow(:,2,1))./sqrt(n))]);
pcc = (1-r1lgpow(:,2,2)./r1lgpow(:,1,2)).*100;
disp(['percent change L2/L1 = ' num2str(mean(pcc)) '+-' num2str(std(pcc)./sqrt(n))])

% power L0 vs L1 connected columns
figure
plot(1,r1lgpow(:,1,2),'bo','markerfacecolor','b')
hold on
plot(2,r1lgpow(:,2,2),'ro','markerfacecolor','r')
axis([.5,2.5,0,300])
for i = 1:size(r1lgpow,1)    
    plot([1,2],[r1lgpow(:,1,2),r1lgpow(:,2,2)],'k')
end
set(gca,'xtick',[1,2]);
set(gca,'xticklabel',{'L0','L1'})
title('change in gamma power')
ylabel('gamma power')
set(gcf,'OuterPosition',[186   300   240   513])


% % coherence iso vs cross, trial by trial standard error
% for i =1:size(r1filltbtcohery,1)
%     figure
%     fill(fillx,squeeze(r1filltbtcohery(i,1,2,:)),[.3,.3,1])
%     hold on
%     fill(fillx,squeeze(r1filltbtcohery(i,1,3,:)),[1,.3,.3])
%     legend('iso','cross')
% end

% coherence iso vs cross, bi-deriv, jack knife confidence intervals
for i =1:size(r1fillcoherybd,1)
    figure
    fill(fillx,squeeze(r1fillcoherybd(i,1,2,:)),[.5,.5,1])
    hold on
    fill(fillx,squeeze(r1fillcoherybd(i,1,3,:)),[1,.5,.5])
    plot(chf(12:103),squeeze(r1_lfpcoherencebd(i,di(i),1,2,12:103)),'b')
    plot(chf(12:103),squeeze(r1_lfpcoherencebd(i,di(i),1,3,12:103)),'r')
%     fill(fillx,squeeze(r1fillcoherybd(i,1,3,:)),[.3,1,.3])
    legend('iso','cross')
end


% coherence iso vs cross scatter
[p,s] = signrank(r1lglfpcoherbd(:,1,2),r1lglfpcoherbd(:,1,3));
n = size(r1lglfpcoherbd,1); 
figure
plot(r1lglfpcoherbd(:,1,2),r1lglfpcoherbd(:,1,3),'ko','markerfacecolor','k')
line([0,1],[0,1],'color','k')
axis square
xlabel('gamma coherence same orientation')
ylabel('gamma coherence 90 off orientation')
title(['signrank p = ' num2str(p)]);
disp(['same ori coherency: ' num2str(nanmean(r1lglfpcoherbd(:,1,2))) '+-' num2str(std(r1lglfpcoherbd(:,1,2))./sqrt(n))]);
disp(['diff ori coherency: ' num2str(nanmean(r1lglfpcoherbd(:,1,3))) '+-' num2str(std(r1lglfpcoherbd(:,1,3))./sqrt(n))]);

%coherence iso vs. cross connected columns
figure
plot(1,r1tbtlglfpcoher(:,1,2),'ko','markerfacecolor','k')
hold on
plot(2,r1tbtlglfpcoher(:,1,3),'o','color',[.5,.5,.5],'markerfacecolor',[.5,.5,.5])
axis([.5,2.5,.25,.9])
for i = 1:size(r1tbtlglfpcoher,1)    
    plot([1,2],[r1tbtlglfpcoher(:,1,2),r1tbtlglfpcoher(:,1,3)],'k')
end
set(gca,'xtick',[1,2]);
set(gca,'xticklabel',{'iso','cross'})
title('coherence between sites')
ylabel('coherence')
set(gcf,'OuterPosition',[186   438   240   513])


% coherence light off vs on, trial by trial standard error
for i =1:size(r1filltbtcohery,1)
    figure
    fill(fillx,squeeze(r1filltbtcohery(i,1,2,:)),[.3,.3,1])
    hold on
    fill(fillx,squeeze(r1filltbtcohery(i,2,2,:)),[1,.3,.3])
    legend('L0','L1')
end

% coherence L0 vs L1 scatter
[p,s] = signrank(r1tbtlglfpcoher(:,1,2),r1tbtlglfpcoher(:,2,2));
n = size(r1tbtlglfpcoher,1); 
figure
plot(r1tbtlglfpcoher(:,1,2),r1tbtlglfpcoher(:,2,2),'ko','markerfacecolor','k')
line([0,1],[0,1],'color','k')
axis square
xlabel('gamma power same orientation')
ylabel('gamma power 90 off orientation')
title(['signrank p = ' num2str(p)]);
disp(['same ori L0 coherency: ' num2str(nanmean(r1tbtlglfpcoher(:,1,2))) '+-' num2str(std(r1tbtlglfpcoher(:,1,2))./sqrt(n))]);
disp(['same ori L1 coherency: ' num2str(nanmean(r1tbtlglfpcoher(:,2,2))) '+-' num2str(std(r1tbtlglfpcoher(:,2,2))./sqrt(n))]);

% coherence L0 vs L1 connected columns
figure
plot(1,r1tbtlglfpcoher(:,1,2),'ko','markerfacecolor','k')
hold on
plot(2,r1tbtlglfpcoher(:,2,2),'ro','markerfacecolor','r')
axis([.5,2.5,.25,.9])
for i = 1:size(r1tbtlglfpcoher,1)    
    plot([1,2],[r1tbtlglfpcoher(:,1,2),r1tbtlglfpcoher(:,2,2)],'k')
end
set(gca,'xtick',[1,2]);
set(gca,'xticklabel',{'L0','L1'})
title('coherence between sites')
ylabel('coherence')
set(gcf,'OuterPosition',[186   438   240   513])


%same for running
% 1 center only, 2 same phase 3 offset phase 4 no center
for i = 1:length(blocks)
    figure('name',animalids{i});
    subplot(2,2,1)
    semilogy(chf,squeeze(r1_lfpspect_c(i,3,1,2,:)),'g')
    hold on    
    semilogy(chf,squeeze(r1_lfpspect_c(i,3,1,3,:)),'r')
    semilogy(chf,squeeze(r1_lfpspect_c(i,3,1,1,:)),'b')
    semilogy(chf,squeeze(r1_lfpspect_c(i,3,1,4,:)),'c')
    ax = axis;
    axis([0,120,1,ax(4)])
    legend('same','90deg','on field only','off field only')
    xlabel('frequency')
    ylabel('psd')
    
    subplot(2,2,2)
    plot(chf,squeeze(r1_lfpcoherence(i,3,1,2,:)),'g');
    hold on
    plot(chf,squeeze(r1_lfpcoherence(i,3,1,3,:)),'r');
    plot(chf,squeeze(r1_lfpcoherence(i,3,1,1,:)),'b');
    plot(chf,squeeze(r1_lfpcoherence(i,3,1,4,:)),'c');
    axis([0,120,0,1]);
    xlabel('frequency')
    ylabel('site coherence')
    
%     subplot(2,3,3)
%     plot(chf,squeeze(r1_sfcoher_c(i,3,1,2,:)),'g');
%     hold on
%     plot(chf,squeeze(r1_sfcoher_c(i,3,1,3,:)),'r');
%     plot(chf,squeeze(r1_sfcoher_c(i,3,1,1,:)),'b');
%     plot(chf,squeeze(r1_sfcoher_c(i,3,1,4,:)),'c');
%     axis([0,120,0,1]);
%     xlabel('frequency')
%     ylabel('spike field coherence')
    
    subplot(2,2,3)
    semilogy(chf,squeeze(r1_lfpspect_c(i,3,1,1,:)),'b')
    hold on    
    semilogy(chf,squeeze(r1_lfpspect_c(i,3,1,2,:)),'c')
    semilogy(chf,squeeze(r1_lfpspect_c(i,3,2,1,:)),'r')
    semilogy(chf,squeeze(r1_lfpspect_c(i,3,2,2,:)),'m')
    ax = axis;
    axis([0,120,1,ax(4)])
    legend('small l0','same l0','small l1','same l1')
    xlabel('frequency')
    ylabel('psd')
    
    subplot(2,2,4)
    plot(chf,squeeze(r1_lfpcoherence(i,3,1,1,:)),'b');
    hold on
    plot(chf,squeeze(r1_lfpcoherence(i,3,1,2,:)),'c');
    plot(chf,squeeze(r1_lfpcoherence(i,3,2,1,:)),'r');
    plot(chf,squeeze(r1_lfpcoherence(i,3,2,2,:)),'m');
    axis([0,120,0,1]);
    xlabel('frequency')
    ylabel('site coherence')
    
%     subplot(2,3,6)
%     plot(chf,squeeze(r1_sfcoher_c(i,3,1,1,:)),'b');
%     hold on
%     plot(chf,squeeze(r1_sfcoher_c(i,3,1,2,:)),'c');
%     plot(chf,squeeze(r1_sfcoher_c(i,3,2,1,:)),'r');
%     plot(chf,squeeze(r1_sfcoher_c(i,3,2,2,:)),'m');
%     axis([0,120,0,1]);
%     xlabel('frequency')
%     ylabel('spike field coherence')
    
    set(gcf,'OuterPosition',[50,100,1200,700])
end

% 1 center only, 2 same phase 3 offset phase 4 no center
for i = 1:length(blocks)
    
    depthax = centdepth(i):-spacing(i):centdepth(i)-(length(centelectrodes(i,1):centelectrodes(i,2))-1)*spacing(i);
    
    figure
    subplot(3,3,1)
    imagesc(chf,depthax,log(squeeze(lfpspect_c(i,:,1,1,:))))
    axis([1,120,depthax(end)-10,depthax(1)+10])
    caxis([-2,9])
    title('power center')
    
    subplot(3,3,2)
    imagesc(chf,depthax,log(squeeze(lfpspect_c(i,:,1,2,:))))
    axis([1,120,depthax(end)-10,depthax(1)+10])
    caxis([-2,9])
    title('power iso')
    
    subplot(3,3,3)
    imagesc(chf,depthax,log(squeeze(lfpspect_c(i,:,1,3,:))))
    axis([1,120,depthax(end)-10,depthax(1)+10])
    caxis([-2,9])
    title('power cross')
    
    subplot(3,3,4)
    imagesc(chf,depthax,log(squeeze(lfpspect_c(i,:,1,2,:)))-log(squeeze(lfpspect_c(i,:,1,1,:))));
    axis([1,120,depthax(end)-10,depthax(1)+10])
    caxis([-2,2])
    title('large - small')
    
    subplot(3,3,5)
    imagesc(chf,depthax,log(squeeze(lfpspect_c(i,:,1,3,:)))-log(squeeze(lfpspect_c(i,:,1,2,:))));
    axis([1,120,depthax(end)-10,depthax(1)+10])
    caxis([-2,2])
    title('iso - cross')    
    
    subplot(3,3,7)
    imagesc(chf,depthax,log(squeeze(lfpspect_c(i,:,2,1,:)))-log(squeeze(lfpspect_c(i,:,1,1,:))));
    axis([1,120,depthax(end)-10,depthax(1)+10])
    caxis([-2,2])
    title('small L1-L0')
    
    subplot(3,3,8)
    imagesc(chf,depthax,log(squeeze(lfpspect_c(i,:,2,2,:)))-log(squeeze(lfpspect_c(i,:,1,2,:))));
    axis([1,120,depthax(end)-10,depthax(1)+10])
    caxis([-2,2])
    title('iso L1-L0')
    
    subplot(3,3,9)
    imagesc(chf,depthax,log(squeeze(lfpspect_c(i,:,2,3,:)))-log(squeeze(lfpspect_c(i,:,1,3,:))));
    axis([1,120,depthax(end)-10,depthax(1)+10])
    caxis([-2,2])
    title('cross L1-L0')
    
    
    figure
    subplot(3,3,1)
    imagesc(chf,depthax,squeeze(lfpcoherence(i,:,1,1,:)))
    axis([1,120,depthax(end)-10,depthax(1)+10])
    caxis([0,1])
    title('coherence RF only')
    
    subplot(3,3,2)
    imagesc(chf,depthax,squeeze(lfpcoherence(i,:,1,2,:)))
    axis([1,120,depthax(end)-10,depthax(1)+10])
    caxis([0,1])
    title('match surround')    
    
    subplot(3,3,3)
    imagesc(chf,depthax,squeeze(lfpcoherence(i,:,1,3,:)))
    axis([1,120,depthax(end)-10,depthax(1)+10])
    caxis([0,1])
    title('rotated surround')    
    
    subplot(3,3,4)
    imagesc(chf,depthax,squeeze(lfpcoherence(i,:,1,2,:))-squeeze(lfpcoherence(i,:,1,1,:)))
    axis([1,120,depthax(end)-10,depthax(1)+10])
    caxis([-.5,.5])
    title('iso large - small')    
    
    subplot(3,3,5)
    imagesc(chf,depthax,squeeze(lfpcoherence(i,:,1,3,:))-squeeze(lfpcoherence(i,:,1,2,:)))
    axis([1,120,depthax(end)-10,depthax(1)+10])
    caxis([-.3,.3])
    title('cross surround - iso surround')
    
    subplot(3,3,7)
    imagesc(chf,depthax,squeeze(lfpcoherence(i,:,2,1,:))-squeeze(lfpcoherence(i,:,1,1,:)))
    axis([1,120,depthax(end)-10,depthax(1)+10])
    caxis([-.3,.3])
    title('small L1-L0')       
    
    subplot(3,3,8)
    imagesc(chf,depthax,squeeze(lfpcoherence(i,:,2,2,:))-squeeze(lfpcoherence(i,:,1,2,:)))
    axis([1,120,depthax(end)-10,depthax(1)+10])
    caxis([-.3,.3])
    title('iso L1-L0')    
    
    subplot(3,3,9)
    imagesc(chf,depthax,squeeze(lfpcoherence(i,:,2,3,:))-squeeze(lfpcoherence(i,:,1,3,:)))
    axis([1,120,depthax(end)-10,depthax(1)+10])
    caxis([-.3,.3])
    title('cross L1-L0')    
    
%     figure
%     subplot(3,3,1)
%     imagesc(chf,depthax,squeeze(sfcoher_c(i,:,1,1,:)))
%     axis([1,120,depthax(end)-10,depthax(1)+10])
%     caxis([0,.6])
%     title('SFC small')
%     
%     subplot(3,3,2)
%     imagesc(chf,depthax,squeeze(sfcoher_c(i,:,1,2,:)))
%     axis([1,120,depthax(end)-10,depthax(1)+10])
%     caxis([0,.6])
%     title('SFC iso')
%     
%     subplot(3,3,3)
%     imagesc(chf,depthax,squeeze(sfcoher_c(i,:,1,3,:)))
%     axis([1,120,depthax(end)-10,depthax(1)+10])
%     caxis([0,.6])
%     title('SFC cross')    
%     
%     subplot(3,3,4)
%     imagesc(chf,depthax,squeeze(sfcoher_c(i,:,1,2,:))-squeeze(sfcoher_c(i,:,1,1,:)))
%     axis([1,120,depthax(end)-10,depthax(1)+10])
%     caxis([-.2,.2])
%     title('SFC large - small')
%     
%     subplot(3,3,5)
%     imagesc(chf,depthax,squeeze(sfcoher_c(i,:,1,3,:))-squeeze(sfcoher_c(i,:,1,2,:)))
%     axis([1,120,depthax(end)-10,depthax(1)+10])
%     caxis([-.2,.2])
%     title('SFC cross - iso')
%     
%     subplot(3,3,7)    
%     imagesc(chf,depthax,squeeze(sfcoher_c(i,:,2,1,:))-squeeze(sfcoher_c(i,:,1,1,:)))
%     axis([1,120,depthax(end)-10,depthax(1)+10])
%     caxis([-.2,.2])
%     title('small L1-L0')
%     
%     subplot(3,3,8)    
%     imagesc(chf,depthax,squeeze(sfcoher_c(i,:,2,2,:))-squeeze(sfcoher_c(i,:,1,2,:)))
%     axis([1,120,depthax(end)-10,depthax(1)+10])
%     caxis([-.2,.2])
%     title('iso L1-L0')
%     
%     subplot(3,3,9)    
%     imagesc(chf,depthax,squeeze(sfcoher_c(i,:,2,3,:))-squeeze(sfcoher_c(i,:,1,3,:)))
%     axis([1,120,depthax(end)-10,depthax(1)+10])
%     caxis([-.2,.2])
%     title('cross L1-L0')

end

 disp('');