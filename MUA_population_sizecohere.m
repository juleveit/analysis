function MUA_population_sizecohere


% % PV Halo population
% animalids =    {'150730', '150731', '150804', '150818', '150820', '150824', '151104'};
% blocks    =     [4,        3,         4,       3,        3,        3,        5];
% animal    =     [1,        2,         3,       4,        5,        6,        7];
% centelectrodes =[[1,16];   [1,16];    [1,16];  [1,16];   [1,16];   [17,32];  [17,32]];
% centdepth    =  [400,      400,       400,     400,      400,      400,      400];
% surrelectrodes =[[17,32];  [17,32];   [17,32]; [17,32];  [17,32];  [1,16];   [1,16]];
% surrdepth    =  [400,       400,     400,      400,      400,      400,      400];
% penangle =      [25,        25,      25,       25,       25,       25,       25];
% spacing =       [25,        25,      25,       25,       25,       25,       25];
% popfile = 'C:\Users\Julia\work\data\populations\PV_Halo\size\MUA_population.mat';

% % SOM Halo population
% animalids =    {'150825', '150831', '150902', '150909', '150915', '151023', '151027', '151109', '151110'};
% blocks    =     [5,         4,       3,        4,        3,        14,       3,        11,       13];
% animal    =     [1,         2,       3,        4,        5,        6,        7,        8,        8];
% centelectrodes =[[17,32];   [1,16];  [1,16];   [1,16];   [17,32];  [17,32];  [17,32];  [17,32];  [1,16]];
% centdepth    =  [400,       400,     400,      400,      400,      400,      500,      500,      500];
% surrelectrodes =[[1,16];    [17,32]; [17,32];  [17,32];  [1,16];   [1,16];   [1,16];   [1,16];   [17,32]];
% surrdepth    =  [400,       400,     400,      400,      400,      400,      500,      500,      500];
% penangle =      [25,        25,      25,       25,       25,       25,       25,       25,       25];
% spacing =       [25,        25,      25,       25,       25,       25,       25,       25,       25];
% popfile = 'C:\Users\Julia\work\data\populations\SOM_Halo_later\size\MUA_population.mat';

% % SOM Halo anesthetized pop
% animalids =    {'150826', '151230', '151231'};
% blocks    =     [3,        3,        11];
% animal    =     [1,        2,        3];
% centelectrodes =[[1,16];  [1,16];   [1,16]];
% centdepth    =  [500,      500,      500];
% surrelectrodes =[[17,32]; [17,32];  [17,32]];
% surrdepth    =  [500,      500,      500];
% penangle =      [25,       25,       25];
% spacing =       [25,       25,       25];
% popfile = 'C:\Users\Julia\work\data\populations\SOM_Halo_anesth\size\MUA_population.mat';

% % Scnn eArch pop
% animalids =    {'160219', '160222', '160223'};
% blocks    =     [4,        2,        7];
% animal    =     [1,        2,        3];
% centelectrodes =[[1,16];  [1,16];   [1,16]];
% centdepth    =  [650,      550,      400];
% surrelectrodes =[[17,32]; [17,32];  [17,32]];
% surrdepth    =  [500,      500,      500];
% penangle =      [25,       25,       25];
% spacing =       [25,       25,       25];
% popfile = 'C:\Users\Julia\work\data\populations\Scnn_eArch\size\MUA_population.mat';

% single thing
animalids =    {'171205'};
blocks    =     [7];
animal    =     [1];
centelectrodes =[[1,16]];
centdepth    =  [500];
surrelectrodes =[[33,48]];
surrdepth    =  [500];
penangle =      [30];
spacing =       [25];
popfile = 'C:\Users\Julia\work\data\populations\temp2.mat';



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
        
        % fix so it is usable with new multi-purpose grating stim script
        if isfield(centresult, 'sizeconds')
            allinds = sort(getSpecificIndices(centresult, 'sizeconds'));
            msstamps = centresult.msstamps(allinds);
            light = centresult.light(allinds);
            gratingInfo.Orientation = centresult.gratingInfo.Orientation(allinds);
            gratingInfo.size = centresult.gratingInfo.size(allinds);
            gratingInfo.Contrast = centresult.gratingInfo.Contrast(allinds);
            gratingInfo.tFreq = centresult.gratingInfo.tFreq(allinds);
        else
            msstamps = centresult.msstamps;
            light = centresult.light;
            gratingInfo = centresult.gratingInfo;
        end
          
        for ch = 1: length(centelectrodes(blck,1):centelectrodes(blck,2))
            
            disp(['Block ' int2str(blck) '/' int2str(length(blocks)) '   channel ' int2str(ch) '/' int2str(length(centelectrodes(blck,1):centelectrodes(blck,2)))]);
            
            depthc(ch) = centdepth(blck)-(ch-1)*spacing(blck);
            depths(ch) = surrdepth(blck)-(ch-1)*spacing(blck);   
            
            trialdur = centresult.stimduration*1000;
%             msstamps = centresult.msstamps;
            if length(msstamps)~=length(light)
%                 msstamps([193,291]) = []; % for 151023 block 5
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
            okspeed = find(meanspeed>( mean(meanspeed(notstill))-(1.5*std(meanspeed(notstill))) ) & meanspeed > 1);
            okvar = find(stdspeed<( mean(stdspeed(notstill))+(1.5*std(stdspeed(notstill)))) & stdspeed>.5);
            oktrials = intersect(okspeed,okvar);
            nonoktrials = 1:size(speed,1); nonoktrials(oktrials) = [];
            stilltrials = 1:size(speed,1); stilltrials(notstill) = [];
            
            % find low and high gamma peaks in spectra
            beta = [20,40];
            gamma = [50,70];
            sizes = unique(gratingInfo.size); sizes(sizes==0) = [];
            large = find(gratingInfo.size == max(sizes) & light == 0);
            small = find(gratingInfo.size == min(sizes) & light == 0);
            lr1 = intersect(large,oktrials);
            sr1 = intersect(small,oktrials);
            for i = 1:length(large)
                [pl(i,:),f] = pmtm(centresult.lfp(ch,msstamps(large(i))+300:msstamps(large(i))+1300),3,[],1000);
                [ps(i,:),f] = pmtm(centresult.lfp(ch,msstamps(small(i))+300:msstamps(small(i))+1300),3,[],1000);
            end            
            plr1 = nan(1,size(pl,2)); psr1 = nan(1,size(pl,2));
            if length(lr1>=5)
                for i = 1:length(lr1)
                    [plr1(i,:),f] = mtspectrumc(centresult.lfp(ch,msstamps(lr1(i))+700:msstamps(lr1(i))+1500),params);
%                     [plr1(i,:),f] = pmtm(centresult.lfp(ch,centresult.msstamps(lr1(i))+300:centresult.msstamps(lr1(i))+1300),3,[],1000);
                end
            end
            if length(sr1>=5)
                for i = 1:length(sr1)
                    [psr1(i,:),f] = mtspectrumc(centresult.lfp(ch,msstamps(sr1(i))+700:msstamps(sr1(i))+1500),params);
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
                
            end
           
            frs_c = sum(resp_c(:,respwin),2)./(length(respwin)/1000);
            frs_s = sum(resp_s(:,respwin),2)./(length(respwin)/1000);
            bl_c = sum(resp_c(:,1:prestim),2)./(prestim/1000);
            bl_s = sum(resp_s(:,1:prestim),2)./(prestim/1000);
            muacmean_c = mean(muacresp_c(:,respwin),2);
            muacmean_s = mean(muacresp_s(:,respwin),2);

            for l = 1:length(unique(light))
                for d = 1:length(sizes)
                    
                    thisinds = find(gratingInfo.size == sizes(d) & light == l-1);
                    thisruninds = intersect(thisinds,oktrials); thisstillinds = intersect(thisinds,stilltrials);
                  
                    condresp_c(blck,ch,l,d,:) = nanmean(resp_c(thisinds,:),1);
                    condresp_s(blck,ch,l,d,:) = nanmean(resp_s(thisinds,:),1);
                    condfr_c(blck,ch,l,d) = nanmean(frs_c(thisinds));
                    condfr_s(blck,ch,l,d) = nanmean(frs_s(thisinds));
                    condmuacresp_c(blck,ch,l,d,:) = nanmean(muacresp_c(thisinds,:));
                    condmuacresp_s(blck,ch,l,d,:) = nanmean(muacresp_s(thisinds,:));
                    condmua_c(blck,ch,l,d) = nanmean(muacmean_c(thisinds));
                    condmua_s(blck,ch,l,d) = nanmean(muacmean_s(thisinds));
                    condlfpresp_c(blck,ch,l,d,:) = nanmean(lfpresp_c(thisinds,:));
                    condlfpresp_s(blck,ch,l,d,:) = nanmean(lfpresp_s(thisinds,:));
                    
                    [lfpspect_c(blck,ch,l,d,:),chf,lfpspecterr_c(blck,ch,l,d,:,:)] = mtspectrumc(squeeze(lfpresp_c(thisinds,1001:1800))',params);
                    [lfpspect_s(blck,ch,l,d,:),chf,lfpspecterr_s(blck,ch,l,d,:,:)] = mtspectrumc(squeeze(lfpresp_s(thisinds,1001:1800))',params);
                    tbtlfpspect_c(blck,ch,l,d,:) = nanmean(lfpspecttbt_c(thisinds,:));
                    tbtlfpspect_s(blck,ch,l,d,:) = nanmean(lfpspecttbt_s(thisinds,:));
                    tbtlfpspecterr_c(blck,ch,l,d,:) = nanstd(lfpspecttbt_c(thisinds,:))./sqrt(length(thisinds));
                    tbtlfpspecterr_s(blck,ch,l,d,:) = nanstd(lfpspecttbt_s(thisinds,:))./sqrt(length(thisinds));
                    [muacspect_c(blck,ch,l,d,:),chf,muacspecterr_c(blck,ch,l,d,:,:)] = mtspectrumc(squeeze(muacresp_c(thisinds,1001:1800))',params);
                    [muacspect_s(blck,ch,l,d,:),chf,muacspecterr_s(blck,ch,l,d,:,:)] = mtspectrumc(squeeze(muacresp_s(thisinds,1001:1800))',params);
                    
                    
                    [lfpcoherence(blck,ch,l,d,:),a,b,c,de,cfx,e,f,lfpcoherr(blck,ch,l,d,:,:)] = coherencyc(lfpresp_c(thisinds,respwin)',lfpresp_s(thisinds,respwin)',params);        
                    condtbtlfpcoherence(blck,ch,l,d,:) = nanmean(tbtlfpcoher(thisinds,:));
                    condtbtlfpcoherr(blck,ch,l,d,:) = nanstd(tbtlfpcoher(thisinds,:))./sqrt(length(thisinds));     
                    condxc(blck,ch,l,d,:) = squeeze(nanmean(xc(thisinds,:)));
                    condmcxc(blck,ch,l,d,:) = squeeze(nanmean(mcxc(thisinds,:)));
%                     [sfcoher_c(blck,ch,l,d,:),a,b,c,de,sfcfx,e,f,g,sfcoherr_c(blck,ch,l,d,:,:)] = coherencycpt(lfpresp_c(thisinds,1001:1800)',ptresp_c(thisinds),params);
%                     [sfcoher_s(blck,ch,l,d,:),a,b,c,de,sfcfx,e,f,g,sfcoherr_s(blck,ch,l,d,:,:)] = coherencycpt(lfpresp_s(thisinds,1001:1800)',ptresp_s(thisinds),params);
                    condcorstrength(blck,ch,l,d,:) = nanmean(corstrength(thisinds));
                    
                    if ~isempty(thisruninds)
                        [r1_lfpspect_c(blck,ch,l,d,:),chf,r1_lfpspecterr_c(blck,ch,l,d,:,:)] = mtspectrumc(squeeze(lfpresp_c(thisruninds,1001:1800))',params);
                        [r1_lfpspect_s(blck,ch,l,d,:),chf,r1_lfpspecterr_s(blck,ch,l,d,:,:)] = mtspectrumc(squeeze(lfpresp_s(thisruninds,1001:1800))',params);
                        r1_tbtlfpspect_c(blck,ch,l,d,:) = nanmean(lfpspecttbt_c(thisruninds,:));
                        r1_tbtlfpspect_s(blck,ch,l,d,:) = nanmean(lfpspecttbt_s(thisruninds,:));
                        r1_tbtlfpspecterr_c(blck,ch,l,d,:) = nanstd(lfpspecttbt_c(thisruninds,:))./sqrt(length(thisruninds));
                        r1_tbtlfpspecterr_s(blck,ch,l,d,:) = nanstd(lfpspecttbt_s(thisruninds,:))./sqrt(length(thisruninds));
                        [r1_lfpcoherence(blck,ch,l,d,:),a,b,c,de,cfx,e,f,r1_lfpcoherr(blck,ch,l,d,:,:)] = coherencyc(lfpresp_c(thisruninds,respwin)',lfpresp_s(thisruninds,respwin)',params);    
                        r1_tbtlfpcoherence(blck,ch,l,d,:) = nanmean(tbtlfpcoher(thisruninds,:));
                        r1_tbtlfpcoherr(blck,ch,l,d,:) = nanstd(tbtlfpcoher(thisruninds,:))./sqrt(length(thisruninds));                   
                        r1condxc(blck,ch,l,d,:) = squeeze(nanmean(xc(thisruninds,:)));
%                         [r1_sfcoher_c(blck,ch,l,d,:),a,b,c,de,sfcfx,e,f,g,r1_sfcoherr_c(blck,ch,l,d,:,:)] = coherencycpt(lfpresp_c(thisruninds,1001:1800)',ptresp_c(thisruninds),params);
%                         [r1_sfcoher_s(blck,ch,l,d,:),a,b,c,de,sfcfx,e,f,g,r1_sfcoherr_s(blck,ch,l,d,:,:)] = coherencycpt(lfpresp_s(thisruninds,1001:1800)',ptresp_s(thisruninds),params);
                        r1_condcorstrength(blck,ch,l,d) = nanmean(corstrength(thisruninds));
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
                        r0_tbtlfpspect_c(blck,ch,l,d,:) = nanmean(lfpspecttbt_c(thisstillinds,:));
                        r0_tbtlfpspect_s(blck,ch,l,d,:) = nanmean(lfpspecttbt_s(thisstillinds,:));
                        r0_tbtlfpspecterr_c(blck,ch,l,d,:) = nanstd(lfpspecttbt_c(thisstillinds,:))./sqrt(length(thisstillinds));
                        r0_tbtlfpspecterr_s(blck,ch,l,d,:) = nanstd(lfpspecttbt_s(thisstillinds,:))./sqrt(length(thisstillinds));
                        [r0_lfpcoherence(blck,ch,l,d,:),a,b,c,de,cfx,e,f,r0_lfpcoherr(blck,ch,l,d,:,:)] = coherencyc(lfpresp_c(thisstillinds,respwin)',lfpresp_s(thisstillinds,respwin)',params);
                        r0_tbtlfpcoherence(blck,ch,l,d,:) = nanmean(tbtlfpcoher(thisstillinds,:));
                        r0_tbtlfpcoherr(blck,ch,l,d,:) = nanstd(tbtlfpcoher(thisstillinds,:))./sqrt(length(thisstillinds));
                        r0condxc(blck,ch,l,d,:) = squeeze(nanmean(xc(thisstillinds,:)));
%                         [r0_sfcoher_c(blck,ch,l,d,:),a,b,c,de,sfcfx,e,f,g,r0_sfcoherr_c(blck,ch,l,d,:,:)] = coherencycpt(lfpresp_c(thisstillinds,1001:1800)',ptresp_c(thisstillinds),params);
%                         [r0_sfcoher_s(blck,ch,l,d,:),a,b,c,de,sfcfx,e,f,g,r0_sfcoherr_s(blck,ch,l,d,:,:)] = coherencycpt(lfpresp_s(thisstillinds,1001:1800)',ptresp_s(thisstillinds),params);
                        r0_condcorstrength(blck,ch,l,d,:) = nanmean(corstrength(thisstillinds));
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
                        r0condxc(blck,ch,l,d,:) = nan(1,size(xc,2));
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
%     valid(i) = leffect(i,di(i));              % valid if light has a significant effect on large size stimuli
end
% valid(isnan(valid)) = 0;
valid = ones(size(di));
valid = logical(valid);
vi = find(valid);

for i = 1:length(vi)
%     lgpow(i,:,:) = squeeze(lfpspect_c(vi(i),di(vi(i)),:,:,lgi(vi(i),di(vi(i))))); % at second contact only
    r1lgpow(i,:,:) = squeeze(r1_lfpspect_c(vi(i),di(vi(i)),:,:,lgi(vi(i),di(vi(i)))));
    r0lgpow(i,:,:) = squeeze(r0_lfpspect_c(vi(i),di(vi(i)),:,:,lgi(vi(i),di(vi(i)))));
%     lglfpcoher(i,:,:) = squeeze(lfpcoherence(vi(i),di(vi(i)),:,:,lgi(vi(i),di(vi(i)))));
%     r1lglfpcoher(i,:,:) = squeeze(r1_lfpcoherence(vi(i),di(vi(i)),:,:,lgi(vi(i),di(vi(i)))));
    r1tbtlglfpcoher(i,:,:) = squeeze(r1_tbtlfpcoherence(vi(i),di(vi(i)),:,:,lgi(vi(i),di(vi(i)))));
    r0tbtlglfpcoher(i,:,:) = squeeze(r0_tbtlfpcoherence(vi(i),di(vi(i)),:,:,lgi(vi(i),di(vi(i)))));
%     lgsfcoher(i,:,:) = squeeze(sfcoher_c(vi(i),di(vi(i)),:,:,lgi(vi(i),di(vi(i)))));
%     r1lgsfcoher(i,:,:) = squeeze(r1_sfcoher_c(vi(i),di(vi(i)),:,:,lgi(vi(i),di(vi(i)))));
    
    for l = 1:2
        for d = 1:5
            r1fillspecy(i,l,d,:) = [squeeze(r1_lfpspecterr_c(vi(i),di(vi(i)),l,d,1,1:104))',fliplr(squeeze(r1_lfpspecterr_c(vi(i),2,l,d,2,1:104))')];            
            r1filltbtspecy(i,l,d,:) = [(squeeze(r1_tbtlfpspect_c(i,di(i),l,d,1:104))+squeeze(r1_tbtlfpspecterr_c(i,di(i),l,d,1:104)))',fliplr((squeeze(r1_tbtlfpspect_c(i,di(i),l,d,1:104))-squeeze(r1_tbtlfpspecterr_c(i,di(i),l,d,1:104)))')]';
%             r1fillcohery(i,l,d,:) = [squeeze(r1_lfpcoherr(vi(i),di(vi(i)),l,d,1,1:104))',fliplr(squeeze(r1_lfpcoherr(vi(i),2,l,d,2,1:104))')]';
            r1filltbtcohery(i,l,d,:) = [(squeeze(r1_tbtlfpcoherence(vi(i),di(vi(i)),l,d,1:104))+squeeze(r1_tbtlfpcoherr(vi(i),di(vi(i)),l,d,1:104)))',fliplr((squeeze(r1_tbtlfpcoherence(vi(i),di(vi(i)),l,d,1:104))-squeeze(r1_tbtlfpcoherr(vi(i),di(vi(i)),l,d,1:104)))')]';
            r1fillspecy(i,l,d,:) = [squeeze(r1_lfpspecterr_c(vi(i),di(vi(i)),l,d,1,1:104))',fliplr(squeeze(r1_lfpspecterr_c(vi(i),2,l,d,2,1:104))')];            
            r0filltbtspecy(i,l,d,:) = [(squeeze(r0_tbtlfpspect_c(i,di(i),l,d,1:104))+squeeze(r0_tbtlfpspecterr_c(i,di(i),l,d,1:104)))',fliplr((squeeze(r0_tbtlfpspect_c(i,di(i),l,d,1:104))-squeeze(r0_tbtlfpspecterr_c(i,di(i),l,d,1:104)))')]';
            r0filltbtcohery(i,l,d,:) = [(squeeze(r0_tbtlfpcoherence(vi(i),di(vi(i)),l,d,1:104))+squeeze(r0_tbtlfpcoherr(vi(i),di(vi(i)),l,d,1:104)))',fliplr((squeeze(r0_tbtlfpcoherence(vi(i),di(vi(i)),l,d,1:104))-squeeze(r0_tbtlfpcoherr(vi(i),di(vi(i)),l,d,1:104)))')]';
        end
    end    
end
fillx = [chf(1:104),fliplr(chf(1:104))];


% spectrum small vs large
for i =1:size(r1filltbtspecy,1)
    figure
    fill(fillx,squeeze(r1filltbtspecy(i,1,1,:)),[.7,.7,1])
    hold on
    fill(fillx,squeeze(r1filltbtspecy(i,1,5,:)),[.3,.3,1])
    fill(fillx,squeeze(r1filltbtspecy(i,2,1,:)),[1,.7,.7])
    fill(fillx,squeeze(r1filltbtspecy(i,2,5,:)),[1,.3,.3])
    set(gca,'yscale','log')
    legend('small l0','large l0','small l1','large l1')
    title(animalids{i})
end

% coherence small vs large
for i =1:size(r1filltbtcohery,1)
    figure
    fill(fillx,squeeze(r1filltbtcohery(i,1,1,:)),[.7,.7,1])
    hold on
    fill(fillx,squeeze(r1filltbtcohery(i,1,5,:)),[.3,.3,1])
    fill(fillx,squeeze(r1filltbtcohery(i,2,1,:)),[1,.7,.7])
    fill(fillx,squeeze(r1filltbtcohery(i,2,5,:)),[1,.3,.3])
    set(gca,'yscale','log')
    legend('small l0','large l0','small l1','large l1')
    title(animalids{i})
end

figure
plot(1,lgpow(:,1,1),'co')
hold on
plot(2,lgpow(:,1,5),'bo')
plot(3,lgpow(:,2,5),'ro')
axis([.5,3.5,0,300])
for i = 1:7
    plot([1,2],[lgpow(:,1,1),lgpow(:,1,5)],'k')
    plot([2,3],[lgpow(:,1,5),lgpow(:,2,5)],'k')
end
set(gca,'xtick',[1,2,3]);
set(gca,'xticklabel',{'small','large','large+light'})
title('low gamma power')
ylabel('psd at peak')

figure
plot(1,lglfpcoher(:,1,1),'co')
hold on
plot(2,lglfpcoher(:,1,5),'bo')
plot(3,lglfpcoher(:,2,5),'ro')
axis([.5,3.5,.4,.9])
for i = 1:7
    plot([1,2],[lglfpcoher(:,1,1),lglfpcoher(:,1,5)],'k')
    plot([2,3],[lglfpcoher(:,1,5),lglfpcoher(:,2,5)],'k')
end
set(gca,'xtick',[1,2,3]);
set(gca,'xticklabel',{'small','large','largelight'})
title('coherence between sites')
ylabel('coherence')

figure
plot(1,lgsfcoher(:,1,1),'co')
hold on
plot(2,lgsfcoher(:,1,5),'bo')
plot(3,lgsfcoher(:,2,5),'ro')
axis([.5,3.5,.2,.9])
for i = 1:7
    plot([1,2],[lgsfcoher(:,1,1),lgsfcoher(:,1,5)],'k')
    plot([2,3],[lgsfcoher(:,1,5),lgsfcoher(:,2,5)],'k')
end
set(gca,'xtick',[1,2,3]);
set(gca,'xticklabel',{'small','large','large+light'})
title('spike field coherence')
ylabel('coherence')

% runing only
figure
plot(1,r1lgpow(:,1,1),'co')
hold on
plot(2,r1lgpow(:,1,5),'bo')
plot(3,r1lgpow(:,2,5),'ro')
axis([.5,3.5,0,300])
for i = 1:7
    plot([1,2],[r1lgpow(:,1,1),r1lgpow(:,1,5)],'k')
    plot([2,3],[r1lgpow(:,1,5),r1lgpow(:,2,5)],'k')
end
set(gca,'xtick',[1,2,3]);
set(gca,'xticklabel',{'small','large','large+light'})
title('low gamma power')
ylabel('psd at peak')

figure
plot(1,r1lglfpcoher(:,1,1),'co')
hold on
plot(2,r1lglfpcoher(:,1,5),'bo')
plot(3,r1lglfpcoher(:,2,5),'ro')
axis([.5,3.5,.4,.9])
for i = 1:7
    plot([1,2],[r1lglfpcoher(:,1,1),r1lglfpcoher(:,1,5)],'k')
    plot([2,3],[r1lglfpcoher(:,1,5),r1lglfpcoher(:,2,5)],'k')
end
set(gca,'xtick',[1,2,3]);
set(gca,'xticklabel',{'small','large','large+light'})
title('coherence between sites')
ylabel('coherence')

figure
plot(1,r1lgsfcoher(:,1,1),'co')
hold on
plot(2,r1lgsfcoher(:,1,5),'bo')
plot(3,r1lgsfcoher(:,2,5),'ro')
axis([.5,3.5,.2,.9])
for i = 1:7
    plot([1,2],[r1lgsfcoher(:,1,1),r1lgsfcoher(:,1,5)],'k')
    plot([2,3],[r1lgsfcoher(:,1,5),r1lgsfcoher(:,2,5)],'k')
end
set(gca,'xtick',[1,2,3]);
set(gca,'xticklabel',{'small','large','large+light'})
title('spike field coherence')
ylabel('coherence')


for i = 1:length(blocks)
    figure('name',animalids{i});
    subplot(2,3,1)
    semilogy(chf,squeeze(lfpspect_c(i,3,1,1,:)),'b')
    hold on    
    semilogy(chf,squeeze(lfpspect_c(i,3,1,2,:)),'c')
    semilogy(chf,squeeze(lfpspect_c(i,3,1,3,:)),'g')
    semilogy(chf,squeeze(lfpspect_c(i,3,1,4,:)),'r')
    semilogy(chf,squeeze(lfpspect_c(i,3,1,5,:)),'m')
    ax = axis;
    axis([0,120,1,ax(4)])
    legend('smallest','-','-','-','largest')
    xlabel('frequency')
    ylabel('psd')
    
    subplot(2,3,2)
    plot(chf,squeeze(lfpcoherence(i,3,1,1,:)),'b');
    hold on
    plot(chf,squeeze(lfpcoherence(i,3,1,2,:)),'c');
    plot(chf,squeeze(lfpcoherence(i,3,1,3,:)),'g');
    plot(chf,squeeze(lfpcoherence(i,3,1,4,:)),'r');
    plot(chf,squeeze(lfpcoherence(i,3,1,5,:)),'m');
    axis([0,120,0,1]);
    xlabel('frequency')
    ylabel('site coherence')
    
    
    subplot(2,3,3)
    plot(chf,squeeze(sfcoher_c(i,3,1,1,:)),'b');
    hold on
    plot(chf,squeeze(sfcoher_c(i,3,1,2,:)),'c');
    plot(chf,squeeze(sfcoher_c(i,3,1,3,:)),'g');
    plot(chf,squeeze(sfcoher_c(i,3,1,4,:)),'r');
    plot(chf,squeeze(sfcoher_c(i,3,1,5,:)),'m');
    axis([0,120,0,1]);
    xlabel('frequency')
    ylabel('spike field coherence')
    
    subplot(2,3,4)
    semilogy(chf,squeeze(lfpspect_c(i,3,1,1,:)),'b')
    hold on    
    semilogy(chf,squeeze(lfpspect_c(i,3,1,5,:)),'c')
    semilogy(chf,squeeze(lfpspect_c(i,3,2,1,:)),'r')
    semilogy(chf,squeeze(lfpspect_c(i,3,2,5,:)),'m')
    ax = axis;
    axis([0,120,1,ax(4)])
    legend('small l0','large l0','small l1','large l1')
    xlabel('frequency')
    ylabel('psd')
    
    subplot(2,3,5)
    plot(chf,squeeze(lfpcoherence(i,3,1,1,:)),'b');
    hold on
    plot(chf,squeeze(lfpcoherence(i,3,1,5,:)),'c');
    plot(chf,squeeze(lfpcoherence(i,3,2,1,:)),'r');
    plot(chf,squeeze(lfpcoherence(i,3,2,5,:)),'m');
    axis([0,120,0,1]);
    xlabel('frequency')
    ylabel('site coherence')
    
    subplot(2,3,6)
    plot(chf,squeeze(sfcoher_c(i,3,1,1,:)),'b');
    hold on
    plot(chf,squeeze(sfcoher_c(i,3,1,5,:)),'c');
    plot(chf,squeeze(sfcoher_c(i,3,2,1,:)),'r');
    plot(chf,squeeze(sfcoher_c(i,3,2,5,:)),'m');
    axis([0,120,0,1]);
    xlabel('frequency')
    ylabel('spike field coherence')
    set(gcf,'OuterPosition',[50,300,1200,700])
end
    
% same running only
for i = 1:length(blocks)
    figure('name',animalids{i});
    subplot(2,3,1)
    semilogy(chf,squeeze(r1_lfpspect_c(i,3,1,1,:)),'b')
    hold on    
    semilogy(chf,squeeze(r1_lfpspect_c(i,3,1,2,:)),'c')
    semilogy(chf,squeeze(r1_lfpspect_c(i,3,1,3,:)),'g')
    semilogy(chf,squeeze(r1_lfpspect_c(i,3,1,4,:)),'r')
    semilogy(chf,squeeze(r1_lfpspect_c(i,3,1,5,:)),'m')
    ax = axis;
    axis([0,120,1,ax(4)+0.0001])
    legend('smallest','-','-','-','largest')
    xlabel('frequency')
    ylabel('psd')
    
    subplot(2,3,2)
    plot(chf,squeeze(r1_lfpcoherence(i,3,1,1,:)),'b');
    hold on
    plot(chf,squeeze(r1_lfpcoherence(i,3,1,2,:)),'c');
    plot(chf,squeeze(r1_lfpcoherence(i,3,1,3,:)),'g');
    plot(chf,squeeze(r1_lfpcoherence(i,3,1,4,:)),'r');
    plot(chf,squeeze(r1_lfpcoherence(i,3,1,5,:)),'m');
    axis([0,120,0,1]);
    xlabel('frequency')
    ylabel('site coherence')
    
    
    subplot(2,3,3)
    plot(chf,squeeze(r1_sfcoher_c(i,3,1,1,:)),'b');
    hold on
    plot(chf,squeeze(r1_sfcoher_c(i,3,1,2,:)),'c');
    plot(chf,squeeze(r1_sfcoher_c(i,3,1,3,:)),'g');
    plot(chf,squeeze(r1_sfcoher_c(i,3,1,4,:)),'r');
    plot(chf,squeeze(r1_sfcoher_c(i,3,1,5,:)),'m');
    axis([0,120,0,1]);
    xlabel('frequency')
    ylabel('spike field coherence')
    
    subplot(2,3,4)
    semilogy(chf,squeeze(r1_lfpspect_c(i,3,1,1,:)),'b')
    hold on    
    semilogy(chf,squeeze(r1_lfpspect_c(i,3,1,5,:)),'c')
    semilogy(chf,squeeze(r1_lfpspect_c(i,3,2,1,:)),'r')
    semilogy(chf,squeeze(r1_lfpspect_c(i,3,2,5,:)),'m')
    ax = axis;
    axis([0,120,1,ax(4)+0.0001])
    legend('small l0','large l0','small l1','large l1')
    xlabel('frequency')
    ylabel('psd')
    
    subplot(2,3,5)
    plot(chf,squeeze(r1_lfpcoherence(i,3,1,1,:)),'b');
    hold on
    plot(chf,squeeze(r1_lfpcoherence(i,3,1,5,:)),'c');
    plot(chf,squeeze(r1_lfpcoherence(i,3,2,1,:)),'r');
    plot(chf,squeeze(r1_lfpcoherence(i,3,2,5,:)),'m');
    axis([0,120,0,1]);
    xlabel('frequency')
    ylabel('site coherence')
    
    subplot(2,3,6)
    plot(chf,squeeze(r1_sfcoher_c(i,3,1,1,:)),'b');
    hold on
    plot(chf,squeeze(r1_sfcoher_c(i,3,1,5,:)),'c');
    plot(chf,squeeze(r1_sfcoher_c(i,3,2,1,:)),'r');
    plot(chf,squeeze(r1_sfcoher_c(i,3,2,5,:)),'m');
    axis([0,120,0,1]);
    xlabel('frequency')
    ylabel('spike field coherence')
    set(gcf,'OuterPosition',[50,100,1200,700])
end
    

for i = 1:length(blocks)
    
    depthax = centdepth(i):-spacing(i):centdepth(i)-(length(centelectrodes(i,1):centelectrodes(i,2))-1)*spacing(i);
    
    figure
    subplot(3,2,1)
    imagesc(chf,depthax,log(squeeze(lfpspect_c(i,:,1,1,:))))
    axis([1,120,depthax(end)-10,depthax(1)+10])
    caxis([-2,9])
    title('power small')
    
    subplot(3,2,2)
    imagesc(chf,depthax,log(squeeze(lfpspect_c(i,:,1,5,:))))
    axis([1,120,depthax(end)-10,depthax(1)+10])
    caxis([-2,9])
    title('power large')
    
    subplot(3,2,3)
    imagesc(chf,depthax,log(squeeze(lfpspect_c(i,:,1,5,:)))-log(squeeze(lfpspect_c(i,:,1,1,:))));
    axis([1,120,depthax(end)-10,depthax(1)+10])
    caxis([-2,2])
    title('large - small') 
    
    subplot(3,2,5)
    imagesc(chf,depthax,log(squeeze(lfpspect_c(i,:,2,1,:)))-log(squeeze(lfpspect_c(i,:,1,1,:))));
    axis([1,120,depthax(end)-10,depthax(1)+10])
    caxis([-2,2])
    title('small L1-L0')
    
    subplot(3,2,6)
    imagesc(chf,depthax,log(squeeze(lfpspect_c(i,:,2,5,:)))-log(squeeze(lfpspect_c(i,:,1,5,:))));
    axis([1,120,depthax(end)-10,depthax(1)+10])
    caxis([-2,2])
    title('large L1-L0')
    
    
    
    figure
    subplot(3,2,1)
    imagesc(chf,depthax,squeeze(lfpcoherence(i,:,1,1,:)))
    axis([1,120,depthax(end)-10,depthax(1)+10])
    caxis([0,1])
    title('coherence smally')
    
    subplot(3,2,2)
    imagesc(chf,depthax,squeeze(lfpcoherence(i,:,1,5,:)))
    axis([1,120,depthax(end)-10,depthax(1)+10])
    caxis([0,1])
    title('coherence large') 
    
    subplot(3,2,3)
    imagesc(chf,depthax,squeeze(lfpcoherence(i,:,1,5,:))-squeeze(lfpcoherence(i,:,1,1,:)))
    axis([1,120,depthax(end)-10,depthax(1)+10])
    caxis([-.5,.5])
    title('large - small')    
      
    subplot(3,2,5)
    imagesc(chf,depthax,squeeze(lfpcoherence(i,:,2,1,:))-squeeze(lfpcoherence(i,:,1,1,:)))
    axis([1,120,depthax(end)-10,depthax(1)+10])
    caxis([-.3,.3])
    title('small L1-L0')       
    
    subplot(3,2,6)
    imagesc(chf,depthax,squeeze(lfpcoherence(i,:,2,5,:))-squeeze(lfpcoherence(i,:,1,5,:)))
    axis([1,120,depthax(end)-10,depthax(1)+10])
    caxis([-.3,.3])
    title('large L1-L0')    
    
    
    figure
    subplot(3,2,1)
    imagesc(chf,depthax,squeeze(sfcoher_c(i,:,1,1,:)))
    axis([1,120,depthax(end)-10,depthax(1)+10])
    caxis([0,.6])
    title('SFC small')
    
    subplot(3,2,2)
    imagesc(chf,depthax,squeeze(sfcoher_c(i,:,1,5,:)))
    axis([1,120,depthax(end)-10,depthax(1)+10])
    caxis([0,.6])
    title('SFC large') 
    
    subplot(3,2,3)
    imagesc(chf,depthax,squeeze(sfcoher_c(i,:,1,5,:))-squeeze(sfcoher_c(i,:,1,1,:)))
    axis([1,120,depthax(end)-10,depthax(1)+10])
    caxis([-.2,.2])
    title('SFC large - small')
    
    subplot(3,2,5)    
    imagesc(chf,depthax,squeeze(sfcoher_c(i,:,2,1,:))-squeeze(sfcoher_c(i,:,1,1,:)))
    axis([1,120,depthax(end)-10,depthax(1)+10])
    caxis([-.2,.2])
    title('small L1-L0')
    
    subplot(3,2,6)    
    imagesc(chf,depthax,squeeze(sfcoher_c(i,:,2,5,:))-squeeze(sfcoher_c(i,:,1,5,:)))
    axis([1,120,depthax(end)-10,depthax(1)+10])
    caxis([-.2,.2])
    title('large L1-L0')

end


 disp('');