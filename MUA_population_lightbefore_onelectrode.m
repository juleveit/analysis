function MUA_population_lightbefore_onelectrode
% 
% % PV Halo or eArch pop
% animalids  = {'151104', '151211', '160204', '160205', '160217', '160320', '160324', '160328'};
% lfblocks    = [6,        4,        13,       12,       9,        5,        3,        3];
% lablocks    = [5,        3,        14,       4,        11,       4,        4,        4];
% animal     =  [1,        3,        4,        5,        6,        7,        8,        9];
% electrodes = [[17,32];  [1,16];   [1,16];   [17,32];  [1,16];   [1,16];   [1,16];   [1,16]];
% bldepth     = [500,      500,      500,      500,      500,      550,      500,      500];
% penangle    = [25,       25,       25,       25,       25,       25,       25,       25];
% spacing     = [25,       25,       25,       25,       25,       25,       25,       25];
% popfile     = 'C:\Users\Julia\work\data\populations\PV_HaloeArch\lightbefore\MUA_population_singleel.mat';

% SOM Halo
animalids =    {'151022', '151023', '151027', '151109', '151110', '160419', '160422'};
lfblocks    =   [8,        13,       10,       10,       12,       2,        2];
lablocks    =   [6,        12,       3,        11,       13,       3,        3];
animal    =     [1,        2,        3,        4,        5,        6,        8];
electrodes =   [[1,16];   [1,16];   [17,32];  [17,32];  [1,16];   [1,16];   [1,16]];
bldepth    =    [400,      400,      500,      500,      500,      500,      500];
penangle =      [25,       25,       25,       25,       25,       25,       25];
spacing =       [25,       25,       25,       25,       25,       25,       25];
popfile = 'C:\Users\Julia\work\data\populations\SOM_Halo_later\lightbefore\MUA_population_singleel.mat';

% % single thing
% animalids =    {'160217'};
% lfblocks    =   [9];
% lablocks    =   [11];
% animal    =     [1];
% electrodes =[[1,16]];
% bldepth    =  [500];
% penangle =      [25];
% spacing =       [25];
% popfile = 'C:\Users\Julia\work\data\populations\temp2.mat';

depth = bldepth.*(cosd(penangle)*cosd(22));
spacing = spacing.*(cosd(penangle)*cosd(22));
evaldepth = 330;
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

% determine analysis window for light before blocks
% anawin = 'same'; % same window as light after = 700ms after vis onset to exclude light artefact
anawin = 'immediate'; % window right after onset of visual stimulus 

if ~exist(popfile) || recalculate_pop

    cll = 1;
    for blck = 1:length(lfblocks)
        
        basepath = strcat('C:\Users\Julia\work\data\', animalids{blck}, '\');
        lffile = strcat(basepath, 'muaresult_', int2str(lfblocks(blck)), '_', int2str(electrodes(blck,1)), '_', int2str(electrodes(blck,2)), '.mat');
        lafile = strcat(basepath, 'muaresult_', int2str(lablocks(blck)), '_', int2str(electrodes(blck,1)), '_', int2str(electrodes(blck,2)), '.mat');
        altlffile = strcat(basepath, 'muaresult_', int2str(lfblocks(blck)), '_', int2str(electrodes(blck,1)), '-', int2str(electrodes(blck,2)), '.mat');
        altlafile = strcat(basepath, 'muaresult_', int2str(lablocks(blck)), '_', int2str(electrodes(blck,1)), '-', int2str(electrodes(blck,2)), '.mat');
        
        clear result;
        if ~exist(lffile) || recalculate_muafile
            if exist(altlffile)
                load(altlfffile);
                save(lffile,'result');
            else
                result = MUAdataprepare(basepath,animalids{blck},lfblocks(blck),electrodes(blck,1):electrodes(blck,2));
                save(lffile,'result');
            end
        else
            load(lffile);            
        end              
        lfresult = result;        
        clear result;
        if ~exist(lafile) || recalculate_muafile
            if exist(altlafile)
                load(altlafile);
                result = centresult;
                save(lafile,'result');
            else
                result = MUAdataprepare(basepath,animalids{blck},lablocks(blck),electrodes(blck,1):electrodes(blck,2));
                save(lafile,'result');
            end
        else
            load(lafile);
        end        
        laresult = result;  
        clear result;
        
        prestim = 700;
        poststim = 700;
        respwin = 501:1500; % after stimulus onset
        offsetwin = 1501:2500;
        respwin = respwin+prestim;
                
        ch = di(blck);
            
        disp(['Block ' int2str(blck) '/' int2str(length(lfblocks)) '   channel ' int2str(ch) '/' int2str(length(electrodes(blck,1):electrodes(blck,2)))]);
        
        trialdur = lfresult.stimduration*1000;
        lfmsstamps = lfresult.msstamps;
        lamsstamps = laresult.msstamps;
        if length(lfmsstamps)~=length(lfresult.light)
%             lfmsstamps([161,287]) = [];
%             lfmsstamps([33]) = [];
%             lfresult.msstamps = lfmsstamps;
%             result = lfresult;
%             save(lffile,'result');
            pause;
        end
        if length(lamsstamps)~=length(laresult.light)
%             lamsstamps([161,287]) = [];
%             laresult.msstamps = lamsstamps;
%             result = laresult;
%             save(lafile,'result');
            pause;
        end
        
        for i = 1:length(lfmsstamps)
            lfspeed(i,:) = lfresult.runspeed(lfmsstamps(i) - prestim+1:lfmsstamps(i) + trialdur + poststim);
        end
        for i = 1:length(lamsstamps)
            laspeed(i,:) = laresult.runspeed(lamsstamps(i) - prestim+1:lamsstamps(i) + trialdur + poststim);
        end
        
        % figure out sufficiently high and nonvariable runspeed trials
        lfmeanspeed = mean(lfspeed(:,respwin),2);
        lfstdspeed = std(lfspeed(:,respwin),1,2);
        lfnotstill = find(lfmeanspeed>1);
        lfokspeed = find(lfmeanspeed>( mean(lfmeanspeed(lfnotstill))-(1.5*std(lfmeanspeed(lfnotstill))) ) & lfmeanspeed > 1);
        lfokvar = find(lfstdspeed<( mean(lfstdspeed(lfnotstill))+(1.5*std(lfstdspeed(lfnotstill)))) & lfstdspeed>.5);
        lfoktrials = intersect(lfokspeed,lfokvar);
        lfnonoktrials = 1:size(lfspeed,1); lfnonoktrials(lfoktrials) = [];
        lfstilltrials = 1:size(lfspeed,1); lfstilltrials(lfnotstill) = [];
        
        lameanspeed = mean(laspeed(:,respwin),2);
        lastdspeed = std(laspeed(:,respwin),1,2);
        lanotstill = find(lameanspeed>1);
        laokspeed = find(lameanspeed>( mean(lameanspeed(lanotstill))-(1.5*std(lameanspeed(lanotstill))) ) & lameanspeed > 1);
        laokvar = find(lastdspeed<( mean(lastdspeed(lanotstill))+(1.5*std(lastdspeed(lanotstill)))) & lastdspeed>.5);
        laoktrials = intersect(laokspeed,laokvar);
        lanonoktrials = 1:size(laspeed,1); lanonoktrials(laoktrials) = [];
        lastilltrials = 1:size(laspeed,1); lastilltrials(lanotstill) = [];
        
        lfsize = unique(lfresult.gratingInfo.size);
        lasizes = unique(laresult.gratingInfo.size);
        [m,lasizeind] = min(abs(lasizes-lfsize));
        
        % find low and high gamma peaks in spectra
        beta = [24,40];        
        lalarge = find(laresult.gratingInfo.size == lasizes(lasizeind) & laresult.light == 0);
        lflarge = find(lfresult.light == 0);
        lflr1 = intersect(lflarge,lfoktrials);
        lalr1 = intersect(lalarge,laoktrials);
        for i = 1:length(lflarge)
            [lfpl(i,:),f] = pmtm(lfresult.lfp(ch,lfresult.msstamps(lflarge(i))+300:lfresult.msstamps(lflarge(i))+1300),3,[],1000);
            [lapl(i,:),f] = pmtm(laresult.lfp(ch,laresult.msstamps(lalarge(i))+300:laresult.msstamps(lalarge(i))+1300),3,[],1000);
        end
        lfplr1 = nan(1,size(lfpl,2)); laplr1 = nan(1,size(lapl,2));
        if length(lflr1>=5)
            for i = 1:length(lflr1)
                [lfplr1(i,:),f] = mtspectrumc(lfresult.lfp(ch,lfresult.msstamps(lflr1(i))+700:lfresult.msstamps(lflr1(i))+1500),params);
                %                     [plr1(i,:),f] = pmtm(result.lfp(ch,result.msstamps(lr1(i))+300:result.msstamps(lr1(i))+1300),3,[],1000);
            end
        end
        if length(lalr1>=5)
            for i = 1:length(lalr1)
                [laplr1(i,:),f] = mtspectrumc(laresult.lfp(ch,laresult.msstamps(lalr1(i))+700:laresult.msstamps(lalr1(i))+1500),params);
                %                     [psr1(i,:),f] = pmtm(result.lfp(ch,result.msstamps(sr1(i))+300:result.msstamps(sr1(i))+1300),3,[],1000);
            end
        end
        b1 = find(f>beta(1),1); b2 = find(f>beta(2),1);
        if ~isnan(lfplr1(1)) % take the peak from the running spectra if there is more than a couple of tunning trials
            lfbsig = nanmean(lfplr1(:,b1:b2));
        else
            lfbsig = nanmean(lfpl(:,b1:b2));
        end
        if ~isnan(laplr1(1))
            labsig = nanmean(laplr1(:,b1:b2));
        else
            labsig = nanmean(lapl(:,b1:b2));
        end
        if isempty(find(diff(lfbsig)>0)) % there is no clear beta peak
            bpi = round((b1+b2)/2);
            lfgpeak(blck) = 0;
        else
            peaks = find(diff(lfbsig)>0)+1;
            pvs = lfbsig(peaks);
            bpi = peaks(pvs == max(pvs));
            bpi = bpi+b1-1;
            lfgpeak(blck) = 1;
        end
        if isempty(find(diff(labsig)>0)) % there is no clear beta peak
            gpi = round((b1+b2)/2);
            lagpeak(blck) = 0;
        else
            peaks = find(diff(labsig)>0)+1;
            pvs = labsig(peaks);
            gpi = peaks(pvs == max(pvs));
            gpi = gpi+b1-1;
            lagpeak(blck) = 1;
        end
        lflgi(blck) = bpi; lalgi(blck) = gpi;
        
        sr = 1000;
        lflfp = lfresult.lfp(ch,:); lalfp = laresult.lfp(ch,:);
        lfgamma = eegfilt(lflfp,sr,f(lflgi(blck))-5,f(lflgi(blck))+5);  
        h = hilbert(lfgamma); lfgpow = abs(h); 
        lagamma = eegfilt(lalfp,sr,f(lalgi(blck))-5,f(lalgi(blck))+5);  
        h = hilbert(lagamma); lagpow = abs(h); 
        
        lfmsStimes = round(lfresult.msStimes{ch});
        if ~isempty(lfmsStimes) & lfmsStimes(1) == 0, lfmsStimes(1) = 1; end
        lamsStimes = round(laresult.msStimes{ch});
        if ~isempty(lamsStimes) & lamsStimes(1) == 0, lamsStimes(1) = 1; end
        
        lfchan = zeros(1,size(lfresult.lfp,2));
        lfchan(lfmsStimes) = 1;        
        lachan = zeros(1,size(laresult.lfp,2));
        lachan(lamsStimes) = 1;        
        
        for i = 1:length(lfmsstamps)
            lfresp(i,:) = lfchan(lfmsstamps(i) - prestim+1:lfmsstamps(i) + trialdur + poststim);
            lflfpresp(i,:) = lfresult.lfp(ch,lfmsstamps(i) - prestim+1:lfmsstamps(i) + trialdur + poststim);
            if strcmp(anawin,'same')
                lflfpspecttbt(i,:) = mtspectrumc(squeeze(lflfpresp(i,prestim+500+200+1:prestim+500+1000-10))',params); % same window as light after
            else
                lflfpspecttbt(i,:) = mtspectrumc(squeeze(lflfpresp(i,prestim+1:prestim+800-10))',params); % first window after vis on
            end
            lfgpowresp(i,:) = lfgpow(lfmsstamps(i) - prestim+1:lfmsstamps(i) + trialdur + poststim);
        end
        for i = 1:length(lamsstamps)
            laresp(i,:) = lachan(lamsstamps(i) - prestim+1:lamsstamps(i) + trialdur + poststim);
            lalfpresp(i,:) = laresult.lfp(ch,lamsstamps(i) - prestim+1:lamsstamps(i) + trialdur + poststim);
            lalfpspecttbt(i,:) = mtspectrumc(squeeze(lalfpresp(i,prestim+500+200+1:prestim+500+1000-10))',params);
            lagpowresp(i,:) = lagpow(lamsstamps(i) - prestim+1:lamsstamps(i) + trialdur + poststim);
        end
        
        lffrs = sum(lfresp(:,respwin),2)./(length(respwin)/1000);
        lafrs = sum(laresp(:,respwin),2)./(length(respwin)/1000);
        
        xfill = [-699:2300,fliplr(-699:2300)];
        lfl0s = find(lfresult.light == 0); lfl1s = find(lfresult.light == 1);
        lal0s = find(laresult.light == 0); lal1s = find(laresult.light == 1);
        lfyl0(blck,:) = [mean(lfgpowresp(lfl0s,1:3000))+std(lfgpowresp(lfl0s,1:3000))./sqrt(length(lfl0s)),fliplr(mean(lfgpowresp(lfl0s,1:3000))-std(lfgpowresp(lfl0s,1:3000))./sqrt(length(lfl0s)))];
        lfyl1(blck,:) = [mean(lfgpowresp(lfl1s,1:3000))+std(lfgpowresp(lfl1s,1:3000))./sqrt(length(lfl1s)),fliplr(mean(lfgpowresp(lfl1s,1:3000))-std(lfgpowresp(lfl1s,1:3000))./sqrt(length(lfl1s)))];
        layl0(blck,:) = [mean(lagpowresp(lal0s,1:3000))+std(lagpowresp(lal0s,1:3000))./sqrt(length(lal0s)),fliplr(mean(lagpowresp(lal0s,1:3000))-std(lagpowresp(lal0s,1:3000))./sqrt(length(lal0s)))];
        layl1(blck,:) = [mean(lagpowresp(lal1s,1:3000))+std(lagpowresp(lal1s,1:3000))./sqrt(length(lal1s)),fliplr(mean(lagpowresp(lal1s,1:3000))-std(lagpowresp(lal1s,1:3000))./sqrt(length(lal1s)))];
        
        for l = 1:length(unique(lfresult.light))
            
            lfthisinds = find(lfresult.light == l-1);
            lathisinds = find(laresult.light == l-1 & laresult.gratingInfo.size == lasizes(lasizeind));
            lfthisruninds = intersect(lfthisinds,lfoktrials); lfthisstillinds = intersect(lfthisinds,lfstilltrials);
            lathisruninds = intersect(lathisinds,laoktrials); lathisstillinds = intersect(lathisinds,lastilltrials);
            
            lfcondresp(blck,l,:) = nanmean(lfresp(lfthisinds,:),1);
            lfcondfr(blck,l) = nanmean(lffrs(lfthisinds));
            lfcondlfpresp(blck,l,:) = nanmean(lflfpresp(lfthisinds,:),1);
            lacondresp(blck,l,:) = nanmean(laresp(lathisinds,:),1);
            lacondfr(blck,l) = nanmean(lafrs(lathisinds));
            lacondlfpresp(blck,l,:) = nanmean(lalfpresp(lathisinds,:),1);
            
            if strcmp(anawin,'same')
                [lflfpspect(blck,l,:),chf,lflfpspecterr(blck,l,:,:)] = mtspectrumc(squeeze(lflfpresp(lfthisinds,prestim+500+200+1:prestim+500+1000-10))',params);
            else
                [lflfpspect(blck,l,:),chf,lflfpspecterr(blck,l,:,:)] = mtspectrumc(squeeze(lflfpresp(lfthisinds,prestim+1:prestim+800-10))',params);
            end
            lftbtlfpspect(blck,l,:) = nanmean(lflfpspecttbt(lfthisinds,:),1);
            lftbtlfpspecterr(blck,l,:) = nanstd(lflfpspecttbt(lfthisinds,:),1,1)./sqrt(length(lfthisinds));
            [lalfpspect(blck,l,:),chf,lalfpspecterr(blck,l,:,:)] = mtspectrumc(squeeze(lalfpresp(lathisinds,prestim+500+200+1:prestim+500+1000-10))',params);
            latbtlfpspect(blck,l,:) = nanmean(lalfpspecttbt(lathisinds,:),1);
            latbtlfpspecterr(blck,l,:) = nanstd(lalfpspecttbt(lathisinds,:),1,1)./sqrt(length(lathisinds));
            
            if ~isempty(lfthisruninds)
                if strcmp(anawin,'same')
                    [lfr1_lfpspect(blck,l,:),chf,lfr1_lfpspecterr(blck,l,:,:)] = mtspectrumc(squeeze(lflfpresp(lfthisruninds,prestim+500+200+1:prestim+500+1000-10))',params);
                else
                    [lfr1_lfpspect(blck,l,:),chf,lfr1_lfpspecterr(blck,l,:,:)] = mtspectrumc(squeeze(lflfpresp(lfthisruninds,prestim+1:prestim+800-10))',params);
                end
                lfr1_tbtlfpspect(blck,l,:) = nanmean(lflfpspecttbt(lfthisruninds,:),1);
                lfr1_tbtlfpspecterr(blck,l,:) = nanstd(lflfpspecttbt(lfthisruninds,:),1,1)./sqrt(length(lfthisruninds));
                lfr1_ntrials(blck,l) = length(lfthisruninds);
            else
                lfr1_lfpspect(blck,l,:) = nan(1,length(lflfpspect(blck,l,:)));
                lfr1_lfpspecterr(blck,l,:,:) = nan(2,length(lflfpspect(blck,l,:)));
                lfr1_tbtlfpspect(blck,l,:) = nan(1,length(lflfpspect(blck,l,:)));
                lfr1_tbtlfpspecterr(blck,l,:) = nan(1,length(lflfpspect(blck,l,:)));
                lfr1_ntrials(blck,l) = 0;
            end
            if ~isempty(lathisruninds)                
                [lar1_lfpspect(blck,l,:),chf,lar1_lfpspecterr(blck,l,:,:)] = mtspectrumc(squeeze(lalfpresp(lathisruninds,prestim+500+200+1:prestim+500+1000-10))',params);
                lar1_tbtlfpspect(blck,l,:) = nanmean(lalfpspecttbt(lathisruninds,:),1);
                lar1_tbtlfpspecterr(blck,l,:) = nanstd(lalfpspecttbt(lathisruninds,:),1,1)./sqrt(length(lathisruninds));
                lar1_ntrials(blck,l) = length(lathisruninds);
            else
                lar1_lfpspect(blck,l,:) = nan(1,length(lalfpspect(blck,l,:)));
                lar1_lfpspecterr(blck,l,:,:) = nan(2,length(lalfpspect(blck,l,:)));
                lar1_tbtlfpspect(blck,l,:) = nan(1,length(lalfpspect(blck,l,:)));
                lar1_tbtlfpspecterr(blck,l,:) = nan(1,length(lalfpspect(blck,l,:)));
                lar1_ntrials(blck,l) = 0;
            end                
            
            if ~isempty(lfthisstillinds)
                if strcmp(anawin,'same')
                    [lfr0_lfpspect(blck,l,:),chf,lfr0_lfpspecterr(blck,l,:,:)] = mtspectrumc(squeeze(lflfpresp(lfthisstillinds,prestim+500+200+1:prestim+500+1000-10))',params);
                else                    
                    [lfr0_lfpspect(blck,l,:),chf,lfr0_lfpspecterr(blck,l,:,:)] = mtspectrumc(squeeze(lflfpresp(lfthisstillinds,prestim+1:prestim+800-10))',params);
                end
                lfr0_tbtlfpspect(blck,l,:) = nanmean(lflfpspecttbt(lfthisstillinds,:),1);
                lfr0_tbtlfpspecterr(blck,l,:) = nanstd(lflfpspecttbt(lfthisstillinds,:),1,1)./sqrt(length(lfthisstillinds));
                lfr0_ntrials(blck,l) = length(lfthisstillinds);
            else
                lfr0_lfpspect(blck,l,:) = nan(1,length(lflfpspect(blck,l,:)));
                lfr0_lfpspecterr(blck,l,:,:) = nan(2,length(lflfpspect(blck,l,:)));
                lfr0_tbtlfpspect(blck,l,:) = nan(1,length(lflfpspect(blck,l,:)));
                lfr0_tbtlfpspecterr(blck,l,:) = nan(1,length(lflfpspect(blck,l,:)));
                lfr0_ntrials(blck,l) = 0;
            end
            if ~isempty(lathisstillinds)                
                [lar0_lfpspect(blck,l,:),chf,lar0_lfpspecterr(blck,l,:,:)] = mtspectrumc(squeeze(lalfpresp(lathisstillinds,prestim+500+200+1:prestim+500+1000-10))',params);
                lar0_tbtlfpspect(blck,l,:) = nanmean(lalfpspecttbt(lathisstillinds,:),1);
                lar0_tbtlfpspecterr(blck,l,:) = nanstd(lalfpspecttbt(lathisstillinds,:),1,1)./sqrt(length(lathisstillinds));
                lar0_ntrials(blck,l) = length(lathisstillinds);
            else
                lar0_lfpspect(blck,l,:) = nan(1,length(lalfpspect(blck,l,:)));
                lar0_lfpspecterr(blck,l,:,:) = nan(2,length(lalfpspect(blck,l,:)));
                lar0_tbtlfpspect(blck,l,:) = nan(1,length(lalfpspect(blck,l,:)));
                lar0_tbtlfpspecterr(blck,l,:) = nan(1,length(lalfpspect(blck,l,:)));
                lar0_ntrials(blck,l) = 0;
            end                
        end
    end
    save(popfile, '-v7.3');
else
    load(popfile);
end

% vi = ones(1,size(lfr1_lfpspect,1)); 
vi = 1:length(depth);
for i = 1:length(vi)
    lflgpow(i,:) = squeeze(lflfpspect(vi(i),:,lflgi(vi(i))));
    lalgpow(i,:) = squeeze(lalfpspect(vi(i),:,lalgi(vi(i))));
    lfr1lgpow(i,:) = squeeze(lfr1_lfpspect(vi(i),:,lflgi(vi(i))));
    lfr0lgpow(i,:) = squeeze(lfr0_lfpspect(vi(i),:,lflgi(vi(i))));
    lar1lgpow(i,:) = squeeze(lar1_lfpspect(vi(i),:,lalgi(vi(i))));
    lar0lgpow(i,:) = squeeze(lar0_lfpspect(vi(i),:,lalgi(vi(i))));
    rellfpow(i,:) = squeeze(lflfpspect(vi(i),:,lflgi(vi(i))))./squeeze(nanmean(lflfpspect(vi(i),:,12:103),3));
    rellapow(i,:) = squeeze(lalfpspect(vi(i),:,lalgi(vi(i))))./squeeze(nanmean(lalfpspect(vi(i),:,12:103),3));
    rellfr1pow(i,:) = squeeze(lfr1_lfpspect(vi(i),:,lflgi(vi(i))))./squeeze(nanmean(lfr1_lfpspect(vi(i),:,12:103),3));
    rellar1pow(i,:) = squeeze(lar1_lfpspect(vi(i),:,lflgi(vi(i))))./squeeze(nanmean(lar1_lfpspect(vi(i),:,12:103),3));
    
    for l = 1:2
        lffillspecy(i,l,:) = [squeeze(lflfpspecterr(vi(i),l,1,1:104))',fliplr(squeeze(lflfpspecterr(vi(i),l,2,1:104))')];
        lffilltbtspecy(i,l,:) = [(squeeze(lftbtlfpspect(i,l,1:104))+squeeze(lftbtlfpspecterr(i,l,1:104)))',fliplr((squeeze(lftbtlfpspect(i,l,1:104))-squeeze(lftbtlfpspecterr(i,l,1:104)))')]';
        lfr1fillspecy(i,l,:) = [squeeze(lfr1_lfpspecterr(vi(i),l,1,1:104))',fliplr(squeeze(lfr1_lfpspecterr(vi(i),l,2,1:104))')];
        lfr1filltbtspecy(i,l,:) = [(squeeze(lfr1_tbtlfpspect(i,l,1:104))+squeeze(lfr1_tbtlfpspecterr(i,l,1:104)))',fliplr((squeeze(lfr1_tbtlfpspect(i,l,1:104))-squeeze(lfr1_tbtlfpspecterr(i,l,1:104)))')]';
        lfr1fillspecy(i,l,:) = [squeeze(lfr1_lfpspecterr(vi(i),l,1,1:104))',fliplr(squeeze(lfr1_lfpspecterr(vi(i),l,2,1:104))')];
        lfr0filltbtspecy(i,l,:) = [(squeeze(lfr0_tbtlfpspect(i,l,1:104))+squeeze(lfr0_tbtlfpspecterr(i,l,1:104)))',fliplr((squeeze(lfr0_tbtlfpspect(i,l,1:104))-squeeze(lfr0_tbtlfpspecterr(i,l,1:104)))')]';
        lafillspecy(i,l,:) = [squeeze(lalfpspecterr(vi(i),l,1,1:104))',fliplr(squeeze(lalfpspecterr(vi(i),l,2,1:104))')];
        lafilltbtspecy(i,l,:) = [(squeeze(latbtlfpspect(i,l,1:104))+squeeze(latbtlfpspecterr(i,l,1:104)))',fliplr((squeeze(latbtlfpspect(i,l,1:104))-squeeze(latbtlfpspecterr(i,l,1:104)))')]';
        lar1fillspecy(i,l,:) = [squeeze(lar1_lfpspecterr(vi(i),l,1,1:104))',fliplr(squeeze(lar1_lfpspecterr(vi(i),l,2,1:104))')];
        lar1filltbtspecy(i,l,:) = [(squeeze(lar1_tbtlfpspect(i,l,1:104))+squeeze(lar1_tbtlfpspecterr(i,l,1:104)))',fliplr((squeeze(lar1_tbtlfpspect(i,l,1:104))-squeeze(lar1_tbtlfpspecterr(i,l,1:104)))')]';
        lar1fillspecy(i,l,:) = [squeeze(lar1_lfpspecterr(vi(i),l,1,1:104))',fliplr(squeeze(lar1_lfpspecterr(vi(i),l,2,1:104))')];
        lar0filltbtspecy(i,l,:) = [(squeeze(lar0_tbtlfpspect(i,l,1:104))+squeeze(lar0_tbtlfpspecterr(i,l,1:104)))',fliplr((squeeze(lar0_tbtlfpspect(i,l,1:104))-squeeze(lar0_tbtlfpspecterr(i,l,1:104)))')]';
    end    
end
fillx = [chf(1:104),fliplr(chf(1:104))];


for i =1:size(lffilltbtspecy,1)
    figure
    
    subplot(2,2,1)
    fill(fillx,squeeze(lffilltbtspecy(i,1,:)),[.3,.3,1])
    hold on
    fill(fillx,squeeze(lffilltbtspecy(i,2,:)),[1,.3,.3])
    set(gca,'yscale','log')
    legend('light first l0','light first l1')    
    title(animalids{i})
    
    subplot(2,2,2)
    fill(fillx,squeeze(lafilltbtspecy(i,1,:)),[.3,.3,1])
    hold on
    fill(fillx,squeeze(lafilltbtspecy(i,2,:)),[1,.3,.3])
    set(gca,'yscale','log')
    legend('light after l0','light after l1')
    
    subplot(2,2,3)
    fill(xfill,lfyl0(i,:),[.3,.3,1]);
    hold on
    fill(xfill,lfyl1(i,:),[1,.3,.3]);
    
    subplot(2,2,4)
    fill(xfill,layl0(i,:),[.3,.3,1]);
    hold on
    fill(xfill,layl1(i,:),[1,.3,.3]);
end


% spectrum iso vs cross running only
for i =1:size(lfr1filltbtspecy,1)
    figure
    
    subplot(2,2,1)
    fill(fillx,squeeze(lfr1filltbtspecy(i,1,:)),[.3,.3,1])
    hold on
    fill(fillx,squeeze(lfr1filltbtspecy(i,2,:)),[1,.3,.3])
    set(gca,'yscale','log')
    legend('light first l0','light first l1')    
    title(animalids{i})
    
    subplot(2,2,2)
    fill(fillx,squeeze(lar1filltbtspecy(i,1,:)),[.3,.3,1])
    hold on
    fill(fillx,squeeze(lar1filltbtspecy(i,2,:)),[1,.3,.3])
    set(gca,'yscale','log')
    legend('light after l0','light after l1')
    
    subplot(2,2,3)
    fill(xfill,lfyl0(i,:),[.3,.3,1]);
    hold on
    fill(xfill,lfyl1(i,:),[1,.3,.3]);
    
    subplot(2,2,4)
    fill(xfill,layl0(i,:),[.3,.3,1]);
    hold on
    fill(xfill,layl1(i,:),[1,.3,.3]);
end

figure
plot(1,lalgpow(:,1),'ko','markerfacecolor','k')
hold on
plot(2,lalgpow(:,2),'ro','markerfacecolor','r')
for i = 1:length(vi)
    line([1,2],[lalgpow(:,1),lalgpow(:,2)],'color','k')
end
set(gca,'xtick',[1,2]);
set(gca,'xticklabel',{'control','SOM Halo'})
title('PV: light after')
ylabel('fold gamma change')
set(gca,'xtick',[1,2]);
axis([0,3,0,800])
set(gcf,'OuterPosition',[573   504   242   513])


lared = lalgpow(:,2)./lalgpow(:,1);
lfred = lflgpow(:,2)./lflgpow(:,1);
rellared = rellapow(:,2)./rellapow(:,1);
rellfred = rellfpow(:,2)./rellfpow(:,1);

figure
plot(1,lared,'ro','markerfacecolor','r')
hold on
plot(2,lfred,'o','color',[1,.8,.8],'markerfacecolor',[1,.8,.8])
axis([0,3,min([lared;lfred])-.05*min([lared;lfred]),max([lfred;lared])+0.05*max([lfred;lared])])
for i = 1:length(vi)
    line([1,2],[lared(i),lfred(i)],'color','k')
end
set(gca,'xtick',[1,2]);
set(gca,'xticklabel',{'light after','light before'})
title('low gamma power change')
ylabel('fold gamma change')
set(gca,'xtick',[1,2]);
set(gcf,'OuterPosition',[573   504   242   513])

axis([0,3,0,1])
set(gca,'xticklabel',{'running','quiescent'})
set(gca,'xtick',[1,2]);
set(gca,'ytick',[0:.25:1]);
ylabel('percent peak gamma reduction')
title(['running vs quiescence p: ' num2str(p)]);
set(gcf,'OuterPosition',[573   504   242   513])

figure
plot(1,rellared,'ro')
hold on
plot(2,rellfred,'bo')
axis([0,3,min([rellared;rellfred])-.05*min([rellared;rellfred]),max([rellfred;rellared])+0.05*max([rellfred;rellared])])
for i = 1:length(vi)
    line([1,2],[rellared(i),rellfred(i)])
end
set(gca,'xtick',[1,2]);
set(gca,'xticklabel',{'light after','light before'})
title('relative low gamma power change')
ylabel('fold change')

plot(lfred,lared,'ro','markerfacecolor','r')
plot(lfred,lared,'bo','markerfacecolor','b')
legend('PV Halo','SOM Halo')
xlabel('gamma change light first')
ylabel('gamma change light after')
title('absolute low gamma power change')
axis([0,3.5,0,3.5])
hold on
plot([1,1],[0,3.5],'k')
plot([0,3.5],[1,1],'k')
axis square
set(gca,'xtick',[0,1,2,3]);
set(gca,'ytick',[0,1,2,3]);

plot(rellfred,rellared,'ro','markerfacecolor','r')
plot(rellfred,rellared,'bo','markerfacecolor','b')
legend('SOM Halo','PV Halo')
xlabel('relative gamma change light first')
ylabel('relative gamma change light after')
title('relative low gamma power change')
axis([0.4,1.5,0.4,1.5])
hold on
plot([0.4,1.5],[1,1],'k')
plot([1,1],[0.4,1.5],'k')
axis square

 disp('');