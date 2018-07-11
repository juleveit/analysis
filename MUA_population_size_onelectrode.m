function MUA_population_size_onelectrode

% % Scnn eArch pop
% animalids  = {'160219', '160222', '160223'};
% blocks     = [4,        2,        7];
% animal     = [1,        2,        3];
% electrodes = [[1,16];  [1,16];   [1,16]];
% bldepth      = [650,      550,      400];
% penangle   = [25,       25,       25];
% spacing    = [25,       25,       25];
% popfile    = 'C:\Users\Julia\work\data\populations\Scnn_eArch\size\MUA_population_singleel.mat';

% single thing
animalids =    {'141015'};
blocks    =     [5];
animal    =     [1];
electrodes =[[1,16]];
bldepth    =  [550];
penangle =      [25];
spacing =       [25];
popfile = 'C:\Users\Julia\work\data\populations\temp2.mat';

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
params.tapers = [2,3]; params.Fs = 1000; params.err = [2, 0.05]; params.trialave = 1;


if ~exist(popfile) || recalculate_pop

    cll = 1;
    for blck = 1:length(blocks)
        
        basepath = strcat('C:\Users\Julia\work\data\', animalids{blck}, '\');
        file = strcat(basepath, 'muaresult_', int2str(blocks(blck)), '_', int2str(electrodes(blck,1)), '_', int2str(electrodes(blck,2)), '.mat');
        
        clear result;
        if ~exist(file) || recalculate_muafile
            result = MUAdataprepare(basepath,animalids{blck},blocks(blck),electrodes(blck,1):electrodes(blck,2));
            save(file,'result')
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
            %                 msstamps([193,291]) = []; % for 151023 block 5
%                             msstamps([302,340]) = []; % for 160328 block 2
%                             result.msstamps = msstamps;
%                             surrresult.msstamps = msstamps;
%                             save(file,'result');
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
            speed(i,:) = result.runspeed(msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim);
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
            [pl(i,:),f] = pmtm(result.lfp(ch,msstamps(large(i))+300:msstamps(large(i))+1300),3,[],1000);
            [ps(i,:),f] = pmtm(result.lfp(ch,msstamps(small(i))+300:msstamps(small(i))+1300),3,[],1000);
        end
        plr1 = nan(1,size(pl,2)); psr1 = nan(1,size(pl,2));
        if length(lr1>=5)
            for i = 1:length(lr1)
                [plr1(i,:),f] = mtspectrumc(result.lfp(ch,msstamps(lr1(i))+700:msstamps(lr1(i))+1500),params);
                %                     [plr1(i,:),f] = pmtm(result.lfp(ch,msstamps(lr1(i))+300:msstamps(lr1(i))+1300),3,[],1000);
            end
        end
        if length(sr1>=5)
            for i = 1:length(sr1)
                [psr1(i,:),f] = mtspectrumc(result.lfp(ch,msstamps(sr1(i))+700:msstamps(sr1(i))+1500),params);
                %                     [psr1(i,:),f] = pmtm(result.lfp(ch,msstamps(sr1(i))+300:msstamps(sr1(i))+1300),3,[],1000);
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
        
        msStimes = round(result.msStimes{ch});
        if ~isempty(msStimes) & msStimes(1) == 0, msStimes(1) = 1; end
        
        chan_c = zeros(1,size(result.lfp,2));
        chan_c(msStimes) = 1;        
        
        for i = 1:length(msstamps)
            resp(i,:) = chan_c(msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim);
            hh = find(resp(i,1001:1800))'; ptresp(i).times = hh./1000;
            lfpresp(i,:) = result.lfp(ch,msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim);
            lfpspecttbt(i,:) = mtspectrumc(squeeze(lfpresp(i,1001:1800))',params);
            muacresp(i,:) =  result.muac(ch,msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim).^2;            
        end
        
        frs = sum(resp(:,respwin),2)./(length(respwin)/1000);
        bl = sum(resp(:,1:prestim),2)./(prestim/1000);
        muacmean = mean(muacresp(:,respwin),2);
        
        for l = 1:length(unique(light))
            for d = 1:length(sizes)
                
                thisinds = find(gratingInfo.size == sizes(d) & light == l-1);
                thisruninds = intersect(thisinds,oktrials); thisstillinds = intersect(thisinds,stilltrials);
                
                condresp(blck,l,d,:) = nanmean(resp(thisinds,:),1);
                condfr(blck,l,d) = nanmean(frs(thisinds));
                condmuacresp(blck,l,d,:) = nanmean(muacresp(thisinds,:));
                condmua(blck,l,d) = nanmean(muacmean(thisinds));
                condlfpresp(blck,l,d,:) = nanmean(lfpresp(thisinds,:));
                
                [lfpspect(blck,l,d,:),chf,lfpspecterr(blck,l,d,:,:)] = mtspectrumc(squeeze(lfpresp(thisinds,1001:1800))',params);
                tbtlfpspect(blck,l,d,:) = nanmean(lfpspecttbt(thisinds,:));
                tbtlfpspecterr(blck,l,d,:) = nanstd(lfpspecttbt(thisinds,:))./sqrt(length(thisinds));
                [muaspect(blck,l,d,:),chf,muaspecterr(blck,l,d,:,:)] = mtspectrumc(squeeze(muacresp(thisinds,1001:1800))',params);
                
                if ~isempty(thisruninds)
                    [r1_lfpspect(blck,l,d,:),chf,r1_lfpspecterr(blck,l,d,:,:)] = mtspectrumc(squeeze(lfpresp(thisruninds,1001:1800))',params);
                    r1_tbtlfpspect(blck,l,d,:) = nanmean(lfpspecttbt(thisruninds,:));
                    r1_tbtlfpspecterr(blck,l,d,:) = nanstd(lfpspecttbt(thisruninds,:))./sqrt(length(thisruninds));
                    r1_ntrials(blck,l,d) = length(thisruninds);
                else
                    r1_lfpspect(blck,l,d,:) = nan(1,length(lfpspect(blck,l,d,:)));
                    r1_lfpspecterr(blck,l,d,:,:) = nan(2,length(lfpspect(blck,l,d,:)));
                    r1_tbtlfpspect(blck,l,d,:) = nan(1,length(lfpspect(blck,l,d,:)));
                    r1_tbtlfpspecterr(blck,l,d,:) = nan(1,length(lfpspect(blck,l,d,:)));
                    r1_ntrials(blck,l,d) = 0;
                end
                
                if ~isempty(thisstillinds)
                    [r0_lfpspect(blck,l,d,:),chf,r0_lfpspecterr(blck,l,d,:,:)] = mtspectrumc(squeeze(lfpresp(thisstillinds,1001:1800))',params);
                    r0_tbtlfpspect(blck,l,d,:) = nanmean(lfpspecttbt(thisstillinds,:));
                    r0_tbtlfpspecterr(blck,l,d,:) = nanstd(lfpspecttbt(thisstillinds,:))./sqrt(length(thisstillinds));
                    r0_ntrials(blck,l,d) = length(thisstillinds);
                else
                    r0_lfpspect(blck,l,d,:) = nan(1,length(lfpspect(blck,l,d,:)));
                    r0_lfpspecterr(blck,l,d,:,:) = nan(2,length(lfpspect(blck,l,d,:)));
                    r0_tbtlfpspect(blck,l,d,:) = nan(1,length(lfpspect(blck,l,d,:)));
                    r0_tbtlfpspecterr(blck,l,d,:) = nan(1,length(lfpspect(blck,l,d,:)));
                    r0_ntrials(blck,l,d) = 0;
                end
            end
        end
    end
    save(popfile, '-v7.3');
else
    load(popfile);
end

vi = ones(size(r1_lfpspect,1))
for i = 1:length(vi)
%     lgpow(i,:,:) = squeeze(lfpspect(vi(i),di(vi(i)),:,:,lgi(vi(i),di(vi(i))))); % at second contact only
    r1lgpow(i,:,:) = squeeze(r1_lfpspect(vi(i),:,:,lgi(vi(i),di(vi(i)))));
    r0lgpow(i,:,:) = squeeze(r0_lfpspect(vi(i),:,:,lgi(vi(i),di(vi(i)))));
%     lglfpcoher(i,:,:) = squeeze(lfpcoherence(vi(i),di(vi(i)),:,:,lgi(vi(i),di(vi(i)))));
%     r1lglfpcoher(i,:,:) = squeeze(r1_lfpcoherence(vi(i),di(vi(i)),:,:,lgi(vi(i),di(vi(i)))));
%     r1tbtlglfpcoher(i,:,:) = squeeze(r1_tbtlfpcoherence(vi(i),:,:,lgi(vi(i),di(vi(i)))));
%     r0tbtlglfpcoher(i,:,:) = squeeze(r0_tbtlfpcoherence(vi(i),:,:,lgi(vi(i),di(vi(i)))));
%     lgsfcoher(i,:,:) = squeeze(sfcoher_c(vi(i),di(vi(i)),:,:,lgi(vi(i),di(vi(i)))));
%     r1lgsfcoher(i,:,:) = squeeze(r1_sfcoher_c(vi(i),di(vi(i)),:,:,lgi(vi(i),di(vi(i)))));
    
    for l = 1:2
        for d = 1:5
            r1fillspecy(i,l,d,:) = [squeeze(r1_lfpspecterr(vi(i),l,d,1,1:104))',fliplr(squeeze(r1_lfpspecterr(vi(i),l,d,2,1:104))')];            
            r1filltbtspecy(i,l,d,:) = [(squeeze(r1_tbtlfpspect(i,l,d,1:104))+squeeze(r1_tbtlfpspecterr(i,l,d,1:104)))',fliplr((squeeze(r1_tbtlfpspect(i,l,d,1:104))-squeeze(r1_tbtlfpspecterr(i,l,d,1:104)))')]';
%             r1fillcohery(i,l,d,:) = [squeeze(r1_lfpcoherr(vi(i),di(vi(i)),l,d,1,1:104))',fliplr(squeeze(r1_lfpcoherr(vi(i),2,l,d,2,1:104))')]';
%             r1filltbtcohery(i,l,d,:) = [(squeeze(r1_tbtlfpcoherence(vi(i),l,d,1:104))+squeeze(r1_tbtlfpcoherr(vi(i),di(vi(i)),l,d,1:104)))',fliplr((squeeze(r1_tbtlfpcoherence(vi(i),di(vi(i)),l,d,1:104))-squeeze(r1_tbtlfpcoherr(vi(i),di(vi(i)),l,d,1:104)))')]';
            r1fillspecy(i,l,d,:) = [squeeze(r1_lfpspecterr(vi(i),l,d,1,1:104))',fliplr(squeeze(r1_lfpspecterr(vi(i),l,d,2,1:104))')];            
            r0filltbtspecy(i,l,d,:) = [(squeeze(r0_tbtlfpspect(i,l,d,1:104))+squeeze(r0_tbtlfpspecterr(i,l,d,1:104)))',fliplr((squeeze(r0_tbtlfpspect(i,l,d,1:104))-squeeze(r0_tbtlfpspecterr(i,l,d,1:104)))')]';
%             r0filltbtcohery(i,l,d,:) = [(squeeze(r0_tbtlfpcoherence(vi(i),l,d,1:104))+squeeze(r0_tbtlfpcoherr(vi(i),di(vi(i)),l,d,1:104)))',fliplr((squeeze(r0_tbtlfpcoherence(vi(i),di(vi(i)),l,d,1:104))-squeeze(r0_tbtlfpcoherr(vi(i),di(vi(i)),l,d,1:104)))')]';
        end
    end    
end
fillx = [chf(1:104),fliplr(chf(1:104))];


% spectrum iso vs cross
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
    semilogy(chf,squeeze(lfpspect(i,3,1,1,:)),'b')
    hold on    
    semilogy(chf,squeeze(lfpspect(i,3,1,2,:)),'c')
    semilogy(chf,squeeze(lfpspect(i,3,1,3,:)),'g')
    semilogy(chf,squeeze(lfpspect(i,3,1,4,:)),'r')
    semilogy(chf,squeeze(lfpspect(i,3,1,5,:)),'m')
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
    semilogy(chf,squeeze(lfpspect(i,3,1,1,:)),'b')
    hold on    
    semilogy(chf,squeeze(lfpspect(i,3,1,5,:)),'c')
    semilogy(chf,squeeze(lfpspect(i,3,2,1,:)),'r')
    semilogy(chf,squeeze(lfpspect(i,3,2,5,:)),'m')
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
    semilogy(chf,squeeze(r1_lfpspect(i,3,1,1,:)),'b')
    hold on    
    semilogy(chf,squeeze(r1_lfpspect(i,3,1,2,:)),'c')
    semilogy(chf,squeeze(r1_lfpspect(i,3,1,3,:)),'g')
    semilogy(chf,squeeze(r1_lfpspect(i,3,1,4,:)),'r')
    semilogy(chf,squeeze(r1_lfpspect(i,3,1,5,:)),'m')
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
    semilogy(chf,squeeze(r1_lfpspect(i,3,1,1,:)),'b')
    hold on    
    semilogy(chf,squeeze(r1_lfpspect(i,3,1,5,:)),'c')
    semilogy(chf,squeeze(r1_lfpspect(i,3,2,1,:)),'r')
    semilogy(chf,squeeze(r1_lfpspect(i,3,2,5,:)),'m')
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
    imagesc(chf,depthax,log(squeeze(lfpspect(i,:,1,1,:))))
    axis([1,120,depthax(end)-10,depthax(1)+10])
    caxis([-2,9])
    title('power small')
    
    subplot(3,2,2)
    imagesc(chf,depthax,log(squeeze(lfpspect(i,:,1,5,:))))
    axis([1,120,depthax(end)-10,depthax(1)+10])
    caxis([-2,9])
    title('power large')
    
    subplot(3,2,3)
    imagesc(chf,depthax,log(squeeze(lfpspect(i,:,1,5,:)))-log(squeeze(lfpspect(i,:,1,1,:))));
    axis([1,120,depthax(end)-10,depthax(1)+10])
    caxis([-2,2])
    title('large - small') 
    
    subplot(3,2,5)
    imagesc(chf,depthax,log(squeeze(lfpspect(i,:,2,1,:)))-log(squeeze(lfpspect(i,:,1,1,:))));
    axis([1,120,depthax(end)-10,depthax(1)+10])
    caxis([-2,2])
    title('small L1-L0')
    
    subplot(3,2,6)
    imagesc(chf,depthax,log(squeeze(lfpspect(i,:,2,5,:)))-log(squeeze(lfpspect(i,:,1,5,:))));
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