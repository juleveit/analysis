function MUA_population_size_nolight


% SOM and PV
animalids = {'150331', '150401','150527','150529','150602','150603','150625','150825','150831','150902','150907','150909', '150915', '151023', '151027', '151109', '151110', '150629', '150730', '150731', '150804', '150818', '150820', '150823', '150824', '151104'};
blocks    = [3,         5,       11,      4,       5,       3,       6,       5,       4,       3,       3,       4,        3,        14,       3,        11,       13,       4,        4,        3,        4,        3,        3,        5,        3,        5];
animal    = [1,         2,       3,       4,       5,       6,       7,       8,       9,       10,      11,      12,       13,       14,       15,       16,       16,       17,       18,       19,       20,       21,       22,       23,       24,       25];
electrodes =[[1,32];    [1,32];  [1,32];  [1,32];  [1,32];  [1,16];  [1,16];  [17,32]; [1,16];  [1,16];  [1,16];  [1,16];   [17,32];  [17,32]; [17,23];  [17,23];  [1,16];   [1,16];   [1,16];   [1,16];   [1,16];   [1,16];   [1,16];   [1,16];   [17,32];  [17,32]];
bldepth    =[1050,      1000,    950,     750,     1000,    550,     450,     400,     400,     400,     400,     400,      400,      400,      500,      500,      500,      450,      300,      300,      400,      400,      400,      400,      400,      500];
penangle =  [25,        25,      25,      25,      25,      10,      10,      25,      25,      25,      25,      25,       25,       25,       25,       25,       25,       10,       10,       25,       25,       25,       25,       25,       25,       25];
spacing =   [25,        25,      25,      25,      25,      25,      25,      25,      25,      25,      25,      25,       25,       25,       25,       25,       25,       25,       25,       25,       25,       25,       25,       25,       25,       25];
popfile = 'C:\Users\Julia\work\data\populations\SOMPVsize_population.mat'; 

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
params.tapers = [5,9]; params.Fs = 1000; params.err = [2, 0.05]; params.trialave = 1;


if ~exist(popfile) || recalculate_pop

    cll = 1;
    for blck = 1:length(blocks)
        
        basepath = strcat('C:\Users\Julia\work\data\', animalids{blck}, '\');
        file = strcat(basepath, 'muaresult_', int2str(blocks(blck)), '_', int2str(electrodes(blck,1)), ':', int2str(electrodes(blck,2)), '.mat');
        
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
        
        for ch = 1: length(electrodes(blck,1):electrodes(blck,2))
            
            disp(['Block ' int2str(blck) '/' int2str(length(blocks)) '   channel ' int2str(ch) '/' int2str(length(electrodes(blck,1):electrodes(blck,2)))]);
            
            depth(ch) = bldepth(blck)-(ch-1)*spacing(blck);
            
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
            sizes = unique(result.gratingInfo.size); sizes(sizes==0) = [];
            large = find(result.gratingInfo.size == max(sizes) & result.light == 0);
            small = find(result.gratingInfo.size == min(sizes) & result.light == 0);
            lr1 = intersect(large,oktrials);
            sr1 = intersect(small,oktrials);
            for i = 1:length(large)
                [pl(i,:),f] = pmtm(result.lfp(ch,result.msstamps(large(i))+300:result.msstamps(large(i))+1300),3,[],1000);
                [ps(i,:),f] = pmtm(result.lfp(ch,result.msstamps(small(i))+300:result.msstamps(small(i))+1300),3,[],1000);
            end            
            plr1 = nan(1,size(pl,2)); psr1 = nan(1,size(pl,2));
            if length(lr1>=5)
                for i = 1:length(lr1)                    
                    [plr1(i,:),f] = pmtm(result.lfp(ch,result.msstamps(lr1(i))+300:result.msstamps(lr1(i))+1300),3,[],1000);
                end
            end
            if length(sr1>=5)
                for i = 1:length(sr1)                    
                    [psr1(i,:),f] = pmtm(result.lfp(ch,result.msstamps(sr1(i))+300:result.msstamps(sr1(i))+1300),3,[],1000);
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

            chan = zeros(1,size(result.lfp,2));
            chan(msStimes) = 1;
    
            for i = 1:length(msstamps)
                resp(i,:) = chan(msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim); 
                hh = find(resp(i,1001:1800))'; ptresp(i).times = hh./1000;                
                lfpresp(i,:) = result.lfp(ch,msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim);
                muacresp(i,:) =  result.muac(ch,msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim).^2;                
            end
           
            frs = sum(resp(:,respwin),2)./(length(respwin)/1000);
            bl = sum(resp(:,1:prestim),2)./(prestim/1000);
            muacmean = mean(muacresp(:,respwin),2);

            for l = 1:length(unique(result.light))
                for d = 1:length(sizes)
                    
                    thisinds = find(result.gratingInfo.size == sizes(d) & result.light == l-1);
                    thisruninds = intersect(thisinds,oktrials); thisstillinds = intersect(thisinds,stilltrials);
                  
                    condresp(blck,ch,l,d,:) = nanmean(resp(thisinds,:),1);
                    condfr(blck,ch,l,d) = nanmean(frs(thisinds));
                    condmuacresp(blck,ch,l,d,:) = nanmean(muacresp(thisinds,:));
                    condmua(blck,ch,l,d) = nanmean(muacmean(thisinds));
                    condlfpresp(blck,ch,l,d,:) = nanmean(lfpresp(thisinds,:));
                    
                    [lfpspect(blck,ch,l,d,:),chf,lfpspecterr(blck,ch,l,d,:,:)] = mtspectrumc(squeeze(lfpresp(thisinds,1001:1800))',params);
                    [muacspect(blck,ch,l,d,:),chf,muacspecterr_c(blck,ch,l,d,:,:)] = mtspectrumc(squeeze(muacresp(thisinds,1001:1800))',params);
                    
                    [sfcoher(blck,ch,l,d,:),a,b,c,de,sfcfx,e,f,g,sfcoherr_c(blck,ch,l,d,:,:)] = coherencycpt(lfpresp(thisinds,1001:1800)',ptresp(thisinds),params);
                    
                    if ~isempty(thisruninds)
                        [r1_lfpspect(blck,ch,l,d,:),chf,r1_lfpspecterr(blck,ch,l,d,:,:)] = mtspectrumc(squeeze(lfpresp(thisruninds,1001:1800))',params);
                        [r1_sfcoher(blck,ch,l,d,:),a,b,c,de,sfcfx,e,f,g,r1_sfcoherr_c(blck,ch,l,d,:,:)] = coherencycpt(lfpresp(thisruninds,1001:1800)',ptresp(thisruninds),params);
                        r1_ntrials(blck,ch,l,d) = length(thisruninds);
                    else
                        r1_lfpspect(blck,ch,l,d,:) = nan(1,length(lfpspect(blck,ch,l,d,:)));
                        r1_lfpspecterr(blck,ch,l,d,:,:) = nan(2,length(lfpspect(blck,ch,l,d,:)));
                        r1_sfcoher(blck,ch,l,d,:) = nan(1,length(sfcoher(blck,ch,l,d,:)));
                        r1_sfcohererr(blck,ch,l,d,:,:) = nan(2,length(sfcoher(blck,ch,l,d,:)));
                        r1_ntrials(blck,ch,l,d) = 0;
                    end
                    
                    if ~isempty(thisstillinds)                    
                        [r0_lfpspect(blck,ch,l,d,:),chf,r0_lfpspecterr(blck,ch,l,d,:,:)] = mtspectrumc(squeeze(lfpresp(thisstillinds,1001:1800))',params);
                        [r0_sfcoher(blck,ch,l,d,:),a,b,c,de,sfcfx,e,f,g,r0_sfcoherr_c(blck,ch,l,d,:,:)] = coherencycpt(lfpresp(thisstillinds,1001:1800)',ptresp(thisstillinds),params);
                        r0_ntrials(blck,ch,l,d) = length(thisstillinds);
                    else
                        r0_lfpspect(blck,ch,l,d,:) = nan(1,length(lfpspect(blck,ch,l,d,:)));
                        r0_lfpspecterr(blck,ch,l,d,:,:) = nan(2,length(lfpspect(blck,ch,l,d,:)));
                        r0_sfcoher(blck,ch,l,d,:) = nan(1,length(sfcoher(blck,ch,l,d,:)));
                        r0_sfcohererr(blck,ch,l,d,:,:) = nan(2,length(sfcoher(blck,ch,l,d,:)));
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


for i = 1:length(blocks)
    lgpow(i,:,:) = squeeze(lfpspect(i,di(i),:,:,lgi(i,di(i)))); % at second contact only
    r1lgpow(i,:,:) = squeeze(r1_lfpspect(i,di(i),:,:,lgi(i,di(i))));
    r1lgsfcoher(i,:,:) = squeeze(r1_sfcoher(i,di(i),:,:,lgi(i,di(i))));
    for l = 1:2
        for d = 1:5
            r1fillspecy(i,l,d,:) = [squeeze(r1_lfpspecterr(i,di(i),l,d,1,1:104))',fliplr(squeeze(r1_lfpspecterr(i,di(i),l,d,2,1:104))')];
        end
    end    
end
fillx = [chf(1:104),fliplr(chf(1:104))];

% spectrum small vs large
for i =1:size(r1fillspecy,1)
    figure
    fill(fillx,squeeze(r1fillspecy(i,1,1,:)),[.8,.8,1])
    hold on
    fill(fillx,squeeze(r1fillspecy(i,1,5,:)),[.3,.3,1])
    set(gca,'yscale','log')
    legend('iso','cross')
end

figure
plot(1,lgpow(:,1,1),'co')
hold on
plot(2,lgpow(:,1,5),'bo')
plot(3,lgpow(:,2,5),'ro')
axis([.5,3.5,0,550])
for i = 1:7
    plot([1,2],[lgpow(:,1,1),lgpow(:,1,5)],'k')
    plot([2,3],[lgpow(:,1,5),lgpow(:,2,5)],'k')
end
set(gca,'xtick',[1,2,3]);
set(gca,'xticklabel',{'small','large','large+light'})
title('low gamma power')
ylabel('psd at peak')

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
    plot(chf,squeeze(sfcoher(i,3,1,1,:)),'b');
    hold on
    plot(chf,squeeze(sfcoher(i,3,1,2,:)),'c');
    plot(chf,squeeze(sfcoher(i,3,1,3,:)),'g');
    plot(chf,squeeze(sfcoher(i,3,1,4,:)),'r');
    plot(chf,squeeze(sfcoher(i,3,1,5,:)),'m');
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
    plot(chf,squeeze(sfcoher(i,3,1,1,:)),'b');
    hold on
    plot(chf,squeeze(sfcoher(i,3,1,5,:)),'c');
    plot(chf,squeeze(sfcoher(i,3,2,1,:)),'r');
    plot(chf,squeeze(sfcoher(i,3,2,5,:)),'m');
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
    plot(chf,squeeze(r1_sfcoher(i,3,1,1,:)),'b');
    hold on
    plot(chf,squeeze(r1_sfcoher(i,3,1,2,:)),'c');
    plot(chf,squeeze(r1_sfcoher(i,3,1,3,:)),'g');
    plot(chf,squeeze(r1_sfcoher(i,3,1,4,:)),'r');
    plot(chf,squeeze(r1_sfcoher(i,3,1,5,:)),'m');
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
    plot(chf,squeeze(r1_sfcoher(i,3,1,1,:)),'b');
    hold on
    plot(chf,squeeze(r1_sfcoher(i,3,1,5,:)),'c');
    plot(chf,squeeze(r1_sfcoher(i,3,2,1,:)),'r');
    plot(chf,squeeze(r1_sfcoher(i,3,2,5,:)),'m');
    axis([0,120,0,1]);
    xlabel('frequency')
    ylabel('spike field coherence')
    set(gcf,'OuterPosition',[50,300,1200,700])
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
    imagesc(chf,depthax,squeeze(sfcoher(i,:,1,1,:)))
    axis([1,120,depthax(end)-10,depthax(1)+10])
    caxis([0,.6])
    title('SFC small')
    
    subplot(3,2,2)
    imagesc(chf,depthax,squeeze(sfcoher(i,:,1,5,:)))
    axis([1,120,depthax(end)-10,depthax(1)+10])
    caxis([0,.6])
    title('SFC large') 
    
    subplot(3,2,3)
    imagesc(chf,depthax,squeeze(sfcoher(i,:,1,5,:))-squeeze(sfcoher(i,:,1,1,:)))
    axis([1,120,depthax(end)-10,depthax(1)+10])
    caxis([-.2,.2])
    title('SFC large - small')
    
    subplot(3,2,5)    
    imagesc(chf,depthax,squeeze(sfcoher(i,:,2,1,:))-squeeze(sfcoher(i,:,1,1,:)))
    axis([1,120,depthax(end)-10,depthax(1)+10])
    caxis([-.2,.2])
    title('small L1-L0')
    
    subplot(3,2,6)    
    imagesc(chf,depthax,squeeze(sfcoher(i,:,2,5,:))-squeeze(sfcoher(i,:,1,5,:)))
    axis([1,120,depthax(end)-10,depthax(1)+10])
    caxis([-.2,.2])
    title('large L1-L0')

end


 disp('');