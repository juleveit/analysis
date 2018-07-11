function orisurround_SU_allcells

animalid = '180601';
block = 6;
rfblock = 5;
lcol = 'r'; %lasercolor

onlymod = 0;
printyn = 1;
sfc = 0;

%Scott Kernel
%%%%%%%%%%%%kernel propeties%%%%%%%%%%%%%%%
kernel_width_eval_s = 0.15;
sdf_freq_hz = 3700;
exp_growth_ms = 2;
exp_decay_ms = 12;
exp_growth_s = exp_growth_ms/1000;
exp_decay_s = exp_decay_ms/1000;
eval_kernel_x_s = 0:1/sdf_freq_hz:kernel_width_eval_s;
exp_kernel = eval_kernel_x_s;
exp_kernel = (1-(exp(-(exp_kernel./exp_growth_s)))).*(exp(-(exp_kernel./exp_decay_s)));
excit_kernel = exp_kernel/sum(exp_kernel); 
%%%%%%%%%%kernal properties%%%%%%%%%%%%%%%

supath = ['C:\Users\Julia\work\data\' animalid '\singleunits\'];
basename = [animalid '_block' int2str(block) '_tet'];
snbasename = [animalid '_block' int2str(rfblock) '_tet'];

files = dir([supath, basename, '*.mat']);
rffiles = dir([supath, snbasename, '*.mat']);

prestim = 300;
poststim = 700;
respwin = 501:1500; % after stimulus onset
respwin = respwin+prestim;
freqbinwidth = 5;

cell = 1;
for fi = 1:length(files)
    
%     if strfind(files(fi).name, 'MU')
%         continue;
%     end
        
    load([supath, files(fi).name]);
%     disp(['now analyzing file: ' files(cell).name]);
    cellname{cell} = files(fi).name;
    
    i = strfind(files(fi).name, 'tet');
    tetno(cell) = strread(files(fi).name(i+3));
    
    wvchan = find(var(result.waveforms) == max(var(result.waveforms)));
    cm = [3,4,1,2]; % confusion matrix? take lfp from two electrodes away to not get too many spike related phase resets
    lfp = result.lfp(:,cm(wvchan))';
        
    sr = 1000;
    
    %find gamma peaks for this animal
    beta = [15,40];
    gamma = [50,70];
    large = find(result.gratingInfo.Orientation_surround ~= -1 & result.light == 0);
    small = find(result.gratingInfo.Orientation_surround == -1 & result.light == 0);
    for i = 1:length(small)
        [pl(i,:),f] = pmtm(lfp(result.msstamps(large(i)):result.msstamps(large(i))+1000),3,[],1000);
        [ps(i,:),f] = pmtm(lfp(result.msstamps(small(i)):result.msstamps(small(i))+1000),3,[],1000);
    end
    b1 = find(f>beta(1),1); b2 = find(f>beta(2),1);
    g1 = find(f>gamma(1),1); g2 = find(f>gamma(2),1);
    bsig = nanmean(pl(:,b1:b2));
    gsig = nanmean(ps(:,g1:g2));
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
        
    msStimes = round(result.spikes);
    if isempty(msStimes), msStimes(1) = 0; end
   if msStimes(1) == 0, msStimes(1) = 1; end  
    
    chan = zeros(1,length(result.lfp));
    chan(msStimes) = 1;
    
    trialdur = result.stimduration*1000;
    msstamps = result.msstamps;
    
     if length(msstamps)~=length(result.light)
%          disp('');
% %         msstamps([62,108,147]) = []; % for 140703 block 8
% %         msstamps([161]) = []; % for 141204 block 3
%         msstamps([303]) = []; % for 150407 block 5
%         msstamps([84,227]) = []; % for 150629 block 8
%         msstamps([100]) = []; % for 150902 block 10
%         msstamps([260]) = []; % for 150902 block 10
%         msstamps([193,291]) = []; % for 151023 block 5
%         msstamps(205) = []; % for 151109 block 9
%         result.msstamps = msstamps;
%         save([supath, files(fi).name],'result');
        pause;
    end
    
    
    for i = 1:length(msstamps)
        resp(i,:) = chan(msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim); 
        lfpresp(i,:) = result.lfp(msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim);
        [lfpspect(i,:),trialfax] = pmtm(lfpresp(i,1001:1800),3,[],sr);
        gamma1resp(i,:) = gpow1(msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim);
        gamma2resp(i,:) = gpow2(msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim);
        gamma3resp(i,:) = gpow3(msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim);
        
        if sfc
            for j = 1:size(phasmat,1)
                allphaseresp(j,i,:) = phasmat(j, msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim);
                allpowresp(j,i,:) = powmat(j, msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim);
            end
        end
        
        speed(i,:) = result.runspeed(msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim);     
    end
    
    
    depth(cell) = result.depth;
    
    spike = result.waveforms(:,wvchan);
    interpspike = spline(1:32,spike,1:.1:32);
    [adiff(cell),swidth(cell)] = spikequant(interpspike);
    
    cellrespl0(cell,:,:) = resp(find(result.light == 0),:);
    cellrespl1(cell,:,:) = resp(find(result.light == 1),:);
    cellresp(cell,:,:) = resp;
    celllfpresp(cell,:,:) = lfpresp;
    
    % figure out sufficiently high and nonvariable runspeed trials
    meanspeed = mean(speed(:,respwin),2);
    stdspeed = std(speed(:,respwin),1,2);
    notstill = find(meanspeed>1);
    okspeed = find(meanspeed>( mean(meanspeed(notstill))-(1.5*std(meanspeed(notstill))) ) );
    okvar = find(stdspeed<( mean(stdspeed(notstill))+(1.5*std(stdspeed(notstill)))) & stdspeed>.5);
    oktrials = intersect(okspeed,okvar);
    nonoktrials = 1:size(resp,1); nonoktrials(oktrials) = [];
    stilltrials = 1:size(resp,1); stilltrials(notstill) = [];
    
    msta = linspace(-prestim,trialdur+poststim,size(resp,2));
    
    lightresp = resp(find(result.light),:);
    nolightresp = resp(find(result.light == 0),:);
    
    lightlfpresp = lfpresp(find(result.light),:);
    nolightlfpresp = lfpresp(find(~result.light),:);
    
    frs = sum(resp(:,respwin),2)./(length(respwin)/1000);
    bl = sum(resp(:,1:prestim),2)./(prestim/1000);
    sc = sum(resp(:,respwin),2);
    gp1 = mean(gamma1resp(:,respwin),2);
    gp2 = mean(gamma2resp(:,respwin),2);
    gp3 = mean(gamma3resp(:,respwin),2);
        
     %determine if cell is visually modulated
    blfr = sum(resp(:,1:prestim),2);
    vrfr = sum(resp(:,prestim+40:2*prestim+40),2);
    vismod(cell) = ttest2(blfr,vrfr);
    
    %determine if cell is modulated by light
    lightmod(cell) = ttest2(frs(find(result.light)),frs(find(result.light == 0)));
    
    lfr(cell) = mean(frs(find(result.light)));
    nlfr(cell) = mean(frs(find(result.light == 0)));
    
   
    
    if sfc
        tmpi = zeros(size(allphaseresp));
        for i = 1:size(allphaseresp,1)
            tmp = zeros(size(gphaseresp));
            tmp(find(resp)) = allphaseresp(i,find(resp));
            tmpi(i,:,:) = tmp;
        end
        alll0phasemat = tmpi(:,find(result.light == 0),:);
        alll1phasemat = tmpi(:,find(result.light == 1),:);
        for i = 1:size(allphaseresp,1)
            allphasesl0{i} = alll0phasemat(i,find(squeeze(alll0phasemat(i,:,:)))); 
            allphasesl1{i} = alll1phasemat(i,find(squeeze(alll1phasemat(i,:,:)))); 
            allrl0(cell,i) = circ_r(allphasesl0{i}');
            allrl1(cell,i) = circ_r(allphasesl1{i}');
            allcmeanl0(cell,i) = circ_mean(allphasesl0{i}');
            allcmeanl1(cell,i) = circ_mean(allphasesl1{i}');
        end
    end
    
    
    binwidth = 25;
    [binnedlight,bta] = binit(mean(lightresp),binwidth);
    binnedlight = binnedlight.*(1000/binwidth);
    [binnednolight,bta] = binit(mean(nolightresp),binwidth);
    binnednolight = binnednolight.*(1000/binwidth);    
    
    %get RF of cell, too
    [on,off,gaussfiton,gaussfitoff,rson(cell),rsoff(cell),xax,yax] = get_rf([supath, rffiles(fi).name]);    
    
    printname = files(fi).name;
    printname(find(printname=='_')) = ' ';
    
    centeroris = unique(result.gratingInfo.Orientation_center);  centeroris(find(centeroris == -1)) = []; %delete control condition
%     help = result.gratingInfo.Orientation_surround - 360; result.gratingInfo.Orientation_surround(result.gratingInfo.Orientation_surround>315) = help(result.gratingInfo.Orientation_surround>315);
    surroundoris = unique(result.gratingInfo.Orientation_surround); surroundoris(find(surroundoris == -1)) = []; 
%     ordf = [0,90,-1,-2];  % 1: iso, 2: cross, 3: center only, 4: surround only
    ordf = [0,30,60,90,-1,-2];  % 1: iso, 2: 30 diff, 3: 60 diff, 4: cross, 5: center only, 6: surround only
    for l = 1:length(unique(result.light))
        for ori = 1:length(centeroris)
            for d = 1:length(ordf)
                if ordf(d) == -1
                    thisinds = find(result.gratingInfo.Orientation_center == centeroris(ori) & ...
                        result.gratingInfo.Orientation_surround == -1 & ...
                        result.light == l-1);
                elseif ordf(d) == -2
                    thisinds = find(result.gratingInfo.Orientation_center == -1 & ...
                        result.gratingInfo.Orientation_surround == centeroris(ori) & ...
                        result.light == l-1);
                else
                    thisinds = find(result.gratingInfo.Orientation_center == centeroris(ori) &...
                        result.gratingInfo.Orientation_surround == centeroris(ori)+ordf(d) & ...
                        result.light == l-1);
                end
                condresp(l,ori,d,:) = mean(resp(thisinds,:),1);
                condlfpspect(l,ori,d,:) = nanmean(lfpspect(thisinds,:));
                condlfpspecterr(l,ori,d,:) = nanstd(lfpspect(thisinds,:))./sqrt(length(thisinds));
                condstspect(l,ori,d,:) = pmtm(mean(resp(thisinds,1001:1800),1),3,[],sr);
                condfr(l,ori,d) = mean(frs(thisinds));%-mean(bl);        
                conderr(l,ori,d) =std(frs(thisinds))./sqrt(length(thisinds));
                condgp1(cell,l,ori,d) = mean(gp1(thisinds));     
                condgp1err(cell,l,ori,d) =std(gp1(thisinds))./sqrt(length(thisinds));
                condgp2(cell,l,ori,d) = mean(gp2(thisinds));     
                condgp2err(cell,l,ori,d) =std(gp2(thisinds))./sqrt(length(thisinds));
                condgp3(cell,l,ori,d) = mean(gp3(thisinds));     
                condgp3err(cell,l,ori,d) =std(gp3(thisinds))./sqrt(length(thisinds));
                
                condfiltresp(cell,l,ori,d,:) = filter(excit_kernel,1,mean(resp(thisinds,:),1));
                
                condz(l,ori,d) = {(sc(thisinds)-mean(sc(thisinds)))/std(sc(thisinds))}; %ecker 2010
                condsc(l,ori,d) = {sc(thisinds)};
                ff(l,ori,d) = var(sc(thisinds))/mean(sc(thisinds));
                
                condlfpresp(l,ori,d,:) = mean(lfpresp(thisinds,:),1);                
                
                condresperr(l,ori,d,:) = nanstd(resp(thisinds,:),1,1)./sqrt(length(thisinds));
                if ~isnan(condresp(l,ori,d,:))
                    [bincondresp(l,ori,d,:),bta] = binit(condresp(l,ori,d,:),binwidth);
                else
                    bincondresp(l,ori,d,:) = binit(condresp(l,ori,d,:),binwidth);
                end
                binconderr(l,ori,d,:) = binit(condresperr(l,ori,d,:),binwidth);
                
                mscc = []; bincc = [];
                for ii = 1:length(thisinds)-1
                    for jj = ii+1:length(thisinds)
                        help = corrcoef(resp(thisinds(ii),:),resp(thisinds(jj),:));
                        mscc = [mscc,help(1,2)];
                        help = corrcoef(binit(resp(thisinds(ii),:),binwidth),binit(resp(thisinds(jj),:),binwidth));
                        bincc = [bincc, help(1,2)];
                    end
                end
                msreliab(l,ori,d) = nanmean(mscc);
                binreliab(l,ori,d) = nanmean(bincc);
                eckerreliability(l,ori,d) = var(frs(thisinds))/var(frs);

                thisruninds = intersect(thisinds,oktrials);
                if ~isempty(thisruninds)
                    runcondresp(l,ori,d,:) = mean(resp(thisruninds,:),1);
                    runcondfr(l,ori,d) = mean(frs(thisruninds));
                    runconderr(l,ori,d) = std(frs(thisruninds))./sqrt(length(thisruninds));
                else
                    runcondresp(l,ori,d,:) = nan(1,size(resp,2));
                    runcondfr(l,ori,d) = NaN;
                    runconderr(l,ori,d) = NaN;
                end  
                
                thisstillinds = intersect(thisinds,stilltrials);
                if ~isempty(thisstillinds)
                    stillcondresp(l,ori,d,:) = mean(resp(thisstillinds,:),1);
                    stillcondfr(l,ori,d) = mean(frs(thisstillinds));
                    stillconderr(l,ori,d) = std(frs(thisstillinds))./sqrt(length(thisstillinds));
                else
                    stillcondresp(l,ori,d,:) = nan(1,size(resp,2));
                    stillcondfr(l,ori,d) = NaN;
                    stillconderr(l,ori,d) = NaN;
                end  
                
            end
        end
        for d = 1:length(ordf)
            if ordf(d) == -1
                thisinds = find(result.gratingInfo.Orientation_center ~= -1 & ...
                    result.gratingInfo.Orientation_surround == -1 & ...
                    result.light == l-1);
            elseif ordf(d) == -2
                thisinds = find(result.gratingInfo.Orientation_center == -1 & ...
                    result.gratingInfo.Orientation_surround ~= -1 & ...
                    result.light == l-1);
            elseif ordf(d) == 0
                thisinds = find(result.gratingInfo.Orientation_center == ...
                    result.gratingInfo.Orientation_surround & ...
                    result.gratingInfo.Orientation_center ~= -1 &...
                    result.light == l-1);
            else
                thisinds = find(abs(result.gratingInfo.Orientation_center - ...
                    result.gratingInfo.Orientation_surround) == ordf(d) & ...
                    result.light == l-1);
            end
            allcondlfpspect(l,d,:,:) = lfpspect(thisinds,:);
            allcondlfpresp(l,d,:,:) = lfpresp(thisinds,:);
            
            params.Fs = 1000; params.trialave = 1; params.err = [2 .05]; params.tapers = [5,9];
            if ~isempty(thisinds)
                [S,chf,Serr]=mtspectrumc(squeeze(lfpresp(thisinds,1001:1800))',params);
                condS(cell,l,d,:) = S(1:150); condSerr(cell,l,d,:,:) = Serr(:,1:150);
            end
        end
    end
    sfx = chf(1:150);
    
    bincondresp = bincondresp.*(1000/binwidth);
    bta = bta-prestim;
%     psthcondplot(bincondresp,binconderr,bta)
    cellbinresp(cell,:,:,:,:) = bincondresp;
    
    contindsnl = find(result.gratingInfo.Orientation_surround == -1 &...
        result.gratingInfo.Orientation_center == -1 & result.light == 0);
    controlresp(1,:) = mean(resp(contindsnl,:),1);
    controlfr(cell,1) = mean(frs(contindsnl));
    controlerr(1) = std(frs(contindsnl))./sqrt(length(contindsnl));
    controlgp1(cell,1) = mean(gp1(contindsnl));
    controlgp1err(cell,1) = std(gp1(contindsnl))./sqrt(length(contindsnl));
    controlgp2(cell,1) = mean(gp2(contindsnl));
    controlgp2err(cell,1) = std(gp2(contindsnl))./sqrt(length(contindsnl));
    controlgp3(cell,1) = mean(gp3(contindsnl));
    controlgp3err(cell,1) = std(gp3(contindsnl))./sqrt(length(contindsnl));
    controllfpspect(cell,1,:) = nanmean(lfpspect(contindsnl,:),1);
    
    contindsl =  find(result.gratingInfo.Orientation_surround == -1 &...
        result.gratingInfo.Orientation_center == -1 & result.light == 1);
    controlresp(2,:) = mean(resp(contindsl,:),1);
    controlfr(cell,2) = mean(frs(contindsl));
    controlerr(2) = std(frs(contindsl))./sqrt(length(contindsl));
    controlgp1(cell,2) = mean(gp1(contindsl));
    controlgp1err(cell,2) = std(gp1(contindsl))./sqrt(length(contindsl));
    controlgp2(cell,2) = mean(gp2(contindsl));
    controlgp2err(cell,2) = std(gp2(contindsl))./sqrt(length(contindsl));
    controlgp3(cell,2) = mean(gp3(contindsl));
    controlgp3err(cell,2) = std(gp3(contindsl))./sqrt(length(contindsl));
    controllfpspect(cell,2,:) = nanmean(lfpspect(contindsl,:),1);

    
    [binnedctrnolight,bta] = binit(controlresp(1,:),binwidth);
    binnedctrlight = binit(controlresp(2,:),binwidth);

    binnedcellrespl0(cell,:) = binnednolight;
    binnedcellrespl1(cell,:) = binnedlight;
    
    nolmaxfr(cell) = max(max(condfr(1,:,:)));
%     lmaxfr(cell) = max(max(condfr(2,:,:)));
    
    cellz(cell,:,:,:) = condz;
    cellsc(cell,:,:,:) = condsc;
    cellff(cell,:,:,:) = ff;
    cellfr(cell,:,:,:) = condfr;
    cellerr(cell,:,:,:) = conderr;
    celleckerrely(cell,:,:,:) = eckerreliability;
    cellmsrely(cell,:,:,:) = msreliab;
    cellbinrely(cell,:,:,:) = binreliab;
    celllfpspect(cell,:,:,:,:) = condlfpspect;
    celllfpspecterr(cell,:,:,:,:) = condlfpspecterr;
    cellalllfpspect(cell,:,:,:,:) = allcondlfpspect;
    cellalllfpresp(cell,:,:,:,:) = allcondlfpresp;
%     figure
%     errorbar(sizes,squeeze(nanmean(condfr(2,:,:),2)),squeeze(nanmean(conderr(2,:,:),2)),'o-','color',lcol,'markersize',8,'linewidth',2)
%     hold on
%     errorbar(sizes,squeeze(nanmean(condfr(1,:,:),2)),squeeze(nanmean(conderr(1,:,:),2)),'ko-','markersize',8,'linewidth',2)
%     xlabel('shown patch size [vd]')
%     ylabel('Firing rate [Hz]')
%     legend({'Light ON','Light OFF'})    
%     set(gca,'xtick',sizes)  
%     title(['cell ' int2str(cell) ' average all orientations'])
    
    prefsize = find(mean(condfr(1,:,:),2) == max(mean(condfr(1,:,:),2)),1);
    
    for l = 1:length(unique(result.light))
        for cond = 1:size(condfr,3)
            if cond == 4, 
                oris = centeroris(1:4); 
                rp = squeeze(condfr(l,1:4,cond));
            else
                oris = centeroris;
                rp = squeeze(condfr(l,:,cond));
            end
            [oriprefratio(cell,l,cond), dirprefratio(cell,l,cond), prefori(cell,l,cond), meanori(cell,l,cond), osi(cell,l,cond), meandir(cell,l,cond), dsi(cell,l,cond)] = getOSI(rp,oris);
        end
    end   
    
    figure
    subplot(2,2,1)
    ta = bta-prestim;
    plot(ta,binnednolight,'k','linewidth',2);
    hold on
    plot(ta,binnedlight,lcol,'linewidth',2);
    mx = max([max(binnednolight),max(binnedlight),.01]);
    axis([-prestim,trialdur+poststim,0,mx]);
    line([0,0],[0,mx],'color','k','linewidth',2);
    line([2000,2000],[0,mx],'color','k','linewidth',2);
    line([500,500],[0,mx],'color','b','linewidth',2)
    line([1500,1500],[0,mx],'color','b','linewidth',2);
    legend({'Light OFF','Light ON'})
    xlabel('time [ms]')
    ylabel('firing rate [Hz]')
    title(['cell ' int2str(cell) ' depth: ' int2str(result.depth), 'cell ' printname ])
    
    colors = ['g','b','m','r'];
    nodiffs = size(condfr,3)-2;
    subplot(2,2,2)
    errorbar(centeroris,squeeze(condfr(1,:,nodiffs+1)),squeeze(conderr(1,:,nodiffs+1)),'ko-','markersize',8,'linewidth',2,'markerfacecolor','k')
    hold on
    for i = 1:nodiffs
        errorbar(centeroris,squeeze(condfr(1,:,i)),squeeze(conderr(1,:,i)),'o-','color',colors(i),'markersize',8,'linewidth',2,'markerfacecolor',colors(i))
    end
    errorbar(centeroris,squeeze(condfr(1,:,nodiffs+2)),squeeze(conderr(1,:,nodiffs+2)),'co-','markersize',8,'linewidth',2,'markerfacecolor','c')
    errorbar(-15,controlfr(1),controlerr(1),'ko','markersize',8,'linewidth',2,'markerfacecolor','k')
    xlabel('shown orientation')
    ylabel('Firing rate [Hz]')
    set(gca,'xtick',oris)
    ordfstr{1} = 'center only'; ordfstr{nodiffs+2} = 'surround only';
    for i = 2:nodiffs+1
        ordfstr{i} = int2str(ordf(i-1));
    end
    legend(ordfstr)
    title([' OSI co: ' num2str(osi(cell,1,3)) ' OSI is: ' num2str(osi(cell,1,1)) ' OSI cs: ' num2str(osi(cell,1,2)) ' OSI so: ' num2str(osi(cell,1,4))])
    
    subplot(2,2,3)
    if length(unique(result.light))>1
        bars = [controlfr(cell,1),nanmean(condfr(1,:,nodiffs+1),2),squeeze(nanmean(condfr(1,:,1:nodiffs),2))',nanmean(condfr(1,:,nodiffs+2),2);...
            controlfr(cell,2),nanmean(condfr(2,:,nodiffs+1),2),squeeze(nanmean(condfr(2,:,1:nodiffs),2))',nanmean(condfr(2,:,nodiffs+2),2)];
        errorbars = [controlerr(1),nanmean(conderr(1,:,nodiffs+1),2),squeeze(nanmean(conderr(1,:,1:nodiffs),2))',nanmean(conderr(1,:,nodiffs+2),2);...
            controlerr(2),nanmean(conderr(2,:,nodiffs+1),2),squeeze(nanmean(conderr(2,:,1:nodiffs),2))',nanmean(conderr(2,:,nodiffs+2),2)];
%         bars = [controlfr(cell,1),nanmean(condfr(1,:,3),2),nanmean(condfr(1,:,1),2),nanmean(condfr(1,:,2),2),nanmean(condfr(1,:,4),2);...
%             controlfr(cell,2),nanmean(condfr(2,:,3),2),nanmean(condfr(2,:,1),2),nanmean(condfr(2,:,2),2),nanmean(condfr(2,:,4),2)];
%         errorbars = [controlerr(1),nanmean(conderr(1,:,3),2),nanmean(conderr(1,:,1),2),nanmean(conderr(1,:,2),2),nanmean(conderr(1,:,4),2);...
%             controlerr(2),nanmean(conderr(2,:,3),2),nanmean(conderr(2,:,1),2),nanmean(conderr(2,:,2),2),nanmean(conderr(2,:,4),2)];
    else
        bars = [controlfr(cell,1),nanmean(condfr(1,:,3),2),nanmean(condfr(1,:,1),2),nanmean(condfr(1,:,2),2),nanmean(condfr(1,:,4),2)];
        errorbars = [controlerr(1),nanmean(conderr(1,:,3),2),nanmean(conderr(1,:,1),2),nanmean(conderr(1,:,2),2),nanmean(conderr(1,:,4),2)];
    end
%     bar(bars')
%     hold on
%     errorbar(bars',errorbars','.');
    barweb(bars',errorbars',[],['gray',ordfstr],'average firing rates','condition','firing rate',[],[],{'no light','light'});
    xlabel('condition')
    ylabel('Firing rate [Hz]')
    set(gca,'xtick',1:length(ordfstr)+1)
    set(gca,'xticklabel',[{'control'}, ordfstr]); %[{'control'},{'co'},{'is'},{'cs'},{'so'}])  
    
    subplot(2,4,7)
    plot(spike)
    axis([0,40,-100,100])
    legend(['width: ' int2str(swidth(cell)) ' adiff: ' num2str(adiff(cell))])
    
    subplot(2,4,8)
    plot(ta,binnedctrnolight)
    hold on
    plot(ta,binnedctrlight,lcol)
    axis([0,2500,0,1])
    
    figure    
    subplot(2,1,1)
    hold on
    plot(bta,squeeze(nanmean(bincondresp(1,:,nodiffs+1,:))),'k')
    for i = 1:nodiffs
        plot(bta,squeeze(nanmean(bincondresp(1,:,i,:))),'color',colors(i));
    end
    plot(bta,squeeze(nanmean(bincondresp(1,:,nodiffs+2,:))),'c')
    legend(ordfstr)
    
    subplot(2,2,3)    
    plot_circle(result.position(1),result.position(2),max(unique(result.gratingInfo.sizesurround))/2,'k',2)
    hold on
    imagesc(xax,yax,on)
    plot_orrf_absdeg(gaussfiton,1,'w',2)
    axis square
    plot_circle(result.position(1),result.position(2),max(unique(result.gratingInfo.sizecenter))/2,'k',2)
    
    subplot(2,2,4)    
    plot_circle(result.position(1),result.position(2),max(unique(result.gratingInfo.sizesurround))/2,'k',2)
    hold on
    imagesc(xax,yax,off)    
    plot_circle(result.position(1),result.position(2),max(unique(result.gratingInfo.sizecenter))/2,'k',2)
    plot_orrf_absdeg(gaussfitoff,1,'w',2)
    axis square
    
%     figure    
%     subplot(2,2,1)
%     bars = [controlfr(cell,:);nanmean(condfr(1,:,3),2),nanmean(condfr(2,:,3),2);...
%         nanmean(condfr(1,:,1),2),nanmean(condfr(2,:,1),2);nanmean(condfr(1,:,2),2),nanmean(condfr(2,:,2),2);...
%         nanmean(condfr(1,:,4),2),nanmean(condfr(2,:,4),2)];
%     errorbars = [controlerr;nanmean(conderr(1,:,3),2),nanmean(conderr(2,:,3),2);...
%         nanmean(conderr(1,:,1),2),nanmean(conderr(2,:,1),2);nanmean(conderr(1,:,2),2),nanmean(conderr(2,:,2),2);...
%         nanmean(conderr(1,:,4),2),nanmean(conderr(2,:,4),2)];
%     barweb(bars',errorbars',[],{'control','light on'},'firing rates','condition','firing rate',[],[],{'cntr','co','is','cs','so'});
%     
%     subplot(2,2,2)
%     bars = [controlgp1(cell,:);nanmean(condgp1(cell,1,:,3),3),nanmean(condgp1(cell,2,:,3),3);...
%         nanmean(condgp1(cell,1,:,1),3),nanmean(condgp1(cell,2,:,1),3);nanmean(condgp1(cell,1,:,2),3),nanmean(condgp1(cell,2,:,2),3);...
%         nanmean(condgp1(cell,1,:,4),3),nanmean(condgp1(cell,2,:,4),3)];
%     errorbars = [controlgp1err(cell,:);nanmean(condgp1err(cell,1,:,3),3),nanmean(condgp1err(cell,2,:,3),3);...
%         nanmean(condgp1err(cell,1,:,1),3),nanmean(condgp1err(cell,2,:,1),3);nanmean(condgp1err(cell,1,:,2),3),nanmean(condgp1err(cell,2,:,2),3);...
%         nanmean(condgp1err(cell,1,:,4),3),nanmean(condgp1err(cell,2,:,4),3)];
%     barweb(bars',errorbars',[],{'control','light on'},'gamma1 power','condition','power',[],[],{'cntr','co','is','cs','so'});
%     
%     subplot(2,2,3)
%     bars = [controlgp2(cell,:);nanmean(condgp2(cell,1,:,3),3),nanmean(condgp2(cell,2,:,3),3);...
%         nanmean(condgp2(cell,1,:,1),3),nanmean(condgp2(cell,2,:,1),3);nanmean(condgp2(cell,1,:,2),3),nanmean(condgp2(cell,2,:,2),3);...
%         nanmean(condgp2(cell,1,:,4),3),nanmean(condgp2(cell,2,:,4),3)];
%     errorbars = [controlgp2err(cell,:);nanmean(condgp2err(cell,1,:,3),3),nanmean(condgp2err(cell,2,:,3),3);...
%         nanmean(condgp2err(cell,1,:,1),3),nanmean(condgp2err(cell,2,:,1),3);nanmean(condgp2err(cell,1,:,2),3),nanmean(condgp2err(cell,2,:,2),3);...
%         nanmean(condgp2err(cell,1,:,4),3),nanmean(condgp2err(cell,2,:,4),3)];
%     barweb(bars',errorbars',[],{'control','light on'},'gamma2 power','condition','power',[],[],{'cntr','co','is','cs','so'});
%     
%     subplot(2,2,4)
%     bars = [controlgp3(cell,:);nanmean(condgp3(cell,1,:,3),3),nanmean(condgp3(cell,2,:,3),3);...
%         nanmean(condgp3(cell,1,:,1),3),nanmean(condgp3(cell,2,:,1),3);nanmean(condgp3(cell,1,:,2),3),nanmean(condgp3(cell,2,:,2),3);...
%         nanmean(condgp3(cell,1,:,4),3),nanmean(condgp3(cell,2,:,4),3)];
%     errorbars = [controlgp3err(cell,:);nanmean(condgp3err(cell,1,:,3),3),nanmean(condgp3err(cell,2,:,3),3);...
%         nanmean(condgp3err(cell,1,:,1),3),nanmean(condgp3err(cell,2,:,1),3);nanmean(condgp3err(cell,1,:,2),3),nanmean(condgp3err(cell,2,:,2),3);...
%         nanmean(condgp3err(cell,1,:,4),3),nanmean(condgp3err(cell,2,:,4),3)];
%     barweb(bars',errorbars',[],{'control','light on'},'gamma2 power','condition','power',[],[],{'cntr','co','is','cs','so'});
    
    
    % running figure
%     runlfr = mean(frs(intersect(find(result.light),oktrials)));
%     runnlfr = mean(frs(intersect(find(~result.light),oktrials)));
%     norunlfr = mean(frs(intersect(find(result.light),stilltrials)));
%     norunnlfr = mean(frs(intersect(find(~result.light),stilltrials)));
%     lfrerr = std(frs(find(result.light)))./sqrt(length(find(result.light)));
%     nlfrerr = std(frs(find(~result.light)))./sqrt(length(find(~result.light)));
%     runlfrerr = std(frs(intersect(find(result.light),oktrials)))./sqrt(length(intersect(find(result.light),oktrials)));
%     runnlfrerr = std(frs(intersect(find(~result.light),oktrials)))./sqrt(length(intersect(find(~result.light),oktrials)));
%     norunlfrerr = std(frs(intersect(find(result.light),stilltrials)))./sqrt(length(intersect(find(result.light),stilltrials)));
%     norunnlfrerr = std(frs(intersect(find(~result.light),stilltrials)))./sqrt(length(intersect(find(~result.light),stilltrials)));
%     
%     l1r1 = frs(intersect(find(result.light),oktrials));
%     l0r1 = frs(intersect(find(~result.light),oktrials));
%     l1r0 = frs(intersect(find(result.light),stilltrials));
%     l0r0 = frs(intersect(find(~result.light),stilltrials));
%     anovavec = [l0r0;l0r1;l1r0;l1r1]; 
%     g1 = [zeros(length(l0r0),1);zeros(length(l0r1),1);ones(length(l1r0),1);ones(length(l1r1),1)]; %light
%     g2 = [zeros(length(l0r0),1);ones(length(l0r1),1);zeros(length(l1r0),1);ones(length(l1r1),1)]; %running
%     [p,table,stats] = anovan(anovavec,{g1 g2},'model','full','display','off');
%     lp(cell) = p(1); rp(cell) = p(2); rlip(cell) = p(3);
%     
%     r0omi(cell) = (norunlfr-norunnlfr)/(norunlfr+norunnlfr);
%     r1omi(cell) = (runlfr-runnlfr)/(runlfr+runnlfr);
%     l0rmi(cell) = (runnlfr-norunnlfr)/(runnlfr+norunnlfr);
%     l1rmi(cell) = (runlfr-norunlfr)/(runlfr+norunlfr);
%     
%     figure
%     subplot(2,2,1)
%     imagesc(speed);
%     colorbar
%     title(['oktrials: ' int2str(length(oktrials)) '/' int2str(size(speed,1))])
%     xlabel('time [ms]')
%     ylabel('trial number')
%     
%     subplot(2,2,2)
%     errorbar(msta,mean(speed(find(~result.light),:)),std(speed(find(~result.light),:))./sqrt(length(find(~result.light))),'b')
%     hold on
%     errorbar(msta,mean(speed(find(result.light),:)),std(speed(find(result.light),:))./sqrt(length(find(result.light))),'r')
%     xlabel('time [ms]')
%     ylabel('average runspeed')
%     legend({'light off' 'light on'})
%     
%     subplot(2,2,3)
%     plot(mean(speed(:,respwin),2),frs,'.')
%     hold on
%     plot(mean(speed(find(result.light),respwin),2),frs(find(result.light)),'r.')
%     xlabel('average runspeed of trial')
%     ylabel('average firing rate of trial')
%     
%     subplot(2,2,4)
%     barweb([nlfr(cell),lfr(cell);runnlfr,runlfr;norunnlfr,norunlfr],...
%         [nlfrerr,lfrerr;runnlfrerr,runlfrerr;norunnlfrerr,norunlfrerr],...
%         [],[{'all'};{'running only'};{'immobile only'}],['ANOVA factor running p: ' num2str(p(2))],...
%         [],'firing rate [Hz]',[],[]);
%     

    cell = cell + 1;
    disp('');
    
end

figure
plot(swidth,adiff,'k.')
xlabel('spike width')
ylabel('amplitude diff')

pfs = find(swidth<125);

fsi = kmeans([swidth',adiff'],2); %kmeans([eslope',ptr',swidth',adiff'],2);
if mean(swidth(find(fsi==1)))<mean(swidth(find(fsi==2)))  %1 is FS
    pfs = find(fsi==1); prs = find(fsi==2);
else
    pfs = find(fsi==2); prs = find(fsi==1);
end

%coherence
params.Fs = 1000; params.trialave = 1; params.err = [2 .05]; params.tapers = [5,9];
cell1 = 19; cell2 = 6; %151022
% cell1 = 4; cell2 = 19; %150909
% cell1 = 5; cell2 = 19; % 150825
% cell1 = 1; cell2 = 16; % 150831
% cell1 = 15; cell2 = 4; % 150902

wo = 60/(1000/2); bw = wo/35;
[b,a] = iirnotch(wo,bw);

for i = 1:size(celllfpresp,2)
    fl1(i,:) = filtfilt(b,a,squeeze(celllfpresp(cell1,i,:)));
    fl2(i,:) = filtfilt(b,a,squeeze(celllfpresp(cell2,i,:)));
%     [coh(i,:),cfxx] = mscohere(squeeze(celllfpresp(cell1,i,respwin)),squeeze(celllfpresp(cell2,i,respwin)),[],[],512,1000);
%     [cpsdy(i,:),cpsdf] = cpsd(squeeze(celllfpresp(cell1,i,respwin)),squeeze(celllfpresp(cell2,i,respwin)),[],[],512,1000);
end
for l = 1:2
    for d = 1:length(ordf)
        if ordf(d) == -1
            thisinds = find(result.gratingInfo.Orientation_center ~= -1 & ...
                result.gratingInfo.Orientation_surround == -1 & ...
                result.light == l-1);
        elseif ordf(d) == -2
            thisinds = find(result.gratingInfo.Orientation_center == -1 & ...
                result.gratingInfo.Orientation_surround ~= -1 & ...
                result.light == l-1);
        elseif ordf(d) == 0
            thisinds = find(result.gratingInfo.Orientation_center == ...
                result.gratingInfo.Orientation_surround & ...
                result.gratingInfo.Orientation_center ~= -1 &...
                result.light == l-1);
        else
            thisinds = find(abs(result.gratingInfo.Orientation_center - ...
                result.gratingInfo.Orientation_surround) == 90 & ...
                result.light == l-1);
        end
        allcondlfpspect(l,d,:,:) = lfpspect(thisinds,:);
%         condcoher(l,d,:) = nanmean(coh(thisinds,:),1);
%         condcoherr(l,d,:) = nanstd(coh(thisinds,:),1)./sqrt(length(thisinds));
        if ~isempty(thisinds)
            [C(l,d,:),phi(l,d,:),S12(l,d,:),S1(l,d,:),S2(l,d,:),...
                cfx,confC(l,d),phistd(l,d,:),Cerr(l,d,:,:)] = coherencyc(squeeze(fl1(thisinds,respwin))',squeeze(fl2(thisinds,respwin))',params);
            
%             for i = 1:length(thisinds)
%                 tbtC(l,d,i,:) = coherencyc(squeeze(celllfpresp(cell1,thisinds(i),respwin))',squeeze(celllfpresp(cell2,thisinds(i),respwin))',params);
%             end
            % this is wrong somehow
            %                 filltbt(l,d,:) = [squeeze(nanmean(tbtC(l,d,:,1:120),6))'+(squeeze(nanstd(tbtC(l,d,:,1:120),1,6))./sqrt(length(thisinds)))',squeeze(nanmean(tbtC(l,d,:,1:120),6))'-(squeeze(nanstd(tbtC(l,d,:,1:120),1,6))./sqrt(length(thisinds)))'];
        else
        end
        fillyc(l,d,:) = [squeeze(Cerr(l,d,1,1:124))',fliplr(squeeze(Cerr(l,d,2,1:124))')];
        fillyspect(l,d,:) = [squeeze(condSerr(cell1,l,d,1,1:124))',fliplr(squeeze(condSerr(cell1,l,d,2,1:124))')];
        
    end
end
contcoher(1,:) = nanmean(coh(contindsnl,:),1);
contcoher(2,:) = nanmean(coh(contindsl,:),1);
fillx = [cfx(1:124),fliplr(cfx(1:124))];

oriratio = (squeeze(nanmean(cellfr(:,1,:,2),3))-squeeze(nanmean(cellfr(:,1,:,1),3)))./(squeeze(nanmean(cellfr(:,1,:,2),3))+squeeze(nanmean(cellfr(:,1,:,1),3)));
oriratiol1 = (squeeze(nanmean(cellfr(:,2,:,2),3))-squeeze(nanmean(cellfr(:,2,:,1),3)))./(squeeze(nanmean(cellfr(:,2,:,2),3))+squeeze(nanmean(cellfr(:,2,:,1),3)));

