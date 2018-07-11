function sizecontanalysis_SU_allcells

% TODO fix for allsorts grating

animalid = '160725';
block = 5;
lcol = 'r'; %lasercolor

onlymod = 0;
printyn = 1;
sfc = 0;

supath = ['C:\Users\Julia\work\data\' animalid '\singleunits\'];
basename = [animalid '_block' int2str(block) '_tet'];

files = dir([supath, basename, '*.mat']);

freqbinwidth = 5;

cell = 1;
for fi = 1:length(files)
    
%     if strfind(files(fi).name, 'MU')
%         continue;
%     end
        
    load([supath, files(fi).name]);
%     disp(['now analyzing file: ' files(cell).name]);

    prestim = 300;
    poststim = 700;
    if result.stimduration == 2
        respwin = 501:1500; % after stimulus onset
    else
        respwin = 1:1000;
    end
    respwin = respwin+prestim;

    cellname{cell} = files(fi).name;
    
    i = strfind(files(fi).name, 'tet');
    tetno = strread(files(fi).name(i+3));
    
    wvchan = find(var(result.waveforms) == max(var(result.waveforms)));
    
    sr = 1000;
    lfp = result.lfp(:,wvchan)';
    nfft = 2^nextpow2(length(lfp));
    fax = sr/2*linspace(0,1,nfft/2+1);
    y = fft(lfp,nfft);
    lfpspectrum = abs(y(1:nfft/2+1));

%     plot(fax,lfpspectrum)
    gamma = eegfilt(lfp,sr,30,90);
    h = hilbert(gamma); gpow = abs(h); gphas = angle(h);
    % %
    if sfc
        for i = 1:100/freqbinwidth
            filtmat(i,:) = eegfilt(lfp,sr,(i-1)*freqbinwidth+1,i*freqbinwidth);
            h = hilbert(filtmat(i,:));
            powmat(i,:) = abs(h); phasmat(i,:) = angle(h);
        end
    end
        
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
%         msstamps([318]) = []; % for 150523 block 11
%         result.msstamps = msstamps;
%         save([supath, files(fi).name],'result');
        pause;
     end    
     
    if isfield(result, 'contconds')
        allinds = sort(getSpecificIndices(result, 'contconds'));
        msstamps = result.msstamps(allinds);
        light = result.light(allinds);
        gratingInfo.Orientation = result.gratingInfo.Orientation(allinds);
        gratingInfo.Contrast = result.gratingInfo.Contrast(allinds);
        gratingInfo.tFreq = result.gratingInfo.tFreq(allinds);
        gratingInfo.size = result.gratingInfo.size(allinds);
    else
        msstamps = result.msstamps;
        light = result.light;
        gratingInfo = result.gratingInfo;
    end    
    
    for i = 1:length(msstamps)
        resp(i,:) = chan(msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim); 
        lfpresp(i,:) = result.lfp(msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim);
        [lfpspect(i,:),trialfax] = pmtm(lfpresp(i,1001:1800),3,[],sr);
        gammaresp(i,:) = gpow(msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim);
        gphaseresp(i,:) = gphas(msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim);
        
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
    
    cellrespl0(cell,:,:) = resp(find(light == 0),:);
    cellrespl1(cell,:,:) = resp(find(light == 1),:);
    cellresp(cell,:,:) = resp;
    
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
    
    lightresp = resp(find(light),:);
    nolightresp = resp(find(light == 0),:);
    
    lightlfpresp = lfpresp(find(light),:);
    nolightlfpresp = lfpresp(find(~light),:);
    
    frs = sum(resp(:,respwin),2)./(length(respwin)/1000);
    bl = sum(resp(:,1:prestim),2)./(prestim/1000);
    sc = sum(resp(:,respwin),2);
        
     %determine if cell is visually modulated
    blfr = sum(resp(:,1:prestim),2);
    vrfr = sum(resp(:,prestim+40:2*prestim+40),2);
    vismod(cell) = ttest2(blfr,vrfr);
    
    %determine if cell is modulated by light
    lightmod(cell) = ttest2(frs(find(light)),frs(find(light == 0)));
    
    lfr(cell) = mean(frs(find(light)));
    nlfr(cell) = mean(frs(find(light == 0)));
    
    % phases
    tmp = zeros(size(gphaseresp));
    tmp(find(resp)) = gphaseresp(find(resp));
    l0phasemat = tmp(find(light == 0),:);
    l1phasemat = tmp(find(light == 1),:);
    l0phases{cell} = l0phasemat(find(l0phasemat));
    l1phases{cell} = l1phasemat(find(l1phasemat));
      
    rl0(cell) = circ_r(l0phases{cell});
    rl1(cell) = circ_r(l1phases{cell});
    cmeanl0(cell) = circ_mean(l0phases{cell});
    cmeanl1(cell) = circ_mean(l1phases{cell});
    
    if sfc
        tmpi = zeros(size(allphaseresp));
        for i = 1:size(allphaseresp,1)
            tmp = zeros(size(gphaseresp));
            tmp(find(resp)) = allphaseresp(i,find(resp));
            tmpi(i,:,:) = tmp;
        end
        alll0phasemat = tmpi(:,find(light == 0),:);
        alll1phasemat = tmpi(:,find(light == 1),:);
        for i = 1:size(allphaseresp,1)
            allphasesl0{i} = alll0phasemat(i,find(squeeze(alll0phasemat(i,:,:)))); 
            allphasesl1{i} = alll1phasemat(i,find(squeeze(alll1phasemat(i,:,:)))); 
            allrl0(cell,i) = circ_r(allphasesl0{i}');
            allrl1(cell,i) = circ_r(allphasesl1{i}');
            allcmeanl0(cell,i) = circ_mean(allphasesl0{i}');
            allcmeanl1(cell,i) = circ_mean(allphasesl1{i}');
        end
    end
    
    
    binwidth = 30;
    [binnedlight,bta] = binit(mean(lightresp),binwidth);
    binnedlight = binnedlight.*(1000/binwidth);
    [binnednolight,bta] = binit(mean(nolightresp),binwidth);
    binnednolight = binnednolight.*(1000/binwidth);
    
    printname = files(fi).name;
    printname(find(printname=='_')) = ' ';
    
    sizes = unique(gratingInfo.size);  sizes(sizes == 0) = []; %!!!!!!!!!! for the one
    oris = unique(gratingInfo.Orientation); oris(oris == -1) = [];
    contrasts = unique(gratingInfo.Contrast); contrasts(find(contrasts == 0)) = [];
    for l = 1:2
        for ori = 1:length(oris)
            for sz = 1:length(sizes)
                for co = 1:length(contrasts)
                    thisinds = find(gratingInfo.Orientation == oris(ori) &...
                        gratingInfo.size == sizes(sz) & gratingInfo.Contrast == contrasts(co) & ...
                        light == l-1);
                    condresp(l,ori,sz,co,:) = mean(resp(thisinds,:),1);
                    condlfpspect(l,ori,sz,co,:) = nanmean(lfpspect(thisinds,:),1);
                    condstspect(l,ori,sz,co,:) = pmtm(mean(resp(thisinds,1001:1800),1),3,[],sr);
                    condfr(l,ori,sz,co) = mean(frs(thisinds));%-mean(bl);        
                    conderr(l,ori,sz,co) =std(frs(thisinds))./sqrt(length(thisinds));
                    
                    condtbtfrvar(cell,l,ori,sz,co) = var(frs(thisinds));
                    condtbtlfpvar(cell,l,ori,sz,co) = nanmean(var(lfpresp(thisinds,respwin),1),2);
                    
                    condz(l,ori,sz,co) = {(sc(thisinds)-mean(sc(thisinds)))/std(sc(thisinds))}; %ecker 2010
                    condsc(l,ori,sz,co) = {sc(thisinds)};
                    ff(l,ori,sz,co) = var(sc(thisinds))/mean(sc(thisinds));

                    condlfpresp(l,ori,sz,co,:) = mean(lfpresp(thisinds,:),1);                

                    condresperr(l,ori,sz,co,:) = nanstd(resp(thisinds,:),1,1)./sqrt(length(thisinds));
                    if ~isnan(condresp(l,ori,sz,co,:))
                        [bincondresp(l,ori,sz,co,:),bta] = binit(condresp(l,ori,sz,co,:),binwidth);
                    else
                        bincondresp(l,ori,sz,co,:) = binit(condresp(l,ori,sz,co,:),binwidth);
                    end
                    binconderr(l,ori,sz,co,:) = binit(condresperr(l,ori,sz,co,:),binwidth);

                    mscc = []; bincc = [];
                    for ii = 1:length(thisinds)-1
                        for jj = ii+1:length(thisinds)
                            help = corrcoef(resp(thisinds(ii),:),resp(thisinds(jj),:));
                            mscc = [mscc,help(1,2)];
                            help = corrcoef(binit(resp(thisinds(ii),:),binwidth),binit(resp(thisinds(jj),:),binwidth));
                            bincc = [bincc, help(1,2)];
                        end
                    end
                    msreliab(l,ori,sz,co) = nanmean(mscc);
                    binreliab(l,ori,sz,co) = nanmean(bincc);
                    eckerreliability(l,ori,sz,co) = var(frs(thisinds))/var(frs);

                    thisruninds = intersect(thisinds,oktrials);
                    if ~isempty(thisruninds)
                        runcondresp(l,ori,sz,co,:) = mean(resp(thisruninds,:),1);
                        runcondfr(l,ori,sz,co) = mean(frs(thisruninds));
                        runconderr(l,ori,sz,co) = std(frs(thisruninds))./sqrt(length(thisruninds));
                    else
                        runcondresp(l,ori,sz,co,:) = nan(1,size(resp,2));
                        runcondfr(l,ori,sz,co) = NaN;
                        runconderr(l,ori,sz,co) = NaN;
                    end  

                    thisstillinds = intersect(thisinds,stilltrials);
                    if ~isempty(thisstillinds)
                        stillcondresp(l,ori,sz,co,:) = mean(resp(thisstillinds,:),1);
                        stillcondfr(l,ori,sz,co) = mean(frs(thisstillinds));
                        stillconderr(l,ori,sz,co) = std(frs(thisstillinds))./sqrt(length(thisstillinds));
                    else
                        stillcondresp(l,ori,sz,co,:) = nan(1,size(resp,2));
                        stillcondfr(l,ori,sz,co) = NaN;
                        stillconderr(l,ori,sz,co) = NaN;
                    end  
                end                
            end
        end
    end
    
    bincondresp = bincondresp.*(1000/binwidth);
    bta = bta-prestim;
%     psthcondplot(bincondresp,binconderr,bta)
    cellbinresp(cell,:,:,:,:,:) = bincondresp;

    binnedcellrespl0(cell,:) = binnednolight;
    binnedcellrespl1(cell,:) = binnedlight;
    
    nolmaxfr(cell) = max(max(max(condfr(1,:,:,:))));
    lmaxfr(cell) = max(max(max(condfr(2,:,:,:))));
    
    cellz(cell,:,:,:,:) = condz;
    cellsc(cell,:,:,:,:) = condsc;
    cellff(cell,:,:,:,:) = ff;
    cellfr(cell,:,:,:,:) = condfr;
    celleckerrely(cell,:,:,:,:) = eckerreliability;
    cellmsrely(cell,:,:,:,:) = msreliab;
    cellbinrely(cell,:,:,:,:) = binreliab;
    celllfpspect(cell,:,:,:,:,:) = condlfpspect;
    
    contindsnl = find(gratingInfo.size == 0 & gratingInfo.Contrast == 0 & light == 0);
    contindsl = find(gratingInfo.size == 0 & gratingInfo.Contrast == 0 & light == 1);
    
    controlresp(1,:) = mean(resp(contindsnl,:),1);
    controlresp(2,:) = mean(resp(contindsl,:),1);
    controllfpspect(1,:) = nanmean(lfpspect(contindsnl,:),1);
    controllfpspect(2,:) = nanmean(lfpspect(contindsl,:),1);
    controlfr(1) = mean(frs(contindsnl));%-mean(bl);
    controlfr(2) = mean(frs(contindsl));%-mean(bl);
    controlerr(1) =std(frs(contindsnl))./sqrt(length(contindsnl));
    controlerr(2) =std(frs(contindsl))./sqrt(length(contindsl));
    
    maxc = length(contrasts);
        
%     [nlsmoriprefratio(cell), nlsmdirprefratio(cell), nlsmprefori, nlsmmeanori, nlsmosi(cell), nlsmmeandir, nlsmdsi(cell)] = getOSI(squeeze(condfr(1,:,1,maxc)),oris);
%     [lsmoriprefratio(cell), lsmdirprefratio(cell), lsmprefori, lsmmeanori, lsmosi(cell), lsmmeandir, lsmdsi(cell)] = getOSI(squeeze(condfr(2,:,1,maxc)),oris);
%     [nllgoriprefratio(cell), nllgdirprefratio(cell), nllgprefori, nllgmeanori, nllgosi(cell), nllgmeandir, nllgdsi(cell)] = getOSI(squeeze(condfr(1,:,2,maxc)),oris);
%     [llgoriprefratio(cell), llgdirprefratio(cell), llgprefori, llgmeanori, llgosi(cell), llgmeandir, llgdsi(cell)] = getOSI(squeeze(condfr(2,:,2,maxc)),oris);
  
    figure
    subplot(2,2,1)
    ta = bta;%-prestim;
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
    
    subplot(2,2,2)
    errorbar(oris,squeeze(condfr(2,:,1,maxc)),squeeze(conderr(2,:,1,maxc)),'o-','color',lcol,'markersize',8,'linewidth',2)
    hold on
    errorbar(oris,squeeze(condfr(1,:,1,maxc)),squeeze(conderr(1,:,1,maxc)),'bo-','markersize',8,'linewidth',2)
    errorbar(oris,squeeze(condfr(2,:,2,maxc)),squeeze(conderr(2,:,2,maxc)),'mo-','color','m','markersize',8,'linewidth',2)
    errorbar(oris,squeeze(condfr(1,:,2,maxc)),squeeze(conderr(1,:,2,maxc)),'co-','markersize',8,'linewidth',2)
    xlabel('shown orientation')
    ylabel('Firing rate [Hz]')
    set(gca,'xtick',oris)
    legend({'Light ON small','Light OFF small','Light ON large','Light OFF large'})    
%     title(['sm OSI: ' num2str(nlsmosi(cell)) '  sm OSI Light: ' num2str(lsmosi(cell)) '  lg OSI: ' num2str(nllgosi(cell)) '  lg OSI Light: ' num2str(llgosi(cell))])
    
%     oneorifr = mean(reshape(condfr(:,:,1,maxc),2,4,2),3);
%     prefori = find(oneorifr(1,:) == max(oneorifr(1,:)),1);
%     ortho = mod(prefori+2,length(oris)/2); if ortho == 0, ortho = length(oris)/2; end
    
%     cellcondresppreforil0(cell,:,:) = squeeze(condresp(1,prefori,:,:));
%     cellcondresppreforil1(cell,:,:) = squeeze(condresp(2,prefori,:,:));
    
%     preffr(cell,:) = oneorifr(:,prefori);

    prefori = 1;
     
    subplot(2,2,3)
    errorbar(contrasts,squeeze(nanmean(condfr(2,[prefori,prefori+(length(oris)/2)],1,:))),...
        squeeze(nanmean(conderr(2,[prefori,prefori+(length(oris)/2)],1,:))),'o-','color',lcol,'markersize',8,'linewidth',2);
    hold on
    errorbar(contrasts,squeeze(nanmean(condfr(1,[prefori,prefori+(length(oris)/2)],1,:))),...
        squeeze(nanmean(conderr(1,[prefori,prefori+(length(oris)/2)],1,:))),'bo-','markersize',8,'linewidth',2);
    errorbar(contrasts,squeeze(nanmean(condfr(2,[prefori,prefori+(length(oris)/2)],2,:))),...
        squeeze(nanmean(conderr(2,[prefori,prefori+(length(oris)/2)],2,:))),'o-','color','m','markersize',8,'linewidth',2);
    errorbar(contrasts,squeeze(nanmean(condfr(1,[prefori,prefori+(length(oris)/2)],2,:))),...
        squeeze(nanmean(conderr(1,[prefori,prefori+(length(oris)/2)],2,:))),'o-','color','c','markersize',8,'linewidth',2);
    xlabel('shown patch size [vd]')
    ylabel('Firing rate [Hz]')
    legend({'Light ON small','Light OFF small','Light ON large','Light OFF large'})    
    set(gca,'xtick',contrasts)  
    title(['cell ' int2str(cell) ' preferred orientations' ' depth: ' int2str(result.depth) '  ' printname])
     
    % alternativ LFP subplot
%     subplot(2,2,4)
%     semilogy(trialfax,squeeze(nanmean(condlfpspect(1,:,1,1,:),2)),'color',[0.8,0.8,1])
%     hold on
%     semilogy(trialfax,squeeze(nanmean(condlfpspect(1,:,1,2,:),2)),'color',[0.6,0.6,1])
%     semilogy(trialfax,squeeze(nanmean(condlfpspect(1,:,1,3,:),2)),'color',[0.4,0.4,1])
%     semilogy(trialfax,squeeze(nanmean(condlfpspect(1,:,1,4,:),2)),'color',[0.2,0.2,1])
% %     semilogy(trialfax,squeeze(nanmean(condlfpspect(1,:,1,5,:),2)),'color',[0,0,1])
%     semilogy(trialfax,squeeze(nanmean(condlfpspect(1,:,2,1,:),2)),'color',[1,0.8,0.8])
%     semilogy(trialfax,squeeze(nanmean(condlfpspect(1,:,2,2,:),2)),'color',[1,0.6,0.6])
%     semilogy(trialfax,squeeze(nanmean(condlfpspect(1,:,2,3,:),2)),'color',[1,0.4,0.4])
%     semilogy(trialfax,squeeze(nanmean(condlfpspect(1,:,2,4,:),2)),'color',[1,0.2,0.2])
% %     semilogy(trialfax,squeeze(nanmean(condlfpspect(1,:,2,5,:),2)),'color',[1,0,0])
%     legend('small lowest','...','...','small full','large lowest','...','...','large full')
%     axis([0,100,2,5000])

    subplot(2,2,4)
    errorbar(sizes,squeeze(nanmean(condfr(1,:,:,1),2)),squeeze(nanmean(conderr(1,:,:,1),2)),'b')
    hold on
    errorbar(sizes,squeeze(nanmean(condfr(1,:,:,3),2)),squeeze(nanmean(conderr(1,:,:,3),2)),'g')
    errorbar(sizes,squeeze(nanmean(condfr(1,:,:,length(contrasts)),2)),squeeze(nanmean(conderr(1,:,:,length(contrasts)),2)),'r')
    xlabel('shown size')
    ylabel('Firing rate [Hz]')
    legend({'low contrast','medium contrast','high contrast'})    
    set(gca,'xtick',sizes)  
%     title(['cell ' int2str(cell) ' preferred orientations' ' depth: ' int2str(result.depth) '  ' printname])

    o = 2; s = 3;
%     condl0 = (gratingInfo.Orientation == oris(o) | gratingInfo.Orientation == oris(o+4)) & (gratingInfo.size == sizes(s) | gratingInfo.size == sizes(s+1)) & light == 0;
%     condl1 = (gratingInfo.Orientation == oris(o) | gratingInfo.Orientation == oris(o+4)) & (gratingInfo.size == sizes(s) | gratingInfo.size == sizes(s+1)) & light == 1;
%     condl0 = gratingInfo.Orientation == oris(o) & gratingInfo.size == sizes(s) & light == 0;
%     condl1 = gratingInfo.Orientation == oris(o) & gratingInfo.size == sizes(s) & light == 1;
    condl0 = light == 0; condl1 = light == 1;
%     figure, rasterplot(resp,condl0,condl1,msta);
%     line([0,2000],[37,37],'color','k','linewidth',2)
%     line([500,1500],[33,33],'color','r','linewidth',2)
    
    % running figure
    runlfr = mean(frs(intersect(find(light),oktrials)));
    runnlfr = mean(frs(intersect(find(~light),oktrials)));
    norunlfr = mean(frs(intersect(find(light),stilltrials)));
    norunnlfr = mean(frs(intersect(find(~light),stilltrials)));
    lfrerr = std(frs(find(light)))./sqrt(length(find(light)));
    nlfrerr = std(frs(find(~light)))./sqrt(length(find(~light)));
    runlfrerr = std(frs(intersect(find(light),oktrials)))./sqrt(length(intersect(find(light),oktrials)));
    runnlfrerr = std(frs(intersect(find(~light),oktrials)))./sqrt(length(intersect(find(~light),oktrials)));
    norunlfrerr = std(frs(intersect(find(light),stilltrials)))./sqrt(length(intersect(find(light),stilltrials)));
    norunnlfrerr = std(frs(intersect(find(~light),stilltrials)))./sqrt(length(intersect(find(~light),stilltrials)));
    
    l1r1 = frs(intersect(find(light),oktrials));
    l0r1 = frs(intersect(find(~light),oktrials));
    l1r0 = frs(intersect(find(light),stilltrials));
    l0r0 = frs(intersect(find(~light),stilltrials));
    anovavec = [l0r0;l0r1;l1r0;l1r1]; 
    g1 = [zeros(length(l0r0),1);zeros(length(l0r1),1);ones(length(l1r0),1);ones(length(l1r1),1)]; %light
    g2 = [zeros(length(l0r0),1);ones(length(l0r1),1);zeros(length(l1r0),1);ones(length(l1r1),1)]; %running
    [p,table,stats] = anovan(anovavec,{g1 g2},'model','full','display','off');
    lp(cell) = p(1); rp(cell) = p(2); rlip(cell) = p(3);
    
    r0omi(cell) = (norunlfr-norunnlfr)/(norunlfr+norunnlfr);
    r1omi(cell) = (runlfr-runnlfr)/(runlfr+runnlfr);
    l0rmi(cell) = (runnlfr-norunnlfr)/(runnlfr+norunnlfr);
    l1rmi(cell) = (runlfr-norunlfr)/(runlfr+norunlfr);
    
%     figure
%     subplot(2,2,1)
%     imagesc(speed);
%     colorbar
%     title(['oktrials: ' int2str(length(oktrials)) '/' int2str(size(speed,1))])
%     xlabel('time [ms]')
%     ylabel('trial number')
%     
%     subplot(2,2,2)
%     errorbar(msta,mean(speed(find(~light),:)),std(speed(find(~light),:))./sqrt(length(find(~light))),'b')
%     hold on
%     errorbar(msta,mean(speed(find(light),:)),std(speed(find(light),:))./sqrt(length(find(light))),'r')
%     xlabel('time [ms]')
%     ylabel('average runspeed')
%     legend({'light off' 'light on'})
%     
%     subplot(2,2,3)
%     plot(mean(speed(:,respwin),2),frs,'.')
%     hold on
%     plot(mean(speed(find(light),respwin),2),frs(find(light)),'r.')
%     xlabel('average runspeed of trial')
%     ylabel('average firing rate of trial')
%     
%     subplot(2,2,4)
%     barweb([nlfr(cell),lfr(cell);runnlfr,runlfr;norunnlfr,norunlfr],...
%         [nlfrerr,lfrerr;runnlfrerr,runlfrerr;norunnlfrerr,norunlfrerr],...
%         [],[{'all'};{'running only'};{'immobile only'}],['ANOVA factor running p: ' num2str(p(2))],...
%         [],'firing rate [Hz]',[],[]);
    
%     %running with LFP
%     for i = 1:20
%         a = squeeze(mean(allpowresp(i,:,respwin),3));
%         b = mean(speed(:,respwin),2);
%         [r,p] = corrcoef(a,b);
%         rr(i) = r(1,2);
%         pp(i) = p(1,2);
%     end
%     fax = 3:5:100
%     
%     figure
%     plot(fax,rr,'o','markerfacecolor','b')
%     hold on
%     ss = find(pp<.05)
%     plot(fax(ss),rr(ss),'ro','markerfacecolor','r')

    figure
    errorbar([0,contrasts],[controlfr(1),squeeze(nanmean(condfr(1,:,1,:),2))'],[controlerr(1),squeeze(nanmean(conderr(1,:,1,:),2))'],'b','linewidth',2)
    hold on    
    errorbar([0,contrasts],[controlfr(1),squeeze(nanmean(condfr(1,:,2,:),2))'],[controlerr(1),squeeze(nanmean(conderr(1,:,2,:),2))'],'g','linewidth',2)
    errorbar([0,contrasts],[controlfr(1),squeeze(nanmean(condfr(1,:,3,:),2))'],[controlerr(1),squeeze(nanmean(conderr(1,:,3,:),2))'],'r','linewidth',2)
    legend(int2str(sizes'))
    
    figure
    errorbar([0,sizes],[controlfr(1),squeeze(nanmean(condfr(1,:,:,1),2))'],[controlerr(1),squeeze(nanmean(conderr(1,:,:,1),2))'],'b','linewidth',2)
    hold on    
    errorbar([0,sizes],[controlfr(1),squeeze(nanmean(condfr(1,:,:,2),2))'],[controlerr(1),squeeze(nanmean(conderr(1,:,:,2),2))'],'c','linewidth',2)
    errorbar([0,sizes],[controlfr(1),squeeze(nanmean(condfr(1,:,:,3),2))'],[controlerr(1),squeeze(nanmean(conderr(1,:,:,3),2))'],'y','linewidth',2)
    errorbar([0,sizes],[controlfr(1),squeeze(nanmean(condfr(1,:,:,4),2))'],[controlerr(1),squeeze(nanmean(conderr(1,:,:,4),2))'],'m','linewidth',2)
    errorbar([0,sizes],[controlfr(1),squeeze(nanmean(condfr(1,:,:,5),2))'],[controlerr(1),squeeze(nanmean(conderr(1,:,:,5),2))'],'r','linewidth',2)
    legend(num2str(contrasts'))
    
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

disp('');