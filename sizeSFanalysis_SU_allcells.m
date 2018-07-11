function sizeSFanalysis_SU_allcells

animalid = '160906';
block = 7;
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
    
    if strfind(files(fi).name, 'MU')
        continue;
    end
        
    load([supath, files(fi).name]);    
    
    prestim = 300;
    poststim = 700;
    if result.stimduration == 2
        respwin = 501:1500; % after stimulus onset
    else
        respwin = 1:1000;
    end
    respwin = respwin+prestim;
    
%     disp(['now analyzing file: ' files(cell).name]);
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
% %         msstamps([169,336]) = []; % for 150523 block 11
%         result.msstamps = msstamps;
%         save([supath, files(fi).name],'result');
        pause;
     end    
     
    % fix so it is usable with new multi-purpose grating stim script
    if isfield(result, 'sfconds')
        allinds = sort(getSpecificIndices(result, 'sfconds'));
        msstamps = result.msstamps(allinds);
        light = result.light(allinds);
        gratingInfo.Orientation = result.gratingInfo.Orientation(allinds);
        gratingInfo.size = result.gratingInfo.size(allinds);
        gratingInfo.Contrast = result.gratingInfo.Contrast(allinds);
        gratingInfo.tFreq = result.gratingInfo.tFreq(allinds);
        gratingInfo.spFreq = result.gratingInfo.spFreq(allinds);
    else
        msstamps = reuslt.msstamps;
        light = result.light;
        gratingInfo = result.gratingInfo;
    end
    
    
    for i = 1:length(msstamps)
        resp(i,:) = chan(msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim); 
        lfpresp(i,:) = result.lfp(msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim);
        [lfpspect(i,:),trialfax] = pmtm(lfpresp(i,respwin(1)+200:respwin(end)),3,[],sr);
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
    spfreqs = unique(gratingInfo.spFreq); spfreqs(spfreqs == 0) = [];
    for l = 1:2
        for ori = 1:length(oris)
            for sz = 1:length(sizes)
                for sf = 1:length(spfreqs)
                    thisinds = find(gratingInfo.Orientation == oris(ori) &...
                        gratingInfo.size == sizes(sz) & gratingInfo.spFreq == spfreqs(sf) & ...
                        light == l-1);
                    condresp(l,ori,sz,sf,:) = mean(resp(thisinds,:),1);
                    condlfpspect(l,ori,sz,sf,:) = nanmean(lfpspect(thisinds,:),1);
                    condstspect(l,ori,sz,sf,:) = pmtm(mean(resp(thisinds,1001:1800),1),3,[],sr);
                    condfr(l,ori,sz,sf) = mean(frs(thisinds));%-mean(bl);        
                    conderr(l,ori,sz,sf) =std(frs(thisinds))./sqrt(length(thisinds));

                    condz(l,ori,sz,sf) = {(sc(thisinds)-mean(sc(thisinds)))/std(sc(thisinds))}; %ecker 2010
                    condsc(l,ori,sz,sf) = {sc(thisinds)};
                    ff(l,ori,sz,sf) = var(sc(thisinds))/mean(sc(thisinds));

                    condlfpresp(l,ori,sz,sf,:) = mean(lfpresp(thisinds,:),1);                

                    condresperr(l,ori,sz,sf,:) = nanstd(resp(thisinds,:),1,1)./sqrt(length(thisinds));
                    if ~isnan(condresp(l,ori,sz,sf,:))
                        [bincondresp(l,ori,sz,sf,:),bta] = binit(condresp(l,ori,sz,sf,:),binwidth);
                    else
                        bincondresp(l,ori,sz,sf,:) = binit(condresp(l,ori,sz,sf,:),binwidth);
                    end
                    binconderr(l,ori,sz,sf,:) = binit(condresperr(l,ori,sz,sf,:),binwidth);

                    mscc = []; bincc = [];
                    for ii = 1:length(thisinds)-1
                        for jj = ii+1:length(thisinds)
                            help = corrcoef(resp(thisinds(ii),:),resp(thisinds(jj),:));
                            mscc = [mscc,help(1,2)];
                            help = corrcoef(binit(resp(thisinds(ii),:),binwidth),binit(resp(thisinds(jj),:),binwidth));
                            bincc = [bincc, help(1,2)];
                        end
                    end
                    msreliab(l,ori,sz,sf) = nanmean(mscc);
                    binreliab(l,ori,sz,sf) = nanmean(bincc);
                    eckerreliability(l,ori,sz,sf) = var(frs(thisinds))/var(frs);

                    thisruninds = intersect(thisinds,oktrials);
                    if ~isempty(thisruninds)
                        runcondresp(l,ori,sz,sf,:) = mean(resp(thisruninds,:),1);
                        runcondfr(l,ori,sz,sf) = mean(frs(thisruninds));
                        runconderr(l,ori,sz,sf) = std(frs(thisruninds))./sqrt(length(thisruninds));
                    else
                        runcondresp(l,ori,sz,sf,:) = nan(1,size(resp,2));
                        runcondfr(l,ori,sz,sf) = NaN;
                        runconderr(l,ori,sz,sf) = NaN;
                    end  

                    thisstillinds = intersect(thisinds,stilltrials);
                    if ~isempty(thisstillinds)
                        stillcondresp(l,ori,sz,sf,:) = mean(resp(thisstillinds,:),1);
                        stillcondfr(l,ori,sz,sf) = mean(frs(thisstillinds));
                        stillconderr(l,ori,sz,sf) = std(frs(thisstillinds))./sqrt(length(thisstillinds));
                    else
                        stillcondresp(l,ori,sz,sf,:) = nan(1,size(resp,2));
                        stillcondfr(l,ori,sz,sf) = NaN;
                        stillconderr(l,ori,sz,sf) = NaN;
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
    cellerr(cell,:,:,:,:) = conderr;
    celleckerrely(cell,:,:,:,:) = eckerreliability;
    cellmsrely(cell,:,:,:,:) = msreliab;
    cellbinrely(cell,:,:,:,:) = binreliab;
    celllfpspect(cell,:,:,:,:,:) = condlfpspect;
    
        
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
    errorbar(oris,squeeze(condfr(1,:,1,1)),squeeze(conderr(1,:,1,1)),'ko-','markerfacecolor','k','linewidth',2)
    hold on
    errorbar(oris,squeeze(condfr(1,:,1,2)),squeeze(conderr(1,:,1,2)),'kd--','markerfacecolor','k','linewidth',2)
    errorbar(oris,squeeze(condfr(1,:,1,3)),squeeze(conderr(1,:,1,3)),'ks:','markerfacecolor','k','linewidth',2)
    errorbar(oris,squeeze(condfr(1,:,2,1)),squeeze(conderr(1,:,2,1)),'go-','markerfacecolor','g','linewidth',2)
    errorbar(oris,squeeze(condfr(1,:,2,2)),squeeze(conderr(1,:,2,2)),'gd--','markerfacecolor','g','linewidth',2)
    errorbar(oris,squeeze(condfr(1,:,2,3)),squeeze(conderr(1,:,2,3)),'gs:','markerfacecolor','g','linewidth',2)
    errorbar(oris,squeeze(condfr(1,:,3,1)),squeeze(conderr(1,:,3,1)),'ro-','markerfacecolor','r','linewidth',2)
    errorbar(oris,squeeze(condfr(1,:,3,2)),squeeze(conderr(1,:,3,2)),'rd--','markerfacecolor','r','linewidth',2)
    errorbar(oris,squeeze(condfr(1,:,3,3)),squeeze(conderr(1,:,3,3)),'rs:','markerfacecolor','r','linewidth',2)
    xlabel('shown orientation')
    ylabel('Firing rate [Hz]')
    set(gca,'xtick',oris)
    legend({'sflow small','sfmed small','sfhigh small','sflow med','sfmed med','sfhigh med','sflow large','sfmed large','sfhigh large'})    
%     title(['sm OSI: ' num2str(nlsmosi(cell)) '  sm OSI Light: ' num2str(lsmosi(cell)) '  lg OSI: ' num2str(nllgosi(cell)) '  lg OSI Light: ' num2str(llgosi(cell))])
    
%     oneorifr = mean(reshape(condfr(:,:,1,maxc),2,4,2),3);
%     prefori = find(oneorifr(1,:) == max(oneorifr(1,:)),1);
%     ortho = mod(prefori+2,length(oris)/2); if ortho == 0, ortho = length(oris)/2; end
    
%     cellcondresppreforil0(cell,:,:) = squeeze(condresp(1,prefori,:,:));
%     cellcondresppreforil1(cell,:,:) = squeeze(condresp(2,prefori,:,:));
    
%     preffr(cell,:) = oneorifr(:,prefori);

    prefori = 1;
     
    subplot(2,2,3)
    errorbar(sizes,squeeze(nanmean(condfr(1,:,:,1),2)),...
        squeeze(nanmean(conderr(1,:,:,1),2)),'ko-','linewidth',2);
    hold on
    errorbar(sizes,squeeze(nanmean(condfr(1,:,:,2),2)),...
        squeeze(nanmean(conderr(1,:,:,2),2)),'go-','linewidth',2);
    errorbar(sizes,squeeze(nanmean(condfr(1,:,:,3),2)),...
        squeeze(nanmean(conderr(1,:,:,3),2)),'ro-','linewidth',2);
    xlabel('shown patch size [vd]')
    ylabel('Firing rate [Hz]')
    legend({'sflow','sfmed','sfhigh'})    
    set(gca,'xtick',sizes)  
    title(['cell ' int2str(cell)  ' depth: ' int2str(result.depth) '  ' printname])
     
    subplot(2,2,4)
    semilogy(trialfax,squeeze(nanmean(condlfpspect(1,:,1,1,:),2)),'color',[.7,.7,.7])
    hold on
    semilogy(trialfax,squeeze(nanmean(condlfpspect(1,:,1,2,:),2)),'color',[.3,.3,.3])
    semilogy(trialfax,squeeze(nanmean(condlfpspect(1,:,1,3,:),2)),'color',[0 0 0])
    semilogy(trialfax,squeeze(nanmean(condlfpspect(1,:,2,1,:),2)),'color',[.7,1,.7])
    semilogy(trialfax,squeeze(nanmean(condlfpspect(1,:,2,2,:),2)),'color',[.3,1,.3])
    semilogy(trialfax,squeeze(nanmean(condlfpspect(1,:,2,3,:),2)),'color',[0 1 0])
    semilogy(trialfax,squeeze(nanmean(condlfpspect(1,:,3,1,:),2)),'color',[1,.7,.7])
    hold on
    semilogy(trialfax,squeeze(nanmean(condlfpspect(1,:,3,2,:),2)),'color',[1,.3,.3])
    semilogy(trialfax,squeeze(nanmean(condlfpspect(1,:,3,3,:),2)),'color',[1 0 0])
    legend('small sflow','...','small sfhigh','med sflow','...','med sfhigh','large sflow','...','large sfhigh')
    axis([0,100,2,5000])

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