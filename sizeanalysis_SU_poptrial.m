function sizeanalysis_SU_poptrial

animalid = '150909';
block = 4;
lcol = 'r'; %lasercolor
ntetrodes = 8;

onlymod = 0;
printyn = 1;
sfc = 0;

supath = ['C:\Users\Julia\work\data\' animalid '\singleunits\'];
basename = [animalid '_block' int2str(block) '_tet'];

files = dir([supath, basename, '*.mat']);

prestim = 300;
poststim = 700;
respwin = 501:1500; % after stimulus onset
respwin = respwin+prestim;
freqbinwidth = 5;

% prestim = 0;
% poststim = 0;
% respwin = 1500:2250;
% respwin = respwin-prestim;
% blwin = 1:1000;

tetdone = zeros(1,ntetrodes);
cell = 1;
for fi = 1:length(files)
    
    if strfind(files(fi).name, 'MU')
        continue;
    end
        
    load([supath, files(fi).name]);
    cellname{cell} = files(fi).name;
    
    % timestamps
    trialdur = result.stimduration*1000;
%     trialdur = result.sweeplength;

    msstamps = result.msstamps;    
     if length(msstamps)~=length(result.light)
% %         msstamps([169,336]) = []; % for 150523 block 11
%         result.msstamps = msstamps;
%         save([supath, files(fi).name],'result');
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

    % LFP if not done yet
    i = strfind(files(fi).name, 'tet');
    tetno = strread(files(fi).name(i+3));
    
    if ~tetdone(tetno)
        for i = 1:4
            lfp((tetno-1)*4+i,:) = result.lfp(:,i);
        end        
        tetdone(tetno) = 1;        
    end
    
    wvchan(cell) = find(var(result.waveforms) == max(var(result.waveforms))); 
  
    msStimes = round(result.spikes);
    if isempty(msStimes), msStimes(1) = 0; end
    if msStimes(1) == 0, msStimes(1) = 1; end  
    channel(cell,:) = zeros(1,length(result.lfp));
    channel(cell,msStimes) = 1;
    
    binwidth = 33.3;
    for i = 1:length(msstamps)
        resp(cell,i,:) = channel(cell,msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim);
        [binresp(cell,i,:),bta] = binit(squeeze(resp(cell,i,:)),binwidth); 
        [acorr(cell,i,:),lags] = xcorr(squeeze(resp(cell,i,respwin)),squeeze(resp(cell,i,respwin)));
    end
    
    depth(cell) = result.depth;
    
    spike = result.waveforms(:,wvchan(cell));
    interpspike = spline(1:32,spike,1:.1:32);
    [adiff(cell),swidth(cell)] = spikequant(interpspike);
    prs(cell) = swidth(cell)>120;
    
    cellrespl0(cell,:,:) = resp(cell,find(light == 0),:);
    cellrespl1(cell,:,:) = resp(cell,find(light == 1),:);
    
    cell = cell+1;
end
    
[mn,ch] = min(abs(depth-200)); % find channel closest to 200mum
chan = wvchan(ch);

% first find gamma peaks for this channel
beta = [15,40];
gamma = [50,70];
large = find(gratingInfo.size == max(unique(gratingInfo.size)) & light == 0);
small = find(gratingInfo.size == min(unique(gratingInfo.size(gratingInfo.size~=0))) & light == 0);
for i = 1:length(large)
    [pl(i,:),f] = pmtm(lfp(chan,msstamps(large(i)):msstamps(large(i))+1000),3,[],1000);
    [ps(i,:),f] = pmtm(lfp(chan,msstamps(small(i)):msstamps(small(i))+1000),3,[],1000);
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
cfg1 = f(bpi); cfg2 = f(gpi); % center frequencies

sr = 1000;
gamma1 = eegfilt(lfp(chan,:),sr,f(bpi)-2.5,f(bpi)+2.5);
gamma2 = eegfilt(lfp(chan,:),sr,f(gpi)-2.5,f(gpi)+2.5);
h1 = hilbert(gamma1); gpow1 = abs(h1); gphas1 = angle(h1);
h2 = hilbert(gamma2); gpow2 = abs(h2); gphas2 = angle(h2);
for cell = 1:size(resp,1)
    g1phases{cell} = gphas1(find(channel(cell,:)));
    g1reighlyp(cell) = circ_rtest(g1phases{cell});
    g1r(cell) = circ_r(g1phases{cell}');
    g1cmean(cell) = circ_mean(g1phases{cell}');
    nspikes(cell) = length(find(channel(cell,:)));
%     g1ppc(cell) = ppc(g1phases{cell}');
    g2phases{cell} = gphas2(find(channel(cell,:)));
    g2reighlyp(cell) = circ_rtest(g2phases{cell});
    g2r(cell) = circ_r(g2phases{cell}');
    g2cmean(cell) = circ_mean(g2phases{cell}');
%     g2ppc(cell) = ppc(g2phases{cell}');
end

% cut trials, too
for i = 1:length(msstamps)
    lfpresp(i,:) = lfp(chan,msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim);
    [lfpspect(i,:),trialfax] = pmtm(squeeze(lfpresp(i,1001:1800)),3,[],sr);
    gamma1resp(i,:) = gamma1(msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim);
    gamma1powresp(i,:) = gpow1(msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim);
    g1phaseresp(i,:) = gphas1(msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim);
    gamma2resp(i,:) = gamma2(msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim);
    gamma2powresp(i,:) = gpow2(msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim);
    g2phaseresp(i,:) = gphas2(msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim);
    for cell = 1:size(resp,1)        
        g1sphases{cell,i} = g1phaseresp(i,find(resp(cell,i,:)));
        g2sphases{cell,i} = g2phaseresp(i,find(resp(cell,i,:)));
    end    
    speed(i,:) = result.runspeed(msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim);
end

% overview figure
for tr = find(result.light == 0)
    clf
    b = g1phaseresp(tr,801:1800);
    plot(gamma1resp(tr,801:1800)./20+7);
    hold on
%     imagesc(repmat(g1phaseresp(tr,801:1800),size(resp,1),1));
%     colormap hsv
    for cell = 1:size(resp,1)
        a = find(squeeze(resp(cell,tr,801:1800)));
        for i = 1:length(a)
            line([a(i),a(i)],[cell-.5,cell+.5],'color','k','linewidth',2)
        end
        if prs(cell), cellstr = ' RS  '; else cellstr = ' FS  '; end
        trppc(tr,cell) = ppc(b(a));
        text(1050,cell,[cellstr int2str(depth(cell)), '  ppc: ' num2str(trppc(tr,cell))])
        title(['trial ' int2str(tr) ' light: ' int2str(result.light(tr)) ' size: ' int2str(result.gratingInfo.size(tr)) ' ori: ' int2str(result.gratingInfo.Orientation(tr))])
    end
    
    pause
end

% figure out sufficiently high and nonvariable runspeed trials
meanspeed = mean(speed(:,respwin),2);
stdspeed = std(speed(:,respwin),1,2);
notstill = find(meanspeed>1);
okspeed = find(meanspeed>( mean(meanspeed(notstill))-(1.5*std(meanspeed(notstill))) ) );
okvar = find(stdspeed<( mean(stdspeed(notstill))+(1.5*std(stdspeed(notstill)))) & stdspeed>.5);
oktrials = intersect(okspeed,okvar);
nonoktrials = 1:size(resp,2); nonoktrials(oktrials) = [];
stilltrials = 1:size(resp,2); stilltrials(notstill) = [];

msta = linspace(-prestim,trialdur+poststim,size(resp,32));

for cell = 1:size(resp,1)
    lightresp(cell,:) = nanmean(resp(cell,find(light),:),2);
    nolightresp(cell,:) = nanmean(resp(cell,find(light == 0),:),2);

    lightlfpresp = nanmean(lfpresp(find(light),:));
    nolightlfpresp = nanmean(lfpresp(find(~light),:));
    
    frs(cell,:) = sum(resp(cell,:,respwin),3)./(length(respwin)/1000);
    bl(cell,:) = sum(resp(cell,:,1:prestim),3)./(prestim/1000);
    sc(cell,:) = sum(resp(cell,:,respwin),3);
        
     %determine if cell is visually modulated
    blfr(cell,:) = sum(resp(cell,:,1:prestim),3);
    vrfr(cell,:) = sum(resp(cell,:,prestim+40:2*prestim+40),3);
    vismod(cell) = ttest2(blfr(cell,:),vrfr(cell,:));
    
    %determine if cell is modulated by light
    lightmod(cell) = ttest2(frs(cell,find(light)),frs(cell,find(light == 0)));
    
    lfr(cell) = mean(frs(cell,find(light)));
    nlfr(cell) = mean(frs(cell,find(light == 0)));
end



sizes = unique(gratingInfo.size);  sizes(find(sizes == 0)) = []; %delete control condition
oris = unique(gratingInfo.Orientation); oris(find(oris == -1)) = [];
for l = 1:2
    for ori = 1:length(oris)
        for sz = 1:length(sizes)
            thisinds = find(gratingInfo.Orientation == oris(ori) &...
                gratingInfo.size == sizes(sz) & ...
                light == l-1);
            condresp(:,l,ori,sz,:) = nanmean(resp(:,thisinds,:),2); % TODO fic dimensions from here on
            condlfpspect(:,l,ori,sz,:) = nanmean(lfpspect(thisinds,:));
            condfr(:,l,ori,sz) = mean(frs(:,thisinds),2);%-mean(bl);
            conderr(:,l,ori,sz) =std(frs(:,thisinds),1,2)./sqrt(length(thisinds));
            
            condz(l,ori,sz) = {(sc(:,thisinds)-mean(sc(:,thisinds),2))/std(sc(:,thisinds),1,2)}; %ecker 2010
            condsc(l,ori,sz) = {sc(thisinds)};
            ff(l,ori,sz) = var(sc(thisinds))/mean(sc(thisinds));
            
            condlfpresp(l,ori,sz,:) = mean(lfpresp(thisinds,:),1);
            
            condresperr(l,ori,sz,:) = nanstd(resp(thisinds,:),1,1)./sqrt(length(thisinds));
            if ~isnan(condresp(l,ori,sz,:))
                [bincondresp(l,ori,sz,:),bta] = binit(condresp(l,ori,sz,:),binwidth);
            else
                bincondresp(l,ori,sz,:) = binit(condresp(l,ori,sz,:),binwidth);
            end
            binconderr(l,ori,sz,:) = binit(condresperr(l,ori,sz,:),binwidth);
            
            mscc = []; bincc = [];
            for ii = 1:length(thisinds)-1
                for jj = ii+1:length(thisinds)
                    help = corrcoef(resp(thisinds(ii),:),resp(thisinds(jj),:));
                    mscc = [mscc,help(1,2)];
                    help = corrcoef(binit(resp(thisinds(ii),:),binwidth),binit(resp(thisinds(jj),:),binwidth));
                    bincc = [bincc, help(1,2)];
                end
            end
            msreliab(l,ori,sz) = nanmean(mscc);
            binreliab(l,ori,sz) = nanmean(bincc);
            eckerreliability(l,ori,sz) = var(frs(thisinds))/var(frs);
            
            thisruninds = intersect(thisinds,oktrials);
            if ~isempty(thisruninds)
                runcondresp(l,ori,sz,:) = mean(resp(thisruninds,:),1);
                runcondfr(l,ori,sz) = mean(frs(thisruninds));
                runconderr(l,ori,sz) = std(frs(thisruninds))./sqrt(length(thisruninds));
            else
                runcondresp(l,ori,sz,:) = nan(1,size(resp,2));
                runcondfr(l,ori,sz) = NaN;
                runconderr(l,ori,sz) = NaN;
            end
            
            thisstillinds = intersect(thisinds,stilltrials);
            if ~isempty(thisstillinds)
                stillcondresp(l,ori,sz,:) = mean(resp(thisstillinds,:),1);
                stillcondfr(l,ori,sz) = mean(frs(thisstillinds));
                stillconderr(l,ori,sz) = std(frs(thisstillinds))./sqrt(length(thisstillinds));
            else
                stillcondresp(l,ori,sz,:) = nan(1,size(resp,2));
                stillcondfr(l,ori,sz) = NaN;
                stillconderr(l,ori,sz) = NaN;
            end
            
        end
    end
end


    
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
    
    sizes = unique(gratingInfo.size);  sizes(find(sizes == 0)) = []; %delete control condition
    oris = unique(gratingInfo.Orientation); oris(find(oris == -1)) = [];
    for l = 1:2
        for ori = 1:length(oris)
            for sz = 1:length(sizes)
                thisinds = find(gratingInfo.Orientation == oris(ori) &...
                    gratingInfo.size == sizes(sz) & ...
                    light == l-1);
                condresp(l,ori,sz,:) = mean(resp(thisinds,:),1);
                condlfpspect(l,ori,sz,:) = nanmean(lfpspect(thisinds,:));
                condstspect(l,ori,sz,:) = pmtm(mean(resp(thisinds,1001:1800),1),3,[],sr);
                condfr(l,ori,sz) = mean(frs(thisinds));%-mean(bl);        
                conderr(l,ori,sz) =std(frs(thisinds))./sqrt(length(thisinds));
                
                condz(l,ori,sz) = {(sc(thisinds)-mean(sc(thisinds)))/std(sc(thisinds))}; %ecker 2010
                condsc(l,ori,sz) = {sc(thisinds)};
                ff(l,ori,sz) = var(sc(thisinds))/mean(sc(thisinds));
                
                condlfpresp(l,ori,sz,:) = mean(lfpresp(thisinds,:),1);                
                
                condresperr(l,ori,sz,:) = nanstd(resp(thisinds,:),1,1)./sqrt(length(thisinds));
                if ~isnan(condresp(l,ori,sz,:))
                    [bincondresp(l,ori,sz,:),bta] = binit(condresp(l,ori,sz,:),binwidth);
                else
                    bincondresp(l,ori,sz,:) = binit(condresp(l,ori,sz,:),binwidth);
                end
                binconderr(l,ori,sz,:) = binit(condresperr(l,ori,sz,:),binwidth);
                
                mscc = []; bincc = [];
                for ii = 1:length(thisinds)-1
                    for jj = ii+1:length(thisinds)
                        help = corrcoef(resp(thisinds(ii),:),resp(thisinds(jj),:));
                        mscc = [mscc,help(1,2)];
                        help = corrcoef(binit(resp(thisinds(ii),:),binwidth),binit(resp(thisinds(jj),:),binwidth));
                        bincc = [bincc, help(1,2)];
                    end
                end
                msreliab(l,ori,sz) = nanmean(mscc);
                binreliab(l,ori,sz) = nanmean(bincc);
                eckerreliability(l,ori,sz) = var(frs(thisinds))/var(frs);

                thisruninds = intersect(thisinds,oktrials);
                if ~isempty(thisruninds)
                    runcondresp(l,ori,sz,:) = mean(resp(thisruninds,:),1);
                    runcondfr(l,ori,sz) = mean(frs(thisruninds));
                    runconderr(l,ori,sz) = std(frs(thisruninds))./sqrt(length(thisruninds));
                else
                    runcondresp(l,ori,sz,:) = nan(1,size(resp,2));
                    runcondfr(l,ori,sz) = NaN;
                    runconderr(l,ori,sz) = NaN;
                end  
                
                thisstillinds = intersect(thisinds,stilltrials);
                if ~isempty(thisstillinds)
                    stillcondresp(l,ori,sz,:) = mean(resp(thisstillinds,:),1);
                    stillcondfr(l,ori,sz) = mean(frs(thisstillinds));
                    stillconderr(l,ori,sz) = std(frs(thisstillinds))./sqrt(length(thisstillinds));
                else
                    stillcondresp(l,ori,sz,:) = nan(1,size(resp,2));
                    stillcondfr(l,ori,sz) = NaN;
                    stillconderr(l,ori,sz) = NaN;
                end  
                
            end
        end
    end
    
    bincondresp = bincondresp.*(1000/binwidth);
    bta = bta-prestim;
%     psthcondplot(bincondresp,binconderr,bta)
    cellbinresp(cell,:,:,:,:) = bincondresp;
    
    contindsnl = find(gratingInfo.size == 0 & light == 0);
    controlresp(1,:) = mean(resp(contindsnl,:),1);
    controlfr(cell,1) = mean(frs(contindsnl));
    controlerr(1) = std(frs(contindsnl))./sqrt(length(contindsnl));
    contindsl = find(gratingInfo.size == 0 & light == 1);
    controlresp(2,:) = mean(resp(contindsl,:),1);
    controlfr(cell,2) = mean(frs(contindsl));
    controlerr(2) = std(frs(contindsl))./sqrt(length(contindsl));
    
    [binnedctrnolight,bta] = binit(controlresp(1,:),binwidth);
    [binnedctrlight,bta] = binit(controlresp(2,:),binwidth);

    binnedcellrespl0(cell,:) = binnednolight;
    binnedcellrespl1(cell,:) = binnedlight;
    
    nolmaxfr(cell) = max(max(condfr(1,:,:)));
    lmaxfr(cell) = max(max(condfr(2,:,:)));
    
    cellz(cell,:,:,:) = condz;
    cellsc(cell,:,:,:) = condsc;
    cellff(cell,:,:,:) = ff;
    cellfr(cell,:,:,:) = condfr;
    celleckerrely(cell,:,:,:) = eckerreliability;
    cellmsrely(cell,:,:,:) = msreliab;
    cellbinrely(cell,:,:,:) = binreliab;
    celllfpspect(cell,:,:,:,:) = condlfpspect;
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
    
    [nloriprefratio(cell), nldirprefratio(cell), nlprefori, nlmeanori, nlosi(cell), nlmeandir, nldsi(cell)] = getOSI(squeeze(condfr(1,:,prefsize)),oris);
    [loriprefratio(cell), ldirprefratio(cell), lprefori, lmeanori, losi(cell), lmeandir, ldsi(cell)] = getOSI(squeeze(condfr(2,:,prefsize)),oris);
    
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
    
    subplot(2,2,2)
    errorbar(oris,squeeze(condfr(2,:,prefsize)),squeeze(conderr(2,:,prefsize)),'o-','color',lcol,'markersize',8,'linewidth',2)
    hold on
    errorbar(oris,squeeze(condfr(1,:,prefsize)),squeeze(conderr(1,:,prefsize)),'ko-','markersize',8,'linewidth',2)
    xlabel('shown orientation')
    ylabel('Firing rate [Hz]')
    set(gca,'xtick',oris)
    legend({'Light ON','Light OFF'})
    title([' OSI: ' num2str(nlosi(cell)) ' OSI Light: ' num2str(losi(cell))])
    
    oneorifr = mean(reshape(condfr(:,:,prefsize),2,4,2),3);
    prefori = find(oneorifr(1,:) == max(oneorifr(1,:)),1);
    ortho = mod(prefori+2,length(oris)/2); if ortho == 0, ortho = length(oris)/2; end
    
    cellcondresppreforil0(cell,:,:) = squeeze(condresp(1,prefori,:,:));
    cellcondresppreforil1(cell,:,:) = squeeze(condresp(2,prefori,:,:));
    
    preffr(cell,:) = oneorifr(:,prefori);
    
    subplot(2,2,3)
    errorbar(sizes,squeeze(nanmean(condfr(2,[prefori,prefori+(length(oris)/2)],:))),...
        squeeze(nanmean(conderr(2,[prefori,prefori+(length(oris)/2)],:))),'o-','color',lcol,'markersize',8,'linewidth',2);
    hold on
    errorbar(sizes,squeeze(nanmean(condfr(1,[prefori,prefori+(length(oris)/2)],:))),...
        squeeze(nanmean(conderr(1,[prefori,prefori+(length(oris)/2)],:))),'ko-','markersize',8,'linewidth',2);
    xlabel('shown patch size [vd]')
    ylabel('Firing rate [Hz]')
    legend({'Light ON','Light OFF'})    
    set(gca,'xtick',sizes)  
    title(['cell ' int2str(cell) ' preferred orientations' ' depth: ' int2str(result.depth) '  ' printname])
    
    sizetunelas = squeeze(nanmean(condfr(2,[prefori,prefori+(length(oris)/2)],:)));
    sizetunenolas = squeeze(nanmean(condfr(1,[prefori,prefori+(length(oris)/2)],:)));
    
    sil(cell) = (sizetunelas(find(sizetunelas == max(sizetunelas),1))-sizetunelas(end))/sizetunelas(find(sizetunelas == max(sizetunelas),1));
    sinl(cell) = (sizetunenolas(find(sizetunenolas == max(sizetunenolas),1))-sizetunenolas(end))/sizetunenolas(find(sizetunenolas == max(sizetunenolas),1));
    
    subplot(2,4,7)
    plot(spike)
    axis([0,40,-100,100])
    legend(['width: ' int2str(swidth(cell)) ' adiff: ' num2str(adiff(cell))])
    
    subplot(2,4,8)
    plot(ta,binnedctrnolight)
    hold on
    plot(ta,binnedctrlight,lcol)
    axis([0,2500,0,1])
    
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
    
    figure
    subplot(2,2,1)
    imagesc(speed);
    colorbar
    title(['oktrials: ' int2str(length(oktrials)) '/' int2str(size(speed,1))])
    xlabel('time [ms]')
    ylabel('trial number')
    
    subplot(2,2,2)
    errorbar(msta,mean(speed(find(~light),:)),std(speed(find(~light),:))./sqrt(length(find(~light))),'b')
    hold on
    errorbar(msta,mean(speed(find(light),:)),std(speed(find(light),:))./sqrt(length(find(light))),'r')
    xlabel('time [ms]')
    ylabel('average runspeed')
    legend({'light off' 'light on'})
    
    subplot(2,2,3)
    plot(mean(speed(:,respwin),2),frs,'.')
    hold on
    plot(mean(speed(find(light),respwin),2),frs(find(light)),'r.')
    xlabel('average runspeed of trial')
    ylabel('average firing rate of trial')
    
    subplot(2,2,4)
    barweb([nlfr(cell),lfr(cell);runnlfr,runlfr;norunnlfr,norunlfr],...
        [nlfrerr,lfrerr;runnlfrerr,runlfrerr;norunnlfrerr,norunlfrerr],...
        [],[{'all'};{'running only'};{'immobile only'}],['ANOVA factor running p: ' num2str(p(2))],...
        [],'firing rate [Hz]',[],[]);
    
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
    
    cell = cell + 1;
    disp('');
    
% end

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

% correlations
infrs = intersect(find(depth>500),prs);
irsresp = cellresp(infrs,:,:);
ncells = size(irsresp,1);

if length(infrs) ~= 0
    
    for i = 1:length(infrs)
        for j = 1:size(irsresp,2)
            binresp(i,j,:) = binit(squeeze(irsresp(i,j,:)),100);
        end
    end
    
    %method 1
    for i = 1:420
        cr(i,:,:) = corr(squeeze(irsresp(:,i,:))');
        bcr(i,:,:) = corr(squeeze(binresp(:,i,:))');
    end

    % method 2
    binsizes = [1,5,10,15,20,25,50,100];;
    for ccbin = 1:length(binsizes)
        pair = 1;
        for i = 1:ncells-1
            for j = i+1:ncells

                % equalize FR of pair
                resp1 = squeeze(irsresp(i,:,:));
                resp2 = squeeze(irsresp(j,:,:));
                s1 = find(resp1); s2 = find(resp2);
                if length(s1)>length(s2)
                    rp = randperm(length(s1)); % delete random spikes from resp1
                    resp1(s1(rp(1:(length(s1)-length(s2))))) = 0;
                else
                    rp = randperm(length(s2)); % delete random spikes from resp2
                    resp2(s2(rp(1:(length(s2)-length(s1))))) = 0;
                end
                % now resp1 and 2 should have equal no of spikes
                %bin
%                 ccbin = 50;
                for ii = 1:size(resp1,2)/binsizes(ccbin);
                    binresp1(:,ii) = sum(resp1(:,(ii-1)*binsizes(ccbin)+1:ii*binsizes(ccbin)),2);
                    binresp2(:,ii) = sum(resp2(:,(ii-1)*binsizes(ccbin)+1:ii*binsizes(ccbin)),2);
                end

                for st = 1:size(irsresp,2)
                    [r,p] = corrcoef(squeeze(binresp1(st,:)),squeeze(binresp2(st,:)));
                    rval(pair,st) = r(1,2);
                    pval(pair,st) = p(1,2);
                end
                pair = pair+1;
            end
        end

        rl0 = rval(:,find(~light));
        rl1 = rval(:,find(light));
        rsl0(ccbin,:) = nanmean(rl0,2);
        rsl1(ccbin,:) = nanmean(rl1,2);
    end

    [s,p] = ttest(nanmean(rl0,2),nanmean(rl1,2));
    figure
    plot(nanmean(rl0,2),nanmean(rl1,2),'k.')
    axis([-.2,1,-.2,1])
    line([-.2,1],[-.2,1],'color','k')
    axis([-.2,1,-.2,1])
    axis square
    line([-.2,1],[0,0],'color','k')
    line([0,0],[-.2,1],'color','k')
    xlabel('mean corrcoef over light OFF conditions')
    ylabel('mean corrcoef over light ON conditions')
    title(['all pairwise correlations between lower layer RS cells p: ' num2str(p)]);
%     if printyn
%         figSize = [30 21];
%         set(gcf,'paperunits','centimeters','papersize',figSize,'paperposition',[0 0 figSize])
%         print(gcf,[popprintpath ,  '18_LLRScorrelations.pdf'], '-dpdf' );
%     end
    
    rl0c0 = nanmean(rval(:,light == 0 & gratingInfo.Contrast == 0),2);
    rl0c1 = nanmean(rval(:,light == 0 & gratingInfo.Contrast == 1),2);
    rl1c0 = nanmean(rval(:,light == 1 & gratingInfo.Contrast == 0),2);
    rl1c1 = nanmean(rval(:,light == 1 & gratingInfo.Contrast == 1),2);
    
    figure
    [s,p] = ttest(rl0c0,rl0c1);
    plot(rl0c0,rl0c1,'.')
    line([-.1,.7],[-.1,.7],'color','k')
    axis([-.1,.7,-.1,.7])
    axis square
    xlabel('average correlation zero contrast (light OFF)')
    ylabel('average correlation full contrast (light OFF)')
    title(['correlation changes with contrast light OFF: p: ' num2str(p)])
    
    figure
    [s,p] = ttest(rl1c0,rl1c1);
    plot(rl1c0,rl1c1,'.')
    line([-.1,.7],[-.1,.7],'color','k')
    axis([-.1,.7,-.1,.7])
    axis square
    xlabel('average correlation zero contrast (light ON)')
    ylabel('average correlation full contrast (light ON)')
    title(['correlation changes with contrast light ON: p: ' num2str(p)])
    
    figure
    [s,p] = ttest(rl0c0,rl1c0);
    plot(rl0c0,rl1c0,'.')
    line([-.1,.7],[-.1,.7],'color','k')
    axis([-.1,.7,-.1,.7])
    axis square
    xlabel('average correlation zero contrast (light OFF)')
    ylabel('average correlation zero contrast (light ON)')
    title(['correlation changes with light, zero contrast: p: ' num2str(p)])
    
    figure
    [s,p] = ttest(rl0c1,rl1c1);
    plot(rl0c1,rl1c1,'.')
    line([-.1,.7],[-.1,.7],'color','k')
    axis([-.1,.7,-.1,.7])
    axis square
    xlabel('average correlation full contrast (light OFF)')
    ylabel('average correlation full contrast (light ON)')
    title(['correlation changes with light, full contrast: p: ' num2str(p)])
    
    [s,p] = ttest(rl0c1-rl0c0,rl1c1-rl1c0)
    
end



%lifetime sparseness
N = 34; %1000ms binned 28:61 on bta (500:1500)
for cl = 1:size(cellbinresp,1)
    for l = 1:2
        for o = 1:8
            for s = 1:5
                slt(cl,l,o,s) = (1 - ((1/N)* ((sum(cellbinresp(cl,l,o,s,28:61),5).^2)/(sum(cellbinresp(cl,l,o,s,28:61).^2,5))) ))...
                 / (1-(1/N));
            end
        end
    end
end

%population sparseness
N = size(cellbinresp,1);
for bn = 1:34
    for l = 1:2
        for o = 1:8
            for s = 1:5
                spop(bn,l,o,s) = (1 - ((1/N)* ((sum(cellbinresp(:,l,o,s,bn+27),1).^2)/(sum(cellbinresp(:,l,o,s,bn+27).^2,1))) ))...
                 / (1-(1/N));
            end
        end
    end
end

ii = 1;
for i = 1:size(cellz,1)-1
    for j = i+1:size(cellz,1)
        for l = 1:2
            for o = 1:8
                for s = 1:5
                    cv = cov(cellz{i,l,o,s},cellz{j,l,o,s});
                    cvs(ii,l,o,s) = cv(1,2);
                    cvcount = cov(cellsc{i,l,o,s},cellsc{j,l,o,s});
                    if size(cvcount) == 1
                        cvcounts(ii,l,o,s) = NaN;
                    else
                        cvcounts(ii,l,o,s) = cvcount(1,2);
                    end
                    cc = corrcoef(cellz{i,l,o,s},cellz{j,l,o,s});
                    if isnan(cc)
                        cc(ii,l,o,s) = NaN;
                    else
                        ccs(ii,l,o,s) = cc(1,2);
                    end
                    pairinds(ii,:) = [i,j];
                end
            end
        end
        ii = ii+1;
    end
end

mscrl0 = nanmean(nanmean(cvs(:,1,:,:),3),4);
mscrl1 = nanmean(nanmean(cvs(:,2,:,:),3),4);
mscountrl0 = nanmean(nanmean(cvcounts(:,1,:,:),3),4);
mscountrl1 = nanmean(nanmean(cvcounts(:,2,:,:),3),4);
mccl0 = nanmean(nanmean(ccs(:,1,:,:),3),4);
mccl1 = nanmean(nanmean(ccs(:,2,:,:),3),4);
mffl0 = nanmean(nanmean(cellff(:,1,:,:),3),4);
mffl1 = nanmean(nanmean(cellff(:,2,:,:),3),4);



figure
[s,p] = ttest(nlosi,losi);
plot(nlosi,losi,'k.')
line([0,1],[0,1],'color','k')
axis square
xlabel('OSI control');
ylabel('OSI light');
title([' p: ' num2str(p)])

figure
[s,p] = ttest(nldsi,ldsi);
plot(nldsi,ldsi,'k.')
line([0,1],[0,1],'color','k')
axis square
xlabel('DSI control');
ylabel('DSI light');
title([' p: ' num2str(p)])

figure
plot(nolmaxfr,lmaxfr,'ko')
% axis([-5,25,-5,25])
axis square
line([-5,25],[-5,25],'color','k')
xlabel('max firing rate no light [Hz]')
ylabel('max firing rate light [Hz]')

figure
plot(nlfr,lfr,'k.')
hold on
plot(nlfr(pfs),lfr(pfs),'ko')
axis square
line([0,30],[0,30],'color','k')
xlabel('mean firing rate no light [Hz]')
ylabel('mean firing rate light [Hz]')

figure
plot((lfr-nlfr)./(lfr+nlfr),depth,'k.')
hold on
plot((lfr(pfs)-nlfr(pfs))./(lfr(pfs)+nlfr(pfs)),depth(pfs),'ko')
line([0,0],[0,1000],'color','k')
axis ij
axis([-1.1,1.1,0,1000])
ylabel('depth[mum]')
xlabel('(light-nolight)/(light+nolight)')
title('mean firing rate changes by layer 4 suppression')

figure
plot((preffr(:,2)-preffr(:,1))./(preffr(:,2)+preffr(:,1)),depth,'k.')
hold on
plot((preffr(pfs,2)-preffr(pfs,1))./(preffr(pfs,2)+preffr(pfs,1)),depth(pfs),'ko')
line([0,0],[0,1000],'color','k')
axis ij
axis([-1.1,1.1,0,1000])
ylabel('depth[mum]')
xlabel('(light-nolight)/(light+nolight)')
title('preferred firing rate changes by layer 4 suppression')

figure
plot((controlfr(:,2)-controlfr(:,1))./(controlfr(:,2)+controlfr(:,1)),depth,'k.')
hold on
plot((controlfr(pfs,2)-controlfr(pfs,1))./(controlfr(pfs,2)+controlfr(pfs,1)),depth(pfs),'ko')
line([0,0],[0,1000],'color','k')
axis ij
axis([-1.1,1.1,0,1000])
ylabel('depth[mum]')
xlabel('(light-nolight)/(light+nolight)')
title('control firing rate changes by layer 4 suppression')

figure
plot(sinl,sil,'k.')
hold on
plot(sinl(pfs),sil(pfs),'ko')
axis([0,1,0,1])
axis square
line([0,1],[0,1],'color','k')
xlabel('SI ((max-largest)/max) no light')
ylabel('SI ((max-largest)/max) light')

figure
plot(sil-sinl,depth,'k.')
hold on
plot(sil(pfs)-sinl(pfs),depth(pfs),'ko')
axis ij
line([0,0],[250,750],'color','k')
ylabel('depth[mum]')
xlabel('SI light - SI no light')

disp('');