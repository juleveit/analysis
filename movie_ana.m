function movie_ana

animalid = '180330';
block = 5;
rfblock = 3;
electrodes = [16];
penangle = [25];
lcol = 'r';

onlymod = 0;
sfc = 0;

basepath = ['C:\Users\Julia\work\data\' animalid '\'];
supath = [basepath 'singleunits\'];
basename = [animalid '_block' int2str(block) '_tet'];
snbasename = [animalid '_block' int2str(rfblock) '_tet'];
load('C:\Users\Julia\work\Matlab\movies\fsmovie.mat');

files = dir([supath, basename, '*.mat']);
rffiles = dir([supath, snbasename, '*.mat']);

prestim = 300;
poststim = 300;
respwin = 1001:2000; % after stimulus onset
respwin = respwin+prestim;
freqbinwidth = 5;

cell = 1;
for fi = 1:length(files)
    
    if strfind(files(fi).name, 'MU')
        continue;
    end
    
    load([supath, files(fi).name]);
    
    % calc spiking to see if includable
    msStimes = round(result.spikes);
    if ~isempty(msStimes) & msStimes(1) == 0, msStimes(1) = 1; end
    
    chan = zeros(1,length(result.lfp));
    chan(msStimes) = 1;
    
    wvchan = find(var(result.waveforms) == max(var(result.waveforms)));
    
    
    % LFP GAMMA
%     lfp = lfp-mean(lfp(1,1:beg),2);
%     stdbl = std(lfp(:,1:beg)');
%     lfp = lfp./stdbl;
    
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
        

    trialdur = result.movieLengthSecs*1000;
    msstamps = result.msstamps;
    
    if length(msstamps)~=length(result.light)
%         msstamps(83) = []; % 151109 - 12
        msstamps(44) = []; % 1160923 - 2
        result.msstamps = msstamps;
        save([supath, files(fi).name],'result');
%         pause;
    end
    if ~isfield(result,'position')
        result.position = [-18,20];
        save([supath, files(fi).name],'result');
    end
    
%     trialnfft = 2^nextpow2(800);
%     trialfax = sr/2*linspace(0,1,trialnfft/2+1);
    for i = 1:length(msstamps)
        resp(i,:) = chan(msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim); 
        lfpresp(i,:) = lfp(msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim);
%         y = fft(lfpresp(i,1001:1800),trialnfft);
%         lfpspect(i,:) = abs(y(1:trialnfft/2+1));
        [lfpspect(i,:),trialfax] = pmtm(lfpresp(i,1501:2300),3,[],sr);
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
    
    
    % figure out sufficiently high and nonvariable runspeed trials
    meanspeed = mean(speed(:,respwin),2);
    stdspeed = std(speed(:,respwin),1,2);
    notstill = find(meanspeed>1);
    okspeed = find(meanspeed>( mean(meanspeed(notstill))-(1.5*std(meanspeed(notstill))) ) );
    okvar = find(stdspeed<( mean(stdspeed(notstill))+(1.5*std(stdspeed(notstill)))) & stdspeed>.5);
    oktrials = intersect(okspeed,okvar);
    nonoktrials = 1:size(resp,1); nonoktrials(oktrials) = [];
    stilltrials = 1:size(resp,1); stilltrials(notstill) = [];

    
    frs = sum(resp(:,respwin),2)./(length(respwin)/1000);
    bl = sum(resp(:,1:prestim),2)./(prestim/1000);
    
    %determine if cell is visually modulated
    blfr = sum(resp(:,1:prestim),2);
    vrfr = sum(resp(:,prestim+40:2*prestim+40-1),2);
    vismod(fi) = ttest2(blfr,vrfr);
    visdriven(fi) = mean(vrfr)>=mean(blfr)+2; % average firing rate is increase at least 2Hz above baseline
    
    if onlymod & ~vismod(fi) %~(lightmod(fi) & vismod(fi))
        continue;
    end
    
    cellname{cell} = files(fi).name;
    
    i = strfind(files(cell).name, 'tet');
    if strcmp(files(cell).name(i+4),'_')
        tetno = strread(files(cell).name(i+3)); % single character number
    else        
        tetno = strread(files(cell).name(i+3:i+4)); % number >10
    end
    tetnos(cell) = tetno;
    if tetno>electrodes(1)/4
        elec(cell) = 2; cellstr = 'elec2';
    else
        elec(cell) = 1; cellstr = 'elec1';
    end
    
    spike = result.waveforms(:,wvchan); 
    interpspike = spline(1:32,spike,1:.1:32);
    [adiff(cell),swidth(cell)] = spikequant(interpspike);
    
    
    l0 = result.light == 0; l1 = result.light == 1;
    
    % phases
    tmp = zeros(size(gphaseresp));
    tmp(find(resp)) = gphaseresp(find(resp));
    l0tmp = tmp(l0,:); l1tmp = tmp(l1,:);
    phasesl0{cell} = l0tmp(find(l0tmp));
    phasesl1{cell} = l1tmp(find(l1tmp));
      
    phaserl0(cell) = circ_r(phasesl0{cell});
    phaserl1(cell) = circ_r(phasesl1{cell});
    cmeanl0(cell) = circ_mean(phasesl0{cell});
    cmeanl1(cell) = circ_mean(phasesl1{cell});
    
    if sfc
        tmpi = zeros(size(allphaseresp));
        for i = 1:size(allphaseresp,1)
            tmp = zeros(size(gphaseresp));
            tmp(find(resp)) = allphaseresp(i,find(resp));
            tmpi(i,:,:) = tmp;
        end
        allphasemat = tmpi;
        for i = 1:size(allphaseresp,1)
            allphases{i} = allphasemat(i,find(squeeze(allphasemat(i,:,:)))); 
            allphaser(cell,i) = circ_r(allphases{i}');
            allcmean(cell,i) = circ_mean(allphases{i}');
        end
    end
    
    depth(cell) = result.depth;
    cellresp(cell,:,:) = resp;
    celllfpresp(cell,:,:) = lfpresp;
    cellchan(cell,:) = chan;
    
    %get RF of cell, too
    [on,off,gaussfiton,gaussfitoff,rson(cell),rsoff(cell),xax,yax] = get_rf([supath, rffiles(fi).name]);
    
    msta = linspace(-prestim,trialdur+poststim,size(resp,2));
    
    baseline = mean(bl);
    baselineerr = std(bl)./(sqrt(size(bl,1)));
    
    binwidth = 33.333;
    clevels = unique(result.contrast); 
    widths = unique(result.aperture); 
    movies = unique(result.movieno); 
        
    for l = 1:2
        for mv = 1:length(movies)
            for co = 1:length(clevels)
                for sz = 1:length(widths)
                    thisinds = find(result.movieno == movies(mv) & ...
                        result.contrast == clevels(co) &...
                        result.aperture == widths(sz) &...
                        result.light == l-1);
                    condn(l,mv,co,sz) = length(thisinds);
                    condresp(l,mv,co,sz,:) = nanmean(resp(thisinds,:),1);
                    condresperr(l,mv,co,sz,:) = nanstd(resp(thisinds,:),1,1)./sqrt(length(thisinds));
                    condlfpresp(l,mv,co,sz,:) = nanmean(lfpresp(thisinds,:),1);
                    if ~isnan(condresp(l,mv,co,sz,:))
                        [bincondresp(l,mv,co,sz,:),bta] = binit(condresp(l,mv,co,sz,:),binwidth);
                    else
                        bincondresp(l,mv,co,sz,:) = binit(condresp(l,mv,co,sz,:),binwidth);
                    end
                    binconderr(l,mv,co,sz,:) = binit(condresperr(l,mv,co,sz,:),binwidth);
                    condfr(l,mv,co,sz) = nanmean(frs(thisinds));
                    trialfr(l,mv,co,sz,:) = frs(thisinds);
                    conderr(l,mv,co,sz) = nanstd(frs(thisinds))./sqrt(length(thisinds));
                    condff(l,mv,co,sz) = var(frs(thisinds))/mean(frs(thisinds));
                    
                    condsparseness(l,mv,co,sz) = sparseness(squeeze(bincondresp(l,mv,co,sz,:)));
                    condsparsenesswin(l,mv,co,sz) = sparseness(squeeze(bincondresp(l,mv,co,sz,find(bta>respwin(1),1):find(bta>respwin(end),1)-1)));
                    
                    %get trial to trial reliability
                    k = 1; clear tc;
                    binrespwin = find(bta>respwin(1),1):find(bta>respwin(end),1)-1;
                    for i = 1:length(thisinds)-1
                        for j = i+1:length(thisinds)
                            a = binit(resp(thisinds(i),:),binwidth);
                            b = binit(resp(thisinds(j),:),binwidth);
                            tc(k) = nancorr(a(binrespwin),b(binrespwin));
                            
                            a = binit(lfpresp(thisinds(i),:),binwidth);
                            b = binit(lfpresp(thisinds(j),:),binwidth);
                            lfptc(k) = nancorr(a(binrespwin),b(binrespwin));
                            
                            k = k+1;
                        end
                    end
                    condrely(l,mv,co,sz) = nanmean(tc);
                    condlfprely(l,mv,co,sz) = nanmean(lfptc);
                            
                end
            end
        end
    end

    bincondresp = bincondresp.*(1000/binwidth);
    bta = bta-prestim;
    
    
    lightmod(cell) = ttest2(frs(find(result.light)),frs(find(result.light == 0)));
    aperturevismod(cell) = ttest(frs(find(result.light == 0 & result.aperture ~= 0 )),bl(find(result.light == 0 & result.aperture ~= 0)));
    contextvismod(cell) = ttest(frs(find(result.light == 0 & result.aperture == 0 )),bl(find(result.light == 0 & result.aperture == 0)));
    
%     figure 
%     rasterplot2(resp,result.aperture == widths(2),result.aperture == widths(1),msta);
%     title(cellname(cell))
 
    figure 
    rasterplot4(resp,result.aperture == widths(2) & result.light == 0 & result.contrast == 1,...
        result.aperture == widths(2) & result.light == 1 & result.contrast == 1,...
        result.aperture == widths(1) & result.light == 0 & result.contrast == 1,...
        result.aperture == widths(1) & result.light == 1 & result.contrast == 1, msta);
    title(cellname(cell))
    
    clts(cell,:,:,:,:) = condsparseness;
    cltswin(cell,:,:,:,:) = condsparsenesswin;
    cfr(cell,:,:,:,:) = condfr;
    cff(cell,:,:,:,:) = condff;
    cellcondrely(cell,:,:,:,:) = condrely;
    celllfprely(cell,:,:,:,:) = condlfprely;
    cellcondresp(cell,:,:,:,:,:) = condresp;
    celltrialfr(cell,:,:,:,:,:) = trialfr;
    bincondcellresp(cell,:,:,:,:,:) = bincondresp;
    
    [mc,mci] = max(clevels);    
    
    figure
    subplot(2,2,1)
    semilogy(trialfax,mean(lfpspect(l0&result.aperture == widths(2),:)),'linewidth',2);
    hold on
    semilogy(trialfax,mean(lfpspect(l1&result.aperture == widths(2),:)),'r','linewidth',2);
    semilogy(trialfax,mean(lfpspect(l0&result.aperture == widths(1),:)),'c','linewidth',2);
    semilogy(trialfax,mean(lfpspect(l1&result.aperture == widths(1),:)),'m','linewidth',2);
%     semilogy(trialfax,mean(lfpspect(l0&result.aperture == widths(2),:))-(std(lfpspect(l0&result.aperture == widths(2),:))./sqrt(size(l0&result.aperture == widths(2),2))));
%     semilogy(trialfax,mean(lfpspect(l0&result.aperture == widths(2),:))+(std(lfpspect(l0&result.aperture == widths(2),:))./sqrt(size(l0&result.aperture == widths(2),2))));
%     semilogy(trialfax,mean(lfpspect(l1&result.aperture == widths(2),:))-(std(lfpspect(l1&result.aperture == widths(2),:))./sqrt(size(l1&result.aperture == widths(2),2))),'r');
%     semilogy(trialfax,mean(lfpspect(l1&result.aperture == widths(2),:))+(std(lfpspect(l1&result.aperture == widths(2),:))./sqrt(size(l1&result.aperture == widths(2),2))),'r');
%     semilogy(trialfax,mean(lfpspect(l0&result.aperture == widths(1),:))-(std(lfpspect(l0&result.aperture == widths(1),:))./sqrt(size(l0&result.aperture == widths(1),2))),'c');
%     semilogy(trialfax,mean(lfpspect(l0&result.aperture == widths(1),:))+(std(lfpspect(l0&result.aperture == widths(1),:))./sqrt(size(l0&result.aperture == widths(1),2))),'c');
%     semilogy(trialfax,mean(lfpspect(l1&result.aperture == widths(1),:))-(std(lfpspect(l1&result.aperture == widths(1),:))./sqrt(size(l1&result.aperture == widths(1),2))),'m');
%     semilogy(trialfax,mean(lfpspect(l1&result.aperture == widths(1),:))+(std(lfpspect(l1&result.aperture == widths(1),:))./sqrt(size(l1&result.aperture == widths(1),2))),'m');
    axis([0,120,...
        min([min(squeeze(mean(lfpspect(:,1:125)))),min(squeeze(mean(lfpspect(:,1:125))))]),...
        max([max(squeeze(mean(lfpspect(:,1:125)))),max(squeeze(mean(lfpspect(:,1:125))))])])
    xlabel('frequency [Hz]')
    ylabel('spectral power')
    legend([{['CRF L0']},{['CRF L1']},{['FS L0']},{['FS L1']}],'location','ne')
    title('LFP spectrum during light on vs off')
    
    subplot(2,2,2)
    plot(mean(gammaresp(l0,:)));
    hold on
    plot(mean(gammaresp(l1,:)),'r')
    title(['gamma power in time depth: ' int2str(depth(cell))])
    
    subplot(2,2,3)
    [to0,ro0] = rose(phasesl0{cell});
    [to1,ro1] = rose(phasesl1{cell});
    polar(to0,ro0,'b')
    hold on    
    polar(to1,ro1,'r')
    title(['gamma phase locking of unit ' int2str(cell) ' spikewidth: ' int2str(swidth(cell))])
    
    subplot(2,2,4)
    plot(mean(lfpresp(l0,:)));
    hold on
    plot(mean(lfpresp(l1,:)),'r')
    
    
    figure
    subplot(2,1,1)
    plot(bta, squeeze(bincondresp(1,1,mci,2,:)),'b') % movie 1 fc aperture ligth off
    hold on
    plot(bta, squeeze(bincondresp(2,1,mci,2,:)),'r') % movie 1 fc aperture ligth on
    plot(bta, squeeze(bincondresp(1,1,mci,1,:)),'c') % movie 1 fc full screen ligth off
    plot(bta, squeeze(bincondresp(2,1,mci,1,:)),'m') % movie 1 fc full screen ligth off
    title(cellname{cell})
    legend({'CRF light OFF','CRF light ON','contextual light OFF','contextual light ON'})
    ax = axis;
    axis([-300,3300,ax(3),ax(4)])
    line([0,0],[ax(3),ax(4)],'color','k')
    line([3000,3000],[ax(3),ax(4)],'color','k')
    line([1000,1000],[ax(3),ax(4)],'color','r')
    line([2000,2000],[ax(3),ax(4)],'color','r')
    
    subplot(2,2,3)
    imagesc(xax,yax,on)
    hold on
    plot_orrf_absdeg(gaussfiton,1,'w',2)
    axis square
    plot_circle(result.position(1),result.position(2),max(unique(result.aperture))/2,'k',2)
    
    subplot(2,2,4)
    imagesc(xax,yax,off)
    hold on
    plot_orrf_absdeg(gaussfitoff,1,'w',2)
    axis square
    plot_circle(result.position(1),result.position(2),max(unique(result.aperture))/2,'k',2)
    

    
%     subplot(2,1,2)
%     plot(bta, squeeze(bincondresp(1,2,mci,2,:)),'b') % movie 2 fc aperture ligth off
%     hold on
%     plot(bta, squeeze(bincondresp(2,2,mci,2,:)),'r') % movie 2 fc aperture ligth on
%     plot(bta, squeeze(bincondresp(1,2,mci,1,:)),'c') % movie 2 fc full screen ligth off
%     plot(bta, squeeze(bincondresp(2,2,mci,1,:)),'m') % movie 2 fc full screen ligth off
%     title(cellname{cell})
%     ax = axis;
%     axis([-300,3300,ax(3),ax(4)])
%     line([0,0],[ax(3),ax(4)],'color','k')
%     line([3000,3000],[ax(3),ax(4)],'color','k')
%     line([1000,1000],[ax(3),ax(4)],'color','r')
%     line([2000,2000],[ax(3),ax(4)],'color','r')
    
    cell = cell+1;
    disp('');
    
end

figure
plot(swidth,adiff,'k.')
xlabel('spike width')
ylabel('amplitude diff')

% adjust depth according to penetration angle
for cl = 1:length(depth)
    depth(cl) = depth(cl).*cosd(penangle(elec(cl)));
end
% l,mv,co,sz

ana_elec = 1;

pfs = find(swidth<125&elec==ana_elec);
prs = find(swidth>=125&elec==ana_elec);
pfsv = swidth<125&elec==ana_elec;
prsv = swidth>=125&elec==ana_elec;


mv = 1;

%firing rates for context
figure
[sr,pr] = ttest(cfr(prs,1,mv,mci,2),cfr(prs,1,mv,mci,1));
[sf,pf] = ttest(cfr(pfs,1,mv,mci,2),cfr(pfs,1,mv,mci,1));
plot(cfr(prs,1,mv,mci,2),cfr(prs,1,mv,mci,1),'o','markerfacecolor','b','markersize',5)
hold on
plot(cfr(pfs,1,mv,mci,2),cfr(pfs,1,mv,mci,1),'ro','markerfacecolor','r','markersize',5)
refline(1,0)
axis square
title(['p RS: ' num2str(pr) '   p FS: ' num2str(pf)]);
xlabel('mean firing rate light OFF CRF')
ylabel('mean firing rate light OFF context')

figure
[sr,pr] = ttest(cfr(prs,2,mv,mci,2),cfr(prs,2,mv,mci,1));
[sf,pf] = ttest(cfr(pfs,2,mv,mci,2),cfr(pfs,2,mv,mci,1));
plot(cfr(prs,2,mv,mci,2),cfr(prs,2,mv,mci,1),'o','markerfacecolor','b','markersize',5)
hold on
plot(cfr(pfs,2,mv,mci,2),cfr(pfs,2,mv,mci,1),'ro','markerfacecolor','r','markersize',5)
refline(1,0)
axis square
title(['p RS: ' num2str(pr) '   p FS: ' num2str(pf)]);
xlabel('mean firing rate light ON CRF')
ylabel('mean firing rate light ON context')

%firing rates for light
figure
[sr,pr] = ttest(cfr(prs,1,mv,mci,2),cfr(prs,2,mv,mci,2));
[sf,pf] = ttest(cfr(pfs,1,mv,mci,2),cfr(pfs,2,mv,mci,2));
plot(cfr(prs,1,mv,mci,2),cfr(prs,2,mv,mci,2),'o','markerfacecolor','b','markersize',5)
hold on
plot(cfr(pfs,1,mv,mci,2),cfr(pfs,2,mv,mci,2),'ro','markerfacecolor','r','markersize',5)
refline(1,0)
axis square
title(['p RS: ' num2str(pr) '   p FS: ' num2str(pf)]);
xlabel('mean firing rate light OFF CRF')
ylabel('mean firing rate light ON CRF')

figure
[sr,pr] = ttest(cfr(prs,1,mv,mci,1),cfr(prs,2,mv,mci,1));
[sf,pf] = ttest(cfr(pfs,1,mv,mci,1),cfr(pfs,2,mv,mci,1));
plot(cfr(prs,1,mv,mci,1),cfr(prs,2,mv,mci,1),'o','markerfacecolor','b','markersize',5)
hold on
plot(cfr(pfs,1,mv,mci,1),cfr(pfs,2,mv,mci,1),'ro','markerfacecolor','r','markersize',5)
refline(1,0)
axis square
title(['p RS: ' num2str(pr) '   p FS: ' num2str(pf)]);
xlabel('mean firing rate light OFF context')
ylabel('mean firing rate light ON context')

%sparseness for context
figure
[sr,pr] = ttest(clts(prs,1,mv,mci,2),clts(prs,1,mv,mci,1));
[sf,pf] = ttest(clts(pfs,1,mv,mci,2),clts(pfs,1,mv,mci,1));
plot(clts(prs,1,mv,mci,2),clts(prs,1,mv,mci,1),'o','markerfacecolor','b','markersize',5)
hold on
plot(clts(pfs,1,mv,mci,2),clts(pfs,1,mv,mci,1),'ro','markerfacecolor','r','markersize',5)
refline(1,0)
axis square
title(['p RS: ' num2str(pr) '   p FS: ' num2str(pf)]);
xlabel('sparseness light OFF CRF')
ylabel('sparseness light OFF context')

figure
[sr,pr] = ttest(clts(prs,2,mv,mci,2),clts(prs,2,mv,mci,1));
[sf,pf] = ttest(clts(pfs,2,mv,mci,2),clts(pfs,2,mv,mci,1));
plot(clts(prs,2,mv,mci,2),clts(prs,2,mv,mci,1),'o','markerfacecolor','b','markersize',5)
hold on
plot(clts(pfs,2,mv,mci,2),clts(pfs,2,mv,mci,1),'ro','markerfacecolor','r','markersize',5)
refline(1,0)
axis square
title(['p RS: ' num2str(pr) '   p FS: ' num2str(pf)]);
xlabel('sparseness light ON CRF')
ylabel('sparseness light ON context')

%sparseness for light
figure
[sr,pr] = ttest(clts(prs,1,mv,mci,2),clts(prs,2,mv,mci,2));
[sf,pf] = ttest(clts(pfs,1,mv,mci,2),clts(pfs,2,mv,mci,2));
plot(clts(prs,1,mv,mci,2),clts(prs,2,mv,mci,2),'o','markerfacecolor','b','markersize',5)
hold on
plot(clts(pfs,1,mv,mci,2),clts(pfs,2,mv,mci,2),'ro','markerfacecolor','r','markersize',5)
refline(1,0)
axis square
title(['p RS: ' num2str(pr) '   p FS: ' num2str(pf)]);
xlabel('sparseness light OFF CRF')
ylabel('sparseness light ON CRF')

figure
[sr,pr] = ttest(clts(prs,1,mv,mci,1),clts(prs,2,mv,mci,1));
[sf,pf] = ttest(clts(pfs,1,mv,mci,1),clts(pfs,2,mv,mci,1));
plot(clts(prs,1,mv,mci,1),clts(prs,2,mv,mci,1),'o','markerfacecolor','b','markersize',5)
hold on
plot(clts(pfs,1,mv,mci,1),clts(pfs,2,mv,mci,1),'ro','markerfacecolor','r','markersize',5)
refline(1,0)
axis square
title(['p RS: ' num2str(pr) '   p FS: ' num2str(pf)]);
xlabel('sparseness light OFF context')
ylabel('sparseness light ON context')

disp('');

end

