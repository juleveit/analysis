function sizeanalysis_broadband

animalid = '150902';
block = 3;
channel = 3;
lcol = 'r'; %lasercolor

onlymod = 0;
printyn = 1;
sfc = 1;

supath = ['C:\Users\Julia\work\data\' animalid '\broadband\'];
% supath = ['C:\Users\Julia\work\data\' animalid '\multiunits\'];
basename = [animalid '_block' int2str(block) '_channel' int2str(channel)];

files = dir([supath, basename, '*.mat']);

prestim = 300;
poststim = 700;
respwin = 501:1500; % after stimulus onset
respwin = respwin+prestim;
freqbinwidth = 5;


load([supath, files.name]);

sr = 1000;
lfp = result.lfp';
nfft = 2^nextpow2(length(lfp));
fax = sr/2*linspace(0,1,nfft/2+1);
y = fft(lfp,nfft);
lfpspectrum = abs(y(1:nfft/2+1));

%     plot(fax,lfpspectrum)
gamma = eegfilt(lfp,sr,25,35);
h = hilbert(gamma); gpow = abs(h); gphas = angle(h);
% %
if sfc
    for i = 1:120/freqbinwidth
        filtmat(i,:) = eegfilt(lfp,sr,(i-1)*freqbinwidth+1,i*freqbinwidth);
        h = hilbert(filtmat(i,:));
        powmat(i,:) = abs(h); phasmat(i,:) = angle(h);
    end
end


trialdur = result.stimduration*1000;
msstamps = result.msstamps;

if length(msstamps)~=length(result.light)
    %          disp('');
    %         msstamps([62,108,147]) = []; % for 140703 block 8
    %         msstamps([161]) = []; % for 141204 block 3
    %         msstamps([303]) = []; % for 150407 block 5
    %         msstamps([169,336]) = []; % for 150523 block 11
    %         msstamps([24]) = []; % for 150730 block 11
    %         msstamps([207]) = []; % for 150730 block 11
    %         result.msstamps = msstamps;
    %         save([supath, files(fi).name],'result');
    pause;
end


for i = 1:length(msstamps)
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

% figure out sufficiently high and nonvariable runspeed trials
meanspeed = mean(speed(:,respwin),2);
stdspeed = std(speed(:,respwin),1,2);
notstill = find(meanspeed>1);
okspeed = find(meanspeed>( mean(meanspeed(notstill))-(1.5*std(meanspeed(notstill))) ) );
okvar = find(stdspeed<( mean(stdspeed(notstill))+(1.5*std(stdspeed(notstill)))) & stdspeed>.5);
oktrials = intersect(okspeed,okvar);
nonoktrials = 1:size(lfpresp,1); nonoktrials(oktrials) = [];
stilltrials = 1:size(lfpresp,1); stilltrials(notstill) = [];


sizes = unique(result.gratingInfo.size);  sizes(find(sizes == 0)) = []; %delete control condition
oris = unique(result.gratingInfo.Orientation); oris(find(oris == -1)) = [];
for l = 1:2
    for ori = 1:length(oris)
        for sz = 1:length(sizes)
            thisinds = find(result.gratingInfo.Orientation == oris(ori) &...0
                result.gratingInfo.size == sizes(sz) & ...
                result.light == l-1);
            condlfpspect(l,ori,sz,:) = nanmean(lfpspect(thisinds,:));
            condlfpresp(l,ori,sz,:) = mean(lfpresp(thisinds,:),1);
            alllfpresp{l,ori,sz} = lfpresp(thisinds,:);
          
            thisruninds = intersect(thisinds,oktrials);
            condruntrials{l,ori,sz} = thisruninds;
            if ~isempty(thisruninds)
                runcondlfpresp(l,ori,sz,:) = mean(lfpresp(thisruninds,:),1);
                runcondlfpspect(l,ori,sz,:) = nanmean(lfpspect(thisruninds,:));
            else
                runcondlfpresp(l,ori,sz,:) = nan(1,size(lfpresp,2));
                runcondlfpspect(l,ori,sz,:) = nan(1,size(lfpspect,2));
            end
            
            thisstillinds = intersect(thisinds,stilltrials);
            condstilltrials{l,ori,sz} = thisstillinds;
            if ~isempty(thisstillinds)
                stillcondlfpresp(l,ori,sz,:) = mean(lfpresp(thisstillinds,:),1);
                stillcondlfpspect(l,ori,sz,:) = nanmean(lfpspect(thisstillinds,:));
            else
                stillcondlfpresp(l,ori,sz,:) = nan(1,size(lfpresp,2));
                stillcondlfpspect(l,ori,sz,:) = nan(1,size(lfpspect,2));
            end
            
        end
    end
end

bincondresp = bincondresp.*(1000/binwidth);
bta = bta-prestim;
%     psthcondplot(bincondresp,binconderr,bta)
cellbinresp(cell,:,:,:,:) = bincondresp;

contindsnl = find(result.gratingInfo.size == 0 & result.light == 0);
controlresp(1,:) = mean(resp(contindsnl,:),1);
controlfr(cell,1) = mean(frs(contindsnl));
controlerr(1) = std(frs(contindsnl))./sqrt(length(contindsnl));
contindsl = find(result.gratingInfo.size == 0 & result.light == 1);
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
%     [loriprefratio(cell), ldirprefratio(cell), lprefori, lmeanori, losi(cell), lmeandir, ldsi(cell)] = getOSI(squeeze(condfr(2,:,prefsize)),oris);

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
%     title([' OSI: ' num2str(nlosi(cell)) ' OSI Light: ' num2str(losi(cell))])

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

%     sil(cell) = (sizetunelas(find(sizetunelas == max(sizetunelas),1))-sizetunelas(end))/sizetunelas(find(sizetunelas == max(sizetunelas),1));
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
%     condl0 = (result.gratingInfo.Orientation == oris(o) | result.gratingInfo.Orientation == oris(o+4)) & (result.gratingInfo.size == sizes(s) | result.gratingInfo.size == sizes(s+1)) & result.light == 0;
%     condl1 = (result.gratingInfo.Orientation == oris(o) | result.gratingInfo.Orientation == oris(o+4)) & (result.gratingInfo.size == sizes(s) | result.gratingInfo.size == sizes(s+1)) & result.light == 1;
%     condl0 = result.gratingInfo.Orientation == oris(o) & result.gratingInfo.size == sizes(s) & result.light == 0;
%     condl1 = result.gratingInfo.Orientation == oris(o) & result.gratingInfo.size == sizes(s) & result.light == 1;
condl0 = result.light == 0; condl1 = result.light == 1;
%     figure, rasterplot(resp,condl0,condl1,msta);
%     line([0,2000],[37,37],'color','k','linewidth',2)
%     line([500,1500],[33,33],'color','r','linewidth',2)

% running figure
runlfr = mean(frs(intersect(find(result.light),oktrials)));
runnlfr = mean(frs(intersect(find(~result.light),oktrials)));
norunlfr = mean(frs(intersect(find(result.light),stilltrials)));
norunnlfr = mean(frs(intersect(find(~result.light),stilltrials)));
lfrerr = std(frs(find(result.light)))./sqrt(length(find(result.light)));
nlfrerr = std(frs(find(~result.light)))./sqrt(length(find(~result.light)));
runlfrerr = std(frs(intersect(find(result.light),oktrials)))./sqrt(length(intersect(find(result.light),oktrials)));
runnlfrerr = std(frs(intersect(find(~result.light),oktrials)))./sqrt(length(intersect(find(~result.light),oktrials)));
norunlfrerr = std(frs(intersect(find(result.light),stilltrials)))./sqrt(length(intersect(find(result.light),stilltrials)));
norunnlfrerr = std(frs(intersect(find(~result.light),stilltrials)))./sqrt(length(intersect(find(~result.light),stilltrials)));

l1r1 = frs(intersect(find(result.light),oktrials));
l0r1 = frs(intersect(find(~result.light),oktrials));
l1r0 = frs(intersect(find(result.light),stilltrials));
l0r0 = frs(intersect(find(~result.light),stilltrials));
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
errorbar(msta,mean(speed(find(~result.light),:)),std(speed(find(~result.light),:))./sqrt(length(find(~result.light))),'b')
hold on
errorbar(msta,mean(speed(find(result.light),:)),std(speed(find(result.light),:))./sqrt(length(find(result.light))),'r')
xlabel('time [ms]')
ylabel('average runspeed')
legend({'light off' 'light on'})

subplot(2,2,3)
plot(mean(speed(:,respwin),2),frs,'.')
hold on
plot(mean(speed(find(result.light),respwin),2),frs(find(result.light)),'r.')
xlabel('average runspeed of trial')
ylabel('average firing rate of trial')

subplot(2,2,4)
barweb([nlfr(cell),lfr(cell);runnlfr,runlfr;norunnlfr,norunlfr],...
    [nlfrerr,lfrerr;runnlfrerr,runlfrerr;norunnlfrerr,norunlfrerr],...
    [],[{'all'};{'running only'};{'immobile only'}],['ANOVA factor running p: ' num2str(p(2))],...
    [],'firing rate [Hz]',[],[]);

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  

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
cell1 = 9; cell2 = 26;
for i = 1:size(celllfpresp,2)
    [coh(i,:),cfx] = mscohere(squeeze(celllfpresp(cell1,i,respwin)),squeeze(celllfpresp(cell2,i,respwin)),[],[],512,1000);
end
for l = 1:2
    for ori = 1:length(oris)
        for sz = 1:length(sizes)
            thisinds = find(result.gratingInfo.Orientation == oris(ori) &...0
                result.gratingInfo.size == sizes(sz) & ...
                result.light == l-1);
            condcoher(l,ori,sz,:) = nanmean(coh(thisinds,:),1);
        end
    end
end
contcoher(1,:) = nanmean(coh(contindsnl,:),1);
contcoher(2,:) = nanmean(coh(contindsl,:),1);
            
