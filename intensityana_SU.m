function stimfreqana_SU

animalid = '160907';
block = 9;

basepath = ['C:\Users\Julia\work\data\' animalid '\'];
% supath = [basepath 'multiunits\'];
supath = [basepath 'singleunits\'];
basename = [animalid '_block' int2str(block) '_tet'];

files = dir([supath, basename, '*.mat']);

prestim = 300;
% poststim = 300;
poststim = 300;
respwin = 501:1500; % after stimulus onset
% respwin = 501:4500; % after stimulus onset
respwin = respwin+prestim;

params.tapers = [2,3]; params.Fs = 1000; params.err = [2, 0.05]; params.trialave = 1;

cell = 1;
for fi = 1:length(files)
    
%     if strfind(files(fi).name, 'MU')
%         continue;
%     end
    
    load([supath, files(fi).name]);
    
    % calc spiking to see if includable
    msStimes = round(result.spikes);
    if ~isempty(msStimes) & msStimes(1) == 0, msStimes(1) = 1; end
    
    chan = zeros(1,length(result.lfp));
    chan(msStimes) = 1;
    
    wvchan = find(var(result.waveforms) == max(var(result.waveforms)));
    
    trialdur = result.stimduration*1000;
    msstamps = result.msstamps;
    if length(msstamps)~=length(result.light)
%         msstamps(385) = []; % for 140703 block 5
%         result.msstamps = msstamps;
%         save([supath, files(fi).name],'result');
        pause;
    end
    
    for i = 1:length(msstamps)
        resp(i,:) = chan(msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim); 
        lfpresp(i,:) = result.lfp(msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim);
        speed(i,:) = result.runspeed(msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim); 
        [lfpspect(i,:),fax] = pmtm(lfpresp(i,respwin),3,[],1000);
        [S(i,:,:),t,f]=mtspecgramc(lfpresp(i,:),[.2,.1],params);
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
    
    cellname{cell} = files(cell).name;
    
    i = strfind(files(cell).name, 'tet');
    if strcmp(files(cell).name(i+4),'_')
        tetno = strread(files(cell).name(i+3)); % single character number
    else        
        tetno = strread(files(cell).name(i+3:i+4)); % number >10
    end
    tetnos(cell) = tetno;
    if tetno>8
        v1(cell) = logical(0); v2(cell) = logical(1); cellstr = 'V2';
    else
        v1(cell) = logical(1); v2(cell) = logical(0); cellstr = 'V1';
    end
    
    spike = result.waveforms(:,wvchan);
    interpspike = spline(1:32,spike,1:.1:32);
    [adiff(cell),swidth(cell),ptr(cell),eslope(cell)] = spikequant(interpspike);
    
    depth(cell) = result.depth;
    cellresp(cell,:,:) = resp;
    celllfpresp(cell,:,:) = lfpresp;
    cellchan(cell,:) = chan;
    
   
    msta = linspace(-prestim,trialdur+poststim,size(resp,2));
    
    printname = files(cell).name;
    printname(find(printname=='_')) = ' ';
    
   
    binwidth = 50;
    oris = unique(result.gratingInfo.Orientation); oris(find(oris == 180)) = [];
    lightlevels = unique(result.light);
    for l = 1:length(lightlevels)
        for ori = 1:length(oris)            
            thisinds = find((result.gratingInfo.Orientation == oris(ori) | result.gratingInfo.Orientation ==oris(ori)+180) &...
                result.light == lightlevels(l));
            condn(l,ori) = length(thisinds);
            condresp(l,ori,:) = nanmean(resp(thisinds,:),1);
            condresperr(l,ori,:) = nanstd(resp(thisinds,:),1,1)./sqrt(length(thisinds));
            condlfpspect(l,ori,:) = nanmean(lfpspect(thisinds,:),1);
            condlfpspecterr(l,ori,:) = nanstd(lfpspect(thisinds,:),1)./sqrt(length(thisinds));
            if ~isnan(condresp(l,ori,:))
                [bincondresp(l,ori,:),bta] = binit(condresp(l,ori,:),binwidth);
            else
                bincondresp(l,ori,:) = binit(condresp(l,ori,:),binwidth);
            end
            binconderr(l,ori,:) = binit(condresperr(l,ori,:),binwidth);
            cfr(l,ori) = nanmean(frs(thisinds));
            cerr(l,ori) =nanstd(frs(thisinds))./sqrt(length(thisinds));
            condspecgram(l,ori,:,:) = nanmean(S(thisinds,:,:),1);
            
            filly(l,ori,:) = [squeeze(condlfpspect(l,ori,1:104)+condlfpspecterr(l,ori,1:104))',fliplr(squeeze(condlfpspect(l,ori,1:104)-condlfpspecterr(l,ori,1:104))')];
        end
    end
    bincondresp = bincondresp.*(1000/binwidth);
    binconderr = binconderr.*(1000/binwidth);
    bta = bta-prestim;
    
    %fill
    fillx = [1:104,fliplr(1:104)];
    
    
    % find peaks    
    beta = [15,40];
    gamma = [50,70];
    b1 = find(fax>beta(1),1); b2 = find(fax>beta(2),1);
    g1 = find(fax>gamma(1),1); g2 = find(fax>gamma(2),1);
    bsig = squeeze(condlfpspect(1,2,b1:b2));
    gsig = squeeze(condlfpspect(1,1,g1:g2));
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
    %%%%%%%%%%%%%%%%%%%%%%
    
    anovavec = [cfr(:,1);cfr(:,2)];
    gr1 = ones(numel(cfr),1); gr1(1:6) = 0;
    gr2 = [lightlevels';lightlevels'];
    [p(cell,:),table,stats] = anovan(anovavec,{gr1 gr2},'display','off');
    
    
%     figure('name' ,['cell: ' int2str(cell), '  p vis: ' num2str(p(cell,1)) '  p light: ' num2str(p(cell,2))])
%     for i = 1:6
%         subplot(2,6,i)
%         boundedline(bta, squeeze(bincondresp(i,1,:)),squeeze(binconderr(i,1,:)));
%         hold on
%         line([500,500],[0,max(max(max(bincondresp)))+max(max(max(binconderr)))],'color','r')
% %         line([5500,5500],[0,max(max(max(bincondresp)))+max(max(max(binconderr)))],'color','r')
%         line([1500,1500],[0,max(max(max(bincondresp)))+max(max(max(binconderr)))],'color','r')
% %         axis([-300,8300,0,max(max(max(bincondresp)))+max(max(max(binconderr)))+1])
%         axis([-300,2300,0,max(max(max(bincondresp)))+max(max(max(binconderr)))+1])
%         title(['level: ' num2str(lightlevels(i))]);
%         
%         subplot(2,6,6+i)
%         boundedline(bta, squeeze(bincondresp(i,2,:)),squeeze(binconderr(i,2,:))); 
%         hold on
%         line([500,500],[0,max(max(max(bincondresp)))+max(max(max(binconderr)))],'color','r')
% %         line([5500,5500],[0,max(max(max(bincondresp)))+max(max(max(binconderr)))],'color','r') 
%         line([1500,1500],[0,max(max(max(bincondresp)))+max(max(max(binconderr)))],'color','r')       
% %         axis([-300,8300,0,max(max(max(bincondresp)))+max(max(max(binconderr)))+1])      
%         axis([-300,2300,0,max(max(max(bincondresp)))+max(max(max(binconderr)))+1])
%     end
    
    figure
    errorbar(lightlevels,cfr(:,1),cerr(:,1),'.')
    hold on
    errorbar(lightlevels,cfr(:,2),cerr(:,2),'r.')
    xlabel('lightintensity % max')
    ylabel('firing rate [Hz]')
    legend([{'no stimulus'},{'drifting grating'}])
    title(['cell: ' cellname{cell}, '  spikewidth: ' num2str(swidth(cell)) '  depth: ' num2str(depth(cell))])
    
    figure
    subplot(2,3,1)
    imagesc(t,f(1:50),log(squeeze(condspecgram(1,1,:,1:50))'));
    axis xy;
    caxis([0,4])
    
    subplot(2,3,2)
    imagesc(t,f(1:50),log(squeeze(condspecgram(2,1,:,1:50))'));
    axis xy;
    title('blank screen')
    caxis([0,4])
    
    subplot(2,3,3)
    imagesc(t,f(1:50),log(squeeze(condspecgram(3,1,:,1:50))'));
    axis xy;
    caxis([0,4])
    
    subplot(2,3,4)
    imagesc(t,f(1:50),log(squeeze(condspecgram(4,1,:,1:50))'));
    axis xy;
    caxis([0,4])
    
    subplot(2,3,5)
    imagesc(t,f(1:50),log(squeeze(condspecgram(5,1,:,1:50))'));
    axis xy;
    caxis([0,4])
    
    subplot(2,3,6)
    imagesc(t,f(1:50),log(squeeze(condspecgram(6,1,:,1:50))'));
    axis xy;
    caxis([0,4])
    
    
    figure
    subplot(2,3,1)
    imagesc(t,f(1:50),log(squeeze(condspecgram(1,2,:,1:50))'));
    axis xy;
    caxis([0,4])
    
    subplot(2,3,2)
    imagesc(t,f(1:50),log(squeeze(condspecgram(2,2,:,1:50))'));
    axis xy;
    caxis([0,4])
    title('drifting grating')
    
    subplot(2,3,3)
    imagesc(t,f(1:50),log(squeeze(condspecgram(3,2,:,1:50))'));
    axis xy;
    caxis([0,4])
    
    subplot(2,3,4)
    imagesc(t,f(1:50),log(squeeze(condspecgram(4,2,:,1:50))'));
    axis xy;
    caxis([0,4])
    
    subplot(2,3,5)
    imagesc(t,f(1:50),log(squeeze(condspecgram(5,2,:,1:50))'));
    axis xy;
    caxis([0,4])
    
    subplot(2,3,6)
    imagesc(t,f(1:50),log(squeeze(condspecgram(6,2,:,1:50))'));
    axis xy;
    caxis([0,4])
    
%     figure
%     semilogy(fax,squeeze(mean(condlfpspect(1,1,:),2)))
%     hold on
%     semilogy(fax,squeeze(mean(condlfpspect(2,1,:),2)),'c')
%     semilogy(fax,squeeze(mean(condlfpspect(3,1,:),2)),'g')
%     semilogy(fax,squeeze(mean(condlfpspect(4,1,:),2)),'y')
%     semilogy(fax,squeeze(mean(condlfpspect(5,1,:),2)),'m')
%     semilogy(fax,squeeze(mean(condlfpspect(6,1,:),2)),'r')
%     axis([0,180,0.3,500]);
% %     legend([{'0'},{'0.08'},{'0.1'},{'0.12'},{'0.16'}])
%     legend([{'0'},{'10'},{'18'},{'32'},{'56'},{'100'}])
    
    figure
    semilogy(fax,squeeze(condlfpspect(1,2,:)),'linewidth',2)
    hold on
    semilogy(fax,squeeze(condlfpspect(2,2,:)),'c','linewidth',2)
    semilogy(fax,squeeze(condlfpspect(3,2,:)),'g','linewidth',2)
    semilogy(fax,squeeze(condlfpspect(4,2,:)),'y','linewidth',2)
    semilogy(fax,squeeze(condlfpspect(5,2,:)),'m','linewidth',2)
    semilogy(fax,squeeze(condlfpspect(6,2,:)),'r','linewidth',2)
    axis([0,180,0.3,500]);
    legend([{'0'},{'0.1'},{'0.2'},{'0.4'},{'0.6'},{'0.8'}])
%     legend([{'0'},{'10'},{'18'},{'32'},{'56'},{'100'}])
%     legend({(num2str(lightlevels))'})
    title('drifting grating')    
    
    figure
    semilogy(fax,squeeze(condlfpspect(1,1,:)),'linewidth',2)
    hold on
    semilogy(fax,squeeze(condlfpspect(2,1,:)),'c','linewidth',2)
    semilogy(fax,squeeze(condlfpspect(3,1,:)),'g','linewidth',2)
    semilogy(fax,squeeze(condlfpspect(4,1,:)),'y','linewidth',2)
    semilogy(fax,squeeze(condlfpspect(5,1,:)),'m','linewidth',2)
    semilogy(fax,squeeze(condlfpspect(6,1,:)),'r','linewidth',2)
    axis([0,180,0.3,500]);
    legend([{'0'},{'0.1'},{'0.2'},{'0.4'},{'0.6'},{'0.8'}])
%     legend([{'0'},{'10'},{'18'},{'32'},{'56'},{'100'}])
%     legend({(num2str(lightlevels))'})
    title('gray blank screen')

    figure
    plot(lightlevels,condlfpspect(:,2,bpi),'o','markerfacecolor','b')
    title(['change in low gamma power at ' num2str(fax(bpi)) ' Hz depth ' num2str(depth(cell))])
    
    condfr(cell,:,:) = cfr; conderr(cell,:,:) = cerr;
    
    cell = cell+1;
    disp('');
    
end

%spike classification
kmeansind = kmeans([eslope',ptr',swidth',adiff'],2);
if mean(swidth(find(kmeansind==1)))<mean(swidth(find(kmeansind==2)))  %1 is FS
    pfs = find(kmeansind==1); prs = find(kmeansind==2); pfsv = kmeansind==1; prsv = kmeansind==2;
else
    pfs = find(kmeansind==2); prs = find(kmeansind==1); pfsv = kmeansind==2; prsv = kmeansind==1;
end

figure
plot(swidth(prs),adiff(prs),'b.')
xlabel('spike width')
ylabel('amplitude diff')
hold on
plot(swidth(pfs),adiff(pfs),'r.')
% axis([5.5,20.5,-.9,.7])

l23 = depth<375;
l4 = depth>=375&depth<=500;
l5 = depth>500&depth<=800;
l23rs = l23&prsv';
l23fs = l23&pfsv';
l4rs = l4&prsv';
l4fs = l4&pfsv';
l5rs = l5&prsv';
l5fs = l5&pfsv';

for i = 1:5
    condomi(:,i,:) = (condfr(:,i+1,:)-condfr(:,1,:))./(condfr(:,i+1,:)+condfr(:,1,:));
end

cond = l5rs;
figure
hold on
for i = find(cond)
    for j = 1:5
        plot(j,condomi(i,j,1),'co')
    end
end
for i = find(cond)
    plot(condomi(i,:,1),'c')
end
axis([.5,5.5,-1.1,1.1])
line([.5,5.5],[0,0],'color','r')

disp('')
