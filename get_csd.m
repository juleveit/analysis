function csd = get_csd(animalid, block, channels)

if nargin == 0
    animalid = '180626';
    block = 12;
    channels = 16;
    allchannels = [32];
    electrode = 1;
    depth = 600;
    penangle = 65;
end

ellength = 25*(channels-1);
highest = depth-ellength;
depthax = depth:-25:highest;
depthax = depthax.*(cosd(90-penangle)*cosd(22));
depth = depth*(cosd(90-penangle)*cosd(22));

printyn = 1;
singleunit = 1;
dolightstuff = 1;
recalculate = 0;

basepath = ['C:\Users\Julia\work\data\' animalid '\'];
if printyn
    if ~exist([basepath, 'pdfs/'],'dir'),mkdir([basepath, 'pdfs/']),end
    if ~exist([basepath, 'pdfs/' 'CSD/'],'dir'),mkdir([basepath, 'pdfs/' 'CSD/']),end
end
printpath = [basepath 'pdfs\CSD\'];

if singleunit
    filepath = ['C:\Users\Julia\work\data\' animalid '\singleunits\'];
    % filepath = ['C:\Users\Julia\work\data\' animalid '\multiunits\'];

    recfiles = dir([filepath, animalid '_block' int2str(block) '*.mat']);

    for tet = 1:channels/4
        filebeg = [animalid '_block' int2str(block) '_tet' int2str(tet)];
        clear file;
        for i = 1:length(recfiles)
            if strfind(recfiles(i).name, filebeg)
                file = recfiles(i).name;
                break;
            else
                continue;
            end
        end
        load([filepath, file]);
        lfp(:,(tet-1)*4+1:tet*4) = result.lfp;
    end
else
    pcs = 0;
    if electrode>1        
        for i = 1:electrode-1
            pcs = allchannels(i)+pcs; % previous channels
        end
    end
    evalchans = pcs+1:pcs+allchannels(electrode);
    if ~exist([basepath 'csdresult_' int2str(block) '.mat']) || recalculate
        result = lfpdataprepare(basepath,animalid,block,evalchans);   
        save([basepath 'csdresult_' int2str(block) '.mat'], 'result')
    else
        load([basepath 'csdresult_' int2str(block) '.mat']);
    end
    lfp = result.lfp';
end

if strcmp(result.protocol, 'grating')
    prestim = 300;
    poststim = 300;
    result.stimduration = result.stimduration.*1000;
    timeax = -prestim+1:result.stimduration+poststim;
else
    if strcmp(result.protocol, 'onoff')|strcmp(result.protocol, 'hartley')
        result.stimduration = 500;
    elseif strcmp(result.protocol, 'motor')
        result.stimduration = 6000;
    else
        result.stimduration = round((1000/60)*result.FramesperStim);;
    end
    prestim = 0; poststim = 0;  
    timeax = 1:result.stimduration;
end


sr = 1000;
gamma = eegfilt(lfp',sr,30,90);
lowgamma = eegfilt(lfp',sr,20,40);
for i = 1:channels
    h = hilbert(gamma(i,:)); gpow(i,:) = abs(h);
    h1 = hilbert(lowgamma(i,:)); g1pow(i,:) = abs(h1);
end

% if length(result.msstamps)~=length(result.light)
%         disp('');
%     %     result.msstamps([54,201,239,316]) = []; % for 140807 block 7
% %         result.msstamps(167) = [];
% %         result.msstamps(92) = [];
% %         result.msstamps([303]) = []; % for 150407 block 5
% %         result.msstamps([190,233]) = []; % for 150219 - 6
% %         result.msstamps([303]) = []; % 150407 - 5
% %         save([basepath 'csdresult_' int2str(block) '.mat'], 'result')
% %         save([filepath, file],'result')
% %         pause;
% end
    
for i = 1:length(result.msstamps)
    lfpresp(i,:,:) = lfp(result.msstamps(i)-prestim+1:result.msstamps(i)+result.stimduration+poststim,:);
    gammaresp(i,:,:) = gpow(:,result.msstamps(i)-prestim+1:result.msstamps(i)+result.stimduration+poststim)';
    lowgammaresp(i,:,:) = g1pow(:,result.msstamps(i)-prestim+1:result.msstamps(i)+result.stimduration+poststim)';
end
% for i = find(result.light)
%     lfpresp(i,:,:) = lfp(result.msstamps(i)+500:result.msstamps(i)+1000,:);
%     gammaresp(i,:,:) = gpow(:,result.msstamps(i)+500:result.msstamps(i)+1000)';
% end

if isfield(result,'light') & find(result.light) & dolightstuff
    
    window = 900+prestim:1411+prestim;
    longwin = prestim:prestim+2000;
    % trialnfft = 2^nextpow2(500);
    % trialfax = sr/2*linspace(0,1,trialnfft/2+1);
    for chan = 1:size(lfpresp,3)
        for i = 1:size(lfpresp,1)
    %         y = fft(lfpresp(i,window,chan),trialnfft);
    %         lfpspect(chan,i,:) = abs(y(1:trialnfft/2+1));
            [lfpspect(chan,i,:),trialfax] = pmtm(lfpresp(i,window,chan),3,[],sr);
            [lfpspectlong(chan,i,:),trialfaxlong] = pmtm(lfpresp(i,longwin,chan),2,[],sr);
        end    
    end
    
    if length(unique(result.gratingInfo.Contrast))  == 1 % size block
        vari = result.gratingInfo.size; contr = 0; condstr = 'SIZE';
    else
        vari = result.gratingInfo.Contrast; contr = 1; condstr = 'CONTRAST'; % contrast block
    end
    if contr, str0 = 'low', str1 = 'high'; else str0 = 'small'; str1 = 'large'; end

    oris = unique(result.gratingInfo.Orientation); oris(find(oris == -1)) = [];
    clevels = unique(vari); clevels(find(clevels == 0)) = []; %delete control condition
    for l = 1:length(unique(result.light))
        for ori = 1:length(oris)
            for co = 1:length(clevels)
                thisinds = find(result.gratingInfo.Orientation == oris(ori) &...
                    vari == clevels(co) & result.light == l-1);
                condlfpspect(:,l,ori,co,:) = nanmean(lfpspect(:,thisinds,:),2);
                condlfpspectlong(:,l,ori,co,:) = nanmean(lfpspectlong(:,thisinds,:),2);
                condlfpspecterr(:,l,ori,co,:) = nanstd(lfpspect(:,thisinds,:),1,2)./sqrt(length(thisinds));
                condlfpresp(:,l,ori,co,:) = squeeze(nanmean(lfpresp(thisinds,:,:),1))';
                condgammaresp(:,l,ori,co,:) = squeeze(nanmean(gammaresp(thisinds,:,:),1))';
                condlowgammaresp(:,l,ori,co,:) = squeeze(nanmean(lowgammaresp(thisinds,:,:),1))';
            end
        end
    end
    
    contindsl0 = find(vari == 0 & result.light == 0);
    contindsl1 = find(vari == 0 & result.light == 1);
    controllfpspect(:,1,:) = nanmean(lfpspect(:,contindsl0,:),2);
    controllfpspect(:,2,:) = nanmean(lfpspect(:,contindsl1,:),2);
    controllfpspecterr(:,1,:) = nanstd(lfpspect(:,contindsl0,:),1,2)./sqrt(length(contindsl0));
    controllfpspecterr(:,2,:) = nanstd(lfpspect(:,contindsl1,:),1,2)./sqrt(length(contindsl0));
    controllfpspectlong(:,1,:) = nanmean(lfpspectlong(:,contindsl0,:),2);
    controllfpspectlong(:,2,:) = nanmean(lfpspectlong(:,contindsl1,:),2);
                
    % l0p = find(result.light == 0); l1p = find(result.light);
    l0inds = find(result.light == 0 & result.gratingInfo.Contrast == 1); 
    l1inds = find(result.light & result.gratingInfo.Contrast == 1); 
    for i = 1:size(lfpspect,1)
        figure
        semilogy(trialfax,squeeze(mean(lfpspect(i,l1inds,:),2)),'r','linewidth',2)
        hold on
        semilogy(trialfax,squeeze(mean(lfpspect(i,l0inds,:),2)),'linewidth',2)
        semilogy(trialfax,squeeze(mean(lfpspect(i,l1inds,:),2))-(squeeze(std(lfpspect(i,l1inds,:),1,2))./sqrt(length(l1inds))),'r')
        semilogy(trialfax,squeeze(mean(lfpspect(i,l1inds,:),2))+(squeeze(std(lfpspect(i,l1inds,:),1,2))./sqrt(length(l1inds))),'r')
        semilogy(trialfax,squeeze(mean(lfpspect(i,l0inds,:),2))-(squeeze(std(lfpspect(i,l0inds,:),1,2))./sqrt(length(l0inds))))
        semilogy(trialfax,squeeze(mean(lfpspect(i,l0inds,:),2))+(squeeze(std(lfpspect(i,l0inds,:),1,2))./sqrt(length(l0inds))))
        axis([0,150,...
            min([min(squeeze(mean(lfpspect(i,l0inds,1:53),2))),min(squeeze(mean(lfpspect(i,l1inds,1:53),2)))]),...
            max([max(squeeze(mean(lfpspect(i,l0inds,1:53),2))),max(squeeze(mean(lfpspect(i,l1inds,1:53),2)))])])
    end
    
    el = 7;
    figure
    semilogy(trialfax,squeeze(mean(condlfpspect(el,1,:,1,:))),'linewidth',2)
    hold on
    semilogy(trialfax,squeeze(mean(condlfpspect(el,2,:,1,:))),'r','linewidth',2)
    semilogy(trialfax,squeeze(mean(condlfpspect(el,1,:,5,:))),'c','linewidth',2)
    semilogy(trialfax,squeeze(mean(condlfpspect(el,2,:,5,:))),'m','linewidth',2)
%     semilogy(trialfax,squeeze(controllfpspect(el,1,:)),'color',[.5,.5,.5],'linewidth',2);
%     semilogy(trialfax,squeeze(controllfpspect(el,2,:)),'k','linewidth',2);
%     semilogy(trialfax,squeeze(mean(condlfpspect(el,1,:,1,:)))+squeeze(mean(condlfpspecterr(el,1,:,1,:))),'b');
%     semilogy(trialfax,squeeze(mean(condlfpspect(el,1,:,1,:)))-squeeze(mean(condlfpspecterr(el,1,:,1,:))),'b');
%     semilogy(trialfax,squeeze(mean(condlfpspect(el,2,:,1,:)))+squeeze(mean(condlfpspecterr(el,2,:,1,:))),'r');
%     semilogy(trialfax,squeeze(mean(condlfpspect(el,2,:,1,:)))-squeeze(mean(condlfpspecterr(el,2,:,1,:))),'r');
%     semilogy(trialfax,squeeze(mean(condlfpspect(el,1,:,5,:)))+squeeze(mean(condlfpspecterr(el,1,:,5,:))),'c');
%     semilogy(trialfax,squeeze(mean(condlfpspect(el,1,:,5,:)))-squeeze(mean(condlfpspecterr(el,1,:,5,:))),'c');
%     semilogy(trialfax,squeeze(mean(condlfpspect(el,2,:,5,:)))+squeeze(mean(condlfpspecterr(el,2,:,5,:))),'m');
%     semilogy(trialfax,squeeze(mean(condlfpspect(el,2,:,5,:)))-squeeze(mean(condlfpspecterr(el,2,:,5,:))),'m');
    axis([0,150,.1,1800])
    legend([{[str0 ' L0']},{[str0 ' L1']},{[str1 ' L0']},{[str1 ' L1']}],'location','ne')
    title('depth: ~250')
    
    %gamma power
    gwin = 600+prestim:800+prestim;
    figure
    plot(squeeze(mean(mean(gammaresp(l1inds,gwin,:),1),2))-squeeze(mean(mean(gammaresp(l0inds,gwin,:),1),2)),depthax,'ko','markerfacecolor','k')
    axis ij
    ax = axis;
    line([ax(1),ax(2)],[375,375],'color','k','linestyle',':');
    line([ax(1),ax(2)],[500,500],'color','k','linestyle',':');
    line([ax(1),ax(2)],[800,800],'color','k','linestyle',':');
    line([0,0],[0,1000],'color','k','linewidth',2);
    
    figure
    subplot(2,2,1)
    imagesc(trialfax,depthax,log(squeeze(nanmean(condlfpspect(:,1,:,1,:),3))))
    axis([-.5,100.5,depthax(end),depthax(1)])
    caxis([0,7])
    colorbar
    title([' spectrum in depth light OFF ' str0])

    subplot(2,2,2)
    imagesc(trialfax,depthax,log(squeeze(nanmean(condlfpspect(:,2,:,1,:),3))))
    axis([-.5,100.5,depthax(end),depthax(1)])
    caxis([0,7])
    colorbar
    title([' spectrum in depth light ON ' str0])
    
    subplot(2,2,3)
    imagesc(trialfax,depthax,log(squeeze(nanmean(condlfpspect(:,1,:,5,:),3))))
    axis([-.5,100.5,depthax(end),depthax(1)])
    caxis([0,7])
    colorbar
    title([' spectrum in depth light OFF ' str1])

    subplot(2,2,4)
    imagesc(trialfax,depthax,log(squeeze(nanmean(condlfpspect(:,2,:,5,:),3))))
    caxis([0,7])
    axis([-.5,100.5,depthax(end),depthax(1)])
    colorbar
    title([' spectrum in depth light ON ' str1])
    
    figure
    subplot(2,2,1)
    imagesc(trialfax,depthax,log(squeeze(nanmean(condlfpspect(:,2,:,1,:),3)))-log(squeeze(nanmean(condlfpspect(:,1,:,1,:),3))));
    caxis([-2,2])
    axis([-.5,100.5,depthax(end),depthax(1)])
    colorbar
    title(['light ON-OFF ' str0]);   
    
    subplot(2,2,2)
    imagesc(trialfax,depthax,log(squeeze(nanmean(condlfpspect(:,2,:,5,:),3)))-log(squeeze(nanmean(condlfpspect(:,1,:,5,:),3))));
    caxis([-2,2])
    axis([-.5,100.5,depthax(end),depthax(1)])
    colorbar
    title(['light ON-OFF ' str1])   
    
    subplot(2,2,3)
    imagesc(trialfax,depthax,log(squeeze(nanmean(condlfpspect(:,1,:,5,:),3)))-log(squeeze(nanmean(condlfpspect(:,1,:,1,:),3))));
    caxis([-2,2])
    axis([-.5,100.5,depthax(end),depthax(1)])
    colorbar
    title(['light OFF ' str1 ' - ' str0]);  
    
    subplot(2,2,4)
    imagesc(trialfax,depthax,log(squeeze(nanmean(condlfpspect(:,2,:,5,:),3)))-log(squeeze(nanmean(condlfpspect(:,2,:,1,:),3))));
    caxis([-2,2])
    axis([-.5,100.5,depthax(end),depthax(1)])
    colorbar
    title(['light ON ' str1 ' - ' str0]);  
    
    %find peak of the beta and gamma rhythms in specified frequency bands
    el = find(depthax<400,1);
    beta = [15,40];
    gamma = [50,70];
    highgamma = [70,120];
    b1 = find(trialfax>beta(1),1); b2 = find(trialfax>beta(2),1);
    g1 = find(trialfax>gamma(1),1); g2 = find(trialfax>gamma(2),1);
    hg1 = find(trialfax>highgamma(1),1); hg2 = find(trialfax>highgamma(2),1);
%     [bp,bpi] = max(squeeze(nanmean(condlfpspect(el,1,:,5,b1:b2),3)),[],2);
%     bpi = bpi+b1-1;    
%     [gp,gpi] = max(squeeze(nanmean(condlfpspect(el,1,:,1,g1:g2),3)),[],2);
% %     [gp,gpi] = max(squeeze(controllfpspect(:,1,g1:g2)),[],2);
%     gpi = gpi+g1-1;
    bsig = squeeze(nanmean(condlfpspect(el,1,:,5,b1:b2),3));
    gsig = squeeze(nanmean(condlfpspect(el,1,:,1,g1:g2),3));
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
        
    figure
    subplot(2,2,1)    
    semilogy(trialfax,squeeze(mean(condlfpspect(el,1,:,1,:))),'linewidth',2)
    hold on
    semilogy(trialfax,squeeze(mean(condlfpspect(el,2,:,1,:))),'r','linewidth',2)
    semilogy(trialfax,squeeze(mean(condlfpspect(el,1,:,5,:))),'c','linewidth',2)
    semilogy(trialfax,squeeze(mean(condlfpspect(el,2,:,5,:))),'m','linewidth',2)
    axis([0,150,.1,1800])
    legend([{[str0 ' L0']},{[str0 ' L1']},{[str1 ' L0']},{[str1 ' L1']}],'location','ne')
    title(['depth: ' int2str(round(depthax(el)))])
    line([beta(1),beta(1)],[.1,1800],'color','k','linestyle','--')
    line([beta(2),beta(2)],[.1,1800],'color','k','linestyle','--')
    line([gamma(1),gamma(1)],[.1,1800],'color','k','linestyle','--')
    line([gamma(2),gamma(2)],[.1,1800],'color','k','linestyle','--')
    line([highgamma(1),highgamma(1)],[.1,1800],'color','k','linestyle','--')
    line([highgamma(2),highgamma(2)],[.1,1800],'color','k','linestyle','--')
    line([trialfax(gpi),trialfax(gpi)],[.1,1800],'color','k')
    line([trialfax(bpi),trialfax(bpi)],[.1,1800],'color','k')
    
    subplot(2,2,2)
    errorbar([0,clevels],[controllfpspect(el,1,bpi), squeeze(nanmean(condlfpspect(el,1,:,:,bpi),3))'],[controllfpspecterr(el,1,bpi), squeeze(nanmean(condlfpspecterr(el,1,:,:,bpi),3))'],'o-')
    hold on
    errorbar([0,clevels],[controllfpspect(el,2,bpi), squeeze(nanmean(condlfpspect(el,2,:,:,bpi),3))'],[controllfpspecterr(el,2,bpi), squeeze(nanmean(condlfpspecterr(el,2,:,:,bpi),3))'],'ro-')
    ylabel(['beta power'])
    title(['beta power: ' int2str(round(trialfax(bpi))) ' Hz'])
    xlabel(condstr)
    
    subplot(2,2,3)
    errorbar([0,clevels], [controllfpspect(el,1,gpi), squeeze(nanmean(condlfpspect(el,1,:,:,gpi),3))'],[controllfpspecterr(el,1,gpi), squeeze(nanmean(condlfpspecterr(el,1,:,:,gpi),3))'],'o-')
    hold on
    errorbar([0,clevels], [controllfpspect(el,1,gpi), squeeze(nanmean(condlfpspect(el,2,:,:,gpi),3))'],[controllfpspecterr(el,1,gpi), squeeze(nanmean(condlfpspecterr(el,2,:,:,gpi),3))'],'ro-')
    ylabel(['gamma power'])
    title(['gamma power: ' int2str(round(trialfax(gpi))) ' Hz'])
    xlabel(condstr)
    
    subplot(2,2,4)
    errorbar([0,clevels], [nanmean(controllfpspect(el,1,hg1:hg2),3), squeeze(nanmean(nanmean(condlfpspect(el,1,:,:,hg1:hg2),3),5))'],[nanmean(controllfpspecterr(el,1,hg1:hg2),3), squeeze(nanmean(nanmean(condlfpspecterr(el,1,:,:,hg1:hg2),3),5))'],'o-');
    hold on
    errorbar([0,clevels], [nanmean(controllfpspect(el,2,hg1:hg2),3), squeeze(nanmean(nanmean(condlfpspect(el,2,:,:,hg1:hg2),3),5))'],[nanmean(controllfpspecterr(el,2,hg1:hg2),3), squeeze(nanmean(nanmean(condlfpspecterr(el,2,:,:,hg1:hg2),3),5))'],'ro-');
    ylabel(['high gamma power'])
    title(['high gamma power: ' int2str(highgamma(1)) '-' int2str(highgamma(2)) ' Hz'])
    xlabel(condstr)
    
    figure
    imagesc(trialfax,depthax,log(squeeze(mean(lfpspect(:,l0inds,:),2))))
%     axis xy
    axis([-.5,100.5,depthax(end),depthax(1)])
    colorbar
    xlabel('frequency [Hz]')
    ylabel(['contactno'])
    title('LFP spectra in depth, light OFF')
    caxis([0,5]);
    if printyn
        figSize = [30 21];
        set(gcf,'paperunits','centimeters','papersize',figSize,'paperposition',[0 0 figSize])
        print(gcf,[printpath ,  '01_L0SpectDepth.pdf'], '-dpdf' );
    end


    figure
    imagesc(trialfax,depthax,log(squeeze(mean(lfpspect(:,l1inds,:),2))))
    axis([-.5,150.5,depthax(end),depthax(1)])
    colorbar
    xlabel('frequency [Hz]')
    ylabel(['contactno'])
    title('LFP spectra in depth, light ON')
    caxis([0,5])
    if printyn
        figSize = [30 21];
        set(gcf,'paperunits','centimeters','papersize',figSize,'paperposition',[0 0 figSize])
        print(gcf,[printpath ,  '02_L1SpectDepth.pdf'], '-dpdf' );
    end

    figure
    imagesc(trialfax,depthax,log(squeeze(mean(lfpspect(:,l1inds,:),2)))-log(squeeze(mean(lfpspect(:,l0inds,:),2))))
    axis([-.5,150.5,depthax(end),depthax(1)])
    colorbar
    xlabel('frequency [Hz]')
    ylabel(['contactno'])
    title('LFP spectra in depth, light ON - light OFF')
    caxis([-1,1])
    if printyn
        figSize = [30 21];
        set(gcf,'paperunits','centimeters','papersize',figSize,'paperposition',[0 0 figSize])
        print(gcf,[printpath ,  '03_L1-L0SpectDepth.pdf'], '-dpdf' );
    end

    mlfpresp = squeeze(mean(lfpresp(l0inds,:,:)));
else
    mlfpresp = squeeze(mean(lfpresp,1));
end

% %spatial smoothing either over one neighbor channel on each side
% %first duplicate uppermost and lowermost channel
% nmlfpresp(:,1) = mlfpresp(:,1);
% nmlfpresp(:,2:channels+1) = mlfpresp;
% nmlfpresp(:,channels+2) = mlfpresp(:,channels);
% %then take weighted average of channel and neighboring channels
% for i = 1:channels
%     smlfpresp(:,i) = (nmlfpresp(:,i)+nmlfpresp(:,i+2)+(2.*nmlfpresp(:,i+1)))./4;
% end

% % or smooth over 2 adjacant
% nmlfpresp(:,1) = mlfpresp(:,1); nmlfpresp(:,2) = mlfpresp(:,1);
% nmlfpresp(:,3:channels+2) = mlfpresp;
% nmlfpresp(:,channels+3) = mlfpresp(:,channels); nmlfpresp(:,channels+4) = mlfpresp(:,channels);
% for i = 1:channels
%     smlfpresp(:,i) = (nmlfpresp(:,i)+nmlfpresp(:,i+4)+(2.*nmlfpresp(:,i+1))+2.*nmlfpresp(:,i+3)...
%         +4.*nmlfpresp(:,i+2))./10;
% end

%or smooth over 3 adjacant
nmlfpresp(:,1) = mlfpresp(:,1); nmlfpresp(:,2) = mlfpresp(:,1); nmlfpresp(:,3) = mlfpresp(:,3);
nmlfpresp(:,4:channels+3) = mlfpresp;
nmlfpresp(:,channels+4) = mlfpresp(:,channels); nmlfpresp(:,channels+5) = mlfpresp(:,channels); nmlfpresp(:,channels+6) = mlfpresp(:,channels);
for i = 1:channels
    smlfpresp(:,i) = (nmlfpresp(:,i)+nmlfpresp(:,i+6)+(2.*nmlfpresp(:,i+1))+2.*nmlfpresp(:,i+5)...
        +(3.*nmlfpresp(:,i+2))+3.*nmlfpresp(:,i+4)+4.*nmlfpresp(:,i+3))./16;
end

%compute actual csd
for i = 1:channels-2
    csd(:,i) = (smlfpresp(:,i)-(2.*smlfpresp(:,i+1))+smlfpresp(:,i+2));
end

% %or over two sites away
% for i = 1:channels-4
%     csd(:,i) = smlfpresp(:,i)+smlfpresp(:,i+4)-(2.*smlfpresp(:,i+2));
% end


if isfield(result,'light') & find(result.light) & dolightstuff
    figure
    imagesc(timeax,depthax,squeeze(mean(lowgammaresp(l1inds,:,:)))')
%     axis([-.5,100.5,.5,32.5])
    colorbar
    xlabel('time [ms]')
    ylabel('contact no.')
    title('Gamma power in time and depth')
    if printyn
        figSize = [30 21];
        set(gcf,'paperunits','centimeters','papersize',figSize,'paperposition',[0 0 figSize])
        print(gcf,[printpath ,  '04_GammaDepth.pdf'], '-dpdf' );
    end
end


figure
imagesc(timeax,depthax,smlfpresp')
colorbar
caxis([-150,150])
axis([-.5,100.5,depthax(end),depthax(1)])
if printyn
    figSize = [30 21];
    set(gcf,'paperunits','centimeters','papersize',figSize,'paperposition',[0 0 figSize])
    print(gcf,[printpath ,  '05_AvgLFPDepth.pdf'], '-dpdf' );
end

% get the depths right:
highestcsd = highest + 25; %because one gets lost for CSD calculation
lowestcsd = depth - 25;
depthaxcsd = lowestcsd:-25:highestcsd;

figure
imagesc(timeax,depthaxcsd,csd')
colorbar
axis([0,200,depthaxcsd(end),depthaxcsd(1)])
if printyn
    figSize = [30 21];
    set(gcf,'paperunits','centimeters','papersize',figSize,'paperposition',[0 0 figSize])
    print(gcf,[printpath ,  '06_CSD.pdf'], '-dpdf' );
end

disp('');