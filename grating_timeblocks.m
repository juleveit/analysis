function grating_timeblocks

animalid = '150422';
blocks = [2:24];
penangle = 25;

lcol = 'c';

onlymod = 0;
sfc = 0;

basepath = ['C:\Users\Julia\work\data\' animalid '\'];
supath = [basepath 'singleunits\'];
basename = [animalid '_block' int2str(blocks(1)) '_tet'];

files = dir([supath, basename, '*.mat']);

prestim = 300;
poststim = 300;
respwin = 1:2000; % after stimulus onset
respwin = respwin+prestim;

cell = 1;
for cl = 1:length(files)
    
%     if strfind(files(cl).name, 'MU')
%         continue;
%     end
    
    for blck = 1:length(blocks)
        
        basename = [animalid '_block' int2str(blocks(blck)) '_tet'];
        SUfiles = dir([supath, basename, '*.mat']);
    
        load([supath, SUfiles(cl).name]);
    
        % calc spiking to see if includable
        msStimes = round(result.spikes);
        if ~isempty(msStimes) & msStimes(1) == 0, msStimes(1) = 1; end

        chan = zeros(1,length(result.lfp));
        chan(msStimes) = 1;

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


        trialdur = result.stimduration*1000;%+1000; %+1000 only for 150413
        msstamps = result.msstamps;
        if length(msstamps)~=length(result.light)
%             disp('');
    %         msstamps([54,201,239,316]) = []; % for 140807 block 7
    %         msstamps(200) = []; % for 140815 block 7
    %         msstamps(123) = []; % for 140815 block 7
    %         msstamps(385) = []; % for 140703 block 5
    %         msstamps(37) = []; % for 150113 block 7
%             msstamps(39) = []; % for 150113 block 11
%             msstamps(18) = []; % for 150414 block 7
%             msstamps(16) = []; % for 150414 block 10
%             msstamps(71) = []; % for 150417 block 2
%             msstamps(14) = []; % for 150417 block 2
%             result.msstamps = msstamps;
%             save([supath, SUfiles(cl).name],'result');
            pause;
        end
        
    %     trialnfft = 2^nextpow2(800);
    %     trialfax = sr/2*linspace(0,1,trialnfft/2+1);
        for i = 1:length(msstamps)
            resp(i,:) = chan(msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim); 
            lfpresp(i,:) = result.lfp(msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim);
    %         y = fft(lfpresp(i,1001:1800),trialnfft);
    %         lfpspect(i,:) = abs(y(1:trialnfft/2+1));
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
        nonoktrials = 1:size(resp,1); nonoktrials(oktrials) = [];
        stilltrials = 1:size(resp,1); stilltrials(notstill) = [];

        frs = sum(resp(:,respwin),2)./(length(respwin)/1000);
        bl = sum(resp(:,1:prestim),2)./(prestim/1000);
        gp = mean(gammaresp(:,respwin),2);
        blgp = mean(gammaresp(:,1:prestim),2);

        %determine if cell is visually modulated
        blfr = sum(resp(:,1:prestim),2);
        vrfr = sum(resp(:,prestim+40:2*prestim+40-1),2);
        vismod(cell,blck) = ttest2(blfr,vrfr);
        visdriven(cell,blck) = mean(vrfr)>=mean(blfr)+2; % average firing rate is increase at least 2Hz above baseline

        cellname{cell,blck} = files(cell).name;

        i = strfind(files(cell).name, 'tet');
        if strcmp(files(cell).name(i+4),'_')
            tetno = strread(files(cell).name(i+3)); % single character number
        else        
            tetno = strread(files(cell).name(i+3:i+4)); % number >10
        end
        tetnos(cell) = tetno;

        spike = result.waveforms(:,wvchan); 
        interpspike = spline(1:32,spike,1:.1:32);
        [adiff(cell,blck),swidth(cell,blck)] = spikequant(interpspike);

        % phases
        tmp = zeros(size(gphaseresp));
        tmp(find(resp)) = gphaseresp(find(resp));
        phases{cell,blck} = tmp(find(tmp));

        phaser(cell,blck) = circ_r(phases{cell});
        cmean(cell,blck) = circ_mean(phases{cell});

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
                allphaser(cell,blck,i) = circ_r(allphases{i}');
                allcmean(cell,blck,i) = circ_mean(allphases{i}');
            end
        end

        depth(cell,blck) = result.depth;
        cellresp(cell,blck,:,:) = resp;
        celllfpresp(cell,blck,:,:) = lfpresp;
%         cellchan(cell,blck,:) = chan;

%         figure
%         subplot(2,2,1)
%         semilogy(trialfax,mean(lfpspect),'linewidth',2);
%         semilogy(trialfax,mean(lfpspect)-(std(lfpspect)./sqrt(size(lfpspect,1))));
%         semilogy(trialfax,mean(lfpspect)+(std(lfpspect)./sqrt(size(lfpspect,1))));
%         axis([0,120,...
%             min([min(squeeze(mean(lfpspect(:,1:125)))),min(squeeze(mean(lfpspect(:,1:125))))]),...
%             max([max(squeeze(mean(lfpspect(:,1:125)))),max(squeeze(mean(lfpspect(:,1:125))))])])
%         xlabel('frequency [Hz]')
%         ylabel('spectral power')
%         title('LFP spectrum during light on vs off')
% 
%         subplot(2,2,2)
%         plot(mean(gammaresp))
%         title(['gamma power in time depth: ' int2str(depth(cell))])
% 
%         subplot(2,2,3)
%         [to,ro] = rose(phases{cell});
%         polar(to,ro,'b')
%         title(['gamma phase locking of unit ' int2str(cell) ' spikewidth: ' int2str(swidth(cell))])
% 
%         subplot(2,2,4)
%         plot(mean(lfpresp));

        msta = linspace(-prestim,trialdur+poststim,size(resp,2));

        baseline(cell,blck) = mean(bl);
        baselineerr = std(bl)./(sqrt(size(bl,1)));

        binwidth = 50;
        oris = unique(result.gratingInfo.Orientation); oris(find(oris == -1)) = [];
        clevels = unique(result.gratingInfo.Contrast); clevels(find(clevels == 0)) = []; %delete control condition
        
        for ori = 1:length(oris)
            thisinds = find(result.gratingInfo.Orientation == oris(ori));
            condn(ori) = length(thisinds);
            condresp(ori,:) = nanmean(resp(thisinds,:),1);
            condresperr(ori,:) = nanstd(resp(thisinds,:),1,1)./sqrt(length(thisinds));
            if ~isnan(condresp(ori,:))
                [bincondresp(ori,:),bta] = binit(condresp(ori,:),binwidth);
            else
                bincondresp(ori,:) = binit(condresp(ori,:),binwidth);
            end
            binconderr(ori,:) = binit(condresperr(ori,:),binwidth);
            condfr(cell,blck,ori) = nanmean(frs(thisinds));
            conderr(cell,blck,ori) =nanstd(frs(thisinds))./sqrt(length(thisinds));
            condgp(cell,blck,ori) = nanmean(gp(thisinds));
            condgperr(cell,blck,ori) = nanstd(gp(thisinds))./sqrt(length(thisinds));
            condtrials(cell,blck,ori,:) = frs(thisinds);
        end

        bincondresp = bincondresp.*(1000/binwidth);
        bta = bta-prestim;

        [oriprefratio(cell,blck), dirprefratio(cell,blck), prefori(cell,blck), meanori(cell,blck), osi(cell,blck), meandir(cell,blck), dsi(cell,blck)] = getOSI(squeeze(condfr(cell,blck,:))',oris);

        oneorifr = mean(reshape(squeeze(condfr(cell,blck,:)),4,2),2);
%         oneorifr = squeeze(condfr(cell,blck,:)); %only for 150413
        prefori = find(oneorifr == max(oneorifr),1);
        [preffr(cell,blck), prefdir(cell,blck)] = max(condfr(cell,blck,:),[],3);
        ortho = mod(prefori+2,length(oris)/2); if ortho == 0, ortho = length(oris)/2; end
        nulldir(cell,blck) = prefdir(cell,blck)-length(oris)/2; 
        if nulldir(cell,blck)<1, nulldir(cell,blck) = nulldir(cell,blck)+length(oris); end
        dosdsi(cell,blck) = (condfr(cell,blck,prefdir(cell,blck))-condfr(cell,blck,nulldir(cell,blck)))./(condfr(cell,blck,prefdir(cell,blck))+condfr(cell,blck,nulldir(cell,blck)));

        [binpref,bta] = binit(squeeze(condresp(prefdir(cell,blck),:)),binwidth); binpref = binpref.*(1000/binwidth);

        [binavg(cell,blck,:),bta] = binit(mean(resp),binwidth); binavg(cell,blck,:) = binavg(cell,blck,:).*(1000/binwidth); % in Hz
        binstd(cell,blck,:) = binit(std(resp)./sqrt(size(resp,1)),binwidth);
        ta = bta-prestim;
        alfpresp(cell,blck,:) = squeeze(mean(lfpresp));

        
%         polar(deg2rad([oris,0]),[squeeze(condfr),condfr(1)])

    %     cresp = [controlfr(cell), squeeze(mean(condfr))];
    %     xlevels = [0,clevels];
    %     
    %     subplot(2,2,3)
    %     errorbar(xlevels,cresp,[controlerr, squeeze(mean(conderr))],'ko','markersize',8)
    %     hold on
    %     xlabel('shown contrast')
    %     ylabel('Firing rate [Hz]')
    %     set(gca,'xtick',clevels); %[0, clevels])
    %     ax = axis;
    %     axis([-.1,1.1,-.1,ax(4)])
    %     title(['contrast response all orientations ']);

    %     params = fit_crf_NR(xlevels,cresp);
    %     plotx = [0:0.01:1];
    %     plot(plotx,NakaRushton(params,plotx),'k','linewidth',2);
    %     
    %     rmax(cell) = params(1); 
    %     c50(cell) = params(3); 
    %     r0(cell) = params(4);
    %     
%         subplot(2,4,7)
%         plot(spike,'linewidth',2)
%         axis([0,40,-100,100])
%         legend(['width: ' int2str(swidth(cell)) ' adiff: ' num2str(adiff(cell))])

    %     subplot(2,4,8)
    %     plot(ta,binctrl,'linewidth',2)
    %     axis([0,2500,min([min(binctrl)]),max([max(binctrl),0.1])])

        % running figure
        runfr = mean(frs(oktrials));
        norunfr = mean(frs(stilltrials));
        runfrerr = std(frs(oktrials))./sqrt(length(oktrials));
        norunfrerr = std(frs(stilltrials))./sqrt(length(stilltrials));

        run1 = frs(oktrials);
        run0 = frs(stilltrials);
        rmi(cell) = (runfr-norunfr)/(runfr+norunfr);

%         figure
%         subplot(2,2,1)
%         imagesc(speed);
%         colorbar
%         title(['oktrials: ' int2str(length(oktrials)) '/' int2str(size(speed,1))])
%         xlabel('time [ms]')
%         ylabel('trial number')
% 
%         subplot(2,2,2)
%         errorbar(msta,mean(speed),std(speed)./sqrt(size(speed,1)),'b')
%         xlabel('time [ms]')
%         ylabel('average runspeed')
% 
%         subplot(2,2,3)
%         plot(mean(speed(:,respwin),2),frs,'.')
%         xlabel('average runspeed of trial')
%         ylabel('average firing rate of trial')
% 
%         subplot(2,2,4)
%         barweb([runfr,norunfr],...
%             [runfrerr,norunfrerr],...
%             [],[{'running only'};{'immobile only'}],[],...
%             [],'firing rate [Hz]',[],[])
    end
    
    cell = cell+1;
    disp('');
    
end

swidth = swidth(:,1);

figure
plot(swidth,adiff,'k.')
xlabel('spike width')
ylabel('amplitude diff')

pfs = find(swidth<120);
prs = find(swidth>=120);
pfsv = swidth<120;
prsv = swidth>=120;

% adjust depth according to penetration angle
depth = depth.*cosd(22).*cosd(penangle);


for cl = 1:size(binavg,1)
    figure
    bamx = max(max(binavg(cl,:,:)));
    cfmx = max(max(condfr(cl,:,:)));
    cgmx = max(max(condgp(cl,:,:)));
    for blck = 1:length(blocks)
        subplot(3,length(blocks),blck)
        boundedline(ta,squeeze(binavg(cl,blck,:)),squeeze(binstd(cl,blck,:)),'k');
        axis([-prestim,trialdur+poststim,-0.05,bamx]);
        line([0,0],[0,bamx],'color','k','linewidth',2);
        line([2000,2000],[0,bamx],'color','k','linewidth',2);
        xlabel('time [ms]')
        ylabel('firingrate [Hz]')
        title(['depth: ' num2str(depth(cl,blck)) '   ' num2str(swidth(cl))])

        subplot(3,length(blocks),length(blocks)+blck)
        errorbar(oris,squeeze(condfr(cl,blck,:)),squeeze(conderr(cl,blck,:)),'ko-','markersize',8,'linewidth',2)
        xlabel('shown orientation')
        ylabel('Firing rate [Hz]')
        set(gca,'xtick',oris)
        title(['DSI: ' num2str(dosdsi(cl,blck))])
        axis([-20,335,-.05,cfmx])

        subplot(3,length(blocks),2*length(blocks)+blck)
        errorbar(oris,squeeze(condgp(cl,blck,:)),squeeze(condgperr(cl,blck,:)),'ko-','markersize',8,'linewidth',2)
        xlabel('shown orientation')
        ylabel('gamma power')
        set(gca,'xtick',oris)
        title(['DSI: ' num2str(dosdsi(cl,blck))])
        axis([-20,335,-.05,cgmx])
    end

end
  
% rs cells
figure
hold on
for i = 1:size(condfr,2)
    plot(i,mean(condfr(prsv,i,:),3),'bo');
end
for i = 1:length(prs)
    errorbar(1:size(condfr,2),mean(condfr(prs(i),:,:),3),mean(conderr(prs(i),:,:),3),'k')
end
title('firing rate RS cells over blocks')
xlabel('block number')
ylabel('firing rate [Hz]')

% fs cells
figure
hold on
for i = 1:size(condfr,2)
    plot(i,mean(condfr(pfsv,i,:),3),'bo');
end
for i = 1:length(pfs)
    errorbar(1:size(condfr,2),mean(condfr(pfs(i),:,:),3),mean(conderr(pfs(i),:,:),3),'k')
end
title('firing rate FS cells over blocks')
xlabel('block number')
ylabel('firing rate [Hz]')

% gp cells
% figure
% hold on
% for i = 1:size(condgp,2)
%     plot(i,mean(condgp(:,i,:),3),'bo');
% end
figure
for i = 1%1:size(condgp,1)
    errorbar(1:size(condgp,2),mean(condgp(i,:,:),3),mean(condgperr(i,:,:),3),'k')
end
title('gamma power over blocks')
xlabel('block number')
ylabel('power [a.u.]')

col = [0,1,2,1,0,0,0,0,0,0,0,0,0,1,1,2,2,1,2,1,2,1,2];
uvuv = zeros(1,length(col)-1);
uvcy = zeros(1,length(col)-1);
cyuv = zeros(1,length(col)-1);
cycy = zeros(1,length(col)-1);
for i = 1:22
    if col(i) == 1 & col(i+1) == 1;
        uvuv(i) = 1;
    elseif col(i) == 1 & col(i+1) == 2;
        uvcy(i) = 1;
    elseif col(i) == 2 & col(i+1) == 1;
        cyuv(i) = 1;
    elseif col(i) == 2 & col(i+1) == 2;
        cycy(i) = 1;
    end
end

figure
plot(dosdsi(:,1),dosdsi(:,3),'.');
line([0,.8],[0,.8])
