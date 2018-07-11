function DS_ana_bar

animalid = '150120';
blocks = [3,6,9,12,17,20];
lcol = 'r';

onlymod = 0;
sfc = 0;

basepath = ['C:\Users\Julia\work\data\' animalid '\'];
supath = [basepath 'singleunits\'];
basename = [animalid '_block' int2str(blocks(1)) '_tet'];

files = dir([supath, basename, '*.mat']);

prestim = 300;
poststim = 300;
respwin = 1:6000; % after stimulus onset
respwin = respwin+prestim;
freqbinwidth = 5;

cl = 1;
for cell = 1:length(files)
    
    if strfind(files(cell).name, 'MU')
        continue;
    end
    
    figure
    
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


        trialdur = 6000; %result.stimduration*1000;
        msstamps = result.msstamps;
        if length(msstamps)~=length(result.light)
    %         msstamps(37) = []; % for 150113 block 7
%             msstamps(39) = []; % for 150113 block 11
%             msstamps(31) = []; % for 150119 block 5
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

        %determine if cell is visually modulated
        blfr = sum(resp(:,1:prestim),2);
        vrfr = sum(resp(:,prestim+40:2*prestim+40-1),2);
        vismod(cl,blck) = ttest2(blfr,vrfr);
        visdriven(cl,blck) = mean(vrfr)>=mean(blfr)+2; % average firing rate is increase at least 2Hz above baseline

        cellname{cl,blck} = files(cl).name;

        i = strfind(files(cl).name, 'tet');
        if strcmp(files(cl).name(i+4),'_')
            tetno = strread(files(cl).name(i+3)); % single character number
        else        
            tetno = strread(files(cl).name(i+3:i+4)); % number >10
        end
        tetnos(cl) = tetno;

        spike = result.waveforms(:,wvchan); 
        interpspike = spline(1:32,spike,1:.1:32);
        [adiff(cl),swidth(cl)] = spikequant(interpspike);

        % phases
        tmp = zeros(size(gphaseresp));
        tmp(find(resp)) = gphaseresp(find(resp));
        phases{cl,blck} = tmp(find(tmp));

        phaser(cl,blck) = circ_r(phases{cl});
        cmean(cl,blck) = circ_mean(phases{cl});

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
                allphaser(cl,blck,i) = circ_r(allphases{i}');
                allcmean(cl,blck,i) = circ_mean(allphases{i}');
            end
        end

        depth(cl) = result.depth;
        cellresp(cl,blck,:,:) = resp;
        celllfpresp(cl,blck,:,:) = lfpresp;
%         cellchan(cl,blck,:) = chan;

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

        baseline = mean(bl);
        baselineerr = std(bl)./(sqrt(size(bl,1)));

        binwidth = 50;
        oris = unique(result.gratingInfo.Orientation); oris(find(oris == -1)) = [];
        
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
            condfr(ori) = nanmean(frs(thisinds));
            conderr(ori) =nanstd(frs(thisinds))./sqrt(length(thisinds));
            condtrials(ori,:) = frs(thisinds);

        end

        bincondresp = bincondresp.*(1000/binwidth);
        bta = bta-prestim;

        [oriprefratio(cl,blck), dirprefratio(cl,blck), prefori(cl,blck), meanori(cl,blck), osi(cl,blck), meandir(cl,blck), dsi(cl,blck)] = getOSI(squeeze(condfr),oris);

        oneorifr = mean(reshape(condfr,4,2),2);
        prefori = find(oneorifr == max(oneorifr),1);
        [preffr(cl,blck), prefdir(cl,blck)] = max(condfr,[],2);
        ortho = mod(prefori+2,length(oris)/2); if ortho == 0, ortho = length(oris)/2; end
        nulldir(cl,blck) = prefdir(cl,blck)-length(oris)/2; 
        if nulldir(cl,blck)<1, nulldir(cl,blck) = nulldir(cl,blck)+length(oris); end
        dosdsi(cl,blck) = (condfr(prefdir(cl,blck))-condfr(nulldir(cl,blck)))./(condfr(prefdir(cl,blck))+condfr(nulldir(cl,blck)));

        [binpref,bta] = binit(squeeze(condresp(prefdir(cl,blck),:)),binwidth); binpref = binpref.*(1000/binwidth);

        [binavg,bta] = binit(mean(resp),binwidth); binavg = binavg.*(1000/binwidth); % in Hz
        binstd = binit(std(resp)./sqrt(size(resp,1)),binwidth);
        ta = bta-prestim;


%         figure
        subplot(3,length(blocks),blck)
        boundedline(ta,binavg,binstd,'k');
        mx = max(binavg);
        axis([-prestim,trialdur+poststim,-0.05,mx]);
        line([0,0],[0,mx],'color','k','linewidth',2);
        line([6000,6000],[0,mx],'color','k','linewidth',2);
        xlabel('time [ms]')
        ylabel('firingrate [Hz]')
        title(['depth: ' num2str(depth(cl)) '  swidth: ' num2str(swidth(cl))])

        subplot(3,length(blocks),length(blocks)+blck)
        errorbar(oris,squeeze(condfr),squeeze(conderr),'ko-','markersize',8,'linewidth',2)
        xlabel('shown orientation')
        ylabel('Firing rate [Hz]')
        set(gca,'xtick',oris)
        title(['DSI 1: ' num2str(dosdsi(cl,blck)) '  DSI 2: ' num2str(dsi(cl,blck))])
        ax = axis;
        axis([-20,335,ax(3),ax(4)])
        
        subplot(3,length(blocks),2*length(blocks)+blck)
        plot_tuning_polar(gca,condtrials,oris,'k')
        title(['DSI: ' num2str(dosdsi(cl,blck))])

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
    
    cl = cl+1;
    disp('');
    
end


% plot(dosdsi(:,1),dosdsi(:,3),'.');


figure
plot(swidth,adiff,'k.')
xlabel('spike width')
ylabel('amplitude diff')

figure
plot(dosdsi(:,1),dosdsi(:,3),'.');
line([0,.8],[0,.8])

figure
plot(1,dosdsi(:,1),'bo')
hold on
plot(2,dosdsi(:,2),'co')
plot(3,dosdsi(:,3),'ro')
plot(4,dosdsi(:,4),'mo')
plot(5,dosdsi(:,5),'bo')
for i = 1:6
    plot(1:5,dosdsi(i,:),'k')
end

figure
plot(1,dsi(:,1),'bo')
hold on
plot(2,dsi(:,2),'co')
plot(3,dsi(:,3),'ro')
plot(4,dsi(:,4),'bo')
for i = 1:16
    plot(1:4,dsi(i,:),'k')
end

disp('');