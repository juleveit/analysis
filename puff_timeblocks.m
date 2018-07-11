function grating_timeblocks

animalid = '150617';
blocks = [1:4];
penangle = 10;

lcol = 'c';

basepath = ['C:\Users\Julia\work\data\' animalid '\'];
supath = [basepath 'singleunits\'];
basename = [animalid '_block' int2str(blocks(1)) '_tet'];

files = dir([supath, basename, '*.mat']);

cell = 1;
for cl = 1:length(files)
    
%     if strfind(files(cl).name, 'MU')
%         continue;
%     end
    
    for blck = 1:length(blocks)
        
        basename = [animalid '_block' int2str(blocks(blck)) '_tet'];
        SUfiles = dir([supath, basename, '*.mat']);
    
        load([supath, SUfiles(cl).name]);        
        
        trialdur = (result.sweeplength-result.initiallag)*1000;%+1000; %+1000 only for 150413
        prestim = result.initiallag*1000;
        respwin = 1:2000; % after stimulus onset
        respwin = respwin+prestim;        
        
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

        if ~isfield(result,'sweeplength')
            result.sweeplength = 4;
            save([supath, SUfiles(cl).name],'result');
        end
        
        
        if length(result.msstamps) ~= length(result.trials)*result.npuffsPerLight
            a = find(diff(result.msstamps)>3000);
            newstamps = [result.msstamps(1),result.msstamps(a+1)'];
            msstamps = newstamps;
        end
        
        % ming chis 150617 experiment
        a = find(diff(result.msstamps)>3500);
        msstamps = [result.msstamps(1);result.msstamps(a+1)];
%         begstamps = result.msstamps(mod(1:length(result.msstamps),result.npuffsPerLight) == 1);
%         msstamps = begstamps;
%         msstamps = result.msstamps;
%         trials = repmat(result.trials',1,result.npuffsPerLight)';
%         trials = trials(:);
        trialdur = 2000;%/result.pufffreq;
        prestim = 3000;
        if length(msstamps)~=length(result.trials)
% %             disp('');
%     %         msstamps([54,201,239,316]) = []; % for 140807 block 7
%     %         msstamps(200) = []; % for 140815 block 7
%     %         msstamps(123) = []; % for 140815 block 7
%     %         msstamps(385) = []; % for 140703 block 5
%     %         msstamps(37) = []; % for 150113 block 7
% %             msstamps(39) = []; % for 150113 block 11
% %             msstamps(18) = []; % for 150414 block 7
% %             msstamps(16) = []; % for 150414 block 10
% %             msstamps(71) = []; % for 150417 block 2
% %             msstamps(14) = []; % for 150417 block 2
% %             result.msstamps = msstamps;
% %             save([supath, SUfiles(cl).name],'result');
            pause;
        end
        
    %     trialnfft = 2^nextpow2(800);
    %     trialfax = sr/2*linspace(0,1,trialnfft/2+1);
        for i = 1:length(msstamps)
            resp(i,:) = chan(msstamps(i) - prestim+1:msstamps(i) + trialdur ); 
            lfpresp(i,:) = result.lfp(msstamps(i) - prestim+1:msstamps(i) + trialdur );
    %         y = fft(lfpresp(i,1001:1800),trialnfft);
    %         lfpspect(i,:) = abs(y(1:trialnfft/2+1));
%             [lfpspect(i,:),trialfax] = pmtm(lfpresp(i,1001:1800),3,[],sr);
            gammaresp(i,:) = gpow(msstamps(i) - prestim+1:msstamps(i) + trialdur);
            gphaseresp(i,:) = gphas(msstamps(i) - prestim+1:msstamps(i) + trialdur );
            [lfpspect(i,:),fax] = pmtm(lfpresp(i,respwin),2,[],1000);
            [lfpspectbefore(i,:),faxb] = pmtm(lfpresp(i,2001:3000),2,[],1000);

            speed(i,:) = result.runspeed(msstamps(i) - prestim+1:msstamps(i) + trialdur);     
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
%         vrfr = sum(resp(:,prestim+40:2*prestim+40-1),2);
%         vismod(cell,blck) = ttest2(blfr,vrfr);
%         visdriven(cell,blck) = mean(vrfr)>=mean(blfr)+2; % average firing rate is increase at least 2Hz above baseline

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

        msta = linspace(-prestim,trialdur,size(resp,2));

        baseline(cell,blck) = mean(bl);
        baselineerr = std(bl)./(sqrt(size(bl,1)));

        binwidth = 50;
        
        lighttypes = unique(result.trials);
        for tr = 1:length(lighttypes)
            thisinds = find(result.trials == tr-1);
            condn(tr) = length(thisinds);
            condresp(tr,:) = nanmean(resp(thisinds,:),1);
            condresperr(tr,:) = nanstd(resp(thisinds,:),1,1)./sqrt(length(thisinds));
            if ~isnan(condresp(tr,:))
                [bincondresp(tr,:),bta] = binit(condresp(tr,:),binwidth);
            else
                bincondresp(tr,:) = binit(condresp(tr,:),binwidth);
            end
            binconderr(tr,:) = binit(condresperr(tr,:),binwidth);
            condfr(cell,blck,tr) = nanmean(frs(thisinds));
            conderr(cell,blck,tr) =nanstd(frs(thisinds))./sqrt(length(thisinds));
            condgp(cell,blck,tr) = nanmean(gp(thisinds));
            condgperr(cell,blck,tr) = nanstd(gp(thisinds))./sqrt(length(thisinds));
            condtrials(cell,blck,tr,:) = frs(thisinds);
        end

        bincondresp = bincondresp.*(1000/binwidth);
        bta = bta-prestim;
       
        [binavg(cell,blck,:,1),bta] = binit(mean(resp(result.trials == 0,:)),binwidth); binavg(cell,blck,:,1) = binavg(cell,blck,:,1).*(1000/binwidth); % in Hz
        [binavg(cell,blck,:,2),bta] = binit(mean(resp(result.trials == 1,:)),binwidth); binavg(cell,blck,:,2) = binavg(cell,blck,:,2).*(1000/binwidth); % in Hz
        binstd(cell,blck,:,1) = binit(std(resp(result.trials == 0,:))./sqrt(length(find(result.trials == 0))),binwidth);
        binstd(cell,blck,:,2) = binit(std(resp(result.trials == 1,:))./sqrt(length(find(result.trials == 1))),binwidth);
        ta = bta-prestim;

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
        boundedline(ta,squeeze(binavg(cl,blck,:,1)),squeeze(binstd(cl,blck,:,1)),'m');
        hold on
        boundedline(ta,squeeze(binavg(cl,blck,:,2)),squeeze(binstd(cl,blck,:,2)),'c');
        axis([-prestim,trialdur,-0.05,bamx]);
        line([0,0],[0,bamx],'color','k','linewidth',2);
        line([2000,2000],[0,bamx],'color','k','linewidth',2);
        xlabel('time [ms]')
        ylabel('firingrate [Hz]')
        title(['depth: ' num2str(depth(cl,blck)) '   ' num2str(swidth(cl))])

        subplot(3,length(blocks),length(blocks)+blck)
        errorbar(squeeze(condfr(cl,blck,:)),squeeze(conderr(cl,blck,:)),'ko-','markersize',8,'linewidth',2)
        xlabel('light color')
        ylabel('Firing rate [Hz]')
        axis([0,3,-.05,cfmx])

        subplot(3,length(blocks),2*length(blocks)+blck)
        errorbar(squeeze(condgp(cl,blck,:)),squeeze(condgperr(cl,blck,:)),'ko-','markersize',8,'linewidth',2)
        xlabel('light color')
        ylabel('gamma power')
        axis([0,3,-.05,cgmx])
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
