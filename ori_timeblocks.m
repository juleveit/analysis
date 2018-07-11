function ori_timeblocks

animalid = '150708';
blocks = [2 4 6];
penangle = 10;

% adjust the paths
basepath = ['C:\Users\Julia\work\data\' animalid '\'];
supath = [basepath 'singleunits\'];
basename = [animalid '_block' int2str(blocks(1)) '_tet'];

files = dir([supath, basename, '*.mat']);

cell = 1;
% loop all cells
for cl = 1:length(files)
    
    % this will skip all multi unit files - if you want them comment it out
    if strfind(files(cl).name, 'MU')
        continue;
    end
    
    % loop all blocks for this cell
    for blck = 1:length(blocks)
        
        basename = [animalid '_block' int2str(blocks(blck)) '_tet'];
        SUfiles = dir([supath, basename, '*.mat']);
    
        load([supath, SUfiles(cl).name]);    
        
        % adjust trial length
        trialdur = result.stimduration*1000;
        prestim = 200;
        respwin = 1:result.stimduration*1000; % after stimulus onset
        respwin = respwin+prestim; 
        fftlen = 500;

        % convert spiketimes to ms
        msStimes = round(result.spikes);
        if ~isempty(msStimes) & msStimes(1) == 0, msStimes(1) = 1; end

        % make spiketrain for entire recording
        chan = zeros(1,length(result.lfp));
        chan(msStimes) = 1;

        % find out which channel of the terode the cell is biggest
        wvchan = find(var(result.waveforms) == max(var(result.waveforms)));
        lfp = result.lfp(:,wvchan)';
        sr = 1000;
    
        % filter in the gamma range - adjust lower and upper bound if you
        % want to restrict the frequency bandwidth
        gamma = eegfilt(lfp,sr,30,90);
        h = hilbert(gamma); gpow = abs(h); gphas = angle(h);
        
        % make trials (0 andn 1 for light color) be one number per visual
        % stimulus presentation
        msstamps = result.msstamps;
        trials = repmat(result.trials',1,20)';
        trials = trials(:);
        
        % important: IF there is faulty timestamps this is going to pause.
        % comment out the pause and in the corrective lines and stop the
        % script after calculation of msstamps
        % Look at plot(diff(msstamps)) and find out which ones to take out.
        % Doublecheck that the result is correct because you are going to
        % overwrite the original files. This has to be done for every cell
        % because each has a different file with the wrong timestamps in it
        if length(msstamps)~=length(trials)
% %             disp('');
% %             msstamps(14) = []; % for 150417 block 2
%             if blck == 1
%                 msstamps([593,765]) = []; % for 150703 block 3
%             elseif blck == 3
%                 msstamps([125,369,640,723]) = []; % for 150703 block 7
%             elseif blck == 4
%                 msstamps(693) = [];
%             end
%             result.msstamps = msstamps;
%             save([supath, SUfiles(cl).name],'result');
            pause;
        end
        
        % loop trials and cut spiketrain, lfp and gammapower into trials of
        % the same length. Also calculate spectrum for each trial
        for i = 1:length(msstamps)
            resp(i,:) = chan(msstamps(i) - prestim+1:msstamps(i) + trialdur ); 
            lfpresp(i,:) = result.lfp(msstamps(i) - prestim+1:msstamps(i) + trialdur );
            gammaresp(i,:) = gpow(msstamps(i) - prestim+1:msstamps(i) + trialdur);
            gphaseresp(i,:) = gphas(msstamps(i) - prestim+1:msstamps(i) + trialdur );
            [lfpspect(i,:),fax] = pmtm(lfpresp(i,respwin(1:fftlen)),2,[],1000);
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

        % average firing, baseline firing, gamma power per trial
        frs = sum(resp(:,respwin),2)./(length(respwin)/1000); % convert to Hz
        bl = sum(resp(:,1:prestim),2)./(prestim/1000);
        gp = mean(gammaresp(:,respwin),2);
        blgp = mean(gammaresp(:,1:prestim),2);

        cellname{cell,blck} = files(cell).name;

        i = strfind(files(cell).name, 'tet');
        if strcmp(files(cell).name(i+4),'_')
            tetno = strread(files(cell).name(i+3)); % single character number
        else        
            tetno = strread(files(cell).name(i+3:i+4)); % number >10
        end
        tetnos(cell) = tetno;

        % save spike waveform and quantify for RS FS diff
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

        msta = linspace(-prestim,trialdur,size(resp,2));

        baseline(cell,blck) = mean(bl);
        baselineerr = std(bl)./(sqrt(size(bl,1)));

        binwidth = 50;
        
        % sort out trials by trial type, average all same conditions
        lighttypes = unique(trials);
        contrasts = unique(result.gratingInfo.Contrast);
        oris = unique(result.gratingInfo.Orientation);
        % loop light type, contrasts and orientations
        for tr = 1:length(lighttypes)
            for co = 1:length(contrasts)
                for ori = 1:length(oris)
                    thisinds = find(trials' == tr-1 & result.gratingInfo.Contrast == contrasts(co) & result.gratingInfo.Orientation == oris(ori));
                    condn(tr,co,ori) = length(thisinds);
                    condresp(tr,co,ori,:) = nanmean(resp(thisinds,:),1);
                    condresperr(tr,co,ori,:) = nanstd(resp(thisinds,:),1,1)./sqrt(length(thisinds));
                    if ~isnan(condresp(tr,co,ori,:))
                        [bincondresp(tr,co,ori,:),bta] = binit(condresp(tr,co,ori,:),binwidth);
                    else
                        bincondresp(tr,co,ori,:) = binit(condresp(tr,co,ori,:),binwidth);
                    end
                    binconderr(tr,co,ori,:) = binit(condresperr(tr,co,ori,:),binwidth);
                    condfr(cell,blck,tr,co,ori) = nanmean(frs(thisinds));
                    conderr(cell,blck,tr,co,ori) =nanstd(frs(thisinds))./sqrt(length(thisinds));
                    condgp(cell,blck,tr,co,ori) = nanmean(gp(thisinds));
                    condgperr(cell,blck,tr,co,ori) = nanstd(gp(thisinds))./sqrt(length(thisinds));
                    condlfpspect(cell,blck,tr,co,ori,:) = nanmean(lfpspect(thisinds,:));
                    condlfpresp(cell,blck,tr,co,ori,:) = nanmean(lfpresp(thisinds,:));
                end
            end
        end

        bincondresp = bincondresp.*(1000/binwidth);
        bta = bta-prestim;
       
        [binavg(cell,blck,:,1),bta] = binit(mean(resp(trials == 0,:)),binwidth); binavg(cell,blck,:,1) = binavg(cell,blck,:,1).*(1000/binwidth); % in Hz
        [binavg(cell,blck,:,2),bta] = binit(mean(resp(trials == 1,:)),binwidth); binavg(cell,blck,:,2) = binavg(cell,blck,:,2).*(1000/binwidth); % in Hz
        binstd(cell,blck,:,1) = binit(std(resp(trials == 0,:))./sqrt(length(find(trials == 0))),binwidth);
        binstd(cell,blck,:,2) = binit(std(resp(trials == 1,:))./sqrt(length(find(trials == 1))),binwidth);
        ta = bta-prestim;
        
        % calculate orientation tuning for both light types
        for tr = 1:length(lighttypes)
            [nlprefratio, nldirpreffratio, nlprefori, nlmeanori, nlosi(cell,blck,tr), nlmeandir, nldsi(cell,blck,tr)] = getOSI(squeeze(condfr(cell,blck,tr,5,:))',oris);
        end

    end
    
    cell = cell+1;
    
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
    cfmx = max(max(max(condfr(cl,:,:,5,:))));
    ccmx = max(max(max(nanmean(condfr(cl,:,:,:,:),5))));
    cgmx = max(max(max(nanmean(condgp(cl,:,:,:,:),5))));
    cgmn = min(min(min(nanmean(condgp(cl,:,:,:,:),5))));
    for blck = 1:length(blocks)
        subplot(3,length(blocks),blck)
        boundedline(ta,squeeze(binavg(cl,blck,:,1)),squeeze(binstd(cl,blck,:,1)),'m');
        hold on
        boundedline(ta,squeeze(binavg(cl,blck,:,2)),squeeze(binstd(cl,blck,:,2)),'c');
        axis([-prestim,trialdur,-0.05,bamx]);
        line([0,0],[0,bamx],'color','k');%,'linewidth',2);
        xlabel('time [ms]')
        ylabel('firingrate [Hz]')
        title(['depth: ' num2str(depth(cl,blck)) '   ' num2str(swidth(cl))])

        % orientation
        subplot(3,length(blocks),length(blocks)+blck)
        errorbar(oris,squeeze(condfr(cl,blck,1,5,:)),squeeze(conderr(cl,blck,1,5,:)),'mo-','markersize',8,'linewidth',2)
        hold on
        errorbar(oris,squeeze(condfr(cl,blck,2,5,:)),squeeze(conderr(cl,blck,2,5,:)),'co-','markersize',8,'linewidth',2)
        xlabel('light color')
        ylabel('Firing rate [Hz]')
        axis([-15,375,-.05,cfmx])        
         
%         % contrast
%         subplot(3,length(blocks),length(blocks)+blck)
%         errorbar(contrasts,squeeze(nanmean(condfr(cl,blck,1,:,:),5)),squeeze(nanmean(conderr(cl,blck,1,:,:),5)),'mo-','markersize',8,'linewidth',2)
%         hold on        
%         errorbar(contrasts,squeeze(nanmean(condfr(cl,blck,2,:,:),5)),squeeze(nanmean(conderr(cl,blck,2,:,:),5)),'co-','markersize',8,'linewidth',2)
%         xlabel('light color')
%         ylabel('Firing rate [Hz]')
%         axis([-.2,1.2,-.05,ccmx])

        subplot(3,length(blocks),2*length(blocks)+blck)
        errorbar(contrasts,squeeze(nanmean(condgp(cl,blck,1,:,:),5)),squeeze(nanmean(condgperr(cl,blck,1,:,:),5)),'mo-','markersize',8,'linewidth',2)
        hold on
        errorbar(contrasts,squeeze(nanmean(condgp(cl,blck,2,:,:),5)),squeeze(nanmean(condgperr(cl,blck,2,:,:),5)),'co-','markersize',8,'linewidth',2)
        xlabel('light color')
        ylabel('gamma power')
        axis([-.2,1.2,cgmn,cgmx])
        title(SUfiles(cl).name)
    end
end

block = 1;
figure
plot(squeeze(nanmean(condfr(:,block,1,5,:),5)),squeeze(nanmean(condfr(:,block,2,5,:),5)),'.')

figure
plot(squeeze(nanmean(condgp(:,block,1,5,:),5)),squeeze(nanmean(condgp(:,block,2,5,:),5)),'.')

  
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

