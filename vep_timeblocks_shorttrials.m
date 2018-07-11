function vep_timeblocks_shorttrials

animalid = '150703';
blocks = [4,6,8,10];
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

        % adjust trial length - short trial version
        trialdur = result.ifi*1000;
        prestim = 200;
        respwin = 1:result.ifi*1000; % after stimulus onset
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
        sr = 1000;
        lfp = result.lfp(:,wvchan)';        
        
        % filter in the gamma range - adjust lower and upper bound if you
        % want to restrict the frequency bandwidth
        gamma = eegfilt(lfp,sr,20,60);
        h = hilbert(gamma); gpow = abs(h); gphas = angle(h);
        % %

      
        % important: IF there is faulty timestamps this is going to pause.
        % comment out the pause and in the corrective lines and stop the
        % script after calculation of msstamps
        % Look at plot(diff(msstamps)) and find out which ones to take out.
        % Doublecheck that the result is correct because you are going to
        % overwrite the original files. This has to be done for every cell
        % because each has a different file with the wrong timestamps in it        
        if length(result.msstamps) ~= length(result.trials)*result.nflashesPerLight
%             a = find(diff(result.msstamps)>3000);
%             newstamps = [result.msstamps(1),result.msstamps(a+1)'];
%             msstamps = newstamps;
%             result.msstamps([144,296]) = [];    % 150701 
%             if blck == 2
%                 result.msstamps([130,282]) = [];    % 150703  block 6    
%             elseif blck == 4
%                 result.msstamps(291) = []; %150703 block 10
%             end
%             save([supath, SUfiles(cl).name],'result');
            pause
        end
        
        % long trials: make trials (light) to be for each stimulus presentation but keep timestamp for
        % beginning of a repetition
        msstamps = result.msstamps;
        trials = repmat(result.trials',1,10)';
        trials = trials(:);
        
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
        frs = sum(resp(:,respwin),2)./(length(respwin)/1000);
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
        
        lighttypes = unique(trials);
        % sort out trials by trial type, average all same conditions
        % loop light type
        for tr = 1:length(lighttypes)
            thisinds = find(trials == tr-1);
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
            condlfpspect(cell,blck,tr,:) = nanmean(lfpspect(thisinds,:));
            condlfpspecterr(cell,blck,tr,:) = nanstd(lfpspect(thisinds,:))./sqrt(length(thisinds));            
            condlfpresp(cell,blck,tr,:) = nanmean(lfpresp(thisinds,:));
%             condtrials(cell,blck,tr,:) = frs(thisinds);
        end

        bincondresp = bincondresp.*(1000/binwidth);
        bta = bta-prestim;
       
        [binavg(cell,blck,:,1),bta] = binit(mean(resp(trials == 0,:)),binwidth); binavg(cell,blck,:,1) = binavg(cell,blck,:,1).*(1000/binwidth); % in Hz
        [binavg(cell,blck,:,2),bta] = binit(mean(resp(trials == 1,:)),binwidth); binavg(cell,blck,:,2) = binavg(cell,blck,:,2).*(1000/binwidth); % in Hz
        binstd(cell,blck,:,1) = binit(std(resp(trials == 0,:))./sqrt(length(find(trials == 0))),binwidth);
        binstd(cell,blck,:,2) = binit(std(resp(trials == 1,:))./sqrt(length(find(trials == 1))),binwidth);
        ta = bta-prestim;
        
        for i = 1:result.repetitions*2
            mfr(cell,blck,i) = mean(frs((i-1)*result.nflashesPerLight+1:i*result.nflashesPerLight));
            mfrerr(cell,blck,i) = std(frs((i-1)*result.nflashesPerLight+1:i*result.nflashesPerLight))./sqrt(result.nflashesPerLight);
            mgp(cell,blck,i) = mean(gp((i-1)*result.nflashesPerLight+1:i*result.nflashesPerLight));
            mgperr(cell,blck,i) = std(gp((i-1)*result.nflashesPerLight+1:i*result.nflashesPerLight))./sqrt(result.nflashesPerLight);
        end

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

%zigzag
for cl = 1:size(binavg,1)
    figure
    for blck = 1:length(blocks)
        subplot(2,length(blocks),blck)
        errorbar(mfr(cl,blck,:),mfrerr(cl,blck,:));
        xlabel('repetition')
        ylabel('firingrate [Hz]')
        title(['depth: ' num2str(depth(cl,blck)) '   ' num2str(swidth(cl))])

        subplot(2,length(blocks),length(blocks)+blck)
        errorbar(mgp(cl,blck,:),mgperr(cl,blck,:));
        xlabel('repetition')
        ylabel('gamma power')

    end

end

%firing rates and gamma
for cl = 1:size(binavg,1)
    figure
    bamx = max(max(binavg(cl,:,:)));
    cfmx = max(max(condfr(cl,:,:)));
    cgmx = max(max(condgp(cl,:,:)));
    cgmn = min(min(condgp(cl,:,:)));
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

        subplot(3,length(blocks),length(blocks)+blck)
        errorbar(squeeze(condfr(cl,blck,:)),squeeze(conderr(cl,blck,:)),'ko-','markersize',8,'linewidth',2)
        xlabel('light color')
        ylabel('Firing rate [Hz]')
        axis([0,3,-.05,cfmx])

        subplot(3,length(blocks),2*length(blocks)+blck)
        errorbar(squeeze(condgp(cl,blck,:)),squeeze(condgperr(cl,blck,:)),'ko-','markersize',8,'linewidth',2)
        xlabel('light color')
        ylabel('gamma power')
        axis([0,3,cgmn,cgmx])
        title(SUfiles(cl).name)
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
