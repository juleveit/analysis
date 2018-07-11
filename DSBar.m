function DSBar

% enter animal id and all blocks from one condition (for example large
% gratings) and adjust the paths to the data

animalid = '150120';
blocks = [3,6,9,12,17,20];

basepath = ['C:\Users\Julia\work\data\' animalid '\'];
supath = [basepath 'singleunits\'];
basename = [animalid '_block' int2str(blocks(1)) '_tet'];

files = dir([supath, basename, '*.mat']);

prestim = 300;
poststim = 300;
trialdur = 6000;
respwin = 1:trialdur; % after stimulus onset
respwin = respwin+prestim;

cl = 1;
for cell = 1:length(files)
    
    if strfind(files(cell).name, 'MU')
        continue;
    end
    
    for blck = 1:length(blocks)
        
        basename = [animalid '_block' int2str(blocks(blck)) '_tet'];
        SUfiles = dir([supath, basename, '*.mat']);
    
        load([supath, SUfiles(cl).name]);
    
        msStimes = round(result.spikes);
        if ~isempty(msStimes) & msStimes(1) == 0, msStimes(1) = 1; end

        chan = zeros(1,length(result.lfp));
        chan(msStimes) = 1;

        wvchan = find(var(result.waveforms) == max(var(result.waveforms)));
        
        msstamps = result.msstamps;
        if length(msstamps)~=length(result.light) % unfortunately sometimes we have a timestamp or 
            % two too many (cheap trigger box)
            % I plot diff(msstimes) to find out which one and delete it
            %and save the file anew so it doesn't happen every time
            
%             msstamps(37) = []; % for 150113 block 7
%             msstamps(39) = []; % for 150113 block 11
%             msstamps(31) = []; % for 150119 block 5
%             result.msstamps = msstamps;
%             save([supath, SUfiles(cl).name],'result');
            pause;
        end
        
        for i = 1:length(msstamps)
            resp(i,:) = chan(msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim); 
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

        depth(cl) = result.depth;
        cellresp(cl,blck,:,:) = resp;
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
            condfr(cl,blck,ori) = nanmean(frs(thisinds));
            conderr(cl,blck,ori) =nanstd(frs(thisinds))./sqrt(length(thisinds));
            condtrials(cl,blck,ori,:) = frs(thisinds);

        end

        bincondresp = bincondresp.*(1000/binwidth);
        bta = bta-prestim;

        [oriprefratio(cl,blck), dirprefratio(cl,blck), prefori(cl,blck), meanori(cl,blck), osi(cl,blck), meandir(cl,blck), dsi(cl,blck)] = getOSI(squeeze(condfr(cl,blck,:))',oris);

        oneorifr = mean(reshape(condfr(cl,blck,:),4,2),2);
        prefori = find(oneorifr == max(oneorifr),1);
        [preffr(cl,blck), prefdir(cl,blck)] = max(condfr(cl,blck,:),[],3);
        ortho = mod(prefori+2,length(oris)/2); if ortho == 0, ortho = length(oris)/2; end
        nulldir(cl,blck) = prefdir(cl,blck)-length(oris)/2; 
        if nulldir(cl,blck)<1, nulldir(cl,blck) = nulldir(cl,blck)+length(oris); end
        dosdsi(cl,blck) = (condfr(cl,blck,prefdir(cl,blck))-condfr(cl,blck,nulldir(cl,blck)))./(condfr(cl,blck,prefdir(cl,blck))+condfr(cl,blck,nulldir(cl,blck)));

        [binpref,bta] = binit(squeeze(condresp(prefdir(cl,blck),:)),binwidth); binpref = binpref.*(1000/binwidth);

        [binavg(cl,blck,:),bta] = binit(mean(resp),binwidth); binavg(cl,blck,:) = binavg(cl,blck,:).*(1000/binwidth); % in Hz
        binstd(cl,blck,:) = binit(std(resp)./sqrt(size(resp,1)),binwidth);
        ta = bta-prestim;

    end
    
    cl = cl+1;
    
end


figure
plot(swidth,adiff,'k.')
xlabel('spike width')
ylabel('amplitude diff')

for cl = 1:size(binavg,1)
    figure
    bamx = max(max(binavg(cl,:,:)));
    cfmx = max(max(condfr(cl,:,:)));
    for blck = 1:length(blocks)
        subplot(3,length(blocks),blck)
        boundedline(ta,squeeze(binavg(cl,blck,:)),squeeze(binstd(cl,blck,:)),'k');
        axis([-prestim,trialdur+poststim,-0.05,bamx]);
        line([0,0],[0,bamx],'color','k','linewidth',2);
        line([6000,6000],[0,bamx],'color','k','linewidth',2);
        xlabel('time [ms]')
        ylabel('firingrate [Hz]')
        title(['depth: ' num2str(depth(cl))])

        subplot(3,length(blocks),length(blocks)+blck)
        errorbar(oris,squeeze(condfr(cl,blck,:)),squeeze(conderr(cl,blck,:)),'ko-','markersize',8,'linewidth',2)
        xlabel('shown orientation')
        ylabel('Firing rate [Hz]')
        set(gca,'xtick',oris)
        title(['DSI: ' num2str(dosdsi(cl,blck))])
        axis([-20,335,-.05,cfmx])

        subplot(3,length(blocks),2*length(blocks)+blck)
        plot_tuning_polar(gca,squeeze(condtrials(cl,blck,:,:)),oris,'k')
        title(['DSI: ' num2str(dosdsi(cl,blck))])
        axis([-cfmx,cfmx,-cfmx,cfmx])
    end

end
        

blblock = 1;
injblock = 3;
recovblock = 5;

figure
plot(dosdsi(:,blblock),dosdsi(:,injblock),'.');
line([0,1],[0,1])
title('OSI change through drug injection')
xlabel('DSI baseline')
ylabel('DSI injection')

figure
plot(1,dosdsi(:,blblock),'bo','markerfacecolor','b')
hold on
plot(2,dosdsi(:,injblock),'ro','markerfacecolor','r')
plot(3,dosdsi(:,recovblock),'co','markerfacecolor','c')
for i = 1:size(dosdsi,1)
    plot(1:3,dosdsi(i,[blblock,injblock,recovblock]),'k')
end
axis([0.8,3.2,0,1])
xlabel('experimental phase')
ylabel('DSI - difference over sum')
title('DSI difference over sum')


figure
plot(1,dsi(:,blblock),'bo','markerfacecolor','b')
hold on
plot(2,dsi(:,injblock),'ro','markerfacecolor','r')
plot(3,dsi(:,recovblock),'co','markerfacecolor','c')
for i = 1:size(dsi,1)
    plot(1:3,dsi(i,[blblock,injblock,recovblock]),'k')
end
axis([0.8,3.2,0,1])
xlabel('experimental phase')
ylabel('DSI - 1-circular variance')
title('DSI - 1-circular variance')
