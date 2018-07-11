function MUA_population_randpulses_oneelectrode


% PV ChR2 population
animalids =    {'160423_2', '160423_2', '160423_2', '160423_2', '160423_2', '160425', '160425', '160425', '160425', '160425', '160503_3', '160503_3', '160503_3', '160503_3'};
blocks    =     [2,          3,          4,          5,          6,          2,        3,        4,        5,        6,        4,          5,          6,          7];
animal    =     [1,          1,          1,          1,          1,          2,        2,        2,        2,        2,        3,          3,          3,          3];
electrodes =   [[1,16];     [1,16];     [1,16];     [1,16];     [1,16];     [1,16];   [1,16];   [1,16];   [1,16];   [1,16];   [1,16];     [1,16];     [1,16];     [1,16]];
bldepth    =    [500,        500,        500,        500,        500,        500,      500,      500,      500,      500,      500,        500,        500,        500];
penangle =      [25,         25,         25,         25,         25,         25,       25,       25,       25,       25,       25,         25,         25,         25];
spacing =       [25,         25,         25,         25,         25,         25,       25,       25,       25,       25,       25,         25,         25,         25];
intensity =     [0.2,        0.3,        0.5,        1,          2,          0.2,      0.3,      0.5,      1,        2,        0.2,        0.3,        0.5,        1];
popfile    = 'C:\Users\Julia\work\data\populations\PV_ChR2\pulses\MUA_population_randpulses_oneelec.mat';
% goodinds = [1,7];
goodinds = 1:length(blocks);
% % 
% % % SOM ChR2 population 
% animalids =    {'160422_2', '160422_2', '160422_2', '160422_2', '160422_2', '160422_2', '160423', '160423', '160423', '160423', '160423', '160425_2', '160425_2', '160425_2', '160425_2', '160503', '160503', '160503', '160503', '160503_2', '160503_2', '160503_2'};
% blocks    =     [2,          7,          3,          4,          5,          6,          4,        5,        6,        7,        8,        3,          4,          5,          6,          2,        3,        4,        5,        3,          4,          5];
% animal    =     [1,          1,          1,          1,          1,          1,          2,        2,        2,        2,        2,        3,          3,          3,          3,          4,        4,        4,        4,        5,          5,          5];
% electrodes =   [[1,16];     [1,16];     [1,16];     [1,16];     [1,16];     [1,16];     [1,16];   [1,16];   [1,16];   [1,16];   [1,16];   [1,16];     [1,16];     [1,16];     [1,16];     [1,16];   [1,16];   [1,16];   [1,16];   [1,16];     [1,16];     [1,16]];
% bldepth    =    [500,        500,        500,        500,        500,        500,        600,      600,      600,      600,      600,      500,        500,        500,        500,        500,      500,      500,      500,      500,        500,        500];
% penangle =      [25,         25,         25,         25,         25,         25,         25,       25,       25,       25,       25,       25,         25,         25,         25,         25,       25,       25,       25,       25,         25,         25];
% spacing =       [25,         25,         25,         25,         25,         25,         25,       25,       25,       25,       25,       25,         25,         25,         25,         25,       25,       25,       25,       25,         25,         25];
% intensity =     [0.1,        0.15,       0.2,        0.3,        0.5,        1,          0.2,      0.3,      0.5,      1,        2,        0.2,        0.3,        0.5,        1,          0.2,      0.3,      0.5,      1,        0.2,        0.3,        0.5];
% popfile    = 'C:\Users\Julia\work\data\populations\PV_ChR2\pulses\MUA_population_randpulses_oneelec_new.mat';
% % goodinds = [3,8,12];
% goodinds = 1:length(blocks);

% % single thing
% animalids =    {'160503_2', '160503_2', '160503_2'};
% blocks    =     [3,        4,        5];
% animal    =     [1,        1,        1];
% electrodes =   [[1,16];   [1,16];   [1,16]];
% bldepth    =    [500,      500,      500];
% penangle =      [25,       25,       25];
% spacing =       [25,       25,       25];
% intensity =     [0.2,      0.3,      0.5];
% popfile = 'C:\Users\Julia\work\data\populations\temp2.mat';
% goodinds = [1,2,3];

depth = bldepth.*(cosd(penangle)*cosd(22));
spacing = spacing.*(cosd(penangle)*cosd(22));
evaldepth = 500;
for i = 1:length(depth)
    for j = 1:length(electrodes(i,1):electrodes(i,2))
        dm(i,j) = depth(i)-((j-1)*spacing(i)); % depth matrix
    end
    [c,di(i)] = min(abs(dm(i,:)-evaldepth)); % get depth index, least distance to evaldepth
end


recalculate_pop = 1;
recalculate_muafile = 0;

% chronux parameters
params.tapers = [2,5]; params.Fs = 1000; params.err = [2, 0.05]; params.trialave = 1;


if ~exist(popfile) || recalculate_pop

    cll = 1;
    for blck = goodinds %1:length(blocks)
        
        basepath = strcat('C:\Users\Julia\work\data\', animalids{blck}, '\');
        file = strcat(basepath, 'muaresult_', int2str(blocks(blck)), '_', int2str(electrodes(blck,1)), '_', int2str(electrodes(blck,2)), '.mat');
        
        clear result;
        if ~exist(file) || recalculate_muafile
            result = MUAdataprepare(basepath,animalids{blck},blocks(blck),electrodes(blck,1):electrodes(blck,2));
            save(file,'result')
        else
            clear centresult; clear surrresult;
            load(file);
            if exist('centresult')
                result = centresult;
            elseif exist('surrresult')
                result = surrresult;
            end            
        end        
        
        
        prestim = 0;
        poststim = 0;
        respwin = 1:result.stimduration*1000; % after stimulus onset
        respwin = respwin+prestim;
                
        ch = di(blck);
        
        lfp = result.lfp(ch,:);
        gamma = eegfilt(lfp,1000,15,40);  % adjust depending on peak
        h = hilbert(gamma); gpow = abs(h); gphas = angle(h);
            
        disp(['Block ' int2str(blck) '/' int2str(length(blocks)) '   channel ' int2str(ch) '/' int2str(length(electrodes(blck,1):electrodes(blck,2)))]);
        
        trialdur = result.stimduration*1000;
        msstamps = result.msstamps;
        if length(msstamps)~=size(result.light,1)
%             msstamps([193,291]) = []; % for 151023 block 5
%             msstamps([62]) = []; % for 16310 block 8
%             msstamps([21]) = []; % for 160404 block 5
%             msstamps([48]) = []; % for 160425 block 
%             result.msstamps = msstamps;
%             save(file,'result');
            pause;
        end
        
        for i = 1:length(msstamps)
            speed(i,:) = result.runspeed(msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim);
        end
        
        % figure out sufficiently high and nonvariable runspeed trials
        meanspeed = mean(speed(:,respwin),2);
        stdspeed = std(speed(:,respwin),1,2);
        notstill = find(meanspeed>1);
        okspeed = find(meanspeed>( mean(meanspeed(notstill))-(1.5*std(meanspeed(notstill))) ) & meanspeed > 1);
        okvar = find(stdspeed<( mean(stdspeed(notstill))+(1.5*std(stdspeed(notstill)))) & stdspeed>.5);
        oktrials = intersect(okspeed,okvar);
        nonoktrials = 1:size(speed,1); nonoktrials(oktrials) = [];
        stilltrials = 1:size(speed,1); stilltrials(notstill) = [];
        
        
        
        msStimes = round(result.msStimes{ch});
        if ~isempty(msStimes) & msStimes(1) == 0, msStimes(1) = 1; end
        
        chan_c = zeros(1,size(result.lfp,2));
        chan_c(msStimes) = 1;        
        
        light = logical(result.light);
        L = length(respwin); nfft = 2^nextpow2(L);
        t = (0:L-1)*(1/1000); fftx = 1000/2*linspace(0,1,nfft/2+1);
        allcyclelen = []; allblcyclelen = [];
        for i = 1:length(msstamps)
            resp(i,:) = chan_c(msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim);
            lfpresp(i,:) = result.lfp(ch,msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim);
            [lfpspect(i,:),fax] = pmtm(lfpresp(i,respwin),3,[],1000);
            hlp = fft(lfpresp(i,respwin),nfft)/L;
            fftspect(i,:) = 2*abs(hlp(1:nfft/2+1));
            bpresp(i,:) = gamma(msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim);
            gphaseresp(i,:) = gphas(msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim);
                        
            blphasediff{i} = find(diff(gphaseresp(i,1:300))<-5);
            allphasediff{i} = find(diff(gphaseresp(i,1:3500))<-5);
            blintervals{i} = diff(blphasediff{i});
            allintervals{i} = diff(allphasediff{i});
            
            meancycle = round(mean(allintervals{i}));
            allcyclelen = [allcyclelen, allintervals{i}];
            allblcyclelen = [allblcyclelen, blintervals{i}];
            
            stims = find(light(i,:));
            for s = 1:length(stims)
                stimphas(i,s) = gphaseresp(i,stims(s)); % find out phase of pulse
                % get all differences of one average phase interval starting 1/3 of an average cycle after stim pulse
                phasdiff = stimphas(i,s)-gphaseresp(i,stims(s)+round(meancycle/3):stims(s)+round(meancycle/3)+meancycle);
                [xx,ind] = min(abs(phasdiff)); % find minimum difference
                nextphas = stims(s)+round(meancycle/3)+ind-1;  % get index of next occurence of stim phase
                stimcyclelen(i,s) = nextphas-stims(s);
                thiscycle = bpresp(i,stims(s)+round(meancycle/3):stims(s)+round(meancycle/3)+meancycle);
%                 plot(thiscycle,[diff(thiscycle),0],'r')                
            end
            
        end
        
        r1allcyclelen = [];
        for i = 1:length(oktrials)
            r1allcyclelen = [r1allcyclelen, allintervals{oktrials(i)}];
        end
        r0allcyclelen = [];
        for i = 1:length(stilltrials)
            r0allcyclelen = [r0allcyclelen, allintervals{stilltrials(i)}];
        end
        
%         intervals = [-pi,-pi/2;-pi/2,0;0,pi/2;pi/2,pi];
        intervals = [-pi,-3*pi/4;-3*pi/4,-pi/2;-pi/2,-pi/4;-pi/4,0;0,pi/4;pi/4,pi/2;pi/2,3*pi/4;3*pi/4,pi];
        for i = 1:size(intervals,1)
            hlp = stimcyclelen(find(stimphas>intervals(i,1)&stimphas<intervals(i,2)));
            meanlen(i) = mean(hlp);
            medianlen(i) = median(hlp);
            meanlenerr(i) = std(hlp)./sqrt(length(hlp));
            
            clear hlp;
            r1stimcyclelen = stimcyclelen(oktrials,:);
            r1stimphas = stimphas(oktrials,:);
            hlp = r1stimcyclelen(find(r1stimphas>intervals(i,1)&r1stimphas<intervals(i,2)));
            r1meanlen(i) = mean(hlp);
            r1medianlen(i) = median(hlp);
            r1meanlenerr(i) = std(hlp)./sqrt(length(hlp));
            
            clear hlp;
            r0stimcyclelen = stimcyclelen(stilltrials,:);
            r0stimphas = stimphas(stilltrials,:);
            hlp = r0stimcyclelen(find(r0stimphas>intervals(i,1)&r0stimphas<intervals(i,2)));
            r0meanlen(i) = mean(hlp);
            r0medianlen(i) = median(hlp);
            r0meanlenerr(i) = std(hlp)./sqrt(length(hlp));            
        end
        
        x = 0:7;
%         figure
%         errorbar(-1,mean(meanblcyclelen),std(meanblcyclelen)./sqrt(length(meanblcyclelen)))
%         hold on
%         errorbar(0,mean(meancyclelen),std(meancyclelen)./sqrt(length(meancyclelen)))
%         errorbar(meanlen,meanlenerr)        
%         plot(x+.5,-cos(x*(pi/4)).*2+mean(meancyclelen),'r')
        
        figure
%         errorbar(-1,median(allblcyclelen),std(allblcyclelen)./sqrt(length(allblcyclelen)))
        hold on
%         errorbar(0,median(allcyclelen),std(allcyclelen)./sqrt(length(allcyclelen)))
        errorbar(medianlen,meanlenerr,'linewidth',2)
        errorbar(r1medianlen,r1meanlenerr,'g','linewidth',2);
        errorbar(r0medianlen,r1meanlenerr,'r','linewidth',2);
        plot([0,9],[median(allcyclelen),median(allcyclelen)],'b--')
%         plot([0,9],[mean(allcyclelen),mean(allcyclelen)],'b--')
        plot([0,9],[median(r1allcyclelen),median(r1allcyclelen)],'g--')
%         plot([0,9],[mean(r1allcyclelen),mean(r1allcyclelen)],'g--')
        plot([0,9],[median(r0allcyclelen),median(r0allcyclelen)],'r--')
%         plot([0,9],[mean(r0allcyclelen),mean(r0allcyclelen)],'r--')
%         plot(x+.5,-cos(x*(pi/4)).*2+median(allcyclelen),'r--')
        set(gca,'xticklabel',{'-pi','-3pi/4','-pi/2','-pi/4','0','pi/4','pi/2','3pi/4'})
        set(gca,'xtick',1:8)
        
        figure
        [a,i] = sort(stimphas(:));
        b = stimcyclelen(:);
        b = stimcyclelen(i);
        y = [b;b;b];
        x = [a;a+2*pi;a+4*pi];
        plot(x,smooth(y,150),'b')
        hold on
        plot([-4,16],[median(allcyclelen),median(allcyclelen)],'k')
        axis([0,2*pi+1,32,45]);
        line([pi,pi],[32,45],'color','k','linestyle','--')
        
%         j = 1;
%         a = find(result.light(j,:));
%         b = a+stimcyclelen(j,:);
%         plot(bpresp(j,:),'b.')
% %         plot(lfpresp(j,:),'b.')
%         hold on
%         plot(result.light(j,:).*50,'r')
%         plot(gphaseresp(j,:).*7,'k.')
%         for i = 1:20
%             line([b(i),b(i)],[0,50],'color','g')
%         end
        
        sclen(blck,:) = stimcyclelen(:);
        sphas(blck,:) = stimphas(:);
        mclen{blck} = allcyclelen;
        mblclen{blck} = allblcyclelen;
        
    end
    save(popfile, '-v7.3');
else
    load(popfile);
end

disp('');