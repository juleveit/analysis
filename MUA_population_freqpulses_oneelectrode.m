function MUA_population_freqpulses_oneelectrode

% % SOM ChR2 population
% animalids  = {'160310', '160310', '160311', '160318', '160318', '160318', '160319', '160319', '160319'};
% blocks     = [7,         8,        5,        9,        10,       11,       3,        4,        5];
% animal     = [1,         1,        1,        2,        2,        2,        3,        3,        3];
% electrodes = [[1,16];  [1,16];    [1,16];   [1,16];   [1,16];   [1,16];   [1,16];   [1,16];   [1,16]];
% bldepth    = [500,      500,       500,      400,      400,      400,      400,      400,      400];
% penangle   = [25,       25,        25,       25,       25,       25,       25,       25,       25];
% spacing    = [25,       25,        25,       25,       25,       25,       25,       25,       25];
% intensity  = [0.2,      0.5,       0.3,      0.2,      0.3,      0.5,      0.2,      0.3,      0.5];
% popfile    = 'C:\Users\Julia\work\data\populations\SOM_ChR2\pulses\MUA_population_freqpulses_oneelec.mat';

% % Scnn ChR2 population
% animalids  = {'160307', '160307', '160321', '160321', '160321'};
% blocks     = [ 5,        6,        3,        4,        5];
% animal     = [ 1,        1,        2,        2,        2];
% electrodes = [[1,16];   [1,16];   [1,16];   [1,16];   [1,16]];
% bldepth    = [500,       500,      550,      550,      550];
% penangle   = [25,        25,       25,       25,       25];
% spacing    = [25,        25,       25,       25,       25];
% intensity  = [0.5,       2,        0.2,      0.3,      0.5];
% popfile    = 'C:\Users\Julia\work\data\populations\scnn_chr2\pulses\MUA_population_freqpulses_oneelec.mat';

% % PV ChR2 population
% animalids =    {'160404', '160404', '160404', '160404', '160405', '160405', '160405', '160406', '160406', '160406', '160407', '160408', '160408', '160408'};
% blocks    =     [2,        3,        4,        5,        2,        3,        4,        2,        3,        4,        2,        2,        3,        4];
% animal    =     [1,        1,        1,        1,        2,        2,        2,        3,        3,        3,        4,        5,        5,        5];
% electrodes =   [[1,16];   [1,16];   [1,16];   [1,16];   [1,16];   [1,16];   [1,16];   [1,16];   [1,16];   [1,16];   [1,16];   [1,16];   [1,16];   [1,16]];
% bldepth    =    [500,      500,      500,      500,      500,      500,      500,      500,      500,      500,      500,      500,      500,      500];
% penangle =      [25,       25,       25,       25,       25,       25,       25,       25,       25,       25,       25,       25,       25,       25];
% spacing =       [25,       25,       25,       25,       25,       25,       25,       25,       25,       25,       25,       25,       25,       25];
% intensity =     [0.3,      0.5,      1,        3,        0.2,      0.3,      0.5,      0.2,      0.3,      0.5,      5,        1,        3,        5];
% popfile    = 'C:\Users\Julia\work\data\populations\PV_ChR2\pulses\MUA_population_freqpulses_oneelec.mat';

% % SOM ChR2 population only new electrode animals
% animalids =    {'160404_2', '160404_2', '160404_2', '160405_2', '160405_2', '160405_2', '160406_2', '160406_2', '160406_2', '160407_2', '160407_2', '160407_2', '160408_2', '160408_2', '160408_2'};
% blocks    =     [5,          6,          7,          6,          7,          8,          4,          5,          6,          5,          6,          7,          4,          5,          6];
% animal    =     [1,          1,          1,          2,          2,          2,          3,          3,          3,          4,          4,          4,          5,          5,          5];
% electrodes =   [[17,32];    [17,32];    [17,32];    [17,32];    [17,32];    [17,32];    [1,16];     [1,16];     [1,16];     [1,16];     [1,16];     [1,16];     [1,16];     [1,16];     [1,16]];
% bldepth    =    [560,        560,        560,        500,        500,        500,        500,        500,        500,        500,        500,        500,        500,        500,        500];
% penangle =      [25,         25,         25,         25,         25,         25,         25,         25,         25,         25,         25,         25,         25,         25,         25];
% spacing =       [25,         25,         25,         25,         25,         25,         25,         25,         25,         25,         25,         25,         25,         25,         25];
% intensity =     [1,          3,          5,          1,          3,          5,          .2,         .5,         1,          0.5,        1,          3,          1,          3,          5];
% popfile    = 'C:\Users\Julia\work\data\populations\PV_ChR2\pulses\MUA_population_freqpulses_oneelec_new.mat';

% single thing
animalids =    {'160907_2', '160907_2', '160907_2', '160907_2'};
blocks    =     [2,          3,          4,          5];
animal    =     [1,          1,          1,          1];
electrodes =   [[1,16];     [1,16];     [1,16];     [1,16]];
bldepth    =    [500,        500,        500,        500];
penangle =      [25,         25,         25,         25];
spacing =       [25,         25,         25,         25];
intensity =     [0.1,        0.2,        0.3,        0.5];
popfile = 'C:\Users\Julia\work\data\populations\temp2.mat';

depth = bldepth.*(cosd(penangle)*cosd(22));
spacing = spacing.*(cosd(penangle)*cosd(22));
evaldepth = 380;
for i = 1:length(depth)
    for j = 1:length(electrodes(i,1):electrodes(i,2))
        dm(i,j) = depth(i)-((j-1)*spacing(i)); % depth matrix
    end
    [c,di(i)] = min(abs(dm(i,:)-evaldepth)); % get depth index, least distance to evaldepth
end


recalculate_pop = 1;
recalculate_muafile = 0;

% chronux parameters
params.tapers = [5,9]; params.Fs = 1000; params.err = [2, 0.05]; params.trialave = 1;


if ~exist(popfile) || recalculate_pop

    cll = 1;
    for blck = 1:length(blocks)
        
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
        
        prestim = 300;
        poststim = 300;
        respwin = 501:3500; % after stimulus onset
        respwin = respwin+prestim;
                
        ch = di(blck);
            
        disp(['Block ' int2str(blck) '/' int2str(length(blocks)) '   channel ' int2str(ch) '/' int2str(length(electrodes(blck,1):electrodes(blck,2)))]);
        
        trialdur = result.stimduration*1000;
        msstamps = result.msstamps;
        if length(msstamps)~=length(result.light)
%             msstamps([193,291]) = []; % for 151023 block 5
%             msstamps([62]) = []; % for 16310 block 8
%             msstamps([21]) = []; % for 160404 block 5
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
        
        
        L = length(respwin); nfft = 2^nextpow2(L);
        t = (0:L-1)*(1/1000); fftx = 1000/2*linspace(0,1,nfft/2+1);
        for i = 1:length(msstamps)
            resp(i,:) = chan_c(msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim);
            lfpresp(i,:) = result.lfp(ch,msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim);
        end        
%         % remove trials with too much noise
%         allvar = std(lfpresp,[],2);
%         badones = find(allvar>median(allvar)+2*std(allvar));
%         lfpresp(badones,:) = nan(length(badones),size(lfpresp,2));
        % only then do the frequency transform
        for i = 1:length(msstamps)            
            [lfpspect(i,:),fax] = pmtm(lfpresp(i,respwin),3,[],1000);
            hlp = fft(lfpresp(i,respwin),nfft)/L;
            fftspect(i,:) = 2*abs(hlp(1:nfft/2+1));
        end
        
        frs = sum(resp(:,respwin),2)./(length(respwin)/1000);
        bl = sum(resp(:,1:prestim),2)./(prestim/1000);
        
        freqs = unique(result.frequencies);
        for f = 1:length(freqs)
            
            thisinds = find(result.frequencies == freqs(f));
            
            condn(blck,f) = length(thisinds);
            condresp(blck,f,:) = nanmean(resp(thisinds,:),1);
            condresperr(blck,f,:) = nanstd(resp(thisinds,:),1,1)./sqrt(length(thisinds));
            condlfpspect(blck,f,:) = nanmean(lfpspect(thisinds,:),1);
            condfftspect(blck,f,:) = nanmean(fftspect(thisinds,:),1);
            condlfpresp(blck,f,:) = nanmean(lfpresp(thisinds,:),1);
            cfr(blck,f) = nanmean(frs(thisinds));
            cerr(blck,f) =nanstd(frs(thisinds))./sqrt(length(thisinds));
            relspect(blck,f,:) = condlfpspect(blck,f,:)./condlfpspect(blck,1,:);
            relfftspect(blck,f,:) = condfftspect(blck,f,:)./condfftspect(blck,1,:);
            
            if freqs(f) == 0
                peakpow(blck,f) = NaN;
                meanpow(blck,f) = NaN;
                powerratio(blck,f) = NaN;
                otherratio(blck,f) = NaN;
                blratio(blck,f) = NaN;
                fftpowerratio(blck,f) = NaN;
                fftotherratio(blck,f) = NaN;
                fftblratio(blck,f) = NaN;
                relpowerratio(blck,f) = NaN;
                relpeakpow(blck,f) = NaN;
                relfftpowerratio(blck,f) = NaN;
                relfftpeakpow(blck,f) = NaN;
            else
                diffs = fax-freqs(f);
                f_ind = find(abs(diffs) == min(abs(diffs)));
                peakpow(blck,f) = mean(condlfpspect(blck,f,f_ind-2:f_ind+2),3);
                blpeakpow(blck,f) = mean(condlfpspect(blck,1,f_ind-2:f_ind+2),3);
                surroundpow(blck,f) = mean([mean(condlfpspect(blck,f,f_ind-7:f_ind-5),3),mean(condlfpspect(blck,f,f_ind+5:f_ind+7),3)]);
                meanpow(blck,f) = mean(condlfpspect(blck,f,1:2000),3);
                powerratio(blck,f) = peakpow(blck,f)./meanpow(blck,f);
                otherratio(blck,f) = peakpow(blck,f)./surroundpow(blck,f);
                blratio(blck,f) = peakpow(blck,f)./blpeakpow(blck,f);
                
                fftpeakpow(blck,f) = condfftspect(blck,f,f_ind);
                fftblpeakpow(blck,f) = condfftspect(blck,1,f_ind);
                fftsurroundpow(blck,f) = mean([mean(condfftspect(blck,f,f_ind-7:f_ind-5),3),mean(condfftspect(blck,f,f_ind+5:f_ind+7),3)]);
                fftmeanpow(blck,f) = mean(condfftspect(blck,f,1:2000),3);
                fftpowerratio(blck,f) = fftpeakpow(blck,f)./fftmeanpow(blck,f);
                fftotherratio(blck,f) = fftpeakpow(blck,f)./fftsurroundpow(blck,f);
                fftblratio(blck,f) = fftpeakpow(blck,f)./fftblpeakpow(blck,f);
                
                relpeakpow(blck,f) = relspect(blck,f,f_ind);
                relmeanpow(blck,f) = mean(relspect(blck,f,1:2000),3);
                relpowerratio(blck,f) = relpeakpow(blck,f)./relmeanpow(blck,f);
                
                relfftpeakpow(blck,f) = relfftspect(blck,f,f_ind);
                relfftmeanpow(blck,f) = mean(relfftspect(blck,f,1:2000));
                relfftpowerratio(blck,f) = relfftpeakpow(blck,f)./relfftmeanpow(blck,f);
            end
            
        end
    end
    save(popfile, '-v7.3');
else
    load(popfile);
end

depth = depth.*(cosd(penangle)*cosd(22));
spacing = spacing.*(cosd(penangle)*cosd(22));

evaldepth = 275;
for i = 1:length(depth)
    for j = 1:length(electrodes(i,1):electrodes(i,2))
        cdm(i,j) = depth(i)-((j-1)*spacing(i)); % center depth matrix
    end
    [c,di(i)] = min(abs(cdm(i,:)-evaldepth)); % get depth index, least distance to evaldepth
%     valid(i) = leffect(i,di(i));              % valid if light has a significant effect on large size stimuli
end
% valid(isnan(valid)) = 0;
valid = ones(size(di));
valid = logical(valid);
vi = find(valid);

pvinds = [1,2,3,4];
% pvinds = [2,6,9,13];
% pvinds = [4,7,10,14]; % strongest
for i = 1:length(pvinds)
    normcurve(i,:) = fftblratio(pvinds(i),1:10)./max(fftblratio(pvinds(i),1:10));
end
errorbar(freqs(1:10),mean(normcurve),std(normcurve)./sqrt(length(pvinds)),'g');

sominds = [1,4,7,10,13];
% sominds = [2,5,8,11,14];
% sominds = [3,6,8,12,15];  %strongest
for i = 1:length(sominds)
    normcurve(i,:) = fftblratio(sominds(i),1:10)./max(fftblratio(sominds(i),1:10));
end
errorbar(freqs(1:10),mean(normcurve),std(normcurve)./sqrt(length(sominds)),'r');

anms = unique(animal);
for i = 1:length(anms)
    thisinds = find(animal == anms(i));

    figure
    plot(freqs(1:10),relfftpeakpow(thisinds,1:10),'o-');
    title(int2str(anms(i)));
end

% ints = unique(intensity);
% for i = 1:length(ints)
%     thisinds = find(intensity == ints(i));
%     
% %     figure    
% %     subplot(2,2,1)
% %     plot(freqs(1:10),powerratio(thisinds,1:10),'o-')
% %     title('peak/all - mt spect')
% %     
% %     subplot(2,2,2)
% %     plot(freqs(1:10),blratio(thisinds,1:10),'o-')
% %     title('peak/baseline - mt spect')
% %     
% %     subplot(2,2,3)
% %     plot(freqs(1:10),fftpowerratio(thisinds,1:10),'o-')
% %     title('peak/all - fft spect')
% %     
% %     subplot(2,2,4)
% %     plot(freqs(1:10),fftblratio(thisinds,1:10),'o-')
% %     title('peak/baseline - fft spect')
%     
%     figure('Name',['intensity: ' num2str(ints(i)) 'V'])
%     subplot(2,2,1)
%     plot(freqs(1:10),relpowerratio(thisinds,1:10),'o-');
%     title('peak/all - mt div by bl')
%     
%     subplot(2,2,2)
%     plot(freqs(1:10),relpeakpow(thisinds,1:10),'o-');
%     title('norm peak - mt div by bl')
%     
%     subplot(2,2,3)
%     plot(freqs(1:10),relfftpowerratio(thisinds,1:10),'o-');
%     title('peak/all - fft div by bl')
%     
%     subplot(2,2,4)
%     plot(freqs(1:10),relfftpeakpow(thisinds,1:10),'o-');
%     title('norm peak - fft div by bl')   
% end


 disp('');