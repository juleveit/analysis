function MUA_population_freqpulses

% % SOM ChR2 population
% animalids  = {'160310', '160310', '160311', '160318', '160318', '160318', '160319', '160319', '160319'};
% blocks     = [7,         8,        5,        9,        10,       11,       3,        4,        5];
% animal     = [1,         1,        1,        2,        2,        2,        3,        3,        3];
% electrodes = [[1,16];  [1,16];    [1,16];   [1,16];   [1,16];   [1,16];   [1,16];   [1,16];   [1,16]];
% bldepth    = [500,      500,       500,      400,      400,      400,      400,      400,      400];
% penangle   = [25,       25,        25,       25,       25,       25,       25,       25,       25];
% spacing    = [25,       25,        25,       25,       25,       25,       25,       25,       25];
% intensity  = [0.2,      0.5,       0.3,      0.2,      0.3,      0.5,      0.2,      0.3,      0.5];
% popfile    = 'C:\Users\Julia\work\data\populations\SOM_ChR2\pulses\MUA_population_freqpulses.mat';

% % Scnn ChR2 population
% animalids  = {'160307', '160307', '160321', '160321', '160321'};
% blocks     = [ 5,        6,        3,        4,        5];
% animal     = [ 1,        1,        2,        2,        2];
% electrodes = [[1,16];   [1,16];   [1,16];   [1,16];   [1,16]];
% bldepth    = [500,       500,      550,      550,      550];
% penangle   = [25,        25,       25,       25,       25];
% spacing    = [25,        25,       25,       25,       25];
% intensity  = [0.5,       2,        0.2,      0.3,      0.5];
% popfile    = 'C:\Users\Julia\work\data\populations\scnn_chr2\pulses\MUA_population_freqpulses.mat';

% single thing
animalids =    {'160907_2'};
blocks    =     [2];
animal    =     [1];
electrodes =[[1,16]];
bldepth    =  [500];
penangle =      [25];
spacing =       [25];
popfile = 'C:\Users\Julia\work\data\populations\temp.mat';


recalculate_pop = 1;
recalculate_muafile = 0;

% chronux parameters
params.tapers = [2,5]; params.Fs = 1000; params.err = [2, 0.05]; params.trialave = 1;


if ~exist(popfile) || recalculate_pop

    cll = 1;
    for blck = 1:length(blocks)
        
        basepath = strcat('C:\Users\Julia\work\data\', animalids{blck}, '\');
        file = strcat(basepath, 'muaresult_', int2str(blocks(blck)), '_', int2str(electrodes(blck,1)), '_', int2str(electrodes(blck,2)), '.mat');
        oldfile = strcat(basepath, 'muaresult_', int2str(blocks(blck)), '_', int2str(electrodes(blck,1)), ':', int2str(electrodes(blck,2)), '.mat');
        
        if ~exist(file) || recalculate_muafile
            if ~exist(oldfile)
                result = MUAdataprepare(basepath,animalids{blck},blocks(blck),electrodes(blck,1):electrodes(blck,2));
                save(file,'result')
            else
                load(oldfile);
                save(file, 'result');
            end
        else
            load(file);
        end       
        
        prestim = 300;
        poststim = 300;
        respwin = 501:3500; % after stimulus onset
        respwin = respwin+prestim;
        
        for ch = 1: length(electrodes(blck,1):electrodes(blck,2))
            
            disp(['Block ' int2str(blck) '/' int2str(length(blocks)) '   channel ' int2str(ch) '/' int2str(length(electrodes(blck,1):electrodes(blck,2)))]);
            
            depth(blck,ch) = bldepth(blck)-(ch-1)*spacing(blck);
            
            trialdur = result.stimduration*1000;
            msstamps = result.msstamps;
            if length(msstamps)~=length(result.light)
%                 msstamps([193,291]) = []; % for 151023 block 5
%                 msstamps([260]) = []; % for 151023 block 5
%                 result.msstamps = msstamps;
%                 surrresult.msstamps = msstamps;
%                 save(surrfile,'surrresult');
%                 save(centfile,'result');
                pause;
            end
            
            for i = 1:length(msstamps)
                speed(i,:) = result.runspeed(msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim);
            end
       
            % figure out sufficiently high and nonvariable runspeed trials
            meanspeed = mean(speed(:,respwin),2);
            stdspeed = std(speed(:,respwin),1,2);
            notstill = find(meanspeed>1);
            okspeed = find(meanspeed>( mean(meanspeed(notstill))-(1.5*std(meanspeed(notstill))) )& meanspeed>1 );
            okvar = find(stdspeed<( mean(stdspeed(notstill))+(1.5*std(stdspeed(notstill)))) & stdspeed>.5);
            oktrials = intersect(okspeed,okvar);
            nonoktrials = 1:size(speed,1); nonoktrials(oktrials) = [];
            stilltrials = 1:size(speed,1); stilltrials(notstill) = [];
            
            msStimes = round(result.msStimes{ch});
            if ~isempty(msStimes) & msStimes(1) == 0, msStimes(1) = 1; end

            chan = zeros(1,size(result.lfp,2));
            chan(msStimes) = 1;
    
            L = length(respwin); nfft = 2^nextpow2(L);
            t = (0:L-1)*(1/1000); fftx = 1000/2*linspace(0,1,nfft/2+1);
            
            for i = 1:length(msstamps)
                resp(i,:) = chan(msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim);
                lfpresp(i,:) = result.lfp(ch,msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim);
                [lfpspect(i,:),fax] = pmtm(lfpresp(i,respwin),3,[],1000);
                hlp = fft(lfpresp(i,respwin),nfft)/L;
                fftspect(i,:) = 2*abs(hlp(1:nfft/2+1));
            end
            
            
            frs = sum(resp(:,respwin),2)./(length(respwin)/1000);
            bl = sum(resp(:,1:prestim),2)./(prestim/1000);
            
            freqs = unique(result.frequencies);
            for f = 1:length(freqs)
                
                thisinds = find(result.frequencies == freqs(f));
                
                condn(blck,ch,f) = length(thisinds);
                condresp(blck,ch,f,:) = nanmean(resp(thisinds,:),1);
                condresperr(blck,ch,f,:) = nanstd(resp(thisinds,:),1,1)./sqrt(length(thisinds));
                condlfpspect(blck,ch,f,:) = nanmean(lfpspect(thisinds,:),1);
                condfftspect(blck,ch,f,:) = nanmean(fftspect(thisinds,:),1);
                condlfpresp(blck,ch,f,:) = nanmean(lfpresp(thisinds,:),1);
                cfr(blck,ch,f) = nanmean(frs(thisinds));
                cerr(blck,ch,f) =nanstd(frs(thisinds))./sqrt(length(thisinds));
                relspect(blck,ch,f,:) = condlfpspect(blck,ch,f,:)./condlfpspect(blck,ch,1,:);
                relfftspect(blck,ch,f,:) = condfftspect(blck,ch,f,:)./condfftspect(blck,ch,1,:);
                
                if freqs(f) == 0
                    peakpow(blck,ch,f) = NaN;
                    meanpow(blck,ch,f) = NaN;
                    powerratio(blck,ch,f) = NaN;
                    otherratio(blck,ch,f) = NaN;
                    blratio(blck,ch,f) = NaN;
                    fftpowerratio(blck,ch,f) = NaN;
                    fftotherratio(blck,ch,f) = NaN;
                    fftblratio(blck,ch,f) = NaN;
                    relpowerratio(blck,ch,f) = NaN;
                    relpeakpow(blck,ch,f) = NaN;
                    relfftpowerratio(blck,ch,f) = NaN;
                    relfftpeakpow(blck,ch,f) = NaN;
                else
                    diffs = fax-freqs(f);
                    f_ind = find(abs(diffs) == min(abs(diffs)));
                    freqind(blck,ch,f) = f_ind;
                    peakpow(blck,ch,f) = mean(condlfpspect(blck,ch,f,f_ind-2:f_ind+2),4);
                    blpeakpow(blck,ch,f) = mean(condlfpspect(blck,ch,1,f_ind-2:f_ind+2),4);
                    surroundpow(blck,ch,f) = mean([mean(condlfpspect(blck,ch,f,f_ind-7:f_ind-5),4),mean(condlfpspect(blck,ch,f,f_ind+5:f_ind+7),4)]);
                    meanpow(blck,ch,f) = mean(condlfpspect(blck,ch,f,1:2000),4);
                    powerratio(blck,ch,f) = peakpow(blck,ch,f)./meanpow(blck,ch,f);
                    otherratio(blck,ch,f) = peakpow(blck,ch,f)./surroundpow(blck,ch,f);
                    blratio(blck,ch,f) = peakpow(blck,ch,f)./blpeakpow(blck,ch,f);
                    
                    fftpeakpow(blck,ch,f) = condfftspect(blck,ch,f,f_ind);
                    fftblpeakpow(blck,ch,f) = condfftspect(blck,ch,1,f_ind);
                    fftsurroundpow(blck,ch,f) = mean([mean(condfftspect(blck,ch,f,f_ind-7:f_ind-5),4),mean(condfftspect(blck,ch,f,f_ind+5:f_ind+7),4)]);
                    fftmeanpow(blck,ch,f) = mean(condfftspect(blck,ch,f,1:2000),4);
                    fftpowerratio(blck,ch,f) = fftpeakpow(blck,ch,f)./fftmeanpow(blck,ch,f);
                    fftotherratio(blck,ch,f) = fftpeakpow(blck,ch,f)./fftsurroundpow(blck,ch,f);
                    fftblratio(blck,ch,f) = fftpeakpow(blck,ch,f)./fftblpeakpow(blck,ch,f);
                    
                    relpeakpow(blck,ch,f) = relspect(blck,ch,f,f_ind);
                    relmeanpow(blck,ch,f) = mean(relspect(blck,ch,f,1:2000),4);
                    relpowerratio(blck,ch,f) = relpeakpow(blck,ch,f)./relmeanpow(blck,ch,f);
                    
                    relfftpeakpow(blck,ch,f) = relfftspect(blck,ch,f,f_ind);
                    relfftmeanpow(blck,ch,f) = mean(relfftspect(blck,ch,f,1:2000));
                    relfftpowerratio(blck,ch,f) = relfftpeakpow(blck,ch,f)./relfftmeanpow(blck,ch,f);
                end
                
            end
        end
        disp('');
    end
    save(popfile, '-v7.3');
else
    load(popfile);
end

depth = depth.*(cosd(penangle)*cosd(22));
spacing = spacing.*(cosd(penangle)*cosd(22));

evaldepth = 275;
for i = 1:length(bldepth)
    for j = 1:length(electrodes(i,1):electrodes(i,2))
        cdm(i,j) = depth(i)-((j-1)*spacing(i)); % center depth matrix
    end
    [c,di(i)] = min(abs(cdm(i,:)-evaldepth)); % get depth index, least distance to evaldepth
    
    lgind(i) = lgi(i,di(i));
%     valid(i) = leffect(i,di(i));              % valid if light has a significant effect on large size stimuli
end
% valid(isnan(valid)) = 0;
% valid = logical(valid);
% vi = find(valid);