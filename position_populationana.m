function position_populationana

% % X98 ChR2
% animalids = {'170404','170501','170501','170501','170502','170502','170504','170504','170504','170531','170531','170531','170601','170601','170601','170602','170602'};
% blocks    = [3,        1,       3,       5,       1,       3,       1,       3,       5,       1,       3,       5,       1,       3,       5,       3,       5];
% stimbl    = [4,        2,       4,       6,       2,       4,       2,       4,       NaN,     2,       4,       6,       2,       4,       6,       4,       6];
% animal    = [1,        2,       2,       2,       3,       3,       4,       4,       4,       5,       5,       5,       6,       6,       6,       7,       7];
% electrodes =[[1,32];  [1,16];  [1,16];  [1,16];  [1,16];  [1,16];  [1,16];  [1,16];  [1,16];  [1,32];  [1,32];  [1,32];  [1,32];  [1,32];  [1,32];  [1,32];  [1,32]];
% penangle =  [30,       30,      30,      30,      30,      30,      30,      30,      30,      30,      30,      30,      30,      30,      30,      30,      30];
% printpath = 'C:\Users\Julia\work\data\populations\X98_ChR2_S1\position\units\';
% runprintpath = 'C:\Users\Julia\work\data\populations\X98_ChR2_S1\position\running\';
% lfpprintpath = 'C:\Users\Julia\work\data\populations\X98_ChR2_S1\position\LFP\';
% popfile = 'C:\Users\Julia\work\data\populations\X98_ChR2_S1\position\pos_population.mat';

% X94 ChR2
animalids = {'170322', '170322','170328','170414','170418','170418','170517','170517','170517','170517','170518','170518','170518','170518','170518','170523','170523','170523','170524','170524','170524'};
blocks    = [1,         2,       1,       2,       1,       3,       1,       3,       5,       7,       1,       3,       5,       7,       9,       1,       3,       5,       1,       3,       5];
stimbl    = [NaN,       NaN,     3,       3,       2,       4,       2,       4,       6,       8,       2,       4,       6,       8,       10,      2,       4,       6,       2,       4,       6];
animal    = [1,         1,       2,       3,       4,       4,       5,       5,       5,       5,       6,       6,       6,       6,       6,       7,       7,       7,       8,       8,       8];
electrodes =[[1,16];   [1,32];  [1,32];  [1,32];  [1,16];  [1,16];  [1,16];  [1,16];  [1,16];  [1,16];  [1,16];  [1,16];  [1,16];  [1,16];  [1,16];  [1,16];  [1,16];  [1,32];  [1,16];  [1,32];  [1,32]];
penangle =  [25,        25,      30,      30,      30,      30,      30,      30,      30,      30,      30,      30,      30,      30,      30,      30,      30,      30,      30,      30,      30];
printpath = 'C:\Users\Julia\work\data\populations\X94_ChR2_S1\position\units\';
runprintpath = 'C:\Users\Julia\work\data\populations\X94_ChR2_S1\position\running\';
lfpprintpath = 'C:\Users\Julia\work\data\populations\X94_ChR2_S1\position\LFP\';
popfile = 'C:\Users\Julia\work\data\populations\X94_ChR2_S1\position\pos_population.mat';

% % GIN ChR2
% animalids = {'170509','170509','170509','170509','170509'};
% blocks    = [1,        3,       5,       7,       9];
% stimbl    = [2,        4,       6,       8,       10];
% animal    = [1,        1,       1,       1,       1];
% electrodes =[[1,16];  [1,16];  [1,16];  [1,16];  [1,16]];
% penangle =  [30,       30,      30,      30,      30];
% printpath = 'C:\Users\Julia\work\data\populations\GIN_ChR2_S1\position\units\';
% runprintpath = 'C:\Users\Julia\work\data\populations\GIN_ChR2_S1\position\running\';
% lfpprintpath = 'C:\Users\Julia\work\data\populations\GIN_ChR2_S1\position\LFP\';
% popfile = 'C:\Users\Julia\work\data\populations\GIN_ChR2_S1\position\pos_population.mat';
% % 
% % X94-X98 combined ChR2 1:8 X94 9:15 X98
% animalids = {'170322', '170322','170328','170414','170418','170418','170517','170517','170517','170517','170518','170518','170518','170518','170518','170523','170523','170523','170524','170524','170524','170404','170501','170501','170501','170502','170502','170504','170504','170504','170531','170531','170531','170601','170601','170601','170602','170602'};
% blocks    = [1,         2,       1,       2,       1,       3,       1,       3,       5,       7,       1,       3,       5,       7,       9,       1,       3,       5,       1,       3,       5,       3,       1,       3,       5,       1,       3,       1,       3,       5,       1,       3,       5,       1,       3,       5,       3,       5];
% stimbl    = [NaN,       NaN,     3,       3,       2,       4,       2,       4,       6,       8,       2,       4,       6,       8,       10,      2,       4,       6,       2,       4,       6,       4,       2,       4,       6,       2,       4,       2,       4,       NaN,     2,       4,       6,       2,       4,       6,       4,       6];
% animal    = [1,         1,       2,       3,       4,       4,       5,       5,       5,       5,       6,       6,       6,       6,       6,       7,       7,       7,       8,       8,       8,       9,       10,      10,      10,      11,      11,      12,      12,      12,      13,      13,      13,      14,      14,      14,      15,      15];
% electrodes =[[1,16];   [1,32];  [1,32];  [1,32];  [1,16];  [1,16];  [1,16];  [1,16];  [1,16];  [1,16];  [1,16];  [1,16];  [1,16];  [1,16];  [1,16];  [1,16];  [1,16];  [1,32];  [1,16];  [1,32];  [1,32];  [1,32];  [1,16];  [1,16];  [1,16];  [1,16];  [1,16];  [1,16];  [1,16];  [1,16];  [1,32];  [1,32];  [1,32];  [1,32];  [1,32];  [1,32];  [1,32];  [1,32]];
% penangle =  [25,        25,      30,      30,      30,      30,      30,      30,     30,       30,      30,      30,      30,      30,      30,      30,      30,      30,      30,      30,      30,      30,      30,      30,      30,      30,      30,      30,      30,      30,      30,      30,      30,      30,      30,      30,      30,      30];
% printpath = 'C:\Users\Julia\work\data\populations\X94X98ChR2S1combined\position\units\';
% runprintpath = 'C:\Users\Julia\work\data\populations\X94X98ChR2S1combined\position\running\';
% lfpprintpath = 'C:\Users\Julia\work\data\populations\X94X98ChR2S1combined\position\LFP\';
% popfile = 'C:\Users\Julia\work\data\populations\X94X98ChR2S1combined\position\pos_population.mat';

% % 5/24 only
% animalids = {'170602','170602','170602'};
% blocks    = [1,        3,       5];
% stimbl    = [2,        4,       6];
% animal    = [1,        1,       1];
% electrodes =[[1,32];  [1,32];  [1,32]];
% penangle =  [30,       30,      30];
% printpath = 'C:\Users\Julia\work\data\populations\X98_ChR2_S1\tmp\units\';
% runprintpath = 'C:\Users\Julia\work\data\populations\X98_ChR2_S1\tmp\running\';
% lfpprintpath = 'C:\Users\Julia\work\data\populations\X98_ChR2_S1\tmp\LFP\';
% popfile = 'C:\Users\Julia\work\data\populations\X98_ChR2_S1\tmp\pos_population.mat';

% 
% % control animals
% animalids = {'170330', '170512', '170512', '170512', '170516', '170516', '170516'};
% blocks    = [1,         1,        3,        5,        1,        3,        5];
% stimbl    = [4,         2,        4,        6,        2,        4,        6];
% animal    = [1,         2,        2,        2,        3,        3,        3];
% electrodes =[[1,32];   [1,16];   [1,16];   [1,16];   [1,16];   [1,16];   [1,16]];
% penangle =  [30,        30,       30,       30,       30,       30,       30];
% printpath = 'C:\Users\Julia\work\data\populations\X94_ChR2_S1\control\units\';
% runprintpath = 'C:\Users\Julia\work\data\populations\X94_ChR2_S1\control\running\';
% lfpprintpath = 'C:\Users\Julia\work\data\populations\X94_ChR2_S1\control\LFP\';
% popfile = 'C:\Users\Julia\work\data\populations\X94_ChR2_S1\control\pos_population.mat';

tic
lcol = 'r'; %lasercolor

recalculate = 0;
printyn = 1
sfc = 1;

%Scott Kernel
%%%%%%%%%%%%kernel propeties%%%%%%%%%%%%%%%
kernel_width_eval_s = 0.15;
sdf_freq_hz = 3700;
exp_growth_ms = 2;
exp_decay_ms = 12;
exp_growth_s = exp_growth_ms/1000;
exp_decay_s = exp_decay_ms/1000;
eval_kernel_x_s = 0:1/sdf_freq_hz:kernel_width_eval_s;
exp_kernel = eval_kernel_x_s;
exp_kernel = (1-(exp(-(exp_kernel./exp_growth_s)))).*(exp(-(exp_kernel./exp_decay_s)));
excit_kernel = exp_kernel/sum(exp_kernel); 
%%%%%%%%%%kernel properties%%%%%%%%%%%%%%%
%%%%% Fieldtrip
addpath C:\Users\Julia\work\Matlab\others\fieldtrip-20160329
ft_defaults
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% chronux parameters
params.tapers = [5,9]; params.Fs = 1000; params.err = [2, 0.05]; params.trialave = 1;

freqbinwidth = 5;
sr = 1000;
prestim = 0;
poststim = 0;
respwin = 1500:2250;
respwin = respwin+prestim;
blwin = 1:1000;

if ~exist(popfile) || recalculate

    cll = 1;
    for blck = 1:length(blocks)

        supath = ['C:\Users\Julia\work\data\' animalids{blck} '\singleunits\'];
        
        
        basename = [animalids{blck} '_block' int2str(blocks(blck)) '_tet'];
        stimblname = [animalids{blck} '_block' int2str(stimbl(blck)) '_tet'];

        files = dir([supath, basename, '*.mat']);
        stimfiles = dir([supath, stimblname, '*.mat']);
        
        for fi = 1:length(files)

            if strfind(files(fi).name, 'MU')
                mu(cll) = 1; su(cll) = 0;
            else
                mu(cll) = 0; su(cll) = 1;
            end
            
            try
                [stimresp, stimlfpspect, stimfax] = get_stimresp([supath, stimfiles(fi).name]);
            catch
            end
            
            load([supath, files(fi).name]);            
            
            i = strfind(files(fi).name, 'tet');
            if strcmp(files(fi).name(i+4),'_')
                tetno = strread(files(fi).name(i+3)); % single character number
            else
                tetno = strread(files(fi).name(i+3:i+4)); % number >10
            end
            if tetno*4<electrodes(blck,1) || tetno*4>electrodes(blck,2) % assure we're only getting V1 in the population
                continue;
            end
            
            cllname{cll} = files(fi).name;
            printname = files(fi).name;
            printname(find(printname=='_')) = ' ';
                        
            % important stuff
            depth(cll) = result.depth;
            pangle(cll) = penangle(blck);
            recording(cll) = blck;
            animalno(cll) = animal(blck);
            
            light = result.randconds(2,:); position = result.randconds(1,:);
            
            wvchan = find(var(result.waveforms) == max(var(result.waveforms)));
            lfp = result.lfp(:,wvchan)';
            spike = result.waveforms(:,wvchan);
            interpspike = spline(1:32,spike,1:.1:32);
            [adiff(cll),swidth(cll)] = spikequant(interpspike);            
            secpersamp = 1/30000;
            interpf = secpersamp/10;
            swidthms(cll) = swidth(cll)*interpf*1000;
            waveform(cll,:) = spike;
            clustqual(cll) = result.clusterquality;

            % get spiketimes
            msStimes = round(result.spikes);
            if ~isempty(msStimes) && msStimes(1) == 0, msStimes(1) = 1; end

            chan = zeros(1,length(result.lfp));
            chan(msStimes) = 1;    
            
            isi = diff(msStimes);
            bursts = legendy_new3(isi,4,1000,3,0,15); %(ISI, fac, sr, min_length_of_burst, local_length, surprise_cutoff)
            surp = [bursts.surprise];
            bursts(surp == 100) = []; % delete probably wrong bursts
            burstbegs = [bursts.begin];
            burstchan = zeros(1,length(result.lfp));
            burstchan(msStimes(burstbegs)) = 1;
            
            gamma = eegfilt(lfp,sr,15,35);  % adjust depending on peak
            h = hilbert(gamma); gpow = abs(h); gphas = angle(h);
            % %
            if sfc
                clear filtmat; clear powmat; clear phasmat;
                for i = 1:100/freqbinwidth
                    filtmat(i,:) = eegfilt(lfp,sr,(i-1)*freqbinwidth+1,i*freqbinwidth);
                    h = hilbert(filtmat(i,:));
                    powmat(i,:) = abs(h); phasmat(i,:) = angle(h);
                end
            end
            
            trialdur = result.sweeplength;
            msstamps = result.msstamps;
            
            if strcmp(animalids{blck}, '170328')
                result.randconds = result.randconds(:,1:396);
                msstamps = msstamps(1:396); light = light(1:396);
                position = position(1:396);
            end
            if length(msstamps)~=size(result.randconds,2)
                disp('');
                %         msstamps([518]) = []; 
                %         result.msstamps = msstamps;
                %         save([supath, files(fi).name],'result');
                pause;
            end
            
            clear speed;
            for i = 1:length(msstamps)
                speed(i,:) = result.runspeed(msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim);
            end
       
            % figure out sufficiently high and nonvariable runspeed trials
            meanspeed = mean(speed(:,respwin),2);
            stdspeed = std(speed(:,respwin),1,2);
            notstill = find(meanspeed>1);
            okspeed = find(meanspeed>( mean(meanspeed(notstill))-(1.5*std(meanspeed(notstill))) ) & meanspeed>1);
            okvar = find(stdspeed<( mean(stdspeed(notstill))+(1.5*std(stdspeed(notstill)))) & stdspeed>.5);
            oktrials = intersect(okspeed,okvar);
            nonoktrials = 1:size(speed,1); nonoktrials(oktrials) = [];
            stilltrials = 1:size(speed,1); stilltrials(notstill) = [];
            
            
            nrun(cll) = length(oktrials);
            nstill(cll) = length(stilltrials);
            if ~isempty(oktrials)
                cellmeanspeed(cll) = mean(meanspeed(oktrials));
                cellminspeed(cll) = min(meanspeed(oktrials));
            else
                cellmeanspeed(cll) = NaN;
                cellminspeed(cll) = NaN;
            end
            
            clear resp; clear lfpresp; clear burstresp; clear msl; clear lfpspect; clear gammaresp;
            clear gammapowresp; clear gphaseresp; clear trialstalfp; clear nspkstrialsta; clear spkxcorr;
            clear allphaseresp; clear allpowresp;
            for i = 1:length(msstamps)
                resp(i,:) = chan(msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim);
                lfpresp(i,:) = lfp(msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim); 
                burstresp(i,:) = burstchan(msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim); 
                msl(i) = mean(find(resp(i,respwin)));
                [lfpspect(i,:),trialfax] = pmtm(lfpresp(i,respwin(1)+200:respwin(end)),3,[],sr);
%                 [lfpsectrogram(i,:,:),ct,cf] = mtspecgramc(lfpresp(i,:)',[.25,.1],params);
                gammaresp(i,:) = gamma(msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim);
                gammapowresp(i,:) = gpow(msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim);
                gphaseresp(i,:) = gphas(msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim);
                
                % Spike triggered LFP
                s = find(resp(i,respwin));
                if ~isempty(s)
                    for si = 1:length(s)
                        lfps(si,:) = lfpresp(i,s(si)+respwin(1)-1-200:s(si)+respwin(1)-1+200);
                    end
                else
                    lfps = nan(1,401);
                end
                trialstalfp(i,:) = mean(lfps,1);
                nspkstrialsta(i) = size(lfps,1);
                clear lfps;
                [spkxcorr(i,:),lags] = xcorr(resp(i,respwin),resp(i,respwin));
                
                % for all frequency bands
                if sfc
                    for j = 1:size(phasmat,1)
                        allphaseresp(j,i,:) = phasmat(j, msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim);
                        allpowresp(j,i,:) = powmat(j, msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim);
                    end
                end
            end            
            
            msta = linspace(-blwin(end),trialdur-blwin(end),size(resp,2));
            
            frs = sum(resp(:,respwin),2)./(length(respwin)/1000);
            bcs = sum(burstresp(:,respwin),2);
            bl = sum(resp(:,blwin),2)./(length(blwin)/1000);
            sc = sum(resp(:,respwin),2); % spike count   
            burstsperspike(cll) = mean(bcs)/mean(sc);
            
%             %determine if cll is touch modulated
%             touchfr = sum(resp(:,blwin(end)+40:blwin(end)+40+200),2);
%             shortbl = sum(resp(:,blwin(end)-40-200:blwin(end)-40),2);
%             touchmod(cll) = ttest2(shortbl,touchfr);
            okcond = zeros(1,size(resp,1)); okcond(oktrials) = 1;
            touchfr = sum(resp(okcond&light == 0,blwin(end)+151:blwin(end)+150+200),2);
            shortbl = sum(resp(okcond&light == 0,blwin(end)-149-200:blwin(end)-150),2);
            if ~isempty(find(okcond))
                [p, touchmod(cll)] = signrank(shortbl,touchfr);
            else
                touchmod(cll) = 0;
            end
            
            %determine if cll is modulated by light
            lightmod(cll) = ranksum(frs(find(light)),frs(find(light == 0)));            
            
            % anova over all things
            [p,table,stats] = anovan(frs,{result.randconds(2,:), result.randconds(1,:)},'model','full','display','off');
            falp(cll) = p(1); fapp(cll) = p(2); falpip(cll) = p(3); % full anova light p, position p, light pos interaction p
            
            [p,table,stats] = anova1(frs(light == 0&position~=0&okcond),result.randconds(1,light == 0&position~=0&okcond),'off')
            tunep(cll) = p;
            
            lfr(cll) = mean(frs(find(light)));
            nlfr(cll) = mean(frs(find(light == 0)));
            
            binwidth = 20;
            [binnedlight,bta] = binit(mean(resp(find(light),:)),binwidth);
            binnedlight = binnedlight.*(1000/binwidth);
            [binnednolight,bta] = binit(mean(resp(find(light==0),:)),binwidth);
            binnednolight = binnednolight.*(1000/binwidth);
            
            clear stalfp; clear nspkssta; clear shufstalfp; clear nspksshufsta;
            for l = 1:length(unique(result.light))
                for pos = 1:length(result.positions)
                    trials = find(position == result.positions(pos) & light == l-1);
                    
                    condresp(cll,l,pos,:) = mean(resp(trials,:),1);
                    condresperr(cll,l,pos,:) = nanstd(resp(trials,:),1,1)./sqrt(length(trials));
                    condlfpresp(cll,l,pos,:) = mean(lfpresp(trials,:),1);
                    
                    condlfpspect(cll,l,pos,:) = nanmean(lfpspect(trials,:));
                    [S,chf,Serr]=mtspectrumc(squeeze(lfpresp(trials,1750:2250))',params);
                    condS(cll,l,pos,:) = S(1:60); condSerr(cll,l,pos,:,:) = Serr(:,1:60);
                    
                    condfr(cll,l,pos) = mean(frs(trials));%-mean(bl);
                    conderr(cll,l,pos) =std(frs(trials))./sqrt(length(trials));
                    condbc(cll,l,pos) = mean(bcs(trials));
                    
                    condmsl(cll,l,pos) = nanmean(msl(trials));
                    condmslerr(cll,l,pos) = nanstd(msl(trials))./sqrt(length(trials));
                    
                    condlfpsta(cll,l,pos,:) = nanmean(trialstalfp(trials,:));
                    condavgnspks(cll,l,pos) = nanmean(nspkstrialsta(trials));
                    condspkxcorr(cll,l,pos,:) = nanmean(spkxcorr(trials,:));
                    
                    condz{cll,l,pos,:} = (sc(trials)-mean(sc(trials)))/std(sc(trials)); %ecker 2010
                    condsc{cll,l,pos,:} = sc(trials);
                    ff(cll,l,pos) = var(sc(trials))/mean(sc(trials));
                    
                    hlp = ones(size(condresp,4),1);
                    [xx,bta] = binit(hlp,binwidth);
                    bincondresp(cll,l,pos,:) = binit(condresp(cll,l,pos,:),binwidth).*(1000/binwidth);
                    binconderr(cll,l,pos,:) = binit(condresperr(cll,l,pos,:),binwidth).*(1000/binwidth);
                    
                    mscc = []; bincc = [];
                    for ii = 1:length(trials)-1
                        for jj = ii+1:length(trials)
                            help = corrcoef(resp(trials(ii),:),resp(trials(jj),:));
                            mscc = [mscc,help(1,2)];
                            help = corrcoef(binit(resp(trials(ii),:),binwidth),binit(resp(trials(jj),:),binwidth));
                            bincc = [bincc, help(1,2)];
                        end
                    end
                    msreliab(cll,l,pos) = nanmean(mscc);
                    binreliab(cll,l,pos) = nanmean(bincc);
                    eckerreliability(cll,l,pos) = var(frs(trials))/var(frs);
                    
                    thisruninds = intersect(trials,oktrials);
                    if ~isempty(thisruninds)
                        thisrunn(cll,l,pos) = length(thisruninds);
                        runcondresp(cll,l,pos,:) = mean(resp(thisruninds,:),1);
                        runcondresperr(cll,l,pos,:) = nanstd(resp(thisruninds,:),1,1)./sqrt(length(thisruninds));
                        runcondfr(cll,l,pos) = mean(frs(thisruninds));
                        runconderr(cll,l,pos) = std(frs(thisruninds))./sqrt(length(thisruninds));
                        runbincondresp(cll,l,pos,:) = binit(runcondresp(cll,l,pos,:),binwidth).*(1000/binwidth);
                        runbinconderr(cll,l,pos,:) = binit(runcondresperr(cll,l,pos,:),binwidth).*(1000/binwidth);
                        [runS,chf,runSerr] = mtspectrumc(squeeze(lfpresp(thisruninds,1750:2250))',params);
                        condrunS(cll,l,pos,:) = runS(1:60); condrunSerr(cll,l,pos,:,:) = runSerr(:,1:60);
                    else
                        thisrunn(cll,l,pos) = 0;
                        runcondresp(cll,l,pos,:) = nan(1,size(resp,2));
                        runcondresperr(cll,l,pos,:) = nan(1,size(resp,2));
                        runcondfr(cll,l,pos) = NaN;
                        runconderr(cll,l,pos) = NaN;
                        runbincondresp(cll,l,pos,:) = nan(1,length(bta));
                        runbinconderr(cll,l,pos,:) = nan(1,length(bta));
                        condrunS(cll,l,pos,:) = nan(1,60); condrunSerr(cll,l,pos,:,:) = nan(2,60);
                    end
                    
                    thisstillinds = intersect(trials,stilltrials);
                    if ~isempty(thisstillinds)
                        thisstilln(cll,l,pos) = length(thisstillinds);
                        stillcondresp(cll,l,pos,:) = mean(resp(thisstillinds,:),1);
                        stillcondresperr(cll,l,pos,:) = nanstd(resp(thisstillinds,:),1,1)./sqrt(length(thisstillinds));
                        stillcondfr(cll,l,pos) = mean(frs(thisstillinds));
                        stillconderr(cll,l,pos) = std(frs(thisstillinds))./sqrt(length(thisstillinds));
                        stillbincondresp(cll,l,pos,:) = binit(stillcondresp(cll,l,pos,:),binwidth).*(1000/binwidth);
                        stillbinconderr(cll,l,pos,:) = binit(stillcondresperr(cll,l,pos,:),binwidth).*(1000/binwidth);
                        [stillS,chf,stillSerr] = mtspectrumc(squeeze(lfpresp(thisstillinds,1750:2250))',params);
                        condstillS(cll,l,pos,:) = stillS(1:60); condstillSerr(cll,l,pos,:,:) = stillSerr(:,1:60);
                    else
                        thisstilln(cll,l,pos) = 0;
                        stillcondresp(cll,l,pos,:) = nan(1,size(resp,2));
                        stillcondresperr(cll,l,pos,:) = nan(1,size(resp,2));
                        stillcondfr(cll,l,pos) = NaN;
                        stillconderr(cll,l,pos) = NaN;
                        stillbincondresp(cll,l,pos,:) = nan(1,length(bta));
                        stillbinconderr(cll,l,pos,:) = nan(1,length(bta));
                        condstillS(cll,l,pos,:) = nan(1,60); condstillSerr(cll,l,pos,:,:) = nan(2,60);
                    end
                    
                    % field trip shitz fo SFC spectra
                    clear data;
                    for i = 1:length(trials)
                        data.trial{i} = [lfpresp(trials(i),:);resp(trials(i),:)];
                        time = 1:result.sweeplength;
                        data.time{i} = time(1)./1000:.001:time(end)./1000;
                        %                 data.time{i} = -.299:.001:2.700;
                    end
                    data.fsample = 1000;
                    data.trialinfo(:,1) = light(trials);
                    data.trialinfo(:,2) = position(trials);
                    cfg = [];
                    cfg.timwin = [-.25 .25];
                    data.label{1,1} = 'lfp';
                    data.label{2,1} = 'spikes';
                    cfg.spikechannel = 'spikes';
                    cfg.channel = 'lfp';
                    cfg.latency = [0.5,1.5];
                    
                    cfg.method = 'mtmfft';
                    cfg.foilim = [5,100];
                    cfg.timwin = [-.15, .15];
                    cfg.taper = 'hanning';
                    cfg.spikechannel = 'spikes';
                    cfg.channel = 'lfp';
                    stsFFT           = ft_spiketriggeredspectrum(cfg, data);
                    ang = angle(stsFFT.fourierspctrm{1});
                    mag = abs(stsFFT.fourierspctrm{1});
                    ftspect(cll,l,pos,:) = squeeze(nanmean(mag(:,1,:)));
                    ftphases{cll,l,pos} = squeeze(ang);
                    ftfax = stsFFT.freq;
                    
                    cfg               = [];
                    cfg.method        = 'ral'; % compute the rayleigh test
                    cfg.spikechannel  = stsFFT.label{1};
                    cfg.channel       = stsFFT.lfplabel; % selected LFP channels
                    cfg.avgoverchan   = 'unweighted'; % weight spike-LFP phases irrespective of LFP power
                    cfg.timwin        = 'all'; % compute over all available spikes in the window
                    cfg.latency       = [0.5 1.5]; % sustained visual stimulation period
                    statSts           = ft_spiketriggeredspectrum_stat(cfg,stsFFT);
                    ftralp(cll,l,pos,:) = statSts.ral;
                    cfg.method = 'ppc0';
                    statSts           = ft_spiketriggeredspectrum_stat(cfg,stsFFT);
                    ftppc(cll,l,pos,:) = statSts.ppc0;
                    cfg.method = 'plv';
                    statSts           = ft_spiketriggeredspectrum_stat(cfg,stsFFT);
                    ftplv(cll,l,pos,:) = statSts.plv;
                end
            end
            sfx = chf(1:60);
            
            ta = bta-blwin(end);
            
            %  figure
            clf;
            subplot(2,2,1)
            plot(ta,squeeze(mean(runbincondresp(cll,1,:,:),3)),'k','linewidth',2);
            hold on
            plot(ta,squeeze(mean(runbincondresp(cll,2,:,:),3)),lcol,'linewidth',2);
            mx = max([max(squeeze(mean(runbincondresp(cll,1,:,:),3))),max(squeeze(mean(runbincondresp(cll,2,:,:),3))),.01]);
            axis([-blwin(end),trialdur-blwin(end),0,mx]);
            line([0,0],[0,mx],'color','k','linewidth',2);
            line([1500,1500],[0,mx],'color','k','linewidth',2);
            line([500,500],[0,mx],'color','b','linewidth',2)
            line([1250,1250],[0,mx],'color','b','linewidth',2);
            legend({'Light OFF','Light ON'})
            xlabel('time [ms]')
            ylabel('firing rate [Hz]')
            title(['cll ' int2str(cll) ' depth: ' int2str(depth(cll)), 'cll ' cllname{cll} ])
            
            subplot(2,2,2)
            errorbar(result.positions,squeeze(runcondfr(cll,1,:)),squeeze(runconderr(cll,1,:)),'ko-','markerfacecolor','k','linewidth',2)
            hold on
            errorbar(result.positions,squeeze(runcondfr(cll,2,:)),squeeze(runconderr(cll,2,:)),'ro-','markerfacecolor','r','linewidth',2)
            ax = axis;
            axis([-1,9,ax(3),ax(4)]);
            legend('light off','light on')
            xlabel('bar position')
            ylabel('firing rate [Hz]')
            
            subplot(2,2,3)
            plot(ta,squeeze(mean(stillbincondresp(cll,1,:,:),3)),'k','linewidth',2);
            hold on
            plot(ta,squeeze(mean(stillbincondresp(cll,2,:,:),3)),lcol,'linewidth',2);
            mx = max([max(squeeze(mean(stillbincondresp(cll,1,:,:),3))),max(squeeze(mean(stillbincondresp(cll,2,:,:),3))),.01]);
            axis([-blwin(end),trialdur-blwin(end),0,mx]);
            line([0,0],[0,mx],'color','k','linewidth',2);
            line([1500,1500],[0,mx],'color','k','linewidth',2);
            line([500,500],[0,mx],'color','b','linewidth',2)
            line([1250,1250],[0,mx],'color','b','linewidth',2);
            legend({'Light OFF','Light ON'})
            xlabel('time [ms]')
            ylabel('firing rate [Hz]')
            title('non-running')
            
            subplot(4,4,11)
            plot(waveform(cll,:))
            axis([0,40,-100,100])
            legend(['width: ' num2str(swidthms(cll))])
            
            try
            subplot(4,4,12)            
            hold on;
            for i = 1:size(stimresp,1)
                if ~isempty(find(stimresp(i,:)))
                    plot(find(stimresp(i,:)),i,'ko','MarkerSize',1.5,'MarkerFaceColor','k')
                end
            end
            % xlim([0 size(resp,2)]),ylim([0 size(resp,1)+1]),set(gca,'visible','off','Xtick',1:1:size(resp,1));
            
            binwidth = 10;            
            subplot(4,4,16)
            [binr,bta] = binit(mean(stimresp,1),binwidth); binr = binr.*(1000/binwidth);
            plot(bta,binr,'k','LineWidth',1.5)
            
            subplot(4,4,15)
            semilogy(stimfax,nanmean(stimlfpspect));
            axis([0,30,min(nanmean(stimlfpspect(:,1:124))),max(nanmean(stimlfpspect))])
            catch
            end
                
            if printyn
                figSize = [30 21];
                set(gcf,'paperunits','centimeters','papersize',figSize,'paperposition',[0 0 figSize])
                if cll<10, printi = ['0', int2str(cll)]; else printi = int2str(cll); end
                print([printpath ,  printi '__' files(fi).name '.pdf'],'-dpdf')
            end
%             
%             clf; 
%             for i = [1,2,4,6,8]
%                 subplot(1,9,i)
%                 plot(ta,squeeze(runbincondresp(cll,1,i,:)),'k','linewidth',2);
%                 hold on                
%                 plot(ta,squeeze(runbincondresp(cll,2,i,:)),'r','linewidth',2);
%             end
            
            % running figure
            runlfr(cll) = mean(frs(intersect(find(result.randconds(2,:)),oktrials)));
            runnlfr(cll) = mean(frs(intersect(find(~result.randconds(2,:)),oktrials)));
            norunlfr(cll) = mean(frs(intersect(find(result.randconds(2,:)),stilltrials)));
            norunnlfr(cll) = mean(frs(intersect(find(~result.randconds(2,:)),stilltrials)));
            lfrerr = std(frs(find(result.randconds(2,:))))./sqrt(length(find(result.randconds(2,:))));
            nlfrerr = std(frs(find(~result.randconds(2,:))))./sqrt(length(find(~result.randconds(2,:))));
            runlfrerr = std(frs(intersect(find(result.randconds(2,:)),oktrials)))./sqrt(length(intersect(find(result.randconds(2,:)),oktrials)));
            runnlfrerr = std(frs(intersect(find(~result.randconds(2,:)),oktrials)))./sqrt(length(intersect(find(~result.randconds(2,:)),oktrials)));
            norunlfrerr = std(frs(intersect(find(result.randconds(2,:)),stilltrials)))./sqrt(length(intersect(find(result.randconds(2,:)),stilltrials)));
            norunnlfrerr = std(frs(intersect(find(~result.randconds(2,:)),stilltrials)))./sqrt(length(intersect(find(~result.randconds(2,:)),stilltrials)));
            
            l1r1 = frs(intersect(find(result.randconds(2,:)),oktrials));
            l0r1 = frs(intersect(find(~result.randconds(2,:)),oktrials));
            l1r0 = frs(intersect(find(result.randconds(2,:)),stilltrials));
            l0r0 = frs(intersect(find(~result.randconds(2,:)),stilltrials));
            anovavec = [l0r0;l0r1;l1r0;l1r1];
            g1 = [zeros(length(l0r0),1);zeros(length(l0r1),1);ones(length(l1r0),1);ones(length(l1r1),1)]; %light
            g2 = [zeros(length(l0r0),1);ones(length(l0r1),1);zeros(length(l1r0),1);ones(length(l1r1),1)]; %running
            [p,table,stats] = anovan(anovavec,{g1 g2},'model','full','display','off');
            lp(cll) = p(1); rp(cll) = p(2); rlip(cll) = p(3);
            
            clf;
            subplot(2,2,1)
            imagesc(speed);
            colorbar
            title(['oktrials: ' int2str(length(oktrials)) '/' int2str(size(speed,1))])
            xlabel('time [ms]')
            ylabel('trial number')
            
            subplot(2,2,2)
            errorbar(msta,mean(speed(find(~light),:)),std(speed(find(~light),:))./sqrt(length(find(~light))),'b')
            hold on
            errorbar(msta,mean(speed(find(light),:)),std(speed(find(light),:))./sqrt(length(find(light))),'r')
            xlabel('time [ms]')
            ylabel('average runspeed')
            legend({'light off' 'light on'})
            
            subplot(2,2,3)
            plot(mean(speed(:,respwin),2),frs,'.')
            hold on
            plot(mean(speed(find(light),respwin),2),frs(find(light)),'r.')
            xlabel('average runspeed of trial')
            ylabel('average firing rate of trial')
            
            subplot(2,2,4)
            barweb([squeeze(mean(condfr(cll,1,:),3)),squeeze(mean(condfr(cll,2,:)));...
                squeeze(mean(runcondfr(cll,1,:),3)),squeeze(mean(runcondfr(cll,2,:)));...
                squeeze(mean(stillcondfr(cll,1,:),3)),squeeze(mean(stillcondfr(cll,2,:)))],...
                [squeeze(mean(conderr(cll,1,:),3)),squeeze(mean(conderr(cll,2,:)));...
                squeeze(mean(runconderr(cll,1,:),3)),squeeze(mean(runconderr(cll,2,:)));...
                squeeze(mean(stillconderr(cll,1,:),3)),squeeze(mean(stillconderr(cll,2,:)))],...
                [],[{'all'};{'running only'};{'immobile only'}],[],...
                [],'firing rate [Hz]',[],[]);
            
            if printyn
                figSize = [30 21];
                set(gcf,'paperunits','centimeters','papersize',figSize,'paperposition',[0 0 figSize])
                if cll<10, printi = ['0', int2str(cll)]; else printi = int2str(cll); end
                print([runprintpath ,  printi '__' files(fi).name '.pdf'],'-dpdf')
            end
            
            clf;% LFP figure
            fillx = [sfx,fliplr(sfx)];
            fillyrun0 = [squeeze(nanmean(condrunSerr(cll,1,:,1,:),3))',fliplr(squeeze(nanmean(condrunSerr(cll,1,:,2,:)))')];
            fillyrun1 = [squeeze(nanmean(condrunSerr(cll,2,:,1,:),3))',fliplr(squeeze(nanmean(condrunSerr(cll,2,:,2,:)))')];
            fillystill0 = [squeeze(nanmean(condstillSerr(cll,1,:,1,:),3))',fliplr(squeeze(nanmean(condstillSerr(cll,1,:,2,:)))')];
                        
            subplot(2,2,1)
            fill(fillx,fillyrun0,'b')
            hold on
            fill(fillx,fillyrun1,'r')
            set(gca,'yscale','log')
            xlabel('frequency [Hz]')
            ylabel('power')
            legend('L0','L1')
            
            subplot(2,2,2)
            fill(fillx,fillyrun0,'c')
            hold on
            fill(fillx,fillystill0,'b')
            set(gca,'yscale','log')
            xlabel('frequency [Hz]')
            ylabel('power')
            legend('running','still')
            
            subplot(2,2,3)
            plot(ftfax,squeeze(nanmean(ftppc(cll,1,:,:),3)),'b')
            hold on
            plot(ftfax,squeeze(nanmean(ftppc(cll,2,:,:),3)),'r')
            xlabel('frequency')
            ylabel('PPC')
            legend('L0','L1')
            
            if printyn
                figSize = [30 21];
                set(gcf,'paperunits','centimeters','papersize',figSize,'paperposition',[0 0 figSize])
                if cll<10, printi = ['0', int2str(cll)]; else printi = int2str(cll); end
                print([lfpprintpath ,  printi '__' files(fi).name '.pdf'],'-dpdf')
            end
            
            disp([files(fi).name '   done'])
            cll = cll + 1;
            
        end
    end
    save(popfile, '-v7.3');
else
    load(popfile);   
end

toc

clear cell;
%spike classification
% kmeansind = kmeans([eslope',ptr',swidth',adiff'],4);
kmeansind = kmeans([swidth',adiff'],3);

m1 = mean(swidth(kmeansind==1)); m2 = mean(swidth(kmeansind==2)); m3 = mean(swidth(kmeansind==3));
fm = min([m1,m2,m3]); fmax = max([m1,m2,m3]);
if fm == m1
    pfs = find(kmeansind==1); pfsv = kmeansind==1;
elseif fm == m2
    pfs = find(kmeansind==2); pfsv = kmeansind == 2;
else
    pfs = find(kmeansind==3); pfsv = kmeansind == 3;
end
if fmax == m1
    prs = find(kmeansind==1); prsv = kmeansind==1;
elseif fmax == m2
    prs = find(kmeansind==2); prsv = kmeansind == 2;
else
    prs = find(kmeansind==3); prsv = kmeansind == 3;
end

% m1 = mean(swidth(kmeansind==1)); m2 = mean(swidth(kmeansind==2)); m3 = mean(swidth(kmeansind==3));
% if find(kmeansind == 4), m4 = mean(swidth(kmeansind == 4)); end
% fm = min([m1,m2,m3,m4]);
% if fm == m1
%     pfs = find(kmeansind==1); prs = find(kmeansind==2|kmeansind==3|kmeansind==4); 
%     pfsv = kmeansind==1; prsv = kmeansind==2|kmeansind==3|kmeansind==4;
% elseif fm == m2
%     pfs = find(kmeansind==2); prs = find(kmeansind==1|kmeansind==3|kmeansind==4);
%     pfsv = kmeansind==2; prsv = kmeansind==1|kmeansind==3|kmeansind==4;
% elseif fm == m3
%     pfs = find(kmeansind==3); prs = find(kmeansind==1|kmeansind==2|kmeansind==4);
%     pfsv = kmeansind==3; prsv = kmeansind==1|kmeansind==2|kmeansind==4;
% else
%     pfs = find(kmeansind==4); prs = find(kmeansind==1|kmeansind==2|kmeansind==3);
%     pfsv = kmeansind==4; prsv = kmeansind==1|kmeansind==2|kmeansind==3;
% end

% if mean(swidth(find(kmeansind==1)))<mean(swidth(find(kmeansind==2)))  %1 is FS
%     pfs = find(kmeansind==1); prs = find(kmeansind==2); pfsv = kmeansind==1;
% else
%     pfs = find(kmeansind==2); prs = find(kmeansind==1); pfsv = kmeansind==2;
% end

figure
plot(swidthms(kmeansind==1),adiff(kmeansind==1),'b.')
xlabel('spike width')
ylabel('amplitude diff')
hold on
plot(swidthms(kmeansind==2),adiff(kmeansind==2),'r.')
if ~isempty(find(kmeansind==3))
    plot(swidthms(kmeansind==3),adiff(kmeansind==3),'g.')
    plot(swidthms(pfsv),adiff(pfsv),'ro');
    plot(swidthms(prsv),adiff(prsv),'o');
end

% figure
% plot(eslope(kmeansind==1),adiff(kmeansind==1),'b.')
% xlabel('end slope')
% ylabel('amplitude diff')
% hold on
% plot(eslope(kmeansind==2),adiff(kmeansind==2),'r.')
% if ~isempty(find(kmeansind==3))
%     plot(eslope(kmeansind==3),adiff(kmeansind==3),'g.')
%     plot(eslope(pfsv),adiff(pfsv),'ro');
%     plot(eslope(prsv),adiff(prsv),'o');
% end
%     % axis([5.5,20.5,-.9,.7])

prsv = swidthms>.36; pfsv = swidthms<=.36;
prs = find(prsv); pfs = find(pfsv);

swamp = max(waveform,[],2)-min(waveform,[],2);
okwv = swamp>42;

for i = 1:length(okwv)
    isodist(i) = clustqual(i).IsolationDistance;
    lratio(i) = clustqual(i).L_Ratio.Lratio;
end

% ok = vismod&okwv'&nlfr>1;
% ok = okwv'&isodist>8;
ok = okwv';

% adjust depth according to penetration angle
% pangle = 10; %TODO temp fix
% depth = depth.*cosd(22).*cosd(pangle);
depth = depth.*cosd(45-pangle);

phe = zeros(1,length(depth));

% %putative halo expressing for X94 ChR2
phe([10, 20, 67,80, 202, 240]) = 1; 

% % putative halo expressing for X98 ChR2
% % phe([2, 10, 38, 71, 90, 129]) = 1; 
% phe([10, 71, 178, 182, 206, 221, 242, 261]) = 1; 

% % putative halo expressing for GIN ChR2
% phe([1,31,38, 59, 90]) = 1; 

% % putative halo for X94 X98 combined
% phe([10, 20, 67,80, 202, 240,     308, 336, 476, 480, 504, 519, 540, 559]) = 1; 

phe = logical(phe);

[sd,si] = sort(depth);
[rsd,rsi] = sort(depth(prs));
[fsd,fsi] = sort(depth(pfs));

l23 = depth<375;
l4 = depth>=375&depth<=550;
% l5 = depth>550&depth<=800;
% l5a = depth>550&depth<=650;
% l5b = depth>650&depth<=800;
l5 = depth>550&depth<=925;
l5a = depth>550&depth<=650;
l5b = depth>650&depth<=925;
l6 = depth>800;
l23rs = l23&prsv&~phe&ok;
l23fs = l23&pfsv&~phe&ok;
l4rs = l4&prsv&~phe&ok;
l4fs = l4&pfsv&~phe&ok;
l5rs = l5&prsv&~phe&ok;
l5fs = l5&pfsv&~phe&ok;
l6rs = l6&prsv&~phe&ok;
l6fs = l6&pfsv&~phe&ok;
okrs = prsv&~phe&ok; 

for i = 1:length(depth)
    mfr(i) = mean(runcondfr(i,1,2:9),3);
    omi(i) = (nanmean(runcondfr(i,2,2:9),3)-nanmean(runcondfr(i,1,2:9),3))./(nanmean(runcondfr(i,2,2:9),3)+nanmean(runcondfr(i,1,2:9),3));
    ominr(i) = (nanmean(stillcondfr(i,2,2:9),3)-nanmean(stillcondfr(i,1,2:9),3))./(nanmean(stillcondfr(i,2,2:9),3)+nanmean(stillcondfr(i,1,2:9),3));
    deltaspikes(i) = nanmean(runcondfr(i,2,2:9),3)-nanmean(runcondfr(i,1,2:9),3);
    dss(i) = squeeze(runcondfr(i,2,1))-squeeze(runcondfr(i,1,1)); % delta spikes spontaneous
    for l = 1:2
        ssi(i,l) = get_ssi(squeeze(runcondfr(i,l,2:8))); % spatial selectivity
    end
    for pos = 1:9
        omicurve(i,pos) = (runcondfr(i,2,pos)-runcondfr(i,1,pos))./(runcondfr(i,2,pos)+runcondfr(i,1,pos));
        omicurvenr(i,pos) = (stillcondfr(i,2,pos)-stillcondfr(i,1,pos))./(stillcondfr(i,2,pos)+stillcondfr(i,1,pos));
        deltacurve(i,pos) = runcondfr(i,2,pos)-runcondfr(i,1,pos);
    end
    [ranked(i,1,:),si] = sort(squeeze(runcondfr(i,1,2:9)),'descend');
    ranked(i,2,:) = squeeze(runcondfr(i,2,si+1));
    normranked(i,:,:) = ranked(i,:,:)./ranked(i,1,1);
    rankedomi(i,:) = (ranked(i,2,:)-ranked(i,1,:))./(ranked(i,2,:)+ranked(i,1,:));
    [linp(i,:),S] = polyfit(squeeze(normranked(i,1,:)),squeeze(normranked(i,2,:)),1);
    normr(i) = S.normr;
    c = corrcoef(squeeze(runcondfr(i,1,:)),squeeze(runcondfr(i,2,:)));
    rsq(i) = c(1,2)^2;
    for l = 1:2
        for pos = 1:9
            efr = (sum(condresp(i,l,pos,1501:2250)).*1000)/750;
            bfr = (sum(condresp(i,l,pos,251:1000)).*1000)/750;
            blsfr(i,l,pos) = efr-bfr;
        end
    end
end


cond4 = l4&~phe&~mu&prsv;%&nlfr>1&lightmod ;
cond5 = l5&~phe&~mu&prsv;%&nlfr>1&lightmod;
[p,s] = ranksum(omicurve(cond4,1),omicurve(cond5,1));
[p,s] = ranksum(omi(cond4),omi(cond5));


% combined population histograms of multiplicative and additive parts
% only for combined pop! adjust as needed
x98 = zeros(1,587); x98(297:end) = 1; % adjust
x94 = zeros(1,587); x94(1:296) = 1;

figure
cond = cond5;
bins = [0:.1:1.1]
[y,x] = histc(linp(cond&x98&rsq>0.25,1),bins);
% [y,x] = histc(linp(cond4&rsq>0.3,1),bins);
plot(bins,y,'m','linewidth',2)
hold on
[y,x] = histc(linp(cond&x94&rsq>0.25,1),bins);
% [y,x] = histc(linp(cond5&rsq>0.3,1),bins);
plot(bins,y,'c','linewidth',2)
legend('X98','X94')
% legend('L4','L5')
xlabel('multiplication factor')
ylabel('count')
title('Layer 5: histogram of multiplicative factors')

[p,s] = ranksum(linp(cond5&x94&rsq>0.25,1),linp(cond5&x98&rsq>0.25,1))

figure
cond = cond5;
bins = [-.6:.1:.6];
[y,x] = histc(linp(cond&x98&rsq>0.25,2),bins);
plot(bins,y,'m','linewidth',2)
hold on
[y,x] = histc(linp(cond&x94&rsq>0.25,2),bins);
plot(bins,y,'c','linewidth',2)
legend('X98','X94')
xlabel('y offset')
ylabel('count')
title('Layer 5: histogram of y offsets')

figure
plot(ssi(cond5&x94,1),ssi(cond5&x94,2),'co','markerfacecolor','c')
hold on
plot(ssi(cond5&x98,1),ssi(cond5&x98,2),'bo','markerfacecolor','b')
axis([0,1,0,1])
refline(1,0)
axis square
xlabel('ssi control')
ylabel('ssi light')
title('spatial selectivity')
legend('X94','X98')

% bursts
bphz = condbc./condfr; % bursts per Hz
plot(squeeze(nanmean(bphz(prsv&~phe&~mu&x94,1,:),3)),squeeze(nanmean(bphz(prsv&~phe&~mu&x94,2,:),3)),'co','markerfacecolor','c')
hold on
plot(squeeze(nanmean(bphz(prsv&~phe&~mu&x98,1,:),3)),squeeze(nanmean(bphz(prsv&~phe&~mu&x98,2,:),3)),'mo','markerfacecolor','m')



figure
plot(squeeze(bphz(prsv&~phe&~mu,1,:)),depth(prsv&~phe&~mu),'ko','markerfacecolor','k')
hold on
plot(squeeze(bphz(pfsv&~phe&~mu,1,:)),depth(pfsv&~phe&~mu),'go','markerfacecolor','g')
axis ij
line([0,.15],[575,575],'color','k','linestyle',':')
line([0,.15],[375,375],'color','k','linestyle',':')
axis([0,0.15,100,1000])

figure
hold on
ls = '-';
% cond = l23rs&~phe& nlfr>1;
% errorbar(nanmean(rankedomi(cond,:)),nanstd(rankedomi(cond,:))./sqrt(length(find(cond))),'r','linestyle',ls)
cond = l4rs&~phe& nlfr>1;
errorbar(nanmean(rankedomi(cond,:)),nanstd(rankedomi(cond,:))./sqrt(length(find(cond))),'m','linestyle',ls)
cond = l5rs&~phe& nlfr>1;
errorbar(nanmean(rankedomi(cond,:)),nanstd(rankedomi(cond,:))./sqrt(length(find(cond))),'b','linestyle',ls)
legend('L4','L5')

%running OMI no touch - scatter
figure
plot(omicurve(prsv& mu == 0,1),depth(prsv& mu == 0),'ko')
hold on
plot(omicurve(pfsv& mu == 0,1),depth(pfsv& mu == 0),'go')
plot(omicurve(phe& mu == 0,1),depth(phe& mu == 0),'ro')
plot(omicurve(find(mu),1),depth(find(mu)),'co','markerfacecolor','c')
plot(omicurve(prsv&lightmod& mu == 0,1),depth(prsv&lightmod& mu == 0),'ko','markerfacecolor','k')
plot(omicurve(pfsv&lightmod& mu == 0,1),depth(pfsv&lightmod& mu == 0),'go','markerfacecolor','g')
plot(omicurve(phe&lightmod& mu == 0,1),depth(phe&lightmod& mu == 0),'ro','markerfacecolor','r')
axis ij
line([0,0],[100,900],'color','k')
line([-1,1],[575,575],'color','k','linestyle',':')
line([-1,1],[375,375],'color','k','linestyle',':')
legend('RS','FS','Chr2','MU')
xlabel('OMI')
ylabel('cortical depth')
title('OMI running - no touch')

cond = ~phe & ~mu & prsv & x94;
figure
plot(omicurve(cond,1),depth(cond),'ko','markerfacecolor','k','markersize',4)
% hold on
% plot(omicurve(cond&lightmod,1),depth(cond&lightmod),'ko','markerfacecolor','k')
axis ij
line([0,0],[0,1200],'color','k')
line([-1,1],[575,575],'color','k','linestyle',':')
line([-1,1],[375,375],'color','k','linestyle',':')
xlabel('OMI')
ylabel('cortical depth')
title('OMI running')

% running OMI no touch - running average
cond = ~phe & ~mu & nlfr>1 & prsv & x94;%   & touchmod;
figure
[x,y,xerr] = runningAverage(depth(cond),omicurve(cond,1),0,15);
plot(x,y,'k','linewidth',2)
hold on
plot(x+xerr,y,'k')
plot(x-xerr,y,'k')
axis ij
line([0,0],[100,900],'color','k')
line([-.8,.1],[575,575],'color','k','linestyle',':')
line([-.8,.1],[375,375],'color','k','linestyle',':')
axis([-.8,.1,200,950])
xlabel('average OMI')
ylabel('cortical depth um')

% running OMI mean touch scatter
figure
plot(omi(prsv& mu == 0),depth(prsv& mu == 0),'ko')
hold on
plot(omi(pfsv& mu == 0),depth(pfsv& mu == 0),'go')
plot(omi(phe& mu == 0),depth(phe& mu == 0),'ro')
plot(omi(find(mu)),depth(find(mu)),'co','markerfacecolor','c')
plot(omi(prsv&lightmod& mu == 0),depth(prsv&lightmod& mu == 0),'ko','markerfacecolor','k')
plot(omi(pfsv&lightmod& mu == 0),depth(pfsv&lightmod& mu == 0),'go','markerfacecolor','g')
plot(omi(phe&lightmod& mu == 0),depth(phe&lightmod& mu == 0),'ro','markerfacecolor','r')
axis ij
line([0,0],[100,900],'color','k')
line([-1,1],[575,575],'color','k','linestyle',':')
line([-1,1],[375,375],'color','k','linestyle',':')
legend('RS','FS','ChR2-cell','MU','RS sig','FS sig','ChR2 sig')
xlabel('OMI')
ylabel('cortical depth')
title('OMI running')

cond = ~phe & ~mu & prsv & x94 & touchmod;
figure
plot(omi(cond),depth(cond),'ko','markerfacecolor','k','markersize',4)
% hold on
% plot(omi(cond&lightmod),depth(cond&lightmod),'ko','markerfacecolor','k')
axis ij
line([0,0],[0,1200],'color','k')
line([-1,1],[575,575],'color','k','linestyle',':')
line([-1,1],[375,375],'color','k','linestyle',':')
xlabel('OMI')
ylabel('cortical depth')
title('OMI running')

% OMI histograms between X94 X98
figure
cond = cond4;
bins = [-1:.1:1];
[y,x] = histc(omi(cond&x98),bins);
plot(bins,y./sum(cond&x98),'m','linewidth',2)
hold on
[y,x] = histc(omi(cond&x94),bins);
plot(bins,y./sum(cond&x94),'c','linewidth',2)
legend('X98','X94')
xlabel('OMI')
ylabel('percent cells')
title('Layer 4: histogram of OMI')


figure
cond = ~phe & ~mu & nlfr>1 & prsv;
[x,y,xerr] = runningAverage(depth(cond),omi(cond),0,17);
plot(x,y,'k','linewidth',2)
hold on
plot(x+xerr,y,'k')
plot(x-xerr,y,'k')
axis ij
line([0,0],[100,900],'color','k')
line([-.5,.1],[575,575],'color','k','linestyle',':')
line([-.5,.1],[375,375],'color','k','linestyle',':')
axis([-.5,.1,200,950])
xlabel('average OMI')
ylabel('cortical depth um')

cond = prsv&~phe&mu == 0;
figure
[x,y,xerr] = runningAverage(depth(cond),omi(cond),0,15);
plot(x,y,'k','linewidth',2)
hold on
plot(x+xerr,y,'k')
plot(x-xerr,y,'k')
axis ij
line([0,0],[100,900],'color','k')
line([-.5,.1],[575,575],'color','k','linestyle',':')
line([-.5,.1],[375,375],'color','k','linestyle',':')
axis([-.5,.1,200,950])
xlabel('average OMI')
ylabel('cortical depth um')

% non-running OMI in depth
figure
plot(ominr(prs),depth(prs),'ko')
hold on
plot(ominr(pfs),depth(pfs),'go')
plot(ominr(phe),depth(phe),'ro')
plot(ominr(find(mu)),depth(find(mu)),'co')
plot(ominr(prsv&lightmod& mu == 0),depth(prsv&lightmod& mu == 0),'ko','markerfacecolor','k')
plot(ominr(pfsv&lightmod& mu == 0),depth(pfsv&lightmod& mu == 0),'go','markerfacecolor','g')
plot(ominr(phe&lightmod& mu == 0),depth(phe&lightmod& mu == 0),'ro','markerfacecolor','r')
plot(ominr(lightmod & mu == 1),depth(lightmod& mu == 1),'co','markerfacecolor','c')
axis ij
line([0,0],[100,900],'color','k')
line([-1,1],[550,550],'color','k','linestyle',':')
line([-1,1],[375,375],'color','k','linestyle',':')
legend('RS','FS','ChR2 cell')
xlabel('OMI')
ylabel('cortical depth')
title('OMI non-running')

% no running but touch
cond = ~phe & ~mu  & pfsv & x94;
figure
plot(ominr(cond),depth(cond),'ko','markerfacecolor','k','markersize',4)
% hold on
% plot(omi(cond&lightmod),depth(cond&lightmod),'ko','markerfacecolor','k')
axis ij
line([0,0],[100,900],'color','k')
line([-1,1],[575,575],'color','k','linestyle',':')
line([-1,1],[375,375],'color','k','linestyle',':')
xlabel('OMI')
ylabel('cortical depth')
title('OMI non-running')

% no running no touch
cond = ~phe & ~mu & prsv & x94;
figure
plot(omicurvenr(cond,1),depth(cond),'ko','markerfacecolor','k','markersize',4)
% hold on
% plot(omicurve(cond&lightmod,1),depth(cond&lightmod),'ko','markerfacecolor','k')
axis ij
line([0,0],[0,1200],'color','k')
line([-1,1],[575,575],'color','k','linestyle',':')
line([-1,1],[375,375],'color','k','linestyle',':')
xlabel('OMI')
ylabel('cortical depth')
title('OMI non-running')

figure
cond = ~phe & ~mu & nlfr>1 & prsv;
[x,y,xerr] = runningAverage(depth(cond),ominr(cond),0,15);
plot(x,y,'k','linewidth',2)
hold on
plot(x+xerr,y,'k')
plot(x-xerr,y,'k')
axis ij
line([0,0],[100,900],'color','k')
line([-.5,.1],[550,550],'color','k','linestyle',':')
line([-.5,.1],[375,375],'color','k','linestyle',':')
axis([-.6,.1,200,950])
xlabel('average OMI (non-runnning')
ylabel('cortical depth um')


stepsize = 50;
a = 100:stepsize:1000;
cond = ~phe & ~mu & prsv & x94;
clear momi; clear binns; clear momierr;
for i = 1:length(a)
    dinds = depth>=a(i)-stepsize/2 & depth<a(i)+stepsize/2;
    momi(i) = nanmean(omi(dinds&cond)); % meanomi for bin
    binns(i) = length(find(dinds&cond));
    momierr(i) = nanstd(omi(dinds&cond))./sqrt(length(find(dinds&cond)));
end

figure
barh(a,momi,'k')
hold on
herrorbar(momi,a,momierr,'k.')
axis ij
axis([-.5,.2,0,1050])
line([-.5,.2],[375,375],'linestyle',':','color','k')
line([-.5,.2],[575,575],'linestyle',':','color','k')
xlabel('OMI')
ylabel('depth (um)')

hold on
for i = 1:length(a)
    text(.15,a(i),int2str(binns(i)))
end


% spatial selectivity
figure
plot(ssi(prs,1),depth(prs),'ko','markerfacecolor','k')
hold on
plot(ssi(pfs,1),depth(pfs),'go','markerfacecolor','g')
axis ij
line([0,1],[550,550],'color','k','linestyle',':')
line([0,1],[375,375],'color','k','linestyle',':')
legend('RS','FS')
xlabel('SSI')
ylabel('cortical depth')
title('SSI running')

% spatial selectivity change
cond = ~phe & ~mu & nlfr>1 & prsv;
figure
plot(ssi(cond,2)-ssi(cond,1),depth(cond),'ko','markerfacecolor','k')
% hold on
% plot(ssi(pfs,2)-ssi(pfs,1),depth(pfs),'go','markerfacecolor','g')
axis ij
line([0,0],[100,900],'color','k')
line([-.6,.6],[550,550],'color','k','linestyle',':')
line([-.6,.6],[375,375],'color','k','linestyle',':')
axis([-.6,.6,100,900])
legend('RS','FS')
xlabel('SSI')
ylabel('cortical depth')
title('change in SSI with light')

cond = ~phe & ~mu & nlfr>1 & prsv;
figure
[x,y,xerr] = runningAverage(depth(cond),ssi(cond,2)-ssi(cond,1),0,15);
plot(x,y,'k','linewidth',2)
hold on
plot(x+xerr,y,'k')
plot(x-xerr,y,'k')
axis ij
line([0,0],[100,900],'color','k')
line([-.5,.1],[575,575],'color','k','linestyle',':')
line([-.5,.1],[375,375],'color','k','linestyle',':')
axis([-.5,.1,200,950])
xlabel('average SSI change')
ylabel('cortical depth um')

function ssi = get_ssi(curve)
ssi  = 1 - (((norm(curve)/max(curve)) - 1)./((sqrt(length(curve)))-1));