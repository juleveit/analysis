function size_populationana

% % SOM Arch
% animalids = {'140618', '140619', '140620'};
% blocks    = [3,         7,        4];
% lcol = 'r'; %lasercolor
% printpath = 'C:\Users\Julia\work\data\ArchSomPlots\';

% % Scnn Halo
% %SCNN Population
% animalids = {'140618', '140619', '140620', '140703', '140703', '140704', '140709', '140709', '140711', '140715', '140717', '140717', '140806', '140806', '141024', '141028', '141103', '141110', '141112', '141112', '141113', '141201', '141201b', '141204', '141209'};
% blocks    = [3,         7,        4,        8,        14,       5,        5,        8,        5,        8,        3,        14        4         8,        10,       8,        7,        3,        4,        9,        4,        2,        2,         3,        2];
% animal =    [1          2         3         4         4         5         6         6         7,        8         9,        9         10        10,       11,       12,       13,       14,       15,       15,       16,       17,       18,        19,       20];
% electrodes= [[1,32];   [1,32];    [1,32];   [1,32];   [1,32];   [1,32];   [1,32];   [1,32];   [1,32];   [1,32];   [1,32];   [1,32];   [1,32];   [1,32];   [1,32];   [1,32];   [1,16];   [1,16];   [1,16];   [1,16];   [1,16];   [1,16];   [1,16];    [1,16];   [1,16];];
% penangle =  [10,        10,       10,       10,       10,       10,       10,       10,       10,       10,       10,       10,       10,       10,       10,       10,       10,       10,       10,       10,       10,       10,       10,        10,       10];
% printpath = 'C:\Users\Julia\work\data\populations\Scnn_Halo_withnew\size\units\';
% runprintpath = 'C:\Users\Julia\work\data\populations\Scnn_Halo_withnew\size\running\';
% l5RSprintpath = 'C:\Users\Julia\work\data\populations\Scnn_Halo_withnew\size\l5rs\'
% popfile = 'C:\Users\Julia\work\data\populations\Scnn_Halo_withnew\size\size_population_withnewanimals.mat';

% % Scnn Halo
% %SCNN Population scotts paper
% animalids = {'140618', '140619', '140620', '140703', '140703', '140704', '140709', '140709', '140711', '140715', '140717', '140717', '140806', '140806'};
% blocks    = [3,         7,        4,        8,        14,       5,        5,        8,        5,        8,         3,        14        4         8];
% animal =    [1          2         3         4         4         5         6         6         7,        8          9,        9         10        10];
% electrodes= [32,        32,       32,       32,       32,       32,       32,       32,       32,       32,        32,       32,       32,       32];
% penangle =  [10,        10,       10,       10,       10,       10,       10,       10,       10,       10,        10,       10,       10,       10];
% printpath = 'C:\Users\Julia\work\data\populations\size\units\';
% runprintpath = 'C:\Users\Julia\work\data\populations\size\running\';
% l5RSprintpath = 'C:\Users\Julia\work\data\populations\size\l5rs\'
% popfile = 'C:\Users\Julia\work\data\populations\size\size_population.mat';

% %SCNN-CHR2 Population 
% animalids = {'150323', '150325'};
% blocks    = [6,         4];
% animal =    [1          2];
% penangle = 10;
% printpath = 'C:\Users\Julia\work\data\populations\scnn_chr2\size\units\';
% runprintpath = 'C:\Users\Julia\work\data\populations\scnn_chr2\size\running\';
% condprintpath = 'C:\Users\Julia\work\data\populations\scnn_chr2\size\conditions\';
% l5RSprintpath = 'C:\Users\Julia\work\data\populations\scnn_chr2\size\l5rs\';
% popfile = 'C:\Users\Julia\work\data\populations\scnn_chr2\size_population.mat';

% %SCNN-eArch Population 
% animalids = {'160219', '160222', '160223'};
% blocks    = [ 4,        2,        7];
% animal =    [ 1         2,        3];
% electrodes =[[1,16];   [1,16];   [1,16]];
% penangle =  [ 25,       25,       25];
% printpath = 'C:\Users\Julia\work\data\populations\scnn_eArch\size\units\';
% runprintpath = 'C:\Users\Julia\work\data\populations\scnn_eArch\size\running\';
% condprintpath = 'C:\Users\Julia\work\data\populations\scnn_eArch\size\conditions\';
% l5RSprintpath = 'C:\Users\Julia\work\data\populations\scnn_eArch\size\l5rs\';
% popfile = 'C:\Users\Julia\work\data\populations\scnn_eArch\size_population.mat';

% % SOM population old
% animalids = {'140807', '140807', '140812', '140812', '140815', '140815', '140818'};
% blocks    = [5,         8,        3,        5,        4,        6,        5];
% animal =    [1          1         2         2         3         3         4];
% printpath = 'C:\Users\Julia\work\data\populations\SOM_Halo\size\units\';
% runprintpath = 'C:\Users\Julia\work\data\populations\SOM_Halo\size\running\';
% popfile = 'C:\Users\Julia\work\data\populations\SOM_Halo\size\contrast_population.mat';
% % popfile = 'C:\Users\Julia\work\data\populations\SOM_Halo\size\contrast_population_onlymale.mat';

% % SOM P0 population
% animalids = {'140807', '140807', '140812', '140812', '140815', '140815', '140818', '141014', '141014', '141015', '141106', '141110b', '141111', '150123', '150213','150523'};
% blocks    = [5,         8,        3,        5,        4,        6,        5,        3,        10,       5,        3,        7,         4,        6,        3,       3];
% animal =    [1          1         2         2         3         3         4,        5,        5,        6,        7,        8,         9,        10,       11,      12];
% electrodes =[32,        32,       32,       32,       32,       32,       32,       32,       32,       32,       32,       32,        32,       32,       16,      32];
% penangle =  [10,        10,       10,       10,       10,       10,       10,       10,       10,       10,       10,       10,        10,       10,       20,      25];
% printpath = 'C:\Users\Julia\work\data\populations\SOM_Halo\size\units\';
% runprintpath = 'C:\Users\Julia\work\data\populations\SOM_Halo\size\running\';
% popfile = 'C:\Users\Julia\work\data\populations\SOM_Halo\size\size_population.mat';
% % popfile = 'C:\Users\Julia\work\data\populations\SOM_Halo\size\contrast_population_onlymale.mat';

% % undecided about 150407 and 150406 (P1 morning injections, not much light effect)
% % % 
% % % % SOM later population
% animalids = {'150331', '150401','150527','150529','150602','150603','150625','150825','150831','150902','150907','150909', '150915', '150916', '151023', '151027','151109', '151110', '151209', '160122'};
% blocks    = [3,         5,       11,      4,       5,       3,       6,       5,       4,       3,       3,       4,        3,        2,       14,       3,        11,       13,       6,        3];
% animal    = [1,         2,       3,       4,       5,       6,       7,       8,       9,       10,      11,      12,       13,       14,      15,       16,       17,       18,       19,       20];
% electrodes =[[1,32];    [1,32];  [1,32];  [1,32];  [1,32];  [1,16];  [1,16];  [17,32]; [1,16];  [1,16];  [1,16];  [1,16];   [17,32];  [1,16];  [17,32];  [17,32]; [17,32];  [1,16];   [17,32];  [1,16]];
% penangle =  [25,        25,      25,      25,      25,      10,      10,      25,      25,      25,      25,      25,       25,       25,      25,       25,       25,       25,       25,       25];
% % age       [P6,        P6,      P2?(P0), P2?(P0), P2?(P0), P1
% printpath = 'C:\Users\Julia\work\data\populations\SOM_Halo_later\size\units\';
% runprintpath = 'C:\Users\Julia\work\data\populations\SOM_Halo_later\size\running\';
% oscillprintpath = 'C:\Users\Julia\work\data\populations\SOM_Halo_later\size\withlfpoverview\';
% popfile = 'C:\Users\Julia\work\data\populations\SOM_Halo_later\size\size_population.mat';

% % SOM Halo+ChR2 % 150914 not all sizes - up to 9 all Halo then CHr2
% animalids = {'150331', '150401','150602','150625','150909', '151023', '151209', '160122', '151229', '160311', '160405_2'};
% blocks    = [3,         5,       5,       6,       4,        14,       6,        3,        10,       7,        3];
% animal    = [1,         2,       3,       4,       5,        6,        7,        8,        9,        10,       11];
% electrodes =[[1,32];    [1,32];  [1,32]; [1,16];  [1,16];   [17,32];  [17,32];  [1,16];   [1,16];   [1,16];   [17,32]];
% penangle =  [25,        25,      25,      10,      25,       25,       25,       25,       25,       25,       25];
% printpath = 'C:\Users\Julia\work\data\populations\SOM_HaloandChR2\size\units\';
% runprintpath = 'C:\Users\Julia\work\data\populations\SOM_HaloandChR2\size\running\';
% oscillprintpath = 'C:\Users\Julia\work\data\populations\SOM_HaloandChR2\size\withlfpoverview\';
% popfile = 'C:\Users\Julia\work\data\populations\SOM_HaloandChR2\size\size_population.mat';

% % PV-earch population
% animalids = {'151117', '151211', '160114', '160115'};
% blocks    = [2,         3,        4,        3];
% animal    = [1,         2,        3,        4];
% electrodes =[[1,16];    [1,16];  [1,16];   [1,16]];
% penangle =  [25,        25,       25,       25];
% printpath = 'C:\Users\Julia\work\data\populations\PV_eArch\size\units\';
% runprintpath = 'C:\Users\Julia\work\data\populations\PV_eArch\size\running\';
% oscillprintpath = 'C:\Users\Julia\work\data\populations\PV_eArch\size\withlfpoverview\';
% popfile = 'C:\Users\Julia\work\data\populations\PV_eArch\size\size_population.mat';

% % control SOM only population
% animalids = {'151214', '151215', '151217', '160112', '160113'};
% blocks    = [2,         7,        3,        3,        6];
% animal    = [1,         2,        3,        4,        5];
% electrodes =[[1,16];   [1,16];   [1,16];   [1,16];   [1,16]];
% penangle =  [25,        25,       25,       25,       25];
% printpath = 'C:\Users\Julia\work\data\populations\control\size\units\';
% runprintpath = 'C:\Users\Julia\work\data\populations\control\size\running\';
% oscillprintpath = 'C:\Users\Julia\work\data\populations\control\size\withlfpoverview\';
% popfile = 'C:\Users\Julia\work\data\populations\control\size\size_population.mat';
% 
% % anesthetized SOM only population
% animalids = {'150826', '151230', '151231'};
% blocks    = [3,         3,        11];
% animal    = [1,         2,        3];
% electrodes =[[1,16];   [1,16];   [1,16]];
% penangle =  [25,        25,       25];
% printpath = 'C:\Users\Julia\work\data\populations\SOM_Halo_anesth\size\units\';
% runprintpath = 'C:\Users\Julia\work\data\populations\SOM_Halo_anesth\size\running\';
% oscillprintpath = 'C:\Users\Julia\work\data\populations\SOM_Halo_anesth\size\withlfpoverview\';
% popfile = 'C:\Users\Julia\work\data\populations\SOM_Halo_anesth\size\size_population.mat';

% % anesthetized SOM+PV combined only population
% animalids = {'150826', '151230', '151231', '150821', '150824'};
% blocks    = [3,         3,        11,       3,        18];
% animal    = [1,         2,        3,        4,        5];
% electrodes =[[1,16];   [1,16];   [1,16];   [1,16];   [17,32]];
% penangle =  [25,        25,       25,       25,       25];
% printpath = 'C:\Users\Julia\work\data\populations\SOMPVcombined_anesth\size\units\';
% runprintpath = 'C:\Users\Julia\work\data\populations\SOMPVcombined_anesth\size\running\';
% oscillprintpath = 'C:\Users\Julia\work\data\populations\SOMPVcombined_anesth\size\withlfpoverview\';
% popfile = 'C:\Users\Julia\work\data\populations\SOMPVcombined_anesth\size\size_population.mat';

% % PV Halo population
% animalids = {'150629', '150730', '150731', '150804','150818','150820','150823','150824', '151104'};
% blocks    = [4,         4,        3,         4,      3,       3,       5,       3,        5];
% animal    = [1,         2,        3,         4,      5,       6,       7,       8,        9];
% electrodes =[[1,16];    [1,16];   [1,16];   [1,16];  [1,16];  [1,16];  [1,16];  [17,32];  [17,32]];
% penangle =  [10,        10,       25,        25,     25,      25,      25,      25,       25];
% printpath = 'C:\Users\Julia\work\data\populations\PV_Halo\size\units\';
% runprintpath = 'C:\Users\Julia\work\data\populations\PV_Halo\size\running\';
% oscillprintpath = 'C:\Users\Julia\work\data\populations\PV_Halo\size\withlfpoverview\';
% popfile = 'C:\Users\Julia\work\data\populations\PV_Halo\size\size_population.mat';
% 
% HUGEass combined population of all animals - 1:20 SOM, 21:32 PV Halo 33:39 PV eArch
animalids = {'150331', '150401','150527','150529','150602','150603','150625','150825','150831','150902','150907','150909', '150915', '150916', '151023', '151027','151109', '151110', '151209', '160122', '150629', '150730', '150731', '150804','150818','150820','150823','150824', '151104', '160217', '160324', '160328', '151117', '151211', '160114', '160115', '160204', '160205', '160210'};
blocks    = [3,         5,       11,      4,       5,       3,       6,       5,       4,       3,       3,       4,        3,        2,       14,       3,        11,       13,       6,        3,        4,        4,        3,        4,       3,       3,       5,       3,        5,        5,        2,        2,        2,        3,        4,        3,        3,        4,        9];
animal    = [1,         2,       3,       4,       5,       6,       7,       8,       9,       10,      11,      12,       13,       14,      15,       16,       17,       18,       19,       20,       21,       22,       23,       24,      25,      26,      27,      28,       29,       30,       31,       32,       33,       34,       35,       36,       37,       38,       39];
electrodes =[[1,32];    [1,32];  [1,32];  [1,32];  [1,32];  [1,16];  [1,16];  [17,32]; [1,16];  [1,16];  [1,16];  [1,16];   [17,32];  [1,16];  [17,32];  [17,32];  [17,32]; [1,16];   [17,32];  [1,16];   [1,16];   [1,16];   [1,16];   [1,16];  [1,16];  [1,16];  [1,16];  [17,32];  [17,32];  [1, 16];  [1, 16];  [1, 16];  [1,16];    [1,16];  [1,16];   [1,16];   [1,16];   [17,32];  [17,32]];
penangle =  [25,        25,      25,      25,      25,      10,      10,      25,      25,      25,      25,      25,       25,       25,      25,       25,       25,       25,       25,       25,       25,       25,       25,       25,      25,      25,      25,      25,       25,       25,       25,       25,       25,        25,      25,       25,       25,       25,       25];
% age       [P6,        P6,      P2?(P0), P2?(P0), P2?(P0), P1
printpath = 'C:\Users\Julia\work\data\populations\SOMPVcombined\size\units\';
runprintpath = 'C:\Users\Julia\work\data\populations\SOMPVcombined\size\running\';
oscillprintpath = 'C:\Users\Julia\work\data\populations\SOMPVcombined\size\withlfpoverview\';
popfile = 'C:\Users\Julia\work\data\populations\SOMPVcombined\size\size_population.mat';

% % PV Halo + PV eArch population, up to animal 12 all Halo then eArch
% animalids = {'150629', '150730', '150731', '150804','150818','150820','150823','150824', '151104', '160217', '160324', '160328', '151117', '151211', '160114', '160115', '160204', '160205', '160210'};
% blocks    = [4,         4,        3,         4,      3,       3,       5,       3,        5,        5,        2,        2,        2,         3,        4,        3,       3,        4,        9];
% animal    = [1,         2,        3,         4,      5,       6,       7,       8,        9,        10,       11,       12,       13,        14,       15,       16,      17,       18,       19];
% electrodes =[[1,16];    [1,16];   [1,16];   [1,16];  [1,16];  [1,16];  [1,16];  [17,32];  [17,32]; [1, 16];  [1, 16];  [1, 16];  [1,16];    [1,16];  [1,16];   [1,16];   [1,16];   [17,32];  [17,32]];
% penangle =  [10,        10,       25,        25,     25,      25,      25,      25,       25,       25,       25,       25,       25,        25,       25,       25,      25,       25,       25];
% printpath = 'C:\Users\Julia\work\data\populations\PV_HaloeArch\size\units\';
% runprintpath = 'C:\Users\Julia\work\data\populations\PV_HaloeArch\size\running\';
% oscillprintpath = 'C:\Users\Julia\work\data\populations\PV_HaloeArch\size\withlfpoverview\';
% popfile = 'C:\Users\Julia\work\data\populations\PV_HaloeArch\size\size_population.mat';

tic
lcol = 'r'; %lasercolor

recalculate = 0;
printyn = 1;

% chronux parameters
params.tapers = [2,5]; params.Fs = 1000; params.err = [2, 0.05]; params.trialave = 1;

%Scott Kernel
%%%%%%%%%%%%kernal propeties%%%%%%%%%%%%%%%
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
%%%%%%%%%%kernal properties%%%%%%%%%%%%%%%

onlymod = 0; % run only for visually modulated
sr = 1000;

if ~exist(popfile) || recalculate

    cll = 1;
    for blck = 1:length(blocks)

        supath = ['C:\Users\Julia\work\data\' animalids{blck} '\singleunits\'];
        basename = [animalids{blck} '_block' int2str(blocks(blck)) '_tet'];

        files = dir([supath, basename, '*.mat']);

        prestim = 300;
        poststim = 700;
        respwin = 501:1500; % after stimulus onset
        offsetwin = 1501:2500;
        respwin = respwin+prestim;
        offsetwin = offsetwin+prestim;
        freqbinwidth = 5;
        
        clear filtmat; clear powmat; clear phasmat;

        for fi = 1:length(files)

            if strfind(files(fi).name, 'MU')
                continue;
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
            
            wvchan = find(var(result.waveforms) == max(var(result.waveforms)));
            cm = [3,4,1,2]; % confusion matrix? take lfp from two electrodes away to not get too many spike related phase resets
            lfp = result.lfp(:,cm(wvchan))';
            
            
            trialdur = result.stimduration*1000;
            msstamps = result.msstamps;            
            
            if length(msstamps)~=length(result.light)
%                 disp('');
%             msstamps(16) = []; % for 150414 block 10
%             result.msstamps = msstamps;
%             save([supath, SUfiles(cl).name],'result');            
                pause;
            end
            
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
            if ~isempty(oktrials)
                cellmeanspeed(cll) = mean(meanspeed(oktrials));
                cellminspeed(cll) = min(meanspeed(oktrials));
            else
                cellmeanspeed(cll) = NaN;
                cellminspeed(cll) = NaN;
            end

            
            %find gamma peaks for this animal
            beta = [15,40];
            gamma = [50,70];
            large = find(result.gratingInfo.size == max(unique(result.gratingInfo.size)) & result.light == 0);
            small = find(result.gratingInfo.size == min(unique(result.gratingInfo.size(result.gratingInfo.size~=0))) & result.light == 0);
            bsl = find(result.gratingInfo.size == 0 & result.light == 0);
            lr1 = intersect(large,oktrials);
            sr1 = intersect(small,oktrials);
            br1 = intersect(bsl,oktrials);
            for i = 1:length(large)
                [pl(i,:),f] = pmtm(lfp(result.msstamps(large(i))+700:result.msstamps(large(i))+1500),3,[],1000);
                [ps(i,:),f] = pmtm(lfp(result.msstamps(small(i))+700:result.msstamps(small(i))+1500),3,[],1000);              
            end     
            plr1 = nan(1,size(pl,2)); psr1 = nan(1,size(pl,2));
            if length(lr1>=5)
                for i = 1:length(lr1)
                    [plr1(i,:),f] = mtspectrumc(lfp(result.msstamps(lr1(i))+700:result.msstamps(lr1(i))+1500),params);
%                     [plr1(i,:),f] = pmtm(lfp(result.msstamps(lr1(i))+300:result.msstamps(lr1(i))+1300),3,[],1000);
                end
            end
            if length(sr1>=5)
                for i = 1:length(sr1)
                    [psr1(i,:),f] = mtspectrumc(lfp(result.msstamps(sr1(i))+700:result.msstamps(sr1(i))+1500),params); 
%                     y(i,:) = fft(lfp(result.msstamps(sr1(i)):result.msstamps(sr1(i))+2000),2048);
%                     [psr1(i,:),f] = pmtm(lfp(result.msstamps(sr1(i))+300:result.msstamps(sr1(i))+1300),3,[],1000);
                end
            end
            if length(br1>=5)
                for i = 1:length(br1)
                    [pblr1(i,:),f] = mtspectrumc(lfp(result.msstamps(br1(i))+700:result.msstamps(br1(i))+1500),params);
                end
            else
                pblr1 = nan(1,size(pl,2));
            end
            sdbl = std(pblr1,1,1);
            ff = 500*linspace(0,1,1025);
            
            b1 = find(f>beta(1),1); b2 = find(f>beta(2),1);
            g1 = find(f>gamma(1),1); g2 = find(f>gamma(2),1);
            
            if ~isnan(plr1(1))
                bsig = nanmean(plr1(:,b1:b2));
            else
                bsig = nanmean(pl(:,b1:b2));
            end
            if ~isnan(psr1(1))
                gsig = nanmean(psr1(:,g1:g2));
            else
                gsig = nanmean(ps(:,g1:g2));
            end
            
            if isempty(find(diff(bsig)>0)) % there is no clear beta peak
                bpi = round((b1+b2)/2);
            else
                peaks = find(diff(bsig)>0)+1;
                pvs = bsig(peaks);
                bpi = peaks(pvs == max(pvs));
                bpi = bpi+b1-1;
            end
            if isempty(find(diff(gsig)>0)) % there is no clear beta peak
                gpi = round((g1+g2)/2);
            else
                peaks = find(diff(gsig)>0)+1;
                pvs = gsig(peaks);
                gpi = peaks(pvs == max(pvs));
                gpi = gpi+g1-1;
            end
            
            gbandwidth = 20; g3b = find(f>70,1); g3e = find(f>100,1);
            gamma1 = eegfilt(lfp,sr,f(bpi)-gbandwidth/2,f(bpi)+gbandwidth/2);
            gamma2 = eegfilt(lfp,sr,f(gpi)-gbandwidth/2,f(gpi)+gbandwidth/2);
%             gamma3 = eegfilt(lfp,sr,f(g3b),f(g3e));
            h1 = hilbert(gamma1); gpow1 = abs(h1); gphas1 = angle(h1);
            h2 = hilbert(gamma2); gpow2 = abs(h2); gphas2 = angle(h2);
%             h3 = hilbert(gamma3); gpow3 = abs(h3); gphas3 = angle(h3);
            
            lgi(cll) = bpi; hgi(cll) = gpi;
            
            lgpeakratio = nanmean(nanmean(plr1(:,bpi-2:bpi+2),2))./...
                nanmean([nanmean(nanmean(plr1(:,bpi-12:bpi-8),2)),nanmean(nanmean(plr1(:,bpi+8:bpi+12),2))]);
            
            for i = 1:100/freqbinwidth
                filtmat(i,:) = eegfilt(lfp,sr,(i-1)*freqbinwidth+1,i*freqbinwidth);
                h = hilbert(filtmat(i,:));
                powmat(i,:) = abs(h); phasmat(i,:) = angle(h);
            end
                            
            trialnfft = 2^nextpow2(800);
            fftfax = sr/2*linspace(0,1,trialnfft/2+1);
            
            gamma1phases = zeros(length(msstamps),trialdur+poststim+prestim);
            gamma2phases = zeros(length(msstamps),trialdur+poststim+prestim);
            vst = zeros(8,length(msstamps),trialdur+poststim+prestim);
            clear lfpspect; clear fax; clear lfpoffsspect;
            for i = 1:length(msstamps)
                resp(i,:) = chan(msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim);  
                hh = find(resp(i,1001:1800))'; ptresp(i).times = hh./1000;                
                burstresp(i,:) = burstchan(msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim); 
                lfpresp(i,:) = lfp(msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim);
                gamma1resp(i,:) = gpow1(msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim); 
                gamma2resp(i,:) = gpow2(msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim); 
                gamma1phasresp(i,:) = gphas1(msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim); 
                gamma2phasresp(i,:) = gphas2(msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim); 
                gamma1phases(i,find(resp(i,:))) = gamma1phasresp(i,find(resp(i,:)));
                gamma2phases(i,find(resp(i,:))) = gamma2phasresp(i,find(resp(i,:)));
                
                allS(i,:)=mtspectrumc(squeeze(lfpresp(i,1001:1800))',params);
                allSz(i,:) = (allS(i,:)-mean(pblr1))./sdbl;
                
                y = fft(lfpresp(i,1001:1800),trialnfft);
                fftspect(i,:) = abs(y(1:trialnfft/2+1));
                
                hhh = smooth(squeeze(allS(i,24:45)),7);             
                if isempty(find(diff(hhh)>0)) % there is no clear beta peak
                    g1peakind(i) = NaN;
                else
                    peaks = find(diff(hhh)>0)+1;
                    hvs = hhh(peaks);
                    pind = peaks(hvs == max(hvs));
                    g1peakind(i) = pind+24-1;
                end
                
                ppi = pi;
                for pb = 1:8 % simplest binning, in fixed bins around the circle
                    vstfrg1(pb,i) = length(find(gamma1phases(i,respwin)> -pi+(pb-1)*pi/4 & gamma1phases(i,respwin)< -pi+pb*pi/4));
                    vstfrg2(pb,i) = length(find(gamma2phases(i,respwin)> -pi+(pb-1)*pi/4 & gamma2phases(i,respwin)< -pi+pb*pi/4));
                end
                [lfpspect(i,:),fax] = pmtm(lfpresp(i,respwin),3,[],1000); 
                lfpoffsspect(i,:) = pmtm(lfpresp(i,offsetwin),3,[],1000);
                spks = find(resp(i,respwin));
                lfpmat = zeros(length(spks),301);
                if ~isempty(spks)
                    for ss = 1:length(spks)
                        lfpmat(ss,:) = lfpresp(i,(spks(ss)+respwin(1)-150):(spks(ss)+respwin(1)+150));
                    end
                    stalfp(i,:) = nanmean(lfpmat,1);
                    stan(i) = size(lfpmat,1);
                else
                    stalfp(i,:) = nan(1,301);
                    stan(i) = 0;
                end
                
                for j = 1:size(phasmat,1)
                    allphaseresp(j,i,:) = phasmat(j, msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim);
                    allpowresp(j,i,:) = powmat(j, msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim);
                end
                
                isiresp{i} = diff(find(resp(i,respwin)));
                msl(i) = mean(find(resp(i,respwin)));

            end
            
            lfpspect = lfpspect(:,1:150); fax = fax(1:150);
            lfpoffsspect = lfpoffsspect(:,1:150);             
            
            msta = linspace(-prestim,trialdur+poststim,size(resp,2));

            lightresp = resp(find(result.light),:);
            nolightresp = resp(find(result.light == 0),:);
            lightlfpresp = lfpresp(find(result.light),:);
            nolightlfpresp = lfpresp(find(~result.light),:);

            frs = sum(resp(:,respwin),2)./(length(respwin)/1000);
            brs = sum(burstresp(:,respwin),2);
            bl = sum(resp(:,1:prestim),2)./(prestim/1000);
            sc = sum(resp(:,respwin),2);
            gpower1 = mean(gamma1resp(:,respwin),2);
            gpower2 = mean(gamma2resp(:,respwin),2);
            gps = lfpspect(:,lgi(cll));
            hgps = lfpspect(:,hgi(cll));
            g3ps = nanmean(lfpspect(:,g3b:g3e),2);
            allpows = mean(allpowresp(:,:,respwin),3);
            l0s = find(result.light == 0);
            
            crudegpow(cll,:) = gps(l0s); crudehgpow(cll,:) = hgps(l0s); crudeg3pow(cll,:) = g3ps(l0s);
            firingrates(cll,:) = frs(l0s);
            allpowers(cll,:,:) = allpows(:,l0s);
            [gammafrcorr(cll),gammafrp(cll)] = nancorr(gps(l0s),frs(l0s));
            [highgammafrcorr(cll),highgammafrp(cll)] = nancorr(hgps(l0s),frs(l0s));
            [gamma3frcorr(cll),gamma3frp(cll)] = nancorr(g3ps(l0s),frs(l0s));
            for i = 1:size(allpows,1)
                [allfrcorr(cll,i), allfrp(cll,i)] = nancorr(allpows(i,l0s),frs(l0s));
                sizes = unique(result.gratingInfo.size);
                for s = 1:length(sizes) % single size condition
                    [allscfrcorr(cll,s,i),allscfrp(cll,s,i)] = nancorr(allpows(i,result.light == 0 & result.gratingInfo.size == sizes(s)),frs(result.light == 0 & result.gratingInfo.size == sizes(s)));
                end
            end

%              %determine if cell is visually modulated
%             blfr = sum(resp(:,1:prestim),2);
%             vrfr = sum(resp(:,prestim+40:2*prestim+40-1),2);
%             vm(fi) = ttest2(blfr,vrfr);
%             if isnan(vm(fi)), vm(fi) = 0; end
            
            tmpi = zeros(size(allphaseresp));
            for i = 1:size(allphaseresp,1)
                tmp = zeros(size(resp));
                tmp(find(resp)) = allphaseresp(i,find(resp));
                tmpi(i,:,:) = tmp;
            end
            g1 = zeros(size(gamma1phasresp));
            g1(find(resp)) = gamma1phasresp(find(resp));
            g2 = zeros(size(gamma2phasresp));
            g2(find(resp)) = gamma2phasresp(find(resp));
            
            oris = unique(result.gratingInfo.Orientation);  oris(find(oris == -1)) = [];
            sizes = unique(result.gratingInfo.size);  sizes(find(sizes == 0)) = []; 
            en_respg1 = zeros(8,size(resp,1),length(respwin));
            en_respg2 = zeros(8,size(resp,1),length(respwin));
            es_respg1 = zeros(8,size(resp,1),length(respwin));
            es_respg2 = zeros(8,size(resp,1),length(respwin));
            shownsizes(cll,:) = sizes;
            
                        
            for l = 1:length(unique(result.light))
                nocontinds = find(result.gratingInfo.size == 0 & result.light == l-1);
                spontspect = nanmean(allS(nocontinds,:));
                for sz = 1:length(sizes)
                    thisinds = find(result.gratingInfo.size == sizes(sz) & result.light == l-1);
                    
%                     [C(cll,l,sz,:),phi(cll,l,sz,:),S12(cll,l,sz,:),S1(cll,l,sz,:),...
%                         S2(cll,l,sz,:),chfx,zerosp(cll,l,sz,:),confC(cll,l,sz),...
%                         phistd(cll,l,sz,:),Cerr(cll,l,sz,:,:)] = coherencycpt(squeeze(lfpresp(thisinds,1001:1800))',...
%                         ptresp(thisinds),params);
                     if sz == length(sizes)
                         largegind(cll,l,:) = g1peakind(thisinds);
                     end
                     
                     szfrs(l,sz) = mean(frs(thisinds));
                                   
                    [S,chf,Serr]=mtspectrumc(squeeze(lfpresp(thisinds,1001:1800))',params);
                    condS(cll,l,sz,:) = S(1:150); condSerr(cll,l,sz,:,:) = Serr(:,1:150);
                    
                    condfftspect(cll,l,sz,:) = nanmean(fftspect(thisinds,:));
                    condfftspecterr(cll,l,sz,:) = nanstd(fftspect(thisinds,:))./sqrt(length(thisinds));
                    
                    for i = 1:length(thisinds)
                        hlp(i,:) = allS(thisinds(i),:)./spontspect;
                    end
                    condindratio(cll,l,sz,:) = nanmean(hlp,1);
                    condindratioerr(cll,l,sz,:) = nanstd(hlp,1,1);
                    
%                     [Call(cll,l,sz,:,:)] = coherencycpt(squeeze(lfpresp(thisinds,1001:1800))',...
%                         ptresp(thisinds),noavgparams);
                    
                    % running
                    thisruninds = intersect(thisinds,oktrials); thisstillinds = intersect(thisinds,stilltrials);
                    spontruninds = intersect(nocontinds,oktrials); spontstillinds = intersect(nocontinds,stilltrials);
                    if ~isempty(thisruninds)
                        thisrunn(cll,l,sz) = length(thisruninds);
                        [runS,chf,runSerr] = mtspectrumc(squeeze(lfpresp(thisruninds,1001:1800))',params);
                        condrunS(cll,l,sz,:) = runS(1:150); condrunSerr(cll,l,sz,:,:) = runSerr(:,1:150);
                        r1condallS(cll,l,sz,:) = nanmean(allS(thisruninds,:),1);
                        r1condallSerr(cll,l,sz,:) = nanstd(allS(thisruninds,:),1,1)./sqrt(length(thisruninds));
                        r1condallSz(cll,l,sz,:) = nanmean(allSz(thisruninds,:),1);
                        r1condallSzerr(cll,l,sz,:) = nanstd(allSz(thisruninds,:),1,1);
                        r1condfftspect(cll,l,sz,:) = nanmean(fftspect(thisruninds,:),1);
                        r1condfftspecterr(cll,l,sz,:) = nanstd(fftspect(thisruninds,:),1,1)./sqrt(length(thisruninds));
                        r1evokedspect(cll,l,sz,:) = nanmean(allS(thisruninds,1:103))./nanmean(allS(spontruninds,1:103));
                        if ~isempty(spontruninds)
                            % bootstrap evoked spectrum confidence intervals
                            for b = 1:10000 % 1000 bootstraps
                                rs = randperm(length(spontruninds)); % random spontaneous trial
                                rv = randperm(length(thisruninds)); % random visual trial
                                ratio(b,:) = allS(thisruninds(rv(1)),:)./allS(spontruninds(rs(1)),:);
                            end
                            for frq = 1:103 % sort all the frequencies
                                sorted = sort(ratio(:,frq));
                                r1evokedconf(cll,l,sz,frq,1) = sorted(500);   % bootstrapped 95% confidence interval
                                r1evokedconf(cll,l,sz,frq,2) = sorted(9500);
                                a = allS(thisruninds,frq);
                                b = allS(spontruninds,frq);
                                varratio(frq) = (nanmean(a)^2*nanvar(b)+nanmean(b)^2*nanvar(a))/nanmean(b)^4;
                                r1evokedspecterr(cll,l,sz,frq) = sqrt(varratio(frq))/sqrt(length(spontruninds)); % the math way
                            end   
                        else
                            r1evokedspecterr(cll,l,sz,:) = nan(1,103);
                            r1evokedconf(cll,l,sz,:,:) = nan(103,2);
                        end
                    else
                        thisrunn(cll,l,sz) = 0;
                        condrunS(cll,l,sz,:) = nan(1,150); condrunSerr(cll,l,sz,:,:) = nan(2,150);
                        r1condallS(cll,l,sz,:) = nan(1,size(allS,2));
                        r1condallSerr(cll,l,sz,:) = nan(1,size(allS,2));
                        r1condallSz(cll,l,sz,:) = nan(1,size(allS,2));
                        r1condallSzerr(cll,l,sz,:) = nan(1,size(allS,2));
                        r1condfftspect(cll,l,sz,:) = nan(1,size(allS,2));
                        r1condfftspecterr(cll,l,sz,:) = nan(1,size(allS,2));
                        r1evokedspect(cll,l,sz,:) = nan(1,103);
                        r1evokedspecterr(cll,l,sz,:) = nan(1,103);
                        r1evokedconf(cll,l,sz,:,:) = nan(103,2);
                    end
                    if ~isempty(thisstillinds)
                        thisstilln(cll,l,sz) = length(thisstillinds);
                        [stillS,chf,stillSerr] = mtspectrumc(squeeze(lfpresp(thisstillinds,1001:1800))',params);
                        condstillS(cll,l,sz,:) = stillS(1:150); condstillSerr(cll,l,sz,:,:) = stillSerr(:,1:150);
                        r0condallS(cll,l,sz,:) = nanmean(allS(thisstillinds,:),1);
                        r0condallSerr(cll,l,sz,:) = nanstd(allS(thisstillinds,:),1,1)./sqrt(length(thisstillinds));
                        r0condfftspect(cll,l,sz,:) = nanmean(fftspect(thisstillinds,:),1);
                        r0condfftspecterr(cll,l,sz,:) = nanstd(fftspect(thisstillinds,:),1,1)./sqrt(length(thisstillinds));
                        r0evokedspect(cll,l,sz,:) = nanmean(allS(thisstillinds,1:103))./nanmean(allS(spontstillinds,1:103));
                        if ~isempty(spontstillinds)
                            % bootstrap evoked spectrum confidence intervals
                            for b = 1:10000 % 1000 bootstraps
                                rs = randperm(length(spontstillinds)); % random spontaneous trial
                                rv = randperm(length(thisstillinds)); % random visual trial
                                ratio(b,:) = allS(thisstillinds(rv(1)),:)./allS(spontstillinds(rs(1)),:);
                            end
                            for frq = 1:103 % sort all the frequencies
                                sorted = sort(ratio(:,frq));
                                r0evokedconf(cll,l,sz,frq,1) = sorted(500);   % bootstrapped 95% confidence interval
                                r0evokedconf(cll,l,sz,frq,2) = sorted(9500);
                                a = allS(thisstillinds,frq);
                                b = allS(spontstillinds,frq);
                                varratio(frq) = (nanmean(a)^2*nanvar(b)+nanmean(b)^2*nanvar(a))/nanmean(b)^4;
                                r0evokedspecterr(cll,l,sz,frq) = sqrt(varratio(frq))/sqrt(length(spontstillinds)); % the math way
                            end
                        else
                            r0evokedspecterr(cll,l,sz,:) = nan(1,103);
                            r0evokedconf(cll,l,sz,:,:) = nan(103,2);
                        end                            
                    else
                        thisstilln(cll,l,sz) = 0;
                        condstillS(cll,l,sz,:) = nan(1,150); condstillSerr(cll,l,sz,:,:) = nan(2,150);
                        r0condallS(cll,l,sz,:) = nan(1,size(allS,2));
                        r0condallSerr(cll,l,sz,:) = nan(1,size(allS,2));
                        r0condfftspect(cll,l,sz,:) = nan(1,size(allS,2));
                        r0condfftspecterr(cll,l,sz,:) = nan(1,size(allS,2));
                        r0evokedspect(cll,l,sz,:) = nan(1,103);
                        r0evokedspecterr(cll,l,sz,:) = nan(1,103);
                        r0evokedconf(cll,l,sz,:,:) = nan(103,2);
                    end
                    
                    %TODO fix to be Hz
                    condfiltresp(cll,l,sz,:) = filter(excit_kernel,1,mean(resp(thisinds,:),1));
                    
                    rp = randperm(length(thisinds));
                    for i = 1:length(thisinds)
                        s = find(resp(thisinds(i),respwin)); % find spikes from the actual trial
                        if ~isempty(s)
                            for si = 1:length(s) % take LFPs from random other smae size trial
                                lfps(si,:) = lfpresp(thisinds(i),s(si)+respwin(1)-1-200:s(si)+respwin(1)-1+200);
                            end
                            for si = 1:length(s) % take LFPs from random other smae size trial
                                shuflfps(si,:) = lfpresp(thisinds(rp(i)),s(si)+respwin(1)-1-200:s(si)+respwin(1)-1+200);
                            end
                        else
                            lfps = nan(1,401);
                            shuflfps = nan(1,401);
                        end
                        
                        stavglfp(l,sz,i,:) = mean(lfps,1);
                        nspkssta(l,sz,i) = size(lfps,1);
                        shufstavglfp(l,sz,i,:) = mean(shuflfps,1);
                        nspksshufsta(l,sz,i) = size(shuflfps,1);
                        clear lfps; clear shuflfps;
                    end
                    
                    g1phasmat = g1(thisinds,respwin);
                    g2phasmat = g2(thisinds,respwin);
                    r1phasmat = g1(thisruninds,respwin);
                    r0phasmat = g1(thisstillinds,respwin);
                    r1g2phasmat = g2(thisruninds,respwin);
                    r0g2phasmat = g2(thisstillinds,respwin);
                    allg1phases = squeeze(g1phasmat(find(g1phasmat)));
                    allg2phases = squeeze(g2phasmat(find(g2phasmat)));
                    allr1phases = squeeze(r1phasmat(find(r1phasmat)));
                    allr0phases = squeeze(r0phasmat(find(r0phasmat)));
                    allr1g2phases = squeeze(r1g2phasmat(find(r1g2phasmat)));
                    allr0g2phases = squeeze(r0g2phasmat(find(r0g2phasmat)));
                    g1inds = find(g1phasmat); g2inds = find(g2phasmat);
                    if isrow(allg1phases) allg1phases = allg1phases'; end
                    if isrow(allg2phases) allg2phases = allg2phases'; end
                    if isrow(allr1phases) allr1phases = allr1phases'; end
                    if isrow(allr0phases) allr0phases = allr0phases'; end
                    if isrow(allr1g2phases) allr1g2phases = allr1g2phases'; end
                    if isrow(allr0g2phases) allr0g2phases = allr0g2phases'; end
                    g1orimeanphases{cll,l,sz} = allg1phases;
                    g2orimeanphases{cll,l,sz} = allg2phases;
                    r1orimeanphases{cll,l,sz} = allr1phases;
                    r0orimeanphases{cll,l,sz} = allr0phases;
                    r1g2orimeanphases{cll,l,sz} = allr1g2phases;
                    r0g2orimeanphases{cll,l,sz} = allr0g2phases;
                    prefg1phase(cll,l,sz) = circ_mean(allg1phases);
                    prefg2phase(cll,l,sz) = circ_mean(allg2phases);
                    prefr1phase(cll,l,sz) = circ_mean(allr1phases);
                    prefr0phase(cll,l,sz) = circ_mean(allr0phases);
                    prefr1g2phase(cll,l,sz) = circ_mean(allr1g2phases);
                    prefr0g2phase(cll,l,sz) = circ_mean(allr0g2phases);
                    g1lockpval(cll,l,sz) = circ_rtest(allg1phases);
                    g2lockpval(cll,l,sz) = circ_rtest(allg2phases);
                    r1lockpval(cll,l,sz) = circ_rtest(allr1phases);
                    r0lockpval(cll,l,sz) = circ_rtest(allr0phases);
                    r1g2lockpval(cll,l,sz) = circ_rtest(allr1g2phases);
                    r0g2lockpval(cll,l,sz) = circ_rtest(allr0g2phases);
                    g1orimeanplv(cll,l,sz) = circ_r(allg1phases);
                    g2orimeanplv(cll,l,sz) = circ_r(allg2phases);
                    r1orimeanplv(cll,l,sz) = circ_r(allr1phases);
                    r0orimeanplv(cll,l,sz) = circ_r(allr0phases);
                    r1g2orimeanplv(cll,l,sz) = circ_r(allr1g2phases);
                    r0g2orimeanplv(cll,l,sz) = circ_r(allr0g2phases);
                    r1orimeanppc(cll,l,sz) = ppc(allr1phases);
                    r0orimeanppc(cll,l,sz) = ppc(allr0phases);
                    r1g2orimeanppc(cll,l,sz) = ppc(allr1g2phases);
                    r0g2orimeanppc(cll,l,sz) = ppc(allr0g2phases);
                    [g1angdev(cll,l,sz),g1cstd(cll,l,sz)] = circ_std(allg1phases);
                    [g2angdev(cll,l,sz),g2cstd(cll,l,sz)] = circ_std(allg2phases);
                    [r1angdev(cll,l,sz),r1cstd(cll,l,sz)] = circ_std(allr1phases);
                    [r0angdev(cll,l,sz),r0cstd(cll,l,sz)] = circ_std(allr0phases);
                    [r1g2angdev(cll,l,sz),r1g2cstd(cll,l,sz)] = circ_std(allr1g2phases);
                    [r0g2angdev(cll,l,sz),r0g2cstd(cll,l,sz)] = circ_std(allr0g2phases);
                    nspikesg(cll,l,sz) = length(allg1phases);
                    nr1spikes(cll,l,sz) = length(allr1phases);
                    nr0spikes(cll,l,sz) = length(allr0phases);
                    r1ntrials(cll,l,sz) = length(thisruninds);
                    r0ntrials(cll,l,sz) = length(thisstillinds);
                    for iii = 1:10000
                        ri = rand(1,nspikesg(cll,l,sz));
                        ri = ri.*(2*pi);
                        randp(iii) = circ_rtest(ri);
                        randr(iii) = circ_r(ri');
                    end
                    sortrandr = sort(randr);
                    bootstrapr(cll,l,sz) = sortrandr(9500);
                    g1bslocksig(cll,l,sz) = g1orimeanplv(cll,l,sz)>bootstrapr(cll,l,sz);
                    g2bslocksig(cll,l,sz) = g2orimeanplv(cll,l,sz)>bootstrapr(cll,l,sz);
                    dfsg1 = circ_dist(allg1phases,prefg1phase(cll,l,sz));
                    dfsg2 = circ_dist(allg2phases,prefg2phase(cll,l,sz));
                    npbg1(cll,l,sz) = floor(length(dfsg1)/8); % n spikes per bin (but leaves up to 7 least locked cells out)
                    npbg2(cll,l,sz) = floor(length(dfsg2)/8); 
                    [sdsg1,srtig1] = sort(abs(dfsg1)); % sort differences in ascending order, remember order
                    [sdsg2,srtig2] = sort(abs(dfsg2));
                    for pb = 1:8
                        %equal n spikes bins
                        hh = zeros(length(thisinds),length(respwin));
                        hh(g1inds(srtig1((pb-1)*npbg1(cll,l,sz)+1:pb*npbg1(cll,l,sz)))) = 1;
                        en_respg1(pb,thisinds,:) = hh;
                        hh = zeros(length(thisinds),length(respwin));
                        hh(g2inds(srtig2((pb-1)*npbg2(cll,l,sz)+1:pb*npbg2(cll,l,sz)))) = 1;
                        en_respg2(pb,thisinds,:) = hh;
                        
                        if ~(npbg1(cll,l,sz) == 0)
                            uptodfg1(cll,pb,l,sz) = sdsg1(pb*npbg1(cll,l,sz)); %up to w hat phase difference                            
                            uptodfg2(cll,pb,l,sz) = sdsg2(pb*npbg2(cll,l,sz));
                        else
                            uptodfg1(cll,pb,l,sz) = NaN; uptodfg2(cll,pb,l,sz) = NaN;
                        end
                        
                        % correct light condition to same number of spikes as no light condition when there are more spikes. If there are less, leave it
                        ff = zeros(length(thisinds),length(respwin));
                        if l == 2 & npbg1(cll,1,sz)<npbg1(cll,2,sz)
                            allthisdist = srtig1((pb-1)*npbg1(cll,l,sz)+1:pb*npbg1(cll,l,sz));
                            arp = randperm(length(allthisdist));
                            ff(g1inds(srtig1(arp(1:npbg1(cll,1,sz))))) = 1; 
                            lightcorrectedg1(cll,sz) = 1;
                        else
                            ff(g1inds(srtig1((pb-1)*npbg1(cll,l,sz)+1:pb*npbg1(cll,l,sz)))) = 1;
                        end
                        en_lc_respg1(pb,thisinds,:) = ff; 
                        ff = zeros(length(thisinds),length(respwin));
                        if l == 2 & npbg2(cll,1,sz)<npbg2(cll,2,sz)                            
                            allthisdist = srtig2((pb-1)*npbg2(cll,l,sz)+1:pb*npbg2(cll,l,sz));
                            arp = randperm(length(allthisdist));
                            ff(g2inds(srtig2(arp(1:npbg2(cll,1,sz))))) = 1;   
                            lightcorrectedg2(cll,sz) = 1;              
                        else
                            ff(g2inds(srtig1((pb-1)*npbg2(cll,l,sz)+1:pb*npbg2(cll,l,sz)))) = 1;
                        end
                        en_lc_respg2(pb,thisinds,:) = ff;
                        
                        % equal spaced bins
                        gg = zeros(length(thisinds),length(respwin));
                        gg(g1inds(abs(dfsg1)<=pb*pi/8 & abs(dfsg1)>(pb-1)*pi/8)) = 1;
                        es_respg1(pb,thisinds,:) = gg;
                        gg = zeros(length(thisinds),length(respwin));
                        gg(g2inds(abs(dfsg2)<=pb*pi/8 & abs(dfsg2)>(pb-1)*pi/8)) = 1;
                        es_respg2(pb,thisinds,:) = gg;                      
                    end
                    
                    % field trip values
                    % field trip messing around
                    [ftspect(cll,l,sz,:),ftphases{cll,l,sz},ftfax,ftralp(cll,l,sz,:),...
                        ftppc(cll,l,sz,:),ftplv(cll,l,sz,:)] = get_ft_spectstats(lfpresp,resp,thisinds);
                    
                    % running only
                    [ftr1spect(cll,l,sz,:),ftr1phases{cll,l,sz},ftr1fax,ftr1ralp(cll,l,sz,:),...
                        ftr1ppc(cll,l,sz,:),ftr1plv(cll,l,sz,:)] = get_ft_spectstats(lfpresp,resp,thisruninds);
                                        
                    % still only
                    [ftr0spect(cll,l,sz,:),ftr0phases{cll,l,sz},ftr0fax,ftr0ralp(cll,l,sz,:),...
                        ftr0ppc(cll,l,sz,:),ftr0plv(cll,l,sz,:)] = get_ft_spectstats(lfpresp,resp,thisstillinds);
                    
                end
                
                thisinds = find(result.gratingInfo.size == 0 & result.light == l-1);
                thisruninds = intersect(thisinds,oktrials); thisstillinds = intersect(thisinds,stilltrials);
                g1phasmat = g1(thisinds,respwin);
                g2phasmat = g2(thisinds,respwin);
                r1phasmat = g1(thisruninds,respwin);
                r0phasmat = g1(thisstillinds,respwin);
                r1g2phasmat = g2(thisruninds,respwin);
                r0g2phasmat = g2(thisstillinds,respwin);
                allg1phases = squeeze(g1phasmat(find(g1phasmat)));
                allg2phases = squeeze(g2phasmat(find(g2phasmat)));
                allr1phases = squeeze(r1phasmat(find(r1phasmat)));
                allr0phases = squeeze(r0phasmat(find(r0phasmat)));
                allr1g2phases = squeeze(r1g2phasmat(find(r1g2phasmat)));
                allr0g2phases = squeeze(r0g2phasmat(find(r0g2phasmat)));
                g1inds = find(g1phasmat); g2inds = find(g2phasmat);
                if isrow(allg1phases) allg1phases = allg1phases'; end
                if isrow(allg2phases) allg2phases = allg2phases'; end
                if isrow(allr1phases) allr1phases = allr1phases'; end
                if isrow(allr0phases) allr0phases = allr0phases'; end
                if isrow(allr1g2phases) allr1g2phases = allr1g2phases'; end
                if isrow(allr0g2phases) allr0g2phases = allr0g2phases'; end
                prefg1phasectrl(cll,l) = circ_mean(allg1phases);
                prefg2phasectrl(cll,l) = circ_mean(allg2phases);
                prefr1phasectrl(cll,l) = circ_mean(allr1phases);
                prefr0phasectrl(cll,l) = circ_mean(allr0phases);
                prefr1g2phasectrl(cll,l) = circ_mean(allr1g2phases);
                prefr0g2phasectrl(cll,l) = circ_mean(allr0g2phases);
                g1lockpvalctrl(cll,l) = circ_rtest(allg1phases);
                g2lockpvalctrl(cll,l) = circ_rtest(allg2phases);
                r1lockpvalctrl(cll,l) = circ_rtest(allr1phases);
                r0lockpvalctrl(cll,l) = circ_rtest(allr0phases);
                r1g2lockpvalctrl(cll,l) = circ_rtest(allr1g2phases);
                r0g2lockpvalctrl(cll,l) = circ_rtest(allr0g2phases);
                g1orimeanplvctrl(cll,l) = circ_r(allg1phases);
                g2orimeanplvctrl(cll,l) = circ_r(allg2phases);
                r1orimeanplvctrl(cll,l) = circ_r(allr1phases);
                r0orimeanplvctrl(cll,l) = circ_r(allr0phases);
                r1g2orimeanplvctrl(cll,l) = circ_r(allr1g2phases);
                r0g2orimeanplvctrl(cll,l) = circ_r(allr0g2phases);
                r1orimeanppcctrl(cll,l) = ppc(allr1phases);
                r0orimeanppcctrl(cll,l) = ppc(allr0phases);
                r1g2orimeanppcctrl(cll,l) = ppc(allr1g2phases);
                r0g2orimeanppcctrl(cll,l) = ppc(allr0g2phases);
                [g1angdevctrl(cll,l),g1cstdctrl(cll,l)] = circ_std(allg1phases);
                [g2angdevctrl(cll,l),g2cstdctrl(cll,l)] = circ_std(allg2phases);
                [r1angdevctrl(cll,l),r1cstdctrl(cll,l)] = circ_std(allr1phases);
                [r0angdevctrl(cll,l),r0cstdctrl(cll,l)] = circ_std(allr0phases);
                [r1g2angdevctrl(cll,l),r1g2cstdctrl(cll,l)] = circ_std(allr1g2phases);
                [r0g2angdevctrl(cll,l),r0g2cstdctrl(cll,l)] = circ_std(allr0g2phases);
                nspikesgctrl(cll,l) = length(allg1phases);
                nspikesr1ctrl(cll,l) = length(allr1phases);
                nspikesr0ctrl(cll,l) = length(allr0phases);
                dfsg1 = circ_dist(allg1phases,prefg1phasectrl(cll,l));
                dfsg2 = circ_dist(allg2phases,prefg2phasectrl(cll,l));
                npbctrlg1(cll,l) = floor(length(dfsg1)/8); % n spikes per bin (but leaves up to 7 least locked cells out)
                npbctrlg2(cll,l) = floor(length(dfsg2)/8);
                [sdsg1,srtig1] = sort(abs(dfsg1)); % sort differences in ascending order, remember order
                [sdsg2,srtig2] = sort(abs(dfsg2));
                for pb = 1:8
                    hh = zeros(length(thisinds),length(respwin));
                    hh(g1inds(srtig1((pb-1)*npbctrlg1(cll,l)+1:pb*npbctrlg1(cll,l)))) = 1;
                    en_respg1(pb,thisinds,:) = hh;
                    hh = zeros(length(thisinds),length(respwin));
                    hh(g2inds(srtig2((pb-1)*npbctrlg2(cll,l)+1:pb*npbctrlg2(cll,l)))) = 1;
                    en_respg2(pb,thisinds,:) = hh;
                    if ~(npbctrlg1(cll,l) == 0)
                        uptodfctrlg1(cll,pb,l) = sdsg1(pb*npbctrlg1(cll,l)); %up to w hat phase difference
                        uptodfctrlg2(cll,pb,l) = sdsg2(pb*npbctrlg2(cll,l));
                    else
                        uptodfctrlg1(cll,pb,l) = NaN; uptodfctrlg2(cll,pb,l) = NaN;
                    end
                    
                    % correct light condition to same number of spikes as no light condition when there are more spikes. If there are less, leave it
                    ff = zeros(length(thisinds),length(respwin));
                    if l == 2 & npbctrlg1(cll,1)<npbctrlg1(cll,2)
                        allthisdist = srtig1((pb-1)*npbctrlg1(cll,l)+1:pb*npbctrlg1(cll,l));
                        arp = randperm(length(allthisdist));
                        ff(g1inds(srtig1(arp(1:npbctrlg1(cll,1))))) = 1;
                        lightcorrectedctrlg1(cll,sz) = 1;
                    else
                        ff(g1inds(srtig1((pb-1)*npbctrlg1(cll,l)+1:pb*npbctrlg1(cll,l)))) = 1;
                    end
                    en_lc_respg1(pb,thisinds,:) = ff;
                    ff = zeros(length(thisinds),length(respwin));
                    if l == 2 & npbctrlg2(cll,1)<npbctrlg2(cll,2)                        
                        allthisdist = srtig2((pb-1)*npbctrlg2(cll,l)+1:pb*npbctrlg2(cll,l));
                        arp = randperm(length(allthisdist));
                        ff(g1inds(srtig2(arp(1:npbctrlg2(cll,1))))) = 1;
                        lightcorrectedctrlg1(cll,sz) = 1;
                    else
                        ff(g2inds(srtig1((pb-1)*npbctrlg2(cll,l)+1:pb*npbctrlg2(cll,l)))) = 1;
                    end
                    en_lc_respg2(pb,thisinds,:) = ff;
                    
                    % equal spaced bins
                    gg = zeros(length(thisinds),length(respwin));
                    gg(g1inds(abs(dfsg1)<=pb*pi/8 & abs(dfsg1)>(pb-1)*pi/8)) = 1;
                    es_respg1(pb,thisinds,:) = gg;
                    gg = zeros(length(thisinds),length(respwin));
                    gg(g2inds(abs(dfsg2)<=pb*pi/8 & abs(dfsg2)>(pb-1)*pi/8)) = 1;
                    es_respg2(pb,thisinds,:) = gg;
                end
            end       
            en_vstfrg1 = sum(en_respg1,3); en_vstfrg2 = sum(en_respg2,3);
            es_vstfrg1 = sum(es_respg1,3); es_vstfrg2 = sum(es_respg2,3);
            en_lc_vstfrg1 = sum(en_lc_respg1,3); en_lc_vstfrg2 = sum(en_lc_respg2,3);
            
            chf = chf(1:150);
            
            condstalfp(cll,:,:,:) = nanmean(stavglfp,3);
            shufcondsizestalfp(cll,:,:,:) = nanmean(shufstavglfp,3);
            condnspkssta(cll,:,:) = nanmean(nspkssta,3);
            
            
            
            %old version, delete if the other works
%             sizes = unique(result.gratingInfo.size); sizes(find(sizes == 0)) = []; 
%             for i = 1:length(msstamps)
%                 hlp = gamma1phases(i,respwin);
%                 l = result.light(i)+1;
%                 if result.gratingInfo.size(i) == 0
%                     pphase = prefg1phasectrl(cll,l);
%                 else
%                     sz = find(sizes == result.gratingInfo.size(i));                
%                     pphase = prefg1phase(cll,l,sz);
%                 end
%                 if ~isempty(find(hlp))
%                     dfs = circ_dist(hlp(find(hlp)),pphase);
%                 else
%                     dfs = NaN;
%                 end
%                 %bin other virtual spiketrains
%                 for pb = 1:8 % 
%                     es_vstfrg1(pb,i) = length(find(abs(dfs)> (pb-1)*pi/8 & abs(dfs)<= pb*pi/8));
%                     es_vstfrg2(pb,i) = length(find(abs(dfs)> (pb-1)*pi/8 & abs(dfs)<= pb*pi/8));
%                 end
%             end
            
            clear cfr; clear condresp; clear cgamma1resp; clear cgamm2resp; clear clfpresp; clear clfpspect;
            clear clfpspecterr; clear clfpoffsspect; clear cgpow1; clear cgpow2; clear cbr; clear allfr; clear cerr;
            clear bincondresp; clear cmsl; clear csta; clear cstan;
            binwidth = 30;            
            sizes = unique(result.gratingInfo.size);  sizes(find(sizes == 0)) = []; %delete control condition
            oris = unique(result.gratingInfo.Orientation);  oris(find(oris == -1)) = [];
            for l = 1:length(unique(result.light))
                for sz = 1:length(sizes)
                    for ori = 1:length(oris)
                        thisinds = find(result.gratingInfo.Orientation == oris(ori) &...
                            result.gratingInfo.size == sizes(sz) & ...
                            result.light == l-1);
                        condresp(l,ori,sz,:) = mean(resp(thisinds,:),1);
                        cgamma1resp(l,ori,sz,:) = mean(gamma1resp(thisinds,:),1);
                        cgamma2resp(l,ori,sz,:) = mean(gamma2resp(thisinds,:),1);
                        clfpresp(l,ori,sz,:) = mean(lfpresp(thisinds,:),1);
                        clfpspect(l,ori,sz,:) = nanmean(lfpspect(thisinds,:),1);
                        clfpspecterr(l,ori,sz,:) = nanstd(lfpspect(thisinds,:),1,1)./sqrt(length(thisinds));
                        clfpoffsspect(l,ori,sz,:) = nanmean(lfpoffsspect(thisinds,:),1);
                        cgpow1(l,ori,sz) = mean(gpower1(thisinds));
                        cgpow2(l,ori,sz) = mean(gpower2(thisinds));
                        cfr(l,ori,sz) = mean(frs(thisinds));%-mean(bl); 
                        cbr(l,ori,sz) = mean(brs(thisinds));
                        allfr{l,ori,sz} = frs(thisinds);
                        cerr(l,ori,sz) =std(frs(thisinds))./sqrt(length(thisinds));
                        [bincondresp(l,ori,sz,:),bta] = binit(squeeze(condresp(l,ori,sz,:)),binwidth); 
                        cmsl(l,ori,sz) = mean(msl(thisinds)); % mean spike latency
                        csta(l,ori,sz,:) = nanmean(stalfp(thisinds,:),1);
                        cstan(l,ori,sz) = nanmean(stan(thisinds));
                        
                        for pb = 1:8
                            cvstfrg1(pb,l,ori,sz) = nanmean(vstfrg1(pb,thisinds),2);
                            cvstfrg2(pb,l,ori,sz) = nanmean(vstfrg2(pb,thisinds),2);
                            cenvstfrg1(pb,l,ori,sz) = nanmean(en_vstfrg1(pb,thisinds),2);
                            cenvstfrg2(pb,l,ori,sz) = nanmean(en_vstfrg2(pb,thisinds),2);
                            cesvstfrg1(pb,l,ori,sz) = nanmean(es_vstfrg1(pb,thisinds),2);
                            cesvstfrg2(pb,l,ori,sz) = nanmean(es_vstfrg2(pb,thisinds),2);
                            cenlcvstfrg1(pb,l,ori,sz) = nanmean(en_lc_vstfrg1(pb,thisinds),2);
                            cenlcvstfrg2(pb,l,ori,sz) = nanmean(en_lc_vstfrg2(pb,thisinds),2);
                        end
                        
                        cisi{l,ori,sz} = [];
                        for i = 1:length(thisinds)
                            cisi{l,ori,sz} = [cisi{l,ori,sz},isiresp{thisinds(i)}];
                        end
                        % test each cell for each condition ranksum
                        if l == 1
                            l1inds = find(result.gratingInfo.Orientation == oris(ori) &...
                            result.gratingInfo.size == sizes(sz) & ...
                            result.light == 1);
                            if ~isempty(l1inds) & ~isempty(thisinds)
                                cdiffp(ori,sz) = ranksum(frs(thisinds),frs(l1inds));
                            else
                                cdiffp(ori,sz) = NaN;
                            end
                        end
                        
                        condz(l,ori,sz) = {(sc(thisinds)-mean(sc(thisinds)))/std(sc(thisinds))}; %ecker 2010
                        condsc(l,ori,sz) = {sc(thisinds)};
                        ff(l,ori,sz) = var(sc(thisinds))/mean(sc(thisinds));
                        
                        mscc = []; bincc = [];
                        for ii = 1:length(thisinds)-1
                            for jj = ii+1:length(thisinds)
                                help = corrcoef(resp(thisinds(ii),:),resp(thisinds(jj),:));
                                mscc = [mscc,help(1,2)];
                                help = corrcoef(binit(resp(thisinds(ii),:),binwidth),binit(resp(thisinds(jj),:),binwidth));
                                bincc = [bincc, help(1,2)];
                            end
                        end
                        msreliab(l,ori,sz) = nanmean(mscc);
                        binreliab(l,ori,sz) = nanmean(bincc);
                        eckerreliability(l,ori,sz) = var(frs(thisinds))/var(frs);
                                                                
                        for b = 1:size(phasmat,1)
                            condpowresp(l,ori,sz,b,:) = squeeze(mean(allpowresp(b,thisinds,:),2));
                            condpow(l,ori,sz,b) = squeeze(mean(mean(allpowresp(b,thisinds,respwin),2),3));
                        end
                        
                        %phase locking to all bands and the two peak bands
                        phasemat = tmpi(:,thisinds,respwin);
                        for b = 1:size(allphaseresp,1)
                            allphases = squeeze(phasemat(b,find(squeeze(phasemat(b,:,:)))));
                            condr(cll,b,l,ori,sz) = circ_r(allphases');
                            condppc(cll,b,l,ori,sz) = ppc(allphases);
                            condcmean(cll,b,l,ori,sz) = circ_mean(allphases');
                            condallphases{cll,b,l,ori,sz} = allphases';
                        end
                        g1phasmat = g1(thisinds,respwin);
                        g2phasmat = g2(thisinds,respwin);
                        allg1phases = squeeze(g1phasmat(find(g1phasmat)));
                        allg2phases = squeeze(g2phasmat(find(g2phasmat)));
                        if isrow(allg1phases) allg1phases = allg1phases'; end
                        if isrow(allg2phases) allg2phases = allg2phases'; end
                        condg1phases{cll,l,ori,sz} = allg1phases;
                        condg2phases{cll,l,ori,sz} = allg2phases;
                        condg1rp(cll,l,ori,sz) = circ_rtest(allg1phases);
                        condg2rp(cll,l,ori,sz) = circ_rtest(allg2phases);
                        condg1r(cll,l,ori,sz) = circ_r(allg1phases);
                        condg1ppc(cll,l,ori,sz) = ppc(allg1phases);
                        condg1cmean(cll,l,ori,sz) = circ_mean(allg1phases);
                        [condg1angdev(cll,l,ori,sz),condg1cstd(cll,l,ori,sz)] = circ_std(allg1phases);
                        condg2r(cll,l,ori,sz) = circ_r(allg2phases);
                        condg2ppc(cll,l,ori,sz) = ppc(allg2phases);
                        condg2cmean(cll,l,ori,sz) = circ_mean(allg2phases);
                        [condg2angdev(cll,l,ori,sz),condg1cstd(cll,l,ori,sz)] = circ_std(allg2phases);

                        % running
                        thisruninds = intersect(thisinds,oktrials);
                        if ~isempty(thisruninds)
                            runcondresp(l,ori,sz,:) = mean(resp(thisruninds,:),1);
                            runcondfr(l,ori,sz) = mean(frs(thisruninds));
                            runconderr(l,ori,sz) = std(frs(thisruninds))./sqrt(length(thisruninds));
                            runclfpspect(l,ori,sz,:) = nanmean(lfpspect(thisruninds,:),1);
                            runcsta(l,ori,sz,:) = nanmean(stalfp(thisruninds,:),1);
                            for b = 1:size(phasmat,1)
                                runcondpow(l,ori,sz,b) = squeeze(mean(mean(allpowresp(b,thisruninds,respwin),2),3));
                            end
                            runcondg1pow(l,ori,sz) = mean(gpower1(thisruninds));
                            runcondg2pow(l,ori,sz) = mean(gpower2(thisruninds));
                            g1runphasmat = g1(thisruninds,respwin);
                            g2runphasmat = g2(thisruninds,respwin);
                            allrung1phases = squeeze(g1runphasmat(find(g1runphasmat)));
                            allrung2phases = squeeze(g2runphasmat(find(g2runphasmat)));
                            if isrow(allrung1phases) allrung1phases = allrung1phases'; end
                            if isrow(allrung2phases) allrung2phases = allrung2phases'; end
                            runcondg1ppc(l,ori,sz) = ppc(allrung1phases);
                            runcondg1cmean(l,ori,sz) = circ_mean(allrung1phases);
                            runcondg2ppc(l,ori,sz) = ppc(allrung2phases);
                            runcondg2cmean(l,ori,sz) = circ_mean(allrung2phases);
                            runcondg1phases{cll,l,ori,sz} = allrung1phases;
                            runcondg2phases{cll,l,ori,sz} = allrung2phases;
                        else
                            runcondresp(l,ori,sz,:) = nan(1,size(resp,2));
                            runcondfr(l,ori,sz) = NaN;
                            runconderr(l,ori,sz) = NaN;
                            runclfpspect(l,ori,sz,:) = nan(1,size(lfpspect,2));
                            runcsta(l,ori,sz,:) = nan(1,301);
                            for b = 1:size(phasmat,1)
                                runcondpow(l,ori,sz,b) = NaN;
                            end
                            runcondg1pow(l,ori,sz) = NaN;
                            runcondg2pow(l,ori,sz) = NaN;
                            runcondg1ppc(l,ori,sz) = NaN;
                            runcondg1cmean(l,ori,sz) = NaN;
                            runcondg2ppc(l,ori,sz) = NaN;
                            runcondg2cmean(l,ori,sz) = NaN;
                            runcondg1phases{cll,l,ori,sz} = NaN;
                            runcondg2phases{cll,l,ori,sz} =NaN;
                        end
                        
                        thisstillinds = intersect(thisinds,stilltrials);
                        if ~isempty(thisstillinds)
                            stillcondresp(l,ori,sz,:) = nanmean(resp(thisstillinds,:),1);
                            stillcondfr(l,ori,sz) = nanmean(frs(thisstillinds));
                            stillconderr(l,ori,sz) = nanstd(frs(thisstillinds))./sqrt(length(thisstillinds));
                            stillclfpspect(l,ori,sz,:) = nanmean(lfpspect(thisstillinds,:),1);
                            stillcsta(l,ori,sz,:) = nanmean(stalfp(thisstillinds,:),1);
                            for b = 1:size(phasmat,1)
                                stillcondpow(l,ori,sz,b) = squeeze(mean(mean(allpowresp(b,thisstillinds,respwin),2),3));
                            end
                            stillcondg1pow(l,ori,sz) = mean(gpower1(thisstillinds));
                            stillcondg2pow(l,ori,sz) = mean(gpower2(thisstillinds));
                            g1stillphasmat = g1(thisstillinds,:);
                            g2stillphasmat = g2(thisstillinds,:);
                            allstillg1phases = g1stillphasmat(find(g1stillphasmat));
                            allstillg2phases = g2stillphasmat(find(g2stillphasmat));
                            if isrow(allstillg1phases) allstillg1phases = allstillg1phases'; end
                            if isrow(allstillg2phases) allstillg2phases = allstillg2phases'; end
                            stillcondg1ppc(l,ori,sz) = ppc(allstillg1phases);
                            stillcondg1cmean(l,ori,sz) = circ_mean(allstillg1phases);
                            stillcondg2ppc(l,ori,sz) = ppc(allstillg2phases);
                            stillcondg2cmean(l,ori,sz) = circ_mean(allstillg2phases);
                            stillcondg1phases{cll,l,ori,sz} = allstillg1phases;
                            stillcondg2phases{cll,l,ori,sz} = allstillg2phases;
                        else
                            stillcondresp(l,ori,sz,:) = nan(1,size(resp,2));
                            stillcondfr(l,ori,sz) = NaN;
                            stillconderr(l,ori,sz) = NaN;
                            stillclfpspect(l,ori,sz,:) = nan(1,size(lfpspect,2));
                            stillcsta(l,ori,sz,:) = nan(1,301);
                            for b = 1:size(phasmat,1)
                                stillcondpow(l,ori,sz,b) = NaN;
                            end
                            stillcondg1pow(l,ori,sz) = NaN;
                            stillcondg2pow(l,ori,sz) = NaN;
                            stillcondg1ppc(l,ori,sz) = NaN;
                            stillcondg1cmean(l,ori,sz) = NaN;
                            stillcondg2ppc(l,ori,sz) = NaN;
                            stillcondg2cmean(l,ori,sz) = NaN;
                            stillcondg1phases{cll,l,ori,sz} = NaN;
                            stillcondg2phases{cll,l,ori,sz} = NaN;
                        end
                    end
                end
            end

            bincondresp = bincondresp.*(1000/binwidth);
            maxdriven = max(max(squeeze(cfr(1,:,:))))-mean(bl);
            
            l0isis = []; l1isis = [];
            for i = find(result.light == 0)
                l0isis = [l0isis, isiresp{i}];
            end
            for i = find(result.light == 1)
                l1isis = [l1isis, isiresp{i}];
            end    
            
            for l = 1:length(unique(result.light))
                % normal case e.g. oris = [0,45,90,135,180,225,270,315] or oris = [0,180];
                if oridist(oris(1),oris(length(oris)/2+1)) == 0
                    for ori = 1:length(oris)./2
                        for sz = 1:length(sizes)
                            orifrs(l,ori,sz,:) = [allfr{l,ori,sz}',allfr{l,ori+length(oris)./2,sz}'];
                        end
                    end
                else
                    % weird rare case where oris = [0,45,90,135]; without shuffling
                    for ori = 1:length(oris)
                        for sz = 1:length(sizes)
                            orifrs(l,ori,sz,:) = allfr{l,ori,sz}';
                        end
                    end
                end
            end
            avl0s1 = reshape(orifrs(1,:,1,:),40,1);
            avl0s5 = reshape(orifrs(1,:,5,:),40,1);
            avl1s1 = reshape(orifrs(2,:,1,:),40,1);
            avl1s5 = reshape(orifrs(2,:,5,:),40,1);            
            h1 = [zeros(40,1);zeros(40,1);ones(40,1);ones(40,1)]; %light
            h2 = [zeros(40,1);ones(40,1);zeros(40,1);ones(40,1)]; %size
            [p,table,stats] = anovan([avl0s1;avl0s5;avl1s1;avl1s5],{h1 h2},'model','full','display','off');
            lightp(cll) = p(1); sizep(cll) = p(2); slip(cll) = p(3);
            
            for sz = 1:length(sizes)
                cvd(sz) = ttest2(frs(find(result.gratingInfo.size == 0 & result.light == 0)),...
                    frs(find(result.gratingInfo.size == sizes(sz) & result.light == 0)));
                if ~isnan(cvd(sz)) & cvd(sz) & mean(frs(find(result.gratingInfo.size == 0 & result.light == 0)))>mean(frs(find(result.gratingInfo.size == sizes(sz) & result.light == 0)))
                    cvd(sz) = -1;
                end
            end           

            if onlymod & (find(cvd) || maxdriven<2) %either not modulated by light or no stim drives more than 2 Hz
                continue;
            end
            
            bincondcllresp(cll,:,:,:,:) = bincondresp;
                        
            nok(cll) = length(oktrials);
            nstill(cll) = length(stilltrials);
            
            g1centerfreq(cll) = f(bpi);
            g2centerfreq(cll) = f(gpi);
            
            condfr(cll,:,:,:) = cfr;
            conderr(cll,:,:,:) = cerr;
            condmsl(cll,:,:,:) = cmsl; % mean spike latency
            condgpow1(cll,:,:,:) = cgpow1;
            condgpow2(cll,:,:,:) = cgpow2;
            condlfpspect(cll,:,:,:,:) = clfpspect;
            condlfpresp(cll,:,:,:,:) = clfpresp;
            condlfpspecterr(cll,:,:,:,:) = clfpspecterr;
            condlfpoffsspect(cll,:,:,:,:) = clfpoffsspect;
            condsta(cll,:,:,:,:) = csta;
            condstan(cll,:,:,:) = cstan;
            
            condbr(cll,:,:,:) = cbr;
            nbursts(cll) = length(bursts);
            nspikes{cll} = [bursts.num_spikes];
            
            condvstfrg1(cll,:,:,:,:) = cvstfrg1;
            condvstfrg2(cll,:,:,:,:) = cvstfrg2;
            condenvstfrg1(cll,:,:,:,:) = cenvstfrg1;
            condenvstfrg2(cll,:,:,:,:) = cenvstfrg2;
            condesvstfrg1(cll,:,:,:,:) = cesvstfrg1;
            condesvstfrg2(cll,:,:,:,:) = cesvstfrg2;
            condenlcvstfrg1(cll,:,:,:,:) = cenlcvstfrg1;
            condenlcvstfrg2(cll,:,:,:,:) = cenlcvstfrg2;
            
            condbandpow(cll,:,:,:,:) = condpow;
%             condbandpowresp(cll,:,:,:,:,:) = condpowresp;
            condgamma1resp(cll,:,:,:,:) = cgamma1resp;
            condgamma2resp(cll,:,:,:,:) = cgamma2resp;
            
             
            r0condfr(cll,:,:,:) = stillcondfr;
            r0conderr(cll,:,:,:) = stillconderr;
            r1condfr(cll,:,:,:) = runcondfr;
            r1conderr(cll,:,:,:) = runconderr; 
            ntrialsr0(cll) = length(stilltrials);
            ntrialsr1(cll) = length(oktrials);     
            r0condlfpspect(cll,:,:,:,:) = stillclfpspect;
            r1condlfpspect(cll,:,:,:,:) = runclfpspect;
            runcondbandpow(cll,:,:,:,:) = runcondpow;
            stillcondbandpow(cll,:,:,:,:) = stillcondpow;
            r1condsta(cll,:,:,:,:) = runcsta;
            r0condsta(cll,:,:,:,:) = stillcsta;
            r1g1pow(cll,:,:,:) = runcondg1pow;
            r1g2pow(cll,:,:,:) = runcondg2pow;
            r0g1pow(cll,:,:,:) = stillcondg1pow;
            r0g2pow(cll,:,:,:) = stillcondg2pow;
            r1g1ppc(cll,:,:,:) = runcondg1ppc;
            r1g2ppc(cll,:,:,:) = runcondg2ppc;
            r0g1ppc(cll,:,:,:) = stillcondg1ppc;
            r0g2ppc(cll,:,:,:) = stillcondg2ppc;
            r1g1cmean(cll,:,:,:) = runcondg1cmean;
            r1g2cmean(cll,:,:,:) = runcondg2cmean;
            r0g1cmean(cll,:,:,:) = stillcondg1cmean;
            r0g2cmean(cll,:,:,:) = stillcondg2cmean;
            
            l0power(cll,:) = squeeze(mean(mean(allpowresp(:,find(result.light == 0),respwin),2),3));
            l1power(cll,:) = squeeze(mean(mean(allpowresp(:,find(result.light == 1),respwin),2),3));
            r1power(cll,:) = squeeze(mean(mean(allpowresp(:,oktrials,respwin),2),3));
            r0power(cll,:) = squeeze(mean(mean(allpowresp(:,stilltrials,respwin),2),3));
            
            l0lfp(cll,:) = squeeze(mean(lfpresp(find(result.light == 0),:),1));
            l1lfp(cll,:) = squeeze(mean(lfpresp(find(result.light == 1),:),1));
            r1lfp(cll,:) = squeeze(mean(lfpresp(oktrials,:),1));
            r0lfp(cll,:) = squeeze(mean(lfpresp(stilltrials,:),1));
            
            l0isi{cll} = l0isis; l1isi{cll} = l1isis;
            
            cllname{cll} = files(fi).name;
            animalno(cll) = animal(blck);
            recording(cll) = blck;
            depth(cll) = result.depth;
            pangle(cll) = penangle(blck);
            
            cellz(cll,:,:,:) = condz;
            cellsc(cll,:,:,:) = condsc;
            cellff(cll,:,:,:) = ff;
            
            celleckerrely(cll,:,:,:) = eckerreliability;
            cellmsrely(cll,:,:,:) = msreliab;
            cellbinrely(cll,:,:,:) = binreliab;
                        
            %determine if cell is modulated by light
            lightmod(cll) = ttest2(frs(find(result.light)),frs(find(result.light == 0)));
            vismod(cll) = ttest(frs(find(result.light == 0 & result.gratingInfo.Contrast ~= 0)),bl(find(result.light == 0 & result.gratingInfo.Contrast ~= 0)));
            if ~isnan(vismod(cll)) & vismod(cll) & mean(frs(find(result.light == 0 &...
                    result.gratingInfo.Contrast ~= 0)))<mean(bl(find(result.light == 0 & result.gratingInfo.Contrast ~= 0)))
                vismod(cll) = -1;
            end
            condvismod(cll,:) = cvd;
            
            wvchan = find(var(result.waveforms) == max(var(result.waveforms)));
            spike = result.waveforms(:,wvchan);
            interpspike = spline(1:32,spike,1:.1:32);
            [adiff(cll),swidth(cll),ptr(cll),eslope(cll)] = spikequant(interpspike);
            
            waveform(cll,:) = spike;
            clustqual(cll) = result.clusterquality;
            
            bfr(cll) = mean(bl);
            lfr(cll) = mean(frs(find(result.light)));
            nlfr(cll) = mean(frs(find(result.light == 0)));
           
            alll0phasemat = tmpi(:,find(result.light == 0),:);
            alll1phasemat = tmpi(:,find(result.light == 1),:);
            for i = 1:size(allphaseresp,1)
                allphasesl0{i} = alll0phasemat(i,find(squeeze(alll0phasemat(i,:,:))));
                allphasesl1{i} = alll1phasemat(i,find(squeeze(alll1phasemat(i,:,:))));
                allrl0(cll,i) = circ_r(allphasesl0{i}');
                allrl1(cll,i) = circ_r(allphasesl1{i}');
                allcmeanl0(cll,i) = circ_mean(allphasesl0{i}');
                allcmeanl1(cll,i) = circ_mean(allphasesl1{i}');
            end
            allg1phasematl0 = g1(result.light == 0,:); 
            allg1phasematl1 = g1(result.light == 1,:); 
            allg2phasematl0 = g2(result.light == 0,:); 
            allg2phasematl1 = g2(result.light == 1,:); 
            allg1phasesl0 = allg1phasematl0(find(squeeze(allg1phasematl0)));
            allg1phasesl1 = allg1phasematl1(find(squeeze(allg1phasematl1)));
            allg2phasesl0 = allg2phasematl0(find(squeeze(allg2phasematl0)));
            allg2phasesl1 = allg2phasematl1(find(squeeze(allg2phasematl1)));
            
            for i = 1:20
                [ccoef,p] = corrcoef(squeeze(mean(allpowresp(i,:,respwin),3)),mean(speed(:,respwin),2));
                speedpowcc(cll,i) = ccoef(1,2);
                speedpowp(cll,i) = p(1,2);                
                
                [l0ccoef,l0p] = corrcoef(squeeze(mean(allpowresp(i,find(result.light ==0),respwin),3)),mean(speed(find(result.light == 0),respwin),2));
                l0speedpowcc(cll,i) = l0ccoef(1,2);
                l0speedpowp(cll,i) = l0p(1,2);
                
                [l1ccoef,l1p] = corrcoef(squeeze(mean(allpowresp(i,find(result.light ==1),respwin),3)),mean(speed(find(result.light == 1),respwin),2));
                l1speedpowcc(cll,i) = l1ccoef(1,2);
                l1speedpowp(cll,i) = l1p(1,2);
            end

%             binwidth = 30;
            [binnedlight,bta] = binit(mean(lightresp),binwidth); binnedlight = binnedlight.*(1000/binwidth);
            [binnednolight,bta] = binit(mean(nolightresp),binwidth); binnednolight = binnednolight.*(1000/binwidth);
            
%             % latency - finer binning
%             [binnedresp,fbta] = binit(mean(resp),5); binnedresp = binnedresp.*(1000/5);
%             strt = find(fbta>prestim,1);
%             threshcross = find(binnedresp(strt+1:end)>(mean(binnedresp(1:strt))+2*std(binnedresp(1:strt))),1);
%             if ~isempty(threshcross), latency(cll) = fbta(threshcross); else latency(cll) = NaN; end
            
            printname = files(fi).name;
            printname(find(printname=='_')) = ' ';

            nolbl = mean(bl(find(~result.light)));
            nolblerr = std(bl(find(~result.light)))./(sqrt(length(find(~result.light))));
            lbl = mean(bl(find(result.light)));
            lblerr = std(bl(find(result.light)))./(sqrt(length(find(result.light))));

            %baselinesubtracted condfr
            blscondfr = squeeze(condfr(cll,:,:,:))-bfr(cll);

            %how much does the firing rate change for each condition
            condchange = squeeze(blscondfr(2,:,:))-squeeze(blscondfr(1,:,:));            
            
            contindsnl = find(result.gratingInfo.size == 0 & result.light == 0);
            controlresp(1,:) = mean(resp(contindsnl,:),1);
            controlfr(cll,1) = mean(frs(contindsnl));
            blscontrolfr(cll,1) = controlfr(cll,1)-bfr(cll);
            controlerr(cll,1) = std(frs(contindsnl))./sqrt(length(contindsnl));            
            controllfpspect(cll,1,:) = nanmean(lfpspect(contindsnl,:));
            controllfpspecterr(cll,1,:) = nanstd(lfpspect(contindsnl,:),1,1)./sqrt(length(contindsnl));
            controlgamma1pow(cll,1) = mean(gpower1(contindsnl,:),1);
            controlgamma2pow(cll,1) = mean(gpower2(contindsnl,:),1);
            controlgamma1powerr(cll,1) = std(gpower1(contindsnl,:),1,1)./sqrt(length(contindsnl));
            controlgamma2powerr(cll,1) = std(gpower2(contindsnl,:),1,1)./sqrt(length(contindsnl)); 
            controlsta(cll,1,:) = nanmean(stalfp(contindsnl,:),1);
            controlstan(cll,1) = nanmean(stan(contindsnl));
            contfftspect(cll,1,:) = nanmean(fftspect(contindsnl,:),1);
            contfftspecterr(cll,1,:) = nanstd(fftspect(contindsnl,:),1,1)./sqrt(length(contindsnl));
            
            phasematl0 = tmpi(:,contindsnl,respwin);
            for b = 1:size(allphaseresp,1)
                controlpowresp(cll,b,1,:) = squeeze(mean(allpowresp(b,contindsnl,:),2));
                controlpow(cll,b,1) = squeeze(mean(mean(allpowresp(b,contindsnl,respwin),2),3));
                allbphasesl0 = squeeze(phasematl0(b,find(squeeze(phasematl0(b,:,:)))));
                controlr(cll,b,1) = circ_r(allbphasesl0');
                controlppc(cll,b,1) = ppc(allbphasesl0);
                controlcmean(cll,b,1) = circ_mean(allbphasesl0');
                controlallphases{cll,b,1} = allbphasesl0';
            end
            g1phasmat = g1(contindsnl,respwin);
            g2phasmat = g2(contindsnl,respwin);
            controlg1phasesl0 = squeeze(g1phasmat(find(g1phasmat)));
            controlg2phasesl0 = squeeze(g2phasmat(find(g2phasmat)));
            controlg1phases{cll,1} = controlg1phasesl0;
            controlg2phases{cll,1} = controlg2phasesl0;
            controlg1rp(cll,1) = circ_rtest(controlg1phasesl0);
            controlg2rp(cll,1) = circ_rtest(controlg2phasesl0);
            controlg1r(cll,1) = circ_r(controlg1phasesl0);
            controlg1ppc(cll,1) = ppc(controlg1phasesl0);
            controlg1cmean(cll,1) = circ_mean(controlg1phasesl0);
            controlg2r(cll,1) = circ_r(controlg2phasesl0);
            controlg2ppc(cll,1) = ppc(controlg2phasesl0);
            controlg2cmean(cll,1) = circ_mean(controlg2phasesl0);              
            
            contindsl = find(result.gratingInfo.size == 0 & result.light == 1);
            controlresp(2,:) = mean(resp(contindsl,:),1);
            controlfr(cll,2) = mean(frs(contindsl));
            blscontrolfr(cll,2) = controlfr(cll,2)-bfr(cll);
            controlerr(cll,2) = std(frs(contindsl))./sqrt(length(contindsl));
            controllfpspect(cll,2,:) = nanmean(lfpspect(contindsl,:),1);
            controllfpspecterr(cll,2,:) = nanstd(lfpspect(contindsl,:),1,1)./sqrt(length(contindsl)); 
            controlgamma1pow(cll,2) = mean(gpower1(contindsl,:),1);
            controlgamma2pow(cll,2) = mean(gpower2(contindsl,:),1);
            controlgamma1powerr(cll,2) = std(gpower1(contindsl,:),1,1)./sqrt(length(contindsl));
            controlgamma2powerr(cll,2) = std(gpower2(contindsl,:),1,1)./sqrt(length(contindsl));  
            controlsta(cll,2,:) = nanmean(stalfp(contindsl,:),1);
            controlstan(cll,2) = nanmean(stan(contindsl));     
            contfftspect(cll,2,:) = nanmean(fftspect(contindsl,:),1);
            contfftspecterr(cll,2,:) = nanstd(fftspect(contindsl,:),1,1)./sqrt(length(contindsl));   
                       
            phasematl1 = tmpi(:,contindsl,respwin);
            for b = 1:size(allphaseresp,1)
                controlpowresp(cll,b,2,:) = squeeze(mean(allpowresp(b,contindsl,:),2));
                controlpow(cll,b,2) = squeeze(mean(mean(allpowresp(b,contindsl,respwin),2),3));
                allbphasesl1 = squeeze(phasematl1(b,find(squeeze(phasematl1(b,:,:)))));
                controlr(cll,b,2) = circ_r(allbphasesl1');
                controlppc(cll,b,2) = ppc(allbphasesl1);
                controlcmean(cll,b,2) = circ_mean(allbphasesl1');
                controlallphases{cll,b,2} = allbphasesl1';
            end
            g1phasmat = g1(contindsl,respwin);
            g2phasmat = g2(contindsl,respwin);
            controlg1phasesl1 = squeeze(g1phasmat(find(g1phasmat)));
            controlg2phasesl1 = squeeze(g2phasmat(find(g2phasmat)));
            controlg1phases{cll,2} = controlg1phasesl1;
            controlg2phases{cll,2} = controlg2phasesl1;
            controlg1rp(cll,2) = circ_rtest(controlg1phasesl1);
            controlg2rp(cll,2) = circ_rtest(controlg2phasesl1);
            controlg1r(cll,2) = circ_r(controlg1phasesl1);
            controlg1ppc(cll,2) = ppc(controlg1phasesl1);
            controlg1cmean(cll,2) = circ_mean(controlg1phasesl1);
            controlg2r(cll,2) = circ_r(controlg2phasesl1);
            controlg2ppc(cll,2) = ppc(controlg2phasesl1);
            controlg2cmean(cll,2) = circ_mean(controlg2phasesl1);              
            
            r0l0continds = intersect(contindsnl,stilltrials);
            r0l1continds = intersect(contindsl,stilltrials);
            r1l0continds = intersect(contindsnl,oktrials);
            r1l1continds = intersect(contindsl,oktrials);
            
            % field trip control
            [contl0ftspect(cll,l,sz,:),contl0ftphases{cll,l,sz},ftfax,contl0ftralp(cll,l,sz,:),...
                contl0ftppc(cll,l,sz,:),contl0ftplv(cll,l,sz,:)] = get_ft_spectstats(lfpresp,resp,contindsnl);
            
            [contl1ftspect(cll,l,sz,:),contl1ftphases{cll,l,sz},ftfax,contl1ftralp(cll,l,sz,:),...
                contl1ftppc(cll,l,sz,:),contl1ftplv(cll,l,sz,:)] = get_ft_spectstats(lfpresp,resp,contindsl);
            
            [r1contl0ftspect(cll,l,sz,:),r1contl0ftphases{cll,l,sz},ftfax,r1contl0ftralp(cll,l,sz,:),...
                r1contl0ftppc(cll,l,sz,:),r1contl0ftplv(cll,l,sz,:)] = get_ft_spectstats(lfpresp,resp,r1l0continds);
            
            [r0contl0ftspect(cll,l,sz,:),r0contl0ftphases{cll,l,sz},ftfax,r0contl0ftralp(cll,l,sz,:),...
                r0contl0ftppc(cll,l,sz,:),r0contl0ftplv(cll,l,sz,:)] = get_ft_spectstats(lfpresp,resp,r0l0continds);
            
            [r1contl1ftspect(cll,l,sz,:),r1contl1ftphases{cll,l,sz},ftfax,r1contl1ftralp(cll,l,sz,:),...
                r1contl1ftppc(cll,l,sz,:),r1contl1ftplv(cll,l,sz,:)] = get_ft_spectstats(lfpresp,resp,r1l1continds);
            
            [r0contl1ftspect(cll,l,sz,:),r0contl1ftphases{cll,l,sz},ftfax,r0contl1ftralp(cll,l,sz,:),...
                r0contl1ftppc(cll,l,sz,:),r0contl1ftplv(cll,l,sz,:)] = get_ft_spectstats(lfpresp,resp,r0l1continds);
            
            % % %
                        
            
            if ~isempty(r0l0continds)
                [contS,chf,contSerr] = mtspectrumc(squeeze(lfpresp(r0l0continds,1001:1800))',params);
                r0l0contS(cll,:) = contS(1:150); r0l0contSerr(cll,:,:) = contSerr(:,1:150);
                r0l0contallS(cll,:) = nanmean(allS(r0l0continds,:),1); 
                r0l0contallSerr(cll,:) = nanstd(allS(r0l0continds,:),1,1)./sqrt(length(r0l0continds));
                r0contfftspect(cll,1,:) = nanmean(fftspect(r0l0continds,:),1);
                r0contfftspect(cll,1,:) = nanstd(fftspect(r0l0continds,:),1,1)./sqrt(length(r0l0continds));
            else
                r0l0contS(cll,:) = nan(1,150); r0l0contSerr(cll,:,:) = nan(2,150);
                r0l0contallS(cll,:) = nan(1,513); r0l0contallSerr(cll,:,:) = nan(1,513);
                r0contfftspect(cll,1,:) = nan(1,513); r0contfftspect(cll,1,:) = nan(1,513);
            end
            if ~isempty(r0l1continds)
                [contS,chf,contSerr] = mtspectrumc(squeeze(lfpresp(r0l1continds,1001:1800))',params);
                r0l1contS(cll,:) = contS(1:150); r0l1contSerr(cll,:,:) = contSerr(:,1:150);
                r0l1contallS(cll,:) = nanmean(allS(r0l1continds,:),1); 
                r0l1contallSerr(cll,:) = nanstd(allS(r0l1continds,:),1,1)./sqrt(length(r0l1continds)); 
                r0contfftspect(cll,2,:) = nanmean(fftspect(r0l1continds,:),1);
                r0contfftspect(cll,2,:) = nanstd(fftspect(r0l1continds,:),1,1)./sqrt(length(r0l1continds));
            else
                r0l1contS(cll,:) = nan(1,150); r0l1contSerr(cll,:,:) = nan(2,150);
                r0l1contallS(cll,:) = nan(1,513); r0l1contallSerr(cll,:,:) = nan(1,513);
                r0contfftspect(cll,2,:) = nan(1,513); r0contfftspect(cll,2,:) = nan(1,513);
            end
            if ~isempty(r1l0continds)
                [contS,chf,contSerr] = mtspectrumc(squeeze(lfpresp(r1l0continds,1001:1800))',params);
                r1l0contS(cll,:) = contS(1:150); r1l0contSerr(cll,:,:) = contSerr(:,1:150);
                r1l0contallS(cll,:) = nanmean(allS(r1l0continds,:),1); 
                r1l0contallSerr(cll,:) = nanstd(allS(r1l0continds,:),1,1)./sqrt(length(r1l0continds));
                r1contfftspect(cll,1,:) = nanmean(fftspect(r1l0continds,:),1);
                r1contfftspect(cll,1,:) = nanstd(fftspect(r1l0continds,:),1,1)./sqrt(length(r0l0continds));
            else
                r1l0contS(cll,:) = nan(1,150); r1l0contSerr(cll,:,:) = nan(2,150);
                r1l0contallS(cll,:) = nan(1,513); r1l0contallSerr(cll,:,:) = nan(1,513);
                r1contfftspect(cll,1,:) = nan(1,513); r1contfftspect(cll,1,:) = nan(1,513);
            end
            if ~isempty(r1l1continds)
                [contS,chf,contSerr] = mtspectrumc(squeeze(lfpresp(r1l1continds,1001:1800))',params);
                r1l1contS(cll,:) = contS(1:150); r1l1contSerr(cll,:,:) = contSerr(:,1:150);   
                r1l1contallS(cll,:) = nanmean(allS(r1l1continds,:),1); 
                r1l1contallSerr(cll,:) = nanstd(allS(r1l1continds,:),1,1)./sqrt(length(r1l1continds));  
                r1contfftspect(cll,2,:) = nanmean(fftspect(r1l1continds,:),1);
                r1contfftspect(cll,2,:) = nanstd(fftspect(r1l1continds,:),1,1)./sqrt(length(r1l1continds));
            else
                r1l1contS(cll,:) = nan(1,150); r1l1contSerr(cll,:,:) = nan(2,150);
                r1l1contallS(cll,:) = nan(1,513); r1l1contallSerr(cll,:,:) = nan(1,513);
                r1contfftspect(cll,2,:) = nan(1,513); r1contfftspect(cll,2,:) = nan(1,513);
            end
            
            r0controllfpspect(cll,1,:) = nan(1,size(lfpspect,2));
            r0controllfpspect(cll,1,:) = nanmean(lfpspect(r0l0continds,:),1);
            r0controllfpspect(cll,2,:) = nan(1,size(lfpspect,2));
            r0controllfpspect(cll,2,:) = nanmean(lfpspect(r0l1continds,:),1);
            r1controllfpspect(cll,1,:) = nan(1,size(lfpspect,2));
            r1controllfpspect(cll,1,:) = nanmean(lfpspect(r1l0continds,:),1);
            r1controllfpspect(cll,2,:) = nan(1,size(lfpspect,2));
            r1controllfpspect(cll,2,:) = nanmean(lfpspect(r1l1continds,:),1);
            r1controlsta(cll,1,:) = nanmean(stalfp(r1l0continds,:),1);
            r1controlsta(cll,2,:) = nanmean(stalfp(r1l1continds,:),1);
            r0controlsta(cll,1,:) = nanmean(stalfp(r0l0continds,:),1);
            r0controlsta(cll,2,:) = nanmean(stalfp(r0l1continds,:),1);
            
            r1controlg1pow(cll,1) = nanmean(gpower1(r1l0continds));
            r0controlg1pow(cll,1) = nanmean(gpower1(r0l0continds));
            r1controlg1powerr(cll,1) = nanstd(gpower1(r1l0continds))./sqrt(length(r1l0continds));
            r0controlg1powerr(cll,1) = nanstd(gpower1(r0l0continds))./sqrt(length(r0l0continds));
            r1controlg2pow(cll,1) = nanmean(gpower2(r1l0continds));
            r0controlg2pow(cll,1) = nanmean(gpower2(r0l0continds));
            r1controlg2powerr(cll,1) = nanmean(gpower2(r1l0continds))./sqrt(length(r1l0continds));
            r0controlg2powerr(cll,1) = nanmean(gpower2(r0l0continds))./sqrt(length(r0l0continds));
            r1controlg1pow(cll,2) = nanmean(gpower1(r1l1continds));
            r0controlg1pow(cll,2) = nanmean(gpower1(r0l1continds));
            r1controlg1powerr(cll,2) = nanstd(gpower1(r1l1continds))./sqrt(length(r1l1continds));
            r0controlg1powerr(cll,2) = nanstd(gpower1(r0l1continds))./sqrt(length(r0l1continds));
            r1controlg2pow(cll,2) = nanmean(gpower2(r1l1continds));
            r0controlg2pow(cll,2) = nanmean(gpower2(r0l1continds));
            r1controlg2powerr(cll,2) = nanmean(gpower2(r1l1continds))./sqrt(length(r1l1continds));
            r0controlg2powerr(cll,2) = nanmean(gpower2(r0l1continds))./sqrt(length(r0l1continds));
            
            if ~isempty(contindsl)
                spontdiffp(cll) = ranksum(frs(contindsnl),frs(contindsl));
            else
                spontdiffp(cll) = NaN;
            end
            conddiffp(cll,:,:) = cdiffp;

            prefsize = find(mean(condfr(cll,1,:,:),3) == max(mean(condfr(cll,1,:,:),3)),1);

            [nloriprefratio(cll), nldirprefratio(cll), nlprefori, meanoril0(cll), nlosi(cll), meandirl0, nldsi(cll)] = getOSI(squeeze(condfr(cll,1,:,prefsize))',oris);
            [loriprefratio(cll), ldirprefratio(cll), lprefori, meanoril1(cll), losi(cll), meandirl1, ldsi(cll)] = getOSI(squeeze(condfr(cll,2,:,prefsize))',oris);

            oneorifr = mean(reshape(blscondfr(:,:,prefsize),2,4,2),3);
            prefori = find(oneorifr(1,:) == max(oneorifr(1,:)),1);
            [preffr(cll,:), prefdir] = max(condfr(cll,:,:,prefsize),[],3); prefdir = prefdir(1);
            ortho = mod(prefori+2,length(oris)/2); if ortho == 0, ortho = length(oris)/2; end
                
            sizetunel1(cll,:) = [controlfr(cll,2), squeeze(nanmean(condfr(cll,2,[prefori,prefori+(length(oris)/2)],:),3))'];
            sizetunel0(cll,:) = [controlfr(cll,1), squeeze(nanmean(condfr(cll,1,[prefori,prefori+(length(oris)/2)],:),3))'];
            sizetuneerrl1(cll,:) = [controlerr(cll,2), squeeze(nanmean(conderr(cll,2,[prefori,prefori+(length(oris)/2)],:),3))'];
            sizetuneerrl0(cll,:) = [controlerr(cll,1), squeeze(nanmean(conderr(cll,1,[prefori,prefori+(length(oris)/2)],:),3))'];
            xsizes = [0,sizes];            
            
            [binprefl0,bta] = binit(squeeze(condresp(1,prefdir,prefsize,:)),binwidth); binprefl0 = binprefl0.*(1000/binwidth);
            [binprefl1,bta] = binit(squeeze(condresp(2,prefdir,prefsize,:)),binwidth); binprefl1 = binprefl1.*(1000/binwidth);
            ta = bta-prestim;            
            
            l0prefresp(cll,:) = binprefl0;
            l1prefresp(cll,:) = binprefl1;
            l0meanresp(cll,:) = binnednolight;
            l1meanresp(cll,:) = binnedlight;
            
            % bin and fit a sine of correct temporal frequency
            binprefl0 = binprefl0-mean(bl);  % subtract baseline
            binprefl1 = binprefl1-mean(bl);
            stwin = find(ta>700&ta<1500);  % only take part of response after transient
            spsigl0 = binprefl0(stwin); spsigl1 = binprefl1(stwin);
            tx = ta(stwin);        % time axis for afer transient
            tempfreq = unique(result.gratingInfo.tFreq);
            [sppl0, uu, rsql0(cll)] = fit_fixedfsin(tx,spsigl0,tempfreq,sr); % the fit parameters ( Amplitude, Phase and Offset)
            [sppl1, uu, rsql1(cll)] = fit_fixedfsin(tx,spsigl1,tempfreq,sr);
            f1f0l0(cll) = (sppl0(1))/sppl0(3); % Amplitude = F1, Offset = F0
            f1f0l1(cll) = (sppl1(1))/sppl1(3);
            f1l0(cll) = sppl0(1); f1l1(cll) = sppl1(1); f0l0(cll) = sppl0(3); f0l1(cll) = sppl1(3);
            
%             figure
            clf;
            subplot(2,2,1)
            ta = bta-prestim;
            plot(ta,binnednolight,'k','linewidth',2);
            hold on
            plot(ta,binnedlight,lcol,'linewidth',2);
            mx = max([max(binnednolight),max(binnedlight),.1]);
            axis([-prestim,trialdur+poststim,0,mx]);
            line([0,0],[0,mx],'color','k','linewidth',2);
            line([2000,2000],[0,mx],'color','k','linewidth',2);
            line([500,500],[0,mx],'color','b','linewidth',2)
            line([1500,1500],[0,mx],'color','b','linewidth',2);
            legend({'Light OFF','Light ON'})
            xlabel('time [ms]')
            ylabel('firing rate [Hz]')
            title(['cell ' int2str(cll) ' depth: ' int2str(result.depth), 'cell ' printname ])

            subplot(2,2,2)
            errorbar(oris,squeeze(condfr(cll,2,:,prefsize)),squeeze(conderr(cll,2,:,prefsize)),'o-','color',lcol,'markersize',8,'linewidth',2)
            hold on
            errorbar(oris,squeeze(condfr(cll,1,:,prefsize)),squeeze(conderr(cll,1,:,prefsize)),'ko-','markersize',8,'linewidth',2)
            xlabel('shown orientation')
            ylabel('Firing rate [Hz]')
            set(gca,'xtick',oris)
            legend({'Light ON','Light OFF'})
            title([' OSI: ' num2str(nlosi(cll)) ' OSI Light: ' num2str(losi(cll))])

            subplot(2,2,3)
            errorbar(xsizes,sizetunel1(cll,:),sizetuneerrl1(cll,:),'o-','color',lcol,'markersize',8,'linewidth',2);
            hold on
            errorbar(xsizes,sizetunel0(cll,:),sizetuneerrl0(cll,:),'ko-','markersize',8,'linewidth',2);
            xlabel('shown patch size [vd]')
            ylabel('Firing rate [Hz]')
            legend({'Light ON','Light OFF'})    
            set(gca,'xtick',sizes)  
            title(['cell ' int2str(cll) ' preferred orientations' ' depth: ' int2str(result.depth) '  ' printname])


            sil(cll) = (sizetunel1(cll,find(sizetunel1(cll,:) == max(sizetunel1(cll,:)),1))-sizetunel1(cll,end))/sizetunel1(cll,find(sizetunel1(cll,:) == max(sizetunel1(cll,:)),1));
            sinl(cll) = (sizetunel0(cll,find(sizetunel0(cll,:) == max(sizetunel0(cll,:)),1))-sizetunel0(cll,end))/sizetunel0(cll,find(sizetunel0(cll,:) == max(sizetunel0(cll,:)),1));

            subplot(2,4,7)
            plot(spike)
            axis([0,40,-100,100])
            legend(['width: ' int2str(swidth(cll)) ' adiff: ' num2str(adiff(cll))])
            
            if printyn
                figSize = [30 21];
                set(gcf,'paperunits','centimeters','papersize',figSize,'paperposition',[0 0 figSize])
                if cll<10, printi = ['0', int2str(cll)]; else printi = int2str(cll); end
                print([printpath ,  printi '__' files(fi).name '.pdf'],'-dpdf')
            end
            
            % running figure
            runlfr(cll) = mean(frs(intersect(find(result.light),oktrials)));
            runnlfr(cll) = mean(frs(intersect(find(~result.light),oktrials)));
            norunlfr(cll) = mean(frs(intersect(find(result.light),stilltrials)));
            norunnlfr(cll) = mean(frs(intersect(find(~result.light),stilltrials)));
            lfrerr = std(frs(find(result.light)))./sqrt(length(find(result.light)));
            nlfrerr = std(frs(find(~result.light)))./sqrt(length(find(~result.light)));
            runlfrerr = std(frs(intersect(find(result.light),oktrials)))./sqrt(length(intersect(find(result.light),oktrials)));
            runnlfrerr = std(frs(intersect(find(~result.light),oktrials)))./sqrt(length(intersect(find(~result.light),oktrials)));
            norunlfrerr = std(frs(intersect(find(result.light),stilltrials)))./sqrt(length(intersect(find(result.light),stilltrials)));
            norunnlfrerr = std(frs(intersect(find(~result.light),stilltrials)))./sqrt(length(intersect(find(~result.light),stilltrials)));
            
            l1r1 = frs(intersect(find(result.light),oktrials));
            l0r1 = frs(intersect(find(~result.light),oktrials));
            l1r0 = frs(intersect(find(result.light),stilltrials));
            l0r0 = frs(intersect(find(~result.light),stilltrials));
            anovavec = [l0r0;l0r1;l1r0;l1r1];
            g1 = [zeros(length(l0r0),1);zeros(length(l0r1),1);ones(length(l1r0),1);ones(length(l1r1),1)]; %light
            g2 = [zeros(length(l0r0),1);ones(length(l0r1),1);zeros(length(l1r0),1);ones(length(l1r1),1)]; %running
            [p,table,stats] = anovan(anovavec,{g1 g2},'model','full','display','off');
            lp(cll) = p(1); rp(cll) = p(2); rlip(cll) = p(3);
            
            r0omi(cll) = (norunlfr(cll)-norunnlfr(cll))/(norunlfr(cll)+norunnlfr(cll));
            r1omi(cll) = (runlfr(cll)-runnlfr(cll))/(runlfr(cll)+runnlfr(cll));
            l0rmi(cll) = (runnlfr(cll)-norunnlfr(cll))/(runnlfr(cll)+norunnlfr(cll));
            l1rmi(cll) = (runlfr(cll)-norunlfr(cll))/(runlfr(cll)+norunlfr(cll));
            
            trialfrl0(cll,:) = frs(find(result.light == 0),:);
            trialfrl1(cll,:) = frs(find(result.light == 1),:);
            trialbl(cll,:) = bl;
%             figure
%             clf
%             subplot(2,2,1)
%             imagesc(speed);
%             colorbar
%             title(['oktrials: ' int2str(length(oktrials)) '/' int2str(size(speed,1))])
%             xlabel('time [ms]')
%             ylabel('trial number')
%             
%             subplot(2,2,2)
%             errorbar(msta,mean(speed(find(~result.light),:)),std(speed(find(~result.light),:))./sqrt(length(find(~result.light))),'b')
%             hold on
%             errorbar(msta,mean(speed(find(result.light),:)),std(speed(find(result.light),:))./sqrt(length(find(result.light))),'r')
%             xlabel('time [ms]')
%             ylabel('average runspeed')
%             legend({'light off' 'light on'})
%             
%             subplot(2,2,3)
%             plot(mean(speed(:,respwin),2),frs,'.')
%             hold on
%             plot(mean(speed(find(result.light),respwin),2),frs(find(result.light)),'r.')
%             xlabel('average runspeed of trial')
%             ylabel('average firing rate of trial')
%             
%             subplot(2,2,4)
%             barweb([nlfr(cll),lfr(cll);runnlfr(cll),runlfr(cll);norunnlfr(cll),norunlfr(cll)],...
%                 [nlfrerr,lfrerr;runnlfrerr,runlfrerr;norunnlfrerr,norunlfrerr],...
%                 [],[{'all'};{'running only'};{'immobile only'}],'firing rates with running',...
%                 [],'firing rate [Hz]',[],[]);
%             
%             if printyn
%                 figSize = [30 21];
%                 set(gcf,'paperunits','centimeters','papersize',figSize,'paperposition',[0 0 figSize])
%                 if cll<10, printi = ['0', int2str(cll)]; else printi = int2str(cll); end
%                 print([runprintpath ,  printi '__' files(fi).name '.pdf'],'-dpdf')
%             end           
            
            disp([files(fi).name '   done'])
            cll = cll + 1;
            
        end
    end
    save(popfile, '-v7.3');
%     , 'cellname','depth','adiff','swidth','ptr','eslope','bfr','lfr','nlfr','oris','xsizes',...
%         'controlfr','controlerr','blscontrolfr', 'nlosi','losi','nldsi','ldsi','nloriprefratio','loriprefratio',...
%         'meanoril0','meanoril1','meandirl0','meandirl1',...
%         'nldirprefratio','ldirprefratio','sizetunel0','sizetunel1','sizetuneerrl0','sizetuneerrl1',...
%         'sil','sinl', 'condfr', 'conderr','preffr','allrl0','allrl1','bta',...
%         'allcmeanl0','allcmeanl1','f1f0l0','f1f0l1','f1l0','f1l1','f0l0','f0l1','xsizes',...
%         'rp','lp','rlip','r0omi','r1omi','l0rmi','l1rmi','nok','nstill','runnlfr','runlfr','norunnlfr','norunlfr',...
%         'r0condfr','r0conderr','r1condfr','r1conderr',...
%         'lightmod','animalno','bincondcellresp',... %,'latency'
%         'l0meanresp','l1meanresp','l0prefresp','l1prefresp',...
%         'condbandpow','runcondbandpow','stillcondbandpow','condr','condppc','condcmean','waveform','clustqual',...
%         'condg1r','condg2r','condg1ppc','condg2ppc','condg1cmean','condg2cmean',...
%         'trialfrl0','trialfrl1','trialbl','vismod','condvismod','cellz','cellsc','cellff','recording',...
%         'celleckerrely','cellbinrely','cellmsrely','condbr','nbursts','nspikes','l0isi','l1isi','conddiffp','spontdiffp',...
%         'ta','ntrialsr0','ntrialsr1','lightp','sizep','slip',...
%         'condmsl','condgpow1','condgpow2','condlfpspect','condlfpspecterr','condlfpoffsspect','controllfpspect','controllfpspecterr','fax',...
%         'r0controllfpspect','r1controllfpspect','r0condlfpspect','r1condlfpspect','condgamma1resp','condgamma2resp',...
%         'controlgamma1pow','controlgamma2pow','controlgamma1powerr','controlgamma2powerr','g1centerfreq','g2centerfreq',...
%         'r1g1pow','r1g2pow','r0g1pow','r0g2pow','r1g1pcc','r1g2pcc','r0g1pcc','r0g2pcc',...
%         'r1g1cmean','r1g2cmean','r0g1cmean','r0g2cmean','pangle',...
%         'r1controlg1pow','r1controlg1powerr','r1controlg2pow','r1controlg2powerr','r0controlg1pow','r0controlg1powerr','r0controlg2pow','r0controlg2powerr');  %,'condbandpowresp' %,'condbandpowresp'
else
    load(popfile);   
end

clear cell;
%spike classification
kmeansind = kmeans([eslope',ptr',swidth',adiff'],4);
% kmeansind = kmeans([swidth',ptr'],3);
secpersamp = 1/30000;
interpf = secpersamp/10;
swidthms = swidth*interpf*1000;

% m1 = mean(swidth(kmeansind==1)); m2 = mean(swidth(kmeansind==2)); m3 = mean(swidth(kmeansind==3));
% fm = min([m1,m2,m3]); fmax = max([m1,m2,m3]);
% if fm == m1
%     pfs = find(kmeansind==1); pfsv = kmeansind==1;
% elseif fm == m2
%     pfs = find(kmeansind==2); pfsv = kmeansind == 2;
% else
%     pfs = find(kmeansind==3); pfsv = kmeansind == 3;
% end
% if fmax == m1
%     prs = find(kmeansind==1); prsv = kmeansind==1;
% elseif fmax == m2
%     prs = find(kmeansind==2); prsv = kmeansind == 2;
% else
%     prs = find(kmeansind==3); prsv = kmeansind == 3;
% end

m1 = mean(swidth(kmeansind==1)); m2 = mean(swidth(kmeansind==2)); m3 = mean(swidth(kmeansind==3));
if find(kmeansind == 4), m4 = mean(swidth(kmeansind == 4)); end
fm = min([m1,m2,m3,m4]);
if fm == m1
    pfs = find(kmeansind==1); prs = find(kmeansind==2|kmeansind==3|kmeansind==4); 
    pfsv = kmeansind==1; prsv = kmeansind==2|kmeansind==3|kmeansind==4;
elseif fm == m2
    pfs = find(kmeansind==2); prs = find(kmeansind==1|kmeansind==3|kmeansind==4);
    pfsv = kmeansind==2; prsv = kmeansind==1|kmeansind==3|kmeansind==4;
elseif fm == m3
    pfs = find(kmeansind==3); prs = find(kmeansind==1|kmeansind==2|kmeansind==4);
    pfsv = kmeansind==3; prsv = kmeansind==1|kmeansind==2|kmeansind==4;
else
    pfs = find(kmeansind==4); prs = find(kmeansind==1|kmeansind==2|kmeansind==3);
    pfsv = kmeansind==4; prsv = kmeansind==1|kmeansind==2|kmeansind==3;
end

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

prsv = swidthms>=.38; pfsv = swidthms<=.36;
prs = find(prsv); pfs = find(pfsv);

swamp = max(waveform,[],2)-min(waveform,[],2);
okwv = swamp>42;
vismod(isnan(vismod)) = 0;

for i = 1:length(okwv)
    isodist(i) = clustqual(i).IsolationDistance;
    lratio(i) = clustqual(i).L_Ratio.Lratio;
end

% ok = vismod&okwv'&nlfr>1;
% ok = okwv'&isodist>8;
ok = okwv';

% adjust depth according to penetration angle
% pangle = 10; %TODO temp fix
depth = depth.*cosd(22).*cosd(pangle);

phe = zeros(1,length(depth));
% %putative halo expressing for Scnn pop responsive only: 
% phe([58,83,196,197]) = 1;

%putative halo expressing for Scnn pop ALL cells: 
% phe([16,71,85,103,106,214,216,260,261,350]) = 1;

% %putative halo expressing for SOM P0 pop: 
% phe([23,33,64,88,132,135,198,287]) = 1;
% phe([22,115,125,126]) = 1; %old

% %putative halo expressing for SOM later pop: 87??? 95???
% phe([7,37,106,119,163,189,232]) = 1;

% % putative halo/ChR2 expressing for SOM Halo/ChR2 pop
% phe([7,37,63,65,81,98,107,  128,133,151,158,166,167,175]) = 1;

% %putative halo expressing for SOM+PV combined pop: (added the first PV Halo 1/7/16)
phe([7,37,106,119,163,189,269,293]) = 1; % took out 330 which is 100 in PV 1/11/16

% %putative halo expressing for PV Halo pop:
% phe([39,63,100]) = 1;   % took 100 out 16/01/11 -> not technically FS

%putative halo expressing for PV HaloeArch pop:
% phe([39,63,100,133,146,148,197,200,214,215,216]) = 1;   % took 100 out 16/01/11 -> not technically FS

%putative halo expressing for PV eArch pop:
% phe([7,9,10]) = 1;

% prsv = zeros(1,length(depth));
% prsv(prs) = 1; 
% prsv = logical(prsv'); 
phe = logical(phe);

[sd,si] = sort(depth);
[rsd,rsi] = sort(depth(prs));
[fsd,fsi] = sort(depth(pfs));

% fsl4 = find(depth(pfs)>=375 & depth(pfs)<=500);
% rsl4 = find(depth(prs)>=375 & depth(prs)<=500);
% fsl23 = find(depth(pfs)<375);
% rsl23 = find(depth(prs)<375);
% fsl5 = find(depth(pfs)>500 & depth(pfs)<=800);
% rsl5 = find(depth(prs)>500 & depth(prs)<=800);
l23 = depth<375;
l4 = depth>=375&depth<=550;
l5 = depth>550&depth<=800;
l5a = depth>550&depth<=650;
l5b = depth>650&depth<=800;
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

% % for SOM Halo later pop - one loc per animal between 350 and 500
% %no right depths for animals 2, 5, 11, 14
% % lfpinds = [31,82,92,115,125,133,146,154,165,178,187,202,214,229]; %previsous pop, changed 4/19/16 - tried 336
% lfpinds = [31,82,93,113,120,134,141,154,164,175,184,206,213,228,237,250];

% for SOM PV mixed pop - shooting for 336
% no right depths for animals 2, 5, 11, 14, 22
lfpinds = [31,82,93,115,125,134,141,154,164,175,184,206,217,228,237,250,252,271,287,300,306,311,323,351,359,370,381,395,410,418,431,447,456,462];

% % for SOM Halo+ChR2 pop - one loc per animal shoot for 336 - no right depth for #2 or 3
% lfpinds = [31,71,82,95,112,125,137,150,176];

% % % for PV Halo pop - one loc per animal between 350 and 500
% % % (animal 2 very supreficial)
% lfpinds = [3,8,21,36,49,55,60,72];
% lfpinds = [6,20,36,49,55,60,72,100]; % close to 330 newer 5/16/16

% % % for PV Halo+eArch pop - one loc per animal between 250+400 - shooting for 336
% % % (animal 2 very supreficial)
% lfpinds = [2,8,20,36,49,55,60,72,92,108,119,130,144,159,167,180,196,205,211];

% % for PV eArch pop - one loc per animal 375 and 350
% lfpinds = [1, 19, 29, 43];

% for SOM control pop - one loc per animal 375 and 350
% lfpinds = [7, 12, 18, 25, 31];

% % for SOM anesth pop - one loc per animal 375 and 350
% lfpinds = [1, 15, 32];
% lgi(animalno == 1) = 52;
% lgi(animalno == 2) = 41;
% lgi(animalno == 3) = 56;

% % for SOM PV combined anesth pop - one loc per animal 375 and 350
% lfpinds = [1, 15, 32, 37, 45];
% lgi(animalno == 1) = 52;
% lgi(animalno == 2) = 41;
% lgi(animalno == 3) = 56;
% lgi(animalno == 4) = 40;
% lgi(animalno == 5) = 43;

% for Scnn populations
% no rigth depths for animals 8, 17
% lfpinds = [32,58,92,123,153,155,233,309,323,387,397,422,424,426,468,476,479,485];

% quick and dirty gamma power mean orientation tuning
for i = 1:468
    [x,prefs(i)] = max(nanmean(condfr(i,1,:,:),3));
    [p, d, o, meanorig1l0(i), g1osi(i), meandirg1l0(i), g1dsi(i)] = getOSI(squeeze(r1g1pow(i,1,:,prefs(i)))',oris);
end
plot(meanorig1l0,meanoril0,'ko','markerfacecolor','k','markersize',4)

bta = bta-300;
hcs = find(phe);
for i = 1:length(hcs)
    figure
%     plot(bta(1:95),squeeze(nanmean(nanmean(bincondcllresp(hcs(i),1,:,:,1:95),3),4)),'k','linewidth',2)
%     hold on
%     plot(bta(1:95),squeeze(nanmean(nanmean(bincondcllresp(hcs(i),2,:,:,1:95),3),4)),'r','linewidth',2)
    plot(-299:2700,squeeze(nanmean(condfiltresp(hcs(i),1,:,:),3)).*1000,'k','linewidth',2)
    hold on
    plot(-299:2700,squeeze(nanmean(condfiltresp(hcs(i),2,:,:),3)).*1000,'r','linewidth',2)
    
end

rsphas = prefr1phase(l23rs,1,5); rsphas(isnan(rsphas)) = [];
fsphas = prefr1phase(l23fs,1,5); fsphas(isnan(fsphas)) = [];
somphas = prefr1phase(phe,1,5); somphas(isnan(somphas)) = [];

figure
bar([rad2deg(uncircle(circ_mean(rsphas))),rad2deg(uncircle(circ_mean(fsphas))),rad2deg(uncircle(circ_mean(somphas)))]);
hold on
errorbar([rad2deg(uncircle(circ_mean(rsphas))),rad2deg(uncircle(circ_mean(fsphas))),rad2deg(uncircle(circ_mean(somphas)))],[rad2deg(circ_std(rsphas)./sqrt(length(rsphas))),rad2deg(circ_std(fsphas)./sqrt(length(fsphas))),rad2deg(circ_std(somphas)./sqrt(length(somphas)))],'k.')

figure
x = 1:360;
plot(x,cosd(x));
hold on
herrorbar(rad2deg(uncircle(circ_mean(rsphas))),...
    cosd(rad2deg(uncircle(circ_mean(rsphas)))),...
    rad2deg(circ_std(rsphas)./sqrt(length(rsphas))),'ko')
herrorbar(rad2deg(uncircle(circ_mean(fsphas))),...
    cosd(rad2deg(uncircle(circ_mean(fsphas)))),...
    rad2deg(circ_std(fsphas)./sqrt(length(fsphas))),'go')
herrorbar(rad2deg(uncircle(circ_mean(somphas))),...
    cosd(rad2deg(uncircle(circ_mean(somphas)))),...
    rad2deg(circ_std(somphas)./sqrt(length(somphas))),'ro')

a = find(l23rs);
for i = 1:length(a)
plot(ftfax,squeeze(ftppc(a(i),1,1,:)),'r')
haxes1 = gca;
set(haxes1,'XColor','r','YColor','r')
haxes1_pos = get(haxes1,'Position');
haxes2 = axes('Position',haxes1_pos,'XaxisLocation','top','YAxisLocation','right','Color','none');
hold on
semilogy(fax(1:103),squeeze(condrunS(a(i),1,1,1:103)),'Parent',haxes2,'Color','k')
legend(int2str(nspikesg(a(i),1,1)))
pause
hold off
end

% for only large sizes
% for i = 1:length(hcs)
%     figure
%     plot(bta(1:95),squeeze(nanmean(bincondcllresp(hcs(i),1,:,5,1:95),3)),'k','linewidth',2)
%     hold on
%     plot(bta(1:95),squeeze(nanmean(bincondcllresp(hcs(i),2,:,5,1:95),3)),'r','linewidth',2)
% end

% PPC spectra
ftppc(isinf(ftppc)) = NaN;
r0contl0ftppc(isinf(r0contl0ftppc)) = NaN;
r0contl1ftppc(isinf(r0contl1ftppc)) = NaN;
r1contl0ftppc(isinf(r1contl0ftppc)) = NaN;
r1contl1ftppc(isinf(r1contl1ftppc)) = NaN;
ftr1ppc(isinf(ftr1ppc)) = NaN;
oknppc = nspikesg(:,1,5)>10;
%!!!!! only for SOM HALO - take out really noisy FS cell
% not done yet

a = find(phe);
for i = 1:length(a)
    figure
    plot(ftfax,squeeze(ftppc(a(i),1,5,:)),'b')
    hold on
    plot(ftfax,squeeze(ftppc(a(i),2,5,:)),'r')
    plot(ftfax(ftralp(a(i),1,5,:)<0.05),squeeze(ftppc(a(i),1,5,ftralp(a(i),1,5,:)<0.05)),'bo')
    plot(ftfax(ftralp(a(i),2,5,:)<0.05),squeeze(ftppc(a(i),2,5,ftralp(a(i),2,5,:)<0.05)),'ro')
end

% high gamma for smalles size
for i = 1:length(a)
    figure
    plot(ftfax,squeeze(ftppc(a(i),1,1,:)),'b')
    hold on
    plot(ftfax,squeeze(ftppc(a(i),2,1,:)),'r')
    plot(ftfax(ftralp(a(i),1,1,:)<0.05),squeeze(ftppc(a(i),1,1,ftralp(a(i),1,1,:)<0.05)),'bo')
    plot(ftfax(ftralp(a(i),2,1,:)<0.05),squeeze(ftppc(a(i),2,1,ftralp(a(i),2,1,:)<0.05)),'ro')
end

% example SOM cell ppc spectrum
figure
plot(ftfax,squeeze(ftppc(189,1,5,:)),'k','linewidth',2)
hold on
plot(ftfax,squeeze(ftppc(189,2,5,:)),'r','linewidth',2)
xlabel('frequency (Hz)')
ylabel('PPC')
legend('control','SOM Halo')

figure
errorbar(ftfax,squeeze(nanmean(ftppc(l23rs,1,5,:))),squeeze(nanstd(ftppc(l23rs,1,5,:)))./sqrt(length(find(~isnan(ftppc(l23rs,1,5,1))))));
hold on
errorbar(ftfax,squeeze(nanmean(ftppc(l23fs,1,5,:))),squeeze(nanstd(ftppc(l23fs,1,5,:)))./sqrt(length(find(~isnan(ftppc(l23fs,1,5,1))))),'g')
errorbar(ftfax,squeeze(nanmean(ftppc(phe,1,5,:))),squeeze(nanstd(ftppc(phe,1,5,:)))./sqrt(length(find(~isnan(ftppc(phe,1,5,1))))),'r')
xlabel('frequency [Hz]')
ylabel('PPC')
legend('RS','FS','SOM')
title('ALL cells')

% PPC spectra with light
a = find(l23rs&oknppc');
fillxppc = [ftfax(3:30),fliplr(ftfax(3:30))];
err = squeeze(nanstd(ftppc(a,1,5,3:30)))./sqrt(length(find(~isnan(ftppc(a,1,5,3)))));
fillyppc0 = [(squeeze(nanmean(ftppc(a,1,5,3:30)))+err)',fliplr((squeeze(nanmean(ftppc(a,1,5,3:30)))-err)')];
err = squeeze(nanstd(ftppc(a,2,5,3:30)))./sqrt(length(find(~isnan(ftppc(a,2,5,3)))));
fillyppc1 = [(squeeze(nanmean(ftppc(a,2,5,3:30)))+err)',fliplr((squeeze(nanmean(ftppc(a,2,5,3:30)))-err)')];

figure
fill(fillxppc,fillyppc1,[1,.3,.3])
hold on
fill(fillxppc,fillyppc0,[0,0,0])
plot(ftfax(3:30),squeeze(nanmean(ftppc(a,1,5,3:30))),'w')
plot(ftfax(3:30),squeeze(nanmean(ftppc(a,2,5,3:30))),'k')
% axis([5,105,-.03,.13])
% legend('SOM Halo','SOM cells')
axis([10,100,-.02,.12])
legend('SOM Halo','FS cells')
% axis([10,100,-.02,.1])
% legend('SOM Halo','RS cells')
% axis([5,105,-.04,.18])
% legend('PV Halo','PV cells')
% axis([10,100,-.02,.14])
% legend('PV Halo','RS cells')
% axis([5,105,-.02,.12])
% legend('PV Halo','FS cells')
xlabel('frequency (Hz)')
ylabel('PPC')
title(['PPC spectrum n = ' int2str(length(a))]);


% norm only max size ppcs
for i = 1:size(ftppc,1)
    [m,ind] = max(ftppc(i,1,5,5:end));
    peakind(i) = ind+4;
    peakamp(i) = ftppc(i,1,5,peakind(i));
    peakampl1(i) = ftppc(i,2,5,peakind(i));
    f30amp(i) = ftppc(i,1,5,9);
    f30ampl1(i) = ftppc(i,2,5,9);
    [m,maxpowfreq(i)] = min(abs(ftfax-fax(lgi(i))));
    mpfamp(i) = ftppc(i,1,5,maxpowfreq(i)); % max pow freq ppc - ppc at peak gamma pow
    mpfampl1(i) = ftppc(i,2,5,maxpowfreq(i));
    mpfsig(i) = ftralp(i,1,5,maxpowfreq(i))<0.05;
    mpfsigl1(i) = ftralp(i,2,5,maxpowfreq(i))<0.05;
    normftppc(i,:) = ftppc(i,1,5,:)./m;
    normftppcl1(i,:) = ftppc(i,2,5,:)./m;
    lppcsig(i) = ftralp(i,1,5,peakind(i))<=0.05;
    lppcsigl1(i) = ftralp(i,2,5,peakind(i))<=0.05;
end
mpfamp(isinf(mpfamp)) = NaN;
mpfampl1(isinf(mpfampl1)) = NaN;

figure
errorbar(ftfax,squeeze(nanmean(normftppc(l23rs,:))),squeeze(nanstd(normftppc(l23rs,:)))./sqrt(length(find(~isnan(normftppc(l23rs,1))))));
hold on
errorbar(ftfax,squeeze(nanmean(normftppc(l23fs,:))),squeeze(nanstd(normftppc(l23fs,:)))./sqrt(length(find(~isnan(normftppc(l23fs,1))))),'g')
errorbar(ftfax,squeeze(nanmean(normftppc(phe,:))),squeeze(nanstd(normftppc(phe,:)))./sqrt(length(find(~isnan(normftppc(phe,1))))),'r')
xlabel('frequency [Hz]')
ylabel('normalized PPC')
legend('RS','FS','SOM')
title('ALL cells')

okn = nr1spikes(:,1,5)>10;
locked = lppcsig;
n = sum(~isnan(ftppc(l23rs&okn'&locked,1,5)));
figure
hold on
for i = find(l23rs&okn'&locked)
    plot(ftppc(i,1,5,peakind(i)),ftppc(i,2,5,peakind(i)),'ko','markerfacecolor','k')
end
for i = find(l23fs&okn'&locked)
    plot(ftppc(i,1,5,peakind(i)),ftppc(i,2,5,peakind(i)),'go','markerfacecolor','g')
end
for i = find(phe)
    plot(ftppc(i,1,5,peakind(i)),ftppc(i,2,5,peakind(i)),'ro','markerfacecolor','r')
end

% running blank screen ppc - l and c are 2 and 5 always - superfluous
% dimensions - fix it or only use 2,5 of r1contl0ftppc and r1contl1ftppc
figure
errorbar(ftfax,squeeze(nanmean(r1contl0ftppc(l23rs,2,5,:))),squeeze(nanstd(r1contl0ftppc(l23rs,2,5,:)))./sqrt(length(find(~isnan(r1contl0ftppc(l23rs,2,5,1))))));
hold on
errorbar(ftfax,squeeze(nanmean(r1contl0ftppc(l23fs,2,5,:))),squeeze(nanstd(r1contl0ftppc(l23fs,2,5,:)))./sqrt(length(find(~isnan(r1contl0ftppc(l23fs,2,5,1))))),'g');
errorbar(ftfax,squeeze(nanmean(r1contl0ftppc(phe,2,5,:))),squeeze(nanstd(r1contl0ftppc(phe,2,5,:)))./sqrt(length(find(~isnan(ftppc(phe,1,5,1))))),'r')
xlabel('frequency [Hz]')
ylabel('PPC')
legend('RS','FS','SOM')
title('ALL cells')

figure
errorbar(ftfax,squeeze(nanmean(ftplv(l23rs,1,5,:))),squeeze(nanstd(ftplv(l23rs,1,5,:)))./sqrt(length(find(~isnan(ftplv(l23rs,1,5,1))))));
hold on
errorbar(ftfax,squeeze(nanmean(ftplv(l23fs,1,5,:))),squeeze(nanstd(ftplv(l23fs,1,5,:)))./sqrt(length(find(~isnan(ftplv(l23fs,1,5,1))))),'g')
errorbar(ftfax,squeeze(nanmean(ftplv(phe,1,5,:))),squeeze(nanstd(ftplv(phe,1,5,:)))./sqrt(length(find(~isnan(ftplv(phe,1,5,1))))),'r')
xlabel('frequency [Hz]')
ylabel('PPC')
legend('RS','FS','SOM')
title('ALL cells')

% find significant ones for each frequency, take mean of those
for i = 1:length(ftfax)
    ppcl23(i) = nanmean(ftppc(l23rs'&ftralp(:,1,5,i)<0.05,1,5,i));
    ppcl23err(i) = nanstd(ftppc(l23rs'&ftralp(:,1,5,i)<0.05,1,5,i))./sqrt(length(find(l23rs'&ftralp(:,1,5,i)<0.05)));
    ppcl23fs(i) = nanmean(ftppc(l23fs'&ftralp(:,1,5,i)<0.05,1,5,i));
    ppcl23fserr(i) = nanstd(ftppc(l23fs'&ftralp(:,1,5,i)<0.05,1,5,i))./sqrt(length(find(l23fs'&ftralp(:,1,5,i)<0.05)));
%     ppcl23phe(i) = nanmean(ftppc(phe'&ftralp(:,1,5,i)<0.05,1,5,i));
%     ppcl23pheerr(i) = nanstd(ftppc(phe'&ftralp(:,1,5,i)<0.05,1,5,i))./sqrt(length(find(phe'&ftralp(:,1,5,i)<0.05)));
end

figure
errorbar(ftfax,ppcl23,ppcl23err);
hold on
errorbar(ftfax,ppcl23fs,ppcl23fserr,'g');
% errorbar(ftfax,ppcl23phe,ppcl23pheerr,'r');
xlabel('frequency [Hz]')
ylabel('PPC')
legend('RS','FS','SOM')
title('only significant at that frequency')

% trial by trial correlation of firing rate with 20 frequency bands
% GLM of gamm power and firing rates for each neuron
x = 1:1000;
for i = 1:length(depth)
    
% %     [glmparams(i,:), dsa, stats] = glmfit(crudegpow(i,:),firingrates(i,:),'poisson');
% %     rsq(i) = 1-(sum(stats.resid.^2)/sum((firingrates(i,:)-mean(firingrates(i,:)).^2)));
% 
%     pmod = fitglm(crudegpow(i,:),firingrates(i,:),'Distribution','poisson');
%     rsq(i) = pmod.Rsquared.Ordinary;
%     glmcoeffs(i,:) = pmod.Coefficients.Estimate;
%     fitline(i,:) = x.*glmcoeffs(i,2)+glmcoeffs(i,1);
%     
% %     [polyparams(i,:),S] = polyfit(crudegpow(i,:),firingrates(i,:),2);
% %     polyfun(i,:) = x.^2.*polyparams(i,1)+x.*polyparams(i,2)+polyparams(i,3);    
%     
%     [polyparams(i,:),S] = polyfit(crudegpow(i,:),firingrates(i,:),1);
%     polyfun(i,:) = x.*polyparams(i,1)+polyparams(i,2);
%     
    for j = 1:20
        [sprho(i,j), sprp(i,j)] = corr(squeeze(allpowers(i,j,:)),firingrates(i,:)','type','Spearman');
    end
end

figure
plot(x,mean(fitline(l23rs,:)),'b','linewidth',2');
hold on
plot(x,mean(fitline(l23fs,:)),'g','linewidth',2');
plot(x,mean(fitline(phe,:)),'r','linewidth',2');
plot(x,mean(fitline(l23rs,:))+std(fitline(l23rs,:))./sqrt(length(find(l23rs))),'b--')
plot(x,mean(fitline(l23rs,:))-std(fitline(l23rs,:))./sqrt(length(find(l23rs))),'b--')
plot(x,mean(fitline(l23fs,:))+std(fitline(l23fs,:))./sqrt(length(find(l23fs))),'g--')
plot(x,mean(fitline(l23fs,:))-std(fitline(l23fs,:))./sqrt(length(find(l23fs))),'g--')
plot(x,mean(fitline(phe,:))+std(fitline(phe,:))./sqrt(length(find(phe))),'r--')
plot(x,mean(fitline(phe,:))-std(fitline(phe,:))./sqrt(length(find(phe))),'r--')
% errorbar(x,mean(fitline(l23rs,:)),std(fitline(l23rs,:))./sqrt(length(find(l23rs))));
% hold on
% errorbar(x,mean(fitline(l23fs,:)),std(fitline(l23fs,:))./sqrt(length(find(l23fs))),'g');
% errorbar(x,mean(fitline(phe,:)),std(fitline(phe,:))./sqrt(length(find(phe))),'r');

figure
errorbar(3:5:98,nanmean(sprho(l23rs,:)),nanstd(sprho(l23rs,:))./sqrt(length(find(l23rs))));
hold on
errorbar(3:5:98,nanmean(sprho(l23fs,:)),nanstd(sprho(l23fs,:))./sqrt(length(find(l23fs))),'g');
errorbar(3:5:98,nanmean(sprho(phe,:)),nanstd(sprho(phe,:))./sqrt(length(find(phe))),'r');
legend('RS','FS','SOM')

fillxcorr = [3:5:98, fliplr(3:5:98)];
fillysom = [(nanmean(sprho(phe,:))+nanstd(sprho(phe,:))./sqrt(length(find(phe)))) fliplr((nanmean(sprho(phe,:))-nanstd(sprho(phe,:))./sqrt(length(find(phe)))))];
fillyfs = [(nanmean(sprho(l23fs,:))+nanstd(sprho(l23fs,:))./sqrt(length(find(l23fs)))) fliplr((nanmean(sprho(l23fs,:))-nanstd(sprho(l23fs,:))./sqrt(length(find(l23fs)))))];

figure
fill(fillxcorr,fillysom,[0,0,0])
hold on
fill(fillxcorr,fillyfs,[0.5,0.5,0.5])
axis([3,98,-.3,.4])
legend('SOM cells','FS cells')
xlabel('frequency (Hz)')
ylabel('Spearman rho')
title('SOM n = 7, FS n = 33')


%average firing rate and increases/reductions
r1nlfr = nanmean(nanmean(r1condfr(:,1,:,:),3),4);
r1lfr = nanmean(nanmean(r1condfr(:,2,:,:),3),4);
% firing rate for largest
l0s5fr = squeeze(nanmean(r1condfr(:,1,:,5),3));
l1s5fr = squeeze(nanmean(r1condfr(:,2,:,5),3));
% OMIs
r1omi = (r1lfr-r1nlfr)./(r1lfr+r1nlfr);
r1s5omi = (l1s5fr-l0s5fr)./(l1s5fr+l0s5fr);

nhalo= length(find(phe));
nl23rs = length(find(l23rs));
nl23fs = length(find(l23fs));

redhalo = (1-(r1lfr(phe)./r1nlfr(phe))).*100;
pci = r1lfr./r1nlfr; pci(isinf(pci)) = NaN; pci = (pci-1).*100; % percent increase
incrs = ((r1lfr(l23rs)./r1nlfr(l23rs))-1).*100; incrs(isinf(incrs)) = NaN;
incfs = ((r1lfr(l23fs&~phe)./r1nlfr(l23fs&~phe))-1).*100;
dlt = l1s5fr-l0s5fr;
pcis5 = l1s5fr./l0s5fr; pcis5(isinf(pcis5)) = NaN; 
disp(['size 5: ' int2str(nhalo) ' Halo cells reduce by ' num2str((1-nanmean(pcis5(phe))).*100) '+-' num2str(nanstd(pcis5(phe).*100)./sqrt(nhalo))]);
pcis5 = (pcis5-1).*100;
disp(['size 5: ' int2str(nl23rs) ' L2/3RS cells increase by ' num2str(nanmean(pcis5(l23rs))) '+-' num2str(nanstd(pcis5(l23rs))./sqrt(nhalo))]);
disp(['size 5: ' int2str(nl23fs) ' L2/3FS cells increase by ' num2str(nanmean(pcis5(l23fs))) '+-' num2str(nanstd(pcis5(l23fs))./sqrt(nhalo))]);

%combined pop only
if length(depth)>300
    som = zeros(1,length(depth));
    pv = zeros(1,length(depth));
    som(1:251) = 1;
    pv(252:end) = 1;
    
    bars = [nanmean(pci(l23rs&~phe&som)),nanmean(pci(l23fs&~phe&som));...
        nanmean(pci(l23rs&~phe&pv)),nanmean(pci(l23fs&~phe&pv))];
    errorbars = [nanstd(pci(l23rs&~phe&som))./sqrt(length(find(~isnan(pci(l23rs&~phe&som))))),...
        nanstd(pci(l23fs&~phe&som))./sqrt(length(find(~isnan(pci(l23fs&~phe&som)))));...
        nanstd(pci(l23rs&~phe&pv))./sqrt(length(find(~isnan(pci(l23rs&~phe&pv))))),...
        nanstd(pci(l23fs&~phe&pv))./sqrt(length(find(~isnan(pci(l23fs&~phe&pv)))))];
    figure
    barweb(bars,errorbars,[],{'SOM Halo','PV Halo'},'firing rate percent changes',[],'percent change',[],[],{'RS','FS'});
    
        
    bars = [nanmean(dlt(l23rs&~phe&som)),nanmean(dlt(l23fs&~phe&som));...
        nanmean(dlt(l23rs&~phe&pv)),nanmean(dlt(l23fs&~phe&pv))];
    errorbars = [nanstd(dlt(l23rs&~phe&som))./sqrt(length(find(~isnan(dlt(l23rs&~phe&som))))),...
        nanstd(dlt(l23fs&~phe&som))./sqrt(length(find(~isnan(dlt(l23fs&~phe&som)))));...
        nanstd(dlt(l23rs&~phe&pv))./sqrt(length(find(~isnan(dlt(l23rs&~phe&pv))))),...
        nanstd(dlt(l23fs&~phe&pv))./sqrt(length(find(~isnan(dlt(l23fs&~phe&pv)))))];
    figure
    barweb(bars,errorbars,[],{'SOM Halo','PV Halo'},'spikes added',[],'# of spikes',[],[],{'RS','FS'});
    

    % firing rate plot spread
    figure
    plotSpread({pci(l23rs&pv&~phe),pci(l23fs&pv&~phe),pci(pv&phe),pci(l23rs&som&~phe),pci(l23fs&som&~phe),pci(som&phe)},[],[],{'PV RS','PV FS','PV Halo','SOM RS','SOM FS','SOM Halo'},0)
    % plotSpread({pci(l23rs&pv&~phe),pci(l23fs&pv&~phe),pci(l23rs&som&~phe),pci(l23fs&som&~phe)},[],[],{'PV RS','PV FS','SOM RS','SOM FS'},0)
    ylabel('percent change')
    hold on
    line([.8,1.2],[nanmedian(pci(l23rs&pv&~phe)),nanmedian(pci(l23rs&pv&~phe))],'color','r','linewidth',2)
    line([1.8,2.2],[nanmedian(pci(l23fs&pv&~phe)),nanmedian(pci(l23fs&pv&~phe))],'color','r','linewidth',2)
    line([2.8,3.2],[nanmedian(pci(pv&phe)),nanmedian(pci(pv&phe))],'color','r','linewidth',2)
    line([3.8,4.2],[nanmedian(pci(l23rs&som&~phe)),nanmedian(pci(l23rs&som&~phe))],'color','r','linewidth',2)
    line([4.8,5.2],[nanmedian(pci(l23fs&som&~phe)),nanmedian(pci(l23fs&som&~phe))],'color','r','linewidth',2)
    line([5.8,6.2],[nanmedian(pci(som&phe)),nanmedian(pci(som&phe))],'color','r','linewidth',2)

    [rsomi,x] = hist(r1s5omi(l23rs&som),-.9:.2:.9);
    [haloomi,x] = hist(r1s5omi(phe&som),-.9:.2:.9);
    figure
    bar(x,rsomi,1,'k')
    hold on
    bar(x,haloomi,1,'r')
    xlabel('OMI')
    ylabel('cell count')
    legend('L2/3 RS','putative SOM Halo cells','location','nw')
    
    [rsomi,x] = hist(r1s5omi(l23rs&pv),-.9:.2:.9);
    [haloomi,x] = hist(r1s5omi(phe&pv),-.9:.2:.9);
    figure
    bar(x,rsomi,1,'k')
    hold on
    bar(x,haloomi,1,'r')
    xlabel('OMI')
    ylabel('cell count')
    legend('L2/3 RS','putative PV Halo cells','location','nw')
    
    
end


% firing rate scatterplot
figure
plot(r1nlfr(l23rs),r1lfr(l23rs),'ko','markerfacecolor','k')
hold on
plot(r1nlfr(l23fs),r1lfr(l23fs),'go','markerfacecolor','g')
plot(r1nlfr(phe),r1lfr(phe),'ro','markerfacecolor','r')
% axis([-.5,22,-.5,22])
axis([-.5,35,-.5,35])
legend('RS cells', 'FS cells', 'putative Halo cells','location','nw')
axis square
% line([-.5,22],[-.5,22],'color','k')
line([-.5,35],[-.5,35],'color','k')
xlabel('firing rate control [Hz]')
% ylabel('firing rate PV Halo [Hz]')
ylabel('firing rate SOM Halo [Hz]')

[yrs,x] = hist(r1omi(l23rs),-.9:.2:.9);
[yfs,x] = hist(r1omi(l23fs),-.9:.2:.9);
[yhe,x] = hist(r1omi(phe),-.9:.2:.9);
[ylrs,x] = hist(r1s5omi(l23rs),-.9:.2:.9);
[ylfs,x] = hist(r1s5omi(l23fs),-.9:.2:.9);
[ylhe,x] = hist(r1s5omi(phe),-.9:.2:.9);
x = x-.1;
stairs(x,ylrs,'b')
hold on
stairs(x+.02,ylfs,'g')
stairs(x+.04,ylhe,'r')


% phase locking stats
okn = nr1spikes(:,1,5)>10; % more than 10 spikes to calculate reighly and ppc
nokl23rs = sum(l23rs&okn');
nokl23fs = sum(l23fs&okn');
nokhalo = sum(phe&okn');
lockedl23rs = length(find(r1lockpval(l23rs&okn',1,5)<0.05));
lockedl23fs = length(find(r1lockpval(l23fs&okn',1,5)<0.05));
lockedhalo = length(find(r1lockpval(phe&okn',1,5)<0.05));
disp([int2str(lockedl23rs) ' of ' int2str(nokl23rs) '(' num2str((lockedl23rs/nokl23rs)*100)  '%) L2/3 RS cells are significantly locked during running'])
disp([int2str(lockedl23fs) ' of ' int2str(nokl23fs) '(' num2str((lockedl23fs/nokl23fs)*100)  '%) L2/3 FS cells are significantly locked during running'])
disp([int2str(lockedhalo) ' of ' int2str(nokhalo) '(' num2str((lockedhalo/nokhalo)*100)  '%) Halo cells are significantly locked during running'])
% for PV also took the significantly locked cell with only 5 spikes so
% dropped the <10 criterion
disp(['average ppc at gamma peak frequency for halo cells: ' num2str(nanmean(r1orimeanppc(phe,1,5))) ...
    ' +- ' num2str(nanstd(r1orimeanppc(phe,1,5))./sqrt(nokhalo))])

locked = r1lockpval(:,1,5)<0.05;
[p,s,stats] = signrank(r1orimeanppc(locked&l23rs'&okn,1,5),r1orimeanppc(locked&l23rs'&okn,2,5))
[p,s,stats] = signrank(r1orimeanppc(locked&l23fs'&okn,1,5),r1orimeanppc(locked&l23fs'&okn,2,5))

% phase locking for high gamma during small stim
okns = nr1spikes(:,1,1)>10;
nokl23rss = sum(l23rs&okns');
nokl23fss = sum(l23fs&okns');
lockedl23rss = length(find(r1g2lockpval(l23rs&okns',1,1)<0.05));
lockedl23fss = length(find(r1g2lockpval(l23fs&okns',1,1)<0.05));

% phase locking for high gamma during control
oknc = nspikesr1ctrl(:,1)>10;
nokl23rsc = sum(l23rs&oknc');
nokl23fsc = sum(l23fs&oknc');
lockedl23rs = length(find(r1g2lockpvalctrl(l23rs&oknc',1)<0.05));
lockedl23fs = length(find(r1g2lockpvalctrl(l23fs&oknc',1)<0.05));

spacing = (2*pi)/16;
bincenters = spacing/2:spacing:(16*spacing)-spacing/2;
% bincenters = bincenters - pi;

examplecell = 170;
a = find(l23rs'&nr1spikes(:,1,5)>30);
for cll = 1:length(a)
    figure
    subplot(2,2,1)
    rose(r1orimeanphases{a(cll),1,5},16)
    subplot(2,2,2)
    rose(r1orimeanphases{a(cll),2,5},16)
    title(int2str(a(cll)))
    
    neleml0 = hist(uncircle(r1orimeanphases{a(cll),1,5}),bincenters);
    neleml1 = hist(uncircle(r1orimeanphases{a(cll),2,5}),bincenters);
    subplot(2,2,3)
    plot(bincenters,neleml0)
    hold on
    plot(bincenters,neleml1,'r')
    set(gca,'xtick',[0,pi/2,pi,3*pi/2])
    set(gca,'xticklabel',{'0','pi/2','pi','3pi/2'})
end

for cll = 1:length(depth)
    for i = 1:2
        for j = 1:5
            [s,cstd(cll,i,j)] = circ_std(r1orimeanphases{cll,i,j});
%             peakgammasfc(cll,i,j) = squeeze(C(cll,i,j,lgi(cll)));
            peakgammapow(cll,i,j) = squeeze(condS(cll,i,j,lgi(cll)));
            intgammapow(cll,i,j) = squeeze(nanmean(condS(cll,i,j,lgi(cll)-5:lgi(cll)+5),4));
            r1peakgammapow(cll,i,j) = squeeze(condrunS(cll,i,j,lgi(cll)));
            r0peakgammapow(cll,i,j) = squeeze(condstillS(cll,i,j,lgi(cll)));
            r1intgammapow(cll,i,j) = squeeze(nanmean(condrunS(cll,i,j,lgi(cll)-5:lgi(cll)+5),4));
            r1peakhighgammapow(cll,i,j) = squeeze(condrunS(cll,i,j,hgi(cll)));
            r0peakhighgammapow(cll,i,j) = squeeze(condstillS(cll,i,j,hgi(cll)));
%             sighighersfclarge(cll,i) = Cerr(cll,i,5,1,lgi(cll))>Cerr(cll,i,1,2,lgi(cll));
%             filly(cll,i,j,:) = [squeeze(Cerr(cll,i,j,1,1:104))',fliplr(squeeze(Cerr(cll,i,j,2,1:104))')];
            fillspecy(cll,i,j,:) = [squeeze(condSerr(cll,i,j,1,1:104))',fliplr(squeeze(condSerr(cll,i,j,2,1:104))')];
            blserrs(cll,i,j,:,2) = squeeze(condSerr(cll,i,j,2,:))-squeeze(condS(cll,i,j,:));
            blserrs(cll,i,j,:,1) = squeeze(condS(cll,i,j,:))-squeeze(condSerr(cll,i,j,1,:));
            r0fillspecy(cll,i,j,:) = [squeeze(condstillSerr(cll,i,j,1,12:103))',fliplr(squeeze(condstillSerr(cll,i,j,2,12:103))')];
            r1fillspecy(cll,i,j,:) = [squeeze(condrunSerr(cll,i,j,1,12:103))',fliplr(squeeze(condrunSerr(cll,i,j,2,12:103))')];
            r1tbtfillspecy(cll,i,j,:) = [(squeeze(r1condallS(cll,i,j,12:103))+squeeze(r1condallSerr(cll,i,j,12:103)))',fliplr((squeeze(r1condallS(cll,i,j,12:103))-squeeze(r1condallSerr(cll,i,j,12:103)))')];
            r0tbtfillspecy(cll,i,j,:) = [(squeeze(r0condallS(cll,i,j,12:103))+squeeze(r0condallSerr(cll,i,j,12:103)))',fliplr((squeeze(r0condallS(cll,i,j,12:103))-squeeze(r0condallSerr(cll,i,j,12:103)))')];
%             r1evokedspecy(cll,i,j,:) = [(squeeze(r1evokedspect(cll,i,j,:))+squeeze(r1evokedspecterr(cll,i,j,1,:)))',fliplr((squeeze(r1evokedspect(cll,i,j,:))-squeeze(r1evokedspecterr(cll,i,j,:)))')];
%             r0evokedspecy(cll,i,j,:) = [(squeeze(r0evokedspect(cll,i,j,:))+squeeze(r0evokedspecterr(cll,i,j,:)))',fliplr((squeeze(r0evokedspect(cll,i,j,:))-squeeze(r0evokedspecterr(cll,i,j,:)))')];
            if i == 1
                r1c0tbtfillspecy(cll,i,:) = [(squeeze(r1l0contallS(cll,12:103))+squeeze(r1l0contallSerr(cll,12:103))),fliplr((squeeze(r1l0contallS(cll,12:103))-squeeze(r1l0contallSerr(cll,12:103))))];
                r0c0tbtfillspecy(cll,i,:) = [(squeeze(r0l0contallS(cll,12:103))+squeeze(r0l0contallSerr(cll,12:103))),fliplr((squeeze(r0l0contallS(cll,12:103))-squeeze(r0l0contallSerr(cll,12:103))))];
                r1c0peakhighgammapow(cll,i,:) = squeeze(r1l0contS(cll,hgi(cll)));
                r0c0peakhighgammapow(cll,i,:) = squeeze(r0l0contS(cll,hgi(cll)));
                r1c0peakgammapow(cll,i,:) = squeeze(r1l0contS(cll,lgi(cll)));
                r0c0peakgammapow(cll,i,:) = squeeze(r0l0contS(cll,lgi(cll)));
            else                
                r1c0tbtfillspecy(cll,i,:) = [(squeeze(r1l1contallS(cll,12:103))+squeeze(r1l1contallSerr(cll,12:103))),fliplr((squeeze(r1l1contallS(cll,12:103))-squeeze(r1l1contallSerr(cll,12:103))))];
                r0c0tbtfillspecy(cll,i,:) = [(squeeze(r0l1contallS(cll,12:103))+squeeze(r0l1contallSerr(cll,12:103))),fliplr((squeeze(r0l1contallS(cll,12:103))-squeeze(r0l1contallSerr(cll,12:103))))];
                r1c0peakhighgammapow(cll,i,:) = squeeze(r1l1contS(cll,hgi(cll)));
                r0c0peakhighgammapow(cll,i,:) = squeeze(r0l1contS(cll,hgi(cll)));
                r1c0peakgammapow(cll,i,:) = squeeze(r1l1contS(cll,lgi(cll)));
                r0c0peakgammapow(cll,i,:) = squeeze(r0l1contS(cll,lgi(cll)));
            end
            relgammapow(cll,i,j) = condrunS(cll,i,j,lgi(cll))./nanmean(condrunS(cll,i,j,12:104));
            relgammapowr0(cll,i,j) = condstillS(cll,i,j,lgi(cll))./nanmean(condstillS(cll,i,j,12:104));
            relhgpow(cll,i,j) = condrunS(cll,i,j,hgi(cll))./nanmean(condrunS(cll,i,j,12:104));
        end
    end
%     sighighersfcnolight(cll) = Cerr(cll,1,5,1,lgi(cll))>Cerr(cll,2,5,2,lgi(cll));
    fillycontSr0l0(cll,:) = [squeeze(r0l0contSerr(cll,1,1:104))',fliplr(squeeze(r0l0contSerr(cll,2,1:104))')];
    fillycontSr1l0(cll,:) = [squeeze(r1l0contSerr(cll,1,1:104))',fliplr(squeeze(r1l0contSerr(cll,2,1:104))')];
    fillycontSr0l1(cll,:) = [squeeze(r0l1contSerr(cll,1,1:104))',fliplr(squeeze(r0l1contSerr(cll,2,1:104))')];
    fillycontSr1l1(cll,:) = [squeeze(r1l1contSerr(cll,1,1:104))',fliplr(squeeze(r1l1contSerr(cll,2,1:104))')];
end
fillx = [chf(12:103),fliplr(chf(12:103))];
efillx = [chf(1:103),fliplr(chf(1:103))];

% Halo vs eArch effect - for PV Halo eArch mix pop
pcgc = r1peakgammapow(lfpinds,2,5)./r1peakgammapow(lfpinds,1,5); % percent gamma change
p = ranksum(pcgc(1:12),pcgc(13:19));
disp(['percent change Halo: ' num2str(mean(pcgc(1:12))) '+-' num2str(std(pcgc(1:12))./sqrt(12)) ...
    '    percent change eArch: ' num2str(mean(pcgc(13:19))) '+-' num2str(std(pcgc(13:19))./sqrt(7)) ...
    '    p = ' num2str(p)]);

% evoked spectra with "math" error bars
for i = 1:length(lfpinds)
    figure
    plot(chf(12:103),squeeze(r1evokedspect(lfpinds(i),1,5,12:103)),'k','linewidth',2)
    hold on
    plot(chf(12:103),squeeze(r1evokedspect(lfpinds(i),2,5,12:103)),'r','linewidth',2)
%     fill(efillx,squeeze(r1evokedspecy(lfpinds(i),1,5,:)),'k')
%     hold on
%     fill(efillx,squeeze(r1evokedspecy(lfpinds(i),2,5,:)),'r')
end

% percent gamma reduction for running not runnign
a = r1peakgammapow(lfpinds(1:16),2,5)./r1peakgammapow(lfpinds(1:16),1,5);
b = r0peakgammapow(lfpinds(1:16),2,5)./r0peakgammapow(lfpinds(1:16),1,5);
erra = nanstd(a)./sqrt(length(find(~isnan(a))));
errb = nanstd(b)./sqrt(length(find(~isnan(b))));

figure
[p,s] = signrank(a,b);
plot(1,a,'go','markerfacecolor','g')
hold on
plot(2,b,'ko','markerfacecolor','k')
for i = 1:16
    plot([1,2],[a(i),b(i)],'k')
end
axis([0,3,0,1])
set(gca,'xticklabel',{'running','quiescent'})
set(gca,'xtick',[1,2]);
set(gca,'ytick',[0:.25:1]);
ylabel('percent peak gamma reduction')
title(['running vs quiescence p: ' num2str(p)]);
set(gcf,'OuterPosition',[573   504   242   513])

% running on chnages in gamma SOM and PV Halo
a = r1peakgammapow(lfpinds,2,5)./r1peakgammapow(lfpinds,1,5);
b = r0peakgammapow(lfpinds,2,5)./r0peakgammapow(lfpinds,1,5);
figure
plot(a(1:16),b(1:16),'ro','markerfacecolor','r')
hold on
plot(a(17:34),b(17:34),'bo','markerfacecolor','b')
axis square
hold on
plot([0,4],[1,1],'k')
plot([1,1],[0,4],'k')
set(gca,'xtick',[0:4]);
set(gca,'ytick',[0:4]);
legend('SOM Halo','PV Halo')
xlabel('change with light during running')
ylabel('change with light during quiescence')
title('effect of running on gamma changes')

% doesn't work because r0 peaks are not in same places sometimes - take
% indivisual peaks further down
% % running on chnages in relative gamma SOM and PV Halo
% a = relgammapow(lfpinds,2,5)./relgammapow(lfpinds,1,5);
% b = relgammapowr0(lfpinds,2,5)./relgammapowr0(lfpinds,1,5);
% figure
% plot(a(1:16),b(1:16),'ro','markerfacecolor','r')
% hold on
% plot(a(17:34),b(17:34),'bo','markerfacecolor','b')
% axis([0,2.5,0,2.5])
% axis square
% hold on
% plot([0,2.5],[1,1],'k')
% plot([1,1],[0,2.5],'k')
% set(gca,'xtick',[0:2]);
% set(gca,'ytick',[0:2]);
% legend('SOM Halo','PV Halo')
% xlabel('change with light during running')
% ylabel('change with light during quiescence')
% title('effect of running on relative gamma changes')

% figure
% p = signrank(relgammapow(lfpinds(17:34),1,5),relgammapow(lfpinds(17:34),2,5));
% subplot(1,2,1)
% plot(1,relgammapow(lfpinds(17:34),1,5),'ko','markerfacecolor','k')
% hold on
% plot(2,relgammapow(lfpinds(17:34),2,5),'ro','markerfacecolor','r')
% plot([1,2],[relgammapow(lfpinds(17:34),1,5),relgammapow(lfpinds(17:34),2,5)], 'k')
% axis([0,3,1,8])
% set(gca,'xtick',[1,2]);
% set(gca,'xticklabel',{'control','PV Halo'})
% ylabel('relative gamma power')
% title(['running p: ' num2str(p)])
% 
% subplot(1,2,2)
% p = signrank(relgammapowr0(lfpinds(17:34),1,5),relgammapowr0(lfpinds(17:34),2,5));
% plot(1,relgammapowr0(lfpinds(17:34),1,5),'ko','markerfacecolor','k')
% hold on
% plot(2,relgammapowr0(lfpinds(17:34),2,5),'ro','markerfacecolor','r')
% for i = 17:34
%     plot([1,2],[relgammapowr0(lfpinds(17:34),1,5),relgammapowr0(lfpinds(17:34),2,5)], 'k')
% end
% axis([0,3,1,8])
% set(gca,'xtick',[1,2]);
% set(gca,'xticklabel',{'control','PV Halo'})
% ylabel('relative gamma power')
% title(['quiescent p: ' num2str(p)])

% SOM Halo on high gamma relative
figure
[p,s] = signrank(relhgpow(lfpinds(1:16),1,1),relhgpow(lfpinds(1:16),2,1));
plot(1,relhgpow(lfpinds(1:16),1,1),'ko','markerfacecolor','k')
hold on
plot(2,relhgpow(lfpinds(1:16),2,1),'ro','markerfacecolor','r')
for i = 1:16
    plot([1,2],[relhgpow(lfpinds(i),1,1),relhgpow(lfpinds(i),2,1)],'k')
end
axis([0,3,0,3])
set(gca,'xticklabel',{'control','light'})
set(gca,'xtick',[1,2]);
set(gca,'ytick',[0:3]);
ylabel('relative high gamma power')
title(['SOM Halo p: ' num2str(p)]);
set(gcf,'OuterPosition',[573   504   242   513])

% SOM Halo on high gamma absolute
figure
[p,s] = signrank(r1peakhighgammapow(lfpinds(1:16),1,1),r1peakhighgammapow(lfpinds(1:16),2,1));
plot(1,r1peakhighgammapow(lfpinds(1:16),1,1),'ko','markerfacecolor','k')
hold on
plot(2,r1peakhighgammapow(lfpinds(1:16),2,1),'ro','markerfacecolor','r')
for i = 1:16
    plot([1,2],[r1peakhighgammapow(lfpinds(i),1,1),r1peakhighgammapow(lfpinds(i),2,1)],'k')
end
axis([0,3,0,100])
set(gca,'xticklabel',{'control','light'})
set(gca,'xtick',[1,2]);
set(gca,'ytick',[0:20:100]);
ylabel('absolute high gamma power')
title(['SOM Halo p: ' num2str(p)]);
set(gcf,'OuterPosition',[573   504   242   513])


% PV Halo+eArch on high gamma relative
figure
[p,s] = signrank(relhgpow(lfpinds(17:34),1,1),relhgpow(lfpinds(17:34),2,1));
plot(1,relhgpow(lfpinds(17:34),1,1),'ko','markerfacecolor','k')
hold on
plot(2,relhgpow(lfpinds(17:34),2,1),'ro','markerfacecolor','r')
for i = 17:34
    plot([1,2],[relhgpow(lfpinds(i),1,1),relhgpow(lfpinds(i),2,1)],'k')
end
axis([0,3,0,3])
set(gca,'xticklabel',{'control','light'})
set(gca,'xtick',[1,2]);
set(gca,'ytick',[0:3]);
ylabel('relative high gamma power')
title(['PV Halo+eArch p: ' num2str(p)]);
set(gcf,'OuterPosition',[573   504   242   513])

% PV Halo+eArch on high gamma absolute
figure
[p,s] = signrank(r1peakhighgammapow(lfpinds(17:34),1,1),r1peakhighgammapow(lfpinds(17:34),2,1));
plot(1,r1peakhighgammapow(lfpinds(17:34),1,1),'ko','markerfacecolor','k')
hold on
plot(2,r1peakhighgammapow(lfpinds(17:34),2,1),'ro','markerfacecolor','r')
for i = 17:34
    plot([1,2],[r1peakhighgammapow(lfpinds(i),1,1),r1peakhighgammapow(lfpinds(i),2,1)],'k')
end
axis([0,3,0,100])
set(gca,'xticklabel',{'control','light'})
set(gca,'xtick',[1,2]);
set(gca,'ytick',[0:20:100]);
ylabel('absolute high gamma power')
title(['PV Halo+eArch p: ' num2str(p)]);
set(gcf,'OuterPosition',[573   504   242   513])


% vis to BL ratio for evoked spectra and power
for i = 1:size(condS,1)
    for s = 1:5
        for l = 1:2
            if l == 1
                v2blratior1(i,l,s,:) = squeeze(condrunS(i,l,s,:))./r1l0contS(i,:)';
                v2blratior0(i,l,s,:) = squeeze(condstillS(i,l,s,:))./r0l0contS(i,:)';
            else
                v2blratior1(i,l,s,:) = squeeze(condrunS(i,l,s,:))./r1l1contS(i,:)';
                v2blratior0(i,l,s,:) = squeeze(condstillS(i,l,s,:))./r0l1contS(i,:)';
            end
%             v2blfftratior1(i,l,s,:) = squeeze(r1condfftspect(i,l,s,:))./squeeze(r1contfftspect(i,l,:));
%             v2blfftratior0(i,l,s,:) = squeeze(r0condfftspect(i,l,s,:))./squeeze(r0contfftspect(i,l,:));
        end
        v2nlblratior1(i,1,s,:) = squeeze(condrunS(i,1,s,:))./r1l0contS(i,:)';  % also light on ratio nromalized to no light spontaeous
        v2nlblratior1(i,2,s,:) = squeeze(condrunS(i,2,s,:))./r1l0contS(i,:)';
        v2nlblratior0(i,1,s,:) = squeeze(condstillS(i,1,s,:))./r0l0contS(i,:)';
        v2nlblratior0(i,2,s,:) = squeeze(condstillS(i,2,s,:))./r0l0contS(i,:)';
    end
end

for i = 1:16
% for i = 17:length(lfpinds)
    figure
%     plot(chf(1:150), log10(squeeze(v2blratior1(lfpinds(i),1,5,:))),'k','linewidth',2);
    plot(chf(1:103), smooth(squeeze(v2blratior1(lfpinds(i),1,5,1:103))),'k','linewidth',2);
    hold on
    plot(chf(1:103), smooth(squeeze(v2blratior1(lfpinds(i),2,5,1:103))),'r','linewidth',2)
    plot(chf(1:103), smooth(squeeze(v2nlblratior1(lfpinds(i),2,5,1:103))),'g','linewidth',2)
%     plot([0,150],[0,0],'k')
    plot([0,100],[1,1],'k')
end


% find peak and get peak amplitude and normalized to peak ratios for each animal
for i = 1:34
% for i = 17:34
    for l = 1:2
        [mx,mxind] = max(squeeze(v2blratior1(lfpinds(i),l,5,25:40)));
        evokedpeakindr1(i,l) = mxind+24;
        ratior1amp(i,l) = v2blratior1(lfpinds(i),l,5,mxind+24);
%         ratior1amp(i,2) = v2blratior1(lfpinds(i),2,5,mxind+24);
        normv2blr1(i,l,:) = squeeze(v2blratior1(lfpinds(i),l,5,:))./ratior1amp(i,1);
%         normv2blr1(i,2,:) = squeeze(v2blratior1(lfpinds(i),2,5,:))./mx;
%         normv2nlblr1(i,2,:) = squeeze(v2nlblratior1(lfpinds(i),2,5,:))./mx;
%         v2blratior1atlgpeak(i,l) = squeeze(v2blratior1(lfpinds(i),l,5,lgi(lfpinds(i)))); % just for control - don't use peaks are not the same
%         v2blratior1atlgpeak(i,2) = squeeze(v2blratior1(lfpinds(i),2,5,lgi(lfpinds(i))));

        for s = 1:5
            ratior1ampall(i,s) = v2blratior1(lfpinds(i),1,s,evokedpeakindr1(i,1));
        end    
        normr1rall(i,:) = ratior1ampall(i,:)./max(ratior1ampall(i,:));

        [mx,mxind] = max(squeeze(v2blratior0(lfpinds(i),l,5,25:40)));
        evokedpeakindr0(i,l) = mxind+24;
        ratior0amp(i,l) = v2blratior0(lfpinds(i),l,5,mxind+24);
%         ratior0amp(i,2) = v2blratior0(lfpinds(i),2,5,mxind+24);
        normv2blr0(i,l,:) = squeeze(v2blratior0(lfpinds(i),l,5,:))./ratior0amp(i,1);
%         normv2blr0(i,2,:) = squeeze(v2blratior0(lfpinds(i),2,5,:))./mx;
%         normv2nlblr0(i,2,:) = squeeze(v2nlblratior0(lfpinds(i),2,5,:))./mx;

    end
end

figure
errorbar(sizes,nanmean(normr1rall),nanstd(normr1rall)./sqrt(length(find(~isnan(normr1rall(:,1))))),'ko-','markerfacecolor','k','linewidth',2)
xlabel('stimulus size (degrees)')
ylabel('normalized evoked gamma power')
title('evoked gamma with size')

% evoked low gamma power SOM
figure
[p,s] = signrank(ratior1amp(1:16,1),ratior1amp(1:16,2));
plot(1,ratior1amp(1:16,1),'ko','markerfacecolor','k')
hold on
plot(2,ratior1amp(1:16,2),'ro','markerfacecolor','r')
for i = 1:16
    plot([1,2],[ratior1amp(i,1),ratior1amp(i,2)],'k')
end
axis([0,3,1,11])
set(gca,'xtick',[1,2])
set(gca,'xticklabel',{'control','SOM Halo'})
ylabel('evoked gamma power')
title(['evoked gamma power p: ' num2str(p)]);
set(gcf,'OuterPosition',[573   504   242   513])

% evoked low gamma power PV
figure
[p,s] = signrank(ratior1amp(17:34,1),ratior1amp(17:34,2));
plot(1,ratior1amp(17:34,1),'ko','markerfacecolor','k')
hold on
plot(2,ratior1amp(17:34,2),'ro','markerfacecolor','r')
for i = 17:34
    plot([1,2],[ratior1amp(i,1),ratior1amp(i,2)],'k')
end
axis([0,3,1,12])
set(gca,'xtick',[1,2])
set(gca,'xticklabel',{'control','PV Halo'})
ylabel('evoked gamma power')
title(['evoked gamma power p: ' num2str(p)]);
set(gcf,'OuterPosition',[573   504   242   513])


inds = 17:34;
fillyv2bl(1,:) = [(squeeze(nanmean(normv2blr1(inds,1,1:104),1))+squeeze(nanstd(normv2blr1(inds,1,1:104),1,1))./sqrt(length(find(~isnan(normv2blr1(inds,1,1))))))',...
    fliplr((squeeze(nanmean(normv2blr1(inds,1,1:104),1))-squeeze(nanstd(normv2blr1(inds,1,1:104),1,1))./sqrt(length(find(~isnan(normv2blr1(inds,1,1))))))')];
fillyv2bl(2,:) = [(squeeze(nanmean(normv2blr1(inds,2,1:104),1))+squeeze(nanstd(normv2blr1(inds,2,1:104),1,1))./sqrt(length(find(~isnan(normv2blr1(inds,2,1))))))',...
    fliplr((squeeze(nanmean(normv2blr1(inds,2,1:104),1))-squeeze(nanstd(normv2blr1(inds,2,1:104),1,1))./sqrt(length(find(~isnan(normv2blr1(inds,2,1))))))')];
fillyv2nlbl(2,:) = [(squeeze(nanmean(normv2nlblr1(inds,2,1:104),1))+squeeze(nanstd(normv2nlblr1(inds,2,1:104),1,1))./sqrt(length(find(~isnan(normv2nlblr1(inds,2,1))))))',...
    fliplr((squeeze(nanmean(normv2nlblr1(inds,2,1:104),1))-squeeze(nanstd(normv2nlblr1(inds,2,1:104),1,1))./sqrt(length(find(~isnan(normv2nlblr1(inds,2,1))))))')];
figure
fill(fillx,fillyv2bl(2,:),[1,.3,.3])
hold on
fill(fillx,fillyv2bl(1,:),'k')
axis([10,100,0,1])
legend('SOM Halo','control')
xlabel('frequency (Hz)')
ylabel('normalized power ratio visual/spontaneous')

% errorbar(chf(1:103),squeeze(nanmean(normv2blr1(:,1,1:103),1)),squeeze(nanstd(normv2blr1(:,1,1:103),1,1))./sqrt(length(find(~isnan(normv2blr1(:,1,1))))))
% hold on
% errorbar(chf(1:103),squeeze(nanmean(normv2blr1(:,2,1:103),1)),squeeze(nanstd(normv2blr1(:,2,1:103),1,1))./sqrt(length(find(~isnan(normv2blr1(:,2,1))))),'r')

figure
[p,s] = signrank(ratior1amp(:,1),ratior1amp(:,2));
plot(1,ratior1amp(:,1),'ko','markerfacecolor','k')
hold on
plot(2,ratior1amp(:,2),'ro','markerfacecolor','r')
for i = 1:16
    plot([1,2],[ratior1amp(:,1),ratior1amp(:,2)],'k')
end
axis([0,3,0,11])
set(gca,'xticklabel',{'control','SOM Halo'})
set(gca,'xtick',[1,2]);
set(gca,'ytick',[1:11]);
ylabel('gamma peak / spontaneous')
title(['vis to spont ratio. p: ' num2str(p)]);
set(gcf,'OuterPosition',[573   504   242   513])

figure
%     plot(chf(1:150), log10(squeeze(v2blratior1(lfpinds(i),1,5,:))),'k','linewidth',2);
plot(chf(1:103), smooth(squeeze(v2blratior1(lfpinds(i),1,5,1:103))),'k','linewidth',2);
hold on
plot(chf(1:103), smooth(squeeze(v2blratior1(lfpinds(i),2,5,1:103))),'r','linewidth',2)
plot(chf(1:103), smooth(squeeze(v2nlblratior1(lfpinds(i),2,5,1:103))),'g','linewidth',2)
%     plot([0,150],[0,0],'k')
plot([0,100],[1,1],'k')

for i = 1:length(lfpinds)
    figure
%     plot(chf(1:150), log10(squeeze(v2blratior1(lfpinds(i),1,5,:))),'k','linewidth',2);
    plot(log10(squeeze(v2blratior1(lfpinds(i),1,5,1:103))),'k','linewidth',2);
    hold on
    plot(log10(squeeze(v2blratior1(lfpinds(i),2,5,1:103))),'r','linewidth',2)
%     plot(chf(1:103), log10(squeeze(v2nlblratior1(lfpinds(i),2,5,1:103))),'g','linewidth',2)
%     plot([0,150],[0,0],'k')
    plot([0,100],[0,0],'k')
end

% trial by trial error of ratio to mean baseline
for i = 1:length(lfpinds)
    figure
    plot(chf(1:103),squeeze(condindratio(lfpinds(i),1,5,1:103)),'b','linewidth',2)
    hold on
    plot(chf(1:103),squeeze(condindratio(lfpinds(i),2,5,1:103)),'r','linewidth',2)
    plot(chf(1:103),squeeze(condindratio(lfpinds(i),1,5,1:103))+squeeze(condindratioerr(lfpinds(i),1,5,1:103)),'b')
    plot(chf(1:103),squeeze(condindratio(lfpinds(i),1,5,1:103))-squeeze(condindratioerr(lfpinds(i),1,5,1:103)),'b')
    plot(chf(1:103),squeeze(condindratio(lfpinds(i),2,5,1:103))+squeeze(condindratioerr(lfpinds(i),2,5,1:103)),'r')
    plot(chf(1:103),squeeze(condindratio(lfpinds(i),2,5,1:103))-squeeze(condindratioerr(lfpinds(i),2,5,1:103)),'r')
end

% full width half max of the high gamam reductions
fwhmhgred = [11,NaN,NaN,7,NaN,7,8,9,NaN,8,6,9,NaN,6,NaN,NaN,6,7,7,NaN,7,7,NaN,8,NaN,NaN,11,6,7,8,6,10,7,7,NaN];

% low gamma running no running all animals
for i = 1:16
    figure
    subplot(2,2,1)
    fill(fillx,squeeze(r0tbtfillspecy(lfpinds(i),1,5,:)),'b')
    hold on
    fill(fillx,squeeze(r0tbtfillspecy(lfpinds(i),2,5,:)),'r')
    set(gca,'yscale','log')
        
    subplot(2,2,2)
    fill(fillx,squeeze(r1tbtfillspecy(lfpinds(i),1,5,:)),'b')
    hold on
    fill(fillx,squeeze(r1tbtfillspecy(lfpinds(i),2,5,:)),'r')
    set(gca,'yscale','log')
%     semilogy(chf(1:103),squeeze(condrunS(lfpinds(i),1,5,1:103)));
%     hold on
%     semilogy(chf(1:103),squeeze(condrunS(lfpinds(i),2,5,1:103)),'r');
%     semilogy(chf(1:103),squeeze(condstillS(lfpinds(i),1,5,1:103)),'c');
%     semilogy(chf(1:103),squeeze(condstillS(lfpinds(i),2,5,1:103)),'m');
end

% relative gamma power
% SOM relative to rest of spectrum
figure
[p,s] = signrank(relgammapow(lfpinds(1:16),1,5),relgammapow(lfpinds(1:16),2,5));
plot(1,squeeze(relgammapow(lfpinds(1:16),1,5)),'ko','markerfacecolor','k')
hold on
plot(2,squeeze(relgammapow(lfpinds(1:16),2,5)),'ro','markerfacecolor','r')
for i = 1:16
    plot([1,2],[squeeze(relgammapow(lfpinds(i),1,5)),squeeze(relgammapow(lfpinds(i),2,5))],'k')
end
axis([0,3,1,6])
set(gca,'xticklabel',{'control','SOM Halo'})
set(gca,'xtick',[1,2]);
set(gca,'ytick',[1:6]);
ylabel('relative gamma power to all other frequencies')
title(['relative gamma power. p: ' num2str(p)]);
set(gcf,'OuterPosition',[573   504   242   513])

[p,s] = signrank(relgammapow(lfpinds,1,5),relgammapow(lfpinds,2,5))
disp(['rel. size gamma control: ' num2str(nanmean(relgammapow(lfpinds,1,5))) '+-' num2str(nanstd(relgammapow(lfpinds,1,5))./sqrt(n))]);
disp(['rel. size gamma light: ' num2str(nanmean(relgammapow(lfpinds,2,5))) '+-' num2str(nanstd(relgammapow(lfpinds,2,5))./sqrt(n))]);;
disp(['rel gamma reduction signrank p: ' num2str(p) ' n = ' int2str(length(find(~isnan(relgammapow(lfpinds,1,5)))))])

% PV
figure
[p,s] = signrank(relgammapow(lfpinds(17:34),1,5),relgammapow(lfpinds(17:34),2,5));
plot(1,squeeze(relgammapow(lfpinds(17:34),1,5)),'ko','markerfacecolor','k')
hold on
plot(2,squeeze(relgammapow(lfpinds(17:34),2,5)),'ro','markerfacecolor','r')
for i = 17:34
    plot([1,2],[relgammapow(lfpinds(17:34),1,5),relgammapow(lfpinds(17:34),2,5)],'k')
end
axis([0,3,1.5,8])
set(gca,'xticklabel',{'control','PV Halo'})
set(gca,'xtick',[1,2]);
set(gca,'ytick',[1:8]);
ylabel('relative gamma power to all other frequencies')
title(['relative gamma pwoer. p: ' num2str(p)]);
set(gcf,'OuterPosition',[573   504   242   513])

% all fill spectra light on off
for i = 1:length(lfpinds)    
    figure
    fill(fillx,squeeze(r1tbtfillspecy(lfpinds(i),1,5,:)),'b')
    hold on
    fill(fillx,squeeze(r1tbtfillspecy(lfpinds(i),2,5,:)),'r')
    set(gca,'yscale','log')
    axis([0,100,1,1000])
    legend('no light','light');
    xlabel('frequency [Hz]')
    ylabel('power')
    title([int2str(i),   cllname{lfpinds(i)}]);
end

% all fft spectra light on off
for i = 1:length(lfpinds)    
    figure
    semilogy(chf,squeeze(r1condfftspect(lfpinds(i),1,5,:)))
    hold on
%     semilogy(chf,squeeze(r1condfftspect(lfpinds(i),1,2,:)),'c')
%     semilogy(chf,squeeze(r1condfftspect(lfpinds(i),1,3,:)),'g')
%     semilogy(chf,squeeze(r1condfftspect(lfpinds(i),1,4,:)),'m')
    semilogy(chf,squeeze(r1condfftspect(lfpinds(i),1,1,:)),'r')
    axis([0,100,500,2500])
%     legend('no light','light');
    xlabel('frequency [Hz]')
    ylabel('power')
    title([int2str(i),   cllname{lfpinds(i)}]);
end

% plot all small spectra
for i = 1:length(lfpinds)    
    figure
    fill(fillx,squeeze(r1c0tbtfillspecy(lfpinds(i),1,:)),'b')
%     semilogy(chf,smooth(squeeze(r1condallS(lfpinds(i),1,5,:))),'b')
    hold on
    fill(fillx,squeeze(r1tbtfillspecy(lfpinds(i),1,1,:)),'r')
%     semilogy(chf,smooth(squeeze(r1condallS(lfpinds(i),2,5,:))),'r')
    set(gca,'yscale','log')
    axis([0,100,1,1000])
    % axis([0,100,.9,2.6])
    legend('no stim','small stim');
    xlabel('frequency [Hz]')
    ylabel('power')
    title([int2str(i),   cllname{lfpinds(i)}]);
end


% plot all large L1 L0 spectra
for i = 1:length(lfpinds)    
    figure
    fill(fillx,squeeze(r1tbtfillspecy(lfpinds(i),1,5,:)),'k')
%     semilogy(chf,smooth(squeeze(r1condallS(lfpinds(i),1,5,:))),'b')
    hold on
    fill(fillx,squeeze(r1tbtfillspecy(lfpinds(i),2,5,:)),'r')
    plot(chf(12:103),squeeze(r1condallS(lfpinds(i),1,5,12:103)),'w')
    plot(chf(12:103),squeeze(r1condallS(lfpinds(i),2,5,12:103)),'color',[.5,.5,.5]);
%     semilogy(chf,smooth(squeeze(r1condallS(lfpinds(i),2,5,:))),'r')
    set(gca,'yscale','log')
    axis([0,100,1,1000])
    % axis([0,100,.9,2.6])
    legend('control','light');
    xlabel('frequency [Hz]')
    ylabel('power')
    title([int2str(i),   cllname{lfpinds(i)}]);
end

% anesthetized
figure
fill(fillx,squeeze(r1tbtfillspecy(lfpinds(4),1,5,:)),'b')
hold on
fill(fillx,squeeze(r1tbtfillspecy(lfpinds(4),2,5,:)),'r')
set(gca,'yscale','log')
axis([0,100,1,1000])
% axis([0,100,.9,2.6])
legend('control','light');
xlabel('frequency [Hz]')
ylabel('power')

% anesthetized peak
figure
plot(1,r1peakgammapow(lfpinds,1,5),'ko','markerfacecolor','k')
hold on
plot(2,r1peakgammapow(lfpinds,2,5),'ko','markerfacecolor','k')
% plot(3,r1peakhighgammapow(lfpinds(cphgpeak),1,5),'ko','markerfacecolor','k')
axis([0,3,30,220])
for i = 1:length(lfpinds)
    line([1,2],[r1peakgammapow(lfpinds(i),1,5),r1peakgammapow(lfpinds(i),2,5)],'color','k')
end
set(gca,'xticklabel',{'light off','light on'})
set(gca,'xtick',[1,2]);
ylabel('peak gamma power')
set(gcf,'OuterPosition',[573   504   242   513])

% TODO fix for double index in SOM 16&17 !!!!!!!!!!!!!!!!!!!!!!
% cphgpeak = [1,4,6,7,8,10,12,15,16,18,19,21]; old
% cpnohgpeak = [5,9,11,13,14,17,20,22];
cphgpeak = [1,4,5,6,7,8,10,11,12,14,16,17,18,20,21,22,23,26,27,28,29,30,31,32,33];
cpnohgpeak = [2,3,9,13,15,19,24,25,34]; % 25 and 34 excluded due to line noise, 2/3 no running

% normalized size gamma
for i = 1:length(lfpinds)
    % normalized gamma power and integrated gamma powe
    ngp(i,:) = squeeze(r1peakgammapow(lfpinds(i),1,:))./squeeze(max(r1peakgammapow(lfpinds(i),1,:)));
    nigp(i,:) = squeeze(r1intgammapow(lfpinds(i),1,:))./squeeze(max(r1intgammapow(lfpinds(i),1,:)));
    
    % normalized gamma power with 0 contrast
    hlp = [r1c0peakgammapow(lfpinds(i),1); squeeze(r1peakgammapow(lfpinds(i),1,:))];
    nc0gp(i,:) = hlp./max(hlp);
    nnc0gp(i,:) = hlp;
    
    % non-normalized peak and integrated gamma power
    nngp(i,:) = squeeze(r1peakgammapow(lfpinds(i),1,:));
    nngpl1(i,:) = squeeze(r1peakgammapow(lfpinds(i),2,:));
    pcred(i,:) = nngpl1(i,:)./nngp(i,:);
    nnigp(i,:) = squeeze(r1intgammapow(lfpinds(i),1,:));
    % not only running gamma power
    nr0gp(i,:) = squeeze(peakgammapow(lfpinds(i),1,:))./squeeze(max(peakgammapow(lfpinds(i),1,:)));
    % percent increase and decrease
    pci(i) = nngp(i,5)./nngp(i,1);
    pcd(i) = nngp(i,1)./nngp(i,5);
    
    % normalized and non-normalized high gamma power and percent increase/decrease (running only)
    nhgp(i,:) = squeeze(r1peakhighgammapow(lfpinds(i),1,:))./squeeze(max(r1peakhighgammapow(lfpinds(i),1,:)));
    nnhgp(i,:) = squeeze(r1peakhighgammapow(lfpinds(i),1,:));
    pcihg(i) = nnhgp(i,5)./nnhgp(i,1);
    pcdhg(i) = nnhgp(i,1)./nnhgp(i,5);    
    
    % normalized gamma power with 0 contrast
    hlp = [r1c0peakhighgammapow(lfpinds(i),1); squeeze(r1peakhighgammapow(lfpinds(i),1,:))];
    nc0hgp(i,:) = hlp./max(hlp);
    nnc0hgp(i,:) = hlp;
    
    % spectral difference and spectral percent change
    specdiff(i,:) = squeeze(condrunS(lfpinds(i),2,5,:)-condrunS(lfpinds(i),1,5,:));
    specpc(i,:) = squeeze(condrunS(lfpinds(i),2,5,:)./condrunS(lfpinds(i),1,5,:));
    [pks(i,:),locs(i,:)] = max(-specdiff(i,20:50));
    cdecf(i) = chf(locs(i)+21); % center decrease frequency
    
    beta = [20,40];
    b1r0 = 18;
    b1 = find(f>beta(1),1); b2 = find(f>beta(2),1);
    g1 = find(f>gamma(1),1); g2 = find(f>gamma(2),1);
    for s = 1:5
        for l = 1:2
            smoothspec(i,s,l,:) = smooth(squeeze(condrunS(lfpinds(i),l,s,:)));
            smoothspecr0(i,s,l,:) = smooth(squeeze(condstillS(lfpinds(i),l,s,:)));

            bsig = squeeze(smoothspec(i,s,l,b1:b2));        
            gsig = squeeze(smoothspec(i,s,l,g1:g2));
            bsigr0 = squeeze(smoothspecr0(i,s,l,b1r0:b2));
            gsigr0 = squeeze(smoothspecr0(i,s,l,g1:g2));
        
            if isempty(find(diff(bsig)>0)) % there is no clear beta peak
                slgi(i,s,l) = NaN;
                pgf(i,s,l) = NaN; % peak gamma frequency
            else
                peaks = find(diff(bsig)>0)+1;
                pvs = bsig(peaks);
                bpi = peaks(pvs == max(pvs));
                bpi = bpi+b1-1;
                slgi(i,s,l) = bpi;
                pgf(i,s,l) = chf(bpi); % peak gamma frequency
            end
        
            if isempty(find(diff(gsig)>0)) % there is no clear beta peak
                shgi(i,s,l) = NaN;
            else
                peaks = find(diff(gsig)>0)+1;
                pvs = gsig(peaks);
                gpi = peaks(pvs == max(pvs));
                gpi = gpi+g1-1;
                shgi(i,s,l) = gpi;
            end            
            
            if isempty(find(diff(bsigr0)>0)) % there is no clear beta peak
                slgir0(i,s,l) = NaN;
            else
                peaks = find(diff(bsigr0)>0)+1;
                pvs = bsigr0(peaks);
                bpi = peaks(pvs == max(pvs));
                bpi = bpi+b1r0-1;
                slgir0(i,s,l) = bpi;
            end
            
            if isempty(find(diff(gsigr0)>0)) % there is no clear beta peak
                shgir0(i,s,l) = NaN;
            else
                peaks = find(diff(gsigr0)>0)+1;
                pvs = gsigr0(peaks);
                gpi = peaks(pvs == max(pvs));
                gpi = gpi+g1-1;
                shgir0(i,s,l) = gpi;
            end       
        
            if ~isnan(slgi(i,s,l))
                indpeakgammapow(i,s,l) = smoothspec(i,s,l,slgi(i,s,l));
                relindpeakgammapow(i,s,l) = smoothspec(i,s,l,slgi(i,s,l))./nanmean(smoothspec(i,s,l,12:103),4);
                dipsig = squeeze(smoothspec(i,s,l,slgi(i,s,l)-17:slgi(i,s,l)-3));
                [xx,ind] = min(abs(diff(dipsig(diff(dipsig)>0))));
                hlp = find(diff(dipsig)>0);
                if ~isempty(ind)
                    dipi(i,s,l) = slgi(i,s,l)-17+hlp(ind)-1;
                    inddipgammapow(i,s,l) = smoothspec(i,s,l,dipi(i,s,l));
                else
                    dipi(i,s,l) = NaN;
                    inddipgammapow(i,s,l) = NaN;
                end
            else
                indpeakgammapow(i,s,l) = NaN;
                inddipgammapow(i,s,l) = NaN;
                dipi(i,s,l) = NaN;
            end
        
            if ~isnan(shgi(i,s,l))
                indhighgammapow(i,s,l) = smoothspec(i,s,l,shgi(i,s,l));
                relindhighgammapow(i,s,l) = smoothspec(i,s,l,shgi(i,s,l))./nanmean(smoothspec(i,s,l,12:103),4);
                dipsig = squeeze(smoothspec(i,s,l,shgi(i,s,l)-17:shgi(i,s,l)-3));
                [xx,ind] = min(abs(diff(dipsig(diff(dipsig)>0))));
                hlp = find(diff(dipsig)>0);
                if ~isempty(ind)
                    dipihg(i,s,l) = shgi(i,s,l)-17+hlp(ind)-1;
                    indhighgammadippow(i,s,l) = smoothspec(i,s,l,dipihg(i,s,l));
                else
                    dipihg(i,s,l) = NaN;
                    indhighgammadippow(i,s,l) = NaN;
                end
            else
                indhighgammapow(i,s,l) = NaN;
                indhighgammadippow(i,s,l) = NaN;
                dipihg(i,s,l) = NaN;
            end
            
            if ~isnan(slgir0(i,s,l))
                indpeakgammapowr0(i,s,l) = smoothspecr0(i,s,l,slgir0(i,s,l));
                relindpeakgammapowr0(i,s,l) = smoothspecr0(i,s,l,slgir0(i,s,l))./nanmean(smoothspecr0(i,s,l,12:103),4);
                dipsig = squeeze(smoothspecr0(i,s,l,slgir0(i,s,l)-17:slgir0(i,s,l)-3));
                [xx,ind] = min(abs(diff(dipsig(diff(dipsig)>0))));
                hlp = find(diff(dipsig)>0);
                if ~isempty(ind)
                    dipir0(i,s,l) = slgir0(i,s,l)-17+hlp(ind)-1;
                    inddipgammapowr0(i,s,l) = smoothspecr0(i,s,l,dipir0(i,s,l));
                else
                    dipir0(i,s,l) = NaN;
                    inddipgammapowr0(i,s,l) = NaN;
                end
            else
                indpeakgammapowr0(i,s,l) = NaN;
                inddipgammapowr0(i,s,l) = NaN;
                dipir0(i,s,l) = NaN;
            end            
            
            if ~isnan(shgir0(i,s,l))
                indhighgammapowr0(i,s,l) = smoothspecr0(i,s,l,shgir0(i,s,l));
                relindhighgammapowr0(i,s,l) = smoothspecr0(i,s,l,shgir0(i,s,l))./nanmean(smoothspecr0(i,s,l,12:103),4);
                dipsig = squeeze(smoothspecr0(i,s,l,shgir0(i,s,l)-17:shgir0(i,s,l)-3));
                [xx,ind] = min(abs(diff(dipsig(diff(dipsig)>0))));
                hlp = find(diff(dipsig)>0);
                if ~isempty(ind)
                    dipihgr0(i,s,l) = shgir0(i,s,l)-17+hlp(ind)-1;
                    inddiphighgammapowr0(i,s,l) = smoothspecr0(i,s,l,dipihgr0(i,s,l));
                else
                    dipihgr0(i,s,l) = NaN;
                    inddiphighgammapowr0(i,s,l) = NaN;
                end
            else
                indhighgammapowr0(i,s,l) = NaN;
                inddiphighgammapowr0(i,s,l) = NaN;
                dipihgr0(i,s,l) = NaN;
            end
        end
    end
    
    for l = 1:2
        if l == 1
            smoothc0spec(i,l,:) = smooth(squeeze(r1l0contS(lfpinds(i),:))); 
            smoothc0specr0(i,l,:) = smooth(squeeze(r0l0contS(lfpinds(i),:))); 
        else
            smoothc0spec(i,l,:) = smooth(squeeze(r1l1contS(lfpinds(i),:)));
            smoothc0specr0(i,l,:) = smooth(squeeze(r0l1contS(lfpinds(i),:)));
        end
        bsig = squeeze(smoothc0spec(i,l,b1:b2));
        gsig = squeeze(smoothc0spec(i,l,g1:g2));
        bsigr0 = squeeze(smoothc0specr0(i,l,b1r0:b2));
        gsigr0 = squeeze(smoothc0specr0(i,l,g1:g2));
        
        if isempty(find(diff(bsig)>0)) % there is no clear beta peak
            slgic0(i,l) = NaN;
        else
            peaks = find(diff(bsig)>0)+1;
            pvs = bsig(peaks);
            bpi = peaks(pvs == max(pvs));
            bpi = bpi+b1-1;
            slgic0(i,l) = bpi;
        end
    
        if isempty(find(diff(gsig)>0)) % there is no clear beta peak
            shgic0(i,l) = NaN;
        else
            peaks = find(diff(gsig)>0)+1;
            pvs = gsig(peaks);
            gpi = peaks(pvs == max(pvs));
            gpi = gpi+g1-1;
            shgic0(i,l) = gpi;
        end
        if isempty(find(diff(bsigr0)>0)) % there is no clear beta peak
            slgic0r0(i,l) = NaN;
        else
            peaks = find(diff(bsigr0)>0)+1;
            pvs = bsigr0(peaks);
            bpi = peaks(pvs == max(pvs));
            bpi = bpi+b1r0-1;
            slgic0r0(i,l) = bpi;
        end
        
        if isempty(find(diff(gsigr0)>0)) % there is no clear beta peak
            shgic0r0(i,l) = NaN;
        else
            peaks = find(diff(gsigr0)>0)+1;
            pvs = gsigr0(peaks);
            gpi = peaks(pvs == max(pvs));
            gpi = gpi+g1-1;
            shgic0r0(i,l) = gpi;
        end
    
        if ~isnan(slgic0(i,l))
            indc0peakgammapow(i,l) = smoothc0spec(i,l,slgic0(i,l));
        else
            indc0peakgammapow(i,l) = NaN;
        end
        % gamma power at individual peaks
        hlp = [indc0peakgammapow(i,l), squeeze(indpeakgammapow(i,:,l))];
        indnc0gp(i,l,:) = hlp./max(hlp);
        indnnc0gp(i,l,:) = hlp;
                
        if ~isnan(shgic0(i,l))
            indc0highgammapow(i,l) = smoothc0spec(i,l,shgic0(i,l));
        else
            indc0highgammapow(i,l) = NaN;
        end
        % high gamma power at individual peaks
        hlp = [indc0highgammapow(i,l), squeeze(indhighgammapow(i,:,l))];
        indnc0hgp(i,l,:) = hlp./max(hlp);
        indnnc0hgp(i,l,:) = hlp;
        
        if ~isnan(slgic0r0(i,l))
            indc0peakgammapowr0(i,l) = smoothc0spec(i,l,slgic0r0(i,l));
        else
            indc0peakgammapowr0(i,l) = NaN;
        end
        if ~isnan(shgic0r0(i,l))
            indc0highgammapowr0(i,l) = smoothc0spec(i,l,shgic0r0(i,l));
        else
            indc0highgammapowr0(i,l) = NaN;
        end
        
        % evoked gamma
        if ~isnan(slgi(i,5,l))
            evokedindgamma(i,l) = smoothspec(i,5,l,slgi(i,5,l))./smoothc0spec(i,l,slgi(i,5,l));
        else
            evokedindgamma(i,l) = NaN;
        end
        if ~isnan(slgir0(i,5,l))
            evokedindgammar0(i,l) = smoothspecr0(i,5,l,slgir0(i,5,l))./smoothc0specr0(i,l,slgir0(i,5,l));
        else
            evokedindgammar0(i,l) = NaN;
        end
        if ~isnan(shgi(i,1,l))
            evokedhighgamma(i,l) = smoothspec(i,1,l,shgi(i,1,l))./smoothc0spec(i,l,shgi(i,1,l));
        else
            evokedhighgamma(i,l) = NaN;
        end
        if ~isnan(shgir0(i,1,l))
            evokedhighgammar0(i,l) = smoothspecr0(i,1,l,shgir0(i,1,l))./smoothc0specr0(i,l,shgir0(i,1,l));
        else
            evokedhighgammar0(i,l) = NaN;
        end        
        for s = 1:5
            if ~isnan(slgi(i,s,l))
                evindgamma(i,l,s) = smoothspec(i,s,l,slgi(i,s,l))./smoothc0spec(i,l,slgi(i,s,l));
            else
                evindgamma(i,l,s) = NaN;
            end
        end
        
    end
    
    
    
    figure
    semilogy(squeeze(smoothspec(i,1,1,:)),'b')
    hold on
    semilogy(squeeze(smoothspec(i,2,1,:)),'c')
    semilogy(squeeze(smoothspec(i,3,1,:)),'g')
    semilogy(squeeze(smoothspec(i,4,1,:)),'m')
    semilogy(squeeze(smoothspec(i,5,1,:)),'r')
    if ~isnan(slgi(i,1,1))
        plot(slgi(i,1,1),smoothspec(i,1,1,slgi(i,1,1)),'b*')
    end
    if ~isnan(shgi(i,1,1))
        plot(shgi(i,1,1),smoothspec(i,1,1,shgi(i,1,1)),'b*')
    end
    if ~isnan(dipi(i,1,1))
        plot(dipi(i,1,1),smoothspec(i,1,1,dipi(i,1,1)),'bo')
    end
    if ~isnan(dipihg(i,1,1))
        plot(dipihg(i,1,1),smoothspec(i,1,1,dipihg(i,1,1)),'bo')
    end
    if ~isnan(slgi(i,2,1))
        plot(slgi(i,2,1),smoothspec(i,2,1,slgi(i,2,1)),'c*')
    end
    if ~isnan(shgi(i,2,1))
        plot(shgi(i,2,1),smoothspec(i,2,1,shgi(i,2,1)),'c*')
    end
    if ~isnan(dipi(i,2,1))
        plot(dipi(i,2,1),smoothspec(i,2,1,dipi(i,2,1)),'co')
    end
    if ~isnan(dipihg(i,2,1))
        plot(dipihg(i,2,1),smoothspec(i,2,1,dipihg(i,2,1)),'co')
    end
    if ~isnan(slgi(i,3,1))
        plot(slgi(i,3,1),smoothspec(i,3,1,slgi(i,3,1)),'g*')
    end
    if ~isnan(shgi(i,3,1))
        plot(shgi(i,3,1),smoothspec(i,3,1,shgi(i,3,1)),'g*')
    end
    if ~isnan(dipi(i,3,1))
        plot(dipi(i,3,1),smoothspec(i,3,1,dipi(i,3,1)),'go')
    end
    if ~isnan(dipihg(i,3,1))
        plot(dipihg(i,3,1),smoothspec(i,3,1,dipihg(i,3,1)),'go')
    end
    if ~isnan(slgi(i,4,1))
        plot(slgi(i,4,1),smoothspec(i,4,1,slgi(i,4,1)),'m*')
    end
    if ~isnan(shgi(i,4,1))
        plot(shgi(i,4,1),smoothspec(i,4,1,shgi(i,4,1)),'m*')
    end
    if ~isnan(dipi(i,4,1))
        plot(dipi(i,4,1),smoothspec(i,4,1,dipi(i,4,1)),'mo')
    end
    if ~isnan(dipihg(i,4,1))
        plot(dipihg(i,4,1),smoothspec(i,4,1,dipihg(i,4,1)),'mo')
    end
    if ~isnan(slgi(i,5,1))
        plot(slgi(i,5,1),smoothspec(i,5,1,slgi(i,5,1)),'r*')
    end
    if ~isnan(shgi(i,5,1))
        plot(shgi(i,5,1),smoothspec(i,5,1,shgi(i,5,1)),'r*')
    end
    if ~isnan(dipi(i,5,1))
        plot(dipi(i,5,1),smoothspec(i,5,1,dipi(i,5,1)),'ro')
    end
    if ~isnan(dipihg(i,5,1))
        plot(dipihg(i,5,1),smoothspec(i,5,1,dipihg(i,5,1)),'ro')
    end
    
%     P = findpeaksG(chf(1:104),squeeze(r1condallS(lfpinds(i),1,1,1:104)),0,3,3,3,1)
end

relindpeakgammapow(relindpeakgammapow == 0) = NaN;
relindpeakgammapowr0(relindpeakgammapowr0 == 0) = NaN;
relindhighgammapow(relindhighgammapow == 0) = NaN;
relindhighgammapowr0(relindhighgammapowr0 == 0) = NaN;
reldiplgpowr0 = indpeakgammapowr0./inddipgammapowr0; 
reldiplgpow = indpeakgammapow./inddipgammapow;% power relative to dip
reldiphgpowr0 = indhighgammapowr0./inddiphighgammapowr0;
reldiphgpow = indhighgammapow./indhighgammadippow;
abspeakheight = indpeakgammapow-inddipgammapow;
abspeakheightr0 = indpeakgammapowr0-inddipgammapowr0;
abshgpeakheight = indhighgammapow-indhighgammadippow;
abshgpeakheightr0 = indhighgammapowr0-inddiphighgammapowr0;

% SOM
% absolute low gamma power
sz = 5;
figure
[p,s] = signrank(indpeakgammapow(1:16,sz,1),indpeakgammapow(1:16,sz,2));
plot(1,indpeakgammapow(1:16,sz,1),'ko','markerfacecolor','k')
hold on
plot(2,indpeakgammapow(1:16,sz,2),'ro','markerfacecolor','r')
for i = 1:16
    plot([1,2],[indpeakgammapow(i,sz,1),indpeakgammapow(i,sz,2)],'k')
end
axis([0,3,0,350])
set(gca,'xtick',[1,2])
set(gca,'xticklabel',{'control','SOM Halo'})
ylabel('peak gamma power')
title(['peak gamma power. p: ' num2str(p)]);
set(gcf,'OuterPosition',[573   504   242   513])

% relative low gamma power
figure
[p,s,stats] = signrank(relindpeakgammapow(1:16,5,1),relindpeakgammapow(1:16,5,2));
plot(1,relindpeakgammapow(1:16,5,1),'ko','markerfacecolor','k')
hold on
plot(2,relindpeakgammapow(1:16,5,2),'ro','markerfacecolor','r')
for i = 1:16
    plot([1,2],[relindpeakgammapow(i,5,1),relindpeakgammapow(i,5,2)],'k')
end
axis([0,3,1.5,5.5])
set(gca,'xtick',[1,2])
set(gca,'ytick',[2:5]);
set(gca,'xticklabel',{'control','SOM Halo'})
ylabel('relative gamma power')
title(['relative gamma power. p: ' num2str(p)]);
set(gcf,'OuterPosition',[573   504   242   513])

% relative to dip low gamma power
figure
[p,s,stats] = signrank(reldiplgpow(1:16,5,1),reldiplgpow(1:16,5,2));
plot(1,reldiplgpow(1:16,5,1),'ko','markerfacecolor','k')
hold on
plot(2,reldiplgpow(1:16,5,2),'ro','markerfacecolor','r')
for i = 1:16
    plot([1,2],[reldiplgpow(i,5,1),reldiplgpow(i,5,2)],'k')
end
axis([0,3,1,4])
set(gca,'xtick',[1,2]);
set(gca,'xticklabel',{'control','SOM Halo'})
ylabel('peak gamma power/trough gamma power')
title(['peak height relative to trough. p: ' num2str(p)]);
set(gcf,'OuterPosition',[573   504   242   513])

% absolute peak height low gamma
figure
[p,s] = signrank(abspeakheight(1:16,5,1),abspeakheight(1:16,5,2));
plot(1,abspeakheight(1:16,5,1),'ko','markerfacecolor','k')
hold on
plot(2,abspeakheight(1:16,5,2),'ro','markerfacecolor','r')
for i = 1:16
    plot([1,2],[abspeakheight(i,5,1),abspeakheight(i,5,2)],'k')
end
axis([0,3,0,250])
set(gca,'xtick',[1,2])
set(gca,'xticklabel',{'control','SOM Halo'})
ylabel('peak height')
title(['peak height  p: ' num2str(p)]);
set(gcf,'OuterPosition',[573   504   242   513])

% evoked low gamma power only for largest size
figure
[p,s] = signrank(evokedindgamma(1:16,1),evokedindgamma(1:16,2));
plot(1,evokedindgamma(1:16,1),'ko','markerfacecolor','k')
hold on
plot(2,evokedindgamma(1:16,2),'ro','markerfacecolor','r')
for i = 1:16
    plot([1,2],[evokedindgamma(i,1),evokedindgamma(i,2)],'k')
end
axis([0,3,0,10])
set(gca,'xtick',[1,2])
set(gca,'xticklabel',{'control','SOM Halo'})
ylabel('evoked gamma power')
title(['evoked gamma power p: ' num2str(p)]);
set(gcf,'OuterPosition',[573   504   242   513])

% evoked low gamma power for all sizes
sz = 1;
figure
[p,s] = signrank(evindgamma(1:16,1,sz),evindgamma(1:16,2,sz));
plot(1,evindgamma(1:16,1,sz),'ko','markerfacecolor','k')
hold on
plot(2,evindgamma(1:16,2,sz),'ro','markerfacecolor','r')
for i = 1:16
    plot([1,2],[evindgamma(i,1,sz),evindgamma(i,2,sz)],'k')
end
axis([0,3,0,10])
set(gca,'xtick',[1,2])
set(gca,'xticklabel',{'control','SOM Halo'})
ylabel('evoked gamma power')
title(['evoked gamma power p: ' num2str(p)]);
set(gcf,'OuterPosition',[573   504   242   513])

% PV
% absolute low gamma power
sz = 5;
figure
[p,s,stats] = signrank(indpeakgammapow(17:34,sz,1),indpeakgammapow(17:34,sz,2));
plot(1,indpeakgammapow(17:34,sz,1),'ko','markerfacecolor','k')
hold on
plot(2,indpeakgammapow(17:34,sz,2),'ro','markerfacecolor','r')
for i = 17:34
    plot([1,2],[indpeakgammapow(i,sz,1),indpeakgammapow(i,sz,2)],'k')
end
axis([0,3,0,2500])
set(gca,'xtick',[1,2])
set(gca,'xticklabel',{'control','PV Halo'})
ylabel('peak gamma power')
title(['peak gamma power. p: ' num2str(p)]);
set(gcf,'OuterPosition',[573   504   242   513])

% relative low gamma power
figure
[p,s,stats] = signrank(relindpeakgammapow(17:34,5,1),relindpeakgammapow(17:34,5,2));
plot(1,relindpeakgammapow(17:34,5,1),'ko','markerfacecolor','k')
hold on
plot(2,relindpeakgammapow(17:34,5,2),'ro','markerfacecolor','r')
for i = 17:34
    plot([1,2],[relindpeakgammapow(i,5,1),relindpeakgammapow(i,5,2)],'k')
end
axis([0,3,1,7])
set(gca,'xtick',[1,2])
set(gca,'xticklabel',{'control','PV Halo'})
ylabel('relative gamma power')
title(['relative gamma power. p: ' num2str(p)]);
set(gcf,'OuterPosition',[573   504   242   513])

% relative to dip low gamma power
figure
[p,s,stats] = signrank(reldiplgpow(17:34,5,1),reldiplgpow(17:34,5,2));
plot(1,reldiplgpow(17:34,5,1),'ko','markerfacecolor','k')
hold on
plot(2,reldiplgpow(17:34,5,2),'ro','markerfacecolor','r')
for i = 17:34
    plot([1,2],[reldiplgpow(i,5,1),reldiplgpow(i,5,2)],'k')
end
axis([0,3,1,5.5])
set(gca,'xtick',[1,2]);
set(gca,'ytick',[1:5]);
set(gca,'xticklabel',{'control','PV Halo'})
ylabel('peak gamma power/trough gamma power')
title(['peak height relative to trough. p: ' num2str(p)]);
set(gcf,'OuterPosition',[573   504   242   513])

% absolute peak height low gamma
figure
[p,s] = signrank(abspeakheight(17:34,5,1),abspeakheight(17:34,5,2));
plot(1,abspeakheight(17:34,5,1),'ko','markerfacecolor','k')
hold on
plot(2,abspeakheight(17:34,5,2),'ro','markerfacecolor','r')
for i = 17:34
    plot([1,2],[abspeakheight(i,5,1),abspeakheight(i,5,2)],'k')
end
axis([0,3,0,2000])
set(gca,'xtick',[1,2])
set(gca,'xticklabel',{'control','PV Halo'})
ylabel('peak height')
title(['peak height  p: ' num2str(p)]);
set(gcf,'OuterPosition',[573   504   242   513])

% evoked low gamma power
figure
[p,s] = signrank(evokedindgamma(17:34,1),evokedindgamma(17:34,2));
plot(1,evokedindgamma(17:34,1),'ko','markerfacecolor','k')
hold on
plot(2,evokedindgamma(17:34,2),'ro','markerfacecolor','r')
for i = 17:34
    plot([1,2],[evokedindgamma(i,1),evokedindgamma(i,2)],'k')
end
axis([0,3,0,12])
set(gca,'xtick',[1,2])
set(gca,'xticklabel',{'control','PV Halo'})
ylabel('evoked gamma power')
title(['evoked gamma power p: ' num2str(p)]);
set(gcf,'OuterPosition',[573   504   242   513])

% evoked low gamma power for all sizes
sz = 5;
figure
[p,s] = signrank(evindgamma(17:34,1,sz),evindgamma(17:34,2,sz));
plot(1,evindgamma(17:34,1,sz),'ko','markerfacecolor','k')
hold on
plot(2,evindgamma(17:34,2,sz),'ro','markerfacecolor','r')
for i = 17:34
    plot([1,2],[evindgamma(i,1,sz),evindgamma(i,2,sz)],'k')
end
axis([0,3,0,10])
set(gca,'xtick',[1,2])
set(gca,'xticklabel',{'control','SOM Halo'})
ylabel('evoked gamma power')
title(['evoked gamma power p: ' num2str(p)]);
set(gcf,'OuterPosition',[573   504   242   513])

% SOM
% absolute high gamma power
figure
[p,s,stats] = signrank(indhighgammapow(1:16,1,1),indhighgammapow(1:16,1,2));
plot(1,indhighgammapow(1:16,1,1),'ko','markerfacecolor','k')
hold on
plot(2,indhighgammapow(1:16,1,2),'ro','markerfacecolor','r')
for i = 1:16
    plot([1,2],[indhighgammapow(i,1,1),indhighgammapow(i,1,2)],'k')
end
axis([0,3,0,120])
set(gca,'xtick',[1,2])
set(gca,'xticklabel',{'control','SOM Halo'})
ylabel('peak gamma power')
title(['peak gamma power. p: ' num2str(p)]);
set(gcf,'OuterPosition',[573   504   242   513])

% relative high gamma power
figure
[p,s] = signrank(relindhighgammapow(1:16,1,1),relindhighgammapow(1:16,1,2));
plot(1,relindhighgammapow(1:16,1,1),'ko','markerfacecolor','k')
hold on
plot(2,relindhighgammapow(1:16,1,2),'ro','markerfacecolor','r')
for i = 1:16
    plot([1,2],[relindhighgammapow(i,1,1),relindhighgammapow(i,1,2)],'k')
end
axis([0,3,.4,2.4])
set(gca,'xtick',[1,2])
set(gca,'xticklabel',{'control','SOM Halo'})
ylabel('relative gamma power')
title(['relative gamma power. p: ' num2str(p)]);
set(gcf,'OuterPosition',[573   504   242   513])

% relative to dip high gamma power
figure
[p,s] = signrank(reldiphgpow(1:16,1,1),reldiphgpow(1:16,1,2));
plot(1,reldiphgpow(1:16,1,1),'ko','markerfacecolor','k')
hold on
plot(2,reldiphgpow(1:16,1,2),'ro','markerfacecolor','r')
for i = 1:16
    plot([1,2],[reldiphgpow(i,1,1),reldiphgpow(i,1,2)],'k')
end
axis([0,3,1,4.5])
set(gca,'xtick',[1,2]);
set(gca,'xticklabel',{'control','SOM Halo'})
ylabel('peak gamma power/trough gamma power')
title(['peak height relative to trough. p: ' num2str(p)]);
set(gcf,'OuterPosition',[573   504   242   513])

% absolute peak height high gamma
figure
[p,s] = signrank(abshgpeakheight(1:16,1,1),abshgpeakheight(1:16,1,2));
plot(1,abshgpeakheight(1:16,1,1),'ko','markerfacecolor','k')
hold on
plot(2,abshgpeakheight(1:16,1,2),'ro','markerfacecolor','r')
for i = 1:16
    plot([1,2],[abshgpeakheight(i,1,1),abshgpeakheight(i,1,2)],'k')
end
axis([0,3,0,60])
set(gca,'xtick',[1,2])
set(gca,'xticklabel',{'control','SOM Halo'})
ylabel('peak height')
title(['peak height  p: ' num2str(p)]);
set(gcf,'OuterPosition',[573   504   242   513])

% evoked high gamma power
figure
[p,s] = signrank(evokedhighgamma(1:16,1),evokedhighgamma(1:16,2));
plot(1,evokedhighgamma(1:16,1),'ko','markerfacecolor','k')
hold on
plot(2,evokedhighgamma(1:16,2),'ro','markerfacecolor','r')
for i = 1:16
    plot([1,2],[evokedhighgamma(i,1),evokedhighgamma(i,2)],'k')
end
axis([0,3,.4,2])
set(gca,'xtick',[1,2])
set(gca,'xticklabel',{'control','SOM Halo'})
ylabel('evoked gamma power')
title(['evoked gamma power p: ' num2str(p)]);
set(gcf,'OuterPosition',[573   504   242   513])


% PV
% absolute high gamma power
figure
[p,s,stats] = signrank(indhighgammapow(17:34,1,1),indhighgammapow(17:34,1,2));
plot(1,indhighgammapow(17:34,1,1),'ko','markerfacecolor','k')
hold on
plot(2,indhighgammapow(17:34,1,2),'ro','markerfacecolor','r')
for i = 17:34
    plot([1,2],[indhighgammapow(i,1,1),indhighgammapow(i,1,2)],'k')
end
axis([0,3,0,90])
set(gca,'xtick',[1,2])
set(gca,'xticklabel',{'control','PV Halo'})
ylabel('peak gamma power')
title(['peak gamma power. p: ' num2str(p)]);
set(gcf,'OuterPosition',[573   504   242   513])

% relative high gamma power
figure
[p,s] = signrank(relindhighgammapow(17:34,1,1),relindhighgammapow(17:34,1,2));
plot(1,relindhighgammapow(17:34,1,1),'ko','markerfacecolor','k')
hold on
plot(2,relindhighgammapow(17:34,1,2),'ro','markerfacecolor','r')
for i = 17:34
    plot([1,2],[relindhighgammapow(i,1,1),relindhighgammapow(i,1,2)],'k')
end
axis([0,3,0,2.5])
set(gca,'xtick',[1,2])
set(gca,'xticklabel',{'control','PV Halo'})
ylabel('relative gamma power')
title(['relative gamma power. p: ' num2str(p)]);
set(gcf,'OuterPosition',[573   504   242   513])

% relative to dip high gamma power
figure
[p,s] = signrank(reldiphgpow(17:34,1,1),reldiphgpow(17:34,1,2));
plot(1,reldiphgpow(17:34,1,1),'ko','markerfacecolor','k')
hold on
plot(2,reldiphgpow(17:34,1,2),'ro','markerfacecolor','r')
for i = 17:34
    plot([1,2],[reldiphgpow(i,1,1),reldiphgpow(i,1,2)],'k')
end
axis([0,3,.5,4.1])
set(gca,'xtick',[1,2]);
set(gca,'xticklabel',{'control','PV Halo'})
ylabel('peak gamma power/trough gamma power')
title(['peak height relative to trough. p: ' num2str(p)]);
set(gcf,'OuterPosition',[573   504   242   513])

% absolute peak height high gamma
figure
[p,s] = signrank(abshgpeakheight(17:34,1,1),abshgpeakheight(17:34,1,2));
plot(1,abshgpeakheight(17:34,1,1),'ko','markerfacecolor','k')
hold on
plot(2,abshgpeakheight(17:34,1,2),'ro','markerfacecolor','r')
for i = 17:34
    plot([1,2],[abshgpeakheight(i,1,1),abshgpeakheight(i,1,2)],'k')
end
axis([0,3,-20,50])
set(gca,'xtick',[1,2])
set(gca,'xticklabel',{'control','PV Halo'})
ylabel('peak height')
title(['peak height  p: ' num2str(p)]);
set(gcf,'OuterPosition',[573   504   242   513])

% evoked high gamma power
figure
[p,s] = signrank(evokedhighgamma(17:34,1),evokedhighgamma(17:34,2));
plot(1,evokedhighgamma(17:34,1),'ko','markerfacecolor','k')
hold on
plot(2,evokedhighgamma(17:34,2),'ro','markerfacecolor','r')
for i = 17:34
    plot([1,2],[evokedhighgamma(i,1),evokedhighgamma(i,2)],'k')
end
axis([0,3,.5,1.4])
set(gca,'xtick',[1,2])
set(gca,'xticklabel',{'control','PV Halo'})
ylabel('evoked gamma power')
title(['evoked gamma power p: ' num2str(p)]);
set(gcf,'OuterPosition',[573   504   242   513])



% plot this with relindpeakgammapow and -r0 SOM reduction 0.007 signrank
% Optogenetics relative to running
% absolute power SOM
figure
subplot(1,2,1)
plot(1,indpeakgammapow(1:16,5,1),'ko','markerfacecolor','k')
hold on
plot(2,indpeakgammapow(1:16,5,2),'ro','markerfacecolor','r')
for i = 1:16
    plot([1,2],[indpeakgammapow(1:16,5,1),indpeakgammapow(1:16,5,2)], 'k')
end
axis([0,3,0,350])
set(gca,'xtick',[1,2]);
set(gca,'xticklabel',{'control','SOM Halo'})
ylabel('peak gamma power')
title('running')

subplot(1,2,2)
plot(1,indpeakgammapowr0(1:16,5,1),'ko','markerfacecolor','k')
hold on
plot(2,indpeakgammapowr0(1:16,5,2),'ro','markerfacecolor','r')
for i = 1:16
    plot([1,2],[indpeakgammapowr0(1:16,5,1),indpeakgammapowr0(1:16,5,2)], 'k')
end
axis([0,3,0,350])
set(gca,'xtick',[1,2]);
set(gca,'xticklabel',{'control','SOM Halo'})
ylabel('peak gamma power')
title('quiescent')

% relative power SOM
figure
p = signrank(relindpeakgammapow(1:16,5,1),relindpeakgammapow(1:16,5,2));
subplot(1,2,1)
plot(1,relindpeakgammapow(1:16,5,1),'ko','markerfacecolor','k')
hold on
plot(2,relindpeakgammapow(1:16,5,2),'ro','markerfacecolor','r')
for i = 1:16
    plot([1,2],[relindpeakgammapow(1:16,5,1),relindpeakgammapow(1:16,5,2)], 'k')
end
axis([0,3,1,6])
set(gca,'xtick',[1,2]);
set(gca,'xticklabel',{'control','SOM Halo'})
ylabel('relative gamma power')
title(['running p: ' num2str(p)])

subplot(1,2,2)
p = signrank(relindpeakgammapowr0(1:16,5,1),relindpeakgammapowr0(1:16,5,2));
plot(1,relindpeakgammapowr0(1:16,5,1),'ko','markerfacecolor','k')
hold on
plot(2,relindpeakgammapowr0(1:16,5,2),'ro','markerfacecolor','r')
for i = 1:16
    plot([1,2],[relindpeakgammapowr0(1:16,5,1),relindpeakgammapowr0(1:16,5,2)], 'k')
end
axis([0,3,1,6])
set(gca,'xtick',[1,2]);
set(gca,'xticklabel',{'control','SOM Halo'})
ylabel('relative gamma power')
title(['quiescent p: ' num2str(p)])

% abolute power PV subplot with running
figure
subplot(1,2,1)
plot(1,indpeakgammapow(17:34,5,1),'ko','markerfacecolor','k')
hold on
plot(2,indpeakgammapow(17:34,5,2),'ro','markerfacecolor','r')
for i = 17:34
    plot([1,2],[indpeakgammapow(17:34,5,1),indpeakgammapow(17:34,5,2)], 'k')
end
axis([0,3,0,1200])
set(gca,'xtick',[1,2]);
set(gca,'xticklabel',{'control','PV Halo'})
ylabel('peak gamma power')
title('running')

subplot(1,2,2)
plot(1,indpeakgammapowr0(17:34,5,1),'ko','markerfacecolor','k')
hold on
plot(2,indpeakgammapowr0(17:34,5,2),'ro','markerfacecolor','r')
for i = 17:34
    plot([1,2],[indpeakgammapowr0(17:34,5,1),indpeakgammapowr0(17:34,5,2)], 'k')
end
axis([0,3,0,1200])
set(gca,'xtick',[1,2]);
set(gca,'xticklabel',{'control','PV Halo'})
ylabel('peak gamma power')
title('quiescent')

% relative power PV subplot with running
figure
p = signrank(relindpeakgammapow(17:34,5,1),relindpeakgammapow(17:34,5,2));
subplot(1,2,1)
plot(1,relindpeakgammapow(17:34,5,1),'ko','markerfacecolor','k')
hold on
plot(2,relindpeakgammapow(17:34,5,2),'ro','markerfacecolor','r')
for i = 17:34
    plot([1,2],[relindpeakgammapow(17:34,5,1),relindpeakgammapow(17:34,5,2)], 'k')
end
axis([0,3,1,7])
set(gca,'xtick',[1,2]);
set(gca,'xticklabel',{'control','PV Halo'})
ylabel('relative gamma power')
title(['running p: ' num2str(p)])

subplot(1,2,2)
p = signrank(relindpeakgammapowr0(17:34,5,1),relindpeakgammapowr0(17:34,5,2));
plot(1,relindpeakgammapowr0(17:34,5,1),'ko','markerfacecolor','k')
hold on
plot(2,relindpeakgammapowr0(17:34,5,2),'ro','markerfacecolor','r')
for i = 17:34
    plot([1,2],[relindpeakgammapowr0(17:34,5,1),relindpeakgammapowr0(17:34,5,2)], 'k')
end
axis([0,3,1,7])
set(gca,'xtick',[1,2]);
set(gca,'xticklabel',{'control','PV Halo'})
ylabel('relative gamma power')
title(['quiescent p: ' num2str(p)])


% running on chnages in relative gamma SOM and PV Halo
a = relindpeakgammapow(:,5,2)./relindpeakgammapow(:,5,1);
b = relindpeakgammapowr0(:,5,2)./relindpeakgammapowr0(:,5,1);
figure
plot(a(1:16),b(1:16),'ro','markerfacecolor','r')
hold on
plot(a(17:34),b(17:34),'bo','markerfacecolor','b')
axis([0.4,1.4,0.4,1.4])
axis square
hold on
plot([0.4,1.4],[1,1],'k')
plot([1,1],[0.4,1.4],'k')
% set(gca,'xtick',[0:2]);
% set(gca,'ytick',[0:2]);
legend('SOM Halo','PV Halo')
xlabel('change with light during running')
ylabel('change with light during quiescence')
title('effect of running on relative gamma changes')



% gamma power relative to dip SOM
figure
[p,s] = signrank(reldiplgpow(1:16,5,1),reldiplgpow(1:16,5,2));
plot(1,reldiplgpow(1:16,5,1),'ko','markerfacecolor','k')
hold on
plot(2,reldiplgpow(1:16,5,2),'ro','markerfacecolor','r')
for i = 1:16
    plot([1,2],[reldiplgpow(i,5,1),reldiplgpow(i,5,2)],'k')
end
axis([0,3,1,4])
set(gca,'xtick',[1,2]);
set(gca,'xticklabel',{'control','SOM Halo'})
ylabel('peak gamma power/trough gamma power')
title(['peak height relative to trough. p: ' num2str(p)]);
set(gcf,'OuterPosition',[573   504   242   513])


% gamma power relative to dip PV
figure
[p,s] = signrank(reldiplgpow(17:34,5,1),reldiplgpow(17:34,5,2));
plot(1,reldiplgpow(17:34,5,1),'ko','markerfacecolor','k')
hold on
plot(2,reldiplgpow(17:34,5,2),'ro','markerfacecolor','r')
for i = 17:34
    plot([1,2],[reldiplgpow(i,5,1),reldiplgpow(i,5,2)],'k')
end
axis([0,3,1,5.5])
set(gca,'xtick',[1,2]);
set(gca,'xticklabel',{'control','PV Halo'})
ylabel('peak gamma power/trough gamma power')
title(['peak height relative to trough. p: ' num2str(p)]);
set(gcf,'OuterPosition',[573   504   242   513])

% high gamma power running non-running
figure
[p,s] = signrank(indhighgammapow(:,1,1),indhighgammapowr0(:,1,1));
plot(1,indhighgammapow(:,1,1),'ko','markerfacecolor','k')
hold on
plot(2,indhighgammapowr0(:,1,1),'ko','markerfacecolor','k')
for i = 1:length(lfpinds)
    plot([1,2],[indhighgammapow(i,1,1),indhighgammapowr0(i,1,1)],'k')
end
axis([0,3,0,90])
set(gca,'xtick',[1,2]);
set(gca,'xticklabel',{'running','quiescent'})
ylabel('high gamma power')
title(['high gamma power with state. p: ' num2str(p)]);
set(gcf,'OuterPosition',[573   504   242   513])

% individual high gamma peak with size
n = length(find(~isnan(indnc0hgp(cphgpeak,1,2))));
figure
errorbar([0,sizes],squeeze(nanmean(indnc0hgp(cphgpeak,1,:))),squeeze(nanstd(indnc0hgp(cphgpeak,1,:)))./sqrt(length(find(~isnan(indnc0hgp(cphgpeak,2,1))))),'ko-','markerfacecolor','k','linewidth',2)
xlabel('size (degrees)')
ylabel('normlized individual high gamma peak')
title(['n: ' int2str(n)]);

% high gamma center frequency
hlp = round(chf(shgi(cphgpeak,1)));
figure
hist(hlp,57:67);
set(gca,'ytick',[1,2,3,4]);
xlabel('frequency (Hz)')
ylabel('bin count')
title('center frequencies of high gamma rhythm')


% find out if center frequency decreases with size
for i = 1:size(slgi,1)
    if isempty(find(diff(slgi(i,2:5,1))>0))
        lastfourdecend(i) = 1;
    else
        lastfourdecend(i) = 0;
    end
    if isempty(find(diff(slgi(i,3:5,1))>0))
        lastthreedecend(i) = 1;
    else
        lastthreedecend(i) = 0;
    end
end
lastfourdecend([2,3]) = NaN; lastthreedecend([2,3]) = NaN;
for i = 1:34
    mx(i) = max(slgi(i,2:5,1));
    mn(i) = min(slgi(i,2:5,1));
end
range = mx-mn;

% Fig 1
disp([int2str(length(find(lastfourdecend))) ' out of ' int2str(length(find(~isnan(lastfourdecend)))) ' decrease center freq with size'])
disp(['range: ' num2str(nanmean(range(find(lastfourdecend)))) ' +- ' num2str(nanstd(range(find(lastfourdecend))))]);

% supplemental ANOVA and error bars
figure
sizes(5) = 60;
errorbar(sizes(2:5),nanmean(pgf(:,2:5,1)),nanstd(pgf(:,2:5,1))./sqrt(length(find(~isnan(pgf(:,5,1))))),'ko-','linewidth',2,'markerfacecolor','k')

[p,tbl,stats] = anova1(pgf(:,2:5,1));
[p,tbl,stats] = kruskalwallis(pgf(:,2:5,1));

% dumme version
a = slgi(:,2,1); b = slgi(:,3,1); c = slgi(:,4,1); d = slgi(:,5,1);
a(isnan(a)) = []; b(isnan(b)) = []; c(isnan(c)) = []; d(isnan(d)) = [];
sizes(5) = 60;
figure
errorbar(sizes(2:5),[mean(chf(a)),mean(chf(b)),mean(chf(c)),mean(chf(d))],[std(chf(a))./sqrt(length(a)),std(chf(b))./sqrt(length(b)),std(chf(c))./sqrt(length(c)),std(chf(d))./sqrt(length(d))],'ko-','linewidth',2,'markerfacecolor','k')




% indivisual peaks
ipwc0 = [indc0peakgammapow(:,1),indpeakgammapow(:,:,1)];
for i = 1:34
    nipgp(i,:) = ipwc0(i,:)./max(ipwc0(i,:)); % norm individual peak gamma power
end
for i = 1:6
    nipgperr(i) = nanstd(nipgp(:,i))./sqrt(length(find(~isnan(nipgp(:,i)))));
end
    
for i = 1:5
    impgp(i) = nanmean(indpeakgammapow(:,i,1)); % individual mean peak gamma pow
    imgperr(i) = nanstd(indpeakgammapow(:,i,1))./sqrt(length(find(~isnan(indpeakgammapow(:,i,1)))));    
end
ic0mpgperr = nanstd(indc0peakgammapow(:,1))./sqrt(length(find(~isnan(indc0peakgammapow(:,1)))));

% not normalized power at individual peaks
figure
errorbar([0,sizes],[nanmean(indc0peakgammapow(:,1)),impgp],[ic0mpgperr,imgperr])

% normalized power at individaul peaks
sizes(5) = 60;
figure
errorbar([0,sizes],nanmean(nipgp),nipgperr,'ko-','markerfacecolor','k','linewidth',2)
xlabel('stimulus size [dva]')
ylabel('normalized individual peak gamma power')

% mean over broad peak (+-10Hz)
for i = 1:34
    for l = 1:2
        for j = 1:5
            if ~isnan(slgi(i,j,l))
                bpgpow(i,j,l) = nanmean(smoothspec(i,j,l,slgi(i,j,l)-10:slgi(i,j,l)+10),4);
            else
                bpgpow(i,j,l) = NaN;
            end
        end
        if ~isnan(slgic0(i,l))
            c0bgpow(i,l) = nanmean(smoothc0spec(i,l,slgic0(i,l)-10:slgic0(i,l)+10),3);
        else
            c0bgpow(i,l) = NaN;
        end
    end
    allbpgpow(i,:) = [c0bgpow(i,1),bpgpow(i,:,1)];
    nbpgpow(i,:) = allbpgpow(i,:,1)./max(allbpgpow(i,:,1));
end
for i = 1:6
    nbpgpowerr(i) = nanstd(nbpgpow(:,i))./sqrt(length(find(~isnan(nbpgpow(:,i)))));
end

figure
errorbar([0,sizes],nanmean(nbpgpow),nbpgpowerr,'ko-','markerfacecolor','k','linewidth',2)
xlabel('stimulus size [dva]')
ylabel('normalized broad peak gamma power')

cdecf(isnan(condrunS(lfpinds,1,5))) = NaN;
specdifferr = nanstd(specdiff,1,1)./sqrt(length(find(~isnan(specdiff(:,1)))));
specpcerr = nanstd(specpc(:,:),1,1)./sqrt(length(find(~isnan(specpc(:,1)))));
n = length(find(~isnan(specpc(:,1))));
fillspecpc = [(squeeze(nanmean(specpc(:,12:103)))+specpcerr(12:103)),fliplr(squeeze(nanmean(specpc(:,12:103)))-specpcerr(12:103))];

% FIGURE suppl 2
figure
fill(fillx,fillspecpc,'k')
line([0,100],[1,1],'color','k')
% axis([0,100,.4,1.6])
axis([0,100,.9,3])
xlabel('frequency [Hz]')
ylabel('fold power change with SOM inactivation')
title(['n = ' int2str(n)])
% ylabel('fold power change with PV inactivation')

figure
sizes(5) = 60;
% plot(sizes,ngp','color',[.7,.7,.7])
% plot([0,sizes],nc0gp','color',[.7,.7,.7])
n = length(find(~isnan(ngp(:,5))));
hold on
% errorbar(sizes,squeeze(nanmean(ngp)),squeeze(nanstd(ngp))./sqrt(n),'ko-','markerfacecolor','k','linewidth',2)
errorbar([0,sizes],squeeze(nanmean(nc0gp)),squeeze(nanstd(nc0gp))./sqrt(n),'ko-','markerfacecolor','k','linewidth',2)
axis([-3,65,.3,1.1])
xlabel('stimulus size [dva]')
ylabel('normalized gamma power')

[p, table, stats] = kruskalwallis(nc0gp);

figure
errorbar([0,sizes],nanmean(indnc0gp),nanstd(indnc0gp)./sqrt(length(find(~isnan(indnc0gp(:,5))))),'ko-','markerfacecolor','k','linewidth',2);
axis([-3,65,.3,1])
xlabel('stimulus size [dva]')
ylabel('normalized gamma power individual peaks')

% Suppl 1
figure
% plot(sizes,nhgp(cphgpeak,:)','color',[.7,.7,.7])
% plot([0,sizes],nc0hgp(cphgpeak,:)','color',[.7,.7,.7])
% n = length(find(~isnan(nhgp(:,5))));
n = length(cphgpeak);
[p,t,s] = kruskalwallis(nnc0hgp(cphgpeak,:));
% errorbar(sizes,squeeze(nanmean(nhgp(cphgpeak,:))),squeeze(nanstd(nhgp(cphgpeak,:)))./sqrt(n),'ko-','markerfacecolor','k','linewidth',2)
errorbar([0,sizes],squeeze(nanmean(nc0hgp(cphgpeak,:))),squeeze(nanstd(nc0hgp(cphgpeak,:)))./sqrt(n),'ko-','markerfacecolor','k','linewidth',2)
axis([-3,65,.3,1])
xlabel('stimulus size [dva]')
ylabel('normalized high gamma power')
title(['n: ' int2str(n) ' p: ' num2str(p)])

% supplemental figure 2h - percent reduction with size
pcred = (1-pcred).*100;
figure
errorbar(sizes,nanmean(pcred,1),nanstd(pcred,1,1)./sqrt(length(find(~isnan(pcred(:,1))))),'ko-','markerfacecolor','k')
xlabel('stimulus size [dva]')
ylabel('percent reduction of gamma power')
[p,table,stats] = kruskalwallis(pcred);

% [p,table,stats] = kruskalwallis(ngp);
% [p,table,stats] = kruskalwallis(nhgp);

% high gamma examples
% stim and state dependence
figure
fill(fillx,squeeze(r0tbtfillspecy(115,1,1,:)),'b')
hold on
fill(fillx,squeeze(r1tbtfillspecy(115,1,1,:)),'c')
fill(fillx,squeeze(r1tbtfillspecy(115,1,5,:)),'r')
set(gca,'yscale','log')
legend('still small stimulus', 'running small stimulus','running large stimulus')
axis([0,100,1,1000])
xlabel('frequency [Hz]')
ylabel('power')

% state dependence new version spontaneous
% 244 pretty little bump, 202 also big bump for still
% Suppl 1
i = 6 % example
for i = 1:length(lfpinds)
    figure
    fill(fillx,squeeze(r1tbtfillspecy(lfpinds(i),1,1,:)),'k')
    hold on
    fill(fillx,squeeze(r0tbtfillspecy(lfpinds(i),1,1,:)),[.5,.5,.5])
    plot(fax(1:104),squeeze(r1condallS(lfpinds(i),1,1,1:104)),'w')
    plot(fax(1:104),squeeze(r0condallS(lfpinds(i),1,1,1:104)),'k')
    set(gca,'yscale','log')
    legend('small running', 'small quiescent')
    axis([0,100,2,2800])
    xlabel('frequency [Hz]')
    ylabel('power')
end

% high gamma with stim and state
figure
[p,s] = signrank(r0peakhighgammapow(lfpinds(cphgpeak),1,1),r1peakhighgammapow(lfpinds(cphgpeak),1,1));
plot(1,r0peakhighgammapow(lfpinds(cphgpeak),1,1),'ko','markerfacecolor','k')
hold on
plot(2,r1peakhighgammapow(lfpinds(cphgpeak),1,1),'ko','markerfacecolor','k')
% plot(3,r1peakhighgammapow(lfpinds(cphgpeak),1,5),'ko','markerfacecolor','k')
axis([0,3,0,100])
for i = 1:length(lfpinds(cphgpeak))
    line([1,2],[r0peakhighgammapow(lfpinds(cphgpeak(i)),1,1),r1peakhighgammapow(lfpinds(cphgpeak(i)),1,1)],'color','k')
%     line([2,3],[r1peakhighgammapow(lfpinds(cphgpeak(i)),1,1),r1peakhighgammapow(lfpinds(cphgpeak(i)),1,5)],'color','k')
end
set(gca,'xticklabel',{'quiescent','running'})
set(gca,'xtick',[1,2]);
ylabel('high gamma power')
set(gcf,'OuterPosition',[573   504   242   513])
title(['n: ' int2str(length(find(~isnan(r0peakhighgammapow(lfpinds(cphgpeak),1,1))))) ' p: ' num2str(p)]);

% high gamma with state - spontaneous
plot(1,r0c0peakhighgammapow(lfpinds(cphgpeak),1),'ko','markerfacecolor','k')
hold on
plot(2,r1c0peakhighgammapow(lfpinds(cphgpeak),1),'ko','markerfacecolor','k')
axis([0,3,0,100])
for i = 1:length(lfpinds(cphgpeak))
    line([1,2],[r0c0peakhighgammapow(lfpinds(cphgpeak(i)),1),r1c0peakhighgammapow(lfpinds(cphgpeak(i)),1)],'color','k')
end
set(gca,'xticklabel',{'still spont','running spont'})
set(gca,'xtick',[1,2]);
ylabel('low gamma power')
set(gcf,'OuterPosition',[573   504   242   513])

% Fig 1
% average SOM gamma
cfreqs = chf(lgi(lfpinds)); cfreqs(isnan(r1peakgammapow(lfpinds,1,5))) = NaN; n = sum(~isnan(cfreqs));
disp(['mean center gamma frequency: ' num2str(nanmean(cfreqs)) '+-' num2str(nanstd(cfreqs)./sqrt(n))]);
disp(['abs. size gamma control: ' num2str(nanmean(r1peakgammapow(lfpinds,1,5))) '+-' num2str(nanstd(r1peakgammapow(lfpinds,1,5))./sqrt(n))]);
disp(['abs. size gamma light: ' num2str(nanmean(r1peakgammapow(lfpinds,2,5))) '+-' num2str(nanstd(r1peakgammapow(lfpinds,2,5))./sqrt(n))]);
[s,p] = ttest(r1peakgammapow(lfpinds,1,5),r1peakgammapow(lfpinds,2,5));
disp(['power reduction ttest: p = ' num2str(p)]);
[p,s] = signrank(r1peakgammapow(lfpinds,1,5),r1peakgammapow(lfpinds,2,5));
disp(['power reduction signrank: p = ' num2str(p)]);
disp(['n = ' int2str(n)])
% [p,table,stats] = kruskalwallis(squeeze(r1peakgammapow(lfpinds,1,:))) % effect of size on low gamma
% [p,table,stats] = kruskalwallis(squeeze(r1peakhighgammapow(lfpinds,1,:)))  % effect of size on high gamma
% [p,table,stats] = kruskalwallis(squeeze(r1peakhighgammapow(lfpinds(cphgpeak),1,:))) % size effect on high gamma animals with peak
% [p,table,stats] = kruskalwallis(squeeze(r1peakhighgammapow(lfpinds(cphgpeak),1,:))) % size effect on high gamma animals with peak
% [p,table,stats] = kruskalwallis(nnc0hgp(cphgpeak,:));
[p,table,stats] = kruskalwallis(nnc0gp); % effect of size on low gamma power

figure
plot(awgn(chf(lgi(lfpinds(1:16))),10),awgn(cdecf(1:16),10),'ro','markerfacecolor','r')
axis([20,35,20,35])
axis square
refline(1,0)
xlabel('center frequency size driven gamma')
ylabel('center of largest decrease band with SOM inactivation')

[r,p] = nancorr(cdecf(1:16), chf(lgi(lfpinds(1:16))));

% higher gamma band
cfreqs = chf(hgi(lfpinds)); cfreqs(isnan(r1peakhighgammapow(lfpinds,1,5))) = NaN; n = sum(~isnan(cfreqs));
% Fig 1 higher gamma frequency range
disp([' high gamma center freq = ' num2str(nanmean(cfreqs)) ' +- ' num2str(nanstd(cfreqs)./sqrt(n)) '  n = ' int2str(n)])


% gamma with size
plot(1,r1peakgammapow(lfpinds,1,1),'bo','markerfacecolor','b')
hold on
plot(2,r1peakgammapow(lfpinds,1,5),'ro','markerfacecolor','r')
axis([0,3,0,260])
for i = 1:length(lfpinds)
    line([1,2],[r1peakgammapow(lfpinds(i),1,1),r1peakgammapow(lfpinds(i),1,5)],'color','k')
end
set(gca,'xticklabel',{'small','large'})
set(gca,'xtick',[1,2]);
ylabel('low gamma power')
set(gcf,'OuterPosition',[573   504   242   513])

%all animals
for i = 1:length(lfpinds)
    figure
    semilogy(chf(1:150),squeeze(condrunS(lfpinds(i),1,1,:)),'b')
    hold on
    semilogy(chf(1:150),squeeze(condrunS(lfpinds(i),1,2,:)),'c')
    semilogy(chf(1:150),squeeze(condrunS(lfpinds(i),1,3,:)),'g')
    semilogy(chf(1:150),squeeze(condrunS(lfpinds(i),1,4,:)),'m')
    semilogy(chf(1:150),squeeze(condrunS(lfpinds(i),1,5,:)),'r')
    semilogy(chf(1:150),r1l0contS(lfpinds(i),:),'k')
    title(cllname{lfpinds(i)})
end


%all animals high gamma running dependence
for i = 1:length(lfpinds)
    figure
    semilogy(chf(1:150),squeeze(condstillS(lfpinds(i),1,1,:)),'b')
    hold on
    semilogy(chf(1:150),squeeze(condrunS(lfpinds(i),1,1,:)),'r')
    plot(chf(hgi(lfpinds(i))),squeeze(condrunS(lfpinds(i),1,1,hgi(lfpinds(i)))),'ro','markerfacecolor','r')
%     semilogy(chf(1:150),r1l0contS(lfpinds(i),:),'m')
%     semilogy(chf(1:150),r0l0contS(lfpinds(i),:),'c')
end

figure
plot(r1peakgammapow(lfpinds,1,1),r1peakgammapow(lfpinds,1,5),'ko','markerfacecolor','k')
line([0,250],[0,250],'color','k')
axis square
xlabel('gamma power small stimulus')
ylabel('gamma power large stimulus')

% FIGURE 2b
% SOM or PV Halo
[p,s] = signrank(r1peakgammapow(lfpinds(17:34),1,5),r1peakgammapow(lfpinds(17:34),2,5));
plot(1,r1peakgammapow(lfpinds(17:34),1,5),'ko','markerfacecolor','k')
hold on
plot(2,r1peakgammapow(lfpinds(17:34),2,5),'ro','markerfacecolor','r')
% axis([0,3,0,550])
axis([0,3,0,900])
% for i = lfpinds(17:34)
    plot([1,2],[r1peakgammapow(lfpinds(17:34),1,5),r1peakgammapow(lfpinds(17:34),2,5)],'color','k')
% end
set(gca,'xticklabel',{'control','PV Halo'})
% set(gca,'xticklabel',{'control','SOM Halo'})
set(gca,'xtick',[1,2]);
ylabel('low gamma power')
title(['n = ' int2str(length(find(~isnan(r1peakgammapow(lfpinds(17:34)))))) ' p: ' num2str(p)]);
set(gcf,'OuterPosition',[573   504   242   513])

% size dependence plot
% %som
% li = 9;
%PV
% li = 5;
% example 
li = 20;

for li = 1:length(lfpinds)
figure
fill(fillx,squeeze(r1c0tbtfillspecy(lfpinds(li),1,:)),[.9,.9,.9]);
hold on
fill(fillx,squeeze(r1tbtfillspecy(lfpinds(li),1,1,:)),[.7,.7,.7])
fill(fillx,squeeze(r1tbtfillspecy(lfpinds(li),1,3,:)),[.4,.4,.4])
fill(fillx,squeeze(r1tbtfillspecy(lfpinds(li),1,5,:)),[0, 0, 0])
plot(chf(12:103), r1l0contallS(lfpinds(li),12:103),'k');
plot(chf(12:103), squeeze(r1condallS(lfpinds(li),1,1,12:103)),'color',[.4,.4,.4])
plot(chf(12:103), squeeze(r1condallS(lfpinds(li),1,3,12:103)),'color',[.7,.7,.7])
plot(chf(12:103), squeeze(r1condallS(lfpinds(li),1,5,12:103)),'color',[.9,.9,.9])

set(gca,'yscale','log')
axis([0,100,1,500])
legend('no stimulus','8 degrees','21 degrees','60 degrees')
xlabel('frequency')
ylabel('power')
title(cllname{lfpinds(li)})
end

% inset
plot(sizes,squeeze(condrunS(lfpinds(li),1,:,lgi(lfpinds(li)))),'k','linewidth',2)
hold on
for i = 1:5
    line([sizes(i),sizes(i)],[squeeze(condrunSerr(lfpinds(li),1,i,1,lgi(lfpinds(li)))),squeeze(condrunSerr(lfpinds(li),1,i,2,lgi(lfpinds(li))))],'color','k','linewidth',2)
end
axis([0,70,20,130])
xlabel('size [dva]')
ylabel('low gamma power')

%SOM dependence plot
for i = 1:length(lfpinds)
    figure
    fill(fillx,squeeze(r1tbtfillspecy(lfpinds(i),1,5,:)),[.3,.3,1])
    hold on
    fill(fillx,squeeze(r1tbtfillspecy(lfpinds(i),2,5,:)),[1,.3,.3])
    set(gca,'yscale','log')
    axis([0,100,2,500])
    legend('no light','SOM Halo')
    xlabel('frequency')
    ylabel('power')
end

% change in phase locking
locked = r1lockpval(:,1,5)<0.05;

% PLV
figure
plot(r1orimeanplv(l23rs&okn',1,5),r1orimeanplv(l23rs&okn',2,5),'ko');%,'markerfacecolor','k')
hold on
plot(r1orimeanplv(l23rs&okn'&locked',1,5),r1orimeanplv(l23rs&okn'&locked',2,5),'ko','markerfacecolor','k')
plot(r1orimeanplv(l23fs&okn',1,5),r1orimeanplv(l23fs&okn',2,5),'go');%,'markerfacecolor','g')
plot(r1orimeanplv(l23fs&okn'&locked',1,5),r1orimeanplv(l23fs&okn'&locked',2,5),'go','markerfacecolor','g')
axis([0,1,0,1]);
axis square;
line([0,1],[0,1],'color','k')
xlabel('phase locking value control')
ylabel('phase locking value SOM Halo')

n = sum(~isnan(r1orimeanplv(l23rs&okn'&locked',2,5)));
disp(['RS PLV control: ' num2str(nanmean(r1orimeanplv(l23rs&okn'&locked',1,5))) '+-' num2str(nanmean(r1orimeanplv(l23rs&okn'&locked',1,5))./sqrt(n))]);
disp(['RS PLV SOM Halo: ' num2str(nanmean(r1orimeanplv(l23rs&okn'&locked',2,5))) '+-' num2str(nanmean(r1orimeanplv(l23rs&okn'&locked',2,5))./sqrt(n))]);
[p,s] = signrank(r1orimeanplv(l23rs&okn'&locked',1,5),r1orimeanplv(l23rs&okn'&locked',2,5));
disp(['signrank plv control vs SOM Halo: p = ' num2str(p)]);
disp(['n = ' int2str(n)]);

n = sum(~isnan(r1orimeanplv(l23fs&okn'&locked',2,5)));
disp(['FS PLV control: ' num2str(nanmean(r1orimeanplv(l23fs&okn'&locked',1,5))) '+-' num2str(nanmean(r1orimeanplv(l23fs&okn'&locked',1,5))./sqrt(n))]);
disp(['FS PLV SOM Halo: ' num2str(nanmean(r1orimeanplv(l23fs&okn'&locked',2,5))) '+-' num2str(nanmean(r1orimeanplv(l23fs&okn'&locked',2,5))./sqrt(n))]);
[p,s] = signrank(r1orimeanplv(l23fs&okn'&locked',1,5),r1orimeanplv(l23fs&okn'&locked',2,5));
disp(['signrank plv control vs SOM Halo: p = ' num2str(p)]);
disp(['n = ' int2str(n)]);

% PPC
figure
plot(r1orimeanppc(l23rs&okn'&locked',1,5),r1orimeanppc(l23rs&okn'&locked',2,5),'ko','markerfacecolor','k','markersize',6)
hold on
plot(r1orimeanppc(l23fs&okn'&locked',1,5),r1orimeanppc(l23fs&okn'&locked',2,5),'go','markerfacecolor','g','markersize',6)
% plot(r1orimeanppc(l23rs&okn'&~locked',1,5),r1orimeanppc(l23rs&okn'&~locked',2,5),'ko','markerfacecolor','k','markersize',3)
% plot(r1orimeanppc(l23fs&okn'&~locked',1,5),r1orimeanppc(l23fs&okn'&~locked',2,5),'go','markerfacecolor','g','markersize',3)
axis([-.1,.6,-.1,.6]);
axis square;
line([-.1,.6],[-.1,.6],'color','k')
xlabel('pairwise phase consistency control')
ylabel('pairwise phase consistency PV Halo')
legend('RS cells','FS cells','location','nw')
% legend('RS cells - locked','FS cells - locked','RS cells - not locked','FS cells - not locked','location','nw')

n = sum(~isnan(r1orimeanppc(l23rs&okn'&locked',2,5)));
disp(['RS PPC control: ' num2str(nanmean(r1orimeanppc(l23rs&okn'&locked',1,5))) '+-' num2str(nanmean(r1orimeanppc(l23rs&okn'&locked',1,5))./sqrt(n))]);
disp(['RS PPC SOM Halo: ' num2str(nanmean(r1orimeanppc(l23rs&okn'&locked',2,5))) '+-' num2str(nanmean(r1orimeanppc(l23rs&okn'&locked',2,5))./sqrt(n))]);
[p,s] = signrank(r1orimeanppc(l23rs&okn'&locked',1,5),r1orimeanppc(l23rs&okn'&locked',2,5));
disp(['RS signrank ppc control vs SOM Halo: p = ' num2str(p)]);
disp(['n = ' int2str(n)])

n = sum(~isnan(r1orimeanppc(l23fs&okn'&locked',2,5)));
disp(['FS PPC control: ' num2str(nanmean(r1orimeanppc(l23fs&okn'&locked',1,5))) '+-' num2str(nanmean(r1orimeanppc(l23fs&okn'&locked',1,5))./sqrt(n))]);
disp(['FS PPC SOM Halo: ' num2str(nanmean(r1orimeanppc(l23fs&okn'&locked',2,5))) '+-' num2str(nanmean(r1orimeanppc(l23fs&okn'&locked',2,5))./sqrt(n))]);
[p,s] = signrank(r1orimeanppc(l23fs&okn'&locked',1,5),r1orimeanppc(l23fs&okn'&locked',2,5));
disp(['FS signrank ppc control vs SOM Halo: p = ' num2str(p)]);
disp(['n = ' int2str(n)])

% circular standard deviation
figure
plot(cstd(l23rs&okn'&locked',1,5),cstd(l23rs&okn'&locked',2,5),'ko','markerfacecolor','k','markersize',6)
hold on
plot(cstd(l23fs&okn'&locked',1,5),cstd(l23fs&okn'&locked',2,5),'go','markerfacecolor','g','markersize',6)
plot(cstd(l23rs&okn'&~locked',1,5),cstd(l23rs&okn'&~locked',2,5),'ko','markerfacecolor','k','markersize',3)
plot(cstd(l23fs&okn'&~locked',1,5),cstd(l23fs&okn'&~locked',2,5),'go','markerfacecolor','g','markersize',3)
axis([.7,3,.7,3]);
axis square;
line([.7,3],[.7,3],'color','k')
xlabel('circular standard deviation control')
ylabel('circular standard deviation SOM Halo')
legend('RS cells - locked','FS cells - locked','RS cells - not locked','FS cells - not locked','location','se')

n = sum(~isnan(cstd(l23rs&okn'&locked',2,5)));
disp(['RS cstd control: ' num2str(nanmean(cstd(l23rs&okn'&locked',1,5))) '+-' num2str(nanmean(cstd(l23rs&okn'&locked',1,5))./sqrt(n))]);
disp(['RS cstd SOM Halo: ' num2str(nanmean(cstd(l23rs&okn'&locked',2,5))) '+-' num2str(nanmean(cstd(l23rs&okn'&locked',2,5))./sqrt(n))]);
[p,s] = signrank(cstd(l23rs&okn'&locked',1,5),cstd(l23rs&okn'&locked',2,5));
disp(['RS signrank cstd control vs SOM Halo: p = ' num2str(p)]);
disp(['n = ' int2str(n)])

n = sum(~isnan(cstd(l23fs&okn'&locked',2,5)));
disp(['FS cstd control: ' num2str(nanmean(cstd(l23fs&okn'&locked',1,5))) '+-' num2str(nanmean(cstd(l23fs&okn'&locked',1,5))./sqrt(n))]);
disp(['FS cstd SOM Halo: ' num2str(nanmean(cstd(l23fs&okn'&locked',2,5))) '+-' num2str(nanmean(cstd(l23fs&okn'&locked',2,5))./sqrt(n))]);
[p,s] = signrank(cstd(l23fs&okn'&locked',1,5),cstd(l23fs&okn'&locked',2,5));
disp(['FS signrank cstd control vs SOM Halo: p = ' num2str(p)]);
disp(['n = ' int2str(n)])

for i = 1:length(lfpinds)
    figure
    subplot(2,2,1)
    fill(fillx,fillycontSr0l0(lfpinds(i),:),'b')
    hold on
    fill(fillx,fillycontSr1l0(lfpinds(i),:),'c')
    fill(fillx,fillycontSr0l1(lfpinds(i),:),'r')
    fill(fillx,fillycontSr1l1(lfpinds(i),:),'m')
    set(gca,'yscale','log')
    legend('R0 L0', 'R1 L0', 'R0 L1', 'R1 L1')
    title('spontaneous')    
    
    subplot(2,2,2)
    fill(fillx,squeeze(fillspecy(lfpinds(i),1,1,:)),'b')
    hold on
    fill(fillx,squeeze(fillspecy(lfpinds(i),1,5,:)),'c')
    fill(fillx,squeeze(fillspecy(lfpinds(i),2,1,:)),'r')
    fill(fillx,squeeze(fillspecy(lfpinds(i),2,5,:)),'m')
    set(gca,'yscale','log')
    legend('small L0', 'large L0', 'small L1', 'large L1')
    title('all')    
    ax = axis;
    line([fillx(lgi(lfpinds(i))),fillx(lgi(lfpinds(i)))],[ax(3),ax(4)],'color','k')
    
    subplot(2,2,3)
    fill(fillx,squeeze(r0fillspecy(lfpinds(i),1,1,:)),'b')
    hold on
    fill(fillx,squeeze(r0fillspecy(lfpinds(i),1,5,:)),'c')
    fill(fillx,squeeze(r0fillspecy(lfpinds(i),2,1,:)),'r')
    fill(fillx,squeeze(r0fillspecy(lfpinds(i),2,5,:)),'m')
    set(gca,'yscale','log')
    legend('small L0', 'large L0', 'small L1', 'large L1')
    title('still')    
    ax = axis;
    line([fillx(lgi(lfpinds(i))),fillx(lgi(lfpinds(i)))],[ax(3),ax(4)],'color','k')
    
    subplot(2,2,4)
    fill(fillx,squeeze(r1fillspecy(lfpinds(i),1,1,:)),'b')
    hold on
    fill(fillx,squeeze(r1fillspecy(lfpinds(i),1,5,:)),'c')
    fill(fillx,squeeze(r1fillspecy(lfpinds(i),2,1,:)),'r')
    fill(fillx,squeeze(r1fillspecy(lfpinds(i),2,5,:)),'m')
    set(gca,'yscale','log')    
    legend('small L0', 'large L0', 'small L1', 'large L1')
    title('running')    
    ax = axis;
    line([fillx(lgi(lfpinds(i))),fillx(lgi(lfpinds(i)))],[ax(3),ax(4)],'color','k')
    
end

% for high gamma running
fill(fillx,fillycontSr0l0(lfpinds(1),:),'b')
hold on
fill(fillx,fillycontSr1l0(lfpinds(1),:),'c')

%rs cells
figure
bar([squeeze(nanmean(peakgammasfc(prsv,1,1))),squeeze(nanmean(peakgammasfc(prsv,1,5))),squeeze(nanmean(peakgammasfc(prsv,2,5)))],'facecolor','k')
hold on
errorbar([squeeze(nanmean(peakgammasfc(prsv,1,1))),squeeze(nanmean(peakgammasfc(prsv,1,5))),squeeze(nanmean(peakgammasfc(prsv,2,5)))],[squeeze(nanstd(peakgammasfc(prsv,1,1)))./sqrt(length(find(prsv))),squeeze(nanstd(peakgammasfc(prsv,1,5)))./sqrt(length(find(prsv))),squeeze(nanstd(peakgammasfc(prsv,2,5)))./sqrt(length(find(prsv)))],'k.')

% inhibitory neurons
bars = [squeeze(nanmean(peakgammasfc(pfsv,1,1))),...
    squeeze(nanmean(peakgammasfc(pfsv,1,5))),...
    squeeze(nanmean(peakgammasfc(pfsv,2,5)));...
    squeeze(nanmean(peakgammasfc(phe,1,1))),...
    squeeze(nanmean(peakgammasfc(phe,1,5))),...
    squeeze(nanmean(peakgammasfc(phe,2,5)))]';
errorbars = [squeeze(nanstd(peakgammasfc(pfsv,1,1)))./sqrt(sum(pfsv)),...
    squeeze(nanstd(peakgammasfc(pfsv,1,5)))./sqrt(sum(pfsv)),...
    squeeze(nanstd(peakgammasfc(pfsv,2,5)))./sqrt(sum(pfsv));...
    squeeze(nanstd(peakgammasfc(phe,1,1)))./sqrt(sum(phe)),...
    squeeze(nanstd(peakgammasfc(phe,1,5)))./sqrt(sum(phe)),...
    squeeze(nanstd(peakgammasfc(phe,2,5)))./sqrt(sum(phe))]';



fillspecx = [fax(1:104)',fliplr(fax(1:104)')];
figure
for i = 1:length(depth)
    clf
    fill(fillx,squeeze(filly(i,1,5,:)),[.3,.3,1]); %'FaceAlpha',.5)
    hold on
    plot(chfx,squeeze(C(i,1,5,:)),'color',[.2,.2,1],'linewidth',2);
    fill(fillx,squeeze(filly(i,1,1,:)),[.8,.8,1]); %,'FaceAlpha',.5)
    plot(chfx,squeeze(C(i,1,1,:)),'color',[.7,.7,1],'linewidth',2);
% plot(sfx,squeeze(condS(cell1,1,1,:)),'color',[.7,.7,1])
% plot(sfx,squeeze(condS(cell1,1,3,:)),'color',[.4,.4,1])
% plot(sfx,squeeze(condS(cell1,1,5,:)),'color',[.2,.2,1])
% set(gca,'yscale','log')
% set(gca,'yscale','lin')
axis([0,120,-.05,.35])
title(int2str(i))
pause
end

figure
for i = 1:length(depth)
    clf
    fill(fillx,squeeze(filly(i,1,5,:)),[.3,.3,1]); %'FaceAlpha',.5)
    hold on
    plot(chfx,squeeze(C(i,1,5,:)),'color',[.2,.2,1],'linewidth',2);
    fill(fillx,squeeze(filly(i,2,5,:)),[1,.3,.3]); %,'FaceAlpha',.5)
    plot(chfx,squeeze(C(i,2,5,:)),'color',[1,.2,.2],'linewidth',2);
% plot(sfx,squeeze(condS(cell1,1,1,:)),'color',[.7,.7,1])
% plot(sfx,squeeze(condS(cell1,1,3,:)),'color',[.4,.4,1])
% plot(sfx,squeeze(condS(cell1,1,5,:)),'color',[.2,.2,1])
% set(gca,'yscale','log')
% set(gca,'yscale','lin')
axis([0,120,-.05,.35])
title(int2str(i))
pause
end

for i = 1:length(lfpinds)
    lgpow(i,:) = squeeze(condS(lfpinds(i),1,:,lgi(lfpinds(i))));
    ntmax(i,:) = squeeze(condS(lfpinds(i),1,:,lgi(lfpinds(i))))./max(squeeze(condS(lfpinds(i),1,:,lgi(lfpinds(i)))));
    ntmin(i,:) = squeeze(condS(lfpinds(i),1,:,lgi(lfpinds(i))))./min(squeeze(condS(lfpinds(i),1,:,lgi(lfpinds(i)))));
end

for i = 1:length(depth)
    for l = 1:2
        for sz = 1:5
            phases = [];
            for ori = 1:8
                phases = [phases;condg1phases{i,l,ori,sz}];
            end
            allphases{i,l,sz} = phases;
        end
    end
end

% bta = bta-300;
% % plot population overview
% for i = 1:length(depth)
%     figure('Position',[122  138  1329  804])
%     
%     subplot(2,3,1)
%     plot(bta,squeeze(nanmean(bincondcllresp(i,1,:,1,:),3)),'color',[.5,.5,.5],'linewidth',2)
%     hold on    
%     plot(bta,squeeze(nanmean(bincondcllresp(i,1,:,5,:),3)),'k','linewidth',2)
%     plot(bta,squeeze(nanmean(bincondcllresp(i,2,:,1,:),3)),'color',[1,.8,.2],'linewidth',2)
%     plot(bta,squeeze(nanmean(bincondcllresp(i,2,:,5,:),3)),'r','linewidth',2)
%     ax = axis;
%     axis([-300,2700,ax(3),ax(4)]);
%     line([0,0],[ax(3),ax(4)],'color','k')
%     line([2000,2000],[ax(3),ax(4)],'color','k')
%     line([500,500],[ax(3),ax(4)],'color','r')
%     line([1500,1500],[ax(3),ax(4)],'color','r')
%     xlabel('time [ms]')
%     ylabel('firing rate [Hz]')
%     if find(prs == i) cellstr = 'RS'; else cellstr = 'FS'; end
%     title(['depth: ' int2str(depth(i)) '  '  cellstr '  swidth: ' num2str(swidthms(i)) 'ms'])
%     
%     subplot(2,3,2)
%     errorbar([0,sizes],[controlfr(i,1),squeeze(nanmean(condfr(i,1,:,:),3))'],[controlerr(i,1),squeeze(nanmean(conderr(i,1,:,:),3))'],'k','linewidth',2)
%     hold on
%     errorbar([0,sizes],[controlfr(i,1),squeeze(nanmean(condfr(i,2,:,:),3))'],[controlerr(i,1),squeeze(nanmean(conderr(i,2,:,:),3))'],'r','linewidth',2)
%     ax = axis;
%     axis([-5,65,ax(3),ax(4)]);
%     xlabel('size')
%     ylabel('firing rate [Hz]')
%     title(['size tuning   ' cllname{i}])
%     
%     ps(i) = find(mean(condfr(i,1,:,:),3) == max(mean(condfr(i,1,:,:),3)),1,'last');
%     subplot(2,3,3)
%     errorbar(oris,squeeze(condfr(i,1,:,ps(i))),squeeze(conderr(i,1,:,ps(i))),'k','linewidth',2)
%     hold on
%     errorbar(oris,squeeze(condfr(i,2,:,ps(i))),squeeze(conderr(i,2,:,ps(i))),'r','linewidth',2)
%     ax = axis;
%     axis([-5,320,ax(3),ax(4)])
%     xlabel('orientation')
%     ylabel('firing rate [Hz]')
%     title('orientation tuning')
%     
%     subplot(2,3,4)
%     semilogy(fax,squeeze(nanmean(condlfpspect(i,1,:,1,:),3)),'color',[.5,.5,.5],'linewidth',2)
%     hold on
%     semilogy(fax,squeeze(nanmean(condlfpspect(i,1,:,5,:),3)),'k','linewidth',2)
%     semilogy(fax,squeeze(nanmean(condlfpspect(i,2,:,1,:),3)),'color',[1,.8,.2],'linewidth',2)
%     semilogy(fax,squeeze(nanmean(condlfpspect(i,2,:,5,:),3)),'r','linewidth',2)
%     xlabel('frequency [Hz]')
%     ylabel('power')
%     ax = axis;
%     axis([0,120,ax(3),ax(4)])
%     legend('small','large','small L1','large L1')
%     title('power spectral density')
%         
%     subplot(2,3,5)
%     plot(-200:200,squeeze(condstalfp(i,1,1,:)),'color',[.5,.5,.5],'linewidth',2)
%     hold on
%     plot(-200:200,squeeze(condstalfp(i,1,5,:)),'k','linewidth',2)
%     plot(-200:200,squeeze(condstalfp(i,2,1,:)),'color',[1,.8,.2],'linewidth',2)
%     plot(-200:200,squeeze(condstalfp(i,2,5,:)),'r','linewidth',2)
%     xlabel('time lag [ms]')
%     ylabel('amplitude')
%     title('spike triggered LFP average')
%     
%     
%     subplot(2,3,6)
%     plot(chfx,squeeze(C(i,1,1,:)),'color',[.5,.5,.5],'linewidth',2)
%     hold on
%     plot(chfx,squeeze(C(i,1,5,:)),'k','linewidth',2)
%     plot(chfx,squeeze(C(i,2,1,:)),'color',[1,.8,.2],'linewidth',2)
%     plot(chfx,squeeze(C(i,2,5,:)),'r','linewidth',2)
%     axis([0,120,0,.35])
%     xlabel('frequency [Hz]')
%     ylabel('coherence');
%     title('spike field coherence')
%         
%     figSize = [30 21];
%     set(gcf,'paperunits','centimeters','papersize',figSize,'paperposition',[0 0 figSize])
%     if i<10, printi = ['0', int2str(i)]; else printi = int2str(i); end
%     print([oscillprintpath ,  printi '__' cllname{i} '.pdf'],'-dpdf')
%     
% end
%     

% orientation tuning with size
for i = 1:length(depth)
    for l = 1:2
        for s = 1:5
            [condoprefratio(i,l,s),conddprefratio(i,l,s),condprefori(i,l,s),condmeanori(i,l,s),condosi(i,l,s),condmeandir(i,l,s),conddsi(i,l,s)] = getOSI(squeeze(condfr(i,l,:,s))',oris);
            for gp = 1:8
                [cpr,cdr,cpo,cmo,condg1osi(i,l,s,gp),cmd,condg1dsi(i,l,s,gp)] = getOSI(squeeze(condvstfrg1(i,gp,l,:,s))',oris);
                [cpr,cdr,cpo,cmo,condeng1osi(i,l,s,gp),cmd,condeng1dsi(i,l,s,gp)] = getOSI(squeeze(condenvstfrg1(i,gp,l,:,s))',oris);
                [cpr,cdr,cpo,cmo,condesg1osi(i,l,s,gp),cmd,condesg1dsi(i,l,s,gp)] = getOSI(squeeze(condesvstfrg1(i,gp,l,:,s))',oris);
                [cpr,cdr,cpo,cmo,condenlcg1osi(i,l,s,gp),cmd,condenlcg1dsi(i,l,s,gp)] = getOSI(squeeze(condenlcvstfrg1(i,gp,l,:,s))',oris);
            end
        end
    end
end
rpax = -pi+pi/8:pi/4:pi-pi/8;

for s = 1:5
    pclasg123rs(s) = length(find(l23rs'&squeeze(g1lockpval(:,1,s))<0.05))./length(find(l23rs));
    pclasg123fs(s) = length(find(l23fs'&squeeze(g1lockpval(:,1,s))<0.05))./length(find(l23fs));
    pclasg14fs(s) = length(find(l4fs'&squeeze(g1lockpval(:,1,s))<0.05))./length(find(l4fs));
    pclasg14rs(s) = length(find(l4rs'&squeeze(g1lockpval(:,1,s))<0.05))./length(find(l4rs));
    pclasg15fs(s) = length(find(l5fs'&squeeze(g1lockpval(:,1,s))<0.05))./length(find(l5fs));
    pclasg15rs(s) = length(find(l5rs'&squeeze(g1lockpval(:,1,s))<0.05))./length(find(l5rs));
    pclasg223rs(s) = length(find(l23rs'&squeeze(g2lockpval(:,1,s))<0.05))./length(find(l23rs));
    pclasg223fs(s) = length(find(l23fs'&squeeze(g2lockpval(:,1,s))<0.05))./length(find(l23fs));
    pclasg24fs(s) = length(find(l4fs'&squeeze(g2lockpval(:,1,s))<0.05))./length(find(l4fs));
    pclasg24rs(s) = length(find(l4rs'&squeeze(g2lockpval(:,1,s))<0.05))./length(find(l4rs));
    pclasg25fs(s) = length(find(l5fs'&squeeze(g2lockpval(:,1,s))<0.05))./length(find(l5fs));
    pclasg25rs(s) = length(find(l5rs'&squeeze(g2lockpval(:,1,s))<0.05))./length(find(l5rs));
    
    pclbsg123rs(s) = length(find(l23rs'&squeeze(g1bslocksig(:,1,s))))./length(find(l23rs));
    pclbsg123fs(s) = length(find(l23fs'&squeeze(g1bslocksig(:,1,s))))./length(find(l23fs));
    pclbsg14fs(s) = length(find(l4fs'&squeeze(g1bslocksig(:,1,s))))./length(find(l4fs));
    pclbsg14rs(s) = length(find(l4rs'&squeeze(g1bslocksig(:,1,s))))./length(find(l4rs));
    pclbsg15fs(s) = length(find(l5fs'&squeeze(g1bslocksig(:,1,s))))./length(find(l5fs));
    pclbsg15rs(s) = length(find(l5rs'&squeeze(g1bslocksig(:,1,s))))./length(find(l5rs));
    pclbsg223rs(s) = length(find(l23rs'&squeeze(g2bslocksig(:,1,s))))./length(find(l23rs));
    pclbsg223fs(s) = length(find(l23fs'&squeeze(g2bslocksig(:,1,s))))./length(find(l23fs));
    pclbsg24fs(s) = length(find(l4fs'&squeeze(g2bslocksig(:,1,s))))./length(find(l4fs));
    pclbsg24rs(s) = length(find(l4rs'&squeeze(g2bslocksig(:,1,s))))./length(find(l4rs));
    pclbsg25fs(s) = length(find(l5fs'&squeeze(g2bslocksig(:,1,s))))./length(find(l5fs));
    pclbsg25rs(s) = length(find(l5rs'&squeeze(g2bslocksig(:,1,s))))./length(find(l5rs));
    
    r1pclasg123rs(s) = length(find(l23rs'&squeeze(r1lockpval(:,1,s))<0.05))./length(find(l23rs));
    r1pclasg123fs(s) = length(find(l23fs'&squeeze(r1lockpval(:,1,s))<0.05))./length(find(l23fs));
    r1pclasg14fs(s) = length(find(l4fs'&squeeze(r1lockpval(:,1,s))<0.05))./length(find(l4fs));
    r1pclasg14rs(s) = length(find(l4rs'&squeeze(r1lockpval(:,1,s))<0.05))./length(find(l4rs));
    r1pclasg15fs(s) = length(find(l5fs'&squeeze(r1lockpval(:,1,s))<0.05))./length(find(l5fs));
    r1pclasg15rs(s) = length(find(l5rs'&squeeze(r1lockpval(:,1,s))<0.05))./length(find(l5rs));
    
end

figure
bars = [pclasg123rs;pclasg123fs;pclasg14rs;pclasg14fs;pclasg15rs;pclasg15fs];
errorbars = zeros(6,5);
barweb(bars',errorbars',[],sizes,'percentage locked cells','sizes','percentage significantly locked',[],[],{'81 L2/3 RS','19 L2/3 FS','71 l4 RS','30 L4 FS','144 L5 RS','50 L5 FS'})

figure
bars = [r1pclasg123rs;r1pclasg123fs;r1pclasg14rs;r1pclasg14fs;r1pclasg15rs;r1pclasg15fs];
errorbars = zeros(6,5);
barweb(bars',errorbars',[],sizes,'percentage locked cells running','sizes','percentage significantly locked',[],[],{'81 L2/3 RS','19 L2/3 FS','71 l4 RS','30 L4 FS','144 L5 RS','50 L5 FS'})

figure
bars = [pclbsg123rs;pclbsg123fs;pclbsg14rs;pclbsg14fs;pclbsg15rs;pclbsg15fs];
errorbars = zeros(6,5);
barweb(bars',errorbars',[],sizes,'percentage locked cells','sizes','percentage significantly locked (bootstrap)',[],[],{'81 L2/3 RS','19 L2/3 FS','71 l4 RS','30 L4 FS','144 L5 RS','50 L5 FS'})

figure
bars = [pclasg123rs;pclasg223rs;pclasg14rs;pclasg24rs;pclasg15rs;pclasg25rs];
errorbars = zeros(6,5);
barweb(bars',errorbars',[],sizes,'percentage locked cells','sizes','percentage significantly locked',[],[],{'G1 L2/3 RS','G2 L2/3 RS','G1 L4 RS','G2 L4 RS','G1 L5 RS','G1 L5 RS'})

figure
bars = [pclbsg123rs;pclbsg223rs;pclbsg14rs;pclbsg24rs;pclbsg15rs;pclbsg25rs];
errorbars = zeros(6,5);
barweb(bars',errorbars',[],sizes,'percentage locked cells','sizes','percentage significantly locked (bootstrap)',[],[],{'G1 L2/3 RS','G2 L2/3 RS','G1 L4 RS','G2 L4 RS','G1 L5 RS','G1 L5 RS'})


figure
frs = squeeze(nanmean(condfr(:,1,:,:),3));
frbars = [nanmean(frs(l23rs,:));nanmean(frs(l23fs,:));nanmean(frs(l4rs,:));nanmean(frs(l4fs,:));nanmean(frs(l5rs,:));nanmean(frs(l5fs,:))];
barweb(frbars',errorbars',[],sizes,'firing rate per size','sizes','firing rates',[],[],{'81 L2/3 RS','19 L2/3 FS','71 l4 RS','30 L4 FS','144 L5 RS','50 L5 FS'})

s = 4;
cond = l23rs'&squeeze(g1lockpval(:,1,s))<0.05;
figure
errorbar(rpax,squeeze(nanmean(condg1osi(cond,1,s,:))),squeeze(nanstd(condg1osi(cond,1,s,:)))./sqrt(length(find(cond))));
hold on
errorbar(rpax,squeeze(nanmean(condg1osi(cond,2,s,:))),squeeze(nanstd(condg1osi(cond,2,s,:)))./sqrt(length(find(cond))),'r');
title(['very basic phase binned OSI of ' int2str(length(find(cond))) ' phase locked L2/3 RS neurons for size ' int2str(sizes(s))])
ylabel('OSI')
xlabel('phase of spikes')
legend('light off','light on')

s = 5;
cond = l23rs'&squeeze(g1lockpval(:,1,s))<0.05;
figure
errorbar(squeeze(nanmean(condesg1osi(cond,1,s,:))),squeeze(nanstd(condesg1osi(cond,1,s,:)))./sqrt(length(find(cond))));
hold on
errorbar(squeeze(nanmean(condesg1osi(cond,2,s,:))),squeeze(nanstd(condesg1osi(cond,2,s,:)))./sqrt(length(find(cond))),'r');
title(['equal spaced phase binned OSI of ' int2str(length(find(cond))) ' phase locked L2/3 RS neurons for size ' int2str(sizes(s))])
ylabel('OSI')
xlabel('distance to preferred phase')
legend('light off','light on')

s = 5;
cond = l23rs'&squeeze(g1lockpval(:,1,s))<0.05;
figure
errorbar(squeeze(nanmean(condeng1osi(cond,1,s,:))),squeeze(nanstd(condeng1osi(cond,1,s,:)))./sqrt(length(find(cond))));
hold on
errorbar(squeeze(nanmean(condeng1osi(cond,2,s,:))),squeeze(nanstd(condeng1osi(cond,2,s,:)))./sqrt(length(find(cond))),'r');
title(['equal n phase binned OSI of ' int2str(length(find(cond))) ' phase locked L2/3 RS neurons for size ' int2str(sizes(s))])
ylabel('OSI')
xlabel('distance to preferred phase')
legend('light off','light on')

cond = l5rs;
anovavec = condosi(cond,:,:); anovavec = anovavec(:);
anovavec = condosi(cond,:,:); anovavec = anovavec(:);
help = repmat(xsizes(2:6)',1,2*sum(cond))'; 
gs = help(:); %size
help = [zeros(1,sum(cond)),ones(1,sum(cond))]'; %light
gl = repmat(help,5,1);
[p,table,stats] = anovan(anovavec,{gs,gl},'model','full');
multcompare(stats)
multcompare(stats,'dimension',2)

for cll = 1:length(depth)
    for b = 1:20
        cppcl0(cll,b,:) = reshape(condppc(cll,b,1,:,:),40,1);
        cppcl1(cll,b,:) = reshape(condppc(cll,b,2,:,:),40,1);
        cbpl0(cll,b,:) = reshape(condbandpow(cll,1,:,:,b),40,1);
        cbpl1(cll,b,:) = reshape(condbandpow(cll,2,:,:,b),40,1);
    end
    cfrl0(cll,:) = reshape(condfr(cll,1,:,:),40,1);
    cfrl1(cll,:) = reshape(condfr(cll,2,:,:),40,1);
end

%which conditions are locked to what frequencies
for i = 1:494
    for j = 1:20
        for k = 1:2
            for l = 1:8
                for m = 1:5
                    condlocksig(i,j,k,l,m) = circ_rtest(condallphases{i,j,k,l,m});
                end
            end
        end
    end
end

nlph = cell(length(depth),5); lph = cell(length(depth),5);
nlphg2 = cell(length(depth),5); lphg2 = cell(length(depth),5);
nlap = cell(length(depth),20,5); lap = cell(length(depth),20,5);
for i = 1:length(depth)
    r1g1powl0(i,:) = reshape(r1g1pow(i,1,:,:),1,40);
    r1g1powl1(i,:) = reshape(r1g1pow(i,2,:,:),1,40);
    r0g1powl0(i,:) = reshape(r0g1pow(i,1,:,:),1,40);
    r0g1powl1(i,:) = reshape(r0g1pow(i,2,:,:),1,40);
    r1g2powl0(i,:) = reshape(r1g2pow(i,1,:,:),1,40);
    r1g2powl1(i,:) = reshape(r1g2pow(i,2,:,:),1,40);
    r0g2powl0(i,:) = reshape(r0g2pow(i,1,:,:),1,40);
    r0g2powl1(i,:) = reshape(r0g2pow(i,2,:,:),1,40);
    r1g1ppcl0(i,:) = reshape(r1g1ppc(i,1,:,:),1,40);
    r1g1ppcl1(i,:) = reshape(r1g1ppc(i,2,:,:),1,40);
    r0g1ppcl0(i,:) = reshape(r0g1ppc(i,1,:,:),1,40);
    r0g1ppcl1(i,:) = reshape(r0g1ppc(i,2,:,:),1,40);
    r1g2ppcl0(i,:) = reshape(r1g2ppc(i,1,:,:),1,40);
    r1g2ppcl1(i,:) = reshape(r1g2ppc(i,2,:,:),1,40);
    r0g2ppcl0(i,:) = reshape(r0g2ppc(i,1,:,:),1,40);
    r0g2ppcl1(i,:) = reshape(r0g2ppc(i,2,:,:),1,40);
    
    r0g1corr(i) = nancorr(r0g1powl0(i,:),r0g1ppcl0(i,:));
    r1g1corr(i) = nancorr(r1g1powl0(i,:),r1g1ppcl0(i,:));
    r0g2corr(i) = nancorr(r0g2powl0(i,:),r0g2ppcl0(i,:));
    r1g2corr(i) = nancorr(r1g2powl0(i,:),r1g2ppcl0(i,:));
    
    g1rl0(i,:) = reshape(condg1r(i,1,:,:),1,40);
    g1rl1(i,:) = reshape(condg1r(i,2,:,:),1,40);
    g2rl0(i,:) = reshape(condg2r(i,1,:,:),1,40);
    g2rl1(i,:) = reshape(condg2r(i,2,:,:),1,40);
    g1ppcl0(i,:) = reshape(condg1ppc(i,1,:,:),1,40);
    g1ppcl1(i,:) = reshape(condg1ppc(i,2,:,:),1,40);
    g2ppcl0(i,:) = reshape(condg2ppc(i,1,:,:),1,40);
    g2ppcl1(i,:) = reshape(condg2ppc(i,2,:,:),1,40);
    
    g1powl0(i,:) = reshape(condgpow1(i,1,:,:),1,40);
    g1powl1(i,:) = reshape(condgpow1(i,2,:,:),1,40);
    g2powl0(i,:) = reshape(condgpow2(i,1,:,:),1,40);
    g2powl1(i,:) = reshape(condgpow2(i,2,:,:),1,40);
    
    frl0(i,:) = reshape(condfr(i,1,:,:),1,40);
    frl1(i,:) = reshape(condfr(i,2,:,:),1,40);
    
    [frpowr(i), frpowp(i)] = nancorr(frl0(i,:),g1powl0(i,:));
    [frlockr(i), frlockp(i)] = nancorr(frl0(i,:),g1rl0(i,:));
    [powlockr(i), powlockp(i)] = nancorr(g1powl0(i,:),g1rl0(i,:));
    [frppcr(i), frppcp(i)] = nancorr(frl0(i,:),g1ppcl0(i,:));
    [powppcr(i), powppcp(i)] = nancorr(g1powl0(i,:),g1ppcl0(i,:));
    
    for s = 1:5
        for o = 1:8
            nlph{i,s} = [nlph{i,s}; condg1phases{i,1,o,s}];
            lph{i,s} = [lph{i,s}; condg1phases{i,2,o,s}];
            nlphg2{i,s} = [nlphg2{i,s}; condg2phases{i,1,o,s}];
            lphg2{i,s} = [lphg2{i,s}; condg2phases{i,2,o,s}];
            for b = 1:20
                nlcap{i,b,s} = [nlap{i,b,s}; condallphases{i,b,1,o,s}];
                lcap{i,b,s} = [lap{i,b,s}; condallphases{i,b,2,o,s}];
            end
        end
        sppcl0(i,s) = ppc(nlph{i,s});
        sppcl1(i,s) = ppc(lph{i,s});
        srpl0(i,s) = circ_rtest(nlph{i,s});
        srpl1(i,s) = circ_rtest(lph{i,s});
        scmeanl0(i,s) = circ_mean(nlph{i,s});
        scmeanl1(i,s) = circ_mean(lph{i,s});
        
        sppcl0g2(i,s) = ppc(nlphg2{i,s});
        sppcl1g2(i,s) = ppc(lphg2{i,s});
        srpl0g2(i,s) = circ_rtest(nlphg2{i,s});
        srpl1g2(i,s) = circ_rtest(lphg2{i,s});
        scmeanl0g2(i,s) = circ_mean(nlphg2{i,s});
        scmeanl1g2(i,s) = circ_mean(lphg2{i,s});
        
        for b = 1:20
            sappcl0(i,b,s) = ppc(nlcap{i,b,s});
            sappcl1(i,b,s) = ppc(lcap{i,b,s});
            sarpl0(i,b,s) = circ_rtest(nlcap{i,b,s});
            sarpl1(i,b,s) = circ_rtest(lcap{i,b,s});
            sacmeanl0(i,b,s) = circ_mean(nlcap{i,b,s});
            sacmeanl1(i,b,s) = circ_mean(lcap{i,b,s});
        end
    end
end
locked = mean(srpl0,2)<.05;

for i = find(l23rs)
    figure('Name',[cllname{i}, ' depth: ' int2str(depth(i)) '  swidth: ' int2str(swidth(i))])
    for s = 1:5
        subplot(2,5,s)
        circ_plot(nlph{i,s},'pretty',[],1);
%         circ_plot(nlph{i,s},'hist',[],[],1,'linewidth',2);
        title(['ppc: ' num2str(sppcl0(i,s)) '  p: ' num2str(srpl0(i,s))])
        subplot(2,5,s+5)
        circ_plot(nlphg2{i,s},'pretty','ro',1);
%         circ_plot(nlphg2{i,s},'hist',[],[],1,'linewidth',2);
        title(['ppc: ' num2str(sppcl0g2(i,s)) '  p: ' num2str(srpl0g2(i,s))])
    end
end

for i = 1:length(depth)
    for b = 1:20
        frpowcorr(i,b) = nancorr(cfrl0(i,:),cbpl0(i,b,:));
        frppccorr(i,b) = nancorr(cfrl0(i,:),squeeze(cppcl0(i,b,:)));
        powppccorr(i,b) = nancorr(squeeze(cbpl0(i,b,:)),squeeze(cppcl0(i,b,:)));
        
    end
end

%phase locking reloaded
% f = 6;
% glsml0 = squeeze(nanmean(condr(:,f,1,:,1).^2,4)); % 12 is 56:60 - around the center of the gamma peak
% glsml1 = squeeze(nanmean(condr(:,f,2,:,1).^2,4));
% glbgl0 = squeeze(nanmean(condr(:,f,1,:,5).^2,4));
% glbgl1 = squeeze(nanmean(condr(:,f,2,:,5).^2,4));

glsml0 = squeeze(nanmean(condg1ppc(:,1,:,1),3)); % 12 is 56:60 - around the center of the gamma peak
glsml1 = squeeze(nanmean(condg1ppc(:,2,:,1),3));
glbgl0 = squeeze(nanmean(condg1ppc(:,1,:,5),3));
glbgl1 = squeeze(nanmean(condg1ppc(:,2,:,5),3));

gp1sml0 = squeeze(nanmean(condgpow1(:,1,:,1),3));
gp1sml1 = squeeze(nanmean(condgpow1(:,2,:,1),3));
gp1bgl0 = squeeze(nanmean(condgpow1(:,1,:,5),3));
gp1bgl1 = squeeze(nanmean(condgpow1(:,2,:,5),3));

gp2sml0 = squeeze(nanmean(condgpow2(:,1,:,1),3));
gp2sml1 = squeeze(nanmean(condgpow2(:,2,:,1),3));
gp2bgl0 = squeeze(nanmean(condgpow2(:,1,:,5),3));
gp2bgl1 = squeeze(nanmean(condgpow2(:,2,:,5),3));

frsml0 = squeeze(nanmean(condfr(:,1,:,1),3));
frsml1 = squeeze(nanmean(condfr(:,2,:,1),3));
frbgl0 = squeeze(nanmean(condfr(:,1,:,5),3));
frbgl1 = squeeze(nanmean(condfr(:,2,:,5),3));

ta = linspace(-299,2700,3000);
cond = l23rs;
gam = condgamma1resp;
%across time
plot(ta,squeeze(nanmean(nanmean(gam(cond,1,:,1,:),3),1)),'linewidth',2)
hold on
plot(ta,squeeze(nanmean(nanmean(gam(cond,1,:,5,:),3),1)),'c','linewidth',2)
plot(ta,squeeze(nanmean(nanmean(gam(cond,2,:,1,:),3),1)),'r','linewidth',2)
plot(ta,squeeze(nanmean(nanmean(gam(cond,2,:,5,:),3),1)),'m','linewidth',2)
mn = min(min(min(nanmean(nanmean(gam(cond,:,:,:,:),1),3))));
mx = max(max(max(nanmean(nanmean(gam(cond,:,:,:,:),1),3))));
line([0,0],[mn,mx],'color','k')
line([2000,2000],[mn,mx],'color','k')
line([500,1500],[mx-.05*mx,mx-.05*mx],'color','r','linewidth',4)
axis([-300,2700,mn,mx])
legend([{['L0 small']},{['L0 large']},{['L1 small']},{['L1 large']}],'location','ne')

%across frequencies
lbgl0 = squeeze(nanmean(condppc(:,:,1,:,5).^2,4));
lbgl1 = squeeze(nanmean(condppc(:,:,2,:,5).^2,4));
figure
errorbar(3:5:98,nanmean(lbgl0(l23rs,:)),nanstd(lbgl0(l23rs,:))./sqrt(length(find(l23rs))))
hold on
errorbar(3:5:98,nanmean(lbgl0(l4rs,:)),nanstd(lbgl0(l4rs,:))./sqrt(length(find(l4rs))),'c')
errorbar(3:5:98,nanmean(lbgl0(l5rs,:)),nanstd(lbgl0(l5rs,:))./sqrt(length(find(l5rs))),'g')
errorbar(3:5:98,nanmean(lbgl1(l23rs,:)),nanstd(lbgl1(l23rs,:))./sqrt(length(find(l23rs))),'r')
errorbar(3:5:98,nanmean(lbgl1(l4rs,:)),nanstd(lbgl1(l4rs,:))./sqrt(length(find(l4rs))),'m')
errorbar(3:5:98,nanmean(lbgl1(l5rs,:)),nanstd(lbgl1(l5rs,:))./sqrt(length(find(l5rs))),'y')
legend('L0: L23 RS','L0: L4 RS','L0: L5 RS','L1: L23 RS','L1: L4 RS','L1: L5 RS')

%across cell types
cond = glsml0;
% cond = squeeze(nanmean(nanmean(condg1ppc(:,1,:,:),3),4));
bars = [nanmean(cond(l23rs)),nanmean(cond(l23fs));...
    nanmean(cond(l4rs)),nanmean(cond(l4fs));...
    nanmean(cond(l5rs)),nanmean(cond(l5fs))];
errorbars = [nanstd(cond(l23rs))./sqrt(sum(l23rs)),...
    nanstd(cond(l23fs))./sqrt(sum(l23fs));...
    nanstd(cond(l4rs))./sqrt(sum(l4rs)),...
    nanstd(cond(l4fs))./sqrt(sum(l4fs));...
    nanstd(cond(l5rs))./sqrt(sum(l5rs)),...
    nanstd(cond(l5fs))./sqrt(sum(l5fs))];
figure
barweb(bars, errorbars, [], [{'L23'},{'L4'},{'L5'}], 'Gamma phase locking', [], 'resultant vector length', [], [], [{'RS'},{'FS'}],[],'axis');

%across size
condsm = glsml0;
condbg = glbgl0;
bars = [nanmean(condsm(l23rs)),nanmean(condbg(l23rs));...
    nanmean(condsm(l4rs)),nanmean(condbg(l4rs));...
    nanmean(condsm(l5rs)),nanmean(condbg(l5rs))];
errorbars = [nanstd(condsm(l23rs))./sqrt(sum(l23rs)),...
    nanstd(condbg(l23rs))./sqrt(sum(l23rs));...
    nanstd(condsm(l4rs))./sqrt(sum(l4rs)),...
    nanstd(condbg(l4rs))./sqrt(sum(l4rs));...
    nanstd(condsm(l5rs))./sqrt(sum(l5rs)),...
    nanstd(condbg(l5rs))./sqrt(sum(l5rs))];
figure
barweb(bars, errorbars, [], [{'L23'},{'L4'},{'L5'}], 'Gamma phase locking', [], 'resultant vector length', [], [], [{'small'},{'large'}],[],'axis');

%across light
cond0 = glbgl0;
cond1 = glbgl1;
bars = [nanmean(cond0(l23rs)),nanmean(cond1(l23rs));...
    nanmean(cond0(l4rs)),nanmean(cond1(l4rs));...
    nanmean(cond0(l5rs)),nanmean(cond1(l5rs))];
errorbars = [nanstd(cond0(l23rs))./sqrt(sum(l23rs)),...
    nanstd(cond1(l23rs))./sqrt(sum(l23rs));...
    nanstd(cond0(l4rs))./sqrt(sum(l4rs)),...
    nanstd(cond1(l4rs))./sqrt(sum(l4rs));...
    nanstd(cond0(l5rs))./sqrt(sum(l5rs)),...
    nanstd(cond1(l5rs))./sqrt(sum(l5rs))];
figure
barweb(bars, errorbars, [], [{'L23'},{'L4'},{'L5'}], 'Gamma phase locking', [], 'resultant vector length', [], [], [{'L0'},{'L1'}],[],'axis');


%bursts
for i = 1:length(l0isi)
    [n1(i,:),c1(i,:)] = hist(l1isi{i}(l1isi{i}<=50),1:50);
    [n0(i,:),c0(i,:)] = hist(l0isi{i}(l0isi{i}<=50),1:50);
    nn0(i,:) = n0(i,:)./max(n0(i,:));
    nn1(i,:) = n1(i,:)./max(n1(i,:));
    cv0(i) = std(l0isi{i})/mean(l0isi{i});
    cv1(i) = std(l1isi{i})/mean(l1isi{i});
end

cond = l5rs;
figure
subplot(2,1,1)
errorbar(nanmean(nn0(cond,:)),nanstd(nn0(cond,:))./sqrt(sum(cond)))
hold on
errorbar(nanmean(nn1(cond,:)),nanstd(nn1(cond,:))./sqrt(sum(cond)),'r')

subplot(2,1,2)
[s,p] = ttest(cv0(cond),cv1(cond));
plot(cv0,cv1,'.')
a = axis;
axis([min(a),max(a),min(a),max(a)]);
refline(1,0)
axis square
xlabel('CV light off')
ylabel('CV light on')
title(['p: ' num2str(p)]);

for i = 1:size(condfr,1)
    lfrs(i,:) = reshape(condfr(i,2,:,:),1,40);
    nlfrs(i,:) = reshape(condfr(i,1,:,:),1,40);
end
nlfs = nlfrs(l4fs,:);
lfs = lfrs(l4fs,:);

for i = 1:1000
    a = randperm(size(nlfs,1));
    for j = 1:40
        fspopsl0(j) = sparseness(nlfs(a(1:10),j));
        fspopsl1(j) = sparseness(lfs(a(1:10),j));
    end
    fspsl0(i) = mean(fspopsl0);
    fspsl1(i) = mean(fspopsl1);
end


nlrs = nlfrs(prsv'&ok,:);
lrs = lfrs(prsv'&ok,:);
for i = 1:1000
    a = randperm(size(nlrs, 1));
    for j = 1:40
        rspopsl0(j) = sparseness(nlrs(a(1:40),j));
        rspopsl1(j) = sparseness(lrs(a(1:40),j));
    end
    rspsl0(i) = mean(rspopsl0);
    rspsl1(i) = mean(rspopsl1);
end


disp('')
% % gammapath = 'C:\Users\Julia\work\data\populations\size\gamma\';
% % gammapath = 'C:\Users\Julia\work\data\populations\SOM_Halo\size\gamma\';
% gammapath = 'C:\Users\Julia\work\data\populations\SOM_Halo\size\gamma\running\';
% bandx = 3:5:98;
% for i = 1:size(condbandpow,1)
%     clf;
% %     plot(bandx,squeeze(nanmean(condbandpow(i,1,:,1,:),3)),'b','linewidth',2);
% %     hold on
% %     plot(bandx,squeeze(nanmean(condbandpow(i,1,:,2,:),3)),'c','linewidth',2,'linewidth',2);
% %     plot(bandx,squeeze(nanmean(condbandpow(i,1,:,3,:),3)),'k','linewidth',2);
% %     plot(bandx,squeeze(nanmean(condbandpow(i,1,:,4,:),3)),'m','linewidth',2);
% %     plot(bandx,squeeze(nanmean(condbandpow(i,1,:,5,:),3)),'r','linewidth',2);
%     semilogy(bandx,squeeze(nanmean(runcondbandpow(i,1,:,1,:),3)),'m','linewidth',2)
%     hold on
%     semilogy(bandx,squeeze(nanmean(runcondbandpow(i,1,:,5,:),3)),'r','linewidth',2)
%     semilogy(bandx,squeeze(nanmean(stillcondbandpow(i,1,:,1,:),3)),'c','linewidth',2)
%     semilogy(bandx,squeeze(nanmean(stillcondbandpow(i,1,:,5,:),3)),'b','linewidth',2)
%     xlabel('frequency band')
%     ylabel('power')
% %     legend([{'smallest'},{''},{''},{''},{'largest'}])
%     legend([{'running small'},{'running large'},{'still small'},{'still large'}])
%     title([' i:  ' int2str(i) '  depth: ' num2str(depth(i)) '  width: ' int2str(swidth(i))])
%     figSize = [30 21];
%     set(gcf,'paperunits','centimeters','papersize',figSize,'paperposition',[0 0 figSize])
%     if i<10, printi = ['0', int2str(i)]; else printi = int2str(i); end
%     print([gammapath ,  printi '__' cellname{i} '.pdf'],'-dpdf')
% %     pause
% end

% ta = linspace(-300,2700,100);
% i = 3
% figure
% exampleplot(squeeze(bincondcellresp(i,:,:,:,:)),squeeze(conderr(i,:,:,:)),controlfr(i,:),controlerr(i,:),ta,xsizes,oris)

for rec = 1: max(recording)
    reccells = length(find(recording == rec));
    recinds = find(recording == rec);
    
    
    %choristers, solists - correlation with what the rest of pop is doing
    fan = condfr(recinds,1,:,:);
    rfan = reshape(fan,reccells,40);
    tl0 = trialfrl0(recinds,:); tl1 = trialfrl1(recinds,:); tbl = trialbl(recinds,:);
    for i = 1:reccells
        nrfan(i,:) = rfan(i,:)-min(rfan(i,:));
        nrfan(i,:) = nrfan(i,:)./max(nrfan(i,:));
        ntl0(i,:) = tl0(i,:)-min(tl0(i,:));
        ntl0(i,:) = ntl0(i,:)./max(ntl0(i,:));
        ntl1(i,:) = tl1(i,:)-min(tl1(i,:));
        ntl1(i,:) = ntl1(i,:)./max(ntl1(i,:));
        ntbl(i,:) = tbl(i,:)-min(tbl(i,:));
        ntbl(i,:) = ntbl(i,:)./max(ntbl(i,:));
    end
    for i = 1:reccells
        tmpnrfan = nrfan; tmpnrfan(i,:) = nan(1,40); %get rid of this cell for caclulation of pop avg
        [pc(i),pp(i)] = nancorr(rfan(i,:),nanmean(tmpnrfan));
        tmpl0 = ntl0; tmpl0(i,:) = nan(1,210);
        [l0pc(i),l0pp(i)] = nancorr(tl0(i,:),nanmean(tmpl0));
        tmpl1 = ntl1; tmpl1(i,:) = nan(1,210);
        [l1pc(i),l1pp(i)] = nancorr(tl1(i,:),nanmean(tmpl1));
        tmpbl = ntbl; tmpbl(i,:) = nan(1,420);
        [blpc(i),blpp(i)] = nancorr(tbl(i,:),nanmean(tmpbl));
    end
    popcorr(recinds) = pc;
    popp(recinds) = pp;
    popl0corr(recinds) = l0pc;
    popl1corr(recinds) = l1pc;
    popblcorr(recinds) = blpc;
    clear pc; clear pp; clear l0pc; clear l0pp; clear l1pc; clear l1pp; clear blpp; clear blpc;
    
    % lifetime and population kurtosis - large kurtosis means sparser code
    l0fr = condfr(recinds,1,:,:);
    l1fr = condfr(recinds,2,:,:);
    l0fr = reshape(l0fr,reccells,40);
    l1fr = reshape(l1fr,reccells,40);
    l0large = squeeze(condfr(recinds,1,:,5));
    l1large = squeeze(condfr(recinds,2,:,5));
    l0small = squeeze(condfr(recinds,1,:,1));
    l1small = squeeze(condfr(recinds,2,:,1));
    ltkurtosisl0(recinds) = kurtosis(l0fr,0,2);
    ltkurtosisl1(recinds) = kurtosis(l1fr,0,2);
    popkurtosisl0(rec,:) = kurtosis(l0fr,0,1);
    popkurtosisl1(rec,:) = kurtosis(l1fr,0,1);
    N = 40;
    sltl0(recinds) = (1 - ((1/N)* ((sum(l0fr,2).^2)./(sum(l0fr.^2,2))) )) / (1-(1/N));
    sltl1(recinds) = (1 - ((1/N)* ((sum(l1fr,2).^2)./(sum(l1fr.^2,2))) )) / (1-(1/N));
    N = 8;
    sltl0large(recinds) = (1 - ((1/N)* ((sum(l0large,2).^2)./(sum(l0large.^2,2))) )) / (1-(1/N));
    sltl1large(recinds) = (1 - ((1/N)* ((sum(l1large,2).^2)./(sum(l1large.^2,2))) )) / (1-(1/N));
    sltl0small(recinds) = (1 - ((1/N)* ((sum(l0small,2).^2)./(sum(l0small.^2,2))) )) / (1-(1/N));
    sltl1small(recinds) = (1 - ((1/N)* ((sum(l1small,2).^2)./(sum(l1small.^2,2))) )) / (1-(1/N));
    N = reccells;
    spopl0(rec,:) = (1 - ((1/N)* ((sum(l0fr,1).^2)./(sum(l0fr.^2,1))) )) / (1-(1/N));
    spopl1(rec,:) = (1 - ((1/N)* ((sum(l1fr,1).^2)./(sum(l1fr.^2,1))) )) / (1-(1/N));
    
    ii = 1;
    for i = 1:reccells-1
        for j = i+1:reccells
            for l = 1:2
                for o = 1:8
                    for s = 1:5
                        cv = cov(cellz{recinds(i),l,o,s},cellz{recinds(j),l,o,s});
                        cvs(ii,l,o,s) = cv(1,2);
                        cvcount = cov(cellsc{recinds(i),l,o,s},cellsc{recinds(j),l,o,s});
                        if size(cvcount) == 1
                            cvcounts(ii,l,o,s) = NaN;
                        else
                            cvcounts(ii,l,o,s) = cvcount(1,2);
                        end
                        cc = corrcoef(cellz{recinds(i),l,o,s},cellz{recinds(j),l,o,s});
                        if isnan(cc)
                            cc(ii,l,o,s) = NaN;
                        else
                            ccs(ii,l,o,s) = cc(1,2);
                        end
                        pairinds(ii,:) = [i,j];
                    end
                end
            end
            ii = ii+1;
        end
    end
    if ii == 1
        mscrl0{rec} = NaN; mscrl1{rec} = NaN;
        mscountrl0{rec} = NaN; mscountrl1{rec} = NaN;
        mccl0{rec} = NaN; mccl1{rec} = NaN;
        pairs{rec} = NaN; recs{rec} = recinds;
    else
        mscrl0{rec} = nanmean(nanmean(cvs(:,1,:,:),3),4); % mean spike count covariance on z scores
        mscrl1{rec} = nanmean(nanmean(cvs(:,2,:,:),3),4);
        mscountrl0{rec} = nanmean(nanmean(cvcounts(:,1,:,:),3),4); %on raw spike counts
        mscountrl1{rec} = nanmean(nanmean(cvcounts(:,2,:,:),3),4);
        mccl0{rec} = nanmean(nanmean(ccs(:,1,:,:),3),4); % correlation coefficients
        mccl1{rec} = nanmean(nanmean(ccs(:,2,:,:),3),4);
        pairs{rec} = recinds(pairinds);
        recs{rec} = recinds;
    end
    clear pairinds ccs cvcounts cvs;
%     %lifetime sparseness
%     N = 34; %1000ms binned 28:61 on bta (500:1500)
%     for cl = 1:reccells
%         for l = 1:2
%             for o = 1:8
%                 for s = 1:5
%                     slt(cl,l,o,s) = (1 - ((1/N)* ((sum(bincondcellresp(recinds(cl),l,o,s,28:61),5).^2)/(sum(bincondcellresp(recinds(cl),l,o,s,28:61).^2,5))) ))...
%                         / (1-(1/N));
%                 end
%             end
%         end
%     end
%     
%     %population sparseness
%     N = size(bincondcellresp,1);
%     for bn = 1:34
%         for l = 1:2
%             for o = 1:8
%                 for s = 1:5
%                     spop(bn,l,o,s) = (1 - ((1/N)* ((sum(bincondcellresp(recinds,l,o,s,bn+27),1).^2)/(sum(bincondcellresp(recinds,l,o,s,bn+27).^2,1))) ))...
%                         / (1-(1/N));
%                 end
%             end
%         end
%     end
%     lifetimesparseness{rec} = slt;
%     popsparseness{rec} = spop;
%     clear slt; clear spop;
    
end

for i = 1:length(depth)
    brl1(i,:) = reshape(cellbinrely(i,2,:,:),1,40);
    brl0(i,:) = reshape(cellbinrely(i,1,:,:),1,40);
end

% trial to trial reliability
bars = [nanmean(brl0(l23rs)),nanmean(brl0(l23fs));...
    nanmean(brl0(l4rs)),nanmean(brl0(l4fs));...
    nanmean(brl0(l5rs)),nanmean(brl0(l5fs))];
errorbars = [nanstd(brl0(l23rs))./sqrt(sum(l23rs)),...
    nanstd(brl0(l23fs))./sqrt(sum(l23fs));...
    nanstd(brl0(l4rs))./sqrt(sum(l4rs)),...
    nanstd(brl0(l4fs))./sqrt(sum(l4fs));...
    nanstd(brl0(l5rs))./sqrt(sum(l5rs)),...
    nanstd(brl0(l5fs))./sqrt(sum(l5fs))];
figure
barweb(bars, errorbars, [], [{'L23'},{'L4'},{'L5'}], 'Reliability light OFF', [], 'corr coef', [], [], [{'RS'},{'FS'}],[],'axis');

% trial to trial reliability
bars = [nanmean(brl0(l23rs)),nanmean(brl1(l23rs));...
    nanmean(brl0(l4rs)),nanmean(brl1(l4rs));...
    nanmean(brl0(l5rs)),nanmean(brl1(l5rs))];
errorbars = [nanstd(brl0(l23rs))./sqrt(sum(l23rs)),...
    nanstd(brl1(l23rs))./sqrt(sum(l23rs));...
    nanstd(brl0(l4rs))./sqrt(sum(l4rs)),...
    nanstd(brl1(l4rs))./sqrt(sum(l4rs));...
    nanstd(brl0(l5rs))./sqrt(sum(l5rs)),...
    nanstd(brl1(l5rs))./sqrt(sum(l5rs))];
figure
barweb(bars, errorbars, [], [{'L23'},{'L4'},{'L5'}], 'Reliability RS cells', [], 'corr coef', [], [], [{'L0'},{'L1'}],[],'axis');

% trial to trial reliability
bars = [nanmean(brl0(l23fs)),nanmean(brl1(l23fs));...
    nanmean(brl0(l4fs)),nanmean(brl1(l4fs));...
    nanmean(brl0(l5fs)),nanmean(brl1(l5fs))];
errorbars = [nanstd(brl0(l23fs))./sqrt(sum(l23fs)),...
    nanstd(brl1(l23fs))./sqrt(sum(l23fs));...
    nanstd(brl0(l4fs))./sqrt(sum(l4fs)),...
    nanstd(brl1(l4fs))./sqrt(sum(l4fs));...
    nanstd(brl0(l5fs))./sqrt(sum(l5fs)),...
    nanstd(brl1(l5fs))./sqrt(sum(l5fs))];
figure
barweb(bars, errorbars, [], [{'L23'},{'L4'},{'L5'}], 'Reliability FS cells', [], 'corr coef', [], [], [{'L0'},{'L1'}],[],'axis');



% population correlation
bars = [nanmean(popl0corr(l23rs)),nanmean(popl0corr(l23fs));...
    nanmean(popl0corr(l4rs)),nanmean(popl0corr(l4fs));...
    nanmean(popl0corr(l5rs)),nanmean(popl0corr(l5fs))];
errorbars = [nanstd(popl0corr(l23rs))./sqrt(sum(l23rs)),...
    nanstd(popl0corr(l23fs))./sqrt(sum(l23fs));...
    nanstd(popl0corr(l4rs))./sqrt(sum(l4rs)),...
    nanstd(popl0corr(l4fs))./sqrt(sum(l4fs));...
    nanstd(popl0corr(l5rs))./sqrt(sum(l5rs)),...
    nanstd(popl0corr(l5fs))./sqrt(sum(l5fs))];
figure
barweb(bars, errorbars, [], [{'L23'},{'L4'},{'L5'}], 'correlation with population light OFF', [], 'corr coef', [], [], [{'RS'},{'FS'}],[],'axis');
    
bars = [nanmean(popl0corr(l23rs)),nanmean(popl1corr(l23rs));...
    nanmean(popl0corr(l4rs)),nanmean(popl1corr(l4rs));...
    nanmean(popl0corr(l5rs)),nanmean(popl1corr(l5rs))];
errorbars = [nanstd(popl0corr(l23rs))./sqrt(sum(l23rs)),...
    nanstd(popl1corr(l23rs))./sqrt(sum(l23rs));...
    nanstd(popl0corr(l4rs))./sqrt(sum(l4rs)),...
    nanstd(popl1corr(l4rs))./sqrt(sum(l4rs));...
    nanstd(popl0corr(l5rs))./sqrt(sum(l5rs)),...
    nanstd(popl1corr(l5rs))./sqrt(sum(l5rs))];
figure
barweb(bars, errorbars, [], [{'L23'},{'L4'},{'L5'}], 'correlation with population RS cells', [], 'corr coef', [], [], [{'L0'},{'L1'}],[],'axis');
  
bars = [nanmean(popl0corr(l23fs)),nanmean(popl1corr(l23fs));...
    nanmean(popl0corr(l4fs)),nanmean(popl1corr(l4fs));...
    nanmean(popl0corr(l5fs)),nanmean(popl1corr(l5fs))];
errorbars = [nanstd(popl0corr(l23fs))./sqrt(sum(l23fs)),...
    nanstd(popl1corr(l23fs))./sqrt(sum(l23fs));...
    nanstd(popl0corr(l4fs))./sqrt(sum(l4fs)),...
    nanstd(popl1corr(l4fs))./sqrt(sum(l4fs));...
    nanstd(popl0corr(l5fs))./sqrt(sum(l5fs)),...
    nanstd(popl1corr(l5fs))./sqrt(sum(l5fs))];
figure
barweb(bars, errorbars, [], [{'L23'},{'L4'},{'L5'}], 'correlation with population FS cells', [], 'corr coef', [], [], [{'L0'},{'L1'}],[],'axis');
  

% lifetime sparseness
bars = [nanmean(sltl0(l23rs)),nanmean(sltl0(l23fs));...
    nanmean(sltl0(l4rs)),nanmean(sltl0(l4fs));...
    nanmean(sltl0(l5rs)),nanmean(sltl0(l5fs))];
errorbars = [nanstd(sltl0(l23rs))./sqrt(sum(l23rs)),...
    nanstd(sltl0(l23fs))./sqrt(sum(l23fs));...
    nanstd(sltl0(l4rs))./sqrt(sum(l4rs)),...
    nanstd(sltl0(l4fs))./sqrt(sum(l4fs));...
    nanstd(sltl0(l5rs))./sqrt(sum(l5rs)),...
    nanstd(sltl0(l5fs))./sqrt(sum(l5fs))];
figure
barweb(bars, errorbars, [], [{'L23'},{'L4'},{'L5'}], 'Lifetime sparseness light OFF', [], 'corr coef', [], [], [{'RS'},{'FS'}],[],'axis');

% lifetime sparseness
bars = [nanmean(sltl0(l23rs)),nanmean(sltl1(l23rs));...
    nanmean(sltl0(l4rs)),nanmean(sltl1(l4rs));...
    nanmean(sltl0(l5rs)),nanmean(sltl1(l5rs))];
errorbars = [nanstd(sltl0(l23rs))./sqrt(sum(l23rs)),...
    nanstd(sltl1(l23rs))./sqrt(sum(l23rs));...
    nanstd(sltl0(l4rs))./sqrt(sum(l4rs)),...
    nanstd(sltl1(l4rs))./sqrt(sum(l4rs));...
    nanstd(sltl0(l5rs))./sqrt(sum(l5rs)),...
    nanstd(sltl1(l5rs))./sqrt(sum(l5rs))];
figure
barweb(bars, errorbars, [], [{'L23'},{'L4'},{'L5'}], 'lifetime sparseness RS cells', [], 'corr coef', [], [], [{'L0'},{'L1'}],[],'axis');

% lifetime sparseness
bars = [nanmean(sltl0(l23fs)),nanmean(sltl1(l23fs));...
    nanmean(sltl0(l4fs)),nanmean(sltl1(l4fs));...
    nanmean(sltl0(l5fs)),nanmean(sltl1(l5fs))];
errorbars = [nanstd(sltl0(l23fs))./sqrt(sum(l23fs)),...
    nanstd(sltl1(l23fs))./sqrt(sum(l23fs));...
    nanstd(sltl0(l4fs))./sqrt(sum(l4fs)),...
    nanstd(sltl1(l4fs))./sqrt(sum(l4fs));...
    nanstd(sltl0(l5fs))./sqrt(sum(l5fs)),...
    nanstd(sltl1(l5fs))./sqrt(sum(l5fs))];
figure
barweb(bars, errorbars, [], [{'L23'},{'L4'},{'L5'}], 'lifetime sparseness FS cells', [], 'corr coef', [], [], [{'L0'},{'L1'}],[],'axis');

% lifetime sparseness
bars = [nanmean(sltl0small(l23rs)),nanmean(sltl0large(l23rs));...
    nanmean(sltl0small(l4rs)),nanmean(sltl0large(l4rs));...
    nanmean(sltl0small(l5rs)),nanmean(sltl0large(l5rs))];
errorbars = [nanstd(sltl0small(l23rs))./sqrt(sum(l23rs)),...
    nanstd(sltl0large(l23rs))./sqrt(sum(l23rs));...
    nanstd(sltl0small(l4rs))./sqrt(sum(l4rs)),...
    nanstd(sltl0large(l4rs))./sqrt(sum(l4rs));...
    nanstd(sltl0small(l5rs))./sqrt(sum(l5rs)),...
    nanstd(sltl0large(l5rs))./sqrt(sum(l5rs))];
figure
barweb(bars, errorbars, [], [{'L23'},{'L4'},{'L5'}], 'lifetime sparseness RS cells light OFF', [], 'sparseness', [], [], [{'small'},{'large'}],[],'axis');

% lifetime sparseness
bars = [nanmean(sltl1small(l23rs)),nanmean(sltl1large(l23rs));...
    nanmean(sltl1small(l4rs)),nanmean(sltl1large(l4rs));...
    nanmean(sltl1small(l5rs)),nanmean(sltl1large(l5rs))];
errorbars = [nanstd(sltl1small(l23rs))./sqrt(sum(l23rs)),...
    nanstd(sltl1large(l23rs))./sqrt(sum(l23rs));...
    nanstd(sltl1small(l4rs))./sqrt(sum(l4rs)),...
    nanstd(sltl1large(l4rs))./sqrt(sum(l4rs));...
    nanstd(sltl1small(l5rs))./sqrt(sum(l5rs)),...
    nanstd(sltl1large(l5rs))./sqrt(sum(l5rs))];
figure
barweb(bars, errorbars, [], [{'L23'},{'L4'},{'L5'}], 'lifetime sparseness RS cells light ON', [], 'sparseness', [], [], [{'small'},{'large'}],[],'axis');

% lifetime sparseness
bars = [nanmean(sltl0small(l23rs)),nanmean(sltl1small(l23rs));...
    nanmean(sltl0large(l23rs)),nanmean(sltl1large(l23rs))];
errorbars = [nanstd(sltl0small(l23rs))./sqrt(sum(l23rs)),...
    nanstd(sltl1small(l23rs))./sqrt(sum(l23rs));...
    nanstd(sltl0large(l23rs))./sqrt(sum(l23rs)),...
    nanstd(sltl1large(l23rs))./sqrt(sum(l23rs))];
figure
barweb(bars, errorbars, [], [{'small'},{'large'}], 'lifetime sparseness RS cells', [], 'sparseness', [], [], [{'L0'},{'L1'}],[],'axis');

cond = l23rs;
anovavec = [sltl0small(cond),sltl0large(cond),sltl1small(cond),sltl1large(cond)];
gs = [zeros(1,sum(cond)),ones(1,sum(cond)),zeros(1,sum(cond)),ones(1,sum(cond))]; %size
gl = [zeros(1,sum(cond)),zeros(1,sum(cond)),ones(1,sum(cond)),ones(1,sum(cond))]; %light
[p,table,stats] = anovan(anovavec,{gs,gl},'model','full');
multcompare(stats)
multcompare(stats,'dimension',2)


%firing rate changes
bars = [mean(nlfr(l23rs)),mean(nlfr(l23fs));...
    mean(nlfr(l4rs)),mean(nlfr(l4fs));...
    mean(nlfr(l5rs)),mean(nlfr(l5fs))];
errorbars = [std(nlfr(l23rs))./sqrt(sum(l23rs)),...
    std(nlfr(l23fs))./sqrt(sum(l23fs));...
    std(nlfr(l4rs))./sqrt(sum(l4rs)),...
    std(nlfr(l4fs))./sqrt(sum(l4fs));...
    std(nlfr(l5rs))./sqrt(sum(l5rs)),...
    std(nlfr(l5fs))./sqrt(sum(l5fs))];
figure
barweb(bars, errorbars, [], [{'L23'},{'L4'},{'L5'}], 'firing rates light OFF', [], 'corr coef', [], [], [{'RS'},{'FS'}],[],'axis');

bars = [mean(nlfr(l23rs)),mean(lfr(l23rs));...
    mean(nlfr(l4rs)),mean(lfr(l4rs));...
    mean(nlfr(l5rs)),mean(lfr(l5rs))];
errorbars = [std(nlfr(l23rs))./sqrt(sum(l23rs)),...
    std(lfr(l23rs))./sqrt(sum(l23rs));...
    std(nlfr(l4rs))./sqrt(sum(l4rs)),...
    std(lfr(l4rs))./sqrt(sum(l4rs));...
    std(nlfr(l5rs))./sqrt(sum(l5rs)),...
    std(lfr(l5rs))./sqrt(sum(l5rs))];
figure
barweb(bars, errorbars, [], [{'L23'},{'L4'},{'L5'}], 'firing rates RS cells', [], 'corr coef', [], [], [{'L0'},{'L1'}],[],'axis');

bars = [mean(nlfr(l23fs)),mean(lfr(l23fs));...
    mean(nlfr(l4fs)),mean(lfr(l4fs));...
    mean(nlfr(l5fs)),mean(lfr(l5fs))];
errorbars = [std(nlfr(l23fs))./sqrt(sum(l23fs)),...
    std(lfr(l23fs))./sqrt(sum(l23fs));...
    std(nlfr(l4fs))./sqrt(sum(l4fs)),...
    std(lfr(l4fs))./sqrt(sum(l4fs));...
    std(nlfr(l5fs))./sqrt(sum(l5fs)),...
    std(lfr(l5fs))./sqrt(sum(l5fs))];
figure
barweb(bars, errorbars, [], [{'L23'},{'L4'},{'L5'}], 'firing rates FS cells', [], 'corr coef', [], [], [{'L0'},{'L1'}],[],'axis');


figure
subplot(2,2,1)
semilogy(3:5:98, squeeze(nanmean(nanmean(runcondbandpow(l4rs,1,:,1,:),1),3)),'linewidth',2)
hold on
semilogy(3:5:98, squeeze(nanmean(nanmean(runcondbandpow(l4rs,1,:,5,:),1),3)),'c','linewidth',2)
semilogy(3:5:98, squeeze(nanmean(nanmean(stillcondbandpow(l4rs,1,:,1,:),1),3)),'r','linewidth',2)
semilogy(3:5:98, squeeze(nanmean(nanmean(stillcondbandpow(l4rs,1,:,5,:),1),3)),'m','linewidth',2)
legend([{'running small'},{'running large'},{'still small'},{'still large'}]);
title('RS L4')

subplot(2,2,2)
semilogy(3:5:98, squeeze(nanmean(nanmean(runcondbandpow(l23rs,1,:,1,:),1),3)),'linewidth',2)
hold on
semilogy(3:5:98, squeeze(nanmean(nanmean(runcondbandpow(l23rs,1,:,5,:),1),3)),'c','linewidth',2)
semilogy(3:5:98, squeeze(nanmean(nanmean(stillcondbandpow(l23rs,1,:,1,:),1),3)),'r','linewidth',2)
semilogy(3:5:98, squeeze(nanmean(nanmean(stillcondbandpow(l23rs,1,:,5,:),1),3)),'m','linewidth',2)
legend([{'running small'},{'running large'},{'still small'},{'still large'}]);
title('RS L2/3')

subplot(2,2,3)
semilogy(3:5:98, squeeze(nanmean(nanmean(runcondbandpow(l23rs,1,:,1,:),1),3)),'linewidth',2)
hold on
semilogy(3:5:98, squeeze(nanmean(nanmean(runcondbandpow(l23rs,1,:,5,:),1),3)),'c','linewidth',2)
semilogy(3:5:98, squeeze(nanmean(nanmean(stillcondbandpow(l23rs,1,:,1,:),1),3)),'r','linewidth',2)
semilogy(3:5:98, squeeze(nanmean(nanmean(stillcondbandpow(l23rs,1,:,5,:),1),3)),'m','linewidth',2)
legend([{'running small'},{'running large'},{'still small'},{'still large'}]);
title('RS L5')

figure
subplot(2,2,1)
semilogy(3:5:98, squeeze(nanmean(nanmean(condbandpow(l4rs,1,:,1,:),1),3)),'linewidth',2)
hold on
semilogy(3:5:98, squeeze(nanmean(nanmean(condbandpow(l4rs,1,:,5,:),1),3)),'c','linewidth',2)
semilogy(3:5:98, squeeze(nanmean(nanmean(condbandpow(l4rs,2,:,1,:),1),3)),'r','linewidth',2)
semilogy(3:5:98, squeeze(nanmean(nanmean(condbandpow(l4rs,2,:,5,:),1),3)),'m','linewidth',2)
legend([{'l0 small'},{'l0 large'},{'l1 small'},{'l1 large'}]);
title('RS L4')

subplot(2,2,2)
semilogy(3:5:98, squeeze(nanmean(nanmean(condbandpow(l23rs,1,:,1,:),1),3)),'linewidth',2)
hold on
semilogy(3:5:98, squeeze(nanmean(nanmean(condbandpow(l23rs,1,:,5,:),1),3)),'c','linewidth',2)
semilogy(3:5:98, squeeze(nanmean(nanmean(condbandpow(l23rs,2,:,1,:),1),3)),'r','linewidth',2)
semilogy(3:5:98, squeeze(nanmean(nanmean(condbandpow(l23rs,2,:,5,:),1),3)),'m','linewidth',2)
legend([{'l0 small'},{'l0 large'},{'l1 small'},{'l1 large'}]);
title('RS L2/3')

subplot(2,2,3)
semilogy(3:5:98, squeeze(nanmean(nanmean(condbandpow(l23rs,1,:,1,:),1),3)),'linewidth',2)
hold on
semilogy(3:5:98, squeeze(nanmean(nanmean(condbandpow(l23rs,1,:,5,:),1),3)),'c','linewidth',2)
semilogy(3:5:98, squeeze(nanmean(nanmean(condbandpow(l23rs,2,:,1,:),1),3)),'r','linewidth',2)
semilogy(3:5:98, squeeze(nanmean(nanmean(condbandpow(l23rs,2,:,5,:),1),3)),'m','linewidth',2)
legend([{'l0 small'},{'l0 large'},{'l1 small'},{'l1 large'}]);
title('RS L5')

figure
subplot(2,2,1)
semilogy(3:5:98, squeeze(nanmean(nanmean(condbandpow(l4fs,1,:,1,:),1),3)),'linewidth',2)
hold on
semilogy(3:5:98, squeeze(nanmean(nanmean(condbandpow(l4fs,1,:,5,:),1),3)),'c','linewidth',2)
semilogy(3:5:98, squeeze(nanmean(nanmean(condbandpow(l4fs,2,:,1,:),1),3)),'r','linewidth',2)
semilogy(3:5:98, squeeze(nanmean(nanmean(condbandpow(l4fs,2,:,5,:),1),3)),'m','linewidth',2)
legend([{'l0 small'},{'l0 large'},{'l1 small'},{'l1 large'}]);
title('FS L4')

subplot(2,2,2)
semilogy(3:5:98, squeeze(nanmean(nanmean(condbandpow(l23fs,1,:,1,:),1),3)),'linewidth',2)
hold on
semilogy(3:5:98, squeeze(nanmean(nanmean(condbandpow(l23fs,1,:,5,:),1),3)),'c','linewidth',2)
semilogy(3:5:98, squeeze(nanmean(nanmean(condbandpow(l23fs,2,:,1,:),1),3)),'r','linewidth',2)
semilogy(3:5:98, squeeze(nanmean(nanmean(condbandpow(l23fs,2,:,5,:),1),3)),'m','linewidth',2)
legend([{'l0 small'},{'l0 large'},{'l1 small'},{'l1 large'}]);
title('FS L2/3')

subplot(2,2,3)
semilogy(3:5:98, squeeze(nanmean(nanmean(condbandpow(l23fs,1,:,1,:),1),3)),'linewidth',2)
hold on
semilogy(3:5:98, squeeze(nanmean(nanmean(condbandpow(l23fs,1,:,5,:),1),3)),'c','linewidth',2)
semilogy(3:5:98, squeeze(nanmean(nanmean(condbandpow(l23fs,2,:,1,:),1),3)),'r','linewidth',2)
semilogy(3:5:98, squeeze(nanmean(nanmean(condbandpow(l23fs,2,:,5,:),1),3)),'m','linewidth',2)
legend([{'l0 small'},{'l0 large'},{'l1 small'},{'l1 large'}]);
title('FS L5')

respta = linspace(-299,2700,100);
normwin = find(respta>500&respta<1500);
pastcl0 = nan(size(condfr,1),9);
pastcl1 = nan(size(condfr,1),9);
for i = 1:size(condfr,1)
    prefs(i) = find(mean(condfr(i,1,:,:),3) == max(mean(condfr(i,1,:,:),3)),1,'last'); %preferred size
    aps(i) = find(mean(condfr(i,1,:,:),3) == min(mean(condfr(i,1,:,:),3)),1); % anti preferred size
    prefori(i) = find(condfr(i,1,:,prefs(i)) == max(condfr(i,1,:,prefs(i))),1);    % preferred ori
    preforil1(i) = find(condfr(i,2,:,prefs(i)) == max(condfr(i,2,:,prefs(i))),1);
    preffrl0(i) = condfr(i,1,prefori(i),prefs(i));
    preffrl1(i) = condfr(i,2,prefori(i),prefs(i));
    
    for sz = 1:5
        [orir(i,sz), dirr(i,sz), x, smeanori(i,sz), osi(i,sz), smeandir(i,sz), dsi(i,sz)] = getOSI(squeeze(condfr(i,1,:,sz))',oris);
        [orir1(i,sz), dirr1(i,sz), x, y, osi1(i,sz), c, dsi1(i,sz)] = getOSI(squeeze(condfr(i,2,:,sz))',oris); 
        sztun = squeeze(condfr(i,1,:,sz))';
        szprefori(i,sz) = find(sztun == max(sztun),1);
        
        oridiffcrv(i,sz,:) = squeeze(condfr(i,2,:,sz))-squeeze(condfr(i,1,:,sz));
        orifoldcrv(i,sz,:) = squeeze(condfr(i,2,:,sz))./squeeze(condfr(i,1,:,sz));
    end
    
    prefsizerangel0(i) = (max(condfr(i,1,prefori(i),:))-min(condfr(i,1,prefori(i),:)));
    prefsizerangel1(i) = (max(condfr(i,2,prefori(i),:))-min(condfr(i,2,prefori(i),:)));
    normprefsizerangel0(i) = (max(condfr(i,1,prefori(i),:))-min(condfr(i,1,prefori(i),:)))./max(condfr(i,1,prefori(i),:));
    normprefsizerangel1(i) = (max(condfr(i,2,prefori(i),:))-min(condfr(i,2,prefori(i),:)))./max(condfr(i,1,prefori(i),:));
    meansizerangel0(i) = (max(mean(condfr(i,1,:,:),3))-min(mean(condfr(i,1,:,:),3)));
    meansizerangel1(i) = (max(mean(condfr(i,2,:,:),3))-min(mean(condfr(i,2,:,:),3)));
    normmeansizerangel0(i) = (max(mean(condfr(i,1,:,:),3))-min(mean(condfr(i,1,:,:),3)))./max(mean(condfr(i,1,:,:),3));
    normmeansizerangel1(i) = (max(mean(condfr(i,2,:,:),3))-min(mean(condfr(i,2,:,:),3)))./max(mean(condfr(i,1,:,:),3));
    
    preforirangel0(i) = (max(condfr(i,1,:,prefs(i)))-min(condfr(i,1,:,prefs(i))));
    preforirangel1(i) = (max(condfr(i,2,:,prefs(i)))-min(condfr(i,2,:,prefs(i))));
    normpreforirangel0(i) = (max(condfr(i,1,:,prefs(i)))-min(condfr(i,1,:,prefs(i))))./max(condfr(i,1,:,prefs(i)));
    normpreforirangel1(i) = (max(condfr(i,2,:,prefs(i)))-min(condfr(i,2,:,prefs(i))))./max(condfr(i,1,:,prefs(i)));
    meanorirangel0(i) = (max(mean(condfr(i,1,:,:),4))-min(mean(condfr(i,1,:,:),4)));
    meanorirangel1(i) = (max(mean(condfr(i,2,:,:),4))-min(mean(condfr(i,2,:,:),4)));
    normmeanorirangel0(i) = (max(mean(condfr(i,1,:,:),4))-min(mean(condfr(i,1,:,:),4)))./max(mean(condfr(i,1,:,:),4));
    normmeanorirangel1(i) = (max(mean(condfr(i,2,:,:),4))-min(mean(condfr(i,2,:,:),4)))./max(mean(condfr(i,1,:,:),4));
    
    cafrl0(i,:) = squeeze(mean(condfr(i,1,:,:),3)); % average over orientations
    cafrl1(i,:) = squeeze(mean(condfr(i,2,:,:),3));   
%     carespl0(i,:,:) = squeeze(mean(bincondcellresp(i,1,:,:,:),3));
%     carespl1(i,:,:) = squeeze(mean(bincondcellresp(i,2,:,:,:),3)); 
    [cas,cai] = sort(cafrl0(i,:));    
    scafrl0(i,:) = cafrl0(i,cai); scafrl1(i,:) = cafrl1(i,cai);  % sorted average over ori fr
%     scarespl0(i,:,:) = carespl0(i,cai,:); scarespl1(i,:,:) = carespl1(i,cai,:);
%     for j = 1:5
%         nscarespl0(i,j,:) = scarespl0(i,j,:)./mean(scarespl0(i,j,normwin),3);
%         nscarespl1(i,j,:) = scarespl1(i,j,:)./mean(scarespl0(i,j,normwin),3);
%         ntomscarespl0(i,j,:) = scarespl0(i,j,:)./mean(scarespl0(i,5,normwin),3);
%         ntomscarespl1(i,j,:) = scarespl1(i,j,:)./mean(scarespl0(i,5,normwin),3);
%     end
    ncafrl0(i,:) = cafrl0(i,:)./max(cafrl0(i,:));
    ncafrl1(i,:) = cafrl1(i,:)./max(cafrl0(i,:));
    cafrangel0(i) = (max(cafrl0(i,:))-min(cafrl0(i,:)))/max(cafrl0(i,:));
    cafrangel1(i) = (max(cafrl1(i,:))-min(cafrl1(i,:)))/max(cafrl1(i,:));
    [cas,cai] = sort(ncafrl0(i,:));
    sncafrl0(i,:) = ncafrl0(i,cai); sncafrl1(i,:) = ncafrl1(i,cai);    
%     ncafrl0(i,:) = cafrl0(i,:)/max(cafrl0(i,:));
%     ncafrl1(i,:) = cafrl1(i,:)/max(cafrl0(i,:));
    siminl0(i) = max(cafrl0(i,:))/mean(cafrl0(i,:),2);
    siminl1(i) = max(cafrl1(i,:))/mean(cafrl1(i,:),2);
    
    sizediffcrv(i,:) = cafrl1(i,:)-cafrl0(i,:);
    sizefoldcrv(i,:) = cafrl1(i,:)./cafrl0(i,:);
    
    
    r1oafrl0(i,:) = squeeze(nanmean(r1condfr(i,1,:,:),3));
    r1oafrl1(i,:) = squeeze(nanmean(r1condfr(i,2,:,:),3));
    r0oafrl0(i,:) = squeeze(nanmean(r0condfr(i,1,:,:),3));
    r0oafrl1(i,:) = squeeze(nanmean(r0condfr(i,2,:,:),3));
    
    nr1oafrl0(i,:) = r1oafrl0(i,:)./max(r1oafrl0(i,:));
    nr1oafrl1(i,:) = r1oafrl1(i,:)./max(r1oafrl0(i,:));
    nr0oafrl0(i,:) = r0oafrl0(i,:)./max(r0oafrl0(i,:));
    nr0oafrl1(i,:) = r0oafrl1(i,:)./max(r0oafrl0(i,:));

    
%     condomi(i,:,:) = (condfr(i,2,:,:)-condfr(i,1,:,:))./(condfr(i,2,:,:)+condfr(i,1,:,:));
    condspikediff(i,:,:) = condfr(i,2,:,:)-condfr(i,1,:,:);
%     help = condfr(i,1,:,:); condspikesl0(i,:) = help(:);
%     help = condfr(i,2,:,:); condspikesl1(i,:) = help(:);
    condspikesl0(i,:) = reshape(condfr(i,1,:,:),1,40);
    condspikesl1(i,:) = reshape(condfr(i,2,:,:),1,40);
    condomi(i,:) = (condspikesl1(i,:)-condspikesl0(i,:))./(condspikesl1(i,:)+condspikesl0(i,:));
    normcondspikesl0(i,:) = condspikesl0(i,:)./max(condspikesl0(i,:));
    normcondspikesl1(i,:) = condspikesl1(i,:)./max(condspikesl1(i,:));
    range(i) = max(condspikesl0(i,:))-min(condspikesl0(i,:));
    rangel1(i) = max(condspikesl1(i,:))-min(condspikesl1(i,:));
    normrange(i) = (max(condspikesl0(i,:))-min(condspikesl0(i,:)))./max(condspikesl0(i,:));
    normrangel1(i) = (max(condspikesl1(i,:))-min(condspikesl1(i,:)))./max(condspikesl1(i,:));
    [r(i),p(i)] = nancorr(condspikesl0(i,:),condspikesl1(i,:)-condspikesl0(i,:));
    fitparams(i,:) = polyfit(condspikesl0(i,:),condspikesl1(i,:)-condspikesl0(i,:),1);
    normfitparams(i,:) = polyfit(normcondspikesl0(i,:),normcondspikesl1(i,:)-normcondspikesl0(i,:),1);
    
%     prefsizetunes(i,:,:) = squeeze(condfr(i,:,prefori(i),:));
%     normprefsizetunes(i,:,:) = prefsizetunes(i,:,:)./prefsizetunes(i,1,prefs(i));
    prefsizetunes(i,:,:) = [controlfr(i,:)',squeeze(condfr(i,:,prefori(i),:))];
    nprefsizetunes(i,:,:) = prefsizetunes(i,:,:)./max(prefsizetunes(i,1,:),[],3);
    [ss,ii] = sort(nprefsizetunes(i,1,:));
    nprefsizerankl0(i,:) = nprefsizetunes(i,1,ii);
    nprefsizerankl1(i,:) = nprefsizetunes(i,2,ii);
    
    meansizetunes(i,:,:) = [controlfr(i,:)', squeeze(nanmean(condfr(i,:,:,:),3))];
    nmeansizetunes(i,:,:) = meansizetunes(i,:,:)./max(meansizetunes(i,1,:),[],3);
    
%     ffl0(i,:) = reshape(cellff(i,1,:,:),1,40);
%     ffl1(i,:) = reshape(cellff(i,2,:,:),1,40);
    
    mincondl0(i) = min(condspikesl0(i,:)); mincondl1(i) = min(condspikesl1(i,:));
    maxcondl0(i) = max(condspikesl0(i,:)); maxcondl1(i) = max(condspikesl1(i,:));
    [s,in] = sort(condspikesl0(i,:));
    sortedl0(i,:) = condspikesl0(i,in);  %sorted 40 conditions per cell
    sortedl1(i,:) = condspikesl1(i,in);
    nsortedl0(i,:) = condspikesl0(i,in)./maxcondl0(i);  %sorted 40 conditions per cell
    nsortedl1(i,:) = condspikesl1(i,in)./maxcondl0(i);
    ncontrolfr(i,:) = controlfr(i,:)./maxcondl0(i);
    
%     sffl0(i,:) = ffl0(i,in); sffl1(i,:) = ffl1(i,in);
    % OMI for each size
    for j = 1:6
        sizeomi(i,j) = (sizetunel1(i,j)-sizetunel0(i,j))./(sizetunel1(i,j)+sizetunel0(i,j));
        sizedelta(i,j) = sizetunel1(i,j)-sizetunel0(i,j);
    end
    
    blssizetunel0(i,:) = sizetunel0(i,:)-sizetunel0(i,1);
    blssizetunel1(i,:) = sizetunel1(i,:)-sizetunel1(i,1);
    blssil(i) = (blssizetunel1(i,find(blssizetunel1(i,:) == max(blssizetunel1(i,:)),1))-blssizetunel1(i,end))/blssizetunel1(i,find(blssizetunel1(i,:) == max(blssizetunel1(i,:)),1));
    blssinl(i) = (blssizetunel0(i,find(blssizetunel0(i,:) == max(blssizetunel0(i,:)),1))-blssizetunel0(i,end))/blssizetunel0(i,find(blssizetunel0(i,:) == max(blssizetunel0(i,:)),1));
    blszeroedstl0(i,:) = blssizetunel0(i,:); blszeroedstl0(i,blszeroedstl0(i,:)<0) = 0; % zero below zero points in the tuning curve
    blszeroedstl1(i,:) = blssizetunel1(i,:); blszeroedstl1(i,blszeroedstl1(i,:)<0) = 0;
    blszsil(i) = (blszeroedstl1(i,find(blszeroedstl1(i,:) == max(blszeroedstl1(i,:)),1))-blszeroedstl1(i,end))/blszeroedstl1(i,find(blszeroedstl1(i,:) == max(blszeroedstl1(i,:)),1));
    blszsinl(i) = (blszeroedstl0(i,find(blszeroedstl0(i,:) == max(blszeroedstl0(i,:)),1))-blszeroedstl0(i,end))/blszeroedstl0(i,find(blszeroedstl0(i,:) == max(blssizetunel0(i,:)),1));
    
    
    %normalize average PSTHs to average of no light response period
    l1meanrespn(i,:) = l1meanresp(i,:)./mean(l0meanresp(1,normwin),2);
    l0meanrespn(i,:) = l0meanresp(i,:)./mean(l0meanresp(1,normwin),2);
    l1prefrespn(i,:) = l1prefresp(i,:)./mean(l0prefresp(1,normwin),2);
    l0prefrespn(i,:) = l0prefresp(i,:)./mean(l0prefresp(1,normwin),2);
    
    % running  
    sizetunel0r0(i,:) = squeeze(nanmean(r0condfr(i,1,:,:),3));
    sizetunel1r0(i,:) = squeeze(nanmean(r0condfr(i,2,:,:),3));
    sizetunel0r1(i,:) = squeeze(nanmean(r1condfr(i,1,:,:),3));
    sizetunel1r1(i,:) = squeeze(nanmean(r1condfr(i,2,:,:),3));
    nsizetunel0r0(i,:)= sizetunel0r0(i,:)./max(sizetunel0r1(i,:));
    nsizetunel0r1(i,:)= sizetunel0r1(i,:)./max(sizetunel0r1(i,:));
    nsizetunel1r0(i,:)= sizetunel1r0(i,:)./max(sizetunel0r1(i,:));
    nsizetunel1r1(i,:)= sizetunel1r1(i,:)./max(sizetunel0r1(i,:));
    nsizetunel1r1n2l1(i,:) = sizetunel1r1(i,:)./max(sizetunel1r1(i,:)); % normalized to light on
    
    
    
    % peak align size tuing curve and pref size during running
    if isempty(find(isnan(nanmean(r1condfr(i,1,:,:),3)),1))
        r1prefs(i) = find(nanmean(r1condfr(i,1,:,:),3) == max(nanmean(r1condfr(i,1,:,:),3)),1,'last');
        r1l1prefs(i) = find(nanmean(r1condfr(i,2,:,:),3) == max(nanmean(r1condfr(i,2,:,:),3)),1,'last');
        pastcl0(i,6-r1prefs(i):10-r1prefs(i)) = nsizetunel0r1(i,:);
        pastcl1(i,6-r1l1prefs(i):10-r1l1prefs(i)) = nsizetunel1r1n2l1(i,:);
    else
        r1prefs(i) = NaN;
        pastcl0(i,:) = nan(1,9);
        pastcl1(i,:) = nan(1,9);
    end
    
%     if ~isnan(max(mean(r0condfr(i,1,:,:),3)))
%         psr0(i) = find(mean(r0condfr(i,1,:,:),3) == max(mean(r0condfr(i,1,:,:),3)),1,'last');
%     else
%         psr0(i) = NaN;
%     end
% %     if ~isnan(
%     preforir0(i) = find(r0condfr(i,1,:,ps(i)) == max(r0condfr(i,1,:,ps(i))),1);
%     psr1(i) = find(mean(r1condfr(i,1,:,:),3) == max(mean(r1condfr(i,1,:,:),3)),1,'last');
%     preforir1(i) = find(r1condfr(i,1,:,ps(i)) == max(r1condfr(i,1,:,ps(i))),1);
     if ps(i) == 1|2
         minds = [1,2,3];
     elseif ps(i) == 3
         minds = [2,3,4];
     else
         minds = [3,4,5];
     end
     oritunel0r0(i,:) = squeeze(nanmean(r0condfr(i,1,:,minds),4));
     oritunel1r0(i,:) = squeeze(nanmean(r0condfr(i,2,:,minds),4));
     oritunel0r1(i,:) = squeeze(nanmean(r1condfr(i,1,:,minds),4));
     oritunel1r1(i,:) = squeeze(nanmean(r1condfr(i,2,:,minds),4));
     
     oritunel0(i,:) = squeeze(nanmean(condfr(i,1,:,minds),4));
     oritunel1(i,:) = squeeze(nanmean(condfr(i,2,:,minds),4));
    
end

%preferred orientations at each size
a = find(l23rs'&nanmean(osi,2)>.3);

figure
plot(smeanori(a,4),smeanori(a,5),'.')
refline(1,0)
axis square
xlabel('preferred orientation size 1')
ylabel('preferred orientation size 5')
r = nancorr(smeanori(a,1),smeanori(a,5));

for i = 1:4
    for j = i+1:5
        ccs(i,j) = nancorr(smeanori(a,i),smeanori(a,j));
    end
end

sizes = xsizes(2:end)
bars = [mean(sizes(ps(l23rs))),mean(sizes(ps(l23fs)));...
    mean(sizes(ps(l4rs))),mean(sizes(ps(l4fs)));...
    mean(sizes(ps(l5rs))),mean(sizes(ps(l5fs)))];
errorbars = [std(sizes(ps(l23rs)))./sqrt(sum(l23rs)),...
    std(sizes(ps(l23fs)))./sqrt(sum(l23fs));...
    std(sizes(ps(l4rs)))./sqrt(sum(l4rs)),...
    std(sizes(ps(l4fs)))./sqrt(sum(l4fs));...
    std(sizes(ps(l5rs)))./sqrt(sum(l5rs)),...
    std(sizes(ps(l5fs)))./sqrt(sum(l5fs))];
figure
barweb(bars, errorbars, [], [{'L23'},{'L4'},{'L5'}], 'preferred sizes', [], 'avg size', [], [], [{'RS'},{'FS'}],[],'axis');

bars = [nanmean(osi(l23rs,5)-osi(l23rs,1)),nanmean(osi1(l23rs,5)-osi1(l23rs,1));...
    nanmean(osi(l4rs,5)-osi(l4rs,1)),nanmean(osi1(l4rs,5)-osi1(l4rs,1));...
    nanmean(osi(l5rs,5)-osi(l5rs,1)),nanmean(osi1(l5rs,5)-osi1(l5rs,1))];
errorbars = [nanstd(osi(l23rs,5)-osi(l23rs,1))./sqrt(sum(l23rs)),...
    nanstd(osi1(l23rs,5)-osi1(l23rs,1))./sqrt(sum(l23rs));...
    nanstd(osi(l4rs,5)-osi(l4rs,1))./sqrt(sum(l4rs)),...
    nanstd(osi1(l4rs,5)-osi1(l4rs,1))./sqrt(sum(l4rs));...
    nanstd(osi(l5rs,5)-osi(l5rs,1))./sqrt(sum(l5rs)),...
    nanstd(osi1(l5rs,5)-osi1(l5rs,1))./sqrt(sum(l5rs))];
figure
barweb(bars, errorbars, [], [{'L23'},{'L4'},{'L5'}], 'OSI RS cells LARGE-SMALL', [], 'OSI change', [], [], [{'L0'},{'L1'}],[],'axis');

bars = [nanmean(osi(l23fs,5)-osi(l23fs,1)),nanmean(osi1(l23fs,5)-osi1(l23fs,1));...
    nanmean(osi(l4fs,5)-osi(l4fs,1)),nanmean(osi1(l4fs,5)-osi1(l4fs,1));...
    nanmean(osi(l5fs,5)-osi(l5fs,1)),nanmean(osi1(l5fs,5)-osi1(l5fs,1))];
errorbars = [nanstd(osi(l23fs,5)-osi(l23fs,1))./sqrt(sum(l23fs)),...
    nanstd(osi1(l23fs,5)-osi1(l23fs,1))./sqrt(sum(l23fs));...
    nanstd(osi(l4fs,5)-osi(l4fs,1))./sqrt(sum(l4fs)),...
    nanstd(osi1(l4fs,5)-osi1(l4fs,1))./sqrt(sum(l4fs));...
    nanstd(osi(l5fs,5)-osi(l5fs,1))./sqrt(sum(l5fs)),...
    nanstd(osi1(l5fs,5)-osi1(l5fs,1))./sqrt(sum(l5fs))];
figure
barweb(bars, errorbars, [], [{'L23'},{'L4'},{'L5'}], 'OSI FS cells LARGE-SMALL', [], 'tuning ratio', [], [], [{'L0'},{'L1'}],[],'axis');


bars = [nanmean(orir(l23fs,1)),nanmean(orir1(l23fs,1));...
    nanmean(orir(l4fs,1)),nanmean(orir1(l4fs,1));...
    nanmean(orir(l5fs,1)),nanmean(orir1(l5fs,1))];
errorbars = [nanstd(orir(l23fs,1))./sqrt(sum(l23fs)),...
    nanstd(orir1(l23fs,1))./sqrt(sum(l23fs));...
    nanstd(orir(l4fs,1))./sqrt(sum(l4fs)),...
    nanstd(orir1(l4fs,1))./sqrt(sum(l4fs));...
    nanstd(orir(l5fs,1))./sqrt(sum(l5fs)),...
    nanstd(orir1(l5fs,1))./sqrt(sum(l5fs))];
figure
barweb(bars, errorbars, [], [{'L23'},{'L4'},{'L5'}], 'tuning ratio FS cells SMALL', [], 'tuning ratio', [], [], [{'L0'},{'L1'}],[],'axis');

bars = [nanmean(orir(l23fs,5)),nanmean(orir1(l23fs,5));...
    nanmean(orir(l4fs,5)),nanmean(orir1(l4fs,5));...
    nanmean(orir(l5fs,5)),nanmean(orir1(l5fs,5))];
errorbars = [nanstd(orir(l23fs,5))./sqrt(sum(l23fs)),...
    nanstd(orir1(l23fs,5))./sqrt(sum(l23fs));...
    nanstd(orir(l4fs,5))./sqrt(sum(l4fs)),...
    nanstd(orir1(l4fs,5))./sqrt(sum(l4fs));...
    nanstd(orir(l5fs,5))./sqrt(sum(l5fs)),...
    nanstd(orir1(l5fs,5))./sqrt(sum(l5fs))];
figure
barweb(bars, errorbars, [], [{'L23'},{'L4'},{'L5'}], 'tuning ratio FS cells LARGE', [], 'tuning ratio', [], [], [{'L0'},{'L1'}],[],'axis');

bars = [nanmean(orir(l23rs,1)),nanmean(orir1(l23rs,1));...
    nanmean(orir(l4rs,1)),nanmean(orir1(l4rs,1));...
    nanmean(orir(l5rs,1)),nanmean(orir1(l5rs,1))];
errorbars = [nanstd(orir(l23rs,1))./sqrt(sum(l23rs)),...
    nanstd(orir1(l23rs,1))./sqrt(sum(l23rs));...
    nanstd(orir(l4rs,1))./sqrt(sum(l4rs)),...
    nanstd(orir1(l4rs,1))./sqrt(sum(l4rs));...
    nanstd(orir(l5rs,1))./sqrt(sum(l5rs)),...
    nanstd(orir1(l5rs,1))./sqrt(sum(l5rs))];
figure
barweb(bars, errorbars, [], [{'L23'},{'L4'},{'L5'}], 'tuning ratio RS cells SMALL', [], 'tuning ratio', [], [], [{'L0'},{'L1'}],[],'axis');

bars = [nanmean(orir(l23rs,5)),nanmean(orir1(l23rs,5));...
    nanmean(orir(l4rs,5)),nanmean(orir1(l4rs,5));...
    nanmean(orir(l5rs,5)),nanmean(orir1(l5rs,5))];
errorbars = [nanstd(orir(l23rs,5))./sqrt(sum(l23rs)),...
    nanstd(orir1(l23rs,5))./sqrt(sum(l23rs));...
    nanstd(orir(l4rs,5))./sqrt(sum(l4rs)),...
    nanstd(orir1(l4rs,5))./sqrt(sum(l4rs));...
    nanstd(orir(l5rs,5))./sqrt(sum(l5rs)),...
    nanstd(orir1(l5rs,5))./sqrt(sum(l5rs))];
figure
barweb(bars, errorbars, [], [{'L23'},{'L4'},{'L5'}], 'tuning ratio RS cells LARGE', [], 'tuning ratio', [], [], [{'L0'},{'L1'}],[],'axis');



mffl0 = nanmean(ffl0,2);
mffl1 = nanmean(ffl1,2);

bars = [nanmean(mffl0(l23rs)),nanmean(mffl1(l23rs));...
    nanmean(mffl0(l4rs)),nanmean(mffl1(l4rs));...
    nanmean(mffl0(l5rs)),nanmean(mffl1(l5rs))];
errorbars = [nanstd(mffl0(l23rs))./sqrt(sum(l23rs)),...
    nanstd(mffl1(l23rs))./sqrt(sum(l23rs));...
    nanstd(mffl0(l4rs))./sqrt(sum(l4rs)),...
    nanstd(mffl1(l4rs))./sqrt(sum(l4rs));...
    nanstd(mffl0(l5rs))./sqrt(sum(l5rs)),...
    nanstd(mffl1(l5rs))./sqrt(sum(l5rs))];
figure
barweb(bars, errorbars, [], [{'L23'},{'L4'},{'L5'}], 'mean fano factors RS cells', [], 'fano factor', [], [], [{'L0'},{'L1'}],[],'axis');
    
bars = [nanmean(mffl0(l23fs)),nanmean(mffl1(l23fs));...
    nanmean(mffl0(l4fs)),nanmean(mffl1(l4fs));...
    nanmean(mffl0(l5fs)),nanmean(mffl1(l5fs))];
errorbars = [nanstd(mffl0(l23fs))./sqrt(sum(l23fs)),...
    nanstd(mffl1(l23fs))./sqrt(sum(l23fs));...
    nanstd(mffl0(l4fs))./sqrt(sum(l4fs)),...
    nanstd(mffl1(l4fs))./sqrt(sum(l4fs));...
    nanstd(mffl0(l5fs))./sqrt(sum(l5fs)),...
    nanstd(mffl1(l5fs))./sqrt(sum(l5fs))];
figure
barweb(bars, errorbars, [], [{'L23'},{'L4'},{'L5'}], 'mean fano factors FS cells', [], 'fano factor', [], [], [{'L0'},{'L1'}],[],'axis');

% delta fr per rank plots
cond = l5fs;
figure
errorbar(mean(sortedl1(cond,:)-sortedl0(cond,:)),std(sortedl1(cond,:)-sortedl0(cond,:))./sqrt(sum(cond)),'.')
hold on
line([0,41],[0,0])
xlabel('rank')
ylabel('delta spikes')

%fr vs fr plots
cond = l23rs;
figure
plot(nanmean(nsortedl0(cond,:)),nanmean(nsortedl1(cond,:)),'ko','markersize',5,'markerfacecolor','k');
axis([0,1.3,0,1.3])
refline(1,0)
lsline
hold on
plot(nanmean(ncontrolfr(cond,1)),nanmean(ncontrolfr(cond,2)),'co','markersize',5,'markerfacecolor','c')
errorbar(nanmean(nsortedl0(cond,:)),nanmean(nsortedl1(cond,:)),nanstd(nsortedl1(cond,:))./sqrt(sum(cond)),'k.')
herrorbar(nanmean(nsortedl0(cond,:)),nanmean(nsortedl1(cond,:)),nanstd(nsortedl0(cond,:))./sqrt(sum(cond)),'k.')
errorbar(nanmean(ncontrolfr(cond,1)),nanmean(ncontrolfr(cond,2)),nanstd(ncontrolfr(cond,2))./sqrt(sum(cond)),'c.')
herrorbar(nanmean(ncontrolfr(cond,1)),nanmean(ncontrolfr(cond,2)),nanstd(ncontrolfr(cond,1))./sqrt(sum(cond)),'c.')
axis([0,1.3,0,1.3])
axis square
xlabel('normalized firing rate light OFF')
ylabel('normalized firing rate light ON')
title(['Layer 4 RS cells n = ' int2str(sum(cond))]);

%fr vs fr plots
cond1 = l4fs&nlfr<lfr;
cond2 = l4fs&nlfr>lfr;
figure
plot(nanmean(nsortedl0(cond1,:)),nanmean(nsortedl1(cond1,:)),'go','markersize',3,'markerfacecolor','g');
hold on
plot(nanmean(nsortedl0(cond2,:)),nanmean(nsortedl1(cond2,:)),'bo','markersize',3,'markerfacecolor','b');
axis([0,1.3,0,1.3])
refline(1,0)
lsline
plot(nanmean(ncontrolfr(cond1,1)),nanmean(ncontrolfr(cond1,2)),'co','markersize',5,'markerfacecolor','c')
errorbar(nanmean(nsortedl0(cond1,:)),nanmean(nsortedl1(cond1,:)),nanstd(nsortedl1(cond1,:))./sqrt(sum(cond1)),'g.')
herrorbar(nanmean(nsortedl0(cond1,:)),nanmean(nsortedl1(cond1,:)),nanstd(nsortedl0(cond1,:))./sqrt(sum(cond1)),'g.')
errorbar(nanmean(ncontrolfr(cond1,1)),nanmean(ncontrolfr(cond1,2)),nanstd(ncontrolfr(cond1,2))./sqrt(sum(cond1)),'c.')
herrorbar(nanmean(ncontrolfr(cond1,1)),nanmean(ncontrolfr(cond1,2)),nanstd(ncontrolfr(cond1,1))./sqrt(sum(cond1)),'c.')
plot(nanmean(ncontrolfr(cond2,1)),nanmean(ncontrolfr(cond2,2)),'yo','markersize',5,'markerfacecolor','y')
errorbar(nanmean(nsortedl0(cond2,:)),nanmean(nsortedl1(cond2,:)),nanstd(nsortedl1(cond2,:))./sqrt(sum(cond2)),'b.')
herrorbar(nanmean(nsortedl0(cond2,:)),nanmean(nsortedl1(cond2,:)),nanstd(nsortedl0(cond2,:))./sqrt(sum(cond2)),'b.')
errorbar(nanmean(ncontrolfr(cond2,1)),nanmean(ncontrolfr(cond2,2)),nanstd(ncontrolfr(cond2,2))./sqrt(sum(cond2)),'y.')
herrorbar(nanmean(ncontrolfr(cond2,1)),nanmean(ncontrolfr(cond2,2)),nanstd(ncontrolfr(cond2,1))./sqrt(sum(cond2)),'y.')
axis([0,1.3,0,1.3])
axis square
xlabel('normalized firing rate light OFF')
ylabel('normalized firing rate light ON')
title(['Layer 4 RS cells n = ' int2str(sum(cond1)) ' m = ' int2str(sum(cond2))]);

% xsizes = [0     8    13    21    36    50];
% oris = [0,45,90,135,180,225,270,315];
% controlerr = ones(size(controlfr,1),2);

% for cell = find(l5&prsv'&~phe')
%     figure
%     subplot(2,2,1)
%     errorbar(oris,squeeze(condfr(cell,2,:,ps(cell))),squeeze(conderr(cell,2,:,ps(cell))),'o-','color',lcol,'markersize',8,'linewidth',2)
%     hold on
%     errorbar(oris,squeeze(condfr(cell,1,:,ps(cell))),squeeze(conderr(cell,1,:,ps(cell))),'ko-','markersize',8,'linewidth',2)
%     xlabel('shown orientation')
%     ylabel('Firing rate [Hz]')
%     set(gca,'xtick',oris)
%     legend({'Light ON','Light OFF'})
%     ax = axis;
%     axis([-10,320,ax(3),ax(4)])
%     title(['orientation tuning, preferred size only cell: ' int2str(cell)])
%     
%     subplot(2,2,2)
%     errorbar(oris,squeeze(mean(condfr(cell,2,:,:),4)),squeeze(mean(conderr(cell,2,:,:),4)),'o-','color',lcol,'markersize',8,'linewidth',2)
%     hold on
%     errorbar(oris,squeeze(mean(condfr(cell,1,:,:),4)),squeeze(mean(conderr(cell,1,:,:),4)),'ko-','markersize',8,'linewidth',2)
%     xlabel('shown orientation')
%     ylabel('Firing rate [Hz]')
%     set(gca,'xtick',oris)
%     legend({'Light ON','Light OFF'})
%     ax = axis;
%     axis([-10,320,ax(3),ax(4)])
%     title(['orientation tuning, mean over sizes, depth: ' int2str(depth(cell))])
%     
%     subplot(2,2,3)
%     %    errorbar(xsizes,sizetunel1(cell,:),sizetuneerrl1(cell,:),'o-','color',lcol,'markersize',8,'linewidth',2);
%     %    hold on
%     %    errorbar(xsizes,sizetunel0(cell,:),sizetuneerrl0(cell,:),'ko-','markersize',8,'linewidth',2);
%     errorbar(xsizes,[controlfr(cell,2);squeeze(condfr(cell,2,prefori(cell),:))],[controlerr(cell,2);squeeze(conderr(cell,2,prefori(cell),:))],'o-','color',lcol,'markersize',8,'linewidth',2);
%     hold on
%     errorbar(xsizes,[controlfr(cell,1);squeeze(condfr(cell,1,prefori(cell),:))],[controlerr(cell,1);squeeze(conderr(cell,1,prefori(cell),:))],'ko-','markersize',8,'linewidth',2);
%     xlabel('shown patch size [vd]')
%     ylabel('Firing rate [Hz]')
%     legend({'Light ON','Light OFF'})
%     ax = axis;
%     axis([-5,50,ax(3),ax(4)])
%     set(gca,'xtick',xsizes)
%     title(['size tuning only preferred orientation  ' cellname{cell}])
%     
%     subplot(2,2,4)
%     errorbar(xsizes,[controlfr(cell,2);squeeze(mean(condfr(cell,2,:,:),3))],[controlerr(cell,2);squeeze(mean(conderr(cell,2,:,:),3))],'o-','color',lcol,'markersize',8,'linewidth',2);
%     hold on
%     errorbar(xsizes,[controlfr(cell,1);squeeze(mean(condfr(cell,1,:,:),3))],[controlerr(cell,1);squeeze(mean(conderr(cell,1,:,:),3))],'ko-','markersize',8,'linewidth',2);
%     xlabel('shown patch size [vd]')
%     ylabel('Firing rate [Hz]')
%     legend({'Light ON','Light OFF'})
%     set(gca,'xtick',xsizes)
%     ax = axis;
%     axis([-5,50,ax(3),ax(4)])
%     title(['size tuning all orientations'])
%     
% %     figSize = [30 21];
% %     set(gcf,'paperunits','centimeters','papersize',figSize,'paperposition',[0 0 figSize])
% %     if cell<10, printi = ['0', int2str(cell)]; else printi = int2str(cell); end
% %     print([l5RSprintpath ,  printi '__' cellname{cell} '.pdf'],'-dpdf')
% 
% end

for i = 1:size(condfr,1)
    l0oricurve = squeeze(condfr(i,1,:,prefs(i)));%-controlfr(i,1);
    l1oricurve = squeeze(condfr(i,2,:,prefs(i)));%-controlfr(i,2);
    
    [l0orir(i), l0dirr(i), x, y, l0osi(i), c, l0dsi(i)] = getOSI(l0oricurve',oris);
    [l1orir(i), l1dirr(i), x, y, l1osi(i), c, l1dsi(i)] = getOSI(l1oricurve',oris);
    
    prefol0(i) = find(l0oricurve == max(l0oricurve),1);
    prefol1(i) = find(l1oricurve == max(l1oricurve),1);
    
    shiftcurvel0(i,:) = circshift(l0oricurve,-(prefol0(i)-1));
    shiftcurvel1(i,:) = circshift(l1oricurve,-(prefol0(i)-1));
    shiftcurvectrll0(i) = controlfr(i,1)./shiftcurvel0(i,1);
    shiftcurvectrll1(i) = controlfr(i,2)./shiftcurvel0(i,1);
%     shiftcurvel1(i,:) = shiftcurvel1(i,:)./shiftcurvel0(i,1);
%     shiftcurvel0(i,:) = shiftcurvel0(i,:)./shiftcurvel0(i,1);
    scl0center(i,:) = circshift(shiftcurvel0(i,:),3,2);
    scl1center(i,:) = circshift(shiftcurvel1(i,:),3,2);
    
    nscl0center(i,:) = scl0center(i,:)./max(scl0center(i,:));
    nscl1center(i,:) = scl1center(i,:)./max(scl0center(i,:));
    
    [s,si] = sort(shiftcurvel0(i,:));
    sortshiftl0(i,:) = shiftcurvel0(i,si);
    sortshiftl1(i,:) = shiftcurvel1(i,si);
    
    %worst size   
    nporicurvel0 = squeeze(condfr(i,1,:,aps(i)));
    nporicurvel1 = squeeze(condfr(i,2,:,aps(i)));
    [s,si] = sort(nporicurvel0);
    nprankcurvel0(i,:) = nporicurvel0(si);
    nprankcurvel1(i,:) = nporicurvel1(si);    
    
    oricondfr(i,:,1,:) = nanmean(condfr(i,:,[1,5],:),3);
    oricondfr(i,:,2,:) = nanmean(condfr(i,:,[2,6],:),3);
    oricondfr(i,:,3,:) = nanmean(condfr(i,:,[3,7],:),3);
    oricondfr(i,:,4,:) = nanmean(condfr(i,:,[4,8],:),3);
    l0s5m = max(oricondfr(i,1,:,5)); if l0s5m == 0, l0s5m = 0.1; end
    l0s1(i,:) = oricondfr(i,1,:,1)./l0s5m;
    l0s5(i,:) = oricondfr(i,1,:,5)./l0s5m;
    l1s1(i,:) = oricondfr(i,2,:,1)./l0s5m;
    l1s5(i,:) = oricondfr(i,2,:,5)./l0s5m;
    prefol0s5(i) = find(l0s5(i,:) == max(l0s5(i,:)),1);
    l0s1shift(i,:) = circshift(l0s1(i,:)',-(prefol0s5(i)-1));
    l0s5shift(i,:) = circshift(l0s5(i,:)',-(prefol0s5(i)-1));
    l1s1shift(i,:) = circshift(l1s1(i,:)',-(prefol0s5(i)-1));
    l1s5shift(i,:) = circshift(l1s5(i,:)',-(prefol0s5(i)-1));
    
end
shiftoris = [-135,-90,-45,0,45,90,135,180];
% cond = l5rs;
% anovavec = [scafrl0(cond,1);scafrl0(cond,2);scafrl0(cond,3);scafrl0(cond,4);scafrl0(cond,5);...
%     scafrl1(cond,1);scafrl1(cond,2);scafrl1(cond,3);scafrl1(cond,4);scafrl1(cond,5)];
% gl = [ones(sum(cond)*5,1);ones(sum(cond)*5,1).*2];
% h = [ones(sum(cond),1)];
% gp = [h.*1;h.*2;h.*3;h.*4;h.*5;h.*1;h.*2;h.*3;h.*4;h.*5];
% [p,table,stats] = anovan(anovavec,{gl,gp},'model','full');
% multcompare(stats)
% multcompare(stats,'dimension',2)
cond = l5rs;
figure
errorbar(oris(1:4),mean(l0s1shift(cond,:)),std(l0s1shift(cond,:))./sqrt(sum(cond)),'b','linewidth',2)
hold on
errorbar(oris(1:4),mean(l0s5shift(cond,:)),std(l0s5shift(cond,:))./sqrt(sum(cond)),'c','linewidth',2)
errorbar(oris(1:4),mean(l1s1shift(cond,:)),std(l1s1shift(cond,:))./sqrt(sum(cond)),'r','linewidth',2)
errorbar(oris(1:4),mean(l1s5shift(cond,:)),std(l1s5shift(cond,:))./sqrt(sum(cond)),'m','linewidth',2)
legend([{'L0 small'},{'L0 large'},{'L1 small'},{'L1 large'}])
xlabel('orientation')
ylabel('normalized firing rate')


cond = l4&prsv'&~phe'&ok;

% figure
% boundedline(respta,mean(l0meanrespn(cond,:)),std(l0meanrespn(cond,:))./sqrt(sum(cond)),'b')
% hold on
% boundedline(respta,mean(l1meanrespn(cond,:)),std(l1meanrespn(cond,:))./sqrt(sum(cond)),'r')
% line([0,0],[.5,2],'color','k')
% line([2000,2000],[.5,2],'color','k')
% line([500,500],[.5,2],'color','k','linestyle',':')
% line([1500,1500],[.5,2],'color','k','linestyle',':')

%ranked oris
cond = l23rs;
errorbar(nanmean(sortshiftl1(cond,:)-sortshiftl0(cond,:)),nanstd(sortshiftl1(cond,:)-sortshiftl0(cond,:))./sqrt(sum(cond)),'.')
hold on
errorbar(nanmean(nprankcurvel1(cond,:)-nprankcurvel0(cond,:)),nanstd(nprankcurvel1(cond,:)-nprankcurvel0(cond,:))./sqrt(sum(cond)),'.')
xlabel('rank according to orientation')
ylabel('average number of spikes added')
legend([{'best size'},{'worst size'}],'location','ne')
line([0,9],[0,0])
title('L2/3 RS')

%ranked psths
cond = l5rs;
rank = 5;
figure
boundedline(respta,squeeze(nanmean(nscarespl0(cond,rank,:))),squeeze(nanstd(nscarespl0(cond,rank,:)))./sqrt(sum(cond)),'k')
hold on
boundedline(respta,squeeze(nanmean(nscarespl1(cond,rank,:))),squeeze(nanstd(nscarespl1(cond,rank,:)))./sqrt(sum(cond)),'r')
line([0,0],[.5,2],'color','k')
line([2000,2000],[.5,2],'color','k')
line([500,500],[.5,2],'color','k','linestyle',':')
line([1500,1500],[.5,2],'color','k','linestyle',':')

figure
for i = 1:5
    subplot(3,5,i)
    cond = l23rs&nlfr>lfr;
    rank = i;
    boundedline(respta,squeeze(nanmean(ntomscarespl0(cond,rank,:))),squeeze(nanstd(ntomscarespl0(cond,rank,:)))./sqrt(sum(cond)),'k')
    hold on
    boundedline(respta,squeeze(nanmean(ntomscarespl1(cond,rank,:))),squeeze(nanstd(ntomscarespl1(cond,rank,:)))./sqrt(sum(cond)),'r')
    line([0,0],[0,2.5],'color','k')
    line([2000,2000],[0,2.5],'color','k')
    line([500,500],[0,2.5],'color','k','linestyle',':')
    line([1500,1500],[0,2.5],'color','k','linestyle',':')
    axis([-300,2300,0,2.5])
    if i == 1, ylabel('normalized FR [Hz]'); end
    
    subplot(3,5,5+i)
    cond = l4rs&nlfr>lfr;
    rank = i;
    boundedline(respta,squeeze(nanmean(ntomscarespl0(cond,rank,:))),squeeze(nanstd(ntomscarespl0(cond,rank,:)))./sqrt(sum(cond)),'k')
    hold on
    boundedline(respta,squeeze(nanmean(ntomscarespl1(cond,rank,:))),squeeze(nanstd(ntomscarespl1(cond,rank,:)))./sqrt(sum(cond)),'r')
    line([0,0],[0,2.5],'color','k')
    line([2000,2000],[0,2.5],'color','k')
    line([500,500],[0,2.5],'color','k','linestyle',':')
    line([1500,1500],[0,2.5],'color','k','linestyle',':')
    axis([-300,2300,0,2.5])
    if i == 1, ylabel('normalized FR [Hz]'); end
    
    subplot(3,5,2*5+i)
    cond = l5rs&nlfr>lfr;
    rank = i;
    boundedline(respta,squeeze(nanmean(ntomscarespl0(cond,rank,:))),squeeze(nanstd(ntomscarespl0(cond,rank,:)))./sqrt(sum(cond)),'k')
    hold on
    boundedline(respta,squeeze(nanmean(ntomscarespl1(cond,rank,:))),squeeze(nanstd(ntomscarespl1(cond,rank,:)))./sqrt(sum(cond)),'r')
    line([0,0],[0,2.5],'color','k')
    line([2000,2000],[0,2.5],'color','k')
    line([500,500],[0,2.5],'color','k','linestyle',':')
    line([1500,1500],[0,2.5],'color','k','linestyle',':')
    axis([-300,2300,0,2.5])
    xlabel('time [ms]')
    if i == 1, ylabel('normalized FR [Hz]'); end
end

%ranked across oris
cond = l23rs;
figure
errorbar(mean(sortshiftl0(cond,:)),std(sortshiftl0(cond,:))./sqrt(sum(cond)),'k','linewidth',2)
hold on
errorbar(mean(sortshiftl1(cond,:)),std(sortshiftl0(cond,:))./sqrt(sum(cond)),'r','linewidth',2)
axis([0,9,0,1.2])
xlabel('rank')
ylabel('normalized firing rate')



cond = l4rs;
figure
errorbar(oris,mean(shiftcurvel0(cond,:)),std(shiftcurvel0(cond,:))./sqrt(length(find(cond))),'k','linewidth',2)
hold on
errorbar(oris,mean(shiftcurvel1(cond,:)),std(shiftcurvel1(cond,:))./sqrt(length(find(cond))),'r','linewidth',2)
xlabel('degrees from optimal')
ylabel('normalized FR')
axis([-10,325,-.2,1.2])

oriscenter = circshift(oris,3,2);
oriscenter(3) = -45; oriscenter(2) = -90; oriscenter(1) = -135;

cond = l5fs & nlfr<lfr;%&nlosi>.2;
figure
errorbar(oriscenter,mean(scl0center(cond,:)),std(scl0center(cond,:))./sqrt(length(find(cond))),'k','linewidth',2)
hold on
errorbar(oriscenter,mean(scl1center(cond,:)),std(scl1center(cond,:))./sqrt(length(find(cond))),'r','linewidth',2)
errorbar(225, mean(shiftcurvectrll0(cond)),std(shiftcurvectrll0(cond))./sqrt(sum(cond)),'k','linewidth',2)
errorbar(225, mean(shiftcurvectrll1(cond)),std(shiftcurvectrll1(cond))./sqrt(sum(cond)),'r','linewidth',2)
xlabel('degrees from optimal')
ylabel('normalized FR')
title(['L4 RS cells n = ' int2str(sum(cond))]);
set(gca,'xtick',[-135,-90,-45,0,45,90,135,180])
axis([-150,240,0,1.3])

figure
cond = l23rs;
errorbar(mean(sncafrl0(cond,:)),mean(sncafrl1(cond,:)),std(sncafrl1(cond,:))./sqrt(sum(cond)),'.-')
hold on
herrorbar(mean(sncafrl0(cond,:)),mean(sncafrl1(cond,:)),std(sncafrl0(cond,:))./sqrt(sum(cond)),'.-')
line([0,1.3],[0,1.3])
axis([0,1.3,0,1.3])
xlabel('normalized ranked firing rate light OFF')
ylabel('normalized ranked firing rate light ON')
title('size ranked')
axis square

cond = l5rs&ps==3;
figure
errorbar(xsizes(2:6),squeeze(mean(nprefsizetunes(cond,1,2:6))),...
    squeeze(std(nprefsizetunes(cond,1,2:6)))./sqrt(sum(cond)),'k','linewidth',2);
hold on
errorbar(xsizes(2:6),squeeze(mean(nprefsizetunes(cond,2,2:6))),...
    squeeze(std(nprefsizetunes(cond,2,2:6)))./sqrt(sum(cond)),'r','linewidth',2);
errorbar(xsizes(1),squeeze(mean(nprefsizetunes(cond,1,1))),...
    squeeze(std(nprefsizetunes(cond,1,1)))./sqrt(sum(cond)),'k','linewidth',2);
errorbar(xsizes(1),squeeze(mean(nprefsizetunes(cond,2,1))),...
    squeeze(std(nprefsizetunes(cond,2,1)))./sqrt(sum(cond)),'r','linewidth',2);
xlabel('shown size')
ylabel('normalized firing rate')
axis([-5,50,0,1.3])

cond = l5rs&(ps==1);
figure
errorbar(xsizes(2:6),nanmean(ncafrl0(cond,:)),nanstd(ncafrl0(cond,:))./sqrt(sum(cond)),'k','linewidth',2)
hold on
errorbar(xsizes(2:6),nanmean(ncafrl1(cond,:)),nanstd(ncafrl1(cond,:))./sqrt(sum(cond)),'r','linewidth',2)
xlabel('size [dva]')
ylabel('normalized firing rate')

cond = l5a&prsv'&~phe'&(ps==3);
figure
errorbar(squeeze(mean(normprefsizetunes(cond,1,:))),...
    squeeze(std(normprefsizetunes(cond,1,:)))./sqrt(length(find(cond))),'k','linewidth',2);
hold on
errorbar(squeeze(mean(normprefsizetunes(cond,2,:))),...
    squeeze(std(normprefsizetunes(cond,2,:)))./sqrt(length(find(cond))),'r','linewidth',2);

cond = l5&prsv'&~phe';
a = find(cond);
figure
hold on
for i = a
    plot([min(condspikesl0(i)),max(condspikesl1(i))],[min(condspikesl0(i))*params(i,1)+params(i,2),max(condspikesl1(i))*params(i,1)+params(i,2)])
end

figure
hold on
for i = a
    plot([min(normcondspikesl0(i)),max(normcondspikesl1(i))],[min(normcondspikesl0(i))*normparams(i,1)+normparams(i,2),max(normcondspikesl1(i))*normparams(i,1)+normparams(i,2)])
end

%PSTHs
cond = l5&pfsv'&~phe'&ok;
figure
boundedline(respta,mean(l0meanrespn(cond,:)),std(l0meanrespn(cond,:))./sqrt(sum(cond)),'b')
hold on
boundedline(respta,mean(l1meanrespn(cond,:)),std(l1meanrespn(cond,:))./sqrt(sum(cond)),'r')
xlabel('time [ms]')
ylabel('firing rate [Hz]')
title('mean response')

figure
boundedline(respta,mean(l0prefrespn(cond,:)),std(l0prefrespn(cond,:))./sqrt(sum(cond)),'b')
hold on
boundedline(respta,mean(l1prefrespn(cond,:)),std(l1prefrespn(cond,:))./sqrt(sum(cond)),'r')
xlabel('time [ms]')
ylabel('firing rate [Hz]')
title('preferred response')

cond = l5&prsv'&~phe'&ok;
rank = 1;
figure
boundedline(respta,squeeze(nanmean(nscarespl0(cond,rank,:))),squeeze(nanstd(nscarespl0(cond,rank,:)))./sqrt(sum(cond)),'k')
hold on
boundedline(respta,squeeze(nanmean(nscarespl1(cond,rank,:))),squeeze(nanstd(nscarespl1(cond,rank,:)))./sqrt(sum(cond)),'r')




% Scott's Paper
figure
rsomi = (lfr(prsv'&~phe'&ok)-nlfr(prsv'&~phe'&ok))./(lfr(prsv'&~phe'&ok)+nlfr(prsv'&~phe'&ok));
plot(rsomi,depth(prsv'&~phe'&ok),'ko','markersize',4,'markerfacecolor','k')
line([0,0],[0,1000],'color','k')
line([-1,1,],[375,375],'color','k','linestyle',':')
line([-1,1,],[500,500],'color','k','linestyle',':')
hold on
[x,y,xerr] = runningMedian(depth(prsv'&~phe'&ok),rsomi,0,12);
plot((lfr(prsv&phe&ok')-nlfr(prsv&phe&ok'))./(lfr(prsv&phe&ok')+nlfr(prsv&phe&ok')),depth(prsv&phe&ok'),'mo','markersize',4,'markerfacecolor','m')
plot(x,y,'k','linewidth',2);
plot(x+xerr,y,'k')
plot(x-xerr,y,'k')
axis ij
axis([-1.1,1.1,0,1000])
ylabel('depth[mum]')
xlabel('(light-nolight)/(light+nolight)')
title('RS cells: average firing rate changes by layer 4 suppression')

figure
fsomi = (lfr(pfsv'&~phe'&ok)-nlfr(pfsv'&~phe'&ok))./(lfr(pfsv'&~phe'&ok)+nlfr(pfsv'&~phe'&ok));
plot(fsomi,depth(pfsv'&~phe'&ok),'go','markersize',4,'markerfacecolor','g')
line([0,0],[0,1000],'color','k')
line([-1,1,],[375,375],'color','k','linestyle',':')
line([-1,1,],[500,500],'color','k','linestyle',':')
hold on
plot((lfr(pfsv&phe&ok')-nlfr(pfsv&phe&ok'))./(lfr(pfsv&phe&ok')+nlfr(pfsv&phe&ok')),depth(pfsv&phe&ok'),'mo','markersize',4,'markerfacecolor','m')
[x,y,xerr] = runningMedian(depth(pfsv'&~phe'&ok),fsomi,0,12);
plot(x,y,'g','linewidth',2);
plot(x+xerr,y,'g')
plot(x-xerr,y,'g')
axis ij
axis([-1.1,1.1,0,1000])
ylabel('depth[mum]')
xlabel('(light-nolight)/(light+nolight)')
title('FS cells: average firing rate changes by layer 4 suppression')

% running
figure
plot(r0omi(prsv),depth(prsv),'o','markersize',4,'markerfacecolor','b','color','b')
line([0,0],[0,1000],'color','k')
hold on
[x,y,xerr] = runningMedian(depth(prsv&~phe),r0omi(prsv&~phe),0,12);
plot(x,y,'linewidth',2,'color','b');
plot(x+xerr,y,'color','b')
plot(x-xerr,y,'color','b')
line([-1,1,],[375,375],'color','k','linestyle',':')
line([-1,1,],[500,500],'color','k','linestyle',':')
plot(r0omi(prsv&phe),depth(prsv&phe),'mo','markersize',4,'markerfacecolor','m')
axis ij
axis([-1.1,1.1,150,950])
ylabel('depth[mum]')
xlabel('OMI of average visual response rate')
title('Quiet: RS cells: average firing rate changes by layer 4 suppression')
set(gca,'box','off')

figure
plot(r1omi(prsv),depth(prsv),'o','markersize',4,'markerfacecolor','b','color','b')
line([0,0],[0,1000],'color','k')
hold on
[x,y,xerr] = runningMedian(depth(prsv&~phe),r1omi(prsv&~phe),0,12);
plot(x,y,'linewidth',2,'color','b');
plot(x+xerr,y,'color','b')
plot(x-xerr,y,'color','b')
line([-1,1,],[375,375],'color','k','linestyle',':')
line([-1,1,],[500,500],'color','k','linestyle',':')
plot(r1omi(prsv&phe),depth(prsv&phe),'mo','markersize',4,'markerfacecolor','m')
axis ij
axis([-1.1,1.1,150,950])
ylabel('depth[mum]')
xlabel('OMI of average visual response rate')
title('Running: RS cells: average firing rate changes by layer 4 suppression')
set(gca,'box','off')

figure
plot(r0omi(pfs),depth(pfs),'o','markersize',4,'markerfacecolor','r','color','r')
line([0,0],[0,1000],'color','k')
hold on
[x,y,xerr] = runningMedian(depth(pfs),r0omi(pfs),0,12);
plot(x,y,'linewidth',2,'color','r');
plot(x+xerr,y,'color','r')
plot(x-xerr,y,'color','r')
line([-1,1,],[375,375],'color','k','linestyle',':')
line([-1,1,],[500,500],'color','k','linestyle',':')
axis ij
axis([-1.1,1.1,150,950])
ylabel('depth[mum]')
xlabel('OMI of average visual response rate')
title('Quiet: FS cells: average firing rate changes by layer 4 suppression')
set(gca,'box','off')

figure
plot(r1omi(pfs),depth(pfs),'o','markersize',4,'markerfacecolor','r','color','r')
line([0,0],[0,1000],'color','k')
hold on
[x,y,xerr] = runningMedian(depth(pfs),r1omi(pfs),0,12);
plot(x,y,'linewidth',2,'color','r');
plot(x+xerr,y,'color','r')
plot(x-xerr,y,'color','r')
line([-1,1,],[375,375],'color','k','linestyle',':')
line([-1,1,],[500,500],'color','k','linestyle',':')
axis ij
axis([-1.1,1.1,150,950])
ylabel('depth[mum]')
xlabel('OMI of average visual response rate')
title('Running: FS cells: average firing rate changes by layer 4 suppression')
set(gca,'box','off')

% RMI

figure
plot(l0rmi(prsv),depth(prsv),'o','markersize',4,'markerfacecolor','b','color','b')
line([0,0],[0,1000],'color','k')
hold on
[x,y,xerr] = runningMedian(depth(prsv&~phe),l0rmi(prsv&~phe),0,12);
plot(x,y,'linewidth',2,'color','b');
plot(x+xerr,y,'color','b')
plot(x-xerr,y,'color','b')
line([-1,1,],[375,375],'color','k','linestyle',':')
line([-1,1,],[500,500],'color','k','linestyle',':')
plot(l0rmi(prsv&phe),depth(prsv&phe),'mo','markersize',4,'markerfacecolor','m')
axis ij
axis([-1.1,1.1,150,950])
ylabel('depth[mum]')
xlabel('RMI of average visual response rate')
title('NO light: RS cells: average firing rate changes by layer 4 suppression')
set(gca,'box','off')

figure
plot(l1rmi(prsv),depth(prsv),'o','markersize',4,'markerfacecolor','b','color','b')
line([0,0],[0,1000],'color','k')
hold on
[x,y,xerr] = runningMedian(depth(prsv&~phe),l1rmi(prsv&~phe),0,12);
plot(x,y,'linewidth',2,'color','b');
plot(x+xerr,y,'color','b')
plot(x-xerr,y,'color','b')
line([-1,1,],[375,375],'color','k','linestyle',':')
line([-1,1,],[500,500],'color','k','linestyle',':')
plot(l1rmi(prsv&phe),depth(prsv&phe),'mo','markersize',4,'markerfacecolor','m')
axis ij
axis([-1.1,1.1,150,950])
ylabel('depth[mum]')
xlabel('RMI of average visual response rate')
title('Light ON: RS cells: average firing rate changes by layer 4 suppression')
set(gca,'box','off')

figure
plot(l0rmi(pfs),depth(pfs),'o','markersize',4,'markerfacecolor','r','color','r')
line([0,0],[0,1000],'color','k')
hold on
[x,y,xerr] = runningMedian(depth(pfs),l0rmi(pfs),0,12);
plot(x,y,'linewidth',2,'color','r');
plot(x+xerr,y,'color','r')
plot(x-xerr,y,'color','r')
line([-1,1,],[375,375],'color','k','linestyle',':')
line([-1,1,],[500,500],'color','k','linestyle',':')
axis ij
axis([-1.1,1.1,150,950])
ylabel('depth[mum]')
xlabel('RMI of average visual response rate')
title('NO Light: FS cells: average firing rate changes by layer 4 suppression')
set(gca,'box','off')

figure
plot(l1rmi(pfs),depth(pfs),'o','markersize',4,'markerfacecolor','r','color','r')
line([0,0],[0,1000],'color','k')
hold on
[x,y,xerr] = runningMedian(depth(pfs),l1rmi(pfs),0,12);
plot(x,y,'linewidth',2,'color','r');
plot(x+xerr,y,'color','r')
plot(x-xerr,y,'color','r')
line([-1,1,],[375,375],'color','k','linestyle',':')
line([-1,1,],[500,500],'color','k','linestyle',':')
axis ij
axis([-1.1,1.1,150,950])
ylabel('depth[mum]')
xlabel('RMI of average visual response rate')
title('Light ON: FS cells: average firing rate changes by layer 4 suppression')
set(gca,'box','off')

% population size tuning
figure
errorbar(mean(sizetunel0(l23rs,:)),mean(sizetuneerrl0(l23rs,:)),'b.-')
hold on
errorbar(mean(sizetunel0(l23fs,:)),mean(sizetuneerrl0(l23fs,:)),'r.-')
title('L23 light OFF')
legend([{'RS'},{'FS'}])
axis([0,7,2,14])

figure
errorbar(mean(sizetunel0(l4rs,:)),mean(sizetuneerrl0(l4rs,:)),'b.-')
hold on
errorbar(mean(sizetunel0(l4fs,:)),mean(sizetuneerrl0(l4fs,:)),'r.-')
title('L4 light OFF')
legend([{'RS'},{'FS'}])
axis([0,7,2,14])

figure
errorbar(mean(sizetunel0(l5rs,:)),mean(sizetuneerrl0(l5rs,:)),'b.-')
hold on
errorbar(mean(sizetunel0(l5fs,:)),mean(sizetuneerrl0(l5fs,:)),'r.-')
title('L56 light OFF')
legend([{'RS'},{'FS'}])
axis([0,7,2,14])

figure
errorbar(mean(sizetunel1(l23rs,:)),mean(sizetuneerrl1(l23rs,:)),'b.-')
hold on
errorbar(mean(sizetunel1(l23fs,:)),mean(sizetuneerrl1(l23fs,:)),'r.-')
title('L23 light ON')
legend([{'RS'},{'FS'}])
axis([0,7,2,14])

figure
errorbar(mean(sizetunel1(l4rs,:)),mean(sizetuneerrl1(l4rs,:)),'b.-')
hold on
errorbar(mean(sizetunel1(l4fs,:)),mean(sizetuneerrl1(l4fs,:)),'r.-')
title('L4 light ON')
legend([{'RS'},{'FS'}])
axis([0,7,2,14])

figure
errorbar(mean(sizetunel1(l5rs,:)),mean(sizetuneerrl1(l5rs,:)),'b.-')
hold on
errorbar(mean(sizetunel1(l5fs,:)),mean(sizetuneerrl1(l5fs,:)),'r.-')
title('L5 light ON')
legend([{'RS'},{'FS'}])
axis([0,7,2,14])



% LFP phas stuff
% phase locking
[sd,si] = sort(depth);

% LFP phas stuff
% RS cells
figure
subplot(2,2,1)
imagesc(3:5:98,rsd,allrl0(prsv,:));
caxis([0,.5])
colorbar
title('RS cells: r of spike phases - light OFF')
xlabel('frequency bands')
ylabel('cell number - sorted by depth')

subplot(2,2,2)
imagesc(3:5:98,rsd,allrl1(prsv,:));
caxis([0,.5])
colorbar
title('RS cells: r of spike phases - light ON')
xlabel('frequency bands')
ylabel('cell number - sorted by depth')

subplot(2,2,3)
imagesc(3:5:98,fsd,allrl0(pfsv,:));
caxis([0,.5])
colorbar
title('FS cells: r of spike phases - light OFF')
xlabel('frequency bands')
ylabel('cell number - sorted by depth')

subplot(2,2,4)
imagesc(3:5:98,fsd,allrl1(pfsv,:));
caxis([0,.5])
colorbar
title('FS cells: r of spike phases - light ON')
xlabel('frequency bands')
ylabel('cell number - sorted by depth')


figure
subplot(2,2,1)
errorbar(3:5:98,nanmean(allrl0(l4fs,:)),nanstd(allrl0(l4fs,:))./sqrt(length(find(l4fs))),'r')
hold on
errorbar(3:5:98,nanmean(allrl0(l4rs,:)),nanstd(allrl0(l4rs,:))./sqrt(length(find(l4rs))),'b')
legend([{'FS'},{'RS'}]);
xlabel('frequency [Hz]')
ylabel('average phase of spike');
title('Light OFF: Layer 4');
axis([0,100,0,.3])

subplot(2,2,2)
errorbar(3:5:98,nanmean(allrl0(l23fs,:)),nanstd(allrl0(l23fs,:))./sqrt(length(find(l23fs))),'r')
hold on
errorbar(3:5:98,nanmean(allrl0(l23rs,:)),nanstd(allrl0(l23rs,:))./sqrt(length(find(l23rs))),'b')
legend([{'FS'},{'RS'}]);
xlabel('frequency [Hz]')
ylabel('average phase of spike');
title('Light OFF: Layer 2/3');
axis([0,100,0,.3])

subplot(2,2,3)
errorbar(3:5:98,nanmean(allrl0(l5fs,:)),nanstd(allrl0(l5fs,:))./sqrt(length(find(l5fs))),'r')
hold on
errorbar(3:5:98,nanmean(allrl0(l5rs,:)),nanstd(allrl0(l5rs,:))./sqrt(length(find(l5rs))),'b')
legend([{'FS'},{'RS'}]);
xlabel('frequency [Hz]')
ylabel('average phase of spike');
title('Light OFF: Layer 5');
axis([0,100,0,.3])


figure
subplot(2,2,1)
errorbar(3:5:98,nanmean(allrl1(l4fs,:)),nanstd(allrl1(l4fs,:))./sqrt(length(find(l4fs))),'r')
hold on
errorbar(3:5:98,nanmean(allrl1(l4rs,:)),nanstd(allrl1(l4rs,:))./sqrt(length(find(l4rs))),'b')
legend([{'FS'},{'RS'}]);
xlabel('frequency [Hz]')
ylabel('average phase of spike');
title('Light ON: Layer 4');
axis([0,100,0,.3])

subplot(2,2,2)
errorbar(3:5:98,nanmean(allrl1(l23fs,:)),nanstd(allrl1(l23fs,:))./sqrt(length(find(l23fs))),'r')
hold on
errorbar(3:5:98,nanmean(allrl1(l23rs,:)),nanstd(allrl1(l23rs,:))./sqrt(length(find(l23rs))),'b')
legend([{'FS'},{'RS'}]);
xlabel('frequency [Hz]')
ylabel('average phase of spike');
title('Light ON: Layer 2/3');
axis([0,100,0,.3])

subplot(2,2,3)
errorbar(3:5:98,nanmean(allrl1(l5fs,:)),nanstd(allrl1(l5fs,:))./sqrt(length(find(l5fs))),'r')
hold on
errorbar(3:5:98,nanmean(allrl1(l5rs,:)),nanstd(allrl1(l5rs,:))./sqrt(length(find(l5rs))),'b')
legend([{'FS'},{'RS'}]);
xlabel('frequency [Hz]')
ylabel('average phase of spike');
title('Light ON: Layer 5');
axis([0,100,0,.3])

%light effect on RS
figure
subplot(2,2,1)
errorbar(3:5:98,nanmean(allrl0(l4rs,:)),nanstd(allrl0(l4rs,:))./sqrt(length(find(l4rs))),'b')
hold on
errorbar(3:5:98,nanmean(allrl1(l4rs,:)),nanstd(allrl1(l4rs,:))./sqrt(length(find(l4rs))),'r')
legend([{'L0'},{'L1'}]);
xlabel('frequency [Hz]')
ylabel('average phase of spike');
title('RS cells: Layer 4');
axis([0,100,0,.3])

subplot(2,2,2)
errorbar(3:5:98,nanmean(allrl0(l23rs,:)),nanstd(allrl0(l23rs,:))./sqrt(length(find(l23rs))),'b')
hold on
errorbar(3:5:98,nanmean(allrl1(l23rs,:)),nanstd(allrl1(l23rs,:))./sqrt(length(find(l23rs))),'r')
legend([{'L0'},{'L1'}]);
xlabel('frequency [Hz]')
ylabel('average phase of spike');
title('RS cells: Layer 2/3');
axis([0,100,0,.3])

subplot(2,2,3)
errorbar(3:5:98,nanmean(allrl0(l5rs,:)),nanstd(allrl0(l5rs,:))./sqrt(length(find(l5rs))),'b')
hold on
errorbar(3:5:98,nanmean(allrl1(l5rs,:)),nanstd(allrl1(l5rs,:))./sqrt(length(find(l5rs))),'r')
legend([{'L0'},{'L1'}]);
xlabel('frequency [Hz]')
ylabel('average phase of spike');
title('RS cells: Layer 5');
axis([0,100,0,.3])

%light effect on FS
figure
subplot(2,2,1)
errorbar(3:5:98,nanmean(allrl0(l4fs,:)),nanstd(allrl0(l4fs,:))./sqrt(length(find(l4fs))),'b')
hold on
errorbar(3:5:98,nanmean(allrl1(l4fs,:)),nanstd(allrl1(l4fs,:))./sqrt(length(find(l4fs))),'r')
legend([{'L0'},{'L1'}]);
xlabel('frequency [Hz]')
ylabel('average phase of spike');
title('FS cells: Layer 4');
axis([0,100,0,.3])

subplot(2,2,2)
errorbar(3:5:98,nanmean(allrl0(l23fs,:)),nanstd(allrl0(l23fs,:))./sqrt(length(find(l23fs))),'b')
hold on
errorbar(3:5:98,nanmean(allrl1(l23fs,:)),nanstd(allrl1(l23fs,:))./sqrt(length(find(l23fs))),'r')
legend([{'L0'},{'L1'}]);
xlabel('frequency [Hz]')
ylabel('average phase of spike');
title('FS cells: Layer 2/3');
axis([0,100,0,.3])

subplot(2,2,3)
errorbar(3:5:98,nanmean(allrl0(l5fs,:)),nanstd(allrl0(l5fs,:))./sqrt(length(find(l5fs))),'b')
hold on
errorbar(3:5:98,nanmean(allrl1(l5fs,:)),nanstd(allrl1(l5fs,:))./sqrt(length(find(l5fs))),'r')
legend([{'L0'},{'L1'}]);
xlabel('frequency [Hz]')
ylabel('average phase of spike');
title('FS cells: Layer 5');
axis([0,100,0,.3])


figure
subplot(2,2,1)
imagesc(3:5:98,rsd,allrl1(prsv,:)-allrl0(prsv,:));
caxis([-.2,.2])
colorbar
title('RS cells: delta r of spike phases - light ON-OFF')
xlabel('frequency bands')
ylabel('cell number - sorted by depth')

subplot(2,2,2)
imagesc(3:5:98,fsd,allrl1(pfsv,:)-allrl0(pfsv,:));
caxis([-.2,.2])
colorbar
title('FS cells: delta r of spike phases - light ON-OFF')
xlabel('frequency bands')
ylabel('cell number - sorted by depth')

subplot(2,2,3)
errorbar(3:5:98,nanmean(allrl1(prs,:)-allrl0(prs,:)),nanstd(allrl1(prs,:)-allrl0(prs,:))./sqrt(size(allrl0(prs,:),1)),'.')
hold on
errorbar(3:5:98,nanmean(allrl1(pfs,:)-allrl0(pfs,:)),nanstd(allrl1(pfs,:)-allrl0(pfs,:))./sqrt(size(allrl0(pfs,:),1)),'r.')
plot(3:5:98,nanmean(allrl1(prs,:)-allrl0(prs,:)),'linewidth',2)
plot(3:5:98,nanmean(allrl1(pfs,:)-allrl0(pfs,:)),'r','linewidth',2)
axis([-1,100,-.04,.04])
line([-1,100],[0,0],'color','k')


figure
subplot(2,2,1)
imagesc(3:5:98,sd,rad2deg(uncircle(allcmeanl0(prsv,:))));
caxis([0,360])
colormap hsv
colorbar
title('RS: average spike phases - light OFF')
xlabel('frequency bands')
ylabel('cell number - sorted by depth')

subplot(2,2,2)
imagesc(3:5:98,sd,rad2deg(uncircle(allcmeanl1(prsv,:))));
caxis([0,360])
colormap hsv
colorbar
title('RS: average spike phases - light ON')
xlabel('frequency bands')
ylabel('cell number - sorted by depth')

subplot(2,2,3)
imagesc(3:5:98,sd,rad2deg(uncircle(allcmeanl0(pfsv,:))));
caxis([0,2*pi])
colormap hsv
colorbar
title('FS: average spike phases - light OFF')
xlabel('frequency bands')
ylabel('cell number - sorted by depth')

subplot(2,2,4)
imagesc(3:5:98,sd,rad2deg(uncircle(allcmeanl1(pfsv,:))));
caxis([0,2*pi])
colormap hsv
colorbar
title('FS: average spike phases - light ON')
xlabel('frequency bands')
ylabel('cell number - sorted by depth')

%phase advance

figure
subplot(2,2,1)
errorbar(3:5:98,mean(rad2deg(uncircle(allcmeanl0(l4fs,:)))),std(rad2deg(uncircle(allcmeanl0(l4fs,:))))./sqrt(length(l4fs)),'r')
hold on
errorbar(3:5:98,mean(rad2deg(uncircle(allcmeanl0(l4rs,:)))),std(rad2deg(uncircle(allcmeanl0(l4rs,:))))./sqrt(length(l4rs)),'b')
line([0,100],[pi,pi],'color','k')
legend([{'FS'},{'RS'}]);
xlabel('frequency [Hz]')
ylabel('average phase of spike');
title('Light OFF: Layer 4');

subplot(2,2,2)
errorbar(3:5:98,mean(rad2deg(uncircle(allcmeanl0(l23fs,:)))),std(rad2deg(uncircle(allcmeanl0(l23fs,:))))./sqrt(length(find(l23fs))),'r')
hold on
errorbar(3:5:98,mean(rad2deg(uncircle(allcmeanl0(l23rs,:)))),std(rad2deg(uncircle(allcmeanl0(l23rs,:))))./sqrt(length(find(l23rs))),'b')
line([0,100],[pi,pi],'color','k')
legend([{'FS'},{'RS'}]);
xlabel('frequency [Hz]')
ylabel('average phase of spike');
title('Light OFF: Layer 2/3');

subplot(2,2,3)
errorbar(3:5:98,mean(rad2deg(uncircle(allcmeanl0(l5fs,:)))),std(rad2deg(uncircle(allcmeanl0(l5fs,:))))./sqrt(length(find(l5fs))),'r')
hold on
errorbar(3:5:98,mean(rad2deg(uncircle(allcmeanl0(l5rs,:)))),std(rad2deg(uncircle(allcmeanl0(l5rs,:))))./sqrt(length(find(l5rs))),'b')
line([0,100],[pi,pi],'color','k')
legend([{'FS'},{'RS'}]);
xlabel('frequency [Hz]')
ylabel('average phase of spike');
title('Light OFF: Layer 5');

figure

subplot(2,2,1)
errorbar(3:5:98,mean(rad2deg(uncircle(allcmeanl1(l4fs,:)))),std(rad2deg(uncircle(allcmeanl1(l4fs,:))))./sqrt(length(find(l4fs))),'r')
hold on
errorbar(3:5:98,mean(rad2deg(uncircle(allcmeanl1(l4rs,:)))),std(rad2deg(uncircle(allcmeanl1(l4rs,:))))./sqrt(length(find(l4rs))),'b')
line([0,100],[pi,pi],'color','k')
legend([{'FS'},{'RS'}]);
xlabel('frequency [Hz]')
ylabel('average phase of spike');
title('Light ON: Layer 4');

subplot(2,2,2)
errorbar(3:5:98,mean(rad2deg(uncircle(allcmeanl1(l23fs,:)))),std(rad2deg(uncircle(allcmeanl1(l23fs,:))))./sqrt(length(find(l23fs))),'r')
hold on
errorbar(3:5:98,mean(rad2deg(uncircle(allcmeanl1(l23rs,:)))),std(rad2deg(uncircle(allcmeanl1(l23rs,:))))./sqrt(length(find(l23rs))),'b')
line([0,100],[pi,pi],'color','k')
legend([{'FS'},{'RS'}]);
xlabel('frequency [Hz]')
ylabel('average phase of spike');
title('Light ON: Layer 2/3');

subplot(2,2,3)
errorbar(3:5:98,mean(rad2deg(uncircle(allcmeanl1(l5fs,:)))),std(rad2deg(uncircle(allcmeanl1(l5fs,:))))./sqrt(length(find(l5fs))),'r')
hold on
errorbar(3:5:98,mean(rad2deg(uncircle(allcmeanl1(l5rs,:)))),std(rad2deg(uncircle(allcmeanl1(l5rs,:))))./sqrt(length(find(l5rs))),'b')
line([0,100],[pi,pi],'color','k')
legend([{'FS'},{'RS'}]);
xlabel('frequency [Hz]')
ylabel('average phase of spike');
title('Light ON: Layer 5');

figure
subplot(2,2,1)
imagesc(3:5:98,sd,circ_dist(allcmeanl1(pfsv,:),allcmeanl0(pfsv,:)));
caxis([-pi,pi])
colorbar
title('FS cells: delta spike phases - light ON-OFF')
xlabel('frequency bands')
ylabel('cell number - sorted by depth')

subplot(2,2,2)
imagesc(3:5:98,sd,circ_dist(allcmeanl1(prsv,:),allcmeanl0(prsv,:)));
caxis([-pi,pi])
colorbar
title('RS cells: delta spike phases - light ON-OFF')
xlabel('frequency bands')
ylabel('cell number - sorted by depth')

subplot(2,2,3)
errorbar(3:5:98,circ_mean(circ_dist(allcmeanl1(prs,:),allcmeanl0(prs,:))),circ_std(circ_dist(allcmeanl1(prs,:),allcmeanl0(prs,:)))./sqrt(size(allcmeanl0(prs,:),1)),'.')
hold on
errorbar(3:5:98,circ_mean(circ_dist(allcmeanl1(pfs,:),allcmeanl0(pfs,:))),circ_std(circ_dist(allcmeanl1(pfs,:),allcmeanl0(pfs,:)))./sqrt(size(allcmeanl0(pfs,:),1)),'r.')
plot(3:5:98,circ_mean(circ_dist(allcmeanl1(prs,:),allcmeanl0(prs,:))),'linewidth',2)
plot(3:5:98,circ_mean(circ_dist(allcmeanl1(pfs,:),allcmeanl0(pfs,:))),'r','linewidth',2)
axis([-1,100,-.2,.2])
line([-1,100],[0,0],'color','k')


% orientation
figure
plot(nloriprefratio(prs),depth(prs),'b.')
axis ij
axis([-.1,1.1,150,1050])
title('RS: OSI in depth')
xlabel('OSI')
ylabel('cortical depth [mum]')

figure
plot(nloriprefratio(prs),depth(prs),'b.')
hold on
[x,y,err] = runningAverage(depth(prs),nloriprefratio(prs),0,15);
plot(x,y,'k','linewidth',2)
plot(x+err,y,'k')
plot(x-err,y,'k')
axis ij
axis([-.1,1.1,150,1050])
title('FS: OSI in depth')
xlabel('OSI')
ylabel('cortical depth [mum]')

figure
plot(nloriprefratio(pfs),depth(pfs),'r.')
hold on
[x,y,err] = runningAverage(depth(pfs),nloriprefratio(pfs),0,15);
plot(x,y,'k','linewidth',2)
plot(x+err,y,'k')
plot(x-err,y,'k')
axis ij
axis([-.1,1.1,150,1050])
title('FS: OSI in depth')
xlabel('OSI')
ylabel('cortical depth [mum]')

figure
plot(nldirprefratio(prs),depth(prs),'b.')
hold on
[x,y,err] = runningAverage(depth(prs),nldirprefratio(prs),0,15);
plot(x,y,'k','linewidth',2)
plot(x+err,y,'k')
plot(x-err,y,'k')
axis ij
axis([-.1,1.1,150,1050])
title('RS: DSI in depth')
xlabel('DSI')
ylabel('cortical depth [mum]')

figure
plot(nldirprefratio(pfs),depth(pfs),'r.')
hold on
[x,y,err] = runningAverage(depth(pfs),nldirprefratio(pfs),0,15);
plot(x,y,'k','linewidth',2)
plot(x+err,y,'k')
plot(x-err,y,'k')
axis ij
axis([-.1,1.1,150,1050])
title('FS: DSI in depth')
xlabel('DSI')
ylabel('cortical depth [mum]')

figure
plot((loriprefratio(prs)-nloriprefratio(prs))./(loriprefratio(prs)+nloriprefratio(prs)),depth(prs),'b.')
axis ij
axis([-1.1,1.1,150,1050])
line([0,0],[150,1050],'color','k')
title('RS: OSI changes in depth')
xlabel('(OSI light - OSI no light)/(OSI light + OSI no light)')
ylabel('cortical depth [mum]')

figure
plot((loriprefratio(pfs)-nloriprefratio(pfs))./(loriprefratio(pfs)+nloriprefratio(pfs)),depth(pfs),'r.')
axis ij
axis([-1.1,1.1,150,1050])
line([0,0],[150,1050],'color','k')
title('FS: OSI changes in depth')
xlabel('(OSI light - OSI no light)/(OSI light + OSI no light)')
ylabel('cortical depth [mum]')

figure
plot((ldirprefratio(prs)-nldirprefratio(prs))./(ldirprefratio(prs)+nldirprefratio(prs)),depth(prs),'b.')
axis ij
axis([-1.1,1.1,150,1050])
line([0,0],[150,1050],'color','k')
title('RS: dsi changes in depth')
xlabel('(dsi light - dsi no light)/(dsi light + dsi no light)')
ylabel('cortical depth [mum]')

figure
plot((ldirprefratio(pfs)-nldirprefratio(pfs))./(ldirprefratio(pfs)+nldirprefratio(pfs)),depth(pfs),'r.')
axis ij
axis([-1.1,1.1,150,1050])
line([0,0],[150,1050],'color','k')
title('FS: dsi changes in depth')
xlabel('(dsi light - dsi no light)/(dsi light + DSI no light)')
ylabel('cortical depth [mum]')

figure
[s,p] = ttest(nloriprefratio(prs),loriprefratio(prs));
plot(nloriprefratio(prs),loriprefratio(prs),'b.')
line([0,1],[0,1],'color','k')
axis square
xlabel('OSI control');
ylabel('OSI light');
title(['RS cells: p: ' num2str(p)])

figure
[s,p] = ttest(nloriprefratio(pfs),loriprefratio(pfs));
plot(nloriprefratio(pfs),loriprefratio(pfs),'r.')
line([0,1],[0,1],'color','k')
axis square
xlabel('OSI control');
ylabel('OSI light');
title(['FS cells p: ' num2str(p)])

figure
[s,p] = ttest(nldirprefratio(prs),ldirprefratio(prs));
plot(nldirprefratio(prs),ldirprefratio(prs),'b.')
line([0,1],[0,1],'color','k')
axis square
xlabel('DSI control');
ylabel('DSI light');
title(['RS cells p: ' num2str(p)])

figure
[s,p] = ttest(nldirprefratio(pfs),ldirprefratio(pfs));
plot(nldirprefratio(pfs),ldirprefratio(pfs),'r.')
line([0,1],[0,1],'color','k')
axis square
xlabel('DSI control');
ylabel('DSI light');
title(['FS cells p: ' num2str(p)])

%%%%%%%%%%%%%%%%%%%%%%%%%% baseline subtracted preffr %%%%%%%%%%%%%%%%%%%%%%%%%% 
blspreffr(:,1) = preffr(:,1)-bfr'; blspreffr(:,2) = preffr(:,2)-bfr'; 
blspreffr(find(blspreffr(:,2)<0),2) = 0;

figure
plot(blspreffr(prs,1),depth(prs),'b.')
axis ij
% axis([-.1,1.1,150,1050])
title('RS: BLS Preferred FR in depth')
xlabel('Firing Rate [Hz]')
ylabel('cortical depth [mum]')

figure
plot(blspreffr(pfs,2),depth(pfs),'r.')
axis ij
% axis([-.1,1.1,150,1050])
title('FS: BLS Preferred FR in depth')
xlabel('Firing Rate [Hz]')
ylabel('cortical depth [mum]')

figure
plot((blspreffr(prs,2)-blspreffr(prs,1))./(blspreffr(prs,2)+blspreffr(prs,1)),depth(prs),'b.')
line([0,0],[0,1200],'color','k')
axis ij
axis([-1.1,1.1,0,1200])
ylabel('depth[mum]')
xlabel('(light-nolight)/(light+nolight)')
title('RS cells: BLS preferred firing rate changes by layer 4 suppression')

figure
plot((blspreffr(pfs,2)-blspreffr(pfs,1))./(blspreffr(pfs,2)+blspreffr(pfs,1)),depth(pfs),'r.')
line([0,0],[0,1000],'color','k')
axis ij
% axis([-1.1,1.1,0,1000])
ylabel('depth[mum]')
xlabel('(light-nolight)/(light+nolight)')
title('FS cells: BLS preferred firing rate changes by layer 4 suppression')

figure
plot(blspreffr(prs,1),blspreffr(prs,2),'b.')
line([0,50],[0,50],'color','k')
axis square
xlabel('preferred firing rate no light')
ylabel('preferred firing rate light on')
title('RS cells: BLS pref firing rate changes');

figure
plot(blspreffr(pfs,1),blspreffr(pfs,2),'r.')
line([0,80],[0,80],'color','k')
axis square
xlabel('preferred firing rate no light')
ylabel('preferred firing rate light on')
title('FS cells: BLS pref firing rate changes');

%%%%%%%%%%%%%%%%%%%%%%%% preferred firing rates %%%%%%%%%%%%%%%%%%%%%%
figure
plot(preffr(prs,1),depth(prs),'b.')
axis ij
% axis([-.1,1.1,150,1050])
title('RS: Preferred FR in depth')
xlabel('Firing Rate [Hz]')
ylabel('cortical depth [mum]')

figure
plot(preffr(pfs,2),depth(pfs),'r.')
axis ij
% axis([-.1,1.1,150,1050])
title('FS: Preferred FR in depth')
xlabel('Firing Rate [Hz]')
ylabel('cortical depth [mum]')

figure
plot((preffr(prs,2)-preffr(prs,1))./(preffr(prs,2)+preffr(prs,1)),depth(prs),'b.')
line([0,0],[0,1000],'color','k')
axis ij
axis([-1.1,1.1,0,1000])
ylabel('depth[mum]')
xlabel('(light-nolight)/(light+nolight)')
title('RS cells: preferred firing rate changes by layer 4 suppression')

figure
plot((preffr(pfs,2)-preffr(pfs,1))./(preffr(pfs,2)+preffr(pfs,1)),depth(pfs),'r.')
line([0,0],[0,1000],'color','k')
axis ij
axis([-1.1,1.1,0,1000])
ylabel('depth[mum]')
xlabel('(light-nolight)/(light+nolight)')
title('FS cells: preferred firing rate changes by layer 4 suppression')

figure
plot(preffr(prs,1),preffr(prs,2),'b.')
line([0,50],[0,50],'color','k')
axis square
xlabel('preferred firing rate no light')
ylabel('preferred firing rate light on')
title('RS cells: pref firing rate changes');

figure
plot(preffr(pfs,1),preffr(pfs,2),'r.')
line([0,80],[0,80],'color','k')
axis square
xlabel('preferred firing rate no light')
ylabel('preferred firing rate light on')
title('FS cells: pref firing rate changes');

%%%%%%%%%%%%%%%%%%%%%%%% average firing rates %%%%%%%%%%%%%%%%%%%%%%
figure
plot(nlfr(prs),depth(prs),'b.')
axis ij
% axis([-.1,1.1,150,1050])
title('RS: Average FR in depth')
xlabel('Firing Rate [Hz]')
ylabel('cortical depth [mum]')

figure
plot(nlfr(pfs),depth(pfs),'r.')
axis ij
% axis([-.1,1.1,150,1050])
title('FS: Average FR in depth')
xlabel('Firing Rate [Hz]')
ylabel('cortical depth [mum]')

figure
plot((lfr(prs)-nlfr(prs))./(lfr(prs)+nlfr(prs)),depth(prs),'b.')
line([0,0],[0,1000],'color','k')
axis ij
axis([-1.1,1.1,0,1000])
ylabel('depth[mum]')
xlabel('(light-nolight)/(light+nolight)')
title('RS cells: average firing rate changes by layer 4 suppression')

figure
plot((lfr(pfs)-nlfr(pfs))./(lfr(pfs)+nlfr(pfs)),depth(pfs),'r.')
line([0,0],[0,1000],'color','k')
axis ij
axis([-1.1,1.1,0,1000])
ylabel('depth[mum]')
xlabel('(light-nolight)/(light+nolight)')
title('FS cells: average firing rate changes by layer 4 suppression')

figure
plot(nlfr(prs),lfr(prs),'b.')
line([0,50],[0,50],'color','k')
axis square
xlabel('preferred firing rate no light')
ylabel('preferred firing rate light on')
title('RS cells:average firing rate changes');

figure
plot(nlfr(pfs),lfr(pfs),'r.')
line([0,80],[0,80],'color','k')
axis square
xlabel('preferred firing rate no light')
ylabel('preferred firing rate light on')
title('FS cells: average firing rate changes');

%%%%%%%%%%%%%%%% Control FR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
plot(controlfr(prs,1),depth(prs),'b.')
axis ij
% axis([-.1,1.1,150,1050])
title('RS: Control FR in depth')
xlabel('Firing Rate [Hz]')
ylabel('cortical depth [mum]')

figure
plot(controlfr(pfs,2),depth(pfs),'r.')
axis ij
% axis([-.1,1.1,150,1050])
title('FS: Control FR in depth')
xlabel('Firing Rate [Hz]')
ylabel('cortical depth [mum]')

figure
plot((controlfr(prs,2)-controlfr(prs,1))./(controlfr(prs,2)+controlfr(prs,1)),depth(prs),'b.')
line([0,0],[0,1000],'color','k')
axis ij
axis([-1.1,1.1,0,1000])
ylabel('depth[mum]')
xlabel('(light-nolight)/(light+nolight)')
title('RS cells: control firing rate changes by layer 4 suppression')

figure
plot((controlfr(pfs,2)-controlfr(pfs,1))./(controlfr(pfs,2)+controlfr(pfs,1)),depth(pfs),'r.')
line([0,0],[0,1000],'color','k')
axis ij
axis([-1.1,1.1,0,1000])
ylabel('depth[mum]')
xlabel('(light-nolight)/(light+nolight)')
title('FS cells: control firing rate changes by layer 4 suppression')

figure
plot(controlfr(prs,1),controlfr(prs,2),'b.')
line([0,30],[0,30],'color','k')
axis square
xlabel('control firing rate no light')
ylabel('control firing rate light on')
title('RS cells: control firing rate changes');

figure
plot(controlfr(pfs,1),controlfr(pfs,2),'r.')
line([0,30],[0,30],'color','k')
axis square
xlabel('control firing rate no light')
ylabel('control firing rate light on')
title('FS cells: control firing rate changes');


figure
[x,y,err] = runningAverage(depth(prs),sinl(prs),0,15);
plot(sinl(prs),depth(prs),'b.')
axis ij
hold on
plot(x,y,'k','linewidth',2)
plot(x+err,y,'k')
plot(x-err,y,'k')
axis([-.05,1.05,200,1050])
xlabel('suppression index ((max-largest)/max)')
ylabel('depth [mum]')
title('Suppression Index RS cells')

figure
[x,y,err] = runningAverage(depth(pfs),sinl(pfs),0,15);
plot(sinl(pfs),depth(pfs),'r.')
axis ij
hold on
plot(x,y,'k','linewidth',2)
plot(x+err,y,'k')
plot(x-err,y,'k')
axis([-.05,1.05,200,1050])
xlabel('suppression index ((max-largest)/max)')
ylabel('depth [mum]')
title('Suppression Index FS cells')

figure
plot(sinl(prs),sil(prs),'b.')
axis([0,1,0,1])
axis square
line([0,1],[0,1],'color','k')
xlabel('SI ((max-largest)/max) no light')
ylabel('SI ((max-largest)/max) light')
[s,p] = ttest(sinl(prs),sil(prs))
title(['SI: RS cells, p: ' num2str(p)])

figure
plot(sinl(pfs),sil(pfs),'r.')
axis([0,1,0,1])
axis square
line([0,1],[0,1],'color','k')
xlabel('SI ((max-largest)/max) no light')
ylabel('SI ((max-largest)/max) light')
[s,p] = ttest(sinl(pfs),sil(pfs))
title(['SI: FS cells, p: ' num2str(p)])

figure
plot((sil(prs)-sinl(prs))./(sil(prs)+sinl(prs)),depth(prs),'bo','markersize',4,'markerfacecolor','b')
hold on
[x,y,xerr] = runningMedian(depth(prs),(sil(prs)-sinl(prs))./(sil(prs)+sinl(prs)),0,12);
plot(x,y,'b','linewidth',2);
plot(x+xerr,y,'b')
plot(x-xerr,y,'b')
line([-1,1,],[375,375],'color','k','linestyle',':')
line([-1,1,],[500,500],'color','k','linestyle',':')
axis ij
line([0,0],[150,1050],'color','k')
axis([-1.1,1.1,150,1050])
line([0,0],[150,1000],'color','k')
ylabel('depth[mum]')
xlabel('(SI light - SI no light)/(SI light + SI no light)')

figure
plot((sil(pfs)-sinl(pfs))./(sil(pfs)+sinl(pfs)),depth(pfs),'ro','markersize',4,'markerfacecolor','r')
hold on
[x,y,xerr] = runningMedian(depth(pfs),(sil(pfs)-sinl(pfs))./(sil(pfs)+sinl(pfs)),0,12);
plot(x,y,'r','linewidth',2);
plot(x+xerr,y,'r')
plot(x-xerr,y,'r')
line([-1,1,],[375,375],'color','k','linestyle',':')
line([-1,1,],[500,500],'color','k','linestyle',':')
axis ij
axis([-1.1,1.1,150,1050])
line([0,0],[150,1000],'color','k')
ylabel('depth[mum]')
xlabel('(SI light - SI no light)/(SI light + SI no light)')

end

function [ftspect,ftphases,ftfax,ftralp,ftppc,ftplv] = get_ft_spectstats(lfpresp,resp,thisinds)
 
     if isempty(thisinds) || isempty(find(resp(thisinds,801:1800)))
         ftspect = nan(1,30); ftphases = NaN; ftfax = nan(1,30);
         ftralp = nan(1,30); ftppc = nan(1,30); ftplv = nan(1,30);
     else
         clear data;
         for i = 1:length(thisinds)
             data.trial{i} = [lfpresp(thisinds(i),:);resp(thisinds(i),:)];
             data.time{i} = -.299:.001:2.700;
         end
         data.fsample = 1000;
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
         ftspect = squeeze(nanmean(mag(:,1,:)));
         ftphases = squeeze(ang);
         ftfax = stsFFT.freq;

         cfg               = [];
         cfg.method        = 'ral'; % compute the rayleigh test
         cfg.spikechannel  = stsFFT.label{1};
         cfg.channel       = stsFFT.lfplabel; % selected LFP channels
         cfg.avgoverchan   = 'unweighted'; % weight spike-LFP phases irrespective of LFP power
         cfg.timwin        = 'all'; % compute over all available spikes in the window
         cfg.latency       = [0.5 1.5]; % sustained visual stimulation period
         statSts           = ft_spiketriggeredspectrum_stat(cfg,stsFFT);
         ftralp = statSts.ral;
         cfg.method = 'ppc0';
         statSts           = ft_spiketriggeredspectrum_stat(cfg,stsFFT);
         ftppc = statSts.ppc0;
         cfg.method = 'plv';
         statSts           = ft_spiketriggeredspectrum_stat(cfg,stsFFT);
         ftplv = statSts.plv;
     end
end

function [p, resnorm, rsquared] = fit_fixedfsin(x,y,f,sr)
    %[A, ph, offs]
    range = max(y)-min(y);
    if range~=0
        p0 = [range/2, pi, (max(y)+min(y))/2]; 
        lb = [.5*(range/2),0,min(y)];
        ub = [1.5*(range/2),2*pi,max(y)];
        [p,resnorm,residual,exitflag] = lsqcurvefit(@(p,x) fixedfsin(p,x,f,sr), p0,x,y,lb,ub,optimset('Display','off'));
        restot = sum((y-mean(y)).^2);
        rsquared = 1 - (resnorm/restot);
    else
        p = [NaN,NaN,NaN];
        resnorm = NaN;
        rsquared = NaN;
    end
end

function val = fixedfsin(p,x,f,sr)
    %[A f ph offs]
    val = p(1) * sin(x*((f*2*pi)/sr) + p(2)) + p(3);
end

function depthplot(value,depth,formatstr,formatstr2)
plot(value,depth,formatstr)
axis ij
line([0,0],[0,1200],'color','k')
hold on
[x,y,xerr] = runningAverage(depth,value,0,7);
plot(x,y,formatstr2,'linewidth',2)
plot(x-xerr,y,formatstr2)
plot(x+xerr,y,formatstr2)
end

function out = uncircle(in)
    in(find(in<0)) = in(find(in<0))+2*pi;
    out = in;
end

% difference in directions, maximally 180
function a = dirdist(x,y)
    a = abs(x-y);
    a(a>180) = abs(a(a>180)-360);
end

% always positive angular difference
function a = oridist(x,y)
    a = abs(x-y);
    a(find(a>90)) = abs(a(find(a>90))-180);
    a(find(a>90)) = abs(a(find(a>90))-180);
end

% directed angular difference for 180 deg
function a = oridiff(x,y)
    a = x-y;
    a(find(a>90)) = -abs(a(find(a>90))-180);
    a(find(a<-90)) = abs(abs(a(find(a<-90)))-180);
end
    

function exampleplot(condresp,conderr,controlfr,controlerr,bta,sizex,orix)

subplot(2,2,1)
plot(bta,squeeze(mean(mean(condresp(1,:,:,:),2),3)),'k','linewidth',2);
hold on
plot(bta,squeeze(mean(mean(condresp(2,:,:,:),2),3)),'r','linewidth',2)
ax = axis;
axis([-300,2600,ax(3),ax(4)])
line([-300,2600],[0,0],'color','k')
line([500,500],[ax(3),ax(4)],'color','k')
line([1500,1500],[ax(3),ax(4)],'color','k')

s1 = find(bta>500,1); s2 = find(bta>1500,1);
frsl0 = squeeze(mean(condresp(1,:,:,s1:s2),4));
frsl1 = squeeze(mean(condresp(2,:,:,s1:s2),4));
ps = find(mean(frsl0) == max(mean(frsl0)));
vfrl0 = reshape(frsl0,1,40);
vfrl1 = reshape(frsl1,1,40);

subplot(2,2,2)
plot(frsl0,frsl1,'o','markersize',3)
refline(1,0);

subplot(2,2,3)
errorbar(sizex,[controlfr(1),mean(frsl0)],[controlerr(1),squeeze(mean(conderr(1,:,:),2))'],'k','linewidth',2)
hold on
errorbar(sizex,[controlfr(2),mean(frsl1)],[controlerr(2),squeeze(mean(conderr(2,:,:),2))'],'r','linewidth',2)
ax = axis;
axis([-5,50,ax(3),ax(4)])

subplot(2,2,4)
errorbar(orix,frsl0(:,ps),squeeze(conderr(1,:,ps)),'k','linewidth',2);
hold on
errorbar(orix,frsl1(:,ps),squeeze(conderr(2,:,ps)),'r','linewidth',2);
ax = axis;
axis([-45,350,ax(3),ax(4)])
end
