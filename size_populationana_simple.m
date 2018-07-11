function size_populationana_simple


% % Scnn Halo
% %SCNN Population
% animalids = {'140618', '140619', '140620', '140703', '140703', '140704', '140709', '140709', '140711', '140715', '140717', '140717', '140806', '140806', '141024', '141028', '141103', '141110', '141112', '141112', '141113', '141201', '141201b', '141204', '141209'};
% blocks    = [3,         7,        4,        8,        14,       5,        5,        8,        5,        8,        3,        14        4         8,        10,       8,        7,        3,        4,        9,        4,        2,        2,         3,        2];
% animal =    [1          2         3         4         4         5         6         6         7,        8         9,        9         10        10,       11,       12,       13,       14,       15,       15,       16,       17,       18,        19,       20];
% electrodes= [[1,32];   [1,32];    [1,32];   [1,32];   [1,32];   [1,32];   [1,32];   [1,32];   [1,32];   [1,32];   [1,32];   [1,32];   [1,32];   [1,32];   [1,32];   [1,32];   [1,16];   [1,16];   [1,16];   [1,16];   [1,16];   [1,16];   [1,16];    [1,16];   [1,16];];
% penangle =  [10,        10,       10,       10,       10,       10,       10,       10,       10,       10,       10,       10,       10,       10,       10,       10,       10,       10,       10,       10,       10,       10,       10,        10,       10];
% printpath = 'C:\Users\Julia\work\data\populations\Scnn_Halo_withnew\size_simple\units\';
% runprintpath = 'C:\Users\Julia\work\data\populations\Scnn_Halo_withnew\size_simple\running\';
% l5RSprintpath = 'C:\Users\Julia\work\data\populations\Scnn_Halo_withnew\size_simple\l5rs\'
% popfile = 'C:\Users\Julia\work\data\populations\Scnn_Halo_withnew\size_simple\size_population_withnewanimals.mat';

% % Scnn Halo
% %SCNN Population scotts paper
% animalids = {'140618', '140619', '140620', '140703', '140703', '140704', '140709', '140709', '140711', '140715', '140717', '140717', '140806', '140806'};
% blocks    = [3,         7,        4,        8,        14,       5,        5,        8,        5,        8,         3,        14        4         8];
% animal =    [1          2         3         4         4         5         6         6         7,        8          9,        9         10        10];
% electrodes= [32,        32,       32,       32,       32,       32,       32,       32,       32,       32,        32,       32,       32,       32];
% penangle =  [10,        10,       10,       10,       10,       10,       10,       10,       10,       10,        10,       10,       10,       10];
% printpath = 'C:\Users\Julia\work\data\populations\size_simple\units\';
% runprintpath = 'C:\Users\Julia\work\data\populations\size_simple\running\';
% l5RSprintpath = 'C:\Users\Julia\work\data\populations\size_simple\l5rs\'
% popfile = 'C:\Users\Julia\work\data\populations\size_simple\size_population.mat';

% %SCNN-CHR2 Population 
% animalids = {'150323', '150325'};
% blocks    = [6,         4];
% animal =    [1          2];
% penangle = 10;
% printpath = 'C:\Users\Julia\work\data\populations\scnn_chr2\size_simple\units\';
% runprintpath = 'C:\Users\Julia\work\data\populations\scnn_chr2\size_simple\running\';
% condprintpath = 'C:\Users\Julia\work\data\populations\scnn_chr2\size_simple\conditions\';
% l5RSprintpath = 'C:\Users\Julia\work\data\populations\scnn_chr2\size_simple\l5rs\';
% popfile = 'C:\Users\Julia\work\data\populations\scnn_chr2\size_simple\size_population.mat';

% % SOM P0 population
% animalids = {'140807', '140807', '140812', '140812', '140815', '140815', '140818', '141014', '141014', '141015', '141106', '141110b', '141111', '150123', '150213','150523'};
% blocks    = [5,         8,        3,        5,        4,        6,        5,        3,        10,       5,        3,        7,         4,        6,        3,       3];
% animal =    [1          1         2         2         3         3         4,        5,        5,        6,        7,        8,         9,        10,       11,      12];
% electrodes =[32,        32,       32,       32,       32,       32,       32,       32,       32,       32,       32,       32,        32,       32,       16,      32];
% penangle =  [10,        10,       10,       10,       10,       10,       10,       10,       10,       10,       10,       10,        10,       10,       20,      25];
% printpath = 'C:\Users\Julia\work\data\populations\SOM_Halo\size_simple\units\';
% runprintpath = 'C:\Users\Julia\work\data\populations\SOM_Halo\size_simple\running\';
% popfile = 'C:\Users\Julia\work\data\populations\SOM_Halo\size_simple\size_population.mat';

% % undecided about 150407 and 150406 (P1 morning injections, not much light effect)
% % % 
% % % SOM later population
% animalids = {'150331', '150401','150527','150529','150602','150603','150625','150825','150831','150902','150907','150909', '150915', '150916', '151023', '151027','151109', '151110', '151209', '160122'};
% blocks    = [3,         5,       11,      4,       5,       3,       6,       5,       4,       3,       3,       4,        3,        2,       14,       3,        11,       13,       6,        3];
% animal    = [1,         2,       3,       4,       5,       6,       7,       8,       9,       10,      11,      12,       13,       14,      15,       16,       17,       18,       19,       20];
% electrodes =[[1,32];    [1,32];  [1,32];  [1,32];  [1,32];  [1,16];  [1,16];  [17,32]; [1,16];  [1,16];  [1,16];  [1,16];   [17,32];  [1,16];  [17,32];  [17,32]; [17,32];  [1,16];   [17,32];  [1,16]];
% penangle =  [25,        25,      25,      25,      25,      10,      10,      25,      25,      25,      25,      25,       25,       25,      25,       25,       25,       25,       25,       25];
% % age       [P6,        P6,      P2?(P0), P2?(P0), P2?(P0), P1
% printpath = 'C:\Users\Julia\work\data\populations\SOM_Halo_later\size_simple\units\';
% runprintpath = 'C:\Users\Julia\work\data\populations\SOM_Halo_later\size_simple\running\';
% oscillprintpath = 'C:\Users\Julia\work\data\populations\SOM_Halo_later\size_simple\withlfpoverview\';
% popfile = 'C:\Users\Julia\work\data\populations\SOM_Halo_later\size_simple\size_population.mat';

% % SOM Halo+ChR2 % 150914 not all sizes - up to 9 all Halo then CHr2
% animalids = {'150331', '150401','150602','150625','150909', '151023', '151209', '160122', '151229', '160311', '160405_2'};
% blocks    = [3,         5,       5,       6,       4,        14,       6,        3,        10,       7,        3];
% animal    = [1,         2,       3,       4,       5,        6,        7,        8,        9,        10,       11];
% electrodes =[[1,32];    [1,32];  [1,32]; [1,16];  [1,16];   [17,32];  [17,32];  [1,16];   [1,16];   [1,16];   [17,32]];
% penangle =  [25,        25,      25,      10,      25,       25,       25,       25,       25,       25,       25];
% printpath = 'C:\Users\Julia\work\data\populations\SOM_HaloandChR2\size_simple\units\';
% runprintpath = 'C:\Users\Julia\work\data\populations\SOM_HaloandChR2\size_simple\running\';
% oscillprintpath = 'C:\Users\Julia\work\data\populations\SOM_HaloandChR2\size_simple\withlfpoverview\';
% popfile = 'C:\Users\Julia\work\data\populations\SOM_HaloandChR2\size_simple\size_population.mat';

% % PV-earch population
% animalids = {'151117', '151211', '160114', '160115'};
% blocks    = [2,         3,        4,        3];
% animal    = [1,         2,        3,        4];
% electrodes =[[1,16];    [1,16];  [1,16];   [1,16]];
% penangle =  [25,        25,       25,       25];
% printpath = 'C:\Users\Julia\work\data\populations\PV_eArch\size_simple\units\';
% runprintpath = 'C:\Users\Julia\work\data\populations\PV_eArch\size_simple\running\';
% oscillprintpath = 'C:\Users\Julia\work\data\populations\PV_eArch\size_simple\withlfpoverview\';
% popfile = 'C:\Users\Julia\work\data\populations\PV_eArch\size_simple\size_population.mat';

% % control SOM only population
% animalids = {'151214', '151215', '151217', '160112', '160113'};
% blocks    = [2,         7,        3,        3,        6];
% animal    = [1,         2,        3,        4,        5];
% electrodes =[[1,16];   [1,16];   [1,16];   [1,16];   [1,16]];
% penangle =  [25,        25,       25,       25,       25];
% printpath = 'C:\Users\Julia\work\data\populations\control\size_simple\units\';
% runprintpath = 'C:\Users\Julia\work\data\populations\control\size_simple\running\';
% oscillprintpath = 'C:\Users\Julia\work\data\populations\control\size_simple\withlfpoverview\';
% popfile = 'C:\Users\Julia\work\data\populations\control\size_simple\size_population.mat';
% 
% % anesthetized SOM only population
% animalids = {'150826', '151230', '151231'};
% blocks    = [3,         3,        11];
% animal    = [1,         2,        3];
% electrodes =[[1,16];   [1,16];   [1,16]];
% penangle =  [25,        25,       25];
% printpath = 'C:\Users\Julia\work\data\populations\SOM_Halo_anesth\size_simple\units\';
% runprintpath = 'C:\Users\Julia\work\data\populations\SOM_Halo_anesth\size_simple\running\';
% oscillprintpath = 'C:\Users\Julia\work\data\populations\SOM_Halo_anesth\size_simple\withlfpoverview\';
% popfile = 'C:\Users\Julia\work\data\populations\SOM_Halo_anesth\size_simple\size_population.mat';

% % anesthetized SOM+PV combined only population
% animalids = {'150826', '151230', '151231', '150821', '150824'};
% blocks    = [3,         3,        11,       3,        18];
% animal    = [1,         2,        3,        4,        5];
% electrodes =[[1,16];   [1,16];   [1,16];   [1,16];   [17,32]];
% penangle =  [25,        25,       25,       25,       25];
% printpath = 'C:\Users\Julia\work\data\populations\SOMPVcombined_anesth\size_simple\units\';
% runprintpath = 'C:\Users\Julia\work\data\populations\SOMPVcombined_anesth\size_simple\running\';
% oscillprintpath = 'C:\Users\Julia\work\data\populations\SOMPVcombined_anesth\size_simple\withlfpoverview\';
% popfile = 'C:\Users\Julia\work\data\populations\SOMPVcombined_anesth\size_simple\size_population.mat';

% % PV Halo population
% animalids = {'150629', '150730', '150731', '150804','150818','150820','150823','150824', '151104'};
% blocks    = [4,         4,        3,         4,      3,       3,       5,       3,        5];
% animal    = [1,         2,        3,         4,      5,       6,       7,       8,        9];
% electrodes =[[1,16];    [1,16];   [1,16];   [1,16];  [1,16];  [1,16];  [1,16];  [17,32];  [17,32]];
% penangle =  [10,        10,       25,        25,     25,      25,      25,      25,       25];
% printpath = 'C:\Users\Julia\work\data\populations\PV_Halo\size_simple\units\';
% runprintpath = 'C:\Users\Julia\work\data\populations\PV_Halo\size_simple\running\';
% oscillprintpath = 'C:\Users\Julia\work\data\populations\PV_Halo\size_simple\withlfpoverview\';
% popfile = 'C:\Users\Julia\work\data\populations\PV_Halo\size_simple\size_population.mat';
% 
% HUGEass combined population of all animals - 1:20 SOM, 21:32 PV Halo 33:39 PV eArch
animalids = {'150331', '150401','150527','150529','150602','150603','150625','150825','150831','150902','150907','150909', '150915', '150916', '151023', '151027','151109', '151110', '151209', '160122', '150629', '150730', '150731', '150804','150818','150820','150823','150824', '151104', '160217', '160324', '160328', '151117', '151211', '160114', '160115', '160204', '160205', '160210'};
blocks    = [3,         5,       11,      4,       5,       3,       6,       5,       4,       3,       3,       4,        3,        2,       14,       3,        11,       13,       6,        3,        4,        4,        3,        4,       3,       3,       5,       3,        5,        5,        2,        2,        2,        3,        4,        3,        3,        4,        9];
animal    = [1,         2,       3,       4,       5,       6,       7,       8,       9,       10,      11,      12,       13,       14,      15,       16,       17,       18,       19,       20,       21,       22,       23,       24,      25,      26,      27,      28,       29,       30,       31,       32,       33,       34,       35,       36,       37,       38,       39];
rfblocks =  [1,         1,       10,      1,       4,       2,       4,       3,       3,       1,       1,       1,        2,        1,       9,        2,        6,        9,        3,        1,        2,        2,        2,        2,       2,       1,       3,       2,        4,        4,        1,        1,        1,        2,        3,        1,        1,        2,        5];
electrodes =[[1,32];    [1,32];  [1,32];  [1,32];  [1,32];  [1,16];  [1,16];  [17,32]; [1,16];  [1,16];  [1,16];  [1,16];   [17,32];  [1,16];  [17,32];  [17,32];  [17,32]; [1,16];   [17,32];  [1,16];   [1,16];   [1,16];   [1,16];   [1,16];  [1,16];  [1,16];  [1,16];  [17,32];  [17,32];  [1, 16];  [1, 16];  [1, 16];  [1,16];    [1,16];  [1,16];   [1,16];   [1,16];   [17,32];  [17,32]];
penangle =  [25,        25,      25,      25,      25,      10,      10,      25,      25,      25,      25,      25,       25,       25,      25,       25,       25,       25,       25,       25,       25,       25,       25,       25,      25,      25,      25,      25,       25,       25,       25,       25,       25,        25,      25,       25,       25,       25,       25];
% age       [P6,        P6,      P2?(P0), P2?(P0), P2?(P0), P1
printpath = 'C:\Users\Julia\work\data\populations\SOMPVcombined\size_simple\units\';
runprintpath = 'C:\Users\Julia\work\data\populations\SOMPVcombined\size_simple\running\';
oscillprintpath = 'C:\Users\Julia\work\data\populations\SOMPVcombined\size_simple\withlfpoverview\';
popfile = 'C:\Users\Julia\work\data\populations\SOMPVcombined\size_simple\size_population.mat';

% % PV Halo + PV eArch population, up to animal 12 all Halo then eArch
% animalids = {'150629', '150730', '150731', '150804','150818','150820','150823','150824', '151104', '160217', '160324', '160328', '151117', '151211', '160114', '160115', '160204', '160205', '160210'};
% blocks    = [4,         4,        3,         4,      3,       3,       5,       3,        5,        5,        2,        2,        2,         3,        4,        3,       3,        4,        9];
% animal    = [1,         2,        3,         4,      5,       6,       7,       8,        9,        10,       11,       12,       13,        14,       15,       16,      17,       18,       19];
% rfblocks =  [2,         2,        2,         2,      2,       1,       3,       2,        2,        4,        1,        1,        1,         2,        3,        1,       1,        5,        5];
% electrodes =[[1,16];    [1,16];   [1,16];   [1,16];  [1,16];  [1,16];  [1,16];  [17,32];  [17,32]; [1, 16];  [1, 16];  [1, 16];  [1,16];    [1,16];   [1,16];   [1,16];   [1,16];   [17,32];  [17,32]];
% penangle =  [10,        10,       25,        25,     25,      25,      25,      25,       25,       25,       25,       25,       25,        25,       25,       25,      25,       25,       25];
% printpath = 'C:\Users\Julia\work\data\populations\PV_HaloeArch\size_simple\units\';
% runprintpath = 'C:\Users\Julia\work\data\populations\PV_HaloeArch\size_simple\running\';
% oscillprintpath = 'C:\Users\Julia\work\data\populations\PV_HaloeArch\size_simple\withlfpoverview\';
% popfile = 'C:\Users\Julia\work\data\populations\PV_HaloeArch\size_simple\size_population.mat';

% % X94 CREDOG Halo ?
% animalids = {'170425', '170426','170427','170711','170712'};
% blocks    = [4,         4,       2,       3,       4];
% animal    = [1,         2,       3,       4,       5];
% rfblocks =  [6,         3,       3,       1,       2];
% electrodes =[[1,16];    [1,32];  [1,32]; [1,16];  [1,16]];
% penangle =  [25,        25,      25,      30,      30];
% % age       [P6,        P6,      P2?(P0), P2?(P0), P2?(P0), P1
% printpath = 'C:\Users\Julia\work\data\populations\X94_Halo_V1\size_simple\units\';
% runprintpath = 'C:\Users\Julia\work\data\populations\X94_Halo_V1\size_simple\running\';
% oscillprintpath = 'C:\Users\Julia\work\data\populations\X94_Halo_V1\size_simple\withlfpoverview\';
% popfile = 'C:\Users\Julia\work\data\populations\X94_Halo_V1\size_simple\size_population.mat';


tic
lcol = 'r'; %lasercolor

recalculate = 1;
printyn = 1

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
        snbasename = [animalids{blck} '_block' int2str(rfblocks(blck)) '_tet'];

        files = dir([supath, basename, '*.mat']);
        rffiles = dir([supath, snbasename, '*.mat']);

        prestim = 300;
        poststim = 700;
        respwin = 501:1500; % after stimulus onset
        offsetwin = 1501:2500;
        respwin = respwin+prestim;
        offsetwin = offsetwin+prestim;
        freqbinwidth = 5;
        
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
            
            wvchan = find(var(result.waveforms) == max(var(result.waveforms)));
            lfp = result.lfp(:,wvchan)';
                        
            trialdur = result.stimduration*1000;
                        
            if length(result.msstamps)~=length(result.light)
%                 disp('');
%             msstamps(16) = []; % for 150414 block 10
%             result.msstamps = msstamps;
%             save([supath, SUfiles(cl).name],'result');            
                pause;
            end
            
            % fix so it is usable with new multi-purpose grating stim script
            if isfield(result, 'sizeconds')
                allinds = sort(getSpecificIndices(result, 'sizeconds'));
                msstamps = result.msstamps(allinds);
                light = result.light(allinds);
                gratingInfo.Orientation = result.gratingInfo.Orientation(allinds);
                gratingInfo.size = result.gratingInfo.size(allinds);
                gratingInfo.Contrast = result.gratingInfo.Contrast(allinds);
                gratingInfo.tFreq = result.gratingInfo.tFreq(allinds);
            else
                msstamps = result.msstamps;
                light = result.light;
                gratingInfo = result.gratingInfo;
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
            
            % not so great running correlations
            a = find(msStimes>500&msStimes<msstamps(1)-500);
            b = find(msStimes>msstamps(end)+2500&msStimes<length(result.runspeed)-500);
            strs = nan(1,1001);
            if ~isempty(a)
                for i = 1:length(a)
                    strs(i,:) = result.runspeed(msStimes(a(i))-500:msStimes(a(i))+500);
                end
            end
            n = length(a);
            if ~isempty(b)
                for i = 1:length(b)
                    strs(n+i,:) = result.runspeed(msStimes(b(i))-500:msStimes(b(i))+500);
                end
            end
            if size(strs,1)>1
                stra(cll,:) = nanmean(strs);
                clear strs; 
            end
            a = filter(excit_kernel,1,chan);
            spksg = [a(1:msstamps(1)),a(msstamps(end)+2000:end-1)];
            runsg = [result.runspeed(1:msstamps(1)),result.runspeed(msstamps(end)+2000:end)];
            srcorr(cll) = nancorr(spksg',runsg');
            
            for i = 1:length(msstamps)
                resp(i,:) = chan(msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim);  
                lfpresp(i,:) = lfp(msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim);               
                msl(i) = mean(find(resp(i,respwin)));
            end            
            msta = linspace(-prestim,trialdur+poststim,size(resp,2));

            frs = sum(resp(:,respwin),2)./(length(respwin)/1000);
            bl = sum(resp(:,1:prestim),2)./(prestim/1000);
            sc = sum(resp(:,respwin),2); % spike count
                        
            l0s = find(result.light == 0);            
            firingrates(cll,:) = frs(l0s);
            speeds(cll,:) = meanspeed(l0s);
            
            %get RF of cell, too
            [onfield(cll,:,:),offfield(cll,:,:),gaussfiton(cll),gaussfitoff(cll),rson(cll),rsoff(cll),xax(cll,:),yax(cll,:)] = get_rf([supath, rffiles(fi).name]);
            
            
            binwidth = 30;            
            sizes = unique(gratingInfo.size);  sizes(find(sizes == 0)) = []; %delete control condition
            oris = unique(gratingInfo.Orientation);  oris(find(oris == -1)) = [];
            shownsizes(cll,:) = sizes;                      
            for l = 1:length(unique(light))
                for sz = 1:length(sizes)
                    for ori = 1:length(oris)
                        thisinds = find(gratingInfo.Orientation == oris(ori) &...
                            gratingInfo.size == sizes(sz) & ...
                            light == l-1);
                        condresp(cll,l,ori,sz,:) = mean(resp(thisinds,:),1);
                        condlfpresp(cll,l,ori,sz,:) = mean(lfpresp(thisinds,:),1);
                        condfr(cll,l,ori,sz) = mean(frs(thisinds));%-mean(bl); 
                        conderr(cll,l,ori,sz) =std(frs(thisinds))./sqrt(length(thisinds));
                        allfr{l,ori,sz} = frs(thisinds);
                        [bincondresp(cll,l,ori,sz,:),bta] = binit(squeeze(condresp(cll,l,ori,sz,:)),binwidth); 
                        condmsl(cll,l,ori,sz) = mean(msl(thisinds)); % mean spike latency
                        %TODO fix to be Hz
                        condfiltresp(cll,l,ori,sz,:) = filter(excit_kernel,1,mean(resp(thisinds,:),1));
                        
                        % test each cell for each condition ranksum
                        if l == 1
                            l1inds = find(gratingInfo.Orientation == oris(ori) &...
                            gratingInfo.size == sizes(sz) & ...
                            light == 1);
                            if ~isempty(l1inds) & ~isempty(thisinds)
                                condlightmodp(cll,ori,sz) = ranksum(frs(thisinds),frs(l1inds));
                            else
                                condlightmodp(cll,ori,sz) = NaN;
                            end
                        end
                        mscc = []; bincc = [];
                        for ii = 1:length(thisinds)-1
                            for jj = ii+1:length(thisinds)
                                help = corrcoef(resp(thisinds(ii),:),resp(thisinds(jj),:));
                                mscc = [mscc,help(1,2)];
                                help = corrcoef(binit(resp(thisinds(ii),:),binwidth),binit(resp(thisinds(jj),:),binwidth));
                                bincc = [bincc, help(1,2)];
                            end
                        end
                        msreliab(cll,l,ori,sz) = nanmean(mscc);
                        binreliab(cll,l,ori,sz) = nanmean(bincc);
                        eckerreliability(cll,l,ori,sz) = var(frs(thisinds))/var(frs);
                             
                        % running
                        thisruninds = intersect(thisinds,oktrials);
                        if ~isempty(thisruninds)
                            runcondresp(cll,l,ori,sz,:) = mean(resp(thisruninds,:),1);
                            runcondfr(cll,l,ori,sz) = mean(frs(thisruninds));
                            runconderr(cll,l,ori,sz) = std(frs(thisruninds))./sqrt(length(thisruninds));
                        else
                            runcondresp(cll,l,ori,sz,:) = nan(1,size(resp,2));
                            runcondfr(cll,l,ori,sz) = NaN;
                            runconderr(cll,l,ori,sz) = NaN;
                        end
                        
                        thisstillinds = intersect(thisinds,stilltrials);
                        if ~isempty(thisstillinds)
                            stillcondresp(cll,l,ori,sz,:) = nanmean(resp(thisstillinds,:),1);
                            stillcondfr(cll,l,ori,sz) = nanmean(frs(thisstillinds));
                            stillconderr(cll,l,ori,sz) = nanstd(frs(thisstillinds))./sqrt(length(thisstillinds));
                        else
                            stillcondresp(cll,l,ori,sz,:) = nan(1,size(resp,2));
                            stillcondfr(cll,l,ori,sz) = NaN;
                            stillconderr(cll,l,ori,sz) = NaN;
                        end
                    end
                end
            end

            bincondresp(cll,:,:,:,:) = bincondresp(cll,:,:,:,:).*(1000/binwidth);
            maxdriven = max(max(squeeze(condfr(cll,1,:,:))))-mean(bl);
            
            nrun(cll) = length(oktrials);
            nstill(cll) = length(stilltrials);
            
            cllname{cll} = files(fi).name;
            printname = files(fi).name;
            printname(find(printname=='_')) = ' ';
            animalno(cll) = animal(blck);
            recording(cll) = blck;
            depth(cll) = result.depth;
            pangle(cll) = penangle(blck);
            
            % fix special cases for old data
            for l = 1:length(unique(light))
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
            
%             % anova between largest smalles and light on off
%             avl0s1 = reshape(orifrs(1,:,1,:),40,1);
%             avl0s5 = reshape(orifrs(1,:,5,:),40,1);
%             avl1s1 = reshape(orifrs(2,:,1,:),40,1);
%             avl1s5 = reshape(orifrs(2,:,5,:),40,1);            
%             h1 = [zeros(40,1);zeros(40,1);ones(40,1);ones(40,1)]; %light
%             h2 = [zeros(40,1);ones(40,1);zeros(40,1);ones(40,1)]; %size
%             [p,table,stats] = anovan([avl0s1;avl0s5;avl1s1;avl1s5],{h1 h2},'model','full','display','off');
%             lightp(cll) = p(1); sizep(cll) = p(2); slip(cll) = p(3);
            
            % ranksum for spontaneous vs each size
            for sz = 1:length(sizes) % find out whether cell is visually driven 
                [svdp(cll,sz),cvd(sz)] = ranksum(frs(find(gratingInfo.size == 0 & light == 0)),...
                    frs(find(gratingInfo.size == sizes(sz) & light == 0)));
                if ~isnan(cvd(sz)) & cvd(sz) & mean(frs(find(gratingInfo.size == 0 & light == 0)))>mean(frs(find(gratingInfo.size == sizes(sz) & light == 0)))
                    cvd(sz) = -1;
                end
            end         
            condvismod(cll,:) = cvd;  
            
            %determine if cell is modulated by light or vis stimulus (from own baseline)
            lightmod(cll) = signrank(frs(find(light)),frs(find(light == 0)));
            vismod(cll) = signrank(frs(find(light == 0 & gratingInfo.Contrast ~= 0)),bl(find(light == 0 & gratingInfo.Contrast ~= 0)));
            if ~isnan(vismod(cll)) & vismod(cll) & mean(frs(find(light == 0 &...
                    gratingInfo.Contrast ~= 0)))<mean(bl(find(light == 0 & gratingInfo.Contrast ~= 0)))
                vismod(cll) = -1;
            end
            
            wvchan = find(var(result.waveforms) == max(var(result.waveforms)));
            spike = result.waveforms(:,wvchan);
            interpspike = spline(1:32,spike,1:.1:32);
            [adiff(cll),swidth(cll),ptr(cll),eslope(cll)] = spikequant(interpspike);
            
            waveform(cll,:) = spike;
            clustqual(cll) = result.clusterquality;
            
            bfr(cll) = mean(bl);

            [binnedlight,bta] = binit(mean(resp(light == 1,:),1),binwidth); binnedlight = binnedlight.*(1000/binwidth);
            [binnednolight,bta] = binit(mean(resp(light == 0,:),1),binwidth); binnednolight = binnednolight.*(1000/binwidth);
            
            %baselinesubtracted condfr
            blscondfr = squeeze(condfr(cll,:,:,:))-bfr(cll);

            %how much does the firing rate change for each condition
            condchange(cll,:,:) = squeeze(blscondfr(2,:,:))-squeeze(blscondfr(1,:,:));            
            
            %control firing rates (also with running)
            contindsnl = find(gratingInfo.size == 0 & light == 0);
            controlfr(cll,1) = mean(frs(contindsnl));
            controlerr(cll,1) = std(frs(contindsnl))./sqrt(length(contindsnl)); 
            
            contindsl = find(gratingInfo.size == 0 & light == 1);
            controlfr(cll,2) = mean(frs(contindsl));
            controlerr(cll,2) = std(frs(contindsl))./sqrt(length(contindsl));
            
            r0l0continds = intersect(contindsnl,stilltrials);
            r0l1continds = intersect(contindsl,stilltrials);
            r1l0continds = intersect(contindsnl,oktrials);
            r1l1continds = intersect(contindsl,oktrials);
            
            r0l0ctrlfr(cll) = mean(frs(r0l0continds));
            r0l0ctrlerr(cll) = std(frs(r0l0continds))./sqrt(length(r0l0continds));
            r0l1ctrlfr(cll) = mean(frs(r0l1continds));
            r0l1ctrlerr(cll) = std(frs(r0l1continds))./sqrt(length(r0l1continds));
            r1l0ctrlfr(cll) = mean(frs(r1l0continds));
            r1l0ctrlerr(cll) = std(frs(r1l0continds))./sqrt(length(r1l0continds));
            r1l1ctrlfr(cll) = mean(frs(r1l1continds));
            r1l1ctrlerr(cll) = std(frs(r1l1continds))./sqrt(length(r1l1continds));
            
            if ~isempty(contindsl)
                contlightmodp(cll) = ranksum(frs(contindsnl),frs(contindsl));
            else
                contlightmod(cll) = NaN;
            end

            prefsize(cll) = find(mean(condfr(cll,1,:,:),3) == max(mean(condfr(cll,1,:,:),3)),1);

            [nloriprefratio(cll), nldirprefratio(cll), nlprefori(cll), meanoril0(cll), nlosi(cll), meandirl0, nldsi(cll)] = getOSI(squeeze(condfr(cll,1,:,prefsize(cll)))',oris);
            [loriprefratio(cll), ldirprefratio(cll), lprefori(cll), meanoril1(cll), losi(cll), meandirl1, ldsi(cll)] = getOSI(squeeze(condfr(cll,2,:,prefsize(cll)))',oris);

            oneorifr = mean(reshape(blscondfr(:,:,prefsize(cll)),2,size(orifrs,2),2),3);
            prefori = find(oneorifr(1,:) == max(oneorifr(1,:)),1);
            [preffr(cll,:), prefdir] = max(condfr(cll,:,:,prefsize(cll)),[],3); prefdir = prefdir(1);
            ortho = mod(prefori+2,length(oris)/2); if ortho == 0, ortho = length(oris)/2; end
                
            sizetunel1(cll,:) = [controlfr(cll,2), squeeze(nanmean(condfr(cll,2,[prefori,prefori+(length(oris)/2)],:),3))'];
            sizetunel0(cll,:) = [controlfr(cll,1), squeeze(nanmean(condfr(cll,1,[prefori,prefori+(length(oris)/2)],:),3))'];
            sizetuneerrl1(cll,:) = [controlerr(cll,2), squeeze(nanmean(conderr(cll,2,[prefori,prefori+(length(oris)/2)],:),3))'];
            sizetuneerrl0(cll,:) = [controlerr(cll,1), squeeze(nanmean(conderr(cll,1,[prefori,prefori+(length(oris)/2)],:),3))'];
            xsizes = [0,sizes];            
            
            [binprefl0,bta] = binit(squeeze(condresp(cll,1,prefdir,prefsize(cll),:)),binwidth); binprefl0 = binprefl0.*(1000/binwidth);
            [binprefl1,bta] = binit(squeeze(condresp(cll,2,prefdir,prefsize(cll),:)),binwidth); binprefl1 = binprefl1.*(1000/binwidth);
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
            tempfreq = unique(gratingInfo.tFreq);
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
            errorbar(oris,squeeze(condfr(cll,2,:,prefsize(cll))),squeeze(conderr(cll,2,:,prefsize(cll))),'o-','color',lcol,'markersize',8,'linewidth',2)
            hold on
            errorbar(oris,squeeze(condfr(cll,1,:,prefsize(cll))),squeeze(conderr(cll,1,:,prefsize(cll))),'ko-','markersize',8,'linewidth',2)
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

            mx = max([max(max(onfield(cll,:,:))),max(max(offfield(cll,:,:)))]);
            mn = min([min(min(onfield(cll,:,:))),min(min(offfield(cll,:,:)))]);
            subplot(2,4,7)            
            imagesc(xax(cll,:),yax(cll,:),squeeze(onfield(cll,:,:)))
            hold on
            plot_orrf_absdeg(gaussfiton(cll),1,'w',2)
            if ~(isnan(mx)||isnan(mn)), caxis([mn,mx]);  end
            axis square
            axis xy
            title(['ON rsq: ' num2str(rson(cll))])            
            
            subplot(2,4,8)            
            imagesc(xax(cll,:),yax(cll,:),squeeze(offfield(cll,:,:)))
            hold on
            plot_orrf_absdeg(gaussfitoff(cll),1,'w',2)
            if ~(isnan(mx)||isnan(mn)), caxis([mn,mx]);  end
            axis square
            axis xy
            title(['OFF rsq: ' num2str(rsoff(cll))])
            
            if printyn
                figSize = [30 21];
                set(gcf,'paperunits','centimeters','papersize',figSize,'paperposition',[0 0 figSize])
                if cll<10, printi = ['0', int2str(cll)]; else printi = int2str(cll); end
                print([printpath ,  printi '__' files(fi).name '.pdf'],'-dpdf')
            end
            
            % running figure
            runlfr(cll) = mean(frs(intersect(find(light),oktrials)));
            runnlfr(cll) = mean(frs(intersect(find(~light),oktrials)));
            norunlfr(cll) = mean(frs(intersect(find(light),stilltrials)));
            norunnlfr(cll) = mean(frs(intersect(find(~light),stilltrials)));
            lfrerr = std(frs(find(light)))./sqrt(length(find(light)));
            nlfrerr = std(frs(find(~light)))./sqrt(length(find(~light)));
            runlfrerr = std(frs(intersect(find(light),oktrials)))./sqrt(length(intersect(find(light),oktrials)));
            runnlfrerr = std(frs(intersect(find(~light),oktrials)))./sqrt(length(intersect(find(~light),oktrials)));
            norunlfrerr = std(frs(intersect(find(light),stilltrials)))./sqrt(length(intersect(find(light),stilltrials)));
            norunnlfrerr = std(frs(intersect(find(~light),stilltrials)))./sqrt(length(intersect(find(~light),stilltrials)));
            
            l1r1 = frs(intersect(find(light),oktrials));
            l0r1 = frs(intersect(find(~light),oktrials));
            l1r0 = frs(intersect(find(light),stilltrials));
            l0r0 = frs(intersect(find(~light),stilltrials));
            anovavec = [l0r0;l0r1;l1r0;l1r1];
            g1 = [zeros(length(l0r0),1);zeros(length(l0r1),1);ones(length(l1r0),1);ones(length(l1r1),1)]; %light
            g2 = [zeros(length(l0r0),1);ones(length(l0r1),1);zeros(length(l1r0),1);ones(length(l1r1),1)]; %running
            [p,table,stats] = anovan(anovavec,{g1 g2},'model','full','display','off');
            lp(cll) = p(1); rp(cll) = p(2); rlip(cll) = p(3);
            
            r0omi(cll) = (norunlfr(cll)-norunnlfr(cll))/(norunlfr(cll)+norunnlfr(cll));
            r1omi(cll) = (runlfr(cll)-runnlfr(cll))/(runlfr(cll)+runnlfr(cll));
            l0rmi(cll) = (runnlfr(cll)-norunnlfr(cll))/(runnlfr(cll)+norunnlfr(cll));
            l1rmi(cll) = (runlfr(cll)-norunlfr(cll))/(runlfr(cll)+norunlfr(cll));
            
            trialfrl0(cll,:) = frs(find(light == 0),:);
            trialfrl1(cll,:) = frs(find(light == 1),:);
            trialbl(cll,:) = bl;
            
            disp([files(fi).name '   done'])
            cll = cll + 1;
            
        end
    end
    save(popfile, '-v7.3');
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

%putative halo/ChR2 expressing for SOM Halo/ChR2 pop
% phe([7,37,63,65,81,98,107,  128,133,151,158,166,167,175]) = 1;

% %putative halo expressing for SOM+PV combined pop: (added the first PV Halo 1/7/16)
phe([7,37,106,119,163,189,232]) = 1; % took out 330 which is 100 in PV 1/11/16

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

% % plot Halo cell PSTHs
% bta = bta-300;
% hcs = find(phe);
% for i = 1:length(hcs)
%     figure
%     plot(-299:2700,squeeze(nanmean(condfiltresp(hcs(i),1,:,prefsize(hcs(i)),:),3)).*1000,'k','linewidth',2)
%     hold on
%     plot(-299:2700,squeeze(nanmean(condfiltresp(hcs(i),2,:,prefsize(hcs(i)),:),3)).*1000,'r','linewidth',2)
%     
% end


for i = 1:length(depth)    
    
    [sprho(i), sprp(i)] = corr(speeds(i,:)',firingrates(i,:)','type','Spearman');
    
    onmax(i) = max(max(onfield(i,:,:)));
    offmax(i) = max(max(offfield(i,:,:)));
    onoffratio(i) = log(onmax(i)./offmax(i));  
    monconind(i) = (onmax(i)-offmax(i))/(onmax(i)+offmax(i)); % like Li Tao
    
    ampon(i) = gaussfiton(i).amp;   % amplitude of fit
    ampoff(i) = gaussfitoff(i).amp;
    amponoffratio(i) = log(ampon(i)./ampoff(i));
    offson(i) = gaussfiton(i).offs;   % offset/DC of fit
    offsoff(i) = gaussfitoff(i).offs;
    
    if gaussfiton(i).xspreadDeg>gaussfiton(i).yspreadDeg
        gaussfiton(i).longedge = gaussfiton(i).xspreadDeg;
        gaussfiton(i).shortedge = gaussfiton(i).yspreadDeg;
        onangle(i) = gaussfiton(i).thetadeg;
    else
        gaussfiton(i).longedge = gaussfiton(i).yspreadDeg;
        gaussfiton(i).shortedge = gaussfiton(i).xspreadDeg;
        onangle(i) = gaussfiton(i).thetadeg-90;
        if(onangle(i)<0), onangle(i) = onangle(i)+180; end
    end
    if gaussfitoff(i).xspreadDeg>gaussfitoff(i).yspreadDeg
        gaussfitoff(i).longedge = gaussfitoff(i).xspreadDeg;
        gaussfitoff(i).shortedge = gaussfitoff(i).yspreadDeg;
        offangle(i) = gaussfitoff(i).thetadeg;
    else
        gaussfitoff(i).longedge = gaussfitoff(i).yspreadDeg;
        gaussfitoff(i).shortedge = gaussfitoff(i).xspreadDeg;
        offangle(i) = gaussfitoff(i).thetadeg-90;
        if(offangle(i)<0), offangle(i) = offangle(i)+180; end
    end
    meanspreadon(i) = (gaussfiton(i).xspreadDeg+gaussfiton(i).yspreadDeg)/2;
    meanspreadoff(i) = (gaussfitoff(i).xspreadDeg+gaussfitoff(i).yspreadDeg)/2;
    onoffd(i) = sqrt((gaussfiton(i).xcenterDeg-gaussfitoff(i).xcenterDeg).^2+(gaussfiton(i).ycenterDeg-gaussfitoff(i).ycenterDeg).^2);
    areaon(i) = gaussfiton(i).xspreadDeg.*gaussfiton(i).yspreadDeg.*pi;
    areaoff(i) = gaussfitoff(i).xspreadDeg.*gaussfitoff(i).yspreadDeg.*pi;
    elongationon(i) = gaussfiton(i).longedge./gaussfiton(i).shortedge;
    elongationoff(i) = gaussfitoff(i).longedge./gaussfitoff(i).shortedge;
    overlap(i) = (meanspreadon(i)+meanspreadoff(i)-onoffd(i))./(meanspreadon(i)+meanspreadoff(i)+onoffd(i));
    xdist(i) = gaussfiton(i).xcenterDeg-gaussfitoff(i).xcenterDeg;
    ydist(i) = gaussfiton(i).ycenterDeg-gaussfitoff(i).ycenterDeg;
    connangle(i) = atand(ydist(i)./xdist(i));
    orthconnangle(i) = connangle(i)+90;
end

% receptive field measures
okfiton = rson>.66 & (ampon./onmax)<3;   % decent rsqared and not a crazy high fit as sometimes for these one pixel RFs
okfitoff = rsoff>.66 & (ampoff./offmax)<3; 

% ta = linspace(-300,2700,100);
% i = 3
% figure
% exampleplot(squeeze(bincondcellresp(i,:,:,:,:)),squeeze(conderr(i,:,:,:)),controlfr(i,:),controlerr(i,:),ta,xsizes,oris)

respta = linspace(-299,2700,100);
normwin = find(respta>500&respta<1500);
for i = 1:size(condfr,1)
    prefs(i) = find(mean(condfr(i,1,:,:),3) == max(mean(condfr(i,1,:,:),3)),1,'last'); %preferred size
    aps(i) = find(mean(condfr(i,1,:,:),3) == min(mean(condfr(i,1,:,:),3)),1); % anti preferred size
    prefori(i) = find(condfr(i,1,:,prefs(i)) == max(condfr(i,1,:,prefs(i))),1);    % preferred ori
    preforil1(i) = find(condfr(i,2,:,prefs(i)) == max(condfr(i,2,:,prefs(i))),1);
    preffrl0(i) = condfr(i,1,prefori(i),prefs(i));
    preffrl1(i) = condfr(i,2,prefori(i),prefs(i));
    
    for sz = 1:5
        [orir(i,sz), dirr(i,sz), x, smeanori(i,sz), osi(i,sz), smeandir(i,sz), dsi(i,sz)] = getOSI(squeeze(condfr(i,1,:,sz))',oris);
        [orir1(i,sz), dirr1(i,sz), x, smeanori1(i,sz), osi1(i,sz), smeandir1(i,sz), dsi1(i,sz)] = getOSI(squeeze(condfr(i,2,:,sz))',oris);     
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
    
    r0l0omst(i,:) = [r0l0ctrlfr(i), squeeze(nanmean(stillcondfr(i,1,:,:),3))'];
    r0l1omst(i,:) = [r0l1ctrlfr(i), squeeze(nanmean(stillcondfr(i,2,:,:),3))'];
    r1l0omst(i,:) = [r1l0ctrlfr(i), squeeze(nanmean(runcondfr(i,1,:,:),3))'];
    r1l1omst(i,:) = [r1l1ctrlfr(i), squeeze(nanmean(runcondfr(i,2,:,:),3))'];
    
    nr0l0omst(i,:) = r0l0omst(i,:)./max(r0l0omst(i,:));  %normed to no light condition each
    nr0l1omst(i,:) = r0l1omst(i,:)./max(r0l0omst(i,:));
    nr1l0omst(i,:) = r1l0omst(i,:)./max(r1l0omst(i,:));
    nr1l1omst(i,:) = r1l1omst(i,:)./max(r1l0omst(i,:));
    
    ntr_r0l0omst(i,:) = r0l0omst(i,:)./max(r1l0omst(i,:)); % all normed to running no light
    ntr_r0l1omst(i,:) = r0l1omst(i,:)./max(r1l0omst(i,:));
    ntr_r1l0omst(i,:) = r1l0omst(i,:)./max(r1l0omst(i,:));
    ntr_r1l1omst(i,:) = r1l1omst(i,:)./max(r1l0omst(i,:));
    
    nts_r0l0omst(i,:) = r0l0omst(i,:)./max(r0l0omst(i,:)); % all normed to still no light
    nts_r0l1omst(i,:) = r0l1omst(i,:)./max(r0l0omst(i,:));
    nts_r1l0omst(i,:) = r1l0omst(i,:)./max(r0l0omst(i,:));
    nts_r1l1omst(i,:) = r1l1omst(i,:)./max(r0l0omst(i,:));

    
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
    sizetunel0r0(i,:) = squeeze(nanmean(stillcondfr(i,1,:,:),3));
    sizetunel1r0(i,:) = squeeze(nanmean(stillcondfr(i,2,:,:),3));
    sizetunel0r1(i,:) = squeeze(nanmean(runcondfr(i,1,:,:),3));
    sizetunel1r1(i,:) = squeeze(nanmean(runcondfr(i,2,:,:),3));
    nsizetunel0r0(i,:)= sizetunel0r0(i,:)./max(sizetunel0r1(i,:));
    nsizetunel0r1(i,:)= sizetunel0r1(i,:)./max(sizetunel0r1(i,:));
    nsizetunel1r0(i,:)= sizetunel1r0(i,:)./max(sizetunel0r1(i,:));
    nsizetunel1r1(i,:)= sizetunel1r1(i,:)./max(sizetunel0r1(i,:));
%     if ~isnan(max(mean(stillcondfr(i,1,:,:),3)))
%         psr0(i) = find(mean(stillcondfr(i,1,:,:),3) == max(mean(stillcondfr(i,1,:,:),3)),1,'last');
%     else
%         psr0(i) = NaN;
%     end
% %     if ~isnan(
%     preforir0(i) = find(stillcondfr(i,1,:,ps(i)) == max(stillcondfr(i,1,:,ps(i))),1);
%     psr1(i) = find(mean(runcondfr(i,1,:,:),3) == max(mean(runcondfr(i,1,:,:),3)),1,'last');
%     preforir1(i) = find(runcondfr(i,1,:,ps(i)) == max(runcondfr(i,1,:,ps(i))),1);
     if prefs(i) == 1|2
         minds = [1,2,3];
     elseif prefs(i) == 3
         minds = [2,3,4];
     else
         minds = [3,4,5];
     end
     oritunel0r0(i,:) = squeeze(nanmean(stillcondfr(i,1,:,minds),4));
     oritunel1r0(i,:) = squeeze(nanmean(stillcondfr(i,2,:,minds),4));
     oritunel0r1(i,:) = squeeze(nanmean(runcondfr(i,1,:,minds),4));
     oritunel1r1(i,:) = squeeze(nanmean(runcondfr(i,2,:,minds),4));
     
     oritunel0(i,:) = squeeze(nanmean(condfr(i,1,:,minds),4));
     oritunel1(i,:) = squeeze(nanmean(condfr(i,2,:,minds),4));
    
end
ntr_r0l0omst(isinf(ntr_r0l0omst)) = NaN;
ntr_r0l1omst(isinf(ntr_r0l1omst)) = NaN;
ntr_r1l0omst(isinf(ntr_r1l0omst)) = NaN;
ntr_r1l1omst(isinf(ntr_r1l1omst)) = NaN;
nts_r0l0omst(isinf(nts_r0l0omst)) = NaN;
nts_r0l1omst(isinf(nts_r0l1omst)) = NaN;
nts_r1l0omst(isinf(nts_r1l0omst)) = NaN;
nts_r1l1omst(isinf(nts_r1l1omst)) = NaN;


sizes = xsizes(2:end)
bars = [mean(sizes(prefs(l23rs))),mean(sizes(prefs(l23fs)));...
    mean(sizes(prefs(l4rs))),mean(sizes(prefs(l4fs)));...
    mean(sizes(prefs(l5rs))),mean(sizes(prefs(l5fs)))];
errorbars = [std(sizes(prefs(l23rs)))./sqrt(sum(l23rs)),...
    std(sizes(prefs(l23fs)))./sqrt(sum(l23fs));...
    std(sizes(prefs(l4rs)))./sqrt(sum(l4rs)),...
    std(sizes(prefs(l4fs)))./sqrt(sum(l4fs));...
    std(sizes(prefs(l5rs)))./sqrt(sum(l5rs)),...
    std(sizes(prefs(l5fs)))./sqrt(sum(l5fs))];
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


% delta fr per rank plots
cond = l23rs;
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
