function all_populationana
clear all;
% 
% % SOM Halo population
% animalids =  {'150331', '150401','150527','150529','150602','150603','150625','150825','150831','150902','150909', '150915', '151023', '151027','151109', '151110', '151209', '160122', '160718', '160721', '160725', '160726', '160728', '160729', '160801', '160802', '160804_2', '160920', '160921'};
% animal    =  [1,         2,       3,       4,       5,       6,       7,       8,       9,       10,      11,       12,       13,       14,      15,       16,       17,       18,       19,       20,       21,       22,       23,       24,       25,       26,       27,         28,       29];
% rfblocks   = [1,         1,       1,       1,       1,       2,       4,       3,       3,       1,       10,       8,        9,        6,       6,        9,        7,        4,        1,        1,        1,        2,        1,        1,        2,        4,        2,          7,        6];
% contblocks = [2,         4,       7,       5,       3,       4,       7,       7,       6,       8,       9,        6,        15,       NaN,     NaN,      NaN,      11,       NaN,      NaN,      4,        4,        6,        7,        3,        4,        5,        6,          4,        4];
% sizeblocks = [3,         5,       11,      4,       5,       3,       6,       5,       4,       3,       4,        3,        14,       3,       11,       13,       6,        3,        NaN,      4,        4,        6,        7,        3,        4,        5,        6,          4,        4];
% movieblock = [NaN,       NaN,     4,       7,       7,       6,       5,       11,      7,       12,      8,        12,       11,       NaN,     12,       15,       13,       NaN,      3,        3,        3,        3,        3,        2,        3,        3,        3,          3,        3];
% surrblock  = [NaN,       NaN,     NaN,     NaN,     NaN,     NaN,     9,       8,       8,       10,      6,        9,        5,        5,       9,        11,       NaN,      NaN,      4,        NaN,      NaN,      4,        NaN,      6,        7,        7,        5,          NaN,      7];
% phasblock =  [NaN,       NaN,     NaN,     NaN,     NaN,     NaN,     NaN,     9,       NaN,     11,      7,        11,       6,        8,       7,        10,       NaN,      NaN,      NaN,      NaN,      NaN,      NaN,      NaN,      5,        8,        8,        4,          6,        NaN];                                                     
% electrodes = [[1,32];   [1,32];  [1,32];  [1,32];  [1,32];  [1,16];  [1,16];  [17,32]; [1,16];  [1,16];  [1,16];   [17,32];  [17,32];  [17,32]; [17,32];  [1,16];   [17,32];  [1,16];   [17,32];  [1,16];   [1,16];   [1,16];   [1,16];   [1,16];   [1,16];   [1,16];   [1,16];     [1,16];   [1,16]];
% penangle =   [25,        25,      25,      25,      25,      10,      10,      25,      25,      25,      25,       25,       25,       25,       25,       25,       25,       25,      25,       25,       25,       25,       25,       25,       25,       25,       25,         25,       25];
% printpath = 'C:\Users\Julia\work\data\populations\SOM_Halo_later\all\units\';
% rfprintpath = 'C:\Users\Julia\work\data\populations\SOM_Halo_later\all\RFs\';
% resultpath = 'C:\Users\Julia\work\data\populations\SOM_Halo_later\all\resultfiles\';
% popfile = 'C:\Users\Julia\work\data\populations\SOM_Halo_later\all\all_population.mat';

% % Pv Halo population
% animalids =  {'150629', '150630', '150731', '150804', '150818', '150819', '150820', '150823', '150824', '151211', '160114', '160115', '160204', '160205', '160210', '160217', '160328', '160901', '160902', '160902_2', '160906', '160907', '160922', '160923'};
% animal    =  [ 1,        2,        3,        4,        5,        6,        7,        8,        9,        10,       11,       12,       13,       14,       15,       16,       17,       18,       19,       20,         21,       22,       23,       24];
% rfblocks   = [ 10,       9,        2,        2,        2,        3,        8,        3,        2,        7,        9,        12,       6,        9,        5,        6,        1,        2,        2,        3,          3,        3,        3,        1];
% contblocks = [ 3,        NaN,      12,       11,       NaN,      NaN,      NaN,      10,       8,        NaN,      NaN,      NaN,      NaN,      NaN,      NaN,      10,       NaN,      3,        4,        5,          7,        5,        5,        3];
% sizeblocks = [ 4,        NaN,      3,        4,        3,        7,        3,        5,        3,        3,        4,        10,       3,        4,        9,        5,        2,        3,        4,        5,          7,        5,        5,        3];
% movieblock = [ 9,        10,       7,        5,        6,        11,       9,        8,        5,        9,        13,       13,       10,       10,       NaN,      8,        5,        4,        3,        4,          5,        4,        4,        2];
% surrblock  = [ 8,        6,        5,        7,        4,        NaN,      10,       6,        6,        NaN,      14,       NaN,      9,        NaN,      8,        NaN,      NaN,      NaN,      NaN,      9,          8,        8,        NaN,      6];
% electrodes = [[1,16];   [1,16];   [1,16];   [1,16];   [1,16];   [1,16];   [1,16];   [1,16];   [17,32];  [1,16];   [1,16];   [1,16];   [1,16];   [17,32];  [17,32];  [1,16];   [1,16];   [1,16];   [1,16];   [1,16];     [1,16];   [1,16];   [1,16];   [1,16]];
% penangle =   [10,        10,       25,       25,       25,       25,       25,       25,       25,       25,       25,       25,       25,       25,       25,       25,       25,       25,       25,       25,         25,       25,       25,       25];
% printpath = 'C:\Users\Julia\work\data\populations\PV_Halo\all\units\';
% rfprintpath = 'C:\Users\Julia\work\data\populations\PV_Halo\all\RFs\';
% resultpath = 'C:\Users\Julia\work\data\populations\PV_Halo\all\resultfiles\';
% popfile = 'C:\Users\Julia\work\data\populations\PV_Halo\all\all_population.mat';

% VIP Halo population  -- -think about whether to include 180404 because of bad sorts
animalids =  {'171201', '171205', '180125', '180404', '180410', '180412', '180413'};
animal    =  [ 1,        1,        2,        3,        4,        5,        5];
rfblocks   = [ 3,        5,        1,        5,        4,        9,        6];
contblocks = [ 2,        7,        4,        9,        6,        6,        10];
sizeblocks = [ 2,        7,        4,        9,        6,        6,        10];
movieblock = [ 5,        6,        3,        3,        10,       11,       4];
surrblock  = [ NaN,      NaN,      7,        7,        5,        10,       8];
phasblock =  [ NaN,      NaN,      NaN,      NaN,      NaN,      NaN,      NaN];
electrodes = [[1,16];   [1,16];   [1,16];   [1,16];   [1,16];   [1,16];   [1,16]];
penangle =   [25,        25,       25,       25,       25,       25,       25];
printpath = 'C:\Users\Julia\work\data\populations\VIP_Halo\all\units\';
rfprintpath = 'C:\Users\Julia\work\data\populations\VIP_Halo\all\RFs\';
coszprintpath = 'C:\Users\Julia\work\data\populations\VIP_Halo\all\cosz\';
resultpath = 'C:\Users\Julia\work\data\populations\VIP_Halo\all\resultfiles\';
popfile = 'C:\Users\Julia\work\data\populations\VIP_Halo\all\all_population.mat';

% % control population
% animalids =  {'151222', '160112', '160113', '160711'};
% animal    =  [ 1,        2,        3,        4];
% rfblocks   = [ 10,       2,        8,        2];
% contblocks = [ NaN,      NaN,      NaN,      5];
% sizeblocks = [ NaN,      3,        6,        NaN];
% movieblock = [ 12,       NaN,      9,        3];
% surrblock  = [ NaN,      9,        11,       NaN];
% electrodes = [[1,16];   [1,16];   [1,16];   [1,16]];
% penangle =   [25,        25,       25,       25];
% printpath = 'C:\Users\Julia\work\data\populations\control\all\units\';
% rfprintpath = 'C:\Users\Julia\work\data\populations\control\all\RFs\';
% resultpath = 'C:\Users\Julia\work\data\populations\control\all\resultfiles\';
% popfile = 'C:\Users\Julia\work\data\populations\control\all\all_population.mat';
% 
% % SOM PV combined population: 1:27 SOM 28:51 PV
% animalids =  {'150331', '150401','150527','150529','150602','150603','150625','150825','150831','150902','150909', '150915', '151023', '151027','151109', '151110', '151209', '160122', '160718', '160721', '160725', '160726', '160728', '160729', '160801', '160802', '160804_2', '150629', '150630', '150731', '150804', '150818', '150819', '150820', '150823', '150824', '151211', '160114', '160115', '160204', '160205', '160210', '160217', '160328', '160901', '160902', '160902_2', '160906', '160907', '160922', '160923'};
% animal    =  [1,         2,       3,       4,       5,       6,       7,       8,       9,       10,      11,       12,       13,       14,      15,       16,       17,       18,       19,       20,       21,       22,       23,       24,       25,       26,       27,         28,       29,       30,       31,       32,       33,       34,       35,       36,       37,       38,       39,       40,       41,       42,       43,       44,       45,       46,       47,         48,       49,       50,       51];
% rfblocks   = [1,         1,       1,       1,       1,       2,       4,       3,       3,       1,       10,       8,        9,        6,       6,        9,        7,        4,        1,        1,        1,        2,        1,        1,        2,        4,        2,          10,       9,        2,        2,        2,        3,        8,        3,        2,        7,        9,        12,       6,        9,        5,        6,        1,        2,        2,        3,          3,        3,        3,        1];
% contblocks = [2,         4,       7,       5,       3,       4,       7,       7,       6,       8,       9,        6,        15,       NaN,     NaN,      NaN,      11,       NaN,      NaN,      4,        4,        6,        7,        3,        4,        5,        6,          3,        NaN,      12,       11,       NaN,      NaN,      NaN,      10,       8,        NaN,      NaN,      NaN,      NaN,      NaN,      NaN,      10,       NaN,      3,        4,        5,          7,        5,        5,        3];
% sizeblocks = [3,         5,       11,      4,       5,       3,       6,       5,       4,       3,       4,        3,        14,       3,       11,       13,       6,        3,        NaN,      4,        4,        6,        7,        3,        4,        5,        6,          4,        NaN,      3,        4,        3,        7,        3,        5,        3,        3,        4,        10,       3,        4,        9,        5,        2,        3,        4,        5,          7,        5,        5,        3];
% movieblock = [NaN,       NaN,     4,       7,       7,       6,       5,       11,      7,       12,      8,        12,       11,       NaN,     12,       15,       13,       NaN,      3,        3,        3,        3,        3,        2,        3,        3,        3,          9,        10,       7,        5,        6,        11,       9,        8,        5,        9,        13,       13,       10,       10,       NaN,      8,        5,        4,        3,        4,          5,        4,        4,        2];
% surrblock  = [NaN,       NaN,     NaN,     NaN,     NaN,     NaN,     9,       8,       8,       10,      6,        9,        5,        5,       9,        11,       NaN,      NaN,      4,        NaN,      NaN,      4,        9,        6,        7,        7,        5,          8,        6,        5,        7,        4,        NaN,      10,       6,        6,        NaN,      14,       NaN,      9,        NaN,      8,        NaN,      NaN,      NaN,      NaN,      9,          8,        8,        NaN,      6];
% electrodes = [[1,32];   [1,32];  [1,32];  [1,32];  [1,32];  [1,16];  [1,16];  [17,32]; [1,16];  [1,16];  [1,16];   [17,32];  [17,32];  [17,32]; [17,32];  [1,16];   [17,32];  [1,16];   [17,32];  [1,16];   [1,16];   [1,16];   [1,16];   [1,16];   [1,16];   [1,16];   [1,16];     [1,16];   [1,16];   [1,16];   [1,16];   [1,16];   [1,16];   [1,16];   [1,16];   [17,32];  [1,16];   [1,16];   [1,16];   [1,16];   [17,32];  [17,32];  [1,16];   [1,16];   [1,16];   [1,16];   [1,16];     [1,16];   [1,16];   [1,16];   [1,16]];
% penangle =   [25,        25,      25,      25,      25,      10,      10,      25,      25,      25,      25,       25,       25,       25,       25,       25,       25,       25,      25,       25,       25,       25,       25,       25,       25,       25,       25,         10,       10,       25,       25,       25,       25,       25,       25,       25,       25,       25,       25,       25,       25,       25,       25,       25,       25,       25,       25,         25,       25,       25,       25];
% printpath = 'C:\Users\Julia\work\data\populations\SOMPVcombined\all\units\';
% rfprintpath = 'C:\Users\Julia\work\data\populations\SOMPVcombined\all\RFs\';
% resultpath = 'C:\Users\Julia\work\data\populations\SOMPVcombined\all\resultfiles\';
% popfile = 'C:\Users\Julia\work\data\populations\SOMPVcombined\all\all_population.mat';
% 
% % SOM PV VIP combined population: 1:27 SOM 28:51 PV 52:
% animalids =  {'150331', '150401','150527','150529','150602','150603','150625','150825','150831','150902','150909', '150915', '151023', '151027','151109', '151110', '151209', '160122', '160718', '160721', '160725', '160726', '160728', '160729', '160801', '160802', '160804_2', '150629', '150630', '150731', '150804', '150818', '150819', '150820', '150823', '150824', '151211', '160114', '160115', '160204', '160205', '160210', '160217', '160328', '160901', '160902', '160902_2', '160906', '160907', '160922', '160923', '171201', '171205', '180125', '180404', '180410', '180412', '180413'};
% animal    =  [1,         2,       3,       4,       5,       6,       7,       8,       9,       10,      11,       12,       13,       14,      15,       16,       17,       18,       19,       20,       21,       22,       23,       24,       25,       26,       27,         28,       29,       30,       31,       32,       33,       34,       35,       36,       37,       38,       39,       40,       41,       42,       43,       44,       45,       46,       47,         48,       49,       50,       51,       52,       52,       53,       54,       55,       56,       56];
% rfblocks   = [1,         1,       1,       1,       1,       2,       4,       3,       3,       1,       10,       8,        9,        6,       6,        9,        7,        4,        1,        1,        1,        2,        1,        1,        2,        4,        2,          10,       9,        2,        2,        2,        3,        8,        3,        2,        7,        9,        12,       6,        9,        5,        6,        1,        2,        2,        3,          3,        3,        3,        1,        3,        5,        1,        5,        4,        9,        6];
% contblocks = [2,         4,       7,       5,       3,       4,       7,       7,       6,       8,       9,        6,        15,       NaN,     NaN,      NaN,      11,       NaN,      NaN,      4,        4,        6,        7,        3,        4,        5,        6,          3,        NaN,      12,       11,       NaN,      NaN,      NaN,      10,       8,        NaN,      NaN,      NaN,      NaN,      NaN,      NaN,      10,       NaN,      3,        4,        5,          7,        5,        5,        3,        2,        7,        4,        9,        6,        6,        10];
% sizeblocks = [3,         5,       11,      4,       5,       3,       6,       5,       4,       3,       4,        3,        14,       3,       11,       13,       6,        3,        NaN,      4,        4,        6,        7,        3,        4,        5,        6,          4,        NaN,      3,        4,        3,        7,        3,        5,        3,        3,        4,        10,       3,        4,        9,        5,        2,        3,        4,        5,          7,        5,        5,        3,        2,        7,        4,        9,        6,        6,        10];
% movieblock = [NaN,       NaN,     4,       7,       7,       6,       5,       11,      7,       12,      8,        12,       11,       NaN,     12,       15,       13,       NaN,      3,        3,        3,        3,        3,        2,        3,        3,        3,          9,        10,       7,        5,        6,        11,       9,        8,        5,        9,        13,       13,       10,       10,       NaN,      8,        5,        4,        3,        4,          5,        4,        4,        2,        5,        6,        3,        3,        10,       11,       4];
% surrblock  = [NaN,       NaN,     NaN,     NaN,     NaN,     NaN,     9,       8,       8,       10,      6,        9,        5,        5,       9,        11,       NaN,      NaN,      4,        NaN,      NaN,      4,        NaN,        6,        7,        7,        5,          8,        6,        5,        7,        4,        NaN,      10,       6,        6,        NaN,      14,       NaN,      9,        NaN,      8,        NaN,      NaN,      NaN,      NaN,      9,          8,        8,        NaN,      6,        NaN,      NaN,      7,        7,        5,        10,       8];
% electrodes = [[1,32];   [1,32];  [1,32];  [1,32];  [1,32];  [1,16];  [1,16];  [17,32]; [1,16];  [1,16];  [1,16];   [17,32];  [17,32];  [17,32]; [17,32];  [1,16];   [17,32];  [1,16];   [17,32];  [1,16];   [1,16];   [1,16];   [1,16];   [1,16];   [1,16];   [1,16];   [1,16];     [1,16];   [1,16];   [1,16];   [1,16];   [1,16];   [1,16];   [1,16];   [1,16];   [17,32];  [1,16];   [1,16];   [1,16];   [1,16];   [17,32];  [17,32];  [1,16];   [1,16];   [1,16];   [1,16];   [1,16];     [1,16];   [1,16];   [1,16];   [1,16];   [1,16];   [1,16];   [1,16];   [1,16];   [1,16];   [1,16];   [1,16]];
% penangle =   [25,        25,      25,      25,      25,      10,      10,      25,      25,      25,      25,       25,       25,       25,       25,       25,       25,       25,      25,       25,       25,       25,       25,       25,       25,       25,       25,         10,       10,       25,       25,       25,       25,       25,       25,       25,       25,       25,       25,       25,       25,       25,       25,       25,       25,       25,       25,         25,       25,       25,       25,       25,       25,       25,       25,       25,       25,       25];
% printpath = 'C:\Users\Julia\work\data\populations\SOMPVVIPcombined\all\units\';
% rfprintpath = 'C:\Users\Julia\work\data\populations\SOMPVVIPcombined\all\RFs\';
% resultpath = 'C:\Users\Julia\work\data\populations\SOMPVVIPcombined\all\resultfiles\';
% popfile = 'C:\Users\Julia\work\data\populations\SOMPVVIPcombined\all\all_population.mat';


tic
printyn = 1;

recalculate = 0;
recalculate_rf = 0;
recalculate_co = 0;
recalculate_sz = 0;
recalculate_mv = 0;
recalculate_sr = 0;
recalculate_ph = 0;

if ~exist(popfile) || recalculate

    rfres = []; szres = []; cores = []; mvres = []; srres = []; phres = []; pangle = []; animalno = [];
    
    for blck = 1:length(rfblocks)

        supath = ['C:\Users\Julia\work\data\' animalids{blck} '\singleunits\'];
        rfbasename = [animalids{blck} '_block' int2str(rfblocks(blck))];
        contbasename = [animalids{blck} '_block' int2str(contblocks(blck))];
        sizebasename = [animalids{blck} '_block' int2str(sizeblocks(blck))];
        moviebasename = [animalids{blck} '_block' int2str(movieblock(blck))];
        surrbasename = [animalids{blck} '_block' int2str(surrblock(blck))];
        phasbasename = [animalids{blck} '_block' int2str(phasblock(blck))];
              
        clear rf_result; clear pa; clear anno;
        ci = 1;
        
        rfresultfile = [resultpath, rfbasename, '_RFresults.mat'];
        if ~exist(rfresultfile) || recalculate_rf
            if ~isnan(rfblocks(blck))
                files = dir([supath, [rfbasename '_'], '*.mat']);
                for fi = 1:length(files)
                    i = strfind(files(fi).name, 'tet');
                    if strcmp(files(fi).name(i+4),'_')
                        tetno = strread(files(fi).name(i+3)); % single character number
                    else
                        tetno = strread(files(fi).name(i+3:i+4)); % number >10
                    end
                    if tetno*4<electrodes(blck,1) || tetno*4>electrodes(blck,2) % assure we're only getting V1 in the population
                        continue;
                    end
                    if strfind(files(fi).name, 'MU'), continue; end
                    
                    rf_result(ci) = get_rfresult([supath, files(fi).name]);
                    pa(ci) = penangle(blck); anno(ci) = animal(blck);
                    ci = ci+1;
                end
            end
            save(rfresultfile, 'rf_result', 'pa', 'anno');
        else
            load(rfresultfile);
        end
        rfres = [rfres, rf_result];
        pangle = [pangle, pa]; animalno = [animalno, anno];
        disp(['done with RF: ' rfbasename])
        
        clear co_result;
        ci = 1;
        coresultfile = [resultpath, contbasename, '_COresults.mat'];
        if ~exist(coresultfile) || recalculate_co
            if ~isnan(contblocks(blck))
                files = dir([supath, [contbasename '_'], '*.mat']);
                for fi = 1:length(files)                
                    i = strfind(files(fi).name, 'tet');
                    if strcmp(files(fi).name(i+4),'_')
                        tetno = strread(files(fi).name(i+3)); % single character number
                    else
                        tetno = strread(files(fi).name(i+3:i+4)); % number >10
                    end
                    if tetno*4<electrodes(blck,1) || tetno*4>electrodes(blck,2) % assure we're only getting V1 in the population
                        continue;
                    end
                    if strfind(files(fi).name, 'MU'), continue; end

                    co_result(ci) = get_contrastresult([supath, files(fi).name]);  
                    ci = ci+1;
                end
            else
                co_result = makeNaNco(length(rf_result),animalids(blck));
            end                 
            save(coresultfile, 'co_result');
        else
            load(coresultfile);
        end
        cores = [cores, co_result];
        disp(['done with contrast: ' contbasename])
        
        clear sz_result;
        ci = 1;
        szresultfile = [resultpath, sizebasename, '_SZresults.mat'];
        if ~exist(szresultfile) || recalculate_sz
            if ~isnan(sizeblocks(blck))
                files = dir([supath, [sizebasename '_'], '*.mat']);
                for fi = 1:length(files)                                
                    i = strfind(files(fi).name, 'tet');
                    if strcmp(files(fi).name(i+4),'_')
                        tetno = strread(files(fi).name(i+3)); % single character number
                    else
                        tetno = strread(files(fi).name(i+3:i+4)); % number >10
                    end
                    if tetno*4<electrodes(blck,1) || tetno*4>electrodes(blck,2) % assure we're only getting V1 in the population
                        continue;
                    end
                    if strfind(files(fi).name, 'MU'), continue; end

                    sz_result(ci) = get_sizeresult([supath, files(fi).name]);
                    ci = ci+1;
                end
            else
                sz_result = makeNaNsz(length(rf_result),animalids(blck));
            end                 
            save(szresultfile, 'sz_result');
        else
            load(szresultfile);
        end     
        szres = [szres, sz_result];  
        disp(['done with size: ' sizebasename])      
        
        clear mv_result;
        ci = 1; % cell index
        mvresultfile = [resultpath, moviebasename, '_MVresults.mat'];
        if ~exist(mvresultfile) || recalculate_mv
            if ~isnan(movieblock(blck))
                files = dir([supath, [moviebasename '_'], '*.mat']);
                for fi = 1:length(files)

                    i = strfind(files(fi).name, 'tet');
                    if strcmp(files(fi).name(i+4),'_')
                        tetno = strread(files(fi).name(i+3)); % single character number
                    else
                        tetno = strread(files(fi).name(i+3:i+4)); % number >10
                    end
                    if tetno*4<electrodes(blck,1) || tetno*4>electrodes(blck,2) % assure we're only getting V1 in the population
                        continue;
                    end
                    if strfind(files(fi).name, 'MU'), continue; end
                    
                    mv_result(ci) = get_movieresult([supath, files(fi).name]);
                    ci = ci+1;
                end
            else
                mv_result = makeNaNmv(length(rf_result),animalids(blck));
            end
            save(mvresultfile, 'mv_result');
        else
            load(mvresultfile);
        end
        mvres = [mvres, mv_result];
        disp(['done with movie: ' moviebasename])        
        
        clear sr_result;
        ci = 1; % cell index
        srresultfile = [resultpath, surrbasename, '_SRresults.mat'];
        if ~exist(srresultfile) || recalculate_sr
            if ~isnan(surrblock(blck))
                files = dir([supath, [surrbasename '_'], '*.mat']);
                for fi = 1:length(files)

                    i = strfind(files(fi).name, 'tet');
                    if strcmp(files(fi).name(i+4),'_')
                        tetno = strread(files(fi).name(i+3)); % single character number
                    else
                        tetno = strread(files(fi).name(i+3:i+4)); % number >10
                    end
                    if tetno*4<electrodes(blck,1) || tetno*4>electrodes(blck,2) % assure we're only getting V1 in the population
                        continue;
                    end
                    if strfind(files(fi).name, 'MU'), continue; end
                    
                    sr_result(ci) = get_surrresult([supath, files(fi).name]);
                    ci = ci+1;
                end
            else
                sr_result = makeNaNsr(length(rf_result),animalids(blck));
            end
            save(srresultfile, 'sr_result');
        else
            load(srresultfile);
        end
        srres = [srres, sr_result];
        disp(['done with ori-surround: ' surrbasename])
        
        clear ph_result;
        ci = 1; % cell index
        phresultfile = [resultpath, surrbasename, '_PHresults.mat'];
        if ~exist(phresultfile) || recalculate_ph
            if ~isnan(phasblock(blck))
                files = dir([supath, [phasbasename '_'], '*.mat']);
                for fi = 1:length(files)

                    i = strfind(files(fi).name, 'tet');
                    if strcmp(files(fi).name(i+4),'_')
                        tetno = strread(files(fi).name(i+3)); % single character number
                    else
                        tetno = strread(files(fi).name(i+3:i+4)); % number >10
                    end
                    if tetno*4<electrodes(blck,1) || tetno*4>electrodes(blck,2) % assure we're only getting V1 in the population
                        continue;
                    end
                    if strfind(files(fi).name, 'MU'), continue; end
                    
                    ph_result(ci) = get_phasresult([supath, files(fi).name]);
                    ci = ci+1;
                end
            else
                ph_result = makeNaNph(length(rf_result),animalids(blck));
            end
            save(phresultfile, 'ph_result');
        else
            load(phresultfile);
        end
        phres = [phres, ph_result];
        disp(['done with ori-surround: ' surrbasename])
        
        disp(['done with ', animalids{blck}]);
    end    
    save(popfile, '-v7.3');
else
    load(popfile);   
end

%% extract info from structs
for i = 1:length(mvres)
    cllname{i} = rfres(i).cellname;
    if ~isnan(szres(i).depth)
        depth(i) = szres(i).depth;
        adiff(i) = szres(i).adiff;
        swidth(i) = szres(i).swidth;
        spikeshape(i,:) = szres(i).interpspike;
    elseif ~isnan(cores(i).depth)
        depth(i) = cores(i).depth;
        adiff(i) = cores(i).adiff;
        swidth(i) = cores(i).swidth;
        spikeshape(i,:) = cores(i).interpspike;
    else
        depth(i) = mvres(i).depth;
        adiff(i) = mvres(i).adiff;
        swidth(i) = mvres(i).swidth;
        spikeshape(i,:) = mvres(i).interpspike;
    end
    
    szvismod(i) = szres(i).vismod;
    szlightmod(i) = szres(i).lightmod;
    szlfr(i) = szres(i).lfr;
    sznlfr(i) = szres(i).nlfr;
    szfr(i,:,:) = nanmean(szres(i).condfr,2);
    szerr(i,:,:) = nanmean(szres(i).conderr,2);
    szfrr1(i,:,:) = nanmean(szres(i).runcondfr,2);
    szfrr0(i,:,:) = nanmean(szres(i).stillcondfr,2);
    szr1n(i,:,:) = nanmean(szres(i).runn,2);
    szr0n(i,:,:) = nanmean(szres(i).stilln,2);
    for j = 1:5
        [vr,rps] = max(szres(i).runcondfr(1,:,j),[],2);
        [vs,sps] = max(szres(i).stillcondfr(1,:,j),[],2);
        szfrr1p(i,1,j) = vr; szfrr1p(i,2,j) = szres(i).runcondfr(2,rps,j);
        szfrr0p(i,1,j) = vs; szfrr0p(i,2,j) = szres(i).stillcondfr(2,sps,j);
    end
    szcondrespr1(i,:,:,:) = squeeze(nanmean(szres(i).runcondresp,2));    
    szomlfpspect(i,:,:,:) = szres(i).omlfpspect;
    szomlfpspecterr(i,:,:,:) = szres(i).omlfpspecterr;
    szctrlfr(i,:) = szres(i).controlfr;
    szctrlerr(i,:) = szres(i).controlerr;
    szr0ctrlfr(i,:) = szres(i).r0controlfr;
    szr0ctrlerr(i,:) = szres(i).r0controlerr;
    szr0ctrln(i,:) = szres(i).r0controln;
    szr1ctrlfr(i,:) = szres(i).r1controlfr;
    szr1ctrlerr(i,:) = szres(i).r1controlerr;
    szr1ctrln(i,:) = szres(i).r1controln;
    szbincondresp(i,:,:,:) = nanmean(szres(i).bincondresp,2);
    szrely(i,:,:) = nanmean(szres(i).condrely,2);
    szrelyn(i,:,:) = nanmean(szres(i).condrelyn,2);
    szsparseness(i,:,:) = nanmean(szres(i).condsparseness,2);
    szsparsenesswin(i,:,:) = nanmean(szres(i).condsparsenesswin,2);
    szff(i,:,:) = nanmean(szres(i).ff,2);
    szn(i,:,:) = nanmean(szres(i).condn,2);
    szorifit(i,:,:,:) = szres(i).gaussparams; %[C Rp Rn theta, sigma] [offset, amp1, amp2, center-ori, bandwidth]
    szorifitrs(i,:,:) = szres(i).gaussrsquared;
    if ~isnan(szres(i).prefsize)
        szprefsfit(i,:,:) = szres(i).gaussparams(:,szres(i).prefsize,:);
        szprefsfitrs(i,:) = szres(i).gaussrsquared(:,szres(i).prefsize);
    else
        szprefsfit(i,:,:) = nan(2,5);
        szprefsfitrs(i,:) = nan(2,1);
    end
    szprefsize(i) = szres(i).prefsize;
    szosi(i,:,:) = szres(i).condosi;
    szdsi(i,:,:) = szres(i).conddsi;
    szop(i) = szres(i).anova_op;
    szsp(i) = szres(i).anova_sp;
    szsoip(i) = szres(i).anova_soip;
    
    for l = 1:2
        for sz = 1:5
            
            [a, b, c, d, szr1osi(i,l,sz), e, szr1dsi(i,l,sz)] = getOSI(squeeze(szres(i).runcondfr(l,:,sz)),szres(i).oris);
            [a, b, c, d, szr0osi(i,l,sz), e, szr0dsi(i,l,sz)] = getOSI(squeeze(szres(i).stillcondfr(l,:,sz)),szres(i).oris);
            shiftfrr1 = szres(i).runcondfr(l,:,sz)-min(szres(i).runcondfr(1,:,sz),[],2);
            shiftfrr0 = szres(i).stillcondfr(l,:,sz)-min(szres(i).stillcondfr(1,:,sz),[],2);
            [a, b, c, d, szshr1osi(i,l,sz), e, szshr1dsi(i,l,sz)] = getOSI(shiftfrr1,szres(i).oris);
            [a, b, c, d, szshr0osi(i,l,sz), e, szshr0dsi(i,l,sz)] = getOSI(shiftfrr0,szres(i).oris);
        end
    end
    
    covismod(i) = cores(i).vismod;
    colightmod(i) = cores(i).lightmod;
    colfr(i) = cores(i).lfr;
    conlfr(i) = cores(i).nlfr;
    if size(cores(i).condfr,3) == 1
        cofr(i,:,:) = squeeze(nanmean(cores(i).condfr,2));
        coerr(i,:,:) = squeeze(nanmean(cores(i).conderr,2));
        cofrr1(i,:,:) = squeeze(nanmean(cores(i).runcondfr,2));
        cofrr0(i,:,:) = squeeze(nanmean(cores(i).stillcondfr,2));
        cor1n(i,:,:) = squeeze(nanmean(cores(i).runn,2));
        cor0n(i,:,:) = squeeze(nanmean(cores(i).stilln,2));
        for j = 1:5
            [vr,rps] = max(cores(i).runcondfr(1,:,1,j),[],2);
            [vs,sps] = max(cores(i).stillcondfr(1,:,1,j),[],2);
            cofrr1p(i,1,j) = vr; cofrr1p(i,2,j) = cores(i).runcondfr(2,rps,1,j);
            cofrr0p(i,1,j) = vs; cofrr0p(i,2,j) = cores(i).stillcondfr(2,sps,1,j);
        end
        cocondrespr1(i,:,:,:) = squeeze(nanmean(cores(i).runcondresp,2));        
        cobincondresp(i,:,:,:) = squeeze(nanmean(cores(i).bincondresp,2));
        [m,coprefcl(i)] = max(squeeze(cofr(i,1,:)));
        corely(i,:,:) = squeeze(nanmean(cores(i).condrely,2));
        corelyn(i,:,:) = squeeze(nanmean(cores(i).condrelyn,2));
        cosparseness(i,:,:) = squeeze(nanmean(cores(i).condsparseness,2));
        cosparsenesswin(i,:,:) = squeeze(nanmean(cores(i).condsparsenesswin,2));
        coff(i,:,:) =squeeze( nanmean(cores(i).ff,2));
        con(i,:,:) = squeeze(nanmean(cores(i).condn,2));
        coorifit(i,:,:,:) = squeeze(cores(i).gaussparams);
        coorifitrs(i,:,:) = squeeze(cores(i).gaussrsquared);
        coosi(i,:,:) = squeeze(cores(i).condosi);
        codsi(i,:,:) = squeeze(cores(i).conddsi);
        if ~isnan(coprefcl(i))
            coprefcfit(i,:,:) = cores(i).gaussparams(:,1,coprefcl(i),:);
            coprefcfitrs(i,:) = cores(i).gaussrsquared(:,1,coprefcl(i));
        else
            coprefcfit(i,:,:) = nan(2,5);
            coprefcfitrs(i,:) = nan(2,1);
        end        
        crfparams(i,1,:) = squeeze(cores(i).paramsl0); %[rmax, n, c50, r0]
        crfparams(i,2,:) = squeeze(cores(i).paramsl1);
        crfrsq(i,1) = squeeze(cores(i).rsql0);
        crfrsq(i,2) = squeeze(cores(i).rsql1);        
        
        coszfr(i,:,:,:) = nan(2,3,5);
        coszerr(i,:,:,:) = nan(2,3,5);
        coszfrr1(i,:,:,:) = nan(2,3,5);
        coszfrr0(i,:,:,:) = nan(2,3,5);
        coszerrr1(i,:,:,:) = nan(2,3,5);
        coszerrr0(i,:,:,:) = nan(2,3,5);
        coszr1n(i,:,:,:) = nan(2,3,5);
        coszr0n(i,:,:,:) = nan(2,3,5);
        coszfrr1p(i,:,:,:) = nan(2,3,5);
        coszfrr0p(i,:,:,:) = nan(2,3,5);
        coszcondrespr1(i,:,:,:,:) = nan(2,3,5,3000);  
        coszbincondresp(i,:,:,:,:) = nan(2,3,5,90);
        coszrely(i,:,:,:) = nan(2,3,5);
        coszrelyn(i,:,:,:) = nan(2,3,5);
        coszsparseness(i,:,:,:) = nan(2,3,5);
        coszsparsenesswin(i,:,:,:) = nan(2,3,5);
        coszff(i,:,:,:) = nan(2,3,5);
        coszn(i,:,:,:) = nan(2,3,5);
        coszorifit(i,:,:,:,:) = nan(2,3,5,5);
        coszorifitrs(i,:,:,:) = nan(2,3,5);
        coszosi(i,:,:,:) = nan(2,3,5);
        coszdsi(i,:,:,:) = nan(2,3,5);
        coszprefcfit(i,:,:,:) = nan(2,3,5);
        coszprefcfitrs(i,:,:) = nan(2,3,1);
        
    else % at preferred size
        if isnan(cores(i).depth) 
            coprefsz(i) = 1; 
        else
            coprefsz(i) = find(nanmean(nanmean(cores(i).condfr(1,:,:,:),2),4) == max(nanmean(nanmean(cores(i).condfr(1,:,:,:),2),4)),1);
        end
        cofr(i,:,:) = squeeze(nanmean(cores(i).condfr(:,:,coprefsz(i),:),2));
        coerr(i,:,:) = squeeze(nanmean(cores(i).conderr(:,:,coprefsz(i),:),2));
        cofrr1(i,:,:) = squeeze(nanmean(cores(i).runcondfr(:,:,coprefsz(i),:),2));
        cofrr0(i,:,:) = squeeze(nanmean(cores(i).stillcondfr(:,:,coprefsz(i),:),2));
        cor1n(i,:,:) = squeeze(nanmean(cores(i).runn(:,:,coprefsz(i),:),2));
        cor0n(i,:,:) = squeeze(nanmean(cores(i).stilln(:,:,coprefsz(i),:),2));
        for j = 1:5
            [vr,rps] = max(cores(i).runcondfr(1,:,coprefsz(i),j),[],2);
            [vs,sps] = max(cores(i).stillcondfr(1,:,coprefsz(i),j),[],2);
            cofrr1p(i,1,j) = vr; cofrr1p(i,2,j) = cores(i).runcondfr(2,rps,coprefsz(i),j);
            cofrr0p(i,1,j) = vs; cofrr0p(i,2,j) = cores(i).stillcondfr(2,sps,coprefsz(i),j);
        end
        cocondrespr1(i,:,:,:) = squeeze(nanmean(cores(i).runcondresp(:,:,coprefsz(i),:,:),2));          
        cobincondresp(i,:,:,:) = squeeze(nanmean(cores(i).bincondresp(:,:,coprefsz(i),:,:),2));
        [m,coprefcl(i)] = max(squeeze(cofr(i,1,:)));
        corely(i,:,:) = squeeze(nanmean(cores(i).condrely(:,:,coprefsz(i),:),2));
        corelyn(i,:,:) = squeeze(nanmean(cores(i).condrelyn(:,:,coprefsz(i),:),2));
        cosparseness(i,:,:) = squeeze(nanmean(cores(i).condsparseness(:,:,coprefsz(i),:),2));
        cosparsenesswin(i,:,:) = squeeze(nanmean(cores(i).condsparsenesswin(:,:,coprefsz(i),:),2));
        coff(i,:,:) =squeeze( nanmean(cores(i).ff(:,:,coprefsz(i),:),2));
        con(i,:,:) = squeeze(nanmean(cores(i).condn(:,:,coprefsz(i),:),2));
        coorifit(i,:,:,:) = squeeze(cores(i).gaussparams(:,coprefsz(i),:,:));
        coorifitrs(i,:,:) = squeeze(cores(i).gaussrsquared(:,coprefsz(i),:));
        coosi(i,:,:) = squeeze(cores(i).condosi(:,coprefsz(i),:));
        codsi(i,:,:) = squeeze(cores(i).conddsi(:,coprefsz(i),:));
        if ~isnan(coprefcl(i))
            coprefcfit(i,:,:) = cores(i).gaussparams(:,coprefsz(i),coprefcl(i),:);
            coprefcfitrs(i,:) = cores(i).gaussrsquared(:,coprefsz(i),coprefcl(i));
        else
            coprefcfit(i,:,:) = nan(2,5);
            coprefcfitrs(i,:) = nan(2,1);
        end        
        
        crfparams(i,1,:) = squeeze(cores(i).paramsl0(coprefsz(i),:)); %[rmax, n, c50, r0]
        crfparams(i,2,:) = squeeze(cores(i).paramsl1(coprefsz(i),:));
        crfrsq(i,1) = squeeze(cores(i).rsql0(coprefsz(i)));
        crfrsq(i,2) = squeeze(cores(i).rsql1(coprefsz(i)));
        
        coszfr(i,:,:,:) = squeeze(nanmean(cores(i).condfr,2));
        coszerr(i,:,:,:) = squeeze(nanmean(cores(i).conderr,2));
        coszfrr1(i,:,:,:) = squeeze(nanmean(cores(i).runcondfr,2));
        coszfrr0(i,:,:,:) = squeeze(nanmean(cores(i).stillcondfr,2));
        coszerrr1(i,:,:,:) = squeeze(nanmean(cores(i).runconderr,2));
        coszerrr0(i,:,:,:) = squeeze(nanmean(cores(i).stillconderr,2));
        coszr1n(i,:,:,:) = squeeze(nanmean(cores(i).runn,2));
        coszr0n(i,:,:,:) = squeeze(nanmean(cores(i).stilln,2));
        % get it at preferred orientation only, not average
        for j = 1:5 % contrasts
            for k = 1:3 %sizes
                [vr,rps] = max(cores(i).runcondfr(1,:,k,j),[],2);
                [vs,sps] = max(cores(i).stillcondfr(1,:,k,j),[],2);
                coszfrr1p(i,1,k,j) = vr; coszfrr1p(i,2,k,j) = cores(i).runcondfr(2,rps,k,j);
                coszfrr0p(i,1,k,j) = vs; coszfrr0p(i,2,k,j) = cores(i).stillcondfr(2,sps,k,j);
                coszerrr1p(i,1,k,j) = cores(i).runconderr(1,rps,k,j); coszerrr1p(i,2,k,j) = cores(i).runconderr(2,rps,k,j);
                coszerrr0p(i,1,k,j) = cores(i).stillconderr(1,rps,k,j); coszerrr0p(i,2,k,j) = cores(i).stillconderr(2,rps,k,j);
            end
        end
        coszcondrespr1(i,:,:,:,:) = squeeze(nanmean(cores(i).runcondresp,2));  
        coszbincondresp(i,:,:,:,:) = squeeze(nanmean(cores(i).bincondresp,2));
        coszrely(i,:,:,:) = squeeze(nanmean(cores(i).condrely,2));
        coszrelyn(i,:,:,:) = squeeze(nanmean(cores(i).condrelyn,2));
        coszsparseness(i,:,:,:) = squeeze(nanmean(cores(i).condsparseness,2));
        coszsparsenesswin(i,:,:,:) = squeeze(nanmean(cores(i).condsparsenesswin,2));
        coszff(i,:,:,:) =squeeze(nanmean(cores(i).ff,2));
        coszn(i,:,:,:) = squeeze(nanmean(cores(i).condn,2));
        coszorifit(i,:,:,:,:) = squeeze(cores(i).gaussparams);
        coszorifitrs(i,:,:,:) = squeeze(cores(i).gaussrsquared);
        coszosi(i,:,:,:) = squeeze(cores(i).condosi);
        coszdsi(i,:,:,:) = squeeze(cores(i).conddsi);
        if ~isnan(coprefcl(i))
            coszprefcfit(i,:,:,:) = cores(i).gaussparams(:,:,coprefcl(i),:);
            coszprefcfitrs(i,:,:) = cores(i).gaussrsquared(:,:,coprefcl(i));
        else
            coszprefcfit(i,:,:,:) = nan(2,3,5);
            coszprefcfitrs(i,:,:) = nan(2,3,1);
        end 
    end
    coctrlfr(i,:) = cores(i).controlfr;
    coctrlerr(i,:) = cores(i).controlerr;
    cor0ctrlfr(i,:) = cores(i).r0controlfr;
    cor0ctrlerr(i,:) = cores(i).r0controlerr;
    cor0ctrln(i,:) = cores(i).r0controln;
    cor1ctrlfr(i,:) = cores(i).r1controlfr;
    cor1ctrlerr(i,:) = cores(i).r1controlerr;
    cor1ctrln(i,:) = cores(i).r1controln;    
    
    mvvismod(i) = mvres(i).vismod;
    mvlightmod(i) = mvres(i).lightmod;
    mvlfr(i) = mvres(i).lfr;
    mvnlfr(i) = mvres(i).nlfr;
    if size(mvres(i).condfr,2)>2
        mvfr(i,:,:) = mvres(i).condfr(:,1:2);
        mverr(i,:,:) = mvres(i).conderr(:,1:2);
        mvfrr1(i,:,:) = nanmean(mvres(i).runcondfr(:,1:2),2);
        mvfrr0(i,:,:) = nanmean(mvres(i).stillcondfr(:,1:2),2);
        mvr1n(i,:,:) = nanmean(mvres(i).runn(:,1:2),2);
        mvr0n(i,:,:) = nanmean(mvres(i).stilln(:,1:2),2);
        mvbincondresp(i,:,:,:) = mvres(i).bincondresp(:,1:2,:);     
        mvrely(i,:,:) = mvres(i).condrely(:,1:2);    
        mvrelyn(i,:,:) = mvres(i).condrelyn(:,1:2);   
        mvsparseness(i,:,:) = mvres(i).condsparseness(:,1:2);
        mvsparsenesswin(i,:,:) = mvres(i).condsparsenesswin(:,1:2);
        mvjitter(i,:,:) = mvres(i).meanjitsd(:,1:2);
        mvff(i,:,:) = mvres(i).ff(:,1:2);
        mvn(i,:,:) = mvres(i).condn(:,1:2);
    else
        mvfr(i,:,:) = mvres(i).condfr;
        mverr(i,:,:) = mvres(i).conderr;
        mvfrr1(i,:,:) = nanmean(mvres(i).runcondfr,2);
        mvfrr0(i,:,:) = nanmean(mvres(i).stillcondfr,2);
        mvr1n(i,:,:) = nanmean(mvres(i).runn,2);
        mvr0n(i,:,:) = nanmean(mvres(i).stilln,2);
        mvbincondresp(i,:,:,:) = mvres(i).bincondresp;
        mvrely(i,:,:) = mvres(i).condrely;    
        mvrelyn(i,:,:) = mvres(i).condrelyn;
        mvsparseness(i,:,:) = mvres(i).condsparseness;
        mvsparsenesswin(i,:,:) = mvres(i).condsparsenesswin;
        mvjitter(i,:,:) = mvres(i).meanjitsd;
        mvff(i,:,:) = mvres(i).ff;
        mvn(i,:,:) = mvres(i).condn;
    end
    
    srvismod(i) = srres(i).vismod;
    srlightmod(i) = srres(i).lightmod;
    srlfr(i) = srres(i).lfr;
    srnlfr(i) = srres(i).nlfr;
    srfr(i,:,:) = nanmean(srres(i).condfr,2);
    srerr(i,:,:) = nanmean(srres(i).conderr,2);
    srfrr1(i,:,:) = nanmean(srres(i).runcondfr,2);
    srfrr0(i,:,:) = nanmean(srres(i).stillcondfr,2);
    srerrr1(i,:,:) = nanmean(srres(i).runconderr,2);
    srerrr0(i,:,:) = nanmean(srres(i).stillconderr,2);
    srr1n(i,:,:) = nanmean(srres(i).runn,2);
    srr0n(i,:,:) = nanmean(srres(i).stilln,2);
    for j = 1:4
        [vr,rps] = max(srres(i).runcondfr(1,:,j),[],2);
        [vs,sps] = max(srres(i).stillcondfr(1,:,j),[],2);
        srfrr1p(i,1,j) = vr; srfrr1p(i,2,j) = srres(i).runcondfr(2,rps,j);
        srfrr0p(i,1,j) = vs; srfrr0p(i,2,j) = srres(i).stillcondfr(2,sps,j);
    end
    srcondrespr1(i,:,:,:) = squeeze(nanmean(srres(i).runcondresp,2));
    srbincondresp(i,:,:,:) = nanmean(srres(i).bincondresp,2);
    srrely(i,:,:) = nanmean(srres(i).condrely,2);
    srrelyn(i,:,:) = nanmean(srres(i).condrelyn,2);
    srsparseness(i,:,:) = nanmean(srres(i).condsparseness,2);
    srsparsenesswin(i,:,:) = nanmean(srres(i).condsparsenesswin,2);
    srff(i,:,:) = nanmean(srres(i).ff,2);
    srn(i,:,:) = nanmean(srres(i).condn,2);
    srlfpspect(i,:,:,:) = squeeze(nanmean(srres(i).condlfpspect,2));
    srlfpspecterr(i,:,:,:) = squeeze(nanmean(srres(i).condlfpspecterr,2));
    sromlfpspect(i,:,:,:) = srres(i).omlfpspect;
    sromlfpspecterr(i,:,:,:) = srres(i).omlfpspecterr;
    fax = srres(i).fax;    
    
    phvismod(i) = phres(i).vismod;
    phlightmod(i) = phres(i).lightmod;
    phlfr(i) = phres(i).lfr;
    phnlfr(i) = phres(i).nlfr;
    phfr(i,:,:) = nanmean(phres(i).condfr,2);
    pherr(i,:,:) = nanmean(phres(i).conderr,2);
    phfrr1(i,:,:) = nanmean(phres(i).runcondfr,2);
    phfrr0(i,:,:) = nanmean(phres(i).stillcondfr,2);
    pherrr1(i,:,:) = nanmean(phres(i).runconderr,2);
    pherrr0(i,:,:) = nanmean(phres(i).stillconderr,2);
    phr1n(i,:,:) = nanmean(phres(i).runn,2);
    phr0n(i,:,:) = nanmean(phres(i).stilln,2);
    for j = 1:4
        [vr,rps] = max(phres(i).runcondfr(1,:,j),[],2);
        [vs,sps] = max(phres(i).stillcondfr(1,:,j),[],2);
        phfrr1p(i,1,j) = vr; phfrr1p(i,2,j) = phres(i).runcondfr(2,rps,j);
        phfrr0p(i,1,j) = vs; phfrr0p(i,2,j) = phres(i).stillcondfr(2,sps,j);
    end
    phcondrespr1(i,:,:,:) = squeeze(nanmean(phres(i).runcondresp,2));
    phbincondresp(i,:,:,:) = nanmean(phres(i).bincondresp,2);
    phrely(i,:,:) = nanmean(phres(i).condrely,2);
    phrelyn(i,:,:) = nanmean(phres(i).condrelyn,2);
    phsparseness(i,:,:) = nanmean(phres(i).condsparseness,2);
    phsparsenesswin(i,:,:) = nanmean(phres(i).condsparsenesswin,2);
    phff(i,:,:) = nanmean(phres(i).ff,2);
    phn(i,:,:) = nanmean(phres(i).condn,2);
    phlfpspect(i,:,:,:) = squeeze(nanmean(phres(i).condlfpspect,2));
    phlfpspecterr(i,:,:,:) = squeeze(nanmean(phres(i).condlfpspecterr,2));
    phomlfpspect(i,:,:,:) = phres(i).omlfpspect;
    phomlfpspecterr(i,:,:,:) = phres(i).omlfpspecterr;
    
    for sz = 1:2
        [mvlinp(i,sz,:),S] = polyfit(squeeze(mvres(i).bincondresp(1,sz,40:69)),squeeze(mvres(i).bincondresp(2,sz,40:69)),1);
        mvnormr(i,sz) = S.normr;
        c = corrcoef(squeeze(mvres(i).bincondresp(1,sz,40:69)),squeeze(mvres(i).bincondresp(2,sz,40:69)));
        mvrsq(i,sz) = c(1,2)^2;
    end
    for sz = 1:5
        a = squeeze(szres(i).bincondresp(1,:,sz,25:54)); a = a(:);
        b = squeeze(szres(i).bincondresp(2,:,sz,25:54)); b = b(:);
        [szlinp(i,sz,:),S] = polyfit(a,b,1);
        sznormr(i,sz) = S.normr;
        c = corrcoef(a,b);
        szrsq(i,sz) = c(1,2)^2;
    end
    a = squeeze(szres(i).condfr(1,:,:)); a = a(:);
    hlp = squeeze(szres(i).condfr(2,:,:)./szres(i).condfr(1,:,:)); % percent change
    b = log10(hlp(:));
    ac = a; bc = b; % copy
    ac(find(isinf(b) | isnan(a) | isnan(b))) = [];
    bc(find(isinf(b) | isnan(a) | isnan(b))) = [];
    acc{i} = ac;
    bcc{i} = bc;
    if length(ac)>2
        [szacvcf(i,:),S] = polyfit(ac,bc,1); % sz all conds vs change fit
        c = corrcoef(ac,bc);
        szacvcfrsq(i) = c(1,2)^2;
    else
        szacvcf(i,:) = nan(1,2);
        szacvcfrsq(i) = NaN;
    end
    
    onfitrs(i) = rfres(i).rson;
    offfitrs(i) = rfres(i).rsoff;
    
    % ON OFF overlap?
    onoffd(i) = sqrt((rfres(i).gaussfiton.shiftxcenterDeg - rfres(i).gaussfitoff.shiftxcenterDeg).^2 + (rfres(i).gaussfiton.shiftycenterDeg - rfres(i).gaussfitoff.shiftycenterDeg).^2);
    meanspreadon(i) = (rfres(i).gaussfiton.xspreadDeg + rfres(i).gaussfiton.yspreadDeg)/2;
    meanspreadoff(i) = (rfres(i).gaussfitoff.xspreadDeg + rfres(i).gaussfitoff.yspreadDeg)/2;
    onoffoverlap(i) = (meanspreadon(i) + meanspreadoff(i) - onoffd(i)) ./ (meanspreadon(i)+meanspreadoff(i)+onoffd(i));
    
    % for overlap of movie aperture to RF
    ond(i) = sqrt((rfres(i).gaussfiton.shiftxcenterDeg - mvres(i).position(1)).^2 + (rfres(i).gaussfiton.shiftycenterDeg - mvres(i).position(2)).^2);
    offd(i) = sqrt((rfres(i).gaussfitoff.shiftxcenterDeg - mvres(i).position(1)).^2 + (rfres(i).gaussfitoff.shiftycenterDeg - mvres(i).position(2)).^2);
    onoverlap(i) = (meanspreadon(i) + max(unique(mvres(i).widths))/2-ond(i)) ./ (meanspreadon(i)+max(unique(mvres(i).widths))/2+ond(i));
    offoverlap(i) = (meanspreadoff(i) + max(unique(mvres(i).widths))/2-offd(i)) ./ (meanspreadoff(i)+max(unique(mvres(i).widths))/2+offd(i));
    
    onamp(i) = rfres(i).gaussfiton.amp; % this is amplitude from gaussfit.
    offamp(i) = rfres(i).gaussfitoff.amp; 
    fronfield(i) = rfres(i).fron; % this is from pixels with 2sd over baseline (inside the RF)
    frofffield(i) = rfres(i).froff;
    rfmfron(i) = nanmean(nanmean(rfres(i).on)); % this is mean firing rate across ALL pixels
    rfmfroff(i) = nanmean(nanmean(rfres(i).off));
    onsnr(i) = rfres(i).snron; % if this is NaN, no pixel or only single ones were over 2sd
    offsnr(i) = rfres(i).snroff;
    
    % overlap of center surround grating with RF    cs = center surround
    ondcs(i) = sqrt((rfres(i).gaussfiton.shiftxcenterDeg - srres(i).position(1)).^2 + (rfres(i).gaussfiton.shiftycenterDeg - srres(i).position(2)).^2);
    offdcs(i) = sqrt((rfres(i).gaussfitoff.shiftxcenterDeg - srres(i).position(1)).^2 + (rfres(i).gaussfitoff.shiftycenterDeg - srres(i).position(2)).^2);
    onoverlapcs(i) = (meanspreadon(i) + 15/2-ond(i)) ./ (meanspreadon(i)+15/2+ond(i)); % TODO put in actual width!!!!
    offoverlapcs(i) = (meanspreadoff(i) + 15/2-offd(i)) ./ (meanspreadoff(i)+15/2+offd(i));
    
    %LFP values - movie
    lfprson(i) = rfres(i).lfprson;
    lfprsoff(i) = rfres(i).lfprsoff;
    % is LFP RF overlapping with movie aperture?
    lfpond(i) = sqrt((rfres(i).lfpfiton.shiftxcenterDeg - mvres(i).position(1)).^2 + (rfres(i).lfpfiton.shiftycenterDeg - mvres(i).position(2)).^2);
    lfpoffd(i) = sqrt((rfres(i).lfpfitoff.shiftxcenterDeg - mvres(i).position(1)).^2 + (rfres(i).lfpfitoff.shiftycenterDeg - mvres(i).position(2)).^2);
    lfpmeanspreadon(i) = (rfres(i).lfpfiton.xspreadDeg + rfres(i).lfpfiton.yspreadDeg)/2;
    lfpmeanspreadoff(i) = (rfres(i).lfpfitoff.xspreadDeg + rfres(i).lfpfitoff.yspreadDeg)/2;    
    lfponoverlap(i) = (lfpmeanspreadon(i) + max(unique(mvres(i).widths))/2-ond(i)) ./ (lfpmeanspreadon(i)+max(unique(mvres(i).widths))/2+ond(i));
    lfpoffoverlap(i) = (lfpmeanspreadoff(i) + max(unique(mvres(i).widths))/2-offd(i)) ./ (lfpmeanspreadoff(i)+max(unique(mvres(i).widths))/2+offd(i));
    
    lfponamp(i) = rfres(i).lfpfiton.amp; % this is amplitude from gaussfit.
    lfpoffamp(i) = rfres(i).lfpfitoff.amp; 
    

end
szvismod(isnan(szvismod)) = 0; covismod(isnan(covismod)) = 0; mvvismod(isnan(mvvismod)) = 0;
szbta = szres(1).bta;

%%
% categorize cells
secpersamp = 1/30000;
interpf = secpersamp/10;
swidthms = swidth*interpf*1000;
prsv = swidthms>=.38; pfsv = swidthms<=.36;
prs = find(prsv); pfs = find(pfsv);
depth = depth.*cosd(22).*cosd(pangle);

anids = animalids; % TODO check why srres(161) is an animal not in the mix!! ok for filename just animalid is off: also rfres is off for animal 21 (280) in the SOM populaiton??
ids = zeros(length(animalids),length(depth));
for i = 1:length(animalids)
    for j = 1:length(depth)
        ids(i,j) = strcmp(animalids{i},rfres(j).animalid);
    end
    ids(i,:)= ids(i,:).*i;
end
aids = sum(ids,1);
anls = unique(animalno);
% anls = unique(aids);
for i = 1:length(anls)
    thisan = find(ids(i,:));
%     thisan = find(animalno == anls(i));
    [a,ind] = min(300-depth(thisan));
    lfpinds(i) = thisan(ind);
end

% % for big population
% som = zeros(1,length(depth));
% pv = zeros(1,length(depth));
% som(1:397) = 1;
% pv(398:end) = 1;

% halo expressing clls
phe = zeros(1,length(depth));

% % putative halo expressing for SOM later pop: 
% phe([4, 7, 37, 64, 66, 87, 95, 106, 119, 163, 201, 244, 283, 302, 329]) = 1;
phe([7, 95, 106, 119, 163, 201, 244, 283, 302, 329, 411]) = 1;  %163 weird cell diff fr size affected, contrast not

% % putative halo expressing for PV Halo pop: 
% phe([44, 85, 159, 174, 196, 244, 259, 274, 278, 281, 285, 287, 288, 298, 301, 311, 314]) = 1;

% %putative halo for SOM + PV combined
% phe([4, 7, 37, 64, 66, 87, 95, 106, 119, 163, 201, 244, 283, 302, 329,   441, 482, 556, 571, 593, 641, 656, 671, 675, 678, 682, 684, 685, 695, 698, 708, 711]) = 1;

phe = logical(phe)';

l23 = depth<375;
l4 = depth>=375&depth<=550;
l5 = depth>550&depth<=800;
l5a = depth>550&depth<=650;
l5b = depth>650&depth<=800;
l6 = depth>800;
l23rs = l23&prsv&~phe';
l23fs = l23&pfsv&~phe';
l4rs = l4&prsv&~phe';
l4fs = l4&pfsv&~phe';
l5rs = l5&prsv&~phe';
l5fs = l5&pfsv&~phe';
l6rs = l6&prsv&~phe';
l6fs = l6&pfsv&~phe';
okrs = prsv&~phe';

% figure out which cells to include for aperture fit, firing rates,...

% figure out if RF mappable and overlapping with movie aperture
onfits = [rfres.gaussfiton]; offfits = [rfres.gaussfitoff];
okrfon = [onfits.rsquared]>.5 & meanspreadon>1.5; %[gaussfiton.amp]>1 & 
okrfoff = [offfits.rsquared]>.5 & meanspreadoff>1.5; %[gaussfitoff.amp]>1 & 
nonmappable = ~okrfon&~okrfoff;
apok = (okrfon&onoverlap>.33) | (okrfoff & offoverlap>.33);

% is it overlapping with center surround? 
csapok = (okrfon&onoverlapcs>.33) | (okrfoff & offoverlapcs>.33);
csapnotok = (okrfon&onoverlapcs<.33) & (okrfoff & offoverlapcs<.33);

% which tuning curves fits are ok? 
for l = 1:2
    for sz = 1:5
        okoriszfits(:,l,sz) = szorifitrs(:,l,sz)>=.75;
        okoricofits(:,l,sz) = coorifitrs(:,l,sz)>=.75;
    end
end

% figure out which had enough spikes for reliability estimates
mvrelyok = mvrelyn>mvn; % need at least as many spikes as number of trials that were cross correlated
szrelyok = szrelyn>szn;
corelyok = corelyn>con;


gamma = [15,40];
g1 = find(fax>gamma(1),1); g2 = find(fax>gamma(2),1);

% LFP spectra
for i = 1:length(lfpinds)
    for l = 1:2
        for c = 1:4            
            % find peaks for all conditions             
            hh = smooth(squeeze(sromlfpspect(lfpinds(i),l,c,1:103)));
            gsig = hh(g1:g2);
            if isempty(find(diff(gsig)>0)) % there is no clear beta peak
                gpi = round((g1+g2)/2);
%                 gpi = NaN;
                haspeak(i,l,c) = 0;
            else
                peaks = find(diff(gsig)>0)+1;
                pvs = gsig(peaks);
                gpi = peaks(pvs == max(pvs));
                gpi = gpi+g1-1;
                haspeak(i,l,c) = 1;
            end            
            gi(i,l,c) = gpi;
            if ~isnan(gi(i,1,1))
                gp(i,l,c) = sromlfpspect(lfpinds(i),l,c,gi(i,1,1));
            else
                gp(i,l,c) = NaN;
            end
            
            filly(i,l,c,:) = [squeeze(sromlfpspect(lfpinds(i),l,c,1:103)+sromlfpspecterr(lfpinds(i),l,c,1:103))',fliplr(squeeze(sromlfpspect(lfpinds(i),l,c,1:103)-sromlfpspecterr(lfpinds(i),l,c,1:103))')];
        end
        for sz = 1:5         
            % find peaks for all conditions             
            hh = smooth(squeeze(szomlfpspect(lfpinds(i),l,sz,1:103)));
            gsig = hh(g1:g2);
            if isempty(find(diff(gsig)>0)) % there is no clear beta peak
                gpi = round((g1+g2)/2);
                haspeaksz(i,l,sz) = 0;
%                 gpi = NaN;
            else
                peaks = find(diff(gsig)>0)+1;
                pvs = gsig(peaks);
                gpi = peaks(pvs == max(pvs));
                gpi = gpi+g1-1;
                haspeaksz(i,l,sz) = 1;
            end            
            gisz(i,l,sz) = gpi;
            szfilly(i,l,sz,:) = [squeeze(szomlfpspect(lfpinds(i),l,sz,1:103)+szomlfpspecterr(lfpinds(i),l,sz,1:103))',fliplr(squeeze(szomlfpspect(lfpinds(i),l,sz,1:103)-szomlfpspecterr(lfpinds(i),l,sz,1:103))')];
        end
        for sz = 1:5
            if ~isnan(gisz(i,1,5))
                gpsz(i,l,sz) = szomlfpspect(lfpinds(i),l,sz,gisz(i,1,5));
            else
                gpsz(i,l,sz) = NaN;
            end
        end           
    end
end
fillx = [fax(1:103),fliplr(fax(1:103))];


figure
plot(1,squeeze(gp(:,1,1)),'ko','markerfacecolor','k')
hold on
plot(2,squeeze(gp(:,1,2)),'o','color',[.5,.5,.5],'markerfacecolor',[.5,.5,.5])
axis([0,3,0,400])
for i = 1:size(gp,1)
    line([1,2],[gp(i,1,1),gp(i,1,2)],'color','k')
end
set(gca,'xticklabel',{'iso','cross'})
set(gca,'xtick',[1,2]);
ylabel('peak gamma power')
set(gcf,'OuterPosition',[573   504   242   513])

figure
plot(1,squeeze(gp(:,1,1)),'ko','markerfacecolor','k')
hold on
plot(2,squeeze(gp(:,2,1)),'o','color',[1,0,0],'markerfacecolor',[1,0,0])
axis([0,3,0,400])
for i = 1:size(gp,1)
    line([1,2],[gp(i,1,1),gp(i,2,1)],'color','k')
end
set(gca,'xticklabel',{'iso control','iso light'})
set(gca,'xtick',[1,2]);
ylabel('peak gamma power')
set(gcf,'OuterPosition',[573   504   242   513])

figure
plot(1,squeeze(gp(:,1,2)),'ko','markerfacecolor','k')
hold on
plot(2,squeeze(gp(:,2,2)),'o','color',[1,0,0],'markerfacecolor',[1,0,0])
axis([0,3,00,400])
for i = 1:size(gp,1)
    line([1,2],[gp(i,1,2),gp(i,2,2)],'color','k')
end
set(gca,'xticklabel',{'cross control','cross light'})
set(gca,'xtick',[1,2]);
ylabel('peak gamma power')
set(gcf,'OuterPosition',[573   504   242   513])

gpdiff = squeeze(gp(:,2,:)-gp(:,1,:));
gpratio = squeeze(gp(:,2,:)./gp(:,1,:));

figure
plot(1,squeeze(gpdiff(:,1)),'ko','markerfacecolor','k')
hold on
plot(2,squeeze(gpdiff(:,2)),'o','color',[.5,.5,.5],'markerfacecolor',[.5,.5,.5])
axis([0,3,0,160])
for i = 1:size(gpdiff,1)
    line([1,2],[gpdiff(i,1),gpdiff(i,2)],'color','k')
end
set(gca,'xticklabel',{'delta power iso','delta power cross'})
set(gca,'xtick',[1,2]);
ylabel('peak gamma power')
set(gcf,'OuterPosition',[573   504   242   513])

figure
plot(1,squeeze(gpratio(:,1)),'ko','markerfacecolor','k')
hold on
plot(2,squeeze(gpratio(:,2)),'o','color',[.5,.5,.5],'markerfacecolor',[.5,.5,.5])
axis([0,3,1,2.1])
for i = 1:size(gpratio,1)
    line([1,2],[gpratio(i,1),gpratio(i,2)],'color','k')
end
set(gca,'xticklabel',{'power ratio iso','power ratio cross'})
set(gca,'xtick',[1,2]);
ylabel('peak gamma power')
set(gcf,'OuterPosition',[573   504   242   513])

% for size
figure
plot(1,squeeze(gpsz(:,1,1)),'ko','markerfacecolor','k')
hold on
plot(2,squeeze(gpsz(:,1,5)),'o','color',[.5,.5,.5],'markerfacecolor',[.5,.5,.5])
axis([0,3,0,400])
for i = 1:size(gpsz,1)
    line([1,2],[gpsz(i,1,1),gpsz(i,1,5)],'color','k')
end
set(gca,'xticklabel',{'small','large'})
set(gca,'xtick',[1,2]);
ylabel('peak gamma power')
set(gcf,'OuterPosition',[573   504   242   513])

figure
plot(1,squeeze(gpsz(:,1,1)),'ko','markerfacecolor','k')
hold on
plot(2,squeeze(gpsz(:,2,1)),'o','color',[1,0,0],'markerfacecolor',[1,0,0])
axis([0,3,0,400])
for i = 1:size(gpsz,1)
    line([1,2],[gpsz(i,1,1),gpsz(i,2,1)],'color','k')
end
set(gca,'xticklabel',{'small control','small light'})
set(gca,'xtick',[1,2]);
ylabel('peak gamma power')
set(gcf,'OuterPosition',[573   504   242   513])

figure
plot(1,squeeze(gpsz(:,1,5)),'ko','markerfacecolor','k')
hold on
plot(2,squeeze(gpsz(:,2,5)),'o','color',[1,0,0],'markerfacecolor',[1,0,0])
axis([0,3,0,400])
for i = 1:size(gpsz,1)
    line([1,2],[gpsz(i,1,5),gpsz(i,2,5)],'color','k')
end
set(gca,'xticklabel',{'large control','large light'})
set(gca,'xtick',[1,2]);
ylabel('peak gamma power')
set(gcf,'OuterPosition',[573   504   242   513])

gpdiff = squeeze(gpsz(:,2,:)-gpsz(:,1,:));
gpratio = squeeze(gpsz(:,2,:)./gpsz(:,1,:));

for i = 1:size(filly,1)
    figure
    fill(fillx,squeeze(filly(i,1,1,:))','b')
    hold on
    fill(fillx,squeeze(filly(i,1,2,:))','c')
    fill(fillx,squeeze(filly(i,2,1,:))','r')
    fill(fillx,squeeze(filly(i,2,2,:))','m')
% 
%     fill(fillx,squeeze(filly(i,1,3,:))','k')
%     fill(fillx,squeeze(filly(i,1,4,:))','g')
    set(gca,'yscale','log')
end

    

% % correlations weird
% an = 2;
% a = find(animalno == an);
% for i = 1:length(a)
%     for j = 1:length(a)
%         for l = 1:2
%             for s = 1:5
%                 for o = 1:size(szres(a(i)).allresp,2)
%                     cc = [];
%                     for rep = 1:min(size(szres(a(i)).allresp{l,o,s},1),size(szres(a(j)).allresp{l,o,s},1))
%                         [cc(rep,:),lags] = xcorr(squeeze(szres(a(i)).allresp{l,o,s}(rep,:)),squeeze(szres(a(j)).allresp{l,o,s}(rep,:)),100);
%                     end
%                     if ~isempty(cc)
%                         xco(i,j,l,s,o,:) = nanmean(cc,1);
%                     else
%                         xco(i,j,l,s,o,:) = nan(1,201);
%                     end
%                 end
%             end
%         end
%     end
% end
            

% figures
% for hillels VIP grant
for i = 1:length(depth)
    [mx,prefs(i)] = max(szfrr1(i,1,:),[],3);
    [mx,prefc(i)] = max(cofrr1(i,1,:),[],3);
end
norml23rs = prefs ~= 5 & prefs ~=4 & (prefc == 4 | prefc == 5) & l23rs;


% center surround

%make raster plots
%get all different orientation in one matrix
ex = 358;
ta = [-299:2700]; prestim = 300;
isol0 = []; isol1 = []; crossl0 = []; crossl1 = [];
for i = 1:size(srres(ex).allresp,2)
    isol0 = [isol0; srres(ex).allresp{1,i,1}];
    isol1 = [isol1; srres(ex).allresp{2,i,1}];
    crossl0 = [crossl0; srres(ex).allresp{1,i,2}];
    crossl1 = [crossl1; srres(ex).allresp{2,i,2}];
end

cond1 = isol0; cond2 = crossl0;
figure
subplot(2,1,1)
hold on;
for i = 1:size(cond1,1)
    if ~isempty(find(cond1(i,:)))
        plot(find(cond1(i,:))-prestim,i,'ko','MarkerSize',1.5,'MarkerFaceColor','k')
    end
end
line([0,2000],[42,42],'color','k','linewidth',3);
% line([0,2000],[22,22],'color','k','linewidth',3);
% line([500,1500],[21,21],'color','r','linewidth',3);
axis([-300,2500,0,43])
subplot(2,1,2)
hold on;
for i = 1:size(cond2,1)
    if ~isempty(find(cond2(i,:)))
        plot(find(cond2(i,:))-prestim,i,'ko','MarkerSize',1.5,'MarkerFaceColor','k')
    end
end

frdiff = squeeze(srfrr1(:,2,:)-srfrr1(:,1,:));
pci = squeeze(srfrr1(:,2,:)./srfrr1(:,1,:));
pci(isinf(pci)) = NaN;
lpci = log2(pci);
lpci(isinf(lpci)) = NaN;
omi = (squeeze(srfrr1(:,2,:))-squeeze(srfrr1(:,1,:)))./(squeeze(srfrr1(:,2,:))+squeeze(srfrr1(:,1,:)));

% iso vs cross
figure
plot(srfrr1(l23rs,1,2),srfrr1(l23rs,1,1),'k.','markersize',15);
axis([0,15,0,15])
axis square
refline(1,0);
xlabel('firing rate cross - control (Hz)')
ylabel('firing rate iso - control')

% iso vs light
figure
plot(srfrr1(l23rs,1,1),srfrr1(l23rs,2,1),'k.');
axis([0,25,0,25])
axis square
refline(1,0);
xlabel('firing rate iso - control (Hz)')
ylabel('firing rate iso - light (Hz)')

% cross vs light
figure
plot(srfrr1(l23rs,1,2),srfrr1(l23rs,2,2),'k.');
axis([0,35,0,35])
axis square
refline(1,0);
xlabel('firing rate cross - control (Hz)')
ylabel('firing rate cross - light (Hz)')

figure
plot(1,lpci(l23rs,1),'k.','markersize',15)
hold on
plot(2,lpci(l23rs,2),'k.','markersize',15)
axis([0,3,-3,3])
a = find(l23rs);
for i = 1:length(a)
    line([1,2],[lpci(a(i),1),lpci(a(i),2)],'color','k')
end
line([0,3],[0,0],'color','k')
set(gca,'xtick',[1,2])
set(gca,'xticklabel',{'iso','cross'})
ylabel('OMR')

meas = frdiff-repmat(spontdiff,1,4);
cllt = l23rs;
figure
errorbar(nanmean(meas(cllt,:)),nanstd(meas(cllt,:))./sqrt(length(find(cllt))),'ko','linewidth',2,'markerfacecolor','k')
set(gca,'xtick',[1,2,3,4])
set(gca,'xticklabel',{'iso','cross','center','aperture'})
ylabel('delta FR')

%single units
binwidth = 33.33;
for i = 1:size(srcondrespr1,1)
    for j = 1:size(srcondrespr1,2)
        for k = 1:size(srcondrespr1,3)
            srbinrespr1(i,j,k,:) = binit(squeeze(srcondrespr1(i,j,k,:)),binwidth);
        end
    end
end
a = find(l23rs);
for i = 1:length(a)
    figure
    plot(szbta,squeeze(srbinrespr1(a(i),1,1,:)),'k','linewidth',2)
    hold on
    plot(szbta,squeeze(srbinrespr1(a(i),1,2,:)),'color',[.5,.5,.5],'linewidth',2)
    plot(szbta,squeeze(srbinrespr1(a(i),2,1,:)),'r','linewidth',2)
    plot(szbta,squeeze(srbinrespr1(a(i),2,2,:)),'color',[1,.5,.5],'linewidth',2)
    legend('iso','cross','iso - L1','cross - L1')
    ax = axis;
    axis([-300,2700,ax(3),ax(4)])
    line([0,0],[ax(3),ax(4)],'color','k')
    line([2000,2000],[ax(3),ax(4)],'color','k')
    line([500,500],[ax(3),ax(4)],'color','r')
    line([1500,1500],[ax(3),ax(4)],'color','r')
    xlabel('time(ms)')
    ylabel('firing rate (Hz)')
end
    

% size tuning
% raster plots
ex = 358;
ta = [-299:2700]; prestim = 300;
smalll0 = []; smalll1 = []; largel0 = []; largel1 = [];
for i = 1:size(srres(ex).allresp,2)
    smalll0 = [smalll0; szres(ex).allresp{1,i,2}];
    smalll1 = [smalll1; szres(ex).allresp{2,i,2}];
    largel0 = [largel0; szres(ex).allresp{1,i,5}];
    largel1 = [largel1; szres(ex).allresp{2,i,5}];
end

cond1 = smalll0; cond2 = largel0;
figure
subplot(2,1,1)
hold on;
for i = 1:size(cond1,1)
    if ~isempty(find(cond1(i,:)))
        plot(find(cond1(i,:))-prestim,i,'ko','MarkerSize',1.5,'MarkerFaceColor','k')
    end
end
line([0,2000],[42,42],'color','k','linewidth',3);
% line([0,2000],[22,22],'color','k','linewidth',3);
% line([500,1500],[21,21],'color','r','linewidth',3);
axis([-300,2500,0,43])
subplot(2,1,2)
hold on;
for i = 1:size(cond2,1)
    if ~isempty(find(cond2(i,:)))
        plot(find(cond2(i,:))-prestim,i,'ko','MarkerSize',1.5,'MarkerFaceColor','k')
    end
end
axis([-300,2500,0,43])

frdiff = squeeze(szfrr1(:,2,:)-szfrr1(:,1,:));
pci = squeeze(szfrr1(:,2,:)./szfrr1(:,1,:));
pci(isinf(pci)) = NaN;
lpci = log2(pci);
lpci(isinf(lpci)) = NaN;
omi = (squeeze(szfrr1(:,2,:))-squeeze(szfrr1(:,1,:)))./(squeeze(szfrr1(:,2,:))+squeeze(szfrr1(:,1,:)));

meas = lpci;
cllt = l23rs;
figure
errorbar(nanmean(meas(cllt,:)),nanstd(meas(cllt,:))./sqrt(length(find(cllt))),'ko-','linewidth',2,'markerfacecolor','k')
set(gca,'xtick',[1,2,3,4,5])
set(gca,'xticklabel',szres(1).sizes)
xlabel('size (degrees)')
ylabel('delta FR')


for i = 1:length(depth)
    [mx, prefsr1(i)] = max(szfrr1(i,1,:));
    nmszfr(i,:,:) = szfrr1(i,:,:)./mx;
    nctrl(i,:) = szr1ctrlfr(i,:)./mx;
    for j = 1:2
        si(i,j) = (max(szfrr1(i,j,:))-szfrr1(i,j,5))/(max(szfrr1(i,j,:))-szr1ctrlfr(i,j));
        sinbl(i,j) = (max(szfrr1(i,j,:))-szfrr1(i,j,5))/(max(szfrr1(i,j,:)));
    end
end
normszfr = cat(3,nctrl,nmszfr);
normszfr(isinf(normszfr)) = NaN;

cllt = l23rs;
figure
errorbar([0,szres(1).sizes], squeeze(nanmean(normszfr(cllt&(prefsr1~=4&prefsr1~=5),1,:))),squeeze(nanstd(normszfr(cllt&(prefsr1~=4&prefsr1~=5),1,:)))./sqrt(length(find(cllt&(prefsr1~=4&prefsr1~=5)))),'ko-','linewidth',2,'markerfacecolor','k')
hold on
errorbar([0,szres(1).sizes], squeeze(nanmean(normszfr(cllt&(prefsr1~=4&prefsr1~=5),2,:))),squeeze(nanstd(normszfr(cllt&(prefsr1~=4&prefsr1~=5),2,:)))./sqrt(length(find(cllt&(prefsr1~=4&prefsr1~=5)))),'ro-','linewidth',2,'markerfacecolor','r')
legend('control','light')
xlabel('size (degrees)')
ylabel('normalized FR')

figure
errorbar([0,szres(1).sizes], squeeze(nanmean(normszfr(cllt,1,:))),squeeze(nanstd(normszfr(cllt,1,:)))./sqrt(length(find(cllt))),'ko-','linewidth',2,'markerfacecolor','k')
hold on
errorbar([0,szres(1).sizes], squeeze(nanmean(normszfr(cllt,2,:))),squeeze(nanstd(normszfr(cllt,2,:)))./sqrt(length(find(cllt))),'ro-','linewidth',2,'markerfacecolor','r')
legend('control','light')
xlabel('size (degrees)')
ylabel('normalized FR')


% contrast@!
frdiff = squeeze(cofrr1(:,2,:)-cofrr1(:,1,:));
pci = squeeze(cofrr1(:,2,:)./cofrr1(:,1,:));
pci(isinf(pci)) = NaN;
lpci = log2(pci);
lpci(isinf(lpci)) = NaN;
omi = (squeeze(cofrr1(:,2,:))-squeeze(cofrr1(:,1,:)))./(squeeze(cofrr1(:,2,:))+squeeze(cofrr1(:,1,:)));

meas = lpci;
cllt = l23rs;
figure
errorbar(cores(1).clevels, squeeze(nanmean(meas(cllt,:))), squeeze(nanstd(meas(cllt,:)))./sqrt(length(find(cllt))),'ko-','linewidth',2,'markerfacecolor','k')
xlabel('contrast');

for i = 1:length(depth)
    crf(i,:) = [cor1ctrlfr(i,1);squeeze(cofrr1(i,1,:))];
    blscrf(i,:) = squeeze(cofrr1(i,1,:))-cor1ctrlfr(i,1);
    ncfr(i,:) = crf(i,:)./max(crf(i,:));
    if max(blscrf(i,:) <=0)
        nblscrf(i,:) = nan(1,5);
    else
        nblscrf(i,:) = blscrf(i,:)./max(blscrf(i,:));
    end
end

figure
plot(1,lpci(l23rs,1),'k.','markersize',15)
hold on
plot(2,lpci(l23rs,2),'k.','markersize',15)
axis([0,3,-4,3])
a = find(l23rs);
for i = 1:length(a)
    line([1,2],[lpci(a(i),1),lpci(a(i),2)],'color','k')
end
line([0,3],[0,0],'color','k')
set(gca,'xtick',[1,2])
set(gca,'xticklabel',{'low','high'})
ylabel('OMR')

%make raster plots
%get all different orientation in one matrix
ex = 11;
ta = [-299:2700]; prestim = 300;
lowl0 = []; lowl1 = []; highl0 = []; highl1 = [];
for i = 1:size(cores(ex).allresp,2)
    lowl0 = [lowl0; cores(ex).allresp{1,i,1,1}];
    lowl1 = [lowl1; cores(ex).allresp{2,i,1,1}];
    highl0 = [highl0; cores(ex).allresp{1,i,1,5}];
    highl1 = [highl1; cores(ex).allresp{2,i,1,5}];
end

cond1 = highl0; cond2 = highl1;
figure
subplot(2,1,1)
hold on;
for i = 1:size(cond1,1)
    if ~isempty(find(cond1(i,:)))
        plot(find(cond1(i,:))-prestim,i,'ko','MarkerSize',1.5,'MarkerFaceColor','k')
    end
end
line([0,2000],[22,22],'color','k','linewidth',3);
line([500,1500],[21,21],'color','r','linewidth',3);
axis([-300,2500,0,23])
subplot(2,1,2)
hold on;
for i = 1:size(cond2,1)
    if ~isempty(find(cond2(i,:)))
        plot(find(cond2(i,:))-prestim,i,'ro','MarkerSize',1.5,'MarkerFaceColor','r')
    end
end
axis([-300,2500,0,20])

% contrast - size interaction
frdiff = squeeze(coszfrr1p(:,2,:,:)-coszfrr1p(:,1,:,:));
pci = squeeze(coszfrr1p(:,2,:,:)./coszfrr1p(:,1,:,:));
pci(isinf(pci)) = NaN;
lpci = log2(pci);
lpci(isinf(lpci)) = NaN;
omi = (squeeze(coszfrr1p(:,2,:,:))-squeeze(coszfrr1p(:,1,:,:)))./(squeeze(coszfrr1p(:,2,:,:))+squeeze(coszfrr1p(:,1,:,:)));

for i = 1:length(depth)
    for j = 1:5
        nmcoszfr(i,:,:,j) = coszfrr1(i,:,:,j)./max(coszfrr1(i,1,:,j),[],3);
        nmcoctrl(i,1,1,j) = cor1ctrlfr(i,1)./max(coszfrr1(i,1,:,j),[],3);
        nmcoctrl(i,2,1,j) = cor1ctrlfr(i,2)./max(coszfrr1(i,1,:,j),[],3);
        si(i,j) = (max(coszfrr1(i,1,:,j))-coszfrr1(i,1,3,j))/max(coszfrr1(i,1,:,j))-cor1ctrlfr(i,1);
        sinbl(i,j) = (max(coszfrr1(i,1,:,j))-coszfrr1(i,1,3,j))/max(coszfrr1(i,1,:,j));
    end
end
normcoszfr = cat(3,nmcoctrl,nmcoszfr);
normcoszfr(isinf(normcoszfr)) = NaN;

figure
plot(1,lpci(l23rs,3,1),'k.','markersize',15)
hold on
plot(2,lpci(l23rs,3,5),'k.','markersize',15)
axis([0,3,-3,2])
a = find(l23rs);
for i = 1:length(a)
    line([1,2],[lpci(a(i),3,1),lpci(a(i),3,5)],'color','k')
end
line([0,3],[0,0],'color','k')
set(gca,'xtick',[1,2])
set(gca,'xticklabel',{'low','high'})
ylabel('OMR')

data = normcoszfr;
cllt = l23rs;
figure
errorbar([0,cores(1).sizes],squeeze(nanmean(data(cllt,1,:,1))),squeeze(nanstd(data(cllt,1,:,1)))./sqrt(length(find(cllt))),'o-','color',[.9,.9,.9],'markerfacecolor',[.9,.9,.9])
hold on
errorbar([0,cores(1).sizes],squeeze(nanmean(data(cllt,1,:,2))),squeeze(nanstd(data(cllt,1,:,2)))./sqrt(length(find(cllt))),'o-','color',[.7,.7,.7],'markerfacecolor',[.7,.7,.7])
errorbar([0,cores(1).sizes],squeeze(nanmean(data(cllt,1,:,3))),squeeze(nanstd(data(cllt,1,:,3)))./sqrt(length(find(cllt))),'o-','color',[.5,.5,.5],'markerfacecolor',[.5,.5,.5])
errorbar([0,cores(1).sizes],squeeze(nanmean(data(cllt,1,:,4))),squeeze(nanstd(data(cllt,1,:,4)))./sqrt(length(find(cllt))),'o-','color',[.3,.3,.3],'markerfacecolor',[.3,.3,.3])
errorbar([0,cores(1).sizes],squeeze(nanmean(data(cllt,1,:,5))),squeeze(nanstd(data(cllt,1,:,5)))./sqrt(length(find(cllt))),'o-','color',[.1,.1,.1],'markerfacecolor',[.1,.1,.1])

cllt = l23rs;
%low contrast large vs small
figure
plot(squeeze(coszfrr1p(cllt,1,1,1)),squeeze(coszfrr1p(cllt,1,3,1)),'k.')
axis([0,16,0,16])
axis square
refline(1,0);
xlabel('firing rate small - low contrast')
ylabel('firing rate large - low contrast')

%high contrast large vs small
figure
plot(squeeze(coszfrr1p(cllt,1,1,5)),squeeze(coszfrr1p(cllt,1,3,5)),'k.')
axis([0,16,0,16])
axis square
refline(1,0);
xlabel('firing rate small - high contrast')
ylabel('firing rate large - high contrast')

%low contrast large light vs no light
figure
plot(squeeze(coszfrr1p(cllt,1,3,1)),squeeze(coszfrr1p(cllt,2,3,1)),'k.')
axis([0,14,0,14])
axis square
refline(1,0);
xlabel('firing rate large - low contrast - control')
ylabel('firing rate large - low contrast - light')

%high contrast large light vs no light
figure
plot(squeeze(coszfrr1p(cllt,1,3,5)),squeeze(coszfrr1p(cllt,2,3,5)),'k.')
axis([0,16,0,16])
axis square
refline(1,0);
xlabel('firing rate large - high contrast - control')
ylabel('firing rate large - high contrast - light')

a = find(cllt);
for i = 1:length(a)
    figure
    errorbar([0,cores(1).clevels],[cor1ctrlfr(a(i),1),squeeze(coszfrr1p(a(i),1,1,:))'],[cor1ctrlerr(a(i),1),squeeze(coszerrr1p(a(i),1,1,:))'],'o-','color',[.7,.7,.7],'markerfacecolor',[.7,.7,.7]);
    hold on
    errorbar([0,cores(1).clevels],[cor1ctrlfr(a(i),1),squeeze(coszfrr1p(a(i),1,2,:))'],[cor1ctrlerr(a(i),1),squeeze(coszerrr1p(a(i),1,2,:))'],'o-','color',[.4,.4,.4],'markerfacecolor',[.4,.4,.4]);
    errorbar([0,cores(1).clevels],[cor1ctrlfr(a(i),1),squeeze(coszfrr1p(a(i),1,3,:))'],[cor1ctrlerr(a(i),1),squeeze(coszerrr1p(a(i),1,3,:))'],'o-','color','k','markerfacecolor','k');
    legend('8 degrees','20 degrees','60 degrees')
    xlabel('contrast level')
    ylabel('firing rate')
end

for i = 1:length(a)
    figure
    errorbar([0,cores(1).sizes],[cor1ctrlfr(a(i),1),squeeze(coszfrr1(a(i),1,:,1))'],[cor1ctrlerr(a(i),1),squeeze(coszerrr1(a(i),1,:,1))'],'o-','color',[.9,.9,.9],'markerfacecolor',[.9,.9,.9]);
    hold on
    errorbar([0,cores(1).sizes],[cor1ctrlfr(a(i),1),squeeze(coszfrr1(a(i),1,:,2))'],[cor1ctrlerr(a(i),1),squeeze(coszerrr1(a(i),1,:,2))'],'o-','color',[.7,.7,.7],'markerfacecolor',[.7,.7,.7]);
    errorbar([0,cores(1).sizes],[cor1ctrlfr(a(i),1),squeeze(coszfrr1(a(i),1,:,3))'],[cor1ctrlerr(a(i),1),squeeze(coszerrr1(a(i),1,:,3))'],'o-','color',[.5,.5,.5],'markerfacecolor',[.5,.5,.5]);
    errorbar([0,cores(1).sizes],[cor1ctrlfr(a(i),1),squeeze(coszfrr1(a(i),1,:,4))'],[cor1ctrlerr(a(i),1),squeeze(coszerrr1(a(i),1,:,4))'],'o-','color',[.3,.3,.3],'markerfacecolor',[.3,.3,.3]);
    errorbar([0,cores(1).sizes],[cor1ctrlfr(a(i),1),squeeze(coszfrr1(a(i),1,:,5))'],[cor1ctrlerr(a(i),1),squeeze(coszerrr1(a(i),1,:,5))'],'o-','color','k','markerfacecolor','k');
    xlabel('stimulus size')
    ylabel('firing rate')
    legend('low','-','-','-','high')
end

for i = 1:length(a)
    clf;
    for j  = 1:5
        subplot(2,5,j)
        errorbar([0,cores(1).sizes],[cor1ctrlfr(a(i),1),squeeze(coszfrr1(a(i),1,:,j))'],[cor1ctrlerr(a(i),1),squeeze(coszerrr1(a(i),1,:,j))'],'ko-','markerfacecolor','k','linewidth',2);
        hold on
        errorbar([0,cores(1).sizes],[cor1ctrlfr(a(i),2),squeeze(coszfrr1(a(i),2,:,j))'],[cor1ctrlerr(a(i),2),squeeze(coszerrr1(a(i),2,:,j))'],'ro-','markerfacecolor','r','linewidth',2);
        ax = axis;
        axis([-10,62,ax(3),ax(4)])
        xlabel('size')
        if j == 1, title('low'); ylabel('firing rate'); elseif j == 5 title('high'); end
    end
    for j = 1:3
        subplot(2,3,j+3)
        errorbar([0,cores(1).clevels],[cor1ctrlfr(a(i),1),squeeze(coszfrr1(a(i),1,j,:))'],[cor1ctrlerr(a(i),1),squeeze(coszerrr1(a(i),1,j,:))'],'ko-','markerfacecolor','k','linewidth',2);
        hold on
        errorbar([0,cores(1).clevels],[cor1ctrlfr(a(i),2),squeeze(coszfrr1(a(i),2,j,:))'],[cor1ctrlerr(a(i),2),squeeze(coszerrr1(a(i),2,j,:))'],'ro-','markerfacecolor','r','linewidth',2);
        ax = axis;
        axis([-.1,.9,ax(3),ax(4)])
        xlabel('contrast')
        if j == 1, title('8'); ylabel('firing rate'); elseif j == 2 title('20'); elseif j == 3 title('60'); end
    end    
    if printyn
        figSize = [30 21];
        set(gcf,'paperunits','centimeters','papersize',figSize,'paperposition',[0 0 figSize])
        if i<10, printi = ['0', int2str(i)]; else printi = int2str(i); end
        print([coszprintpath ,  printi '__' cllname{i} '.pdf'],'-dpdf')
    end
end

for i = 1:length(depth)
    for j = 1:5
        nmcoszfr(i,:,:,j) = coszfrr1(i,:,:,j)./max(coszfrr1(i,1,:,j),[],3);
        nmcoctrl(i,:,1,j) = coctrlfr(i,:)./max(coszfrr1(i,1,:,j),[],3);
    end
end
normcoszfr = cat(3,nmcoctrl,nmcoszfr);
normcoszfr(isinf(normcoszfr)) = NaN;

data = normcoszfr;
cllt = l23rs;
figure
errorbar([0,cores(1).sizes],squeeze(nanmean(data(cllt,1,:,1))),squeeze(nanstd(data(cllt,1,:,1)))./sqrt(length(find(cllt))),'o-','color',[.9,.9,.9],'markerfacecolor',[.9,.9,.9])
hold on
errorbar([0,cores(1).sizes],squeeze(nanmean(data(cllt,1,:,2))),squeeze(nanstd(data(cllt,1,:,2)))./sqrt(length(find(cllt))),'o-','color',[.7,.7,.7],'markerfacecolor',[.7,.7,.7])
errorbar([0,cores(1).sizes],squeeze(nanmean(data(cllt,1,:,3))),squeeze(nanstd(data(cllt,1,:,3)))./sqrt(length(find(cllt))),'o-','color',[.5,.5,.5],'markerfacecolor',[.5,.5,.5])
errorbar([0,cores(1).sizes],squeeze(nanmean(data(cllt,1,:,4))),squeeze(nanstd(data(cllt,1,:,4)))./sqrt(length(find(cllt))),'o-','color',[.3,.3,.3],'markerfacecolor',[.3,.3,.3])
errorbar([0,cores(1).sizes],squeeze(nanmean(data(cllt,1,:,5))),squeeze(nanstd(data(cllt,1,:,5)))./sqrt(length(find(cllt))),'o-','color',[.1,.1,.1],'markerfacecolor',[.1,.1,.1])

cl = 5;
figure
errorbar([0,cores(1).sizes],squeeze(nanmean(data(cllt,1,:,cl))),squeeze(nanstd(data(cllt,1,:,cl)))./sqrt(length(find(cllt))),'ko-','markerfacecolor','k','linewidth',2)
hold on
errorbar([0,cores(1).sizes],squeeze(nanmean(data(cllt,2,:,cl))),squeeze(nanstd(data(cllt,2,:,cl)))./sqrt(length(find(cllt))),'ro-','markerfacecolor','r','linewidth',2)
xlabel('size')
ylabel('normalized FR')
title(['contrast: ' num2str(cores(1).clevels(cl))]);

for i = 1:length(depth)
    for j = 1:3
        nmtofccoszfr(i,:,j,:) = coszfrr1(i,:,j,:)./coszfrr1(i,1,j,5); % normalized to full contrast
        nmtomxcoszfr(i,:,j,:) = coszfrr1(i,:,j,:)./max(coszfrr1(i,1,j,:),[],4); % normalized to max response at that size
        nmxctrl(i,:,j) = coctrlfr(i,:)./max(coszfrr1(i,1,j,:)); % normalize control to max
        nfcctrl(i,:,j) = coctrlfr(i,:)./coszfrr1(i,1,j,5); % normalize control to max
    end
end
normmxcoszfr = cat(4,nmxctrl,nmtomxcoszfr);
normfccoszfr = cat(4,nfcctrl,nmtofccoszfr);
normmxcoszfr(isinf(normmxcoszfr)) = NaN;
normfccoszfr(isinf(normfccoszfr)) = NaN;

data = normmxcoszfr;
cllt = l23rs;
sz = 3;
figure
errorbar([0,cores(1).clevels],squeeze(nanmean(data(cllt,1,sz,:))),squeeze(nanstd(data(cllt,1,sz,:)))./sqrt(length(find(cllt))),'ko-','markerfacecolor','k','linewidth',2)
hold on
errorbar([0,cores(1).clevels],squeeze(nanmean(data(cllt,2,sz,:))),squeeze(nanstd(data(cllt,2,sz,:)))./sqrt(length(find(cllt))),'ro-','markerfacecolor','r','linewidth',2)
xlabel('contrast')
ylabel('normalized FR')
title(['size: ' num2str(cores(1).sizes(sz))]);

figure
errorbar([0,cores(1).clevels],squeeze(nanmean(data(l23rs,1,1,:))),squeeze(nanstd(data(l23rs,1,1,:))),'o-','color',[.7,.7,.7],'markerfacecolor',[.7,.7,.7])
hold on
errorbar([0,cores(1).clevels],squeeze(nanmean(data(l23rs,1,2,:))),squeeze(nanstd(data(l23rs,1,2,:))),'o-','color',[.4,.4,.4],'markerfacecolor',[.4,.4,.4])
errorbar([0,cores(1).clevels],squeeze(nanmean(data(l23rs,1,3,:))),squeeze(nanstd(data(l23rs,1,3,:))),'o-','color',[0,0,0],'markerfacecolor',[0,0,0])
errorbar([0,cores(1).clevels],squeeze(nanmean(data(l23rs,2,1,:))),squeeze(nanstd(data(l23rs,2,1,:))),'o-','color',[1,.7,.7],'markerfacecolor',[1,.7,.7])
errorbar([0,cores(1).clevels],squeeze(nanmean(data(l23rs,2,2,:))),squeeze(nanstd(data(l23rs,2,2,:))),'o-','color',[1,.4,.4],'markerfacecolor',[1,.4,.4])
errorbar([0,cores(1).clevels],squeeze(nanmean(data(l23rs,2,3,:))),squeeze(nanstd(data(l23rs,2,3,:))),'o-','color',[1,0,0],'markerfacecolor',[1,0,0])
legend('8 L0','20 L0','60 L0','8 L1','20 L1','60 L1')

data = lpci;
errorbar(cores(1).clevels,squeeze(nanmean(data(l23rs,1,:))),squeeze(nanstd(data(l23rs,1,:))),'o-','color',[.7,.7,.7],'markerfacecolor',[.7,.7,.7])
hold on
errorbar(cores(1).clevels,squeeze(nanmean(data(l23rs,2,:))),squeeze(nanstd(data(l23rs,2,:))),'o-','color',[.4,.4,.4],'markerfacecolor',[.4,.4,.4])
errorbar(cores(1).clevels,squeeze(nanmean(data(l23rs,3,:))),squeeze(nanstd(data(l23rs,3,:))),'o-','color','k','markerfacecolor','k')

% orientation tuning
% tuning width with optogenetics - pref size
figure
cond = l23rs' & pv' & szprefsfitrs(:,1)>.75 & szprefsfitrs(:,2)>.75;
[p,s] = signrank(squeeze(szprefsfit(cond,1,5)),squeeze(szprefsfit(cond,2,5)));
plot(squeeze(szprefsfit(cond,1,5)),squeeze(szprefsfit(cond,2,5)),'ko','markerfacecolor','k','markersize',4)
refline(1,0)
axis square
title(['p: ' num2str(p,3)])

% tuning width with optogenetics - pref clevel
figure
cond = l23rs' & pv' & coprefcfitrs(:,1)>.75 & coprefcfitrs(:,2)>.75;
[p,s] = signrank(squeeze(coprefcfit(cond,1,5)),squeeze(coprefcfit(cond,2,5)));
plot(squeeze(coprefcfit(cond,1,5)),squeeze(coprefcfit(cond,2,5)),'ko','markerfacecolor','k','markersize',4)
refline(1,0)
axis square
title(['p: ' num2str(p,3)])

% plot all the tuning curves at different contrasts
cond = l23rs;
hlp = find(cond);
for i = 1:length(hlp)
    clf
    hold on
    errorbar(cores(hlp(i)).oris,squeeze(cores(hlp(i)).condfr(1,:,1)),squeeze(cores(hlp(i)).conderr(1,:,1)),'bo','markerfacecolor','b')
    errorbar(cores(hlp(i)).oris,squeeze(cores(hlp(i)).condfr(1,:,2)),squeeze(cores(hlp(i)).conderr(1,:,2)),'co','markerfacecolor','c')
    errorbar(cores(hlp(i)).oris,squeeze(cores(hlp(i)).condfr(1,:,3)),squeeze(cores(hlp(i)).conderr(1,:,3)),'go','markerfacecolor','g')
    errorbar(cores(hlp(i)).oris,squeeze(cores(hlp(i)).condfr(1,:,4)),squeeze(cores(hlp(i)).conderr(1,:,4)),'mo','markerfacecolor','m')
    errorbar(cores(hlp(i)).oris,squeeze(cores(hlp(i)).condfr(1,:,5)),squeeze(cores(hlp(i)).conderr(1,:,5)),'ro','markerfacecolor','r')
    if coorifitrs(hlp(i),1,1)>.75
        plot(0:365,doublegauss(squeeze(cores(hlp(i)).gaussparams(1,1,:)),0:365),'b','linewidth',2);
    else
        plot(0:365,doublegauss(squeeze(cores(hlp(i)).gaussparams(1,1,:)),0:365),'b:');
    end
    if coorifitrs(hlp(i),1,2)>.75
        plot(0:365,doublegauss(squeeze(cores(hlp(i)).gaussparams(1,2,:)),0:365),'c','linewidth',2);
    else        
        plot(0:365,doublegauss(squeeze(cores(hlp(i)).gaussparams(1,2,:)),0:365),'c:');
    end
    if coorifitrs(hlp(i),1,3)>.75
        plot(0:365,doublegauss(squeeze(cores(hlp(i)).gaussparams(1,3,:)),0:365),'g','linewidth',2);
    else        
        plot(0:365,doublegauss(squeeze(cores(hlp(i)).gaussparams(1,3,:)),0:365),'g:');
    end
    if coorifitrs(hlp(i),1,4)>.75
        plot(0:365,doublegauss(squeeze(cores(hlp(i)).gaussparams(1,4,:)),0:365),'m','linewidth',2);
    else
        plot(0:365,doublegauss(squeeze(cores(hlp(i)).gaussparams(1,4,:)),0:365),'m:');
    end
    if coorifitrs(hlp(i),1,5)>.75
        plot(0:365,doublegauss(squeeze(cores(hlp(i)).gaussparams(1,5,:)),0:365),'r','linewidth',2);
    else
        plot(0:365,doublegauss(squeeze(cores(hlp(i)).gaussparams(1,5,:)),0:365),'r:');
    end
    line([0,365],[cores(hlp(i)).cresp(1,1),cores(hlp(i)).cresp(1,1)]);
    title(int2str(hlp(i)))
    pause
end

% plot all the tuning curves at different sizes
cond = l23rs;
hlp = find(cond);
for i = 1:length(hlp)
    clf
    hold on
    errorbar(szres(hlp(i)).oris,squeeze(szres(hlp(i)).condfr(1,:,1)),squeeze(szres(hlp(i)).conderr(1,:,1)),'bo','markerfacecolor','b')
    errorbar(szres(hlp(i)).oris,squeeze(szres(hlp(i)).condfr(1,:,2)),squeeze(szres(hlp(i)).conderr(1,:,2)),'co','markerfacecolor','c')
    errorbar(szres(hlp(i)).oris,squeeze(szres(hlp(i)).condfr(1,:,3)),squeeze(szres(hlp(i)).conderr(1,:,3)),'go','markerfacecolor','g')
    errorbar(szres(hlp(i)).oris,squeeze(szres(hlp(i)).condfr(1,:,4)),squeeze(szres(hlp(i)).conderr(1,:,4)),'mo','markerfacecolor','m')
    errorbar(szres(hlp(i)).oris,squeeze(szres(hlp(i)).condfr(1,:,5)),squeeze(szres(hlp(i)).conderr(1,:,5)),'ro','markerfacecolor','r')
    if szorifitrs(hlp(i),1,1)>.75
        plot(0:365,doublegauss(squeeze(szres(hlp(i)).gaussparams(1,1,:)),0:365),'b','linewidth',2);
    else
        plot(0:365,doublegauss(squeeze(szres(hlp(i)).gaussparams(1,1,:)),0:365),'b:');
    end
    if szorifitrs(hlp(i),1,2)>.75
        plot(0:365,doublegauss(squeeze(szres(hlp(i)).gaussparams(1,2,:)),0:365),'c','linewidth',2);
    else        
        plot(0:365,doublegauss(squeeze(szres(hlp(i)).gaussparams(1,2,:)),0:365),'c:');
    end
    if szorifitrs(hlp(i),1,3)>.75
        plot(0:365,doublegauss(squeeze(szres(hlp(i)).gaussparams(1,3,:)),0:365),'g','linewidth',2);
    else        
        plot(0:365,doublegauss(squeeze(szres(hlp(i)).gaussparams(1,3,:)),0:365),'g:');
    end
    if szorifitrs(hlp(i),1,4)>.75
        plot(0:365,doublegauss(squeeze(szres(hlp(i)).gaussparams(1,4,:)),0:365),'m','linewidth',2);
    else
        plot(0:365,doublegauss(squeeze(szres(hlp(i)).gaussparams(1,4,:)),0:365),'m:');
    end
    if szorifitrs(hlp(i),1,5)>.75
        plot(0:365,doublegauss(squeeze(szres(hlp(i)).gaussparams(1,5,:)),0:365),'r','linewidth',2);
    else
        plot(0:365,doublegauss(squeeze(szres(hlp(i)).gaussparams(1,5,:)),0:365),'r:');
    end
    line([0,365],[szres(hlp(i)).controlfr(1),szres(hlp(i)).controlfr(1)]);
    title(int2str(hlp(i)))
    pause
end

% plot the tuning curves light no light
cond = l23rs'&szorifitrs(:,2,3)>.75&szorifitrs(:,1,3)>.75;
hlp = find(cond);
for i = 1:length(hlp)
    clf
    hold on
    errorbar(szres(hlp(i)).oris,squeeze(szres(hlp(i)).condfr(1,:,3)),squeeze(szres(hlp(i)).conderr(1,:,3)),'ko','markerfacecolor','k')
    errorbar(szres(hlp(i)).oris,squeeze(szres(hlp(i)).condfr(2,:,3)),squeeze(szres(hlp(i)).conderr(2,:,3)),'ro','markerfacecolor','r')
    if szorifitrs(hlp(i),1,3)>.75
        plot(0:365,doublegauss(squeeze(szres(hlp(i)).gaussparams(1,3,:)),0:365),'k','linewidth',2);
    else
        plot(0:365,doublegauss(squeeze(szres(hlp(i)).gaussparams(1,3,:)),0:365),'k:');
    end
    if szorifitrs(hlp(i),1,3)>.75
        plot(0:365,doublegauss(squeeze(szres(hlp(i)).gaussparams(2,3,:)),0:365),'r','linewidth',2);
    else
        plot(0:365,doublegauss(squeeze(szres(hlp(i)).gaussparams(2,3,:)),0:365),'r:');
    end
    line([0,365],[szres(hlp(i)).controlfr(1),szres(hlp(i)).controlfr(1)]);
    title(int2str(hlp(i)))
    pause
end

bars = [squeeze(nanmean(cofr(l23rs,1,:)))',squeeze(nanmean(szfr(l23rs,1,:)))',squeeze(nanmean(mvfr(l23rs,1,:)))'];
errbars = [(squeeze(nanstd(cofr(l23rs,1,:)))./sqrt(sum(l23rs)))',(squeeze(nanstd(szfr(l23rs,1,:)))./sqrt(sum(l23rs)))',(squeeze(nanstd(mvfr(l23rs,1,:)))./sqrt(sum(l23rs)))'];
figure
barweb(bars,errbars,[],[],'firing rates','conditions','firing rate (Hz)',[],[],{'low','.','.','.','high','small','.','.','.','large','full screen','aperture'})

bars = [squeeze(nanmean(cosparseness(l23rs,1,:)))',squeeze(nanmean(szsparseness(l23rs,1,:)))',squeeze(nanmean(mvsparseness(l23rs,1,:)))'];
errbars = [(squeeze(nanstd(cosparseness(l23rs,1,:)))./sqrt(sum(l23rs)))',(squeeze(nanstd(szsparseness(l23rs,1,:)))./sqrt(sum(l23rs)))',(squeeze(nanstd(mvsparseness(l23rs,1,:)))./sqrt(sum(l23rs)))'];
figure
barweb(bars,errbars,[],[],'sparseness','conditions','sparseness')

bars = [squeeze(nanmean(corely(l23rs,1,:)))',squeeze(nanmean(szrely(l23rs,1,:)))',squeeze(nanmean(mvrely(l23rs,1,:)))'];
errbars = [(squeeze(nanstd(corely(l23rs,1,:)))./sqrt(sum(l23rs)))',(squeeze(nanstd(szrely(l23rs,1,:)))./sqrt(sum(l23rs)))',(squeeze(nanstd(mvrely(l23rs,1,:)))./sqrt(sum(l23rs)))'];
figure
barweb(bars,errbars,[],[],'reliability','conditions','reliability')

% firing rate movie FS vs large grating
rsp = signrank(mvfr(l23rs,1,1),szfr(l23rs,1,5));
fsp = signrank(mvfr(l23fs,1,1),szfr(l23fs,1,5));
figure
loglog(mvfr(l23rs,1,1),szfr(l23rs,1,5),'ko','markerfacecolor','k','markersize',5)
hold on
loglog(mvfr(l23fs,1,1),szfr(l23fs,1,5),'go','markerfacecolor','g','markersize',5)
loglog(nanmean(mvfr(l23rs,1,1)),nanmean(szfr(l23rs,1,5)),'k+','markersize',10);
loglog(nanmean(mvfr(l23fs,1,1)),nanmean(szfr(l23fs,1,5)),'g+','markersize',10);
% axis([0,1,0,1])
axis square
refline(1,0);
xlabel('movie fr full screen')
ylabel('grating FR large')
title(['firing rate movie full vs large grating RS: ' num2str(rsp) ' FS: ' num2str(fsp)])

% firing rate movie small vs small grating
rsp = signrank(mvfr(l23rs&apok,1,2),szfr(l23rs&apok,1,2));
fsp = signrank(mvfr(l23fs&apok,1,2),szfr(l23fs&apok,1,2));
figure
loglog(mvfr(l23rs&apok,1,2),szfr(l23rs&apok,1,2),'ko','markerfacecolor','k','markersize',5)
hold on
loglog(mvfr(l23fs&apok,1,2),szfr(l23fs&apok,1,2),'go','markerfacecolor','g','markersize',5)
loglog(nanmean(mvfr(l23rs&apok,1,2)),nanmean(szfr(l23rs&apok,1,2)),'k+','markersize',10);
loglog(nanmean(mvfr(l23fs&apok,1,2)),nanmean(szfr(l23fs&apok,1,2)),'g+','markersize',10);
% axis([0,1,0,1])
axis square
refline(1,0);
xlabel('movie fr aperture')
ylabel('grating FR small')
title(['firing rate movie aperture vs small grating RS: ' num2str(rsp) ' FS: ' num2str(fsp)])

% sparseness movie vs grating LARGE
rsp = signrank(mvsparseness(l23rs,1,1),szsparseness(l23rs,1,5));
fsp = signrank(mvsparseness(l23fs,1,1),szsparseness(l23fs,1,5));
figure
plot(mvsparseness(l23rs,1,1),szsparseness(l23rs,1,5),'ko','markerfacecolor','k','markersize',5)
hold on
plot(mvsparseness(l23fs,1,1),szsparseness(l23fs,1,5),'go','markerfacecolor','g','markersize',5)
plot(nanmean(mvsparseness(l23rs,1,1)),nanmean(szsparseness(l23rs,1,5)),'k+','markersize',10);
plot(nanmean(mvsparseness(l23fs,1,1)),nanmean(szsparseness(l23fs,1,5)),'g+','markersize',10);
axis([0,1,0,1])
axis square
refline(1,0);
xlabel('sparseness movie full screen')
ylabel('sparseness large grating')
title(['sparseness movie vs large grating RS: ' num2str(rsp) ' FS: ' num2str(fsp)])

% sparseness movie vs grating small
rsp = signrank(mvsparseness(l23rs&apok,1,2),szsparseness(l23rs&apok,1,2));
fsp = signrank(mvsparseness(l23fs&apok,1,2),szsparseness(l23fs&apok,1,2));
figure
plot(mvsparseness(l23rs&apok,1,2),szsparseness(l23rs&apok,1,2),'ko','markerfacecolor','k','markersize',5)
hold on
plot(mvsparseness(l23fs&apok,1,2),szsparseness(l23fs&apok,1,2),'go','markerfacecolor','g','markersize',5)
plot(nanmean(mvsparseness(l23rs&apok,1,2)),nanmean(szsparseness(l23rs&apok,1,2)),'k+','markersize',10);
plot(nanmean(mvsparseness(l23fs&apok,1,2)),nanmean(szsparseness(l23fs&apok,1,2)),'g+','markersize',10);
axis([0,1,0,1])
axis square
refline(1,0);
xlabel('sparseness movie aperture')
ylabel('sparseness small grating')
title(['sparseness movie vs small grating RS: ' num2str(rsp) ' FS: ' num2str(fsp)])

% reliability movie vs grating large
rsp = signrank(mvrely(l23rs&mvrelyok(:,1,1)'&szrelyok(:,1,5)',1,1),szrely(l23rs&mvrelyok(:,1,1)'&szrelyok(:,1,5)',1,5));
fsp = signrank(mvrely(l23fs&mvrelyok(:,1,1)'&szrelyok(:,1,5)',1,1),szrely(l23fs&mvrelyok(:,1,1)'&szrelyok(:,1,5)',1,5));
figure
plot(mvrely(l23rs&mvrelyok(:,1,1)'&szrelyok(:,1,5)',1,1),szrely(l23rs&mvrelyok(:,1,1)'&szrelyok(:,1,5)',1,5),'ko','markerfacecolor','k','markersize',5)
hold on
plot(mvrely(l23fs&mvrelyok(:,1,1)'&szrelyok(:,1,5)',1,1),szrely(l23fs&mvrelyok(:,1,1)'&szrelyok(:,1,5)',1,5),'go','markerfacecolor','g','markersize',5)
plot(nanmean(mvrely(l23rs&mvrelyok(:,1,1)'&szrelyok(:,1,5)',1,1)),nanmean(szrely(l23rs&mvrelyok(:,1,1)'&szrelyok(:,1,5)',1,5)),'k+','markersize',10);
plot(nanmean(mvrely(l23fs&mvrelyok(:,1,1)'&szrelyok(:,1,5)',1,1)),nanmean(szrely(l23fs&mvrelyok(:,1,1)'&szrelyok(:,1,5)',1,5)),'g+','markersize',10);
axis([0,1,0,1])
axis square
refline(1,0);
xlabel('reliability movie full screen')
ylabel('reliability large grating')
title(['reliability movie vs large grating RS: ' num2str(rsp) ' FS: ' num2str(fsp)])

% reliability movie vs grating small
rsp = signrank(mvrely(l23rs&mvrelyok(:,1,2)'&szrelyok(:,1,2)'&apok,1,2),szrely(l23rs&mvrelyok(:,1,2)'&szrelyok(:,1,2)'&apok,1,2));
fsp = signrank(mvrely(l23fs&mvrelyok(:,1,2)'&szrelyok(:,1,2)'&apok,1,2),szrely(l23fs&mvrelyok(:,1,2)'&szrelyok(:,1,2)'&apok,1,2));
figure
plot(mvrely(l23rs&mvrelyok(:,1,2)'&szrelyok(:,1,2)'&apok,1,2),szrely(l23rs&mvrelyok(:,1,2)'&szrelyok(:,1,2)'&apok,1,2),'ko','markerfacecolor','k','markersize',5)
hold on
plot(mvrely(l23fs&mvrelyok(:,1,2)'&szrelyok(:,1,2)'&apok,1,2),szrely(l23fs&mvrelyok(:,1,2)'&szrelyok(:,1,2)'&apok,1,2),'go','markerfacecolor','g','markersize',5)
plot(nanmean(mvrely(l23rs&mvrelyok(:,1,2)'&szrelyok(:,1,2)'&apok,1,2)),nanmean(szrely(l23rs&mvrelyok(:,1,2)'&szrelyok(:,1,2)'&apok,1,2)),'k+','markersize',10);
plot(nanmean(mvrely(l23fs&mvrelyok(:,1,2)'&szrelyok(:,1,2)'&apok,1,2)),nanmean(szrely(l23fs&mvrelyok(:,1,2)'&szrelyok(:,1,2)'&apok,1,2)),'g+','markersize',10);
axis([0,1,0,1])
axis square
refline(1,0);
xlabel('reliability movie aperture')
ylabel('reliability small grating')
title(['reliability movie vs small grating RS: ' num2str(rsp) ' FS: ' num2str(fsp)])

%movie jitter
rsp = signrank(mvjitter(l23rs,1,1),mvjitter(l23rs,1,2));
fsp = signrank(mvjitter(l23fs,1,1),mvjitter(l23fs,1,2));
figure
plot(mvjitter(l23rs,1,2),mvjitter(l23rs,1,1),'ko','markerfacecolor','k','markersize',5)
hold on
plot(mvjitter(l23fs,1,2),mvjitter(l23fs,1,1),'go','markerfacecolor','g','markersize',5)
% axis([0,1,0,1])
axis square
refline(1,0);
xlabel('jitter sd small aperture')
ylabel('jitter sd full screen')
title(['jitter sd large vs small RS: ' num2str(rsp) ' FS: ' num2str(fsp)])

% sparseness differences PV and SOM
mvdeltasplg = mvsparsenesswin(:,2,1)-mvsparsenesswin(:,1,1);
mvdeltaspsm = mvsparsenesswin(:,2,2)-mvsparsenesswin(:,1,2);
grdeltasplg = szsparsenesswin(:,2,5)-szsparsenesswin(:,1,5);
grdeltaspsm = szsparsenesswin(:,2,2)-szsparsenesswin(:,1,2);

mvfrchlg = log10(mvfr(:,2,1)./mvfr(:,1,1)); % firing rate change large
mvfrchsm = log10(mvfr(:,2,2)./mvfr(:,1,2));
grfrchlg = log10(szfr(:,2,5)./szfr(:,1,5));
grfrchsm = log10(szfr(:,2,2)./szfr(:,1,2));
mvfrchlg(isinf(mvfrchlg)) = NaN;
mvfrchsm(isinf(mvfrchsm)) = NaN;
grfrchlg(isinf(grfrchlg)) = NaN;
grfrchsm(isinf(grfrchsm)) = NaN;

mvdeltafrlg = mvfr(:,2,1)-mvfr(:,1,1); % delta firing rate
mvdeltafrsm = mvfr(:,2,2)-mvfr(:,1,2);
grdeltafrlg = szfr(:,2,5)-szfr(:,1,5);
grdeltafrsm = szfr(:,2,2)-szfr(:,1,2);

mvfrratiolg = mvfr(:,2,1)./mvfr(:,1,1); % unloged ratio
mvfrratiosm = mvfr(:,2,2)./mvfr(:,1,2);
grfrratiolg = szfr(:,2,5)./szfr(:,1,5);
grfrratiosm = szfr(:,2,2)./szfr(:,1,2);
mvfrratiolg(isinf(mvfrratiolg)) = NaN;
mvfrratiosm(isinf(mvfrratiosm)) = NaN;
grfrratiolg(isinf(grfrratiolg)) = NaN;
grfrratiosm(isinf(grfrratiosm)) = NaN;

% firing change movie large
[ypv,xpv] = hist(mvdeltafrlg(l23rs&pv),-5:1:10);
[ysom,xsom]= hist(mvdeltafrlg(l23rs&som),-5:1:10);
[p,s] = ranksum(mvdeltafrlg(l23rs&som),mvdeltafrlg(l23rs&pv));
figure
stairs(xpv,ypv,'g')
hold on
stairs(xsom,ysom,'r')
legend('PV','SOM')
xlabel('delta firing movie large')
ylabel('count')
title(['p: ' num2str(p,3)]);

% sparseness change movie large
[ypv,xpv] = hist(mvdeltasplg(l23rs&pv),-.2:.01:.2);
[ysom,xsom]= hist(mvdeltasplg(l23rs&som),-.2:.01:.2);
[p,s] = ranksum(mvdeltasplg(l23rs&som),mvdeltasplg(l23rs&pv));
figure
stairs(xpv,ypv,'g')
hold on
stairs(xsom,ysom,'r')
legend('PV','SOM')
xlabel('delta sparseness movie large')
ylabel('count')
title(['p: ' num2str(p,3)]);

% firing change movie small
[ypv,xpv] = hist(mvdeltafrsm(l23rs&pv),-5:1:10);
[ysom,xsom]= hist(mvdeltafrsm(l23rs&som),-5:1:10);
[p,s] = ranksum(mvdeltafrsm(l23rs&som),mvdeltafrsm(l23rs&pv));
figure
stairs(xpv,ypv,'g')
hold on
stairs(xsom,ysom,'r')
legend('PV','SOM')
xlabel('delta firing movie large')
ylabel('count')
title(['p: ' num2str(p,3)]);

%sparseness movie small
[ypv,xpv] = hist(mvdeltaspsm(l23rs&pv&apok),-.2:.01:.2);
[ysom,xsom]= hist(mvdeltaspsm(l23rs&som&apok),-.2:.01:.2);
[p,s] = ranksum(mvdeltaspsm(l23rs&som&apok),mvdeltaspsm(l23rs&pv&apok));
figure
stairs(xpv,ypv,'g')
hold on
stairs(xsom,ysom,'r')
legend('PV','SOM')
xlabel('delta sparseness movie small')
ylabel('count')
title(['p: ' num2str(p,3)]);

% firing change grating large
[ypv,xpv] = hist(grdeltafrlg(l23rs&pv),-5:1:10);
[ysom,xsom]= hist(grdeltafrlg(l23rs&som),-5:1:10);
[p,s] = ranksum(grdeltafrlg(l23rs&som),grdeltafrlg(l23rs&pv));
figure
stairs(xpv,ypv,'g')
hold on
stairs(xsom,ysom,'r')
legend('PV','SOM')
xlabel('delta firing movie large')
ylabel('count')
title(['p: ' num2str(p,3)]);

% sparseness grating large
[ypv,xpv] = hist(grdeltasplg(l23rs&pv),-.2:.01:.2);
[ysom,xsom]= hist(grdeltasplg(l23rs&som),-.2:.01:.2);
[p,s] = ranksum(grdeltasplg(l23rs&som),grdeltasplg(l23rs&pv))
figure
stairs(xpv,ypv,'g')
hold on
stairs(xsom,ysom,'r')
legend('PV','SOM')
xlabel('delta sparseness grating large')
ylabel('count')
title(['p: ' num2str(p,3)]);

% firing change grating small
[ypv,xpv] = hist(grdeltafrsm(l23rs&pv),-5:1:10);
[ysom,xsom]= hist(grdeltafrsm(l23rs&som),-5:1:10);
[p,s] = ranksum(grdeltafrsm(l23rs&som),grdeltafrsm(l23rs&pv));
figure
stairs(xpv,ypv,'g')
hold on
stairs(xsom,ysom,'r')
legend('PV','SOM')
xlabel('delta firing movie large')
ylabel('count')
title(['p: ' num2str(p,3)]);

% sparseness grating small
[ypv,xpv] = hist(grdeltaspsm(l23rs&pv&apok),-.2:.01:.2);
[ysom,xsom]= hist(grdeltaspsm(l23rs&som&apok),-.2:.01:.2);
[p,s] = ranksum(grdeltaspsm(l23rs&som&apok),grdeltaspsm(l23rs&pv&apok))
figure
stairs(xpv,ypv,'g')
hold on
stairs(xsom,ysom,'r')
legend('PV','SOM')
xlabel('delta sparseness grating small')
ylabel('count')
title(['p: ' num2str(p,3)]);

% firing rate movie vs grating SOM
[ypv,xpv] = hist(mvdeltafrlg(l23rs&som),-5:1:10);
[ysom,xsom]= hist(grdeltafrlg(l23rs&som),-5:1:10);
[p,s] = ranksum(mvdeltafrlg(l23rs&som),grdeltafrlg(l23rs&som));
figure
stairs(xpv,ypv,'g')
hold on
stairs(xsom,ysom,'r')
legend('movie','grating')
xlabel('SOM: delta firing movie vs grating large')
ylabel('count')
title(['p: ' num2str(p,3)]);

% firing rate movie vs grating PV
[ypv,xpv] = hist(mvdeltafrlg(l23rs&pv),-5:1:10);
[ysom,xsom]= hist(grdeltafrlg(l23rs&pv),-5:1:10);
[p,s] = ranksum(mvdeltafrlg(l23rs&pv),grdeltafrlg(l23rs&pv));
figure
stairs(xpv,ypv,'g')
hold on
stairs(xsom,ysom,'r')
legend('movie','grating')
xlabel('PV: delta firing movie vs grating large')
ylabel('count')
title(['p: ' num2str(p,3)]);

bars = [nanmean(mvdeltafrlg(l23rs&som)); nanmean(grdeltafrlg(l23rs&som)); nanmean(mvdeltafrlg(l23rs&pv)); nanmean(grdeltafrlg(l23rs&pv))];
errbars = [nanstd(mvdeltafrlg(l23rs&som))./sqrt(length(find(~isnan(mvdeltafrlg(l23rs&som)))));...
    nanstd(grdeltafrlg(l23rs&som))./sqrt(length(find(~isnan(grdeltafrlg(l23rs&som)))));...
    nanstd(mvdeltafrlg(l23rs&pv))./sqrt(length(find(~isnan(mvdeltafrlg(l23rs&pv)))));...
    nanstd(grdeltafrlg(l23rs&pv))./sqrt(length(find(~isnan(grdeltafrlg(l23rs&pv)))))];
% barweb(bars, errbars, [], [],

anovavec = [mvdeltafrlg(l23rs&som); grdeltafrlg(l23rs&som); mvdeltafrlg(l23rs&pv); grdeltafrlg(l23rs&pv)];
hlp1 = ones(size(mvdeltafrlg(l23rs&som)));
hlp2 = ones(size(mvdeltafrlg(l23rs&pv)));
hlptype = [hlp1;hlp1;hlp2.*2;hlp2.*2]; % twos for pv
hlpstim = [hlp1;hlp1.*2;hlp2;hlp2.*2]; % twos for gratings
[p,table,stats] = anovan(anovavec,{hlptype hlpstim},'model','full','display','off');

bars = [nanmean(mvdeltafrlg(l23rs&som)); nanmean(grdeltafrlg(l23rs&som)); nanmean(mvdeltafrlg(l23rs&pv)); nanmean(grdeltafrlg(l23rs&pv))];
errbars = [nanstd(mvdeltafrlg(l23rs&som))./sqrt(length(find(~isnan(mvdeltafrlg(l23rs&som)))));...
    nanstd(grdeltafrlg(l23rs&som))./sqrt(length(find(~isnan(grdeltafrlg(l23rs&som)))));...
    nanstd(mvdeltafrlg(l23rs&pv))./sqrt(length(find(~isnan(mvdeltafrlg(l23rs&pv)))));...
    nanstd(grdeltafrlg(l23rs&pv))./sqrt(length(find(~isnan(grdeltafrlg(l23rs&pv)))))];

anovavec = [mvdeltafrlg(l23rs&som); grdeltafrlg(l23rs&som); mvdeltafrlg(l23rs&pv); grdeltafrlg(l23rs&pv)];
hlp1 = ones(size(mvdeltafrlg(l23rs&som)));
hlp2 = ones(size(mvdeltafrlg(l23rs&pv)));
hlptype = [hlp1;hlp1;hlp2.*2;hlp2.*2]; % twos for pv
hlpstim = [hlp1;hlp1.*2;hlp2;hlp2.*2]; % twos for gratings
[p,table,stats] = anovan(anovavec,{hlptype hlpstim},'model','full','display','off');

% correlation to population /  noise correlation stuff? 
ltstrt = find(szbta>500,1); ltend = find(szbta>1500,1)-1;
anno = 2;
inds = find(animalno == anno);
for i = 1:length(depth)
    for lg = 1:2
        for sz = 1:5
            normszres(i,lg,sz,:) = szbincondresp(i,lg,sz,:)./mean(szbincondresp(i,lg,sz,ltstrt:ltend),4);
        end
    end
end
for i = 1:length(inds)
    for j = 1:length(inds)
        for lg = 1:2
            for sz = 1:5
                cellcorr(i,j,lg,sz) = nancorr(squeeze(normszres(inds(i),lg,sz,ltstrt:ltend)),squeeze(normszres(inds(j),lg,sz,ltstrt:ltend)));
            end
        end
    end
end
for i = 1:length(inds)
    for j = 1:length(inds)
        for lg = 1:2
            tccorr(i,j,lg) = nancorr(szfr(inds(i),lg,:),szfr(inds(j),lg,:));
        end
    end
end

%big figure
for i = 1:length(depth)
    clf
    
    subplot(2,6,1)
    imagesc(rfres(i).xax,rfres(i).yax,rfres(i).on)
    hold on
    plot_orrf_absdeg(rfres(i).gaussfiton,1,'w',2)
    axis square
    axis xy
    plot_circle(mvres(i).position(1),mvres(i).position(2),max(unique(mvres(i).widths))/2,'k',2)
    title(['rs: ' num2str(onfitrs(i),2) '  on overlap: ' num2str(onoverlap(i),2)])
    
    subplot(2,6,2)
    imagesc(rfres(i).xax,rfres(i).yax,rfres(i).off)
    hold on
    plot_orrf_absdeg(rfres(i).gaussfitoff,1,'w',2)
    axis square
    axis xy
    plot_circle(mvres(i).position(1),mvres(i).position(2),max(unique(mvres(i).widths))/2,'k',2)
    title(['rs: ' num2str(offfitrs(i),2) '  off overlap: ' num2str(offoverlap(i),2)])
    
    subplot(2,6,3)
    imagesc(rfres(i).xax,rfres(i).yax,rfres(i).on-rfres(i).off)
    hold on
    plot_orrf_absdeg(rfres(i).gaussfiton,1,'w',2)
    plot_orrf_absdeg(rfres(i).gaussfitoff,1,'k',2)
    axis square
    axis xy
    mx = max(max(rfres(i).on-rfres(i).off)); mn = min(min(rfres(i).on-rfres(i).off)); am = max([abs(mx),abs(mn)]);
    caxis([-am,am]);
       
    if ~isnan(szres(i).depth)
        subplot(2,6,4)
        errorbar(szres(i).oris, squeeze(szres(i).condfr(1,:,szres(i).prefsize)),squeeze(szres(i).conderr(1,:,szres(i).prefsize)),'ko','markerfacecolor','k')
        hold on
        errorbar(szres(i).oris, squeeze(szres(i).condfr(2,:,szres(i).prefsize)),squeeze(szres(i).conderr(2,:,szres(i).prefsize)),'ro','markerfacecolor','r')
        plot(0:365,doublegauss(squeeze(szres(i).gaussparams(1,szres(i).prefsize,:)),0:365),'k','linewidth',2);
        plot(0:365,doublegauss(squeeze(szres(i).gaussparams(2,szres(i).prefsize,:)),0:365),'r','linewidth',2);
        ax = axis;
        axis([-10,340,ax(3),ax(4)]);
        title(['OSI: ' num2str(szres(i).condosi(1,szres(i).prefsize),2) '  L1: '  num2str(szres(i).condosi(2,szres(i).prefsize),2) '  Rs: ' num2str(szres(i).gaussrsquared(1,szres(i).prefsize),2) ' L1: ' num2str(szres(i).gaussrsquared(2,szres(i).prefsize),2)])
        
        subplot(2,6,5)
        errorbar([0,szres(i).sizes], [szres(i).controlfr(1); squeeze(nanmean(szres(i).condfr(1,:,:),2))], [szres(i).controlerr(1); squeeze(nanmean(szres(i).conderr(1,:,:),2))],'ko-','linewidth',2,'markerfacecolor','k')
        hold on    
        errorbar([0,szres(i).sizes], [szres(i).controlfr(2); squeeze(nanmean(szres(i).condfr(2,:,:),2))], [szres(i).controlerr(2); squeeze(nanmean(szres(i).conderr(2,:,:),2))],'ro-','linewidth',2,'markerfacecolor','r')
        if swidthms(i) >.38 str = 'RS'; elseif swidthms(i) <.36 str = 'FS'; else str = 'in betw.'; end
        title(['swidth: ' num2str(swidthms(i),2) '  ' str])
        ax = axis;
        axis([0,62,ax(3),ax(4)]);
    elseif ~isnan(cores(i).depth)
        subplot(2,6,4)
        errorbar(cores(i).oris,squeeze(cores(i).condfr(1,:,5)),squeeze(cores(i).conderr(1,:,5)),'ko','markerfacecolor','k')
        hold on
        errorbar(cores(i).oris,squeeze(cores(i).condfr(2,:,5)),squeeze(cores(i).conderr(2,:,5)),'ro','markerfacecolor','r')
        plot(0:365,doublegauss(squeeze(cores(i).gaussparams(1,5,:)),0:365),'k','linewidth',2);
        plot(0:365,doublegauss(squeeze(cores(i).gaussparams(2,5,:)),0:365),'r','linewidth',2);
        ax = axis;
        axis([-10,340,ax(3),ax(4)]);
        title(['OSI: ' num2str(cores(i).condosi(1,5),2) '  L1: '  num2str(cores(i).condosi(2,5),2) '  Rs: ' num2str(cores(i).gaussrsquared(1,5),2) ' L1: ' num2str(cores(i).gaussrsquared(2,5),2)])
    end

    subplot(2,6,6)
    errorbar(cores(i).xlevels, [cores(i).controlfr(1); squeeze(nanmean(cores(i).condfr(1,:,1,:),2))], [cores(i).controlerr(1); squeeze(nanmean(cores(i).conderr(1,:,1,:),2))],'ko','markerfacecolor','k')
    hold on
    errorbar(cores(i).xlevels, [cores(i).controlfr(2); squeeze(nanmean(cores(i).condfr(2,:,1,:),2))], [cores(i).controlerr(2); squeeze(nanmean(cores(i).conderr(2,:,1,:),2))],'ro','markerfacecolor','r')
    plot(0:.05:1,NakaRushton(cores(i).paramsl0(1,:),0:.05:1),'k','linewidth',2)
    plot(0:.05:1,NakaRushton(cores(i).paramsl1(1,:),0:.05:1),'r','linewidth',2)
    ax = axis;
    axis([-.05,1.05,ax(3),ax(4)]);
    if depth(i)<375 layerstr = 'L2/3'; elseif depth(i)<=550 layerstr = 'L4'; elseif depth(i)<=800  layerstr = 'L5'; else layerstr = 'L6'; end
    title(['depth: ' num2str(depth(i),3) '  ' layerstr])
    
    if ~isnan(szres(i).depth)
        subplot(2,3,4)
        rasterplotall(szres(i).allresp);
        title([szres(i).printname, '   rely: ' num2str(nanmean(szres(i).condrely(1,:,szres(i).prefsize)),2) ' L1: ' num2str(nanmean(szres(i).condrely(2,:,szres(i).prefsize)),2)]);
%         axis([-300,2300,0,400])
    end
    if ~isnan(cores(i).depth)
        subplot(2,3,5)    
        rasterplotall(squeeze(cores(i).allresp(:,:,1,:)));
        title([cores(i).printname, '   rely: ' num2str(nanmean(cores(i).condrely(1,:,5)),2) ' L1: ' num2str(nanmean(cores(i).condrely(2,:,5)),2)]);
%         axis([-300,2300,0,400])
    end
    if ~isnan(mvres(i).depth)
        subplot(2,3,6)
        rasterplotallmv(mvres(i).condallresp);
        title([mvres(i).printname, '   rely: ' num2str(mvres(i).condrely(1,1),2) ' L1: ' num2str(mvres(i).condrely(2,1),2)]);
%         axis([-300,3300,0,120])
    end
    if printyn
        figSize = [30 21];
        set(gcf,'paperunits','centimeters','papersize',figSize,'paperposition',[0 0 figSize])
        if i<10, printi = ['0', int2str(i)]; else printi = int2str(i); end
        print([printpath ,  printi '__' cllname{i} '.pdf'],'-dpdf')
    end
end


% RF figure
for i = 1:length(depth)
    clf
    
    subplot(2,3,1)
    imagesc(rfres(i).xax,rfres(i).yax,rfres(i).on)
    hold on
    plot_orrf_absdeg(rfres(i).gaussfiton,1,'w',2)
    axis square
    axis xy
    plot_circle(mvres(i).position(1),mvres(i).position(2),max(unique(mvres(i).widths))/2,'k',2)
    title(['rs: ' num2str(onfitrs(i),2) '  on overlap: ' num2str(onoverlap(i),2)])
    
    subplot(2,3,2)
    imagesc(rfres(i).xax,rfres(i).yax,rfres(i).off)
    hold on
    plot_orrf_absdeg(rfres(i).gaussfitoff,1,'w',2)
    axis square
    axis xy
    plot_circle(mvres(i).position(1),mvres(i).position(2),max(unique(mvres(i).widths))/2,'k',2)
    title(['rs: ' num2str(offfitrs(i),2) '  off overlap: ' num2str(offoverlap(i),2)])
    
    subplot(2,3,3)
    imagesc(rfres(i).xax,rfres(i).yax,rfres(i).on-rfres(i).off)
    hold on
    plot_orrf_absdeg(rfres(i).gaussfiton,1,'w',2)
    plot_orrf_absdeg(rfres(i).gaussfitoff,1,'k',2)
    axis square
    axis xy
    mx = max(max(rfres(i).on-rfres(i).off)); mn = min(min(rfres(i).on-rfres(i).off)); am = max([abs(mx),abs(mn)]);
    caxis([-am,am]);
    title(['on off overlap: ' num2str(onoffoverlap(i),2)])
    
    subplot(2,3,4)
    imagesc(rfres(i).xax,rfres(i).yax,rfres(i).lfpon)
    hold on
    plot_orrf_absdeg(rfres(i).lfpfiton,1,'w',2)
    axis square
    axis xy
    plot_circle(mvres(i).position(1),mvres(i).position(2),max(unique(mvres(i).widths))/2,'k',2)
    title(['rs: ' num2str(lfprson(i),2) '  on overlap: ' num2str(lfponoverlap(i),2)])
    mx = max(max(rfres(i).lfpoff-rfres(i).lfpon)); mn = min(min(rfres(i).lfpon-rfres(i).lfpoff)); am = max([abs(mx),abs(mn)]);
    caxis([-am,am]);
    
    subplot(2,3,5)
    imagesc(rfres(i).xax,rfres(i).yax,rfres(i).lfpoff)
    hold on
    plot_orrf_absdeg(rfres(i).lfpfitoff,1,'w',2)
    axis square
    axis xy
    plot_circle(mvres(i).position(1),mvres(i).position(2),max(unique(mvres(i).widths))/2,'k',2)
    title(['rs: ' num2str(lfprsoff(i),2) '  off overlap: ' num2str(lfpoffoverlap(i),2)])
    mx = max(max(rfres(i).lfpoff-rfres(i).lfpon)); mn = min(min(rfres(i).lfpon-rfres(i).lfpoff)); am = max([abs(mx),abs(mn)]);
    caxis([-am,am]);
    
    subplot(2,3,6)
    imagesc(rfres(i).xax,rfres(i).yax,rfres(i).lfpon-rfres(i).lfpoff)
    hold on
    plot_orrf_absdeg(rfres(i).lfpfiton,1,'w',2)
    plot_orrf_absdeg(rfres(i).lfpfitoff,1,'k',2)
    axis square
    axis xy
    mx = max(max(rfres(i).lfpoff-rfres(i).lfpon)); mn = min(min(rfres(i).lfpon-rfres(i).lfpoff)); am = max([abs(mx),abs(mn)]);
    caxis([-am,am]);

    if printyn
        figSize = [30 21];
        set(gcf,'paperunits','centimeters','papersize',figSize,'paperposition',[0 0 figSize])
        if i<10, printi = ['0', int2str(i)]; else printi = int2str(i); end
        print([rfprintpath ,  printi '__' cllname{i} '.pdf'],'-dpdf')
    end
end

end


function val=NakaRushton(p,x)
% parameters of Naka-Rushton function as in Disney et al., Neuron, 2007
% [R_max, contrast Exponent n,  50%firing-Contrast, spontaneous rate sFR]

val = p(4)+p(1)*((x.^p(2))./(x.^p(2)+p(3).^p(2)));
end

function val = doublegauss(p,x)
    %[C Rp Rn theta, sigma]
%     val = C + Rp* exp(-(dirdist(x-theta)^2)/(2*(sigma^2)) )+ Rn* exp(-(dirdist(x+pi-theta)^2)/(2*(sigma^2)) )
    val = p(1) + p(2)* exp(-(dirdiff(x-p(4)).^2)/(2*(p(5).^2)) ) + p(3)* exp(-(dirdiff(x+180-p(4)).^2)/(2*(p(5).^2)) );
end

