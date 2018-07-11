function MUA_population_freqpulses_oneelectrode


% % PV Ai32 V1 5342 - day 1
% animalids =     {'160926',  '160926',  '160926'};
% blocks    =     {'092616a', '092616b', '092616c'};
% animal    =     [1,          1,         1];
% depth    =      [150,        150,       150];
% penangle =      [0,          0,         0];
% intensity =     [0.1,        0.3,       1];
% popfile = 'C:\Users\Julia\work\data\populations\temp2.mat';
% 
% % PV Ai32 V1 5342 - day 2
% animalids =     {'160927',  '160927',  '160927'};
% blocks    =     {'092716a', '092716b', '092716c'};
% animal    =     [1,          1,         1];
% depth    =      [200,        200,       200];
% penangle =      [0,          0,         0];
% intensity =     [0.3,        1,         3];
% popfile = 'C:\Users\Julia\work\data\populations\temp2.mat';

% % PV Ai32 V1 - 5341
% animalids =     {'160928',  '160927_2', '160928',  '160927_2', '160927_2'};
% blocks    =     {'092816a', '092716d',  '092816b', '092716e',  '092716f'};
% animal    =     [1,          1,          1,         1,          1];
% depth    =      [150,        150,        150,       150,        150];
% penangle =      [0,          0,          0,         0,          0];
% intensity =     [0.1,        0.2,        0.3,       0.5,        1];
% popfile = 'C:\Users\Julia\work\data\populations\temp2.mat';

% % PV Ai32 S1 5341
% animalids =     {'160928',  '160927_2', '160927_2', '160927_2'};
% blocks    =     {'092816c', '092716g',  '092716h',  '092716i'};
% animal    =     [1,          1,          1,           1];
% depth    =      [150,        150,        150,         150];
% penangle =      [0,          0,          0,           0];
% intensity =     [0.1,        0.2,        0.5,         1];
% popfile = 'C:\Users\Julia\work\data\populations\temp2.mat';

% % SOM Ai32 V1 5556
% animalids =     {'160929',  '160929',  '160929',  '160929'};
% blocks    =     {'092916d','092916a', '092916b', '092916c'};
% animal    =     [1,          1,          1,          1];
% depth    =      [150,        150,        150,        150];
% penangle =      [0,          0,          0,          0];
% intensity =     [0.05,       0.1,        0.2,        0.3];
% popfile = 'C:\Users\Julia\work\data\populations\temp2.mat';

% % SOM Ai32 V1 5555
% animalids =     {'160930',  '160930'};
% blocks    =     {'093016a', '093016b'};
% animal    =     [1,          1];
% depth    =      [150,        150];
% penangle =      [0,          0];
% intensity =     [0.05,       0.1];
% popfile = 'C:\Users\Julia\work\data\populations\temp2.mat';

% % SOM Ai32 V1 5552
% animalids =     {'160930_2', '160930_2', '160930_2'};
% blocks    =     {'093016c',  '093016d',  '093016e'};
% animal    =     [1,          1,           1];
% depth    =      [150,        150,         150];
% penangle =      [0,          0,           0];
% intensity =     [0.05,       0.1,         0.2];
% popfile = 'C:\Users\Julia\work\data\populations\temp2.mat';

% % SOM Ai32 S1 5552
% animalids =     {'160930_2'};
% blocks    =     {'093016f'};
% animal    =     [1];
% depth    =      [150];
% penangle =      [0];
% intensity =     [0.1];
% popfile = 'C:\Users\Julia\work\data\populations\temp2.mat';

% % SOM Ai32 V1 5560
% animalids =     {'161003_2', '161003_2', '161003_2', '161003_2'};
% blocks    =     {'100316f',  '100316d',  '100316e',  '100316g'};
% animal    =      [1,          1,          1,          1];
% depth    =      [150,         150,        150,        150];
% penangle =      [0,           0,          0,          0];
% intensity =     [0.05,        0.1,        0.1,        0.2];
% popfile = 'C:\Users\Julia\work\data\populations\temp2.mat';

% % SOM Ai32 S1 5560
% animalids =     {'161003_2', '161003_2'};
% blocks    =     {'100316h',  '100316i'};
% animal    =      [1,          1];
% depth    =      [150,         150];
% penangle =      [0,           0];
% intensity =     [0.3,         0.5];
% popfile = 'C:\Users\Julia\work\data\populations\temp2.mat';

% % SOM Ai32 V1 5562
% animalids =     {'161004', '161004', '161004'};
% blocks    =     {'100416c','100416a','100416b'};
% animal    =      [1,        1,        1];
% depth    =      [200,       200,       200];
% penangle =      [0,         0,         0];
% intensity =     [0.05,      0.1,       0.2];
% popfile = 'C:\Users\Julia\work\data\populations\temp2.mat';

% % SOM Ai32 S1 5562
% animalids =     {'161004',  '161004'};
% blocks    =     {'100416e', '100416d'};
% animal    =      [1,         1];
% depth    =      [150,        150];
% penangle =      [0,          0];
% intensity =     [0.05,       0.1];
% popfile = 'C:\Users\Julia\work\data\populations\temp2.mat';

% % SOM Ai32 V1 5556
% animalids =     {'161004_2', '161004_2', '161004_2'};
% blocks    =     {'100416g',  '100416f',  '100416h'};
% animal    =     [1,           1,          1];
% depth    =      [150,         150,        150];
% penangle =      [0,           0,          0];
% intensity =     [0.05,        0.1,        0.2];
% popfile = 'C:\Users\Julia\work\data\populations\temp2.mat';

% % SOM Ai32 S1 5556
% animalids =     {'161004_2', '161004_2', '161004_2'};
% blocks    =     {'100416j',  '100416i',  '100416k'};
% animal    =     [1,           1,          1];
% depth    =      [150,         150,        150];
% penangle =      [0,           0,          0];
% intensity =     [0.05,        0.1,        0.2];
% popfile = 'C:\Users\Julia\work\data\populations\temp2.mat';

% % SOM Ai32 V1 5561
% animalids =     {'161005',  '161005'};
% blocks    =     {'100516b', '100516a'};
% animal    =     [1,          1];
% depth    =      [150,        150];
% penangle =      [0,          0];
% intensity =     [0.05,       0.1];
% popfile = 'C:\Users\Julia\work\data\populations\temp2.mat';

% % SOM Ai32 S1 5561
% animalids =     {'161005',  '161005'};
% blocks    =     {'100516d', '100516c'};
% animal    =     [1,          1];
% depth    =      [150,        150];
% penangle =      [0,          0];
% intensity =     [0.05,       0.1];
% popfile = 'C:\Users\Julia\work\data\populations\temp2.mat';

% % PV Ai32 V1 5608
% animalids =     {'161010',  '161010'};
% blocks    =     {'101016b', '101016a'};
% animal    =     [1,         1];
% depth    =      [150,       150];
% penangle =      [0,         0];
% intensity =     [0.05,      0.1];
% popfile = 'C:\Users\Julia\work\data\populations\temp2.mat';

% % PV Ai32 S1 5608
% animalids =     {'161010',  '161010'};
% blocks    =     {'101016d', '101016c'};
% animal    =     [1,         1];
% depth    =      [150,       150];
% penangle =      [0,         0];
% intensity =     [0.05,      0.1];
% popfile = 'C:\Users\Julia\work\data\populations\temp2.mat';

% % PV Ai32 V1 5607
% animalids =     {'161011',  '161011'};
% blocks    =     {'101116b', '101116a'};
% animal    =     [1,         1];
% depth    =      [150,       150];
% penangle =      [0,         0];
% intensity =     [0.05,      0.1];
% popfile = 'C:\Users\Julia\work\data\populations\temp2.mat';

% % PV Ai32 S1 5607
% animalids =     {'161011',  '161011'};
% blocks    =     {'101116d', '101116c'};
% animal    =     [1,         1];
% depth    =      [150,       150];
% penangle =      [0,         0];
% intensity =     [0.05,      0.1];
% popfile = 'C:\Users\Julia\work\data\populations\temp2.mat';

% % PV ChR2 6268 2/17/17 S1
% animalids =     {'170217',  '170217',  '170217',  '170217',  '170217'};
% blocks    =     {'021717d', '021717e', '021717c', '021717f', '021717b'};
% animal    =     [1,          1,         1,         1,         1];
% depth    =      [150,        150,       150,       150,       150];
% penangle =      [0,          0,         0,         0,         0];
% intensity =     [0.05,       0.1,       0.5,       1,         5];
% popfile = 'C:\Users\Julia\work\data\populations\temp2.mat';

% % PV ChR2 6272 2/18/17 V1
% animalids =     {'170218',  '170218',  '170218',  '170218'};
% blocks    =     {'021817k', '021817l', '021817m', '021817n'};
% animal    =     [1,          1,         1,         1];
% depth    =      [150,        150,       150,       150];
% penangle =      [0,          0,         0,         0];
% intensity =     [3,          3,         5,         5];
% pulsewidth =    [1,          1,         1,         3];
% iso         =   [1.5,        2.5,       1.5,       1.5];
% popfile = 'C:\Users\Julia\work\data\populations\temp2.mat';

% % PV ChR2 6269 2/21/17 V1 - pop TODO
% animalids =     {'170221'};
% blocks    =     {'022117i'};
% animal    =     [1];
% depth    =      [150];
% penangle =      [0];
% intensity =     [1];
% pulsewidth =    [3];
% iso         =   [2.5];
% popfile = 'C:\Users\Julia\work\data\populations\temp2.mat';

% % PV ChR2 6273 2/21/17 V1
% animalids =     {'170222'};
% blocks    =     {'022217b'};
% animal    =     [1];
% depth    =      [150];
% penangle =      [0];
% intensity =     [1];
% pulsewidth =    [3];
% iso         =   [2.5];
% popfile = 'C:\Users\Julia\work\data\populations\temp2.mat';

% % PV ChR2 6228 2/27/17 V1
% animalids =     {'170227',  '170227',  '170227',  '170227'};
% blocks    =     {'022717c', '022717d', '022717e', '022717f'};
% animal    =     [1,          1,         1,         1];
% depth    =      [150,        150,       150,       150];
% penangle =      [0,          0,         0,         0];
% intensity =     [0.1,        0.5,       0.3,       0.4];
% pulsewidth =    [1,          1,         1,         1];
% iso         =   [2.5,        2.5,       2.5,       2.5];
% popfile = 'C:\Users\Julia\work\data\populations\temp2.mat';

% % PV ChR2 6228 2/27/17 S1
% animalids =     {'170227',  '170227',  '170227'};
% blocks    =     {'022717g', '022717h', '022717i'};
% animal    =     [1,          1,         1];
% depth    =      [150,        150,       150];
% penangle =      [0,          0,         0];
% intensity =     [0.5,        0.1,       0.05];
% pulsewidth =    [1,          1,         1];
% iso         =   [2.5,        2.5,       2.5];
% popfile = 'C:\Users\Julia\work\data\populations\temp2.mat';

% % PV ChR2 6226 2/28/17 V1
% animalids =     {'170228',  '170228',  '170228'};
% blocks    =     {'022817f', '022817g', '022817h'};
% animal    =     [1,          1,         1];
% depth    =      [150,        150,       150];
% penangle =      [0,          0,         0];
% intensity =     [0.5,        0.1,       0.2];
% pulsewidth =    [1,          1,         1];
% iso         =   [2.5,        2.5,       2.5];
% popfile = 'C:\Users\Julia\work\data\populations\temp2.mat';

% % PV ChR2 6226 2/28/17 S1
% animalids =     {'170228',  '170228',  '170228',  '170228'};
% blocks    =     {'022817a', '022817b', '022817d', '022817e'};
% animal    =     [1,          1,         1,         1];
% depth    =      [150,        150,       150,       150];
% penangle =      [0,          0,         0,         0];
% intensity =     [0.5,        0.1,       0.1,       0.1];
% pulsewidth =    [1,          1,         1,         1];
% iso         =   [2.5,        2.5,       2.5,       2.5];
% popfile = 'C:\Users\Julia\work\data\populations\temp2.mat';

% % PV ChR2 6272 2/18/17 S1
% animalids =     {'170218',  '170218',  '170218',  '170218',  '170218',  '170218',  '170218',  '170218',  '170218',  '170218'};
% blocks    =     {'021817a', '021817b', '021817c', '021817d', '021817e', '021817f', '021817g', '021817h', '021817i', '021817j'};
% animal    =     [1,          1,         1,         1,         1,         1,         1,         1,         1,         1];
% depth    =      [150,        150,       150,       150,       150,       150,       150,       150,       150,       150];
% penangle =      [0,          0,         0,         0,          0,         0,         0,          0,         0,         0];
% intensity =     [0.1,        1,         3,         3,          3,         3,         5,          3,         3,         5];
% pulsewidth =    [3,          3,         3,         3,          3,         3,         3,          1,         1,         1];
% iso         =   [1.5,        1.5,       1.5,       0.5,        2.5,       1.5,       1.5,        1.5,       2.5,       2.5];
% popfile = 'C:\Users\Julia\work\data\populations\temp2.mat';

% % SOM Ai32 V1 population @ 0.1V
% animalids =     {'160929', '160930', '160930_2', '161003_2', '161004',  '161004_2', '161005'};
% blocks    =     {'092916a','093016b','093016d',  '100316e',  '100416a', '100416h',  '100516b'};
% animal    =     [1,         2,        3,          4,          5,         6,          7];
% depth    =      [150,       150,      150,        150,        150,       150,        150];
% penangle =      [0,         0,        0,          0,          0,         0,          0];
% intensity =     [0.1,       0.1,      0.1,        0.1         0.1,       0.2,        0.05];
% popfile = 'C:\Users\Julia\work\data\populations\temp2.mat';

% % SOM Ai32 S1 population @ 0.1V
% animalids =     {'160930_2', '161003_2', '161004',  '161004_2', '161005'};
% blocks    =     {'093016f',  '100316h',  '100416d', '100416j',  '100516c'};
% animal    =     [1,           2,          3,         4,          5];
% depth    =      [150,         150,        150,       150,        150];
% penangle =      [0,           0,          0,         0,          0];
% intensity =     [0.1,         0.3,        0.1,       0.1,        0.1];
% popfile = 'C:\Users\Julia\work\data\populations\temp2.mat';

% % PV Ai32 V1 population @ 0.1V (and 1 0.2V)
% animalids =     {'160926', '160927', '160927_2', '160928'};
% blocks    =     {'092616a','092716b','092716d',  '092816a'};
% animal    =     [1,         1,        1,          1];
% depth    =      [150,       150,      150,        150];
% penangle =      [0,         0,        0,          0];
% intensity =     [0.1,       0.1,      0.2,        0.1];
% popfile = 'C:\Users\Julia\work\data\populations\temp2.mat';

% % SOM Ai32 S1 ALL animals all intensities
% animalids =     {'160930_2', '161003_2', '161003_2', '161004', '161004',  '161004_2', '161004_2', '161004_2', '161005',  '161005'};
% blocks    =     {'093016f',  '100316h',  '100316i',  '100416e','100416d', '100416j',  '100416i',  '100416k',  '100516d', '100516c'};
% animal    =     [1,           2,         2,           3,         3,        4,          4,          4,          5,         5];
% depth    =      [150,         150,       150,         150,       150,      150,        150,        150,        150,       150];
% penangle =      [0,           0,         0,           0,         0,        0,          0,          0,          0,         0];
% intensity =     [0.1,         0.3,       0.5,         0.05,       0.1,     0.05,       0.1,        0.2,        0.05,      0.1];
% popfile = 'C:\Users\Julia\work\data\populations\SOMS1.mat';
% 
% PV Ai32 S1 ALL animals all intensities - 
animalids =     {'160927_2', '160927_2', '160927_2', '160928',  '161010',  '161010',  '161011',  '161011'};
blocks    =     {'092716g',  '092716h',  '092716i',  '092816c', '101016d', '101016c', '101116d', '101116c'};
animal    =     [1,           1,          1,          2,         3,         3,         4,         4];
depth    =      [150,         150,        150,        150,       150,       150,       150,       150];
penangle =      [0,           0,          0,          0,         0,         0,         0,         0];
intensity =     [0.2,         0.5,        1,          0.1,       0.05,      0.1,       0.05,      0.1];
popfile = 'C:\Users\Julia\work\data\populations\PVS1.mat';

% % paper
% % SOM Ai32 V1 ALL animals all intensities
% animalids =     {'160929',  '160929',  '160929',  '160929',  '160930',  '160930',  '160930_2', '160930_2', '160930_2', '161003_2', '161003_2', '161003_2', '161003_2', '161004', '161004', '161004',  '161004_2', '161004_2', '161004_2', '161005',  '161005'};
% blocks    =     {'092916d', '092916a', '092916b', '092916c', '093016a', '093016b', '093016c',  '093016d',  '093016e',  '100316f',  '100316d',  '100316e',  '100316g',  '100416c','100416a','100416b', '100416g',  '100416f',  '100416h',  '100516b', '100516a'};
% animal    =     [1,          1,         1,         1,         2,         2,         3,          3,          3,          4,          4,          4,          4,          5,        5,        5,         6,          6,          6,          7,         7];
% depth    =      [150,        150,       150,       150,       150,       150,       150,        150,        150,        150,        150,        150,        150,        150,      150,      150,       150,        150,        150,        150,       150];
% penangle =      [0,          0,         0,         0,         0,         0,         0,          0,          0,          0,          0,          0,          0,          0,        0,        0,         0,          0,          0,          0,         0];
% intensity =     [0.05,       0.1,       0.2,       0.3,       0.05,      0.1,       0.05,       0.1,        0.2,        0.05,       0.1,        0.1,        0.2,        0.05,     0.1,      0.2,       0.05,       0.1,        0.2,        0.05,      0.1];
% popfile = 'C:\Users\Julia\work\data\populations\SOMAi32awakePaper.mat';
% % popfile = 'C:\Users\Julia\work\data\populations\SOMChR2V1awake.mat';

% % PV Ai32 V1 ALL animals all intensities - onlty 2 animals but 3 locations
% animalids =     {'160926',  '160926',  '160926',  '160927',  '160927',  '160927',  '160927_2', '160927_2', '160927_2', '160928',  '160928',  '161010',  '161010',  '161011',  '161011'};
% blocks    =     {'092616a', '092616b', '092616c', '092716a', '092716b', '092716c', '092716d',  '092716e',  '092716f',  '092816a', '092816b', '101016b', '101016a', '101116b', '101116a'};
% animal    =     [1,          1,         1,         2,         2,         2,         3,         3,          3,           4,         4,         5,         5,         6,         6];
% depth    =      [150,        150,       150,       200,       200,       200,       150,       150,        150,         150,       150,       150,       150,       150,       150];
% penangle =      [0,          0,         0,         0,         0,         0,         0,         0,          0,           0,         0,         0,         0,         0,         0];
% intensity =     [0.1,        0.3,       1,         0.3,       1,         3,         0.2,       0.5,        1,           0.1,       0.3,       0.05,      0.1,       0.05,      0.1];
% % popfile = 'C:\Users\Julia\work\data\populations\PVChR2V1awake.mat';
% popfile = 'C:\Users\Julia\work\data\populations\PVAi32awakePaper.mat';

% % PV ChR2  anesth V1 population
% animalids =     {'170221',  '170221',  '170221',  '170222',  '170227',  '170227',  '170227',  '170227',  '170228',  '170228',  '170228'};
% blocks    =     {'022117a', '022117c', '022117d', '022217b', '022717c', '022717d', '022717e', '022717f', '022817f', '022817g', '022817h'};
% animal    =     [1,          1,         1,         2,         3,         3,         3,         3,         4,         4,         4];
% depth    =      [150,        150,       150,       150,       150,       150,       150,       150,       150,       150,       150];
% penangle =      [0,          0,         0,         0,         0,         0,         0,         0,         0,         0,         0];
% intensity =     [1,          0.3,       0.3,       3,         0.1,       0.5,       0.3,       0.4,       0.5,       0.1,       0.2];
% pulsewidth =    [3,          3,         1,         1,         1,         1,         1,         1,         1,         1,         1];
% iso         =   [2.5,        2.5,       2.5,       2.5,       2.5,       2.5,       2.5,       2.5,       2.5,       2.5,       2.5];
% popfile = 'C:\Users\Julia\work\data\populations\PVChR2V1.mat';

% % PV ChR2  anesth S1 population
% animalids =     {'170221',  '170221',  '170221',  '170221',  '170221',  '170227',  '170227',  '170227',  '170228',  '170228',  '170228',  '170228'};
% blocks    =     {'022117e', '022117f', '022117g', '022117h', '022117i', '022717g', '022717h', '022717i', '022817a', '022817b', '022817d', '022817e'};
% animal    =     [1,          1,         1,         1,         1,         2,         2,         2,         3,         3,         3,         3];
% depth    =      [150,        150,       150,       150,       150,       150,       150,       150,       150,       150,       150,       150];
% penangle =      [0,          0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0];
% intensity =     [0.1,        0.1,       0.1,       0.03,      0.05,      0.5,       0.1,       0.05,      0.5,       0.1,       0.1,       0.1];
% pulsewidth =    [1,          1,         1,         1,         1,         1,         1,         1,         1,         1,         1,         1];
% iso         =   [2.5,        2.5,       2.5,       2.5,       2.5,       2.5,       2.5,       2.5,       2.5,       2.5,       2.5,       2.5];
% popfile = 'C:\Users\Julia\work\data\populations\PVChR2S1.mat';

recalculate_pop = 0;
recalculate_muafile = 0;

% chronux parameters
params.tapers = [2,5]; params.Fs = 1000; params.err = [2, 0.05]; params.trialave = 1;


if ~exist(popfile) || recalculate_pop

    cll = 1;
    for blck = 1:length(blocks)
        
        basepath = strcat('C:\Users\Julia\work\data\', animalids{blck}, '\');
        file = strcat(basepath, blocks{blck}, '.mat');
        load(file);
        
        freqs = ExpStruct.testPulseFrequencies;
        shownfreqs = ExpStruct.stimOrder;
        validtrials = find(ExpStruct.autoPulseFrequency);
        
        sr = 20000;
        
        prestim = 500;
        poststim = 500;
        respwin = 501:3500; % after stimulus onset
              
%         % Todo tomorrow with wheel plugged in
%         for i = 1:length(msstamps)
%             speed(i,:) = result.runspeed(msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim);
%         end
%         
%         % figure out sufficiently high and nonvariable runspeed trials
%         meanspeed = mean(speed(:,respwin),2);
%         stdspeed = std(speed(:,respwin),1,2);
%         notstill = find(meanspeed>1);
%         okspeed = find(meanspeed>( mean(meanspeed(notstill))-(1.5*std(meanspeed(notstill))) )& meanspeed>1 );
%         okvar = find(stdspeed<( mean(stdspeed(notstill))+(1.5*std(stdspeed(notstill)))) & stdspeed>.5);
%         oktrials = intersect(okspeed,okvar);
%         nonoktrials = 1:size(speed,1); nonoktrials(oktrials) = [];
%         stilltrials = 1:size(speed,1); stilltrials(notstill) = [];
            
        L = length(respwin); nfft = 2^nextpow2(L);
        t = (0:L-1)*(1/1000); fftx = 1000/2*linspace(0,1,nfft/2+1);
        
      
        for i = 1:length(validtrials)
            lfpresp(i,:) = resample(sweeps{validtrials(i)}(:,1),1,sr/1000);
            lfpresp(i,:) = lfpresp(i,:)-mean(lfpresp(i,respwin));
            if find(lfpresp(i,respwin)<-0.7), 
                lfpresp(i,:) = nan(1,size(lfpresp,2)); 
            end
        end        
%         % remove trials with too much noise
%         allvar = std(lfpresp,[],2);
%         badones = find(allvar>median(allvar)+2*std(allvar));
%         lfpresp(badones,:) = nan(length(badones),size(lfpresp,2));
        % only then do the frequency transform
        for i = 1:length(validtrials)          
            [lfpspect(i,:),fax] = pmtm(lfpresp(i,respwin),3,[],1000);
            hlp = fft(lfpresp(i,respwin),nfft)/L;
            fftspect(i,:) = 2*abs(hlp(1:nfft/2+1));
            chronuxspect(i,:)=mtspectrumc(squeeze(lfpresp(i,respwin))',params);
        end
        
        for f = 1:length(freqs)
            
            thisinds = find(shownfreqs(validtrials) == f,30,'first');
            
            condn(blck,f) = length(thisinds);
            condlfpspect(blck,f,:) = nanmean(lfpspect(thisinds,:),1);
            condfftspect(blck,f,:) = nanmean(fftspect(thisinds,:),1);
            condchronspect(blck,f,:) = nanmean(chronuxspect(thisinds,:),1);
            condlfpresp(blck,f,:) = nanmean(lfpresp(thisinds,:),1);
            relspect(blck,f,:) = condlfpspect(blck,f,:)./condlfpspect(blck,1,:);
            relfftspect(blck,f,:) = condfftspect(blck,f,:)./condfftspect(blck,1,:);
            relchronspect(blck,f,:) = condchronspect(blck,f,:)./condchronspect(blck,1,:);
            
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
                meanpow(blck,f) = mean(condlfpspect(blck,f,1:500),3);
                powerratio(blck,f) = peakpow(blck,f)./meanpow(blck,f);
                otherratio(blck,f) = peakpow(blck,f)./surroundpow(blck,f);
                blratio(blck,f) = peakpow(blck,f)./blpeakpow(blck,f);
                
                fftpeakpow(blck,f) = condfftspect(blck,f,f_ind);
                fftblpeakpow(blck,f) = condfftspect(blck,1,f_ind);
                fftsurroundpow(blck,f) = mean([mean(condfftspect(blck,f,f_ind-7:f_ind-5),3),mean(condfftspect(blck,f,f_ind+5:f_ind+7),3)]);
                fftmeanpow(blck,f) = mean(condfftspect(blck,f,1:500),3);
                fftpowerratio(blck,f) = fftpeakpow(blck,f)./fftmeanpow(blck,f);
                fftotherratio(blck,f) = fftpeakpow(blck,f)./fftsurroundpow(blck,f);
                fftblratio(blck,f) = fftpeakpow(blck,f)./fftblpeakpow(blck,f);
                
                chrpeakpow(blck,f) = mean(condchronspect(blck,f,f_ind-2:f_ind+2),3);
                chrblpeakpow(blck,f) = mean(condchronspect(blck,1,f_ind-2:f_ind+2),3);
                chrblratio(blck,f) = chrpeakpow(blck,f)./chrblpeakpow(blck,f);
                
                relpeakpow(blck,f) = relspect(blck,f,f_ind);
                relmeanpow(blck,f) = mean(relspect(blck,f,1:500),3);
                relpowerratio(blck,f) = relpeakpow(blck,f)./relmeanpow(blck,f);
                
                relfftpeakpow(blck,f) = relfftspect(blck,f,f_ind);
                relfftmeanpow(blck,f) = mean(relfftspect(blck,f,1:500));
                relfftpowerratio(blck,f) = relfftpeakpow(blck,f)./relfftmeanpow(blck,f);
                
                relchrpeakpow(blck,f) = relchronspect(blck,f,f_ind);
            end
            
        end
    end
    save(popfile, '-v7.3');
else
    load(popfile);
end

blck = 8;
figure
semilogy(fax,squeeze(condlfpspect(blck,1,:)))
hold on
semilogy(fax,squeeze(condlfpspect(blck,2,:)),'c')
semilogy(fax,squeeze(condlfpspect(blck,3,:)),'g')
semilogy(fax,squeeze(condlfpspect(blck,4,:)),'y')
semilogy(fax,squeeze(condlfpspect(blck,5,:)),'m')
semilogy(fax,squeeze(condlfpspect(blck,6,:)),'r')
semilogy(fax,squeeze(condlfpspect(blck,7,:)),'k')
semilogy(fax,squeeze(condlfpspect(blck,8,:)),'b')
semilogy(fax,squeeze(condlfpspect(blck,9,:)),'c')
semilogy(fax,squeeze(condlfpspect(blck,10,:)),'g')
semilogy(fax,squeeze(condlfpspect(blck,11,:)),'y')
semilogy(fax,squeeze(condlfpspect(blck,12,:)),'m')
semilogy(fax,squeeze(condlfpspect(blck,1,:)),'linewidth',2)

blck = 1;
figure
semilogy(fax,squeeze(condfftspect(blck,1,:)))
hold on
semilogy(fax,squeeze(condfftspect(blck,2,:)),'c')
semilogy(fax,squeeze(condfftspect(blck,3,:)),'g')
semilogy(fax,squeeze(condfftspect(blck,4,:)),'y')
semilogy(fax,squeeze(condfftspect(blck,5,:)),'m')
semilogy(fax,squeeze(condfftspect(blck,6,:)),'r')
semilogy(fax,squeeze(condfftspect(blck,7,:)),'k')
semilogy(fax,squeeze(condfftspect(blck,8,:)),'b')
semilogy(fax,squeeze(condfftspect(blck,9,:)),'c')
semilogy(fax,squeeze(condfftspect(blck,10,:)),'g')
semilogy(fax,squeeze(condfftspect(blck,11,:)),'y')
semilogy(fax,squeeze(condfftspect(blck,12,:)),'m')
semilogy(fax,squeeze(condfftspect(blck,1,:)),'linewidth',2)

blck = 1;
figure
semilogy(fax,squeeze(condchronspect(blck,1,:)))
hold on
semilogy(fax,squeeze(condchronspect(blck,2,:)),'c')
semilogy(fax,squeeze(condchronspect(blck,3,:)),'g')
semilogy(fax,squeeze(condchronspect(blck,4,:)),'y')
semilogy(fax,squeeze(condchronspect(blck,5,:)),'m')
semilogy(fax,squeeze(condchronspect(blck,6,:)),'r')
semilogy(fax,squeeze(condchronspect(blck,7,:)),'k')
semilogy(fax,squeeze(condchronspect(blck,8,:)),'b')
semilogy(fax,squeeze(condchronspect(blck,9,:)),'c')
semilogy(fax,squeeze(condchronspect(blck,10,:)),'g')
semilogy(fax,squeeze(condchronspect(blck,11,:)),'y')
semilogy(fax,squeeze(condchronspect(blck,12,:)),'m')
semilogy(fax,squeeze(condchronspect(blck,1,:)),'linewidth',2)

% pvinds = [1,4,7,10]; % paper
% pvinds = [1,4,7,10,13,15]; % more data actual paper
pvinds = [4,6,8]; % PV Ai32 S1
% pvinds = [2,6,9,13];
% pvinds = [4,7,10,14]; % strongest
% pvinds = [2,4,7,10]; % PV ChR2 V1 anesth
% pvinds = [3,6,9]; % PV ChR2 S1 anesth
for i = 1:length(pvinds)
    fftnormcurve(i,:) = fftblratio(pvinds(i),:)./max(fftblratio(pvinds(i),:));
    normcurve(i,:) = blratio(pvinds(i),:)./max(blratio(pvinds(i),:));
    chrnormcurve(i,:) = chrblratio(pvinds(i),:)./max(chrblratio(pvinds(i),:));
    normpowercurve(i,:) = powerratio(pvinds(i),:)./max(powerratio(pvinds(i),:));
end
errorbar(freqs,mean(fftnormcurve,1),std(fftnormcurve,1,1)./sqrt(length(pvinds)),'g');
errorbar(freqs,mean(normcurve,1),std(normcurve,1,1)./sqrt(length(pvinds)),'g');
errorbar(freqs,mean(chrnormcurve,1),std(chrnormcurve,1,1)./sqrt(length(pvinds)),'g');

errorbar(freqs,mean(normcurve),std(normcurve)./sqrt(length(pvinds)),'k','linewidth',2)

[p,tbl,stats] = anova1(normcurve)
[p,tbl,stats] = kruskalwallis(normcurve)

% sominds = [1,4,7,10,13];
% sominds = [2,5,8,11,14];
% sominds = [3,6,8,12,15];  %strongest
sominds = [4,6,9,12,15,19,21]; % paper
% sominds = [2, 5, 7, 10]; % S1
for i = 1:length(sominds)
    fftnormcurve(i,:) = fftblratio(sominds(i),:)./max(fftblratio(sominds(i),:));
    normcurve(i,:) = blratio(sominds(i),:)./max(blratio(sominds(i),:));
    chrnormcurve(i,:) = chrblratio(sominds(i),:)./max(chrblratio(sominds(i),:));
end
errorbar(freqs,mean(fftnormcurve),std(fftnormcurve)./sqrt(length(sominds)),'r');
errorbar(freqs,mean(normcurve),std(normcurve)./sqrt(length(sominds)),'r');
errorbar(freqs,mean(chrnormcurve),std(chrnormcurve)./sqrt(length(sominds)),'r');

errorbar(freqs,nanmean(normcurve),nanstd(normcurve)./sqrt(length(sominds)),'k','linewidth',2)

anms = unique(animal);
for i = 1:length(anms)
    thisinds = find(animal == anms(i));
    tmp = find(intensity(thisinds) == max(intensity(thisinds)),1);
    maxinds(i) = thisinds(tmp);
    tmp = find(intensity(thisinds) == min(intensity(thisinds)),1);
    mininds(i) = thisinds(tmp);

    figure
    plot(freqs,blratio(thisinds,:),'o-');
    title(int2str(anms(i)));
    legend(num2str(intensity(thisinds)'))
    set(gca,'xtick',freqs)
end

% anesthetized vs awake
load('C:\Users\Julia\work\manuscripts\gamma\4th\cardindatas\awakepop.mat')
load('C:\Users\Julia\work\manuscripts\gamma\4th\cardindatas\anesthpop.mat')
load('C:\Users\Julia\work\manuscripts\gamma\4th\cardindatas\ans1.mat')
load('C:\Users\Julia\work\manuscripts\gamma\4th\cardindatas\anv1.mat')
freqs = [0,8,16,24,32,40,48,56,64,72,80,100];
awcurve(:,1) = [];
ancurve(:,1) = [];
awvec = awcurve(:);
anvec = ancurve(:);
anovavec = [awvec;anvec];
ga = [zeros(length(awvec),1);ones(length(anvec),1)];
hlp = repmat(freqs(2:end),size(awcurve,1),1);
fawvec = hlp(:);
hlp = repmat(freqs(2:end),size(ancurve,1),1);
fanvec = hlp(:);
gf = [fawvec;fanvec];
[p,table,stats] = anovan(anovavec,{gf,ga},'model','full')
multcompare(stats)
multcompare(stats,'dimension',2)

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