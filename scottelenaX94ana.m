load('C:\Users\Julia\work\data\others\elenascottX94\ScottElenaCREDOGX94.mat')

mouseno = [1    1    1    1    2    2    3    3    4    4    5    6]
fileid =  [1331 1332 1333 1334 1339 1340 1344 1345 1350 1351 1324 1349]

cll = 1;
for fi = 1:length(fileid)
    structname=strcat('ID',num2str(fileid(fi)));
    tetname = fieldnames(eval(structname));
    
    for tet = 1:length(tetname)
        thisunitname=strcat(structname,'.',tetname{tet},'.');
        nounit=size(eval(strcat(thisunitname,'SingleUnitTrials')),1);
        
        for i = 1:nounit
            
            for l  = 1:2
                for s = 1:2
                    if l == 1 & s == 1, cond = 1;  % no light , bar
                    elseif l == 1 & s == 2, cond = 2; % no light, no bar
                    elseif l == 2 & s == 1, cond = 3; % light, bar
                    elseif l == 2 & s == 2, cond = 4; % light, no bar
                    end
                    condfr(cll,l,s) = strcat(thisunitname,'SingleUnitDriveCounts'){i,cond}
            
            
end