function MUA_RFseparation

% SOM Halo population
animalids =    {'150825', '150831', '150902', '150909', '150915', '151022', '151023', '151027', '151109', '151110'};
blocks_c  =     [3,         3,       1,        1,        8,        2,        4,        6,        6,        5];
blocks_s  =     [4,         2,       2,        2,        10,       3,        3,        6,        8,        6];
animal    =     [1,         2,       3,        4,        5,        6,        7,        8,        9,        9];
centelectrodes =[[17,32];   [1,16];  [1,16];   [1,16];   [17,32];  [1,16];   [17,32];  [17,32];  [17,32];  [1,16]];
centdepth    =  [400,       400,     400,      400,      400,      400,      400,      500,      500,      500];
surrelectrodes =[[1,16];    [17,32]; [17,32];  [17,32];  [1,16];   [17,32];  [1,16];   [1,16];   [1,16];   [17,32]];
surrdepth    =  [400,       400,     400,      400,      400,      400,      400,      500,      500,      500];
penangle =      [25,        25,      25,       25,       25,       25,       25,       25,       25,       25];
spacing =       [25,        25,      25,       25,       25,       25,       25,       25,       25,       25];
popfile = 'C:\Users\Julia\work\data\populations\SOM_Halo_later\RFs\MUA_population.mat';

centdepth = centdepth.*(cosd(penangle)*cosd(22));
surrdepth = surrdepth.*(cosd(penangle)*cosd(22));
spacing = spacing.*(cosd(penangle)*cosd(22));
evaldepth = 300;
for i = 1:length(centdepth)
    for j = 1:length(centelectrodes(i,1):centelectrodes(i,2))
        cdm(i,j) = centdepth(i)-((j-1)*spacing(i)); % depth matrix
        sdm(i,j) = surrdepth(i)-((j-1)*spacing(i)); % depth matrix
    end
    [c,cdi(i)] = min(abs(cdm(i,:)-evaldepth)); % get depth index, least distance to evaldepth
    [c,sdi(i)] = min(abs(sdm(i,:)-evaldepth)); % get depth index, least distance to evaldepth
end
sdi(5) = 4;

recalculate_pop = 1;
recalculate_muafile = 0;


if ~exist(popfile) || recalculate_pop

    for blck = 1:length(blocks_c)
        
        basepath = strcat('C:\Users\Julia\work\data\', animalids{blck}, '\');
        centfile = strcat(basepath, 'muaresult_', int2str(blocks_c(blck)), '_', int2str(centelectrodes(blck,1)), '-', int2str(centelectrodes(blck,2)), '.mat');
        surrfile = strcat(basepath, 'muaresult_', int2str(blocks_s(blck)), '_', int2str(surrelectrodes(blck,1)), '-', int2str(surrelectrodes(blck,2)), '.mat');
        oldcentfile = strcat(basepath, 'muaresult_', int2str(blocks_c(blck)), '_', int2str(centelectrodes(blck,1)), ':', int2str(centelectrodes(blck,2)), '.mat');
        oldsurrfile = strcat(basepath, 'muaresult_', int2str(blocks_s(blck)), '_', int2str(surrelectrodes(blck,1)), ':', int2str(surrelectrodes(blck,2)), '.mat');
                
        if ~exist(centfile) || recalculate_muafile
            if ~exist(oldcentfile)
                centresult = MUAdataprepare(basepath,animalids{blck},blocks(blck),centelectrodes(blck,1):centelectrodes(blck,2));
                save(centfile,'centresult')
            else
                load(oldcentfile);
                save(centfile, 'centresult');
            end
        else
            load(centfile);
        end        
        if ~exist(surrfile) || recalculate_muafile
            if ~exist(oldsurrfile)
                surrresult = MUAdataprepare(basepath,animalids{blck},blocks(blck),surrelectrodes(blck,1):surrelectrodes(blck,2));
                save(surrfile,'surrresult')
            else
                load(oldsurrfile);
                save(surrfile, 'surrresult');
            end
        else
            load(surrfile);
        end
        
        [on,off,gaussfiton,gaussfitoff,rson,rsoff,xax,yax] = get_rf_MUA(centfile);
        centon(:,:,blck) = on(:,:,cdi(blck)); centoff(:,:,blck) = off(:,:,cdi(blck));
        centonfit(blck) = gaussfiton(cdi(blck)); centofffit(blck) = gaussfitoff(cdi(blck));
        centonrse(blck) = rson(cdi(blck)); centoffrse(blck) = rsoff(cdi(blck));
        xax_c(:,blck) = xax; yax_c(:,blck) = yax;
        
        figure
        hold on
        for i = 1:length(gaussfiton)
            if(rson(i)>.75)
                plot_orrf_absdeg(gaussfiton(i),1,'r',1);
            end
            if(rsoff(i)>.75)
                plot_orrf_absdeg(gaussfitoff(i),1,'b',1);
            end
        end
        plot_orrf_absdeg(gaussfiton(cdi(blck)),1,'r',3);
        plot_orrf_absdeg(gaussfitoff(cdi(blck)),1,'b',3);
        
        [on,off,gaussfiton,gaussfitoff,rson,rsoff,xax,yax] = get_rf_MUA(surrfile);
        surron(:,:,blck) = on(:,:,sdi(blck)); surroff(:,:,blck) = off(:,:,sdi(blck));
        surronfit(blck) = gaussfiton(sdi(blck)); surrofffit(blck) = gaussfitoff(sdi(blck));
        surronrse(blck) = rson(sdi(blck)); surroffrse(blck) = rsoff(sdi(blck));
        xax_s(:,blck) = xax; yax_s(:,blck) = yax;
                
        for i = 1:length(gaussfiton)
            if(rson(i)>.75)
                plot_orrf_absdeg(gaussfiton(i),1,'r',1);
            end
            if(rson(i)>.75)
                plot_orrf_absdeg(gaussfitoff(i),1,'b',1);
            end
        end
        plot_orrf_absdeg(gaussfitoff(cdi(blck)),1,'b',3);
        plot_orrf_absdeg(gaussfiton(cdi(blck)),1,'r',3);
        
        if centonrse(blck)<0.75 & centoffrse(blck)<0.75
            disp(['shitty center fields here block ' int2str(block)]);
            centx(blck) = NaN; centy(blck) = NaN;
        elseif centonrse(blck)<0.75
            centx(blck) = centofffit(blck).shiftxcenterDeg; centy(blck) = centofffit(blck).shiftycenterDeg;
        elseif centoffrse(blck)<0.75
            centx(blck) = centonfit(blck).shiftxcenterDeg; centy(blck) = centonfit(blck).shiftycenterDeg;
        else
            centx(blck) = mean([centonfit(blck).shiftxcenterDeg,centofffit(blck).shiftxcenterDeg]);
            centy(blck) = mean([centonfit(blck).shiftycenterDeg,centofffit(blck).shiftycenterDeg]);
        end
            
        if surronrse(blck)<0.75 & surroffrse(blck)<0.75
            disp(['shitty surround fields here block ' int2str(blck)]);
            surrx(blck) = NaN; surry(blck) = NaN;
        elseif surronrse(blck)<0.75
            surrx(blck) = surrofffit(blck).shiftxcenterDeg; surry(blck) = surrofffit(blck).shiftycenterDeg;
        elseif surroffrse(blck)<0.75
            surrx(blck) = surronfit(blck).shiftxcenterDeg; surry(blck) = surronfit(blck).shiftycenterDeg;
        else
            surrx(blck) = mean([surronfit(blck).shiftxcenterDeg,surrofffit(blck).shiftxcenterDeg]);
            surry(blck) = mean([surronfit(blck).shiftycenterDeg,surrofffit(blck).shiftycenterDeg]);
        end
            
        separation(blck) = sqrt(abs(centx(blck)-surrx(blck))^2+abs(centy(blck)-surry(blck))^2);
    end
    save(popfile);
end

disp('');
