function lfp_population_ana

% % SOM contrast population
% animalids = {'140807', '140807', '140812', '140812', '140815', '140815', '140815', '140818', '141014', '141014', '141015', '141106', '141110b', '141111', '150123', '150130', '150406', '150407'};
% blocks    = [4,         7,        2,        6,        3,        7,        8,        4,        2,        9,        7,        6,        2,         2,        7,        4,        3,        3];
% animal =    [1          1         2         2         3         3         3         4,        5,        5,        6,        7,        8,         9,        10,       11,       12,       13];
% channels =  [32,        32,       32,       32,       32,       32,       32,       32,       32,       32,       32,       32,       32,        32,       32,       32,       32,       32];
% depths =    [500,       950,      400,      950,      400,      950,      1150,     950,      800,      1050,     950,      950,      950,       850,      900,      950,      900,      900];
% penangle = 80;
% popfile = 'C:\Users\Julia\work\data\populations\SOM_Halo\contrast\contrast_lfp_population.mat';
% printpath1 = 'C:\Users\Julia\work\data\populations\SOM_Halo\contrast\LFP\spectra\';
% printpath2 = 'C:\Users\Julia\work\data\populations\SOM_Halo\contrast\LFP\spectradiff\';
% printpath3 = 'C:\Users\Julia\work\data\populations\SOM_Halo\contrast\LFP\gamma\';

% % SOM size population
% animalids = {'140807', '140807', '140812', '140812', '140815', '140815', '140818', '141014', '141014', '141015', '141106', '141110b', '141111', '150123',  '150130', '150216', '150219', '150219', '150407'};
% blocks    = [5,         8,        3,        5,        4,        6,        5,        3,        10,       5,        3,        7,         4,        6,         5,        6,        6,        8,        5];
% animal =    [1          1         2         2         3         3         4,        5,        5,        6,        7,        8,         9,        10,        11,       13,       14,       14,       15];
% channels =  [32,        32,       32,       32,       32,       32,       32,       32,       32,       32,       32,       32,       32,        32,        32,       32,       32,       32,       32];
% depths =    [500,       950,      400,      950,      400,      950,      950,      800,      1050,     950,      950,      950,       850,      900,       950,      950,      500,      1050,     900];
% penangle = 80;
% popfile = 'C:\Users\Julia\work\data\populations\SOM_Halo\size\size_lfp_population.mat';
% printpath1 = 'C:\Users\Julia\work\data\populations\SOM_Halo\size\LFP\spectra\';
% printpath2 = 'C:\Users\Julia\work\data\populations\SOM_Halo\size\LFP\spectradiff\';
% printpath3 = 'C:\Users\Julia\work\data\populations\SOM_Halo\size\LFP\gamma\';

% SOM size population
animalids = {'140807', '140812', '140815', '140818', '141015', '141106', '141110b', '150130', '150216'};
blocks    = [8,         5,        6,        5,        5,        3,        7,         5,        6];
animal =    [1          2         3         4,        6,        7,        8,         11,       13];
channels =  [32,        32,       32,       32,       32,       32,       32,        32,       32];
depths =    [950,       950,      950,      950,      950,      950,      950,       950,      950];
penangle = 80;
popfile = 'C:\Users\Julia\work\data\populations\SOM_Halo\size\size_lfp_population_950.mat';
printpath1 = 'C:\Users\Julia\work\data\populations\SOM_Halo\size\LFP\spectra\';
printpath2 = 'C:\Users\Julia\work\data\populations\SOM_Halo\size\LFP\spectradiff\';
printpath3 = 'C:\Users\Julia\work\data\populations\SOM_Halo\size\LFP\gamma\';


printyn = 1;
dolightstuff = 1;
recalculate_csdmat = 0;
recalculate = 1;

if ~exist(popfile) || recalculate

    for blck = 1:length(blocks)

        basepath = ['C:\Users\Julia\work\data\' animalids{blck} '\'];
        
        if animal(blck)<10, panino = ['0', int2str(animal(blck))]; else panino = int2str(animal(blck)); end
        printname = [panino '_' animalids{blck} '_block_' int2str(blocks(blck))];

        if ~exist([basepath 'csdresult_' int2str(blocks(blck)) '.mat']) || recalculate_csdmat
            result = lfpdataprepare(basepath,animalids{blck},blocks(blck),channels(blck));
            save([basepath 'csdresult_' int2str(blocks(blck)) '.mat'], 'result')
        else
            load([basepath 'csdresult_' int2str(blocks(blck)) '.mat']);
        end
        lfp = result.lfp';        
        
        ellength = 25*(channels(blck)-1);
        highest = depths(blck)-ellength;
        depthax(blck,:) = depths(blck):-25:highest;
        depthax(blck,:) = depthax(blck,:).*(cosd(90-penangle)*cosd(22));
        
        if length(result.msstamps)~=length(result.light)
                disp('');
%                 result.msstamps([54,201,239,316]) = []; % for 140807 block 7
        %         result.msstamps(167) = [];
        %         result.msstamps(92) = [];
        %         result.msstamps([303]) = []; % for 150407 block 5
%                 save([basepath 'csdresult_' int2str(blocks(blck)) '.mat'], 'result')
%                 pause;
        end

        if strcmp(result.protocol, 'grating')
            prestim = 300;
            poststim = 700;
            result.stimduration = result.stimduration.*1000;
            timeax = -prestim+1:result.stimduration+poststim;
        else
            prestim = 0; poststim = 0; result.stimduration = round((1000/60)*result.FramesperStim); 
            timeax = 1:result.stimduration;
        end
        sr = 1000;
        gamma = eegfilt(lfp',sr,30,90);
        clear gpow;
        for i = 1:channels
            h = hilbert(gamma(i,:)); gpow(i,:) = abs(h);
        end

        for i = 1:length(result.msstamps)
            lfpresp(i,:,:) = lfp(result.msstamps(i)-prestim+1:result.msstamps(i)+result.stimduration+poststim,:);
            gammaresp(i,:,:) = gpow(:,result.msstamps(i)-prestim+1:result.msstamps(i)+result.stimduration+poststim)';
        end

        if isfield(result,'light') & find(result.light) & dolightstuff

            window = 900+prestim:1411+prestim;
            for chan = 1:size(lfpresp,3)
                for i = 1:size(lfpresp,1)
                    [lfpspect(chan,i,:),trialfax] = pmtm(lfpresp(i,window,chan),3,[],sr);
                end    
            end

            if length(unique(result.gratingInfo.Contrast))  == 1 % size block
                vari = result.gratingInfo.size; contr = 0; condstr = 'SIZE';
            else
                vari = result.gratingInfo.Contrast; contr = 1; condstr = 'CONTRAST'; % contrast block
            end
            if contr
                str0 = 'low'; str1 = 'high'; 
                contindsl0 = find(result.gratingInfo.Contrast == 0 & result.light == 0);
                contindsl1 = find(result.gratingInfo.Contrast == 0 & result.light == 1);
            else
                str0 = 'small'; str1 = 'large'; 
                contindsl0 = find(result.gratingInfo.size == 0 & result.light == 0);
                contindsl1 = find(result.gratingInfo.size == 0 & result.light == 1);
            end

            oris = unique(result.gratingInfo.Orientation); oris(find(oris == -1)) = [];
            clevels = unique(vari); clevels(find(clevels == 0)) = []; %delete control condition
            for l = 1:length(unique(result.light))
                for ori = 1:length(oris)
                    for co = 1:length(clevels)
                        thisinds = find(result.gratingInfo.Orientation == oris(ori) &...
                            vari == clevels(co) & result.light == l-1);
                        condlfpspect(:,l,ori,co,:) = nanmean(lfpspect(:,thisinds,:),2);
                        for chan = 1:size(lfpresp,3)
                            cgresp(chan,l,ori,co,:) = squeeze(nanmean(gammaresp(thisinds,:,chan)))';
                            cgpow(chan,l,ori,co) = nanmean(cgresp(chan,l,ori,co,window),5);
                        end
                    end
                end
            end


            l0inds = find(result.light == 0 & result.gratingInfo.Contrast == 1); 
            l1inds = find(result.light & result.gratingInfo.Contrast == 1); 
            
            controllfpspect(blck,1,:,:) = nanmean(lfpspect(:,contindsl0,:),2);
            controllfpspecterr(blck,1,:,:) = nanstd(lfpspect(:,contindsl0,:),1,2)./sqrt(length(contindsl0));
            controllfpspectstd(blck,1,:,:) = nanstd(lfpspect(:,contindsl0,:),1,2);
            controllfpspect(blck,2,:,:) = nanmean(lfpspect(:,contindsl1,:),2);
            controllfpspecterr(blck,2,:,:) = nanstd(lfpspect(:,contindsl1,:),1,2)./sqrt(length(contindsl1));
            controllfpspectstd(blck,2,:,:) = nanstd(lfpspect(:,contindsl1,:),1,2);
            
            %gamma power
            gwin = 600+prestim:800+prestim;
            figure
            plot(squeeze(mean(mean(gammaresp(l1inds,gwin,:),1),2))-squeeze(mean(mean(gammaresp(l0inds,gwin,:),1),2)),depthax(blck,:),'ko','markerfacecolor','k')
            axis ij
            ax = axis;
            line([ax(1),ax(2)],[375,375],'color','k','linestyle',':');
            line([ax(1),ax(2)],[500,500],'color','k','linestyle',':');
            line([ax(1),ax(2)],[800,800],'color','k','linestyle',':');
            line([0,0],[0,1000],'color','k','linewidth',2);

            figure
            subplot(2,2,1)
            imagesc(trialfax,depthax(blck,:),log(squeeze(nanmean(condlfpspect(:,1,:,1,:),3))))
            axis([-.5,100.5,depthax(blck,end),depthax(blck,1)])
            caxis([0,7])
        %     axis xy;
            title([' spectrum in depth light OFF ' str0])

            subplot(2,2,2)
            imagesc(trialfax,depthax(blck,:),log(squeeze(nanmean(condlfpspect(:,2,:,1,:),3))))
            axis([-.5,100.5,depthax(blck,end),depthax(blck,1)])
            caxis([0,7])
        %     axis xy
            title([' spectrum in depth light ON ' str0])

            subplot(2,2,3)
            imagesc(trialfax,depthax(blck,:),log(squeeze(nanmean(condlfpspect(:,1,:,5,:),3))))
            axis([-.5,100.5,depthax(blck,end),depthax(blck,1)])
            caxis([0,7])
        %     axis xy;
            title([' spectrum in depth light OFF ' str1])

            subplot(2,2,4)
            imagesc(trialfax,depthax(blck,:),log(squeeze(nanmean(condlfpspect(:,2,:,5,:),3))))
            caxis([0,7])
            axis([-.5,100.5,depthax(blck,end),depthax(blck,1)])
        %     axis xy
            title([' spectrum in depth light ON ' str1])
            
            if printyn
                figSize = [30 21];
                set(gcf,'paperunits','centimeters','papersize',figSize,'paperposition',[0 0 figSize])
                print(gcf,[printpath1 ,  printname '_spectra.pdf'], '-dpdf' );
            end

            figure
            subplot(2,2,1)
            imagesc(trialfax,depthax(blck,:),log(squeeze(nanmean(condlfpspect(:,2,:,1,:),3)))-log(squeeze(nanmean(condlfpspect(:,1,:,1,:),3))));
            caxis([-2,2])
            axis([-.5,100.5,depthax(blck,end),depthax(blck,1)])
        %     axis xy
            title(['light ON-OFF ' str0]);   

            subplot(2,2,2)
            imagesc(trialfax,depthax(blck,:),log(squeeze(nanmean(condlfpspect(:,2,:,5,:),3)))-log(squeeze(nanmean(condlfpspect(:,1,:,5,:),3))));
            caxis([-2,2])
            axis([-.5,100.5,depthax(blck,end),depthax(blck,1)])
        %     axis xy
            title(['light ON-OFF ' str1])   

            subplot(2,2,3)
            imagesc(trialfax,depthax(blck,:),log(squeeze(nanmean(condlfpspect(:,1,:,5,:),3)))-log(squeeze(nanmean(condlfpspect(:,1,:,1,:),3))));
            caxis([-2,2])
            axis([-.5,100.5,depthax(blck,end),depthax(blck,1)])
        %     axis xy
            title(['light OFF ' str1 ' - ' str0]);  

            subplot(2,2,4)
            imagesc(trialfax,depthax(blck,:),log(squeeze(nanmean(condlfpspect(:,2,:,5,:),3)))-log(squeeze(nanmean(condlfpspect(:,2,:,1,:),3))));
            caxis([-2,2])
            axis([-.5,100.5,depthax(blck,end),depthax(blck,1)])
        %     axis xy
            title(['light ON ' str1 ' - ' str0]);  
            
            if printyn
                figSize = [30 21];
                set(gcf,'paperunits','centimeters','papersize',figSize,'paperposition',[0 0 figSize])
                print(gcf,[printpath2 ,  printname '_spectradiff.pdf'], '-dpdf' );
            end
            
            figure
            subplot(2,2,1)
            imagesc(timeax,depthax(blck,:),squeeze(nanmean(cgresp(:,1,:,1,:),3)));
            caxis([0,100])
            title(['gamma power in time L0 ' str0]);

            subplot(2,2,2)
            imagesc(timeax,depthax(blck,:),squeeze(nanmean(cgresp(:,1,:,5,:),3)));
            caxis([0,100])
            title(['gamma power in time L0 ' str1])   

            subplot(2,2,3)
            imagesc(timeax,depthax(blck,:),squeeze(nanmean(cgresp(:,2,:,1,:),3)));
            caxis([0,100])
            title(['gamma power in time L1 ' str0]);  

            subplot(2,2,4)
            imagesc(timeax,depthax(blck,:),squeeze(nanmean(cgresp(:,2,:,5,:),3)));
            caxis([0,100])
            title(['gamma power in time L1 ' str1]);  
            
            
            figure
            subplot(2,2,1)
            imagesc(clevels,depthax(blck,:),squeeze(nanmean(cgpow(:,1,:,:),3)));
            caxis([0,75])
            colorbar
            title(['L0 gamma power with ' condstr ': ' str0 ' to ' str1 ]);

            subplot(2,2,2)
            imagesc(clevels,depthax(blck,:),squeeze(nanmean(cgpow(:,2,:,:),3)));
            caxis([0,75])
            colorbar
            title(['L1 gamma power with ' condstr ': ' str0 ' to ' str1 ])   

            subplot(2,2,3)
            imagesc(clevels,depthax(blck,:),squeeze(nanmean(cgpow(:,2,:,:),3))-squeeze(nanmean(cgpow(:,1,:,:),3)));
            caxis([-10,10])
            colorbar
            title(['gamma power difference with light L1-L0' ]);  

            if printyn
                figSize = [30 21];
                set(gcf,'paperunits','centimeters','papersize',figSize,'paperposition',[0 0 figSize])
                print(gcf,[printpath3 ,  printname '_gammadiff.pdf'], '-dpdf' );
            end

            condspect(blck,:,:,:,:,:) = condlfpspect;
            condgammapow(blck,:,:,:,:) = cgpow;
            condgammaresp(blck,:,:,:,:,:) = cgresp;
            
        end
    end
    save(popfile, '-v7.3', 'condspect','condgammapow','condgammaresp','depthax','trialfax','clevels','timeax')
else
    load(popfile);
end

disp('')'