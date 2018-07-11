function sizecont_populationana


% % % SOM later population
% animalids = {'160729', '160801', '160802', '160804_2'};
% blocks    = [ 4,        6,        6,        7];
% animal    = [ 1,        2,        3,        4];
% electrodes =[[1,16];   [1,16];   [1,16];   [1,16]  ];
% penangle =  [ 25,       25,       25,       25];
% printpath = 'C:\Users\Julia\work\data\populations\SOM_Halo_later\sizecont\units\';
% runprintpath = 'C:\Users\Julia\work\data\populations\SOM_Halo_later\sizecont\running\';
% oscillprintpath = 'C:\Users\Julia\work\data\populations\SOM_Halo_later\sizecont\withlfpoverview\';
% coszprintpath = 'C:\Users\Julia\work\data\populations\SOM_Halo_later\sizecont\szcont\';
% popfile = 'C:\Users\Julia\work\data\populations\SOM_Halo_later\sizecont\sizecont_population.mat';

% % % VIP Halo population
% animalids =  {'171201', '171205', '180125', '180404', '180410', '180412', '180413'};
% blocks =     [ 2,        7,        4,        9,        6,        6,        10];
% animal    =  [ 1,        1,        2,        3,        4,        5,        5];
% electrodes = [[1,16];   [1,16];   [1,16];   [1,16];   [1,16];   [1,16];   [1,16]];
% penangle =   [25,        25,       25,       25,       25,       25,       25];
% printpath = 'C:\Users\Julia\work\data\populations\VIP_Halo\sizecont\units\';
% runprintpath = 'C:\Users\Julia\work\data\populations\VIP_Halo\sizecont\running\';
% oscillprintpath = 'C:\Users\Julia\work\data\populations\VIP_Halo\sizecont\withlfpoverview\';
% coszprintpath = 'C:\Users\Julia\work\data\populations\VIP_Halo\sizecont\szcont\';
% popfile = 'C:\Users\Julia\work\data\populations\VIP_Halo\sizecont\sizecont_population.mat';

% % SOM VIP combined population
animalids = {'160729', '160801', '160802', '160804_2', '171201', '171205', '180125', '180404', '180410', '180412', '180413'};
blocks    = [ 4,        6,        6,        7,          2,        7,        4,        9,        6,        6,        10];
animal    = [ 1,        2,        3,        4,          5,        5,        6,        7,        8,        9         9];
electrodes =[[1,16];   [1,16];   [1,16];   [1,16];     [1,16];   [1,16];   [1,16];   [1,16];   [1,16];   [1,16];   [1,16]  ];
penangle =  [ 25,       25,       25,       25,         25,       25,       25,       25,       25,       25,       25];
printpath = 'C:\Users\Julia\work\data\populations\SOMVIPcombined\sizecont\units\';
runprintpath = 'C:\Users\Julia\work\data\populations\SOMVIPcombined\sizecont\running\';
oscillprintpath = 'C:\Users\Julia\work\data\populations\SOMVIPcombined\sizecont\withlfpoverview\';
coszprintpath = 'C:\Users\Julia\work\data\populations\SOMVIPcombined\sizecont\szcont\';
popfile = 'C:\Users\Julia\work\data\populations\SOMVIPcombined\sizecont\sizecont_population.mat';


tic
lcol = 'r'; %lasercolor

recalculate = 0;
printyn = 1;

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
                allinds = sort(getSpecificIndices(result, 'contconds'));
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
            
            
            for i = 1:length(msstamps)
                resp(i,:) = chan(msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim);  
                lfpresp(i,:) = lfp(msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim);               
                msl(i) = mean(find(resp(i,respwin)));
            end            
            msta = linspace(-prestim,trialdur+poststim,size(resp,2));

            frs = sum(resp(:,respwin),2)./(length(respwin)/1000);
            bl = sum(resp(:,1:prestim),2)./(prestim/1000);
            sc = sum(resp(:,respwin),2); % spike count
                        
            l0s = find(light == 0);            
            firingrates(cll,:) = frs(l0s);
            speeds(cll,:) = meanspeed(l0s);
            
            
            binwidth = 30;              
            sizes = unique(gratingInfo.size);  sizes(sizes == 0) = []; sizes(sizes==10) =[];%!!!!!!!!!! for the one
            oris = unique(gratingInfo.Orientation); oris(oris == -1) = [];
            contrasts = unique(gratingInfo.Contrast); contrasts(find(contrasts == 0)) = [];
            shownsizes(cll,:) = sizes;   shownclevels(cll,:) = contrasts; shownoris(cll,:) = oris;
            for l = 1:2
                for ori = 1:length(oris)
                    for sz = 1:length(sizes)
                        for co = 1:length(contrasts)
                            thisinds = find(gratingInfo.Orientation == oris(ori) &...
                                gratingInfo.size == sizes(sz) & gratingInfo.Contrast == contrasts(co) & ...
                                light == l-1);                            
                            
                            condresp(cll,l,ori,sz,co,:) = mean(resp(thisinds,:),1);
                            allresp{cll,l,ori,sz,co,:,:} = resp(thisinds,:);
                            condlfpresp(cll,l,ori,sz,co,:) = mean(lfpresp(thisinds,:),1);
                            condfr(cll,l,ori,sz,co) = mean(frs(thisinds));%-mean(bl); 
                            conderr(cll,l,ori,sz,co) =std(frs(thisinds))./sqrt(length(thisinds));
                            allfr{l,ori,sz,co} = frs(thisinds);
                            [bincondresp(cll,l,ori,sz,co,:),bta] = binit(squeeze(condresp(cll,l,ori,sz,co,:)),binwidth); 
                            condmsl(cll,l,ori,sz,co) = mean(msl(thisinds)); % mean spike latency
                            %TODO fix to be Hz
                            condfiltresp(cll,l,ori,sz,co,:) = filter(excit_kernel,1,mean(resp(thisinds,:),1));

                            % test each cell for each condition ranksum
                            if l == 1
                                l1inds = find(gratingInfo.Orientation == oris(ori) &...
                                gratingInfo.size == sizes(sz) & ...
                                light == 1);
                                if ~isempty(l1inds) & ~isempty(thisinds)
                                    condlightmodp(cll,ori,sz,co) = ranksum(frs(thisinds),frs(l1inds));
                                else
                                    condlightmodp(cll,ori,sz,co) = NaN;
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
                            msreliab(cll,l,ori,sz,co) = nanmean(mscc);
                            binreliab(cll,l,ori,sz,co) = nanmean(bincc);
                            eckerreliability(cll,l,ori,sz,co) = var(frs(thisinds))/var(frs);

                            % running
                            thisruninds = intersect(thisinds,oktrials);
                            if ~isempty(thisruninds)
                                runcondresp(cll,l,ori,sz,co,:) = mean(resp(thisruninds,:),1);
                                runcondfr(cll,l,ori,sz,co) = mean(frs(thisruninds));
                                runconderr(cll,l,ori,sz,co) = std(frs(thisruninds))./sqrt(length(thisruninds));
                            else
                                runcondresp(cll,l,ori,sz,co,:) = nan(1,size(resp,2));
                                runcondfr(cll,l,ori,sz,co) = NaN;
                                runconderr(cll,l,ori,sz,co) = NaN;
                            end

                            thisstillinds = intersect(thisinds,stilltrials);
                            if ~isempty(thisstillinds)
                                stillcondresp(cll,l,ori,sz,co,:) = nanmean(resp(thisstillinds,:),1);
                                stillcondfr(cll,l,ori,sz,co) = nanmean(frs(thisstillinds));
                                stillconderr(cll,l,ori,sz,co) = nanstd(frs(thisstillinds))./sqrt(length(thisstillinds));
                            else
                                stillcondresp(cll,l,ori,sz,co,:) = nan(1,size(resp,2));
                                stillcondfr(cll,l,ori,sz,co) = NaN;
                                stillconderr(cll,l,ori,sz,co) = NaN;
                            end
                            
                        end
                    end
                end
            end
            
            bincondresp(cll,:,:,:,:,:) = bincondresp(cll,:,:,:,:,:).*(1000/binwidth);
            maxdriven = max(max(max(squeeze(condfr(cll,1,:,:,:)))))-mean(bl);
            
            nrun(cll) = length(oktrials);
            nstill(cll) = length(stilltrials);
            
            cllname{cll} = files(fi).name;
            printname = files(fi).name;
            printname(find(printname=='_')) = ' ';
            animalno(cll) = animal(blck);
            recording(cll) = blck;
            depth(cll) = result.depth;
            pangle(cll) = penangle(blck);            
            
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
            blscondfr = squeeze(condfr(cll,:,:,:,:))-bfr(cll);

            %how much does the firing rate change for each condition
            condchange(cll,:,:,:) = squeeze(blscondfr(2,:,:,:))-squeeze(blscondfr(1,:,:,:));            
            
            %control firing rates (also with running)
            contindsnl = find(gratingInfo.Contrast == 0 & light == 0);
            controlfr(cll,1) = mean(frs(contindsnl));
            controlerr(cll,1) = std(frs(contindsnl))./sqrt(length(contindsnl)); 
            
            contindsl = find(gratingInfo.Contrast == 0 & light == 1);
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

            prefsize(cll) = find(nanmean(nanmean(condfr(cll,1,:,:,:),3),5) == max(nanmean(nanmean(condfr(cll,1,:,:,:),3),5)),1);
            prefcont(cll) = find(nanmean(nanmean(condfr(cll,1,:,:,:),3),4) == max(nanmean(nanmean(condfr(cll,1,:,:,:),3),4)),1);


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
            errorbar(oris,squeeze(condfr(cll,2,:,prefsize(cll),prefcont(cll))),squeeze(conderr(cll,2,:,prefsize(cll),prefcont(cll))),'o-','color',lcol,'markersize',8,'linewidth',2)
            hold on
            errorbar(oris,squeeze(condfr(cll,1,:,prefsize(cll),prefcont(cll))),squeeze(conderr(cll,1,:,prefsize(cll),prefcont(cll))),'ko-','markersize',8,'linewidth',2)
            xlabel('shown orientation')
            ylabel('Firing rate [Hz]')
            set(gca,'xtick',oris)
            legend({'Light ON','Light OFF'})

            subplot(2,2,3)
            errorbar([0,sizes],[controlfr(cll,1),squeeze(nanmean(condfr(cll,1,:,:,prefcont(cll)),3))'],[controlerr(cll,1), squeeze(nanmean(conderr(cll,1,:,:,prefcont(cll)),3))'],'ko-','color','k','markersize',8,'linewidth',2);
            hold on
            errorbar([0,sizes],[controlfr(cll,2),squeeze(nanmean(condfr(cll,2,:,:,prefcont(cll)),3))'],[controlerr(cll,2), squeeze(nanmean(conderr(cll,2,:,:,prefcont(cll)),3))'],'o-','color',lcol,'markersize',8,'linewidth',2);
            xlabel('shown patch size [vd]')
            ylabel('Firing rate [Hz]')
            legend({'Light ON','Light OFF'})    
            set(gca,'xtick',sizes)  
            title(['cell ' int2str(cll) ' preferred orientations' ' depth: ' int2str(result.depth) '  ' printname])

            subplot(2,2,4)
            errorbar([0,contrasts],[controlfr(cll,1),squeeze(nanmean(condfr(cll,1,:,prefsize(cll),:),3))'],[controlerr(cll,1), squeeze(nanmean(conderr(cll,1,:,prefsize(cll),:),3))'],'ko-','color','k','markersize',8,'linewidth',2);
            hold on
            errorbar([0,contrasts],[controlfr(cll,2),squeeze(nanmean(condfr(cll,2,:,prefsize(cll),:),3))'],[controlerr(cll,2), squeeze(nanmean(conderr(cll,2,:,prefsize(cll),:),3))'],'o-','color',lcol,'markersize',8,'linewidth',2);
            xlabel('contrast')
            ylabel('Firing rate [Hz]')
            legend({'Light ON','Light OFF'})    
            set(gca,'xtick',contrasts)  
            title(['cell ' int2str(cll) ' preferred orientations' ' depth: ' int2str(result.depth) '  ' printname])

            if printyn
                figSize = [30 21];
                set(gcf,'paperunits','centimeters','papersize',figSize,'paperposition',[0 0 figSize])
                if cll<10, printi = ['0', int2str(cll)]; else printi = int2str(cll); end
                print([printpath ,  printi '__' files(fi).name '.pdf'],'-dpdf')
            end
            
            disp([files(fi).name '   done'])
            cll = cll + 1;            
        end
    end
    save(popfile, '-v7.3');
else
    load(popfile);   
end


secpersamp = 1/30000;
interpf = secpersamp/10;
swidthms = swidth*interpf*1000;
prsv = swidthms>=.38; pfsv = swidthms<=.36;
prs = find(prsv); pfs = find(pfsv);
depth = depth.*cosd(22).*cosd(pangle);

% halo expressing clls
phe = zeros(1,length(depth));

% % putative halo expressing for SOM later pop: 
phe([1]) = 1;

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

for i = 1:length(depth)
    for l = 1:2
        for sz = 1:3
            for cl = 1:5
                [mx,mxi] = max(runcondfr(i,l,:,sz,cl),[],3);
                pofrr1(i,l,sz,cl) = runcondfr(i,l,mxi,sz,cl);
                poerrr1(i,l,sz,cl) = runconderr(i,l,mxi,sz,cl);
            end
        end
    end
end
mofrr1 = squeeze(nanmean(runcondfr,3));
moerrr1 = squeeze(nanmean(runconderr,3));

for i = 1:length(depth)
    for j = 1:5
        nmcoszfr(i,:,:,j) = mofrr1(i,:,:,j)./max(mofrr1(i,1,:,j),[],3);
        nmcoctrl(i,1,1,j) = r1l0ctrlfr(i)./max(mofrr1(i,1,:,j),[],3);
        nmcoctrl(i,2,1,j) = r1l1ctrlfr(i)./max(mofrr1(i,1,:,j),[],3);
        si(i,j) = (max(mofrr1(i,1,:,j))-mofrr1(i,1,3,j))/max(mofrr1(i,1,:,j))-r1l0ctrlfr(i);
        sinbl(i,j) = (max(mofrr1(i,1,:,j))-mofrr1(i,1,3,j))/max(mofrr1(i,1,:,j));
    end
end
normcoszfr = cat(3,nmcoctrl,nmcoszfr);
normcoszfr(isinf(normcoszfr)) = NaN;

% TODO change to allresp
ex = 30;
ta = [-299:2700]; prestim = 300;
smalll0 = []; smalll1 = []; largel0 = []; largel1 = [];
for i = 1:size(allresp,3)
    smalll0 = [smalll0; allresp{ex,1,i,1,1}];
    smalll1 = [smalll1; allresp{ex,2,i,1,1}];
    largel0 = [largel0; allresp{ex,1,i,3,1}];
    largel1 = [largel1; allresp{ex,2,i,3,1}];
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
line([0,2000],[22,22],'color','k','linewidth',3);
% line([0,2000],[22,22],'color','k','linewidth',3);
% line([500,1500],[21,21],'color','r','linewidth',3);
axis([-300,2500,0,23])
subplot(2,1,2)
hold on;
for i = 1:size(cond2,1)
    if ~isempty(find(cond2(i,:)))
        plot(find(cond2(i,:))-prestim,i,'ko','MarkerSize',1.5,'MarkerFaceColor','k')
    end
end
axis([-300,2500,0,23])


data = normcoszfr;
cllt = l23rs;
figure
errorbar([0,sizes],squeeze(nanmean(data(cllt,1,:,1))),squeeze(nanstd(data(cllt,1,:,1)))./sqrt(length(find(cllt))),'o-','color',[.9,.9,.9],'markerfacecolor',[.9,.9,.9])
hold on
errorbar([0,sizes],squeeze(nanmean(data(cllt,1,:,2))),squeeze(nanstd(data(cllt,1,:,2)))./sqrt(length(find(cllt))),'o-','color',[.7,.7,.7],'markerfacecolor',[.7,.7,.7])
errorbar([0,sizes],squeeze(nanmean(data(cllt,1,:,3))),squeeze(nanstd(data(cllt,1,:,3)))./sqrt(length(find(cllt))),'o-','color',[.5,.5,.5],'markerfacecolor',[.5,.5,.5])
errorbar([0,sizes],squeeze(nanmean(data(cllt,1,:,4))),squeeze(nanstd(data(cllt,1,:,4)))./sqrt(length(find(cllt))),'o-','color',[.3,.3,.3],'markerfacecolor',[.3,.3,.3])
errorbar([0,sizes],squeeze(nanmean(data(cllt,1,:,5))),squeeze(nanstd(data(cllt,1,:,5)))./sqrt(length(find(cllt))),'o-','color',[.1,.1,.1],'markerfacecolor',[.1,.1,.1])


frdiff = squeeze(mofrr1(:,2,:,:)-mofrr1(:,1,:,:));
pci = squeeze(mofrr1(:,2,:,:)./mofrr1(:,1,:,:));
pci(isinf(pci)) = NaN;
lpci = log2(pci);
lpci(isinf(lpci)) = NaN;
omi = (squeeze(mofrr1(:,2,:,:))-squeeze(mofrr1(:,1,:,:)))./(squeeze(mofrr1(:,2,:,:))+squeeze(mofrr1(:,1,:,:)));

figure
plot(1,lpci(l23rs,3,1),'k.','markersize',15)
hold on
plot(2,lpci(l23rs,3,5),'k.','markersize',15)
axis([0,3,-5,4])
a = find(l23rs);
for i = 1:length(a)
    line([1,2],[lpci(a(i),3,1),lpci(a(i),3,5)],'color','k')
end
line([0,3],[0,0],'color','k')
set(gca,'xtick',[1,2])
set(gca,'xticklabel',{'low','high'})
ylabel('OMR')

cllt = l23rs;
%low contrast large vs small
figure
plot(squeeze(mofrr1(cllt,1,1,1)),squeeze(mofrr1(cllt,1,3,1)),'k.')
axis([0,7,0,7])
axis square
refline(1,0);
xlabel('firing rate small - low contrast')
ylabel('firing rate large - low contrast')

%high contrast large vs small
figure
plot(squeeze(mofrr1(cllt,1,1,5)),squeeze(mofrr1(cllt,1,3,5)),'k.')
axis([0,14,0,14])
axis square
refline(1,0);
xlabel('firing rate small - high contrast')
ylabel('firing rate large - high contrast')

%low contrast large light vs no light
figure
plot(squeeze(mofrr1(cllt,1,3,1)),squeeze(mofrr1(cllt,2,3,1)),'k.')
axis([0,14,0,14])
axis square
refline(1,0);
xlabel('firing rate large - low contrast - control')
ylabel('firing rate large - low contrast - light')

%high contrast large light vs no light
figure
plot(squeeze(mofrr1(cllt,1,3,5)),squeeze(mofrr1(cllt,2,3,5)),'k.')
axis([0,26,0,26])
axis square
refline(1,0);
xlabel('firing rate large - high contrast - control')
ylabel('firing rate large - high contrast - light')

a = find(cllt);
for i = 1:length(a)
    figure
    errorbar([0,contrasts],[r1l0ctrlfr(a(i)),squeeze(mofrr1(a(i),1,1,:))'],[r1l0ctrlerr(a(i)),squeeze(moerrr1(a(i),1,1,:))'],'o-','color',[.7,.7,.7],'markerfacecolor',[.7,.7,.7]);
    hold on
    errorbar([0,contrasts],[r1l0ctrlfr(a(i)),squeeze(mofrr1(a(i),1,2,:))'],[r1l0ctrlerr(a(i)),squeeze(moerrr1(a(i),1,2,:))'],'o-','color',[.4,.4,.4],'markerfacecolor',[.4,.4,.4]);
    errorbar([0,contrasts],[r1l0ctrlfr(a(i)),squeeze(mofrr1(a(i),1,3,:))'],[r1l0ctrlerr(a(i)),squeeze(moerrr1(a(i),1,3,:))'],'o-','color','k','markerfacecolor','k');
    legend('8 degrees','20 degrees','60 degrees')
    xlabel('contrast level')
    ylabel('firing rate')
end

for i = 1:length(a)
    figure
    errorbar([0,sizes],[r1l0ctrlfr(a(i)),squeeze(mofrr1(a(i),1,:,1))'],[r1l0ctrlerr(a(i)),squeeze(moerrr1(a(i),1,:,1))'],'o-','color',[.9,.9,.9],'markerfacecolor',[.9,.9,.9]);
    hold on
    errorbar([0,sizes],[r1l0ctrlfr(a(i)),squeeze(mofrr1(a(i),1,:,2))'],[r1l0ctrlerr(a(i)),squeeze(moerrr1(a(i),1,:,2))'],'o-','color',[.7,.7,.7],'markerfacecolor',[.7,.7,.7]);
    errorbar([0,sizes],[r1l0ctrlfr(a(i)),squeeze(mofrr1(a(i),1,:,3))'],[r1l0ctrlerr(a(i)),squeeze(moerrr1(a(i),1,:,3))'],'o-','color',[.5,.5,.5],'markerfacecolor',[.5,.5,.5]);
    errorbar([0,sizes],[r1l0ctrlfr(a(i)),squeeze(mofrr1(a(i),1,:,4))'],[r1l0ctrlerr(a(i)),squeeze(moerrr1(a(i),1,:,4))'],'o-','color',[.3,.3,.3],'markerfacecolor',[.3,.3,.3]);
    errorbar([0,sizes],[r1l0ctrlfr(a(i)),squeeze(mofrr1(a(i),1,:,5))'],[r1l0ctrlerr(a(i)),squeeze(moerrr1(a(i),1,:,5))'],'o-','color','k','markerfacecolor','k');
    xlabel('stimulus size')
    ylabel('firing rate')
    legend('low','-','-','-','high')
end

for i = 1:length(a)
    clf;
    for j  = 1:5
        subplot(2,5,j)
        errorbar([0,sizes],[r1l0ctrlfr(a(i)),squeeze(mofrr1(a(i),1,:,j))'],[r1l0ctrlerr(a(i)),squeeze(moerrr1(a(i),1,:,j))'],'ko-','markerfacecolor','k','linewidth',2);
        hold on
        errorbar([0,sizes],[r1l1ctrlfr(a(i)),squeeze(mofrr1(a(i),2,:,j))'],[r1l1ctrlerr(a(i)),squeeze(moerrr1(a(i),2,:,j))'],'ro-','markerfacecolor','r','linewidth',2);
        ax = axis;
        axis([-10,62,ax(3),ax(4)])
        xlabel('size')
        if j == 1, title('low'); ylabel('firing rate'); elseif j == 5 title('high'); end
    end
    for j = 1:3
        subplot(2,3,j+3)
        errorbar([0,contrasts],[r1l0ctrlfr(a(i)),squeeze(mofrr1(a(i),1,j,:))'],[r1l0ctrlerr(a(i)),squeeze(moerrr1(a(i),1,j,:))'],'ko-','markerfacecolor','k','linewidth',2);
        hold on
        errorbar([0,contrasts],[r1l1ctrlfr(a(i)),squeeze(mofrr1(a(i),2,j,:))'],[r1l1ctrlerr(a(i)),squeeze(moerrr1(a(i),2,j,:))'],'ro-','markerfacecolor','r','linewidth',2);
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
        nmcoszfr(i,:,:,j) = mofrr1(i,:,:,j)./max(mofrr1(i,1,:,j),[],3);
        nmcoctrl(i,1,1,j) = r1l0ctrlfr(i)./max(mofrr1(i,1,:,j),[],3);
        nmcoctrl(i,2,1,j) = r1l1ctrlfr(i)./max(mofrr1(i,1,:,j),[],3);
        si(i,j) = (max(mofrr1(i,1,:,j))-mofrr1(i,1,3,j))/max(mofrr1(i,1,:,j))-r1l0ctrlfr(i);
        sinbl(i,j) = (max(mofrr1(i,1,:,j))-mofrr1(i,1,3,j))/max(mofrr1(i,1,:,j));
    end
end
normcoszfr = cat(3,nmcoctrl,nmcoszfr);
normcoszfr(isinf(normcoszfr)) = NaN;

data = normcoszfr;
cllt = l23rs;
figure
errorbar([0,sizes],squeeze(nanmean(data(cllt,1,:,1))),squeeze(nanstd(data(cllt,1,:,1)))./sqrt(length(find(cllt))),'o-','color',[.9,.9,.9],'markerfacecolor',[.9,.9,.9])
hold on
errorbar([0,sizes],squeeze(nanmean(data(cllt,1,:,2))),squeeze(nanstd(data(cllt,1,:,2)))./sqrt(length(find(cllt))),'o-','color',[.7,.7,.7],'markerfacecolor',[.7,.7,.7])
errorbar([0,sizes],squeeze(nanmean(data(cllt,1,:,3))),squeeze(nanstd(data(cllt,1,:,3)))./sqrt(length(find(cllt))),'o-','color',[.5,.5,.5],'markerfacecolor',[.5,.5,.5])
errorbar([0,sizes],squeeze(nanmean(data(cllt,1,:,4))),squeeze(nanstd(data(cllt,1,:,4)))./sqrt(length(find(cllt))),'o-','color',[.3,.3,.3],'markerfacecolor',[.3,.3,.3])
errorbar([0,sizes],squeeze(nanmean(data(cllt,1,:,5))),squeeze(nanstd(data(cllt,1,:,5)))./sqrt(length(find(cllt))),'o-','color',[.1,.1,.1],'markerfacecolor',[.1,.1,.1])

cl = 1;
figure
errorbar([0,sizes],squeeze(nanmean(data(cllt,1,:,cl))),squeeze(nanstd(data(cllt,1,:,cl)))./sqrt(length(find(cllt))),'ko-','markerfacecolor','k','linewidth',2)
hold on
errorbar([0,sizes],squeeze(nanmean(data(cllt,2,:,cl))),squeeze(nanstd(data(cllt,2,:,cl)))./sqrt(length(find(cllt))),'ro-','markerfacecolor','r','linewidth',2)
xlabel('size')
ylabel('normalized FR')
title(['contrast: ' num2str(contrasts(cl))]);


for i = 1:length(depth)
    for j = 1:3
        nmtofccoszfr(i,:,j,:) = mofrr1(i,:,j,:)./mofrr1(i,1,j,5); % normalized to full contrast
        nmtomxcoszfr(i,:,j,:) = mofrr1(i,:,j,:)./max(mofrr1(i,1,j,:),[],4); % normalized to max response at that size
        nmxctrl(i,1,j) = r1l0ctrlfr(i)./max(mofrr1(i,1,j,:)); % normalize control to max
        nmxctrl(i,2,j) = r1l1ctrlfr(i)./max(mofrr1(i,1,j,:)); % normalize control to max
        nfcctrl(i,1,j) = r1l0ctrlfr(i)./mofrr1(i,1,j,5); % normalize control to max clevel
        nfcctrl(i,2,j) = r1l1ctrlfr(i)./mofrr1(i,1,j,5); % normalize control to max clevel
    end
end
normmxcoszfr = cat(4,nmxctrl,nmtomxcoszfr);
normfccoszfr = cat(4,nfcctrl,nmtofccoszfr);
normmxcoszfr(isinf(normmxcoszfr)) = NaN;
normfccoszfr(isinf(normfccoszfr)) = NaN;


data = normmxcoszfr;
sz = 1;
figure
errorbar([0,contrasts],squeeze(nanmean(data(cllt,1,sz,:))),squeeze(nanstd(data(cllt,1,sz,:)))./sqrt(length(find(cllt))),'ko-','markerfacecolor','k','linewidth',2)
hold on
errorbar([0,contrasts],squeeze(nanmean(data(cllt,2,sz,:))),squeeze(nanstd(data(cllt,2,sz,:)))./sqrt(length(find(cllt))),'ro-','markerfacecolor','r','linewidth',2)
xlabel('contrast')
ylabel('normalized FR')
title(['size: ' num2str(sizes(sz))]);

figure
errorbar([0,contrasts],squeeze(nanmean(data(cllt,1,1,:))),squeeze(nanstd(data(cllt,1,1,:)))./sqrt(length(find(cllt))),'o-','color',[.7,.7,.7],'markerfacecolor',[.7,.7,.7])
hold on
errorbar([0,contrasts],squeeze(nanmean(data(cllt,1,2,:))),squeeze(nanstd(data(cllt,1,2,:)))./sqrt(length(find(cllt))),'o-','color',[.4,.4,.4],'markerfacecolor',[.4,.4,.4])
errorbar([0,contrasts],squeeze(nanmean(data(cllt,1,3,:))),squeeze(nanstd(data(cllt,1,3,:)))./sqrt(length(find(cllt))),'o-','color',[0,0,0],'markerfacecolor',[0,0,0])
errorbar([0,contrasts],squeeze(nanmean(data(cllt,2,1,:))),squeeze(nanstd(data(cllt,2,1,:)))./sqrt(length(find(cllt))),'o-','color',[1,.7,.7],'markerfacecolor',[1,.7,.7])
errorbar([0,contrasts],squeeze(nanmean(data(cllt,2,2,:))),squeeze(nanstd(data(cllt,2,2,:)))./sqrt(length(find(cllt))),'o-','color',[1,.4,.4],'markerfacecolor',[1,.4,.4])
errorbar([0,contrasts],squeeze(nanmean(data(cllt,2,3,:))),squeeze(nanstd(data(cllt,2,3,:)))./sqrt(length(find(cllt))),'o-','color',[1,0,0],'markerfacecolor',[1,0,0])
legend('8 L0','20 L0','60 L0','8 L1','20 L1','60 L1')
