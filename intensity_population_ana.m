function intensity_population_ana

% %SCNN-CHR2 Population 
% animalids = {'150323', '150325','150326'};
% blocks    = [7,         7,       5];
% animal =    [1          2,       3];
% popfile = 'C:\Users\Julia\work\data\populations\scnn_chr2\intensity_population.mat';


% SOM later population
animalids = {'150331', '150401','150527','150529','150602','150603','150825','150831','150902', '151023', '151109', '151110'};
blocks    = [5,         9,       13,      9,       9,       8,       16,      10,      15,       17,       11,       16];
animal    = [1,         2,       3,       4,       5,       6,       7,       8,       9,        11,       12,       13];
electrodes =[[1,32];    [1,32];  [1,32];  [1,32];  [1,32];  [1,16]; [17,32]; [1,16];  [1,16];   [17,32];  [17,32];  [1,16]];
penangle =  [25,        25,      25,      25,      25,      10,      25,      25,      25,       25,       25,       25];
% age       [P6,        P6,      P2?(P0), P2?(P0), P2?(P0), P1
popfile = 'C:\Users\Julia\work\data\populations\SOM_Halo_later\intensity\intensity_population.mat';

recalculate = 0;

prestim = 300;
% poststim = 300;
poststim = 300;
respwin = 501:1500; % after stimulus onset
% respwin = 501:4500; % after stimulus onset
respwin = respwin+prestim;


if ~exist(popfile) || recalculate
    
    cll = 1;
    for blck = 1: length(blocks)

        supath = ['C:\Users\Julia\work\data\' animalids{blck} '\singleunits\'];
        basename = [animalids{blck} '_block' int2str(blocks(blck)) '_tet'];

        files = dir([supath, basename, '*.mat']);

        prestim = 300;
        poststim = 700;
        respwin = 501:1500; % after stimulus onset
        respwin = respwin+prestim;
        
        for fi = 1:length(files)

            if strfind(files(fi).name, 'MU')
                continue;
            end

            load([supath, files(fi).name]);


            % calc spiking to see if includable
            msStimes = round(result.spikes);
            if ~isempty(msStimes) & msStimes(1) == 0, msStimes(1) = 1; end

            chan = zeros(1,length(result.lfp));
            chan(msStimes) = 1;

            wvchan = find(var(result.waveforms) == max(var(result.waveforms)));

            trialdur = result.stimduration*1000;
            msstamps = result.msstamps;
            if length(msstamps)~=length(result.light)
        %         msstamps(385) = []; % for 140703 block 5
        %         result.msstamps = msstamps;
        %         save([supath, files(fi).name],'result');
                pause;
            end

            for i = 1:length(msstamps)
                resp(i,:) = chan(msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim); 
                lfpresp(i,:) = result.lfp(msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim);
                speed(i,:) = result.runspeed(msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim);     
            end


            % figure out sufficiently high and nonvariable runspeed trials
            meanspeed = mean(speed(:,respwin),2);
            stdspeed = std(speed(:,respwin),1,2);
            notstill = find(meanspeed>1);
            okspeed = find(meanspeed>( mean(meanspeed(notstill))-(1.5*std(meanspeed(notstill))) ) );
            okvar = find(stdspeed<( mean(stdspeed(notstill))+(1.5*std(stdspeed(notstill)))) & stdspeed>.5);
            oktrials = intersect(okspeed,okvar);
            nonoktrials = 1:size(resp,1); nonoktrials(oktrials) = [];
            stilltrials = 1:size(resp,1); stilltrials(notstill) = [];


            frs = sum(resp(:,respwin),2)./(length(respwin)/1000);
            bl = sum(resp(:,1:prestim),2)./(prestim/1000);

            %determine if cell is visually modulated
            blfr = sum(resp(:,1:prestim),2);

            cellname{cll} = files(fi).name;

            i = strfind(files(fi).name, 'tet');
            if strcmp(files(fi).name(i+4),'_')
                tetno = strread(files(fi).name(i+3)); % single character number
            else        
                tetno = strread(files(fi).name(i+3:i+4)); % number >10
            end
            tetnos(cll) = tetno;
            if tetno>8
                v1(cll) = logical(0); v2(cll) = logical(1); cellstr = 'V2';
            else
                v1(cll) = logical(1); v2(cll) = logical(0); cellstr = 'V1';
            end

            spike = result.waveforms(:,wvchan);
            interpspike = spline(1:32,spike,1:.1:32);
            [adiff(cll),swidth(cll),ptr(cll),eslope(cll)] = spikequant(interpspike);

            depth(cll) = result.depth;
            cellresp(cll,:,:) = resp;
            celllfpresp(cll,:,:) = lfpresp;

            msta = linspace(-prestim,trialdur+poststim,size(resp,2));

            printname = files(fi).name;
            printname(find(printname=='_')) = ' ';


            binwidth = 50;
            oris = unique(result.gratingInfo.Orientation); oris(find(oris == 180)) = [];
            lightlevels = unique(result.light);
            for l = 1:length(lightlevels)
                for ori = 1:length(oris)            
                    thisinds = find((result.gratingInfo.Orientation == oris(ori) | result.gratingInfo.Orientation ==oris(ori)+180) &...
                        result.light == lightlevels(l));
                    condn(l,ori) = length(thisinds);
                    condresp(l,ori,:) = nanmean(resp(thisinds,:),1);
                    condresperr(l,ori,:) = nanstd(resp(thisinds,:),1,1)./sqrt(length(thisinds));
                    if ~isnan(condresp(l,ori,:))
                        [bincondresp(l,ori,:),bta] = binit(condresp(l,ori,:),binwidth);
                    else
                        bincondresp(l,ori,:) = binit(condresp(l,ori,:),binwidth);
                    end
                    binconderr(l,ori,:) = binit(condresperr(l,ori,:),binwidth);
                    cfr(l,ori) = nanmean(frs(thisinds));
                    cerr(l,ori) =nanstd(frs(thisinds))./sqrt(length(thisinds));

                end
            end
            bincondresp = bincondresp.*(1000/binwidth);
            binconderr = binconderr.*(1000/binwidth);
            bta = bta-prestim;

        %     anovavec = [cfr(:,1);cfr(:,2)];
        %     g1 = ones(numel(cfr),1); g1(1:6) = 0;
        %     g2 = [lightlevels';lightlevels'];
        %     [p(cll,:),table,stats] = anovan(anovavec,{g1 g2},'display','off');


        %     figure('name' ,['cll: ' int2str(cll), '  p vis: ' num2str(p(cll,1)) '  p light: ' num2str(p(cll,2))])
        %     for i = 1:6
        %         subplot(2,6,i)
        %         boundedline(bta, squeeze(bincondresp(i,1,:)),squeeze(binconderr(i,1,:)));
        %         hold on
        %         line([500,500],[0,max(max(max(bincondresp)))+max(max(max(binconderr)))],'color','r')
        % %         line([5500,5500],[0,max(max(max(bincondresp)))+max(max(max(binconderr)))],'color','r')
        %         line([1500,1500],[0,max(max(max(bincondresp)))+max(max(max(binconderr)))],'color','r')
        % %         axis([-300,8300,0,max(max(max(bincondresp)))+max(max(max(binconderr)))+1])
        %         axis([-300,2300,0,max(max(max(bincondresp)))+max(max(max(binconderr)))+1])
        %         title(['level: ' num2str(lightlevels(i))]);
        %         
        %         subplot(2,6,6+i)
        %         boundedline(bta, squeeze(bincondresp(i,2,:)),squeeze(binconderr(i,2,:))); 
        %         hold on
        %         line([500,500],[0,max(max(max(bincondresp)))+max(max(max(binconderr)))],'color','r')
        % %         line([5500,5500],[0,max(max(max(bincondresp)))+max(max(max(binconderr)))],'color','r') 
        %         line([1500,1500],[0,max(max(max(bincondresp)))+max(max(max(binconderr)))],'color','r')       
        % %         axis([-300,8300,0,max(max(max(bincondresp)))+max(max(max(binconderr)))+1])      
        %         axis([-300,2300,0,max(max(max(bincondresp)))+max(max(max(binconderr)))+1])
        %     end
        %     
        %     figure
        %     errorbar(lightlevels,cfr(:,1),cerr(:,1),'.')
        %     hold on
        %     errorbar(lightlevels,cfr(:,2),cerr(:,2),'r.')
        %     xlabel('lightintensity % max')
        %     ylabel('firing rate [Hz]')
        %     legend([{'no stimulus'},{'drifting grating'}])
        %     title(['cell: ' cellname{cll}, '  spikewidth: ' num2str(swidth(cll)) '  depth: ' num2str(depth(cll))])

            condfr(cll,:,:) = cfr; conderr(cll,:,:) = cerr;

            cll = cll+1;
            disp('');

        end
    end
    save(popfile, '-v7.3');
else
    load(popfile);
end

%spike classification
kmeansind = kmeans([eslope',ptr',swidth',adiff'],2);
if mean(swidth(find(kmeansind==1)))<mean(swidth(find(kmeansind==2)))  %1 is FS
    pfs = find(kmeansind==1); prs = find(kmeansind==2); pfsv = kmeansind==1; prsv = kmeansind==2;
else
    pfs = find(kmeansind==2); prs = find(kmeansind==1); pfsv = kmeansind==2; prsv = kmeansind==1;
end

secpersamp = 1/30000;
interpf = secpersamp/10;
swidthms = swidth*interpf*1000;

figure
plot(swidthms(prs),adiff(prs),'b.')
xlabel('spike width')
ylabel('amplitude diff')
hold on
plot(swidthms(pfs),adiff(pfs),'r.')
% axis([5.5,20.5,-.9,.7])

l23 = depth<375;
l4 = depth>=375&depth<=500;
l5 = depth>500&depth<=800;
l23rs = l23&prsv';
l23fs = l23&pfsv';
l4rs = l4&prsv';
l4fs = l4&pfsv';
l5rs = l5&prsv';
l5fs = l5&pfsv';

for i = 1:5
    condomi(:,i,:) = (condfr(:,i+1,:)-condfr(:,1,:))./(condfr(:,i+1,:)+condfr(:,1,:));
end

cond = l5rs;
figure
hold on
for i = find(cond)
    for j = 1:5
        plot(j,condomi(i,j,1),'o')
    end
end
for i = find(cond)
    plot(condomi(i,:,1))
end
axis([.5,5.5,-1.1,1.1])
line([.5,5.5],[0,0],'color','r')

disp('')
