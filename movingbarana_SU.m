function movingbarana_SU

animalid = '150127';
block = 5;

basepath = ['C:\Users\Julia\work\data\' animalid '\'];
% supath = [basepath 'multiunits\'];
supath = [basepath 'singleunits\'];
basename = [animalid '_block' int2str(block) '_tet'];

files = dir([supath, basename, '*.mat']);

prestim = 300;
poststim = 300;
stimduration = 6;
respwin = 1000:6000;

cell = 1;
for fi = 1:length(files)
    
%     if strfind(files(fi).name, 'MU')
%         continue;
%     end
    
    load([supath, files(fi).name]);
    
    % calc spiking to see if includable
    msStimes = round(result.spikes);
    if ~isempty(msStimes) & msStimes(1) == 0, msStimes(1) = 1; end
    
    chan = zeros(1,length(result.lfp));
    chan(msStimes) = 1;
    
    wvchan = find(var(result.waveforms) == max(var(result.waveforms)));
    
    trialdur = stimduration*1000;
    msstamps = result.msstamps;
%     if length(msstamps)~=length(result.light)
% %         msstamps(385) = []; % for 140703 block 5
% %         result.msstamps = msstamps;
% %         save([supath, files(fi).name],'result');
%         pause;
%     end
    
    for i = 1:length(msstamps)
        resp(i,:) = chan(msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim); 
        lfpresp(i,:) = result.lfp(msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim);
        [lfpspect(i,:), fax] = pmtm(lfpresp(i,2001:3000),3,1:100,1000);
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
    
    cellname{cell} = files(cell).name;
    
    i = strfind(files(cell).name, 'tet');
    if strcmp(files(cell).name(i+4),'_')
        tetno = strread(files(cell).name(i+3)); % single character number
    else        
        tetno = strread(files(cell).name(i+3:i+4)); % number >10
    end
    tetnos(cell) = tetno;
    if tetno>8
        v1(cell) = logical(0); v2(cell) = logical(1); cellstr = 'V2';
    else
        v1(cell) = logical(1); v2(cell) = logical(0); cellstr = 'V1';
    end
    
    spike = result.waveforms(:,wvchan);
    interpspike = spline(1:32,spike,1:.1:32);
    [adiff(cell),swidth(cell),ptr(cell),eslope(cell)] = spikequant(interpspike);
    
    depth(cell) = result.depth;
    cellresp(cell,:,:) = resp;
    celllfpresp(cell,:,:) = lfpresp;
    cellchan(cell,:) = chan;
    
   
    msta = linspace(-prestim,trialdur+poststim,size(resp,2));
    
    printname = files(cell).name;
    printname(find(printname=='_')) = ' ';
    
   
    binwidth = 50;
    oris = unique(result.gratingInfo.Orientation); %oris(find(oris == 180)) = [];
    lightlevels = unique(result.light);
    for ori = 1:length(oris)
        thisinds = find(result.gratingInfo.Orientation == oris(ori));
        condn(ori) = length(thisinds);
        condresp(ori,:) = nanmean(resp(thisinds,:),1);
        condresperr(ori,:) = nanstd(resp(thisinds,:),1,1)./sqrt(length(thisinds));
        if ~isnan(condresp(ori,:))
            [bincondresp(ori,:),bta] = binit(condresp(ori,:),binwidth);
        else
            bincondresp(ori,:) = binit(condresp(ori,:),binwidth);
        end
        binconderr(ori,:) = binit(condresperr(ori,:),binwidth);
        cfr(ori) = nanmean(frs(thisinds));
        cerr(ori) =nanstd(frs(thisinds))./sqrt(length(thisinds));
        
    end
    bincondresp = bincondresp.*(1000/binwidth);
    binconderr = binconderr.*(1000/binwidth);
    bta = bta-prestim;
    
    figure('name',['cell: ' int2str(cell), '  spikewidth: ' num2str(swidth(cell)) '  depth: ' num2str(depth(cell))])

    for i = 1:4
        subplot(2,4,i)
%         boundedline(bta, squeeze(bincondresp(1,:)),squeeze(binconderr(1,:)));
%         hold on
%         boundedline(bta, squeeze(bincondresp(i+1,:)),squeeze(binconderr(i+1,:)),'r');
        plot(bta, squeeze(bincondresp(1,:)));
        hold on
        plot(bta, squeeze(bincondresp(i+1,:)),'r');
        title(int2str(oris(i+1)));
        axis([-300,6300,0,max(max(max(bincondresp)))+max(max(max(binconderr)))+1])
        
        subplot(2,4,4+i)
%         boundedline(bta, squeeze(bincondresp(1,:)),squeeze(binconderr(1,:)));
%         hold on
%         boundedline(bta, squeeze(bincondresp(i+5,:)),squeeze(binconderr(i+5,:)),'r'); 
        plot(bta, squeeze(bincondresp(1,:)));
        hold on
        plot(bta, squeeze(bincondresp(i+5,:)),'r');
        title(int2str(oris(i+5)));
        axis([-300,6300,0,max(max(max(bincondresp)))+max(max(max(binconderr)))+1])
    end
    
    xoris = [-45,oris(2:9)];
    figure
    errorbar(xoris,cfr,cerr,'.')
    xlabel('orienation of bar')
    ylabel('firing rate [Hz]')
    title(['cell: ' int2str(cell), '  spikewidth: ' num2str(swidth(cell)) '  depth: ' num2str(depth(cell))])
    
    condfr(cell,:,:) = cfr; conderr(cell,:,:) = cerr;
    
    cell = cell+1;
    disp('');
    
end

%spike classification
kmeansind = kmeans([eslope',ptr',swidth',adiff'],2);
if mean(swidth(find(kmeansind==1)))<mean(swidth(find(kmeansind==2)))  %1 is FS
    pfs = find(kmeansind==1); prs = find(kmeansind==2); pfsv = kmeansind==1; prsv = kmeansind==2;
else
    pfs = find(kmeansind==2); prs = find(kmeansind==1); pfsv = kmeansind==2; prsv = kmeansind==1;
end

figure
plot(swidth(prs),adiff(prs),'b.')
xlabel('spike width')
ylabel('amplitude diff')
hold on
plot(swidth(pfs),adiff(pfs),'r.')
% axis([5.5,20.5,-.9,.7])

fsl4 = find(depth(pfs)>=375 & depth(pfs)<=500);
rsl4 = find(depth(prs)>=375 & depth(prs)<=500);
fsl23 = find(depth(pfs)<375);
rsl23 = find(depth(prs)<375);
fsl5 = find(depth(pfs)>500 & depth(pfs)<=800);
rsl5 = find(depth(prs)>500 & depth(prs)<=800);
l23 = depth<375;
l4 = depth>=375&depth<=500;
l5 = depth>500&depth<=800;


