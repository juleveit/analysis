function ori_population_decode

animalid = '150909';
block = 6;
lcol = 'r'; %lasercolor

onlymod = 0;
printyn = 1;
sfc = 0;

supath = ['C:\Users\Julia\work\data\' animalid '\singleunits\'];
basename = [animalid '_block' int2str(block) '_tet'];

files = dir([supath, basename, '*.mat']);

prestim = 300;
poststim = 700;
respwin = 501:1500; % after stimulus onset
respwin = respwin+prestim;
freqbinwidth = 5;

cell = 1;
for fi = 1:length(files)
    
    if strfind(files(fi).name, 'MU')
        continue;
    end
        
    load([supath, files(fi).name]);
    cellname{cell} = files(fi).name;
    
    i = strfind(files(fi).name, 'tet');
    if strcmp(files(fi).name(i+4),'_')
        tetno(cell) = strread(files(fi).name(i+3)); % single character number
    else
        tetno(cell) = strread(files(fi).name(i+3:i+4)); % number >10
    end

    % timestamps
    trialdur = result.stimduration*1000;
%     trialdur = result.sweeplength;

    msstamps = result.msstamps;    
     if length(msstamps)~=length(result.light)
% %         msstamps([169,336]) = []; % for 150523 block 11
%         result.msstamps = msstamps;
%         save([supath, files(fi).name],'result');
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

    wvchan(cell) = find(var(result.waveforms) == max(var(result.waveforms))); 
  
    msStimes = round(result.spikes);
    if isempty(msStimes), msStimes(1) = 0; end
    if msStimes(1) == 0, msStimes(1) = 1; end  
    channel(cell,:) = zeros(1,length(result.lfp));
    channel(cell,msStimes) = 1;
    
    binwidth = 33.3;
    for i = 1:length(msstamps)
        resp(cell,i,:) = channel(cell,msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim);
        [binresp(cell,i,:),bta] = binit(squeeze(resp(cell,i,:)),binwidth); 
    end
    
    depth(cell) = result.depth;
    
    spike = result.waveforms(:,wvchan(cell));
    interpspike = spline(1:32,spike,1:.1:32);
    [adiff(cell),swidth(cell)] = spikequant(interpspike);
    prs(cell) = swidth(cell)>120;
    
    cell = cell+1;
end
    
frs = sum(resp(:,:,respwin),3)./(length(respwin)/1000);
bl = sum(resp(:,:,1:prestim),3)./(prestim/1000);
sc = sum(resp(:,:,respwin),3);

for i = 1:length(msstamps)
    speed(i,:) = result.runspeed(msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim);
end
% figure out sufficiently high and nonvariable runspeed trials
meanspeed = mean(speed(:,respwin),2);
stdspeed = std(speed(:,respwin),1,2);
notstill = find(meanspeed>1);
speedcutoff = ( mean(meanspeed(notstill))-(1.5*std(meanspeed(notstill))) );
if speedcutoff<1 % too much variance
    speedcutoff = 1.5;
end
okspeed = find(meanspeed> speedcutoff );
okvar = find(stdspeed<( mean(stdspeed(notstill))+(1.5*std(stdspeed(notstill)))) & stdspeed>.5);
oktrials = intersect(okspeed,okvar);
nonoktrials = 1:size(resp,2); nonoktrials(oktrials) = [];
stilltrials = 1:size(resp,2); stilltrials(notstill) = [];

msta = linspace(-prestim,trialdur+poststim,size(resp,3));

for cell = 1:size(resp,1)
     %determine if cell is visually modulated
    vrfr(cell,:) = sum(resp(cell,:,prestim+60:2*prestim+60),3);
    [p,vismod(cell)] = signrank(bl(cell,:),vrfr(cell,:));
    
    %determine if cell is modulated by light
    [p,lightmod(cell)] = signrank(frs(cell,find(light)),frs(cell,find(light == 0)));
    
    lfr(cell) = mean(frs(cell,find(light)));
    nlfr(cell) = mean(frs(cell,find(light == 0)));
end

% for crosscorr concatenate all non-light trials
celltrials = reshape(resp(:,light == 0,:),size(resp,1),size(resp,3)*(size(resp,2)/2));
for i = 1:size(resp,1)
    celltrials5ms(i,:) = binit(celltrials(i,:),5);
    celltrials50ms(i,:) = binit(celltrials(i,:),50);
end
for i = 1:size(resp,1)
    for j = 2:size(resp,1)
        [c1,lags1] = xcorr(celltrials(i,:),celltrials(j,:),50);
        [c5,lags5] = xcorr(celltrials5ms(i,:),celltrials5ms(j,:),500);
        [c50,lags50] = xcorr(celltrials50ms(i,:),celltrials50ms(j,:),1000);
        subplot(2,2,1)
        plot(lags1,c1)
        subplot(2,2,2)
        plot(lags5,c5)
        subplot(2,2,3)
        plot(lags50,c50)
        title([cellname{i}, '  ' cellname{j}]);
        pause
        clf
    end
end

timebin = 100; corrbin = 33.3;
oricond = []; szcond = []; lcond = []; respmat = [];
timeoricond = []; timeszcond = []; timelcond = []; timerespmat = [];
sizes = unique(gratingInfo.size);  sizes(find(sizes == 0)) = []; %delete control condition
oris = unique(gratingInfo.Orientation); oris(find(oris == -1)) = [];
for l = 1:2
    for ori = 1:length(oris)
        for sz = 1:length(sizes)
            thisinds = find(gratingInfo.Orientation == oris(ori) &...
                gratingInfo.size == sizes(sz) & ...
                light == l-1);
            condresp(:,l,ori,sz,:) = nanmean(resp(:,thisinds,:),2); 
            condfr(:,l,ori,sz) = mean(frs(:,thisinds),2);%-mean(bl);
            conderr(:,l,ori,sz) =std(frs(:,thisinds),1,2)./sqrt(length(thisinds));
            
            respmat = [respmat, frs(:,thisinds)];
            oricond = [oricond, ones(1,length(thisinds)).*ori];
            szcond = [szcond, ones(1,length(thisinds)).*sz];
            lcond = [lcond, ones(1,length(thisinds)).*l];
            
            for i = 1:length(thisinds)
                for t = 1:floor(1000/timebin)
                    timerespmat = [timerespmat, sum(resp(:,thisinds(i),respwin(1)+(t-1)*timebin:respwin(1)+t*timebin),3)./(timebin/1000)];
                    timeoricond = [timeoricond, ori];
                    timeszcond = [timeszcond, sz];
                    timelcond = [timelcond, l];
                end
                for t = 1:floor(1000/corrbin)
                    corrbinresp(:,i,t) = sum(resp(:,thisinds(i),respwin(1)+round((t-1)*corrbin):respwin(1)+round(t*corrbin)),3)./(corrbin/1000);
                end
                for n1 = 1:size(resp,1)
                    for n2 = 1:size(resp,1)
                        cc(i,n1,n2) = nancorr(squeeze(corrbinresp(n1,i,:)),squeeze(corrbinresp(n2,i,:)));
                    end
                end                
            end
            condcorrs(l,ori,sz,:,:) = nanmean(cc,1);
            
            condz(:,l,ori,sz) = {(sc(:,thisinds)-repmat(mean(sc(:,thisinds),2),1,length(thisinds)))./repmat(std(sc(:,thisinds),1,2),1,length(thisinds))}; %ecker 2010
            condsc(:,l,ori,sz) = {sc(thisinds)};
            ff(:,l,ori,sz) = var(sc(:,thisinds),1,2)./mean(sc(:,thisinds),2);
            
            thisruninds = intersect(thisinds,oktrials);
            if ~isempty(thisruninds)
                runcondresp(:,l,ori,sz,:) = mean(resp(:,thisruninds,:),2);
                runcondfr(:,l,ori,sz) = mean(frs(:,thisruninds),2);
                runconderr(:,l,ori,sz) = std(frs(:,thisruninds),1,2)./sqrt(length(thisruninds));
            else
                runcondresp(:,l,ori,sz,:) = nan(size(resp,1),size(resp,3));
                runcondfr(:,l,ori,sz) = nan(1,size(resp,1));
                runconderr(:,l,ori,sz) = nan(1,size(resp,1));
            end
            
            thisstillinds = intersect(thisinds,stilltrials);
            if ~isempty(thisstillinds)
                stillcondresp(:,l,ori,sz,:) = mean(resp(:,thisstillinds,:),2);
                stillcondfr(:,l,ori,sz) = mean(frs(:,thisstillinds),2);
                stillconderr(:,l,ori,sz) = std(frs(:,thisstillinds),1,2)./sqrt(length(thisstillinds));
            else
                stillcondresp(:,l,ori,sz,:) = nan(size(resp,1),size(resp,3));
                stillcondfr(:,l,ori,sz) = nan(1,size(resp,1));
                stillconderr(:,l,ori,sz) = nan(1,size(resp,1));
            end
            
        end
    end
end

cond = szcond == 5 & lcond == 1;
badlines = [];
valmat = respmat(:,cond); conds = oricond(cond);
for i = 1:size(valmat,1)
    if ~(max(valmat(i,:)) == 0)
        valmat(i,:) = valmat(i,:); %./max(valmat(i,:));
    else
        badlines = [badlines,i];
    end
end
valmat(badlines,:) = [];
% if orientation do
conds = 1+mod(conds-1,max(conds)/2);
try
    [B,dev,stats] = mnrfit(valmat',conds');
    pihat = mnrval(B,valmat');


    [m,mpo] = max(pihat,[],2);
    for i = 1:size(pihat,2)
        for j = 1:size(pihat,2)
            cm(i,j) = length(find(conds == i & mpo' == j));
        end
        cmn(i,:) = cm(i,:)./length(find(conds == i));
    end
catch
    cmn = nan(length(unique(conds)));
end

clear cmn; clear cm;
for l = 1:2
    for sz = 1:5
        cond = szcond == sz & lcond == l;
        valmat = respmat(:,cond); conds = oricond(cond);
        badlines = [];
        for i = 1:size(valmat,1)
            if max(valmat(i,:)) == 0
                badlines = [badlines,i];
            end
        end
        if ~isempty(badlines)
            valmat(badlines,:) = [];
        end
        % for orientation, not direction        
        conds = 1+mod(conds-1,max(conds)/2);
        
        % cross validate
        tenp = floor(size(valmat,2)/10);
        rp = randperm(size(valmat,2));
        valmat = valmat(:,rp); conds = conds(rp);
        for f = 1:10 %size(valmat,2)
            testinds = (f-1)*tenp+1:(f-1)*tenp+tenp;
            traindat = valmat; traindat(:,testinds) = [];
            trainconds = conds; trainconds(testinds) = [];
            testdat = valmat(:,testinds);
            testconds = conds(testinds);

            try
                [B,dev,stats] = mnrfit(traindat',trainconds');

                pihat = mnrval(B,testdat');
                [m,mpo] = max(pihat,[],2);
                for i = 1:size(pihat,2)
                    for j = 1:size(pihat,2)
                        cm(f,i,j) = length(find(testconds == i & mpo' == j));
                    end
                    cmn(l,sz,f,i,:) = cm(f,i,:)./length(find(testconds == i));
                end
            catch
                cmn(l,sz,f,:,:) = nan(length(unique(conds)));
            end
        end
        percentcorrect(l,sz) = nanmean(diag(squeeze(nanmean(cmn(l,sz,:,:,:),3))));
        percentcorrecterr(l,sz) = nanmean(diag(squeeze(nanstd(cmn(l,sz,:,:,:),1,3))));
    end
end

% with time features
cond = timeszcond == 5 & timelcond == 1;
badlines = [];
valmat = timerespmat(:,cond); conds = timeoricond(cond);
for i = 1:size(valmat,1)
    if ~(max(valmat(i,:)) == 0)
        valmat(i,:) = valmat(i,:); %./max(valmat(i,:));
    else
        badlines = [badlines,i];
    end
end
valmat(badlines,:) = [];
% if orientation do
conds = 1+mod(conds-1,max(conds)/2);
try
    [B,dev,stats] = mnrfit(valmat',conds');
    pihat = mnrval(B,valmat');


    [m,mpo] = max(pihat,[],2);
    for i = 1:size(pihat,2)
        for j = 1:size(pihat,2)
            cm(i,j) = length(find(conds == i & mpo' == j));
        end
        cmn(i,:) = cm(i,:)./length(find(conds == i));
    end
catch
    cmn = nan(length(unique(conds)));
end

% all with time features

clear cmn; clear cm;
for l = 1:2
    for sz = 1:5
        cond = timeszcond == sz & timelcond == l;
        valmat = timerespmat(:,cond); conds = timeoricond(cond);
        badlines = [];
        for i = 1:size(valmat,1)
            if max(valmat(i,:)) == 0
                badlines = [badlines,i];
            end
        end
        if ~isempty(badlines)
            valmat(badlines,:) = [];
        end
        % for orientation, not direction        
%         conds = 1+mod(conds-1,max(conds)/2);
        
        % cross validate
        tenp = floor(size(valmat,2)/10);
        rp = randperm(size(valmat,2));
        valmat = valmat(:,rp); conds = conds(rp);
        for f = 1:10 %size(valmat,2)
            testinds = (f-1)*tenp+1:(f-1)*tenp+tenp;
            traindat = valmat; traindat(:,testinds) = [];
            trainconds = conds; trainconds(testinds) = [];
            testdat = valmat(:,testinds);
            testconds = conds(testinds);

            try
                [B,dev,stats] = mnrfit(traindat',trainconds');

                pihat = mnrval(B,testdat');
                [m,mpo] = max(pihat,[],2);
                for i = 1:size(pihat,2)
                    for j = 1:size(pihat,2)
                        cm(f,i,j) = length(find(testconds == i & mpo' == j));
                    end
                    cmn(l,sz,f,i,:) = cm(f,i,:)./length(find(testconds == i));
                end
            catch
                cmn(l,sz,f,:,:) = nan(length(unique(conds)));
            end
        end
        percentcorrect(l,sz) = nanmean(diag(squeeze(nanmean(cmn(l,sz,:,:,:),3))));
        percentcorrecterr(l,sz) = nanmean(diag(squeeze(nanstd(cmn(l,sz,:,:,:),1,3))));
    end
end


disp('');
