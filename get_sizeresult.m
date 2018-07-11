function sz_result = get_sizeresult(filename)

% chronux parameters
params.tapers = [5,9]; params.Fs = 1000; params.err = [2, 0.05]; params.trialave = 1;

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

load([filename]);

sz_result.animalid = result.animalid;
sz_result.position = result.position;

prestim = 300;
poststim = 700;
if result.stimduration == 2
    respwin = 501:1500; % after stimulus onset
else
    respwin = 1:1000;
end
respwin = respwin+prestim;

[p,sz_result.cellname] = fileparts(filename);

i = strfind(filename, 'tet');
sz_result.tetno = strread(filename(i+3));

wvchan = find(var(result.waveforms) == max(var(result.waveforms)));

sr = 1000;
lfp = result.lfp(:,wvchan)';

msStimes = round(result.spikes);
if isempty(msStimes), msStimes(1) = 0; end
if msStimes(1) == 0, msStimes(1) = 1; end

chan = zeros(1,length(result.lfp));
chan(msStimes) = 1;

trialdur = result.stimduration*1000;
msstamps = result.msstamps;

if length(msstamps)~=length(result.light)
    %          disp('');
    %         msstamps([62,108,147]) = []; % for 140703 block 8
    %         msstamps([161]) = []; % for 141204 block 3
    %         msstamps([303]) = []; % for 150407 block 5
    %         msstamps([169,336]) = []; % for 150523 block 11
    %         msstamps([24]) = []; % for 150730 block 11
    %         msstamps([207]) = []; % for 150730 block 11
    %         msstamps([242]) = []; % for 151210 block 4
    %         msstamps([412]) = []; % for 151210 block 3
    %         msstamps([83]) = []; % for 160125 block 7
    %         msstamps([302,340]) = []; % for 160328 block 2
    %         msstamps([518]) = []; % for 160726 block 6
    %         msstamps([390,631,875]) = []; % for 170426 block 2
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

for i = 1:length(msstamps)
    resp(i,:) = chan(msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim);
    lfpresp(i,:) = result.lfp(msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim, wvchan);
    [lfpspect(i,:),trialfax] = mtspectrumc(squeeze(lfpresp(i,1001:1800))',params);
    
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

sz_result.depth = result.depth;

spike = result.waveforms(:,wvchan);
sz_result.interpspike = spline(1:32,spike,1:.1:32);
[sz_result.adiff,sz_result.swidth] = spikequant(sz_result.interpspike);

msta = linspace(-prestim,trialdur+poststim,size(resp,2));

frs = sum(resp(:,respwin),2)./(length(respwin)/1000);
bl = sum(resp(:,1:prestim),2)./(prestim/1000);
sc = sum(resp(:,respwin),2);

%determine if cell is visually modulated
blfr = sum(resp(:,1:prestim),2);
vrfr = sum(resp(:,prestim+40:2*prestim+40),2);
sz_result.vismod = signrank(blfr,vrfr);

%determine if cll is modulated by light
sz_result.lightmod = signrank(frs(find(light)),frs(find(light == 0)));

sz_result.lfr = mean(frs(find(light)));
sz_result.nlfr = mean(frs(find(light == 0)));

binwidth = 33.33;

[p,printname] = fileparts(filename);
printname(find(printname=='_')) = ' ';
sz_result.printname = printname;

anovavals = []; anovaori = []; anovasz = [];
sizes = unique(gratingInfo.size);  sizes(find(sizes == 0)) = []; %delete control condition
oris = unique(gratingInfo.Orientation); oris(find(oris == -1)) = [];
for l = 1:length(unique(light))
    for sz = 1:length(sizes)
        for ori = 1:length(oris)
            thisinds = find(gratingInfo.Orientation == oris(ori) &...0
                gratingInfo.size == sizes(sz) & ...
                light == l-1);
            sz_result.condn(l,ori,sz) = length(thisinds);
            sz_result.condresp(l,ori,sz,:) = mean(resp(thisinds,:),1);
            sz_result.condresperr(l,ori,sz,:) = nanstd(resp(thisinds,:),1,1)./sqrt(length(thisinds));
            sz_result.allresp{l,ori,sz} = resp(thisinds,:);
            sz_result.condlfpresp(l,ori,sz,:) = mean(lfpresp(thisinds,:),1);
            sz_result.alllfpresp{l,ori,sz} = lfpresp(thisinds,:);
            sz_result.condlfpspect(l,ori,sz,:) = nanmean(lfpspect(thisinds,:),1);
            
            sz_result.condfr(l,ori,sz) = mean(frs(thisinds));%-mean(bl);
            sz_result.conderr(l,ori,sz) =std(frs(thisinds))./sqrt(length(thisinds));
            
            if l == 1
                anovavals = [anovavals; frs(thisinds)];
                anovaori = [anovaori, ones(1,length(thisinds)).*ori];
                anovasz = [anovasz, ones(1,length(thisinds)).*sz];
            end
            
            sz_result.condz(l,ori,sz) = {(sc(thisinds)-mean(sc(thisinds)))/std(sc(thisinds))}; %ecker 2010
            sz_result.condsc(l,ori,sz) = {sc(thisinds)};
            sz_result.ff(l,ori,sz) = var(sc(thisinds))/mean(sc(thisinds));
                           
            binn = floor(size(sz_result.condresp(l,ori,sz,:),4)/binwidth);
            sz_result.bta = linspace(1,binn*binwidth,binn);
            if ~isnan(sz_result.condresp(l,ori,sz,:))
                [sz_result.bincondresp(l,ori,sz,:)] = binit(sz_result.condresp(l,ori,sz,:),binwidth);
            else
                sz_result.bincondresp(l,ori,sz,:) = binit(sz_result.condresp(l,ori,sz,:),binwidth);
            end
            sz_result.binconderr(l,ori,sz,:) = binit(sz_result.condresperr(l,ori,sz,:),binwidth);            
            
            sz_result.condsparseness(l,ori,sz) = sparseness(squeeze(sz_result.bincondresp(l,ori,sz,:)));                    
            sz_result.condsparsenesswin(l,ori,sz) = sparseness(squeeze(sz_result.bincondresp(l,ori,sz,find(sz_result.bta-prestim>respwin(1),1):find(sz_result.bta-prestim>respwin(end),1)-1)));                        
            
            thisruninds = intersect(thisinds,oktrials);
            if ~isempty(thisruninds)
                sz_result.runcondresp(l,ori,sz,:) = mean(resp(thisruninds,:),1);
                sz_result.runcondfr(l,ori,sz) = mean(frs(thisruninds));
                sz_result.runconderr(l,ori,sz) = std(frs(thisruninds))./sqrt(length(thisruninds));
                sz_result.runn(l,ori,sz) = length(thisruninds);
            else
                sz_result.runcondresp(l,ori,sz,:) = nan(1,size(resp,2));
                sz_result.runcondfr(l,ori,sz) = NaN;
                sz_result.runconderr(l,ori,sz) = NaN;
                sz_result.runn(l,ori,sz) = 0;
            end
            
            thisstillinds = intersect(thisinds,stilltrials);
            if ~isempty(thisstillinds)
                sz_result.stillcondresp(l,ori,sz,:) = mean(resp(thisstillinds,:),1);
                sz_result.stillcondfr(l,ori,sz) = mean(frs(thisstillinds));
                sz_result.stillconderr(l,ori,sz) = std(frs(thisstillinds))./sqrt(length(thisstillinds));
                sz_result.stilln(l,ori,sz) = length(thisstillinds);
            else
                sz_result.stillcondresp(l,ori,sz,:) = nan(1,size(resp,2));
                sz_result.stillcondfr(l,ori,sz) = NaN;
                sz_result.stillconderr(l,ori,sz) = NaN;
                sz_result.stilln(l,ori,sz) = 0;
            end
            
            %get trial to trial reliability
            clear tc; clear lfptc; clear ftc;
            k = 1;
            if length(thisinds)>1
                binrespwin = find(sz_result.bta>respwin(1),1):find(sz_result.bta>respwin(end),1)-1;
                for i = 1:length(thisinds)-1
                    for j = i+1:length(thisinds)
                        a = binit(resp(thisinds(i),:),binwidth);
                        b = binit(resp(thisinds(j),:),binwidth);
                        tc(k) = nancorr(a(binrespwin),b(binrespwin));

                        a = binit(lfpresp(thisinds(i),:),binwidth);
                        b = binit(lfpresp(thisinds(j),:),binwidth);
                        lfptc(k) = nancorr(a(binrespwin),b(binrespwin));

                        a = filter(excit_kernel,1,resp(thisinds(i),:));
                        b = filter(excit_kernel,1,resp(thisinds(j),:));
                        ftc(k) = nancorr(a(respwin),b(respwin));

                        k = k+1;
                    end
                end
                sz_result.condrely(l,ori,sz) = nanmean(tc);
                sz_result.condrelyn(l,ori,sz) = length(find(resp(thisinds,:)));
                sz_result.condfiltrely(l,ori,sz) = nanmean(ftc);
                sz_result.condlfprely(l,ori,sz) = nanmean(lfptc);  
            else
                sz_result.condrely(l,ori,sz) = NaN;
                sz_result.condrelyn(l,ori,sz) = NaN;
                sz_result.condfiltrely(l,ori,sz) = NaN;
                sz_result.condlfprely(l,ori,sz) = NaN;  
            end
        end
        [sz_result.condoriprefratio(l,sz), sz_result.conddirprefratio(l,sz), prefori, sz_result.condmeanori(l,sz), sz_result.condosi(l,sz), sz_result.condmeandir(l,sz), sz_result.conddsi(l,sz)] = getOSI(squeeze(sz_result.condfr(l,:,sz)),oris);
        [sz_result.gaussparams(l,sz,:), resnorm, sz_result.gaussrsquared(l,sz)] = fit_doublegauss(oris,squeeze(sz_result.condfr(l,:,sz)));
    end
end
sz_result.sizes = sizes; sz_result.oris = oris;
sz_result.bincondresp = sz_result.bincondresp.*(1000/binwidth);
sz_result.bta = sz_result.bta-prestim;

%only for LFP
for l = 1:length(unique(light))
    for sz = 1:length(sizes)        
        thisinds = find(gratingInfo.size == sizes(sz) & light == l-1);
        sz_result.omlfpspect(l,sz,:) = nanmean(lfpspect(thisinds,:),1);
        sz_result.omlfpspecterr(l,sz,:) = nanstd(lfpspect(thisinds,:),1)./sqrt(length(thisinds));
    end
end

[p,table,stats] = anovan(anovavals,{anovasz anovaori},'model','full','display','off');
sz_result.anova_sp = p(1); sz_result.anova_op = p(2); sz_result.anova_soip = p(3);

contindsnl = find(gratingInfo.size == 0 & light == 0);
sz_result.controlresp(1,:) = mean(resp(contindsnl,:),1);
sz_result.allcontrolresp = resp(contindsnl,:);
sz_result.controlfr(1) = mean(frs(contindsnl));
sz_result.controlerr(1) = std(frs(contindsnl))./sqrt(length(contindsnl));
sz_result.controllfpspect(1,:) = nanmean(lfpspect(contindsnl,:));

contindsl = find(gratingInfo.size == 0 & light == 1);
sz_result.controlresp(2,:) = mean(resp(contindsl,:),1);
sz_result.controlfr(2) = mean(frs(contindsl));
sz_result.controlerr(2) = std(frs(contindsl))./sqrt(length(contindsl));
sz_result.controllfpspect(2,:) = nanmean(lfpspect(contindsl,:));

[sz_result.binnedctr(1,:),bta] = binit(sz_result.controlresp(1,:),binwidth);
[sz_result.binnedctr(2,:),bta] = binit(sz_result.controlresp(2,:),binwidth);

l0r0continds = intersect(contindsnl,stilltrials);
l0r1continds = intersect(contindsnl,oktrials);
l1r0continds = intersect(contindsl, stilltrials);
l1r1continds = intersect(contindsl, oktrials);
if ~isempty(l0r0continds) % first dim light
    sz_result.r0contolresp(1,:) = mean(resp(l0r0continds,:),1);
    sz_result.r0controlfr(1) = mean(frs(l0r0continds));
    sz_result.r0controlerr(1) = std(frs(l0r0continds))./sqrt(length(l0r0continds));
    sz_result.r0controln(1) = length(l0r0continds);
else
    sz_result.r0contolresp(1,:) = nan(1,size(resp,2));
    sz_result.r0controlfr(1) = NaN;
    sz_result.r0controlerr(1) = NaN;
    sz_result.r0controln(1) = 0;
end
if ~isempty(l1r0continds) % first dim light
    sz_result.r0contolresp(2,:) = mean(resp(l1r0continds,:),1);
    sz_result.r0controlfr(2) = mean(frs(l1r0continds));
    sz_result.r0controlerr(2) = std(frs(l1r0continds))./sqrt(length(l1r0continds));
    sz_result.r0controln(2) = length(l1r0continds);
else
    sz_result.r0contolresp(2,:) = nan(1,size(resp,2));
    sz_result.r0controlfr(2) = NaN;
    sz_result.r0controlerr(2) = NaN;
    sz_result.r0controln(2) = 0;
end
if ~isempty(l0r1continds) % first dim light
    sz_result.r1contolresp(1,:) = mean(resp(l0r1continds,:),1);
    sz_result.r1controlfr(1) = mean(frs(l0r1continds));
    sz_result.r1controlerr(1) = std(frs(l0r1continds))./sqrt(length(l0r1continds));
    sz_result.r1controln(1) = length(l0r1continds);
else
    sz_result.r1contolresp(1,:) = nan(1,size(resp,2));
    sz_result.r1controlfr(1) = NaN;
    sz_result.r1controlerr(1) = NaN;
    sz_result.r1controln(1) = 0;
end
if ~isempty(l1r1continds) % first dim light
    sz_result.r1contolresp(2,:) = mean(resp(l1r1continds,:),1);
    sz_result.r1controlfr(2) = mean(frs(l1r1continds));
    sz_result.r1controlerr(2) = std(frs(l1r1continds))./sqrt(length(l1r1continds));
    sz_result.r1controln(2) = length(l1r1continds);
else
    sz_result.r1contolresp(2,:) = nan(1,size(resp,2));
    sz_result.r1controlfr(2) = NaN;
    sz_result.r1controlerr(2) = NaN;
    sz_result.r1controln(2) = 0;
end

sz_result.prefsize = find(mean(sz_result.condfr(1,:,:),2) == max(mean(sz_result.condfr(1,:,:),2)),1);

% test for running modulation
l1r1 = frs(intersect(find(light),oktrials));
l0r1 = frs(intersect(find(~light),oktrials));
l1r0 = frs(intersect(find(light),stilltrials));
l0r0 = frs(intersect(find(~light),stilltrials));
anovavec = [l0r0;l0r1;l1r0;l1r1];
g1 = [zeros(length(l0r0),1);zeros(length(l0r1),1);ones(length(l1r0),1);ones(length(l1r1),1)]; %light
g2 = [zeros(length(l0r0),1);ones(length(l0r1),1);zeros(length(l1r0),1);ones(length(l1r1),1)]; %running
[p,table,stats] = anovan(anovavec,{g1 g2},'model','full','display','off');
sz_result.anova_lp = p(1); sz_result.anova_rp = p(2); sz_result.anova_rlip = p(3);


end
