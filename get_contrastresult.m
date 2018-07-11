function co_result = get_contrastresult(filename)

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

co_result.animalid = result.animalid;
co_result.position = result.position;

prestim = 300;
poststim = 700;
if result.stimduration == 2
    respwin = 501:1500; % after stimulus onset
else
    respwin = 1:1000;
end
respwin = respwin+prestim;

[p,co_result.cellname] = fileparts(filename);

i = strfind(filename, 'tet');
co_result.tetno = strread(filename(i+3));

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
if isfield(result, 'contconds')
    allinds = sort(getSpecificIndices(result, 'contconds'));
    msstamps = result.msstamps(allinds);
    light = result.light(allinds);
    gratingInfo.Orientation = result.gratingInfo.Orientation(allinds);
    gratingInfo.Contrast = result.gratingInfo.Contrast(allinds);
    gratingInfo.tFreq = result.gratingInfo.tFreq(allinds);
    gratingInfo.size = result.gratingInfo.size(allinds);
else
    msstamps = result.msstamps;
    light = result.light;
    gratingInfo = result.gratingInfo;
    gratingInfo.size = ones(size(gratingInfo.Orientation)).*result.stimulusSize;
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

co_result.depth = result.depth;

spike = result.waveforms(:,wvchan);
co_result.interpspike = spline(1:32,spike,1:.1:32);
[co_result.adiff,co_result.swidth] = spikequant(co_result.interpspike);

msta = linspace(-prestim,trialdur+poststim,size(resp,2));

frs = sum(resp(:,respwin),2)./(length(respwin)/1000);
bl = sum(resp(:,1:prestim),2)./(prestim/1000);
sc = sum(resp(:,respwin),2);

%determine if cell is visually modulated
blfr = sum(resp(:,1:prestim),2);
vrfr = sum(resp(:,prestim+40:2*prestim+40),2);
co_result.vismod = signrank(blfr,vrfr);

%determine if cll is modulated by light
co_result.lightmod = signrank(frs(find(light)),frs(find(light == 0)));

co_result.lfr = mean(frs(find(light)));
co_result.nlfr = mean(frs(find(light == 0)));

binwidth = 33.33;

[p,printname] = fileparts(filename);
printname(find(printname=='_')) = ' ';
co_result.printname = printname;

anovavals = []; anovaori = []; anovaco = []; 
clevels = unique(gratingInfo.Contrast); clevels(find(clevels == 0)) = []; %delete control condition
oris = unique(gratingInfo.Orientation); oris(find(oris == -1)) = [];
sizes = unique(gratingInfo.size); sizes(find(sizes == 0)) = [];
for l = 1:length(unique(light))
    for co = 1:length(clevels)
        for sz = 1:length(sizes)
            for ori = 1:length(oris)
                thisinds = find(gratingInfo.Orientation == oris(ori) &...0
                    gratingInfo.Contrast == clevels(co) & ...
                    gratingInfo.size == sizes(sz) & ...
                    light == l-1);
                co_result.condn(l,ori,sz,co) = length(thisinds);
                co_result.condresp(l,ori,sz,co,:) = mean(resp(thisinds,:),1);
                co_result.condresperr(l,ori,sz,co,:) = nanstd(resp(thisinds,:),1,1)./sqrt(length(thisinds));
                co_result.allresp{l,ori,sz,co} = resp(thisinds,:);
                co_result.condlfpresp(l,ori,sz,co,:) = mean(lfpresp(thisinds,:),1);
                co_result.alllfpresp{l,ori,sz,co} = lfpresp(thisinds,:);
                co_result.condlfpspect(l,ori,sz,co,:) = nanmean(lfpspect(thisinds,:),1);

                co_result.condfr(l,ori,sz,co) = mean(frs(thisinds));%-mean(bl);
                co_result.conderr(l,ori,sz,co) =std(frs(thisinds))./sqrt(length(thisinds));
                if l == 1
                    anovavals = [anovavals; frs(thisinds)];
                    anovaori = [anovaori, ones(1,length(thisinds)).*ori];
                    anovaco = [anovaco, ones(1,length(thisinds)).*co];
                end

                co_result.condz(l,ori,sz,co) = {(sc(thisinds)-mean(sc(thisinds)))/std(sc(thisinds))}; %ecker 2010
                co_result.condsc(l,ori,sz,co) = {sc(thisinds)};
                co_result.ff(l,ori,sz,co) = var(sc(thisinds))/mean(sc(thisinds));

                binn = floor(size(co_result.condresp(l,ori,sz,co,:),4)/binwidth);
                co_result.bta = linspace(1,binn*binwidth,binn);
                if ~isnan(co_result.condresp(l,ori,sz,co,:))
                    [co_result.bincondresp(l,ori,sz,co,:)] = binit(co_result.condresp(l,ori,sz,co,:),binwidth);
                else
                    co_result.bincondresp(l,ori,sz,co,:) = binit(co_result.condresp(l,ori,sz,co,:),binwidth);
                end
                co_result.binconderr(l,ori,sz,co,:) = binit(co_result.condresperr(l,ori,sz,co,:),binwidth);            

                co_result.condsparseness(l,ori,sz,co) = sparseness(squeeze(co_result.bincondresp(l,ori,sz,co,:)));                    
                co_result.condsparsenesswin(l,ori,sz,co) = sparseness(squeeze(co_result.bincondresp(l,ori,sz,co,find(co_result.bta-prestim>respwin(1),1):find(co_result.bta-prestim>respwin(end),1)-1)));                        

                thisruninds = intersect(thisinds,oktrials);
                if ~isempty(thisruninds)
                    co_result.runcondresp(l,ori,sz,co,:) = mean(resp(thisruninds,:),1);
                    co_result.runcondfr(l,ori,sz,co) = mean(frs(thisruninds));
                    co_result.runconderr(l,ori,sz,co) = std(frs(thisruninds))./sqrt(length(thisruninds));
                    co_result.runn(l,ori,sz,co) = length(thisruninds);
                else
                    co_result.runcondresp(l,ori,sz,co,:) = nan(1,size(resp,2));
                    co_result.runcondfr(l,ori,sz,co) = NaN;
                    co_result.runconderr(l,ori,sz,co) = NaN;
                    co_result.runn(l,ori,sz,co) = 0;
                end

                thisstillinds = intersect(thisinds,stilltrials);
                if ~isempty(thisstillinds)
                    co_result.stillcondresp(l,ori,sz,co,:) = mean(resp(thisstillinds,:),1);
                    co_result.stillcondfr(l,ori,sz,co) = mean(frs(thisstillinds));
                    co_result.stillconderr(l,ori,sz,co) = std(frs(thisstillinds))./sqrt(length(thisstillinds));
                    co_result.stilln(l,ori,sz,co) = length(thisstillinds);
                else
                    co_result.stillcondresp(l,ori,sz,co,:) = nan(1,size(resp,2));
                    co_result.stillcondfr(l,ori,sz,co) = NaN;
                    co_result.stillconderr(l,ori,sz,co) = NaN;
                    co_result.stilln(l,ori,sz,co) = 0;
                end

                %get trial to trial reliability
                clear tc; clear lfptc; clear ftc;
                k = 1;
                if length(thisinds)>1
                    binrespwin = find(co_result.bta>respwin(1),1):find(co_result.bta>respwin(end),1)-1;
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
                    co_result.condrely(l,ori,sz,co) = nanmean(tc);
                    co_result.condrelyn(l,ori,sz,co) = length(find(resp(thisinds,:)));
                    co_result.condfiltrely(l,ori,sz,co) = nanmean(ftc);
                    co_result.condlfprely(l,ori,sz,co) = nanmean(lfptc);  
                else
                    co_result.condrely(l,ori,sz,co) = NaN;
                    co_result.condrelyn(l,ori,sz,co) = NaN;
                    co_result.condfiltrely(l,ori,sz,co) = NaN;
                    co_result.condlfprely(l,ori,sz,co) = NaN; 
                end
            end
        end
        [co_result.condoriprefratio(l,sz,co), co_result.conddirprefratio(l,sz,co), prefori, co_result.condmeanori(l,sz,co), co_result.condosi(l,sz,co), co_result.condmeandir(l,sz,co), co_result.conddsi(l,sz,co)] = getOSI(squeeze(co_result.condfr(l,:,sz,co)),oris);
        [co_result.gaussparams(l,sz,co,:), resnorm, co_result.gaussrsquared(l,sz,co)] = fit_doublegauss(oris,squeeze(co_result.condfr(l,:,sz,co)));
    end
end
co_result.clevels = clevels; co_result.oris = oris; co_result.sizes = sizes;
co_result.bincondresp = co_result.bincondresp.*(1000/binwidth);
co_result.bta = co_result.bta-prestim;
[p,table,stats] = anovan(anovavals,{anovaco anovaori},'model','full','display','off');
co_result.anova_cp = p(1); co_result.anova_op = p(2); co_result.anova_coip = p(3);


contindsnl = find(gratingInfo.Contrast == 0 & light == 0);
co_result.controlresp(1,:) = mean(resp(contindsnl,:),1);
co_result.allcontrolresp = resp(contindsnl,:);
co_result.controlfr(1) = mean(frs(contindsnl));
co_result.controlerr(1) = std(frs(contindsnl))./sqrt(length(contindsnl));
co_result.controllfpspect(1,:) = nanmean(lfpspect(contindsnl,:));

contindsl = find(gratingInfo.Contrast == 0 & light == 1);
co_result.controlresp(2,:) = mean(resp(contindsl,:),1);
co_result.controlfr(2) = mean(frs(contindsl));
co_result.controlerr(2) = std(frs(contindsl))./sqrt(length(contindsl));
co_result.controllfpspect(2,:) = nanmean(lfpspect(contindsl,:));

[co_result.binnedctr(1,:),bta] = binit(co_result.controlresp(1,:),binwidth);
[co_result.binnedctr(2,:),bta] = binit(co_result.controlresp(2,:),binwidth);

l0r0continds = intersect(contindsnl,stilltrials);
l0r1continds = intersect(contindsnl,oktrials);
l1r0continds = intersect(contindsl, stilltrials);
l1r1continds = intersect(contindsl, oktrials);
if ~isempty(l0r0continds) % first dim light
    co_result.r0contolresp(1,:) = mean(resp(l0r0continds,:),1);
    co_result.r0controlfr(1) = mean(frs(l0r0continds));
    co_result.r0controlerr(1) = std(frs(l0r0continds))./sqrt(length(l0r0continds));
    co_result.r0controln(1) = length(l0r0continds);
else
    co_result.r0contolresp(1,:) = nan(1,size(resp,2));
    co_result.r0controlfr(1) = NaN;
    co_result.r0controlerr(1) = NaN;
    co_result.r0controln(1) = 0;
end
if ~isempty(l1r0continds) % first dim light
    co_result.r0contolresp(2,:) = mean(resp(l1r0continds,:),1);
    co_result.r0controlfr(2) = mean(frs(l1r0continds));
    co_result.r0controlerr(2) = std(frs(l1r0continds))./sqrt(length(l1r0continds));
    co_result.r0controln(2) = length(l1r0continds);
else
    co_result.r0contolresp(2,:) = nan(1,size(resp,2));
    co_result.r0controlfr(2) = NaN;
    co_result.r0controlerr(2) = NaN;
    co_result.r0controln(2) = 0;
end
if ~isempty(l0r1continds) % first dim light
    co_result.r1contolresp(1,:) = mean(resp(l0r1continds,:),1);
    co_result.r1controlfr(1) = mean(frs(l0r1continds));
    co_result.r1controlerr(1) = std(frs(l0r1continds))./sqrt(length(l0r1continds));
    co_result.r1controln(1) = length(l0r1continds);
else
    co_result.r1contolresp(1,:) = nan(1,size(resp,2));
    co_result.r1controlfr(1) = NaN;
    co_result.r1controlerr(1) = NaN;
    co_result.r1controln(1) = 0;
end
if ~isempty(l1r1continds) % first dim light
    co_result.r1contolresp(2,:) = mean(resp(l1r1continds,:),1);
    co_result.r1controlfr(2) = mean(frs(l1r1continds));
    co_result.r1controlerr(2) = std(frs(l1r1continds))./sqrt(length(l1r1continds));
    co_result.r1controln(2) = length(l1r1continds);
else
    co_result.r1contolresp(2,:) = nan(1,size(resp,2));
    co_result.r1controlfr(2) = NaN;
    co_result.r1controlerr(2) = NaN;
    co_result.r1controln(2) = 0;
end

co_result.prefcontrast = find(mean(mean(co_result.condfr(1,:,:,:),2),3) == max(mean(mean(co_result.condfr(1,:,:,:),2),3)),1);

co_result.xlevels = [0,clevels];
for sz = 1:length(sizes)
    co_result.cresp(1,sz,:) = [co_result.controlfr(1), squeeze(nanmean(co_result.condfr(1,:,sz,:),2))'];
    co_result.cresp(2,sz,:) = [co_result.controlfr(2), squeeze(nanmean(co_result.condfr(2,:,sz,:),2))'];
    [co_result.paramsl0(sz,:), bu, co_result.rsql0(sz,:)] = fit_crf_NR(co_result.xlevels,co_result.cresp(1,sz,:));
    [co_result.paramsl1(sz,:), bu, co_result.rsql1(sz,:)] = fit_crf_NR(co_result.xlevels,co_result.cresp(2,sz,:));
end


% test for running modulation
l1r1 = frs(intersect(find(light),oktrials));
l0r1 = frs(intersect(find(~light),oktrials));
l1r0 = frs(intersect(find(light),stilltrials));
l0r0 = frs(intersect(find(~light),stilltrials));
anovavec = [l0r0;l0r1;l1r0;l1r1];
g1 = [zeros(length(l0r0),1);zeros(length(l0r1),1);ones(length(l1r0),1);ones(length(l1r1),1)]; %light
g2 = [zeros(length(l0r0),1);ones(length(l0r1),1);zeros(length(l1r0),1);ones(length(l1r1),1)]; %running
[p,table,stats] = anovan(anovavec,{g1 g2},'model','full','display','off');
co_result.anova_lp = p(1); co_result.anova_rp = p(2); co_result.anova_rlip = p(3);
end

function [params, resnorm, rsquared] = fit_crf_NR(x,y)
% [rmax, n, c50, r0]

nrfitmargin = .01; % .5
range = max(y)-min(y);
if range ~= 0
    p0 = [range 2 .5 min(y)];
    lb = [(1-nrfitmargin)*range, .5, 0, min(y)-.1*range];
    ub = [(1+nrfitmargin)*range, 10, 1, min(y)+.1*range];
    %         lb = [-max(y), .5, 0, -10*max(y)];
    %         ub = [10*max(y), 10, 1, max(y)];
    warning off
    [params,resnorm,residual,exitflag] = lsqcurvefit(@(p,x) NakaRushton(p,x),p0,x(:),y(:),lb,ub,optimset('Display','off'));
    warning on
    restot = sum((y-mean(y)).^2);
    rsquared = 1 - (resnorm/restot);
else
    params = [NaN,NaN,NaN,NaN];
    resnorm = NaN;
    residual = NaN;
    exitflag = NaN;
    rsquared = NaN;
end
end

function val=NakaRushton(p,x)
% parameters of Naka-Rushton function as in Disney et al., Neuron, 2007
% [R_max, contrast Exponent n,  50%firing-Contrast, spontaneous rate sFR]

val = p(4)+p(1)*((x.^p(2))./(x.^p(2)+p(3).^p(2)));
end

