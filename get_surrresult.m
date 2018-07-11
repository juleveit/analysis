function sr_result = get_surrresult(filename)

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

sr_result.animalid = result.animalid;
sr_result.position = result.position;

prestim = 300;
poststim = 700;
if result.stimduration == 2
    respwin = 501:1500; % after stimulus onset
else
    respwin = 1:1000;
end
respwin = respwin+prestim;

[p,sr_result.cellname] = fileparts(filename);

i = strfind(filename, 'tet');
sr_result.tetno = strread(filename(i+3));

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

msstamps = result.msstamps;
light = result.light;
gratingInfo = result.gratingInfo;

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

sr_result.depth = result.depth;

spike = result.waveforms(:,wvchan);
sr_result.interpspike = spline(1:32,spike,1:.1:32);
[sr_result.adiff,sr_result.swidth] = spikequant(sr_result.interpspike);

msta = linspace(-prestim,trialdur+poststim,size(resp,2));

frs = sum(resp(:,respwin),2)./(length(respwin)/1000);
bl = sum(resp(:,1:prestim),2)./(prestim/1000);
sc = sum(resp(:,respwin),2);

%determine if cell is visually modulated
blfr = sum(resp(:,1:prestim),2);
vrfr = sum(resp(:,prestim+40:2*prestim+40),2);
sr_result.vismod = signrank(blfr,vrfr);

%determine if cll is modulated by light
sr_result.lightmod = signrank(frs(find(light)),frs(find(light == 0)));

sr_result.lfr = mean(frs(find(light)));
sr_result.nlfr = mean(frs(find(light == 0)));

binwidth = 33.33;

[p,printname] = fileparts(filename);
printname(find(printname=='_')) = ' ';
sr_result.printname = printname;

anovavals = []; anovaori = []; anovad = []; 
centeroris = unique(result.gratingInfo.Orientation_center);  centeroris(find(centeroris == -1)) = []; %delete control condition
surroundoris = unique(result.gratingInfo.Orientation_surround); surroundoris(find(surroundoris == -1)) = []; 
ordf = [0,90,-1,-2];  % 1: iso, 2: cross, 3: center only, 4: surround only
for l = 1:length(unique(result.light))
    for ori = 1:length(centeroris)
        for d = 1:length(ordf)
            if ordf(d) == -1
                thisinds = find(result.gratingInfo.Orientation_center == centeroris(ori) & ...
                    result.gratingInfo.Orientation_surround == -1 & ...
                    result.light == l-1);
            elseif ordf(d) == -2
                thisinds = find(result.gratingInfo.Orientation_center == -1 & ...
                    result.gratingInfo.Orientation_surround == centeroris(ori) & ...
                    result.light == l-1);
            else
                thisinds = find(result.gratingInfo.Orientation_center == centeroris(ori) &...
                    result.gratingInfo.Orientation_surround == centeroris(ori)+ordf(d) & ...
                    result.light == l-1);
            end
            sr_result.condn(l,ori,d) = length(thisinds);
            sr_result.condresp(l,ori,d,:) = mean(resp(thisinds,:),1);
            sr_result.condresperr(l,ori,d,:) = nanstd(resp(thisinds,:),1,1)./sqrt(length(thisinds));
            sr_result.allresp{l,ori,d} = resp(thisinds,:);
            sr_result.condlfpresp(l,ori,d,:) = mean(lfpresp(thisinds,:),1);
            sr_result.alllfpresp{l,ori,d} = lfpresp(thisinds,:);
            sr_result.condlfpspect(l,ori,d,:) = nanmean(lfpspect(thisinds,:),1);
            sr_result.condlfpspecterr(l,ori,d,:) = nanstd(lfpspect(thisinds,:),1)./sqrt(length(thisinds));            
            sr_result.condfr(l,ori,d) = mean(frs(thisinds));%-mean(bl);
            sr_result.conderr(l,ori,d) =std(frs(thisinds))./sqrt(length(thisinds));
            if l == 1
                anovavals = [anovavals; frs(thisinds)];
                anovaori = [anovaori, ones(1,length(thisinds)).*ori];
                anovad = [anovad, ones(1,length(thisinds)).*d];
            end
            
            sr_result.condz(l,ori,d) = {(sc(thisinds)-mean(sc(thisinds)))/std(sc(thisinds))}; %ecker 2010
            sr_result.condsc(l,ori,d) = {sc(thisinds)};
            sr_result.ff(l,ori,d) = var(sc(thisinds))/mean(sc(thisinds));
                            
            binn = floor(size(sr_result.condresp(l,ori,d,:),4)/binwidth);
            sr_result.bta = linspace(1,binn*binwidth,binn);
            if ~isnan(sr_result.condresp(l,ori,d,:))
                [sr_result.bincondresp(l,ori,d,:)] = binit(sr_result.condresp(l,ori,d,:),binwidth);
            else
                sr_result.bincondresp(l,ori,d,:) = binit(sr_result.condresp(l,ori,d,:),binwidth);
            end
            sr_result.binconderr(l,ori,d,:) = binit(sr_result.condresperr(l,ori,d,:),binwidth);            
            
            sr_result.condsparseness(l,ori,d) = sparseness(squeeze(sr_result.bincondresp(l,ori,d,:)));                    
            sr_result.condsparsenesswin(l,ori,d) = sparseness(squeeze(sr_result.bincondresp(l,ori,d,find(sr_result.bta-prestim>respwin(1),1):find(sr_result.bta-prestim>respwin(end),1)-1)));                        
            
            thisruninds = intersect(thisinds,oktrials);
            if ~isempty(thisruninds)
                sr_result.runcondresp(l,ori,d,:) = mean(resp(thisruninds,:),1);
                sr_result.runcondfr(l,ori,d) = mean(frs(thisruninds));
                sr_result.runconderr(l,ori,d) = std(frs(thisruninds))./sqrt(length(thisruninds));
                sr_result.runn(l,ori,d) = length(thisruninds);
            else
                sr_result.runcondresp(l,ori,d,:) = nan(1,size(resp,2));
                sr_result.runcondfr(l,ori,d) = NaN;
                sr_result.runconderr(l,ori,d) = NaN;
                sr_result.runn(l,ori,d) = 0;
            end
            
            thisstillinds = intersect(thisinds,stilltrials);
            if ~isempty(thisstillinds)
                sr_result.stillcondresp(l,ori,d,:) = mean(resp(thisstillinds,:),1);
                sr_result.stillcondfr(l,ori,d) = mean(frs(thisstillinds));
                sr_result.stillconderr(l,ori,d) = std(frs(thisstillinds))./sqrt(length(thisstillinds));
                sr_result.stilln(l,ori,d) = length(thisstillinds);
            else
                sr_result.stillcondresp(l,ori,d,:) = nan(1,size(resp,2));
                sr_result.stillcondfr(l,ori,d) = NaN;
                sr_result.stillconderr(l,ori,d) = NaN;
                sr_result.stilln(l,ori,d) = 0;
            end
            
            %get trial to trial reliability
            clear tc; clear lfptc; clear ftc;
            k = 1;
            if length(thisinds)>1
                binrespwin = find(sr_result.bta>respwin(1),1):find(sr_result.bta>respwin(end),1)-1;
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
                sr_result.condrely(l,ori,d) = nanmean(tc);
                sr_result.condrelyn(l,ori,d) = length(find(resp(thisinds,:)));
                sr_result.condfiltrely(l,ori,d) = nanmean(ftc);
                sr_result.condlfprely(l,ori,d) = nanmean(lfptc);  
            else
                sr_result.condrely(l,ori,d) = NaN;
                sr_result.condrelyn(l,ori,d) = NaN;
                sr_result.condfiltrely(l,ori,d) = NaN;
                sr_result.condlfprely(l,ori,d) = NaN; 
            end
        end
    end
end
sr_result.ordf = ordf; sr_result.centeroris = centeroris;
sr_result.bincondresp = sr_result.bincondresp.*(1000/binwidth);
sr_result.bta = sr_result.bta-prestim;
sr_result.fax = trialfax;
[p,table,stats] = anovan(anovavals,{anovad anovaori},'model','full','display','off');
sr_result.anova_dp = p(1); sr_result.anova_op = p(2); sr_result.anova_doip = p(3);

% just for lfp
for l = 1:length(unique(result.light))
    for d = 1:length(ordf)
        if ordf(d) == -1
            thisinds = find(result.gratingInfo.Orientation_surround == -1 & ...
                result.light == l-1);
        elseif ordf(d) == -2
            thisinds = find(result.gratingInfo.Orientation_center == -1 & ...
                result.light == l-1);
        else
            thisinds = find(result.gratingInfo.Orientation_center ~= -1 &...
                result.gratingInfo.Orientation_surround == gratingInfo.Orientation_center+ordf(d) & ...
                result.light == l-1);
        end        
        sr_result.omlfpspect(l,d,:) = nanmean(lfpspect(thisinds,:),1); % orientation mean
        sr_result.omlfpspecterr(l,d,:) = nanstd(lfpspect(thisinds,:),1)./sqrt(length(thisinds));
    end
end

contindsnl = find(result.gratingInfo.Orientation_center == -1 & result.gratingInfo.Orientation_surround == -1 & light == 0);
sr_result.controlresp(1,:) = mean(resp(contindsnl,:),1);
sr_result.allcontrolresp = resp(contindsnl,:);
sr_result.controlfr(1) = mean(frs(contindsnl));
sr_result.controlerr(1) = std(frs(contindsnl))./sqrt(length(contindsnl));
sr_result.controllfpspect(1,:) = nanmean(lfpspect(contindsnl,:));

contindsl = find(result.gratingInfo.Orientation_center == -1 & result.gratingInfo.Orientation_surround == -1 & light == 1);
sr_result.controlresp(2,:) = mean(resp(contindsl,:),1);
sr_result.controlfr(2) = mean(frs(contindsl));
sr_result.controlerr(2) = std(frs(contindsl))./sqrt(length(contindsl));
sr_result.controllfpspect(2,:) = nanmean(lfpspect(contindsl,:));

[sr_result.binnedctr(1,:),bta] = binit(sr_result.controlresp(1,:),binwidth);
[sr_result.binnedctr(2,:),bta] = binit(sr_result.controlresp(2,:),binwidth);

sr_result.sizecenter = unique(result.gratingInfo.sizecenter);
sr_result.sizesurround = unique(result.gratingInfo.sizesurround);

% test for running modulation
l1r1 = frs(intersect(find(light),oktrials));
l0r1 = frs(intersect(find(~light),oktrials));
l1r0 = frs(intersect(find(light),stilltrials));
l0r0 = frs(intersect(find(~light),stilltrials));
anovavec = [l0r0;l0r1;l1r0;l1r1];
g1 = [zeros(length(l0r0),1);zeros(length(l0r1),1);ones(length(l1r0),1);ones(length(l1r1),1)]; %light
g2 = [zeros(length(l0r0),1);ones(length(l0r1),1);zeros(length(l1r0),1);ones(length(l1r1),1)]; %running
[p,table,stats] = anovan(anovavec,{g1 g2},'model','full','display','off');
sr_result.anova_lp = p(1); sr_result.anova_rp = p(2); sr_result.anova_rlip = p(3);

end

