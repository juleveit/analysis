function ph_result = get_phasresult(filename)

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

ph_result.animalid = result.animalid;
ph_result.position = result.position;

prestim = 300;
poststim = 700;
if result.stimduration == 2
    respwin = 501:1500; % after stimulus onset
else
    respwin = 1:1000;
end
respwin = respwin+prestim;

[p,ph_result.cellname] = fileparts(filename);

i = strfind(filename, 'tet');
ph_result.tetno = strread(filename(i+3));

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

ph_result.depth = result.depth;

spike = result.waveforms(:,wvchan);
ph_result.interpspike = spline(1:32,spike,1:.1:32);
[ph_result.adiff,ph_result.swidth] = spikequant(ph_result.interpspike);

msta = linspace(-prestim,trialdur+poststim,size(resp,2));

frs = sum(resp(:,respwin),2)./(length(respwin)/1000);
bl = sum(resp(:,1:prestim),2)./(prestim/1000);
sc = sum(resp(:,respwin),2);

%determine if cell is visually modulated
blfr = sum(resp(:,1:prestim),2);
vrfr = sum(resp(:,prestim+40:2*prestim+40),2);
ph_result.vismod = signrank(blfr,vrfr);

%determine if cll is modulated by light
ph_result.lightmod = signrank(frs(find(light)),frs(find(light == 0)));

ph_result.lfr = mean(frs(find(light)));
ph_result.nlfr = mean(frs(find(light == 0)));

binwidth = 33.33;

[p,printname] = fileparts(filename);
printname(find(printname=='_')) = ' ';
ph_result.printname = printname;

anovavals = []; anovaori = []; anovad = []; 

oris = unique(result.gratingInfo.Orientation);  oris(find(oris == -1)) = []; %delete control condition
surroundphases = unique(result.gratingInfo.Phase_surround); surroundphases(find(surroundphases == -1)) = []; 
ordf = [0,180,-1,-2];  % 1: iso, 2: 180off, 3: center only, 4: surround only
for l = 1:length(unique(result.light))
    for ori = 1:length(oris)
        for d = 1:length(ordf)
            if ordf(d) == -1
                thisinds = find(result.gratingInfo.Orientation == oris(ori) & ...
                    result.gratingInfo.Phase_surround == -1 & ...
                    result.light == l-1);
            elseif ordf(d) == -2
                thisinds = find(result.gratingInfo.Orientation == -1 & ...
                    result.gratingInfo.Phase_surround == oris(ori) & ...
                    result.light == l-1);
            else
                thisinds = find(result.gratingInfo.Orientation == oris(ori) &...
                    result.gratingInfo.Phase_surround == ordf(d) & ...
                    result.light == l-1);
            end
            ph_result.condn(l,ori,d) = length(thisinds);
            ph_result.condresp(l,ori,d,:) = nanmean(resp(thisinds,:),1);
            ph_result.condresperr(l,ori,d,:) = nanstd(resp(thisinds,:),1,1)./sqrt(length(thisinds));
            ph_result.allresp{l,ori,d} = resp(thisinds,:);
            ph_result.condlfpresp(l,ori,d,:) = nanmean(lfpresp(thisinds,:),1);
            ph_result.alllfpresp{l,ori,d} = lfpresp(thisinds,:);
            ph_result.condlfpspect(l,ori,d,:) = nanmean(lfpspect(thisinds,:),1);
            ph_result.condlfpspecterr(l,ori,d,:) = nanstd(lfpspect(thisinds,:),1)./sqrt(length(thisinds));            
            ph_result.condfr(l,ori,d) = nanmean(frs(thisinds));%-mean(bl);
            ph_result.conderr(l,ori,d) =nanstd(frs(thisinds))./sqrt(length(thisinds));
            if l == 1
                anovavals = [anovavals; frs(thisinds)];
                anovaori = [anovaori, ones(1,length(thisinds)).*ori];
                anovad = [anovad, ones(1,length(thisinds)).*d];
            end
            
            ph_result.condz(l,ori,d) = {(sc(thisinds)-mean(sc(thisinds)))/std(sc(thisinds))}; %ecker 2010
            ph_result.condsc(l,ori,d) = {sc(thisinds)};
            ph_result.ff(l,ori,d) = var(sc(thisinds))/mean(sc(thisinds));
                            
            binn = floor(size(ph_result.condresp(l,ori,d,:),4)/binwidth);
            ph_result.bta = linspace(1,binn*binwidth,binn);
            if ~isnan(ph_result.condresp(l,ori,d,:))
                [ph_result.bincondresp(l,ori,d,:)] = binit(ph_result.condresp(l,ori,d,:),binwidth);
            else
                ph_result.bincondresp(l,ori,d,:) = binit(ph_result.condresp(l,ori,d,:),binwidth);
            end
            ph_result.binconderr(l,ori,d,:) = binit(ph_result.condresperr(l,ori,d,:),binwidth);            
            
            ph_result.condsparseness(l,ori,d) = sparseness(squeeze(ph_result.bincondresp(l,ori,d,:)));                    
            ph_result.condsparsenesswin(l,ori,d) = sparseness(squeeze(ph_result.bincondresp(l,ori,d,find(ph_result.bta-prestim>respwin(1),1):find(ph_result.bta-prestim>respwin(end),1)-1)));                        
            
            thisruninds = intersect(thisinds,oktrials);
            if ~isempty(thisruninds)
                ph_result.runcondresp(l,ori,d,:) = mean(resp(thisruninds,:),1);
                ph_result.runcondfr(l,ori,d) = mean(frs(thisruninds));
                ph_result.runconderr(l,ori,d) = std(frs(thisruninds))./sqrt(length(thisruninds));
                ph_result.runn(l,ori,d) = length(thisruninds);
            else
                ph_result.runcondresp(l,ori,d,:) = nan(1,size(resp,2));
                ph_result.runcondfr(l,ori,d) = NaN;
                ph_result.runconderr(l,ori,d) = NaN;
                ph_result.runn(l,ori,d) = 0;
            end
            
            thisstillinds = intersect(thisinds,stilltrials);
            if ~isempty(thisstillinds)
                ph_result.stillcondresp(l,ori,d,:) = mean(resp(thisstillinds,:),1);
                ph_result.stillcondfr(l,ori,d) = mean(frs(thisstillinds));
                ph_result.stillconderr(l,ori,d) = std(frs(thisstillinds))./sqrt(length(thisstillinds));
                ph_result.stilln(l,ori,d) = length(thisstillinds);
            else
                ph_result.stillcondresp(l,ori,d,:) = nan(1,size(resp,2));
                ph_result.stillcondfr(l,ori,d) = NaN;
                ph_result.stillconderr(l,ori,d) = NaN;
                ph_result.stilln(l,ori,d) = 0;
            end
            
            %get trial to trial reliability
            clear tc; clear lfptc; clear ftc;
            k = 1;
            if length(thisinds)>1
                binrespwin = find(ph_result.bta>respwin(1),1):find(ph_result.bta>respwin(end),1)-1;
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
                ph_result.condrely(l,ori,d) = nanmean(tc);
                ph_result.condrelyn(l,ori,d) = length(find(resp(thisinds,:)));
                ph_result.condfiltrely(l,ori,d) = nanmean(ftc);
                ph_result.condlfprely(l,ori,d) = nanmean(lfptc);  
            else
                ph_result.condrely(l,ori,d) = NaN;
                ph_result.condrelyn(l,ori,d) = NaN;
                ph_result.condfiltrely(l,ori,d) = NaN;
                ph_result.condlfprely(l,ori,d) = NaN; 
            end
        end
    end
end
ph_result.ordf = ordf; ph_result.oris = oris;
ph_result.bincondresp = ph_result.bincondresp.*(1000/binwidth);
ph_result.bta = ph_result.bta-prestim;
ph_result.fax = trialfax;
[p,table,stats] = anovan(anovavals,{anovad anovaori},'model','full','display','off');
ph_result.anova_dp = p(1); ph_result.anova_op = p(2); ph_result.anova_doip = p(3);

% just for lfp
for l = 1:length(unique(result.light))
    for d = 1:length(ordf)
        if ordf(d) == -1
            thisinds = find(result.gratingInfo.Phase_surround == -1 & ...
                result.light == l-1);
        elseif ordf(d) == -2
            thisinds = find(result.gratingInfo.Orientation == -1 & ...
                result.light == l-1);
        else
            thisinds = find(result.gratingInfo.Orientation ~= -1 &...
                result.gratingInfo.Phase_surround == ordf(d) & ...
                result.light == l-1);
        end        
        ph_result.omlfpspect(l,d,:) = nanmean(lfpspect(thisinds,:),1); % orientation mean
        ph_result.omlfpspecterr(l,d,:) = nanstd(lfpspect(thisinds,:),1)./sqrt(length(thisinds));
    end
end

contindsnl = find(result.gratingInfo.Orientation == -1 & result.gratingInfo.Phase_surround == -1 & light == 0);
ph_result.controlresp(1,:) = mean(resp(contindsnl,:),1);
ph_result.allcontrolresp = resp(contindsnl,:);
ph_result.controlfr(1) = mean(frs(contindsnl));
ph_result.controlerr(1) = std(frs(contindsnl))./sqrt(length(contindsnl));
ph_result.controllfpspect(1,:) = nanmean(lfpspect(contindsnl,:));

contindsl = find(result.gratingInfo.Orientation == -1 & result.gratingInfo.Phase_surround == -1 & light == 1);
ph_result.controlresp(2,:) = mean(resp(contindsl,:),1);
ph_result.controlfr(2) = mean(frs(contindsl));
ph_result.controlerr(2) = std(frs(contindsl))./sqrt(length(contindsl));
ph_result.controllfpspect(2,:) = nanmean(lfpspect(contindsl,:));

[ph_result.binnedctr(1,:),bta] = binit(ph_result.controlresp(1,:),binwidth);
[ph_result.binnedctr(2,:),bta] = binit(ph_result.controlresp(2,:),binwidth);

ph_result.sizecenter = unique(result.gratingInfo.sizecenter);
ph_result.sizesurround = unique(result.gratingInfo.sizesurround);

% test for running modulation
l1r1 = frs(intersect(find(light),oktrials));
l0r1 = frs(intersect(find(~light),oktrials));
l1r0 = frs(intersect(find(light),stilltrials));
l0r0 = frs(intersect(find(~light),stilltrials));
anovavec = [l0r0;l0r1;l1r0;l1r1];
g1 = [zeros(length(l0r0),1);zeros(length(l0r1),1);ones(length(l1r0),1);ones(length(l1r1),1)]; %light
g2 = [zeros(length(l0r0),1);ones(length(l0r1),1);zeros(length(l1r0),1);ones(length(l1r1),1)]; %running
[p,table,stats] = anovan(anovavec,{g1 g2},'model','full','display','off');
ph_result.anova_lp = p(1); ph_result.anova_rp = p(2); ph_result.anova_rlip = p(3);

end

