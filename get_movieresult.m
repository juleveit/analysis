function mv_result = get_movieresult(filename)

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


if ~isfield(result,'position')
    result.position = [10,0];
    save(filename,'result');
end
mv_result.animalid = result.animalid;
mv_result.position = result.position;

prestim = 300;
poststim = 300;
respwin = 1001:2000; % after stimulus onset
respwin = respwin+prestim;

[p,mv_result.cellname] = fileparts(filename);

i = strfind(filename, 'tet');
mv_result.tetno = strread(filename(i+3));

wvchan = find(var(result.waveforms) == max(var(result.waveforms)));

sr = 1000;
lfp = result.lfp(:,wvchan)';

msStimes = round(result.spikes);
if isempty(msStimes), msStimes(1) = 0; end
if msStimes(1) == 0, msStimes(1) = 1; end

chan = zeros(1,length(result.lfp));
chan(msStimes) = 1;

trialdur = result.movieLengthSecs*1000;
msstamps = result.msstamps;
light = result.light;

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

mv_result.depth = result.depth;

spike = result.waveforms(:,wvchan);
mv_result.interpspike = spline(1:32,spike,1:.1:32);
[mv_result.adiff,mv_result.swidth] = spikequant(mv_result.interpspike);

msta = linspace(-prestim,trialdur+poststim,size(resp,2));

frs = sum(resp(:,respwin),2)./(length(respwin)/1000);
bl = sum(resp(:,1:prestim),2)./(prestim/1000);
sc = sum(resp(:,respwin),2);

%determine if cell is visually modulated
blfr = sum(resp(:,1:prestim),2);
vrfr = sum(resp(:,prestim+40:2*prestim+40),2);
mv_result.vismod = signrank(blfr,vrfr);

%determine if cll is modulated by light
mv_result.lightmod = signrank(frs(find(light)),frs(find(light == 0)));

mv_result.lfr = mean(frs(find(light)));
mv_result.nlfr = mean(frs(find(light == 0)));

msta = linspace(-prestim,trialdur+poststim,size(resp,2));

binwidth = 33.33;

[p,printname] = fileparts(filename);
printname(find(printname=='_')) = ' ';
mv_result.printname = printname;

widths = unique(result.aperture);
for l = 1:2
    for sz = 1:length(widths)
        thisinds = find(result.aperture == widths(sz) &...
            result.light == l-1 & result.contrast == 1);
        mv_result.condn(l,sz) = length(thisinds);
        mv_result.condresp(l,sz,:) = nanmean(resp(thisinds,:),1);
        mv_result.condresperr(l,sz,:) = nanstd(resp(thisinds,:),1,1)./sqrt(length(thisinds));
        mv_result.condallresp{l,sz} = resp(thisinds,:);
        mv_result.condlfpresp(l,sz,:) = nanmean(lfpresp(thisinds,:),1);
        mv_result.condalllfpresp{l,sz} = lfpresp(thisinds,:);
        mv_result.condlfpspect(l,sz,:) = nanmean(lfpspect(thisinds,:));
        
        
        mv_result.condfr(l,sz) = mean(frs(thisinds));%-mean(bl);
        mv_result.conderr(l,sz) =std(frs(thisinds))./sqrt(length(thisinds));
        condfrstd(l,sz) = nanstd(frs(thisinds));
        
        mv_result.condz(l,sz) = {(sc(thisinds)-mean(sc(thisinds)))/std(sc(thisinds))}; %ecker 2010
        mv_result.condsc(l,sz) = {sc(thisinds)};
        mv_result.ff(l,sz) = var(sc(thisinds))/mean(sc(thisinds));
        
        if ~isnan(mv_result.condresp(l,sz,:))
            [mv_result.bincondresp(l,sz,:),mv_result.bta] = binit(mv_result.condresp(l,sz,:),binwidth);
        else
            mv_result.bincondresp(l,sz,:) = binit(mv_result.condresp(l,sz,:),binwidth);
        end
        mv_result.binconderr(l,sz,:) = binit(mv_result.condresperr(l,sz,:),binwidth);
        mv_result.bincondresp(l,sz,:) = mv_result.bincondresp(l,sz,:).*(1000/binwidth);
        mv_result.binconderr(l,sz,:) = mv_result.binconderr(l,sz,:).*(1000/binwidth);
        
        mv_result.condsparseness(l,sz) = sparseness(squeeze(mv_result.bincondresp(l,sz,:)));
        mv_result.condsparsenesswin(l,sz) = sparseness(squeeze(mv_result.bincondresp(l,sz,find(mv_result.bta-prestim>respwin(1),1):find(mv_result.bta-prestim>respwin(end),1)-1)));
        
        thisruninds = intersect(thisinds,oktrials);
        if ~isempty(thisruninds)
            mv_result.runcondresp(l,sz,:) = mean(resp(thisruninds,:),1);
            mv_result.runcondfr(l,sz) = mean(frs(thisruninds));
            mv_result.runconderr(l,sz) = std(frs(thisruninds))./sqrt(length(thisruninds));
            mv_result.runn(l,sz) = length(thisruninds);
        else
            mv_result.runcondresp(l,sz,:) = nan(1,size(resp,2));
            mv_result.runcondfr(l,sz) = NaN;
            mv_result.runconderr(l,sz) = NaN;
            mv_result.runn(l,sz) = 0;
        end
        
        thisstillinds = intersect(thisinds,stilltrials);
        if ~isempty(thisstillinds)
            mv_result.stillcondresp(l,sz,:) = mean(resp(thisstillinds,:),1);
            mv_result.stillcondfr(l,sz) = mean(frs(thisstillinds));
            mv_result.stillconderr(l,sz) = std(frs(thisstillinds))./sqrt(length(thisstillinds));
            mv_result.stilln(l,sz) = length(thisstillinds);
        else
            mv_result.stillcondresp(l,sz,:) = nan(1,size(resp,2));
            mv_result.stillcondfr(l,sz) = NaN;
            mv_result.stillconderr(l,sz) = NaN;
            mv_result.stilln(l,sz) = 0;
        end
        
        %get trial to trial reliability
        clear tc; clear lfptc; clear ftc;
        k = 1;
        binrespwin = find(mv_result.bta>respwin(1),1):find(mv_result.bta>respwin(end),1)-1;
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
        mv_result.condrely(l,sz) = nanmean(tc);
        mv_result.condrelyn(l,sz) = length(find(resp(thisinds,:)));
        mv_result.condfiltrely(l,sz) = nanmean(ftc);
        mv_result.condlfprely(l,sz) = nanmean(lfptc);
    
        % find events and determine jitter
        clear msi; clear mat; clear jitter; clear jittersd;
        rbta =mv_result.bta-prestim;
        events = find(mv_result.bincondresp(l,sz,:)>mv_result.condfr(l,sz)+3*condfrstd(l,sz));
        if find(events == 1), events(1) = []; end
        if find(events == length(mv_result.bta)), events(end) = []; end
        for i = 1:length(events)
            [msm,msi(i)] = min(abs(msta-rbta(events(i))));
            mat = resp(thisinds,msi(i)-round(binwidth):msi(i)+round(binwidth));
            lats = [];
            for j = 1:size(mat,1)
                lats = [lats,find(mat(j,:))];
            end
            jitter(i) = var(lats);
            jittersd(i) = std(lats);
        end
        if ~isempty(events)
            if ~(length(find(mv_result.condallresp{l,sz})) == length(events))
                mv_result.condjit{l,sz} = jitter;
                mv_result.eventtimes{l,sz} = msta(msi);
                mv_result.meanjit(l,sz) = mean(jitter);
                mv_result.meanjitsd(l,sz) = mean(jittersd);
                mv_result.nevents(l,sz) = length(events);
            else
                mv_result.condjit{l,sz} = [];
                mv_result.eventtimes{l,sz} = [];
                mv_result.meanjit(l,sz) = NaN;
                mv_result.meanjitsd(l,sz) = NaN;
                mv_result.nevents(l,sz) = 0;
            end
        else
            mv_result.condjit{l,sz} = [];
            mv_result.eventtimes{l,sz} = [];
            mv_result.meanjit(l,sz) = NaN;
            mv_result.meanjitsd(l,sz) = NaN;
            mv_result.nevents(l,sz) = 0;
        end
    end
end
mv_result.widths = widths;
mv_result.bta = mv_result.bta-prestim;


% test for running modulation
l1r1 = frs(intersect(find(light),oktrials));
l0r1 = frs(intersect(find(~light),oktrials));
l1r0 = frs(intersect(find(light),stilltrials));
l0r0 = frs(intersect(find(~light),stilltrials));
anovavec = [l0r0;l0r1;l1r0;l1r1];
g1 = [zeros(length(l0r0),1);zeros(length(l0r1),1);ones(length(l1r0),1);ones(length(l1r1),1)]; %light
g2 = [zeros(length(l0r0),1);ones(length(l0r1),1);zeros(length(l1r0),1);ones(length(l1r1),1)]; %running
[p,table,stats] = anovan(anovavec,{g1 g2},'model','full','display','off');
mv_result.anova_lp = p(1); mv_result.anova_rp = p(2); mv_result.anova_rlip = p(3);



