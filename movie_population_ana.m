function movie_population_ana

% SOM later population
animalids = {'150401','150527','150529','150602','150603','150623','150625','150825','150831','150902','150909', '150915', '150916', '150916', '151022', '151023', '151109', '151110', '151209', '151221','160718','160721','160725','160726','160726','160728','160728','160729','160801','160802','160804_2'};
blocks    = [7,        4,       7,       6,       6,       9,       5,       11,      7,       12,      8,        12,       4,        4,        10,       11,       12,       15,       13,       5,       3,       3,       3,       3,       3,       3,       4,       2,       3,       3,       3];
rfblocks  = [1,        6,       1,       1,       2,       11,      4,       14,      3,       13,      10,       13,       1,        1,        2,        10,       6,        5,        7,        4,       1,       1,       1,       1,       2,       1,       2,       1,       2,       2,       2];
animal    = [1,        2,       3,       4,       5,       6,       7,       8,       9,       10,      11,       12,       13,       13,       14,       15,       16,       17,       18,       19,      20,      21,      22,      23,      24,      24,      25,      25,      26,      27,      28];
electrodes =[[1,32];  [1,32];  [1,32];  [1,32];  [1,16];  [1,16];  [1,16];  [17,32]; [1,16];  [1,16];  [1,16];   [17,32];  [1,16];   [17,32];  [1,16];   [17,32];  [17,32];  [1,16];   [17,32];  [1,16];  [17,32]; [1,16];  [1,16];  [1,16];  [17,32]; [1,16];  [17,32]; [1,16];  [1,16];  [1,16];  [1,16]];
penangle =  [25,       25,      25,      25,      10,      10,      10,      25,      25,      25,      25,       25,       25,       25,       25,       25,       25,       25,       25,       25,      25,      25,      25,      25,      25,      25,      25,      25,      25,      25,      25];
% age       [P6,       P2?(P0), P2?(P0), P2?(P0), P1       P1
printpath = 'C:\Users\Julia\work\data\populations\SOM_Halo_later\movie\units\';
lfpprintpath = 'C:\Users\Julia\work\data\populations\SOM_Halo_later\movie\lfp\';
rfprintpath = 'C:\Users\Julia\work\data\populations\SOM_Halo_later\movie\rfs\';
popfile = 'C:\Users\Julia\work\data\populations\SOM_Halo_later\movie\movie_population.mat';
% % 
% % PV Halo population
% animalids = {'150629','150630','150731','150731','150804','150804','150818','150819','150820','150823','150824','151211','160114','160115','160204','160205','160217','160328','160901','160902','160902_2','160906','160907','160922','160923'};
% blocks    = [9,        10,      7,       8,       5,       6,       6,       11,      9,       8,       5,       9,       13,      13,      10,      10,      8,       5,       4,       3,       4,         5,       4,       4,       2];
% rfblocks  = [10,       9,       10,      10,      2,       3,       2,       3,       8,       3,       2,       7,       9,       12,      6,       9,       6,       1,       2,       2,       3,         3,       3,       3,       1];
% animal    = [1,        2,       3,       3,       4,       4,       5,       6,       7,       8,       9,       10,      11,      12,      13,      14,      15,      16,      17,      18,      19,        20,      21,      22,      23];
% electrodes =[[1,16];  [1,16];  [1,16];  [17,32]; [1,16];  [17,32]; [1,16];  [1,16];  [1,16];  [1,16];  [17,32]; [1,16];  [1,16];  [1,16];  [1,16];  [17,32]; [1,16];  [1,16];  [1,16];  [1,16];  [1,16];    [1,16];  [1,32];  [1,16];  [1,16]];
% penangle =  [10,       10,      25,      25,      25,      25,      25,      25,      25,      25,      25,      25,      25,      25,      25,      25,      25,      25,      25,      25,      25,        25,      25,      25,      25];
% % age       [P6,       P2?(P0), P2?(P0), P2?(P0), P1       P1
% printpath = 'C:\Users\Julia\work\data\populations\PV_Halo\movie\units\';
% lfpprintpath = 'C:\Users\Julia\work\data\populations\PV_Halo\movie\lfp\';
% rfprintpath = 'C:\Users\Julia\work\data\populations\PV_Halo\movie\rfs\';
% popfile = 'C:\Users\Julia\work\data\populations\PV_Halo\movie\movie_population.mat';

% % controlpop
% animalids = {'151222', '160113', '160711'};
% blocks    = [12,        9,        3];
% rfblocks  = [10,        8,        2];
% animal    = [1,         2,        3];
% electrodes =[[1,16];   [1,16];   [1,16]];
% penangle =  [25,        25,       25];
% % age       [P6,       P2?(P0), P2?(P0), P2?(P0), P1       P1
% printpath = 'C:\Users\Julia\work\data\populations\control\movie\units\';
% lfpprintpath = 'C:\Users\Julia\work\data\populations\control\movie\lfp\';
% rfprintpath = 'C:\Users\Julia\work\data\populations\control\movie\rfs\';
% popfile = 'C:\Users\Julia\work\data\populations\control\movie\movie_population.mat';



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

tic
lcol = 'r'; %lasercolor

recalculate = 0;
printyn = 1;
sfc = 1;

% chronux parameters
params.tapers = [2,5]; params.Fs = 1000; params.err = [2, 0.05]; params.trialave = 1;

sr = 1000;

if ~exist(popfile) || recalculate

    cll = 1;
    for blck = 1:length(blocks)

        supath = ['C:\Users\Julia\work\data\' animalids{blck} '\singleunits\'];
        basename = [animalids{blck} '_block' int2str(blocks(blck)) '_tet'];
        snbasename = [animalids{blck} '_block' int2str(rfblocks(blck)) '_tet'];

        files = dir([supath, basename, '*.mat']);
        rffiles = dir([supath, snbasename, '*.mat']);
        
        prestim = 300;
        poststim = 300;
        respwin = 1001:2000; % after stimulus onset
        respwin = respwin+prestim;
        freqbinwidth = 5;
        
        clear filtmat; clear powmat; clear phasmat;

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
            tetnos(cll) = tetno;
%             if tetno>8
%                 v1(cll) = logical(0); v2(cll) = logical(1); cllstr = 'V2';
%             else
%                 v1(cll) = logical(1); v2(cll) = logical(0); cllstr = 'V1';
%             end
            
            msStimes = round(result.spikes);
            if ~isempty(msStimes) & msStimes(1) == 0, msStimes(1) = 1; end
            
            chan = zeros(1,length(result.lfp));
            chan(msStimes) = 1;
            
            wvchan = find(var(result.waveforms) == max(var(result.waveforms)));
            cm = [3,4,1,2]; % confusion matrix? take lfp from two electrodes away to not get too many spike related phase resets
            lfp = result.lfp(:,cm(wvchan))';
            nfft = 2^nextpow2(length(lfp));
            fax = sr/2*linspace(0,1,nfft/2+1);
            y = fft(lfp,nfft);
            lfpspectrum = abs(y(1:nfft/2+1));
            
            beta = [15,40];
            gamma = [50,70];
            large = find(result.aperture == 0 & result.light == 0);
            small = find(result.aperture == max(unique(result.aperture)) & result.light == 0);
            for i = 1:length(large)
                [pl(i,:),f] = pmtm(lfp(result.msstamps(large(i)):result.msstamps(large(i))+1000),3,[],1000);
                [ps(i,:),f] = pmtm(lfp(result.msstamps(small(i)):result.msstamps(small(i))+1000),3,[],1000);
            end            
            b1 = find(f>beta(1),1); b2 = find(f>beta(2),1);
            g1 = find(f>gamma(1),1); g2 = find(f>gamma(2),1);
            bsig = nanmean(pl(:,b1:b2));
            gsig = nanmean(ps(:,g1:g2));
            if isempty(find(diff(bsig)>0)) % there is no clear beta peak
                bpi = round((b1+b2)/2);
            else
                peaks = find(diff(bsig)>0)+1;
                pvs = bsig(peaks);
                bpi = peaks(pvs == max(pvs));
                bpi = bpi+b1-1;
            end
            if isempty(find(diff(gsig)>0)) % there is no clear beta peak
                gpi = round((g1+g2)/2);
            else
                peaks = find(diff(gsig)>0)+1;
                pvs = gsig(peaks);
                gpi = peaks(pvs == max(pvs));
                gpi = gpi+g1-1;
            end
            
            gamma1 = eegfilt(lfp,sr,f(bpi)-2.5,f(bpi)+2.5);
            gamma2 = eegfilt(lfp,sr,f(gpi)-2.5,f(gpi)+2.5);
            h1 = hilbert(gamma1); gpow1 = abs(h1); gphas1 = angle(h1);
            h2 = hilbert(gamma2); gpow2 = abs(h2); gphas2 = angle(h2);
            lgi(cll) = bpi; hgi(cll) = gpi;
            
            if sfc
                for i = 1:100/freqbinwidth
                    filtmat(i,:) = eegfilt(lfp,sr,(i-1)*freqbinwidth+1,i*freqbinwidth);
                    h = hilbert(filtmat(i,:));
                    powmat(i,:) = abs(h); phasmat(i,:) = angle(h);
                end
            end           
            
            trialdur = result.movieLengthSecs*1000;
            msstamps = result.msstamps;
            
            if length(msstamps)~=length(result.light)
%                         disp('');
                %         msstamps(38) = [];
                %         result.msstamps = msstamps;
                %         save([supath, files(fi).name],'result');
                pause;
            end
            if ~isfield(result,'position')
                result.position = [10,0];
                save([supath, files(fi).name],'result');
            end
            
            %     trialnfft = 2^nextpow2(800);
            %     trialfax = sr/2*linspace(0,1,trialnfft/2+1);
            for i = 1:length(msstamps)
                resp(i,:) = chan(msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim);
                lfpresp(i,:) = lfp(msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim);  
                hh = find(resp(i,respwin(201:end)))'; ptresp(i).times = hh./1000;                
                %         y = fft(lfpresp(i,1001:1800),trialnfft);
                %         lfpspect(i,:) = abs(y(1:trialnfft/2+1));
                [lfpspect(i,:),trialfax] = pmtm(lfpresp(i,1501:2300),3,[],sr);
                gamma1resp(i,:) = gpow1(msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim);
                gamma2resp(i,:) = gpow2(msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim);
                gphase1resp(i,:) = gphas1(msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim);
                gphase2resp(i,:) = gphas2(msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim);
                
                if sfc
                    for j = 1:size(phasmat,1)
                        allphaseresp(j,i,:) = phasmat(j, msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim);
                        allpowresp(j,i,:) = powmat(j, msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim);
                    end
                end
                
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
            vrfr = sum(resp(:,prestim+40:2*prestim+40-1),2);
            vismod(cll) = ttest2(blfr,vrfr);
            visdriven(cll) = mean(vrfr)>=mean(blfr)+2; % average firing rate is increase at least 2Hz above baseline
            
            cllname{cll} = files(fi).name;
            animalno(cll) = animal(blck);
            pangle(cll) = penangle(blck);
            
            spike = result.waveforms(:,wvchan);
            interpspike = spline(1:32,spike,1:.1:32);
            [adiff(cll),swidth(cll),ptr(cll),eslope(cll)] = spikequant(interpspike);
            
            waveform(cll,:) = spike;
            clustqual(cll) = result.clusterquality;
            
            l0 = result.light == 0; l1 = result.light == 1;
            
            % phases
            tmp1 = zeros(size(gphase1resp));
            tmp1(find(resp)) = gphase1resp(find(resp));
            tmp2 = zeros(size(gphase2resp));
            tmp2(find(resp)) = gphase2resp(find(resp));
            l0tmp1 = tmp1(l0,:); l1tmp1 = tmp1(l1,:);
            l0tmp2 = tmp2(l0,:); l1tmp2 = tmp2(l1,:);
            g1phasesl0{cll} = l0tmp1(find(l0tmp1));
            g1phasesl1{cll} = l1tmp1(find(l1tmp1));
            g2phasesl0{cll} = l0tmp2(find(l0tmp2));
            g2phasesl1{cll} = l1tmp2(find(l1tmp2));
            
            g1phaserl0(cll) = circ_r(g1phasesl0{cll});
            g1phaserl1(cll) = circ_r(g1phasesl1{cll});
            g2phaserl0(cll) = circ_r(g2phasesl0{cll});
            g2phaserl1(cll) = circ_r(g2phasesl1{cll});
            g1cmeanl0(cll) = circ_mean(g1phasesl0{cll});
            g1cmeanl1(cll) = circ_mean(g1phasesl1{cll});
            g2cmeanl0(cll) = circ_mean(g2phasesl0{cll});
            g2cmeanl1(cll) = circ_mean(g2phasesl1{cll});
            g1ppcl0(cll) = ppc(g1phasesl0{cll});
            g1ppcl1(cll) = ppc(g1phasesl1{cll});
            g2ppcl0(cll) = ppc(g2phasesl0{cll});
            g2ppcl1(cll) = ppc(g2phasesl1{cll});
            
            if sfc
                tmpi = zeros(size(allphaseresp));
                for i = 1:size(allphaseresp,1)
                    tmp = zeros(size(allphaseresp,2),size(allphaseresp,3));
                    tmp(find(resp)) = allphaseresp(i,find(resp));
                    tmpi(i,:,:) = tmp;
                end
                allphasemat = tmpi;
                for i = 1:size(allphaseresp,1)
                    allphases{i} = allphasemat(i,find(squeeze(allphasemat(i,:,:))));
                    allphaser(cll,i) = circ_r(allphases{i}');
                    allcmean(cll,i) = circ_mean(allphases{i}');
                end
            end
            
            depth(cll) = result.depth;
%             cllresp(cll,:,:) = resp;
%             clllfpresp(cll,:,:) = lfpresp;
            
            %get RF of cll, too
            [on(cll,:,:),off(cll,:,:),gaussfiton(cll),gaussfitoff(cll),rson(cll),rsoff(cll),xax(cll,:),yax(cll,:)] = get_rf([supath, rffiles(fi).name]);
            
            msta = linspace(-prestim,trialdur+poststim,size(resp,2));
            
            baseline = mean(bl);
            baselineerr = std(bl)./(sqrt(size(bl,1)));
            
            binwidth = 33.333;
            clevels = unique(result.contrast);
            widths = unique(result.aperture);
            movies = unique(result.movieno);
            
            for l = 1:2
                for sz = 1:length(widths)
                    thisinds = find(result.aperture == widths(sz) &...
                        result.light == l-1 & result.contrast == 1);
                    condn(cll,l,sz) = length(thisinds);
                    condallresp{cll,l,sz} = resp(thisinds,:);
                    condalllfpresp{cll,l,sz} = lfpresp(thisinds,:);
                    condresp(l,sz,:) = nanmean(resp(thisinds,:),1);
                    condlfpresp(cll,l,sz,:) = nanmean(lfpresp(thisinds,:),1);
                    condresperr(l,sz,:) = nanstd(resp(thisinds,:),1,1)./sqrt(length(thisinds));
                    condfiltresp(cll,l,sz,:) = filter(excit_kernel,1,nanmean(resp(thisinds,:),1));
                    if ~isnan(condresp(l,sz,:))
                        [bincondresp(l,sz,:),bta] = binit(condresp(l,sz,:),binwidth);
                    else
                        bincondresp(l,sz,:) = binit(condresp(l,sz,:),binwidth);
                    end                    
                    bta = bta-prestim;
                    binconderr(l,sz,:) = binit(condresperr(l,sz,:),binwidth); 
                    bincondresp(l,sz,:) = bincondresp(l,sz,:).*(1000/binwidth);
                    binconderr(l,sz,:) = binconderr(l,sz,:).*(1000/binwidth);
                    condfr(cll,l,sz) = nanmean(frs(thisinds));
                    conderr(cll,l,sz) = nanstd(frs(thisinds))./sqrt(length(thisinds));
                    condfrstd(cll,l,sz) = nanstd(frs(thisinds));
                    condff(cll,l,sz) = var(frs(thisinds))/mean(frs(thisinds));
                    
                    [C(cll,l,sz,:),phi(cll,l,sz,:),S12(cll,l,sz,:),S1(cll,l,sz,:),...
                        S2(cll,l,sz,:),chfx,zerosp,confC,...
                        phistd,Cerr(cll,l,sz,:,:)] = coherencycpt(squeeze(lfpresp(thisinds,respwin(201:end)))',...
                        ptresp(thisinds),params);
                                   
                    [S,chf,Serr]=mtspectrumc(squeeze(lfpresp(thisinds,respwin(201:end)))',params);
                    condS(cll,l,sz,:) = S(1:124); condSerr(cll,l,sz,:,:) = Serr(:,1:124);
                    
                    condsparseness(cll,l,sz) = sparseness(squeeze(bincondresp(l,sz,:)));
                    condsparsenesswin(cll,l,sz) = sparseness(squeeze(bincondresp(l,sz,find(bta>respwin(1),1):find(bta>respwin(end),1)-1)));
                    
                    %get trial to trial reliability
                    clear tc; clear lfptc; clear ftc;
                    k = 1;
                    binrespwin = find(bta+prestim>respwin(1),1):find(bta+prestim>respwin(end),1)-1;
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
                    condrely(cll,l,sz) = nanmean(tc);
                    condrelyn(cll,l,sz) = length(find(resp(thisinds,:)));
                    condfiltrely(cll,l,sz) = nanmean(ftc);
                    condlfprely(cll,l,sz) = nanmean(lfptc);
                    
                    
                    % find events and determine jitter
                    clear msi; clear mat; clear jitter; clear jittersd;
                    events = find(bincondresp(l,sz,:)>condfr(cll,l,sz)+3*condfrstd(cll,l,sz));
                    if find(events == 1), events(1) = []; end
                    if find(events == length(bta)), events(end) = []; end
                    for i = 1:length(events)
                        [msm,msi(i)] = min(abs(msta-bta(events(i))));
                        mat = resp(thisinds,msi(i)-round(binwidth):msi(i)+round(binwidth));
                        lats = []; 
                        for j = 1:size(mat,1)
                            lats = [lats,find(mat(j,:))];
                        end
                        jitter(i) = var(lats);
                        jittersd(i) = std(lats);
                    end
                    % if there are events and an event is more than one spike
                    if ~isempty(events) 
                        if ~(length(find(condallresp{cll,l,sz})) == length(events))
                            condjit{cll,l,sz} = jitter;
                            eventtimes{cll,l,sz} = msta(msi);
                            meanjit(cll,l,sz) = mean(jitter);
                            meanjitsd(cll,l,sz) = mean(jittersd);
                            nevents(cll,l,sz) = length(events);
                        else
                            condjit{cll,l,sz} = [];
                            eventtimes{cll,l,sz} = [];
                            meanjit(cll,l,sz) = NaN;
                            meanjitsd(cll,l,sz) = NaN;
                            nevents(cll,l,sz) = 0;
                        end
                    else
                        condjit{cll,l,sz} = [];
                        eventtimes{cll,l,sz} = [];
                        meanjit(cll,l,sz) = NaN;
                        meanjitsd(cll,l,sz) = NaN;
                        nevents(cll,l,sz) = 0;
                    end
                end
            end
            sfx = chf(1:124);
            
            bcresp(cll,:,:,:) = bincondresp;
            bcerr(cll,:,:,:) = binconderr;
            
            lightmod(cll) = ttest2(frs(find(result.light)),frs(find(result.light == 0)));
            aperturevismod(cll) = ttest(frs(find(result.light == 0 & result.aperture ~= 0 )),bl(find(result.light == 0 & result.aperture ~= 0)));
            contextvismod(cll) = ttest(frs(find(result.light == 0 & result.aperture == 0 )),bl(find(result.light == 0 & result.aperture == 0)));
            
%             figure
            clf;
            rasterplot4(resp,result.aperture == widths(2) & result.light == 0,...
                result.aperture == widths(2) & result.light == 1,...
                result.aperture == widths(1) & result.light == 0,...
                result.aperture == widths(1) & result.light == 1,msta);
            title(strcat(cllname(cll), ' depth: ', num2str(depth(cll)), ' swidth: ', num2str(swidth(cll))))
            if printyn
                figSize = [30 21];
                set(gcf,'paperunits','centimeters','papersize',figSize,'paperposition',[0 0 figSize])
                if cll<10, printi = ['0', int2str(cll)]; else printi = int2str(cll); end
                print([printpath ,  printi '__' files(fi).name '.pdf'],'-dpdf')
            end
            
%             figure
            clf;
            subplot(2,2,1)
            semilogy(trialfax,mean(lfpspect(l0&result.aperture == widths(2),:)),'linewidth',2);
            hold on
            semilogy(trialfax,mean(lfpspect(l1&result.aperture == widths(2),:)),'r','linewidth',2);
            semilogy(trialfax,mean(lfpspect(l0&result.aperture == widths(1),:)),'c','linewidth',2);
            semilogy(trialfax,mean(lfpspect(l1&result.aperture == widths(1),:)),'m','linewidth',2);
            %     semilogy(trialfax,mean(lfpspect(l0&result.aperture == widths(2),:))-(std(lfpspect(l0&result.aperture == widths(2),:))./sqrt(size(l0&result.aperture == widths(2),2))));
            %     semilogy(trialfax,mean(lfpspect(l0&result.aperture == widths(2),:))+(std(lfpspect(l0&result.aperture == widths(2),:))./sqrt(size(l0&result.aperture == widths(2),2))));
            %     semilogy(trialfax,mean(lfpspect(l1&result.aperture == widths(2),:))-(std(lfpspect(l1&result.aperture == widths(2),:))./sqrt(size(l1&result.aperture == widths(2),2))),'r');
            %     semilogy(trialfax,mean(lfpspect(l1&result.aperture == widths(2),:))+(std(lfpspect(l1&result.aperture == widths(2),:))./sqrt(size(l1&result.aperture == widths(2),2))),'r');
            %     semilogy(trialfax,mean(lfpspect(l0&result.aperture == widths(1),:))-(std(lfpspect(l0&result.aperture == widths(1),:))./sqrt(size(l0&result.aperture == widths(1),2))),'c');
            %     semilogy(trialfax,mean(lfpspect(l0&result.aperture == widths(1),:))+(std(lfpspect(l0&result.aperture == widths(1),:))./sqrt(size(l0&result.aperture == widths(1),2))),'c');
            %     semilogy(trialfax,mean(lfpspect(l1&result.aperture == widths(1),:))-(std(lfpspect(l1&result.aperture == widths(1),:))./sqrt(size(l1&result.aperture == widths(1),2))),'m');
            %     semilogy(trialfax,mean(lfpspect(l1&result.aperture == widths(1),:))+(std(lfpspect(l1&result.aperture == widths(1),:))./sqrt(size(l1&result.aperture == widths(1),2))),'m');
            axis([0,120,...
                min([min(squeeze(mean(lfpspect(:,1:125)))),min(squeeze(mean(lfpspect(:,1:125))))]),...
                max([max(squeeze(mean(lfpspect(:,1:125)))),max(squeeze(mean(lfpspect(:,1:125))))])])
            xlabel('frequency [Hz]')
            ylabel('spectral power')
            legend([{['CRF L0']},{['CRF L1']},{['FS L0']},{['FS L1']}],'location','ne')
            title(['LFP spectra depth: ' num2str(depth(cll))])
            
            subplot(2,2,2)
            plot(mean(gamma1resp(l0,:)));
            hold on
            plot(mean(gamma1resp(l1,:)),'r')
            title(['gamma1 power in time depth: ' int2str(depth(cll))])
            
            subplot(2,2,3)
            [to0,ro0] = rose(g1phasesl0{cll});
            [to1,ro1] = rose(g1phasesl1{cll});
            polar(to0,ro0,'b')
            hold on
            polar(to1,ro1,'r')
            title(['gamma phase locking of unit ' int2str(cll) ' spikewidth: ' int2str(swidth(cll))])
            
            subplot(2,2,4)
            plot(mean(lfpresp(l0,:)));
            hold on
            plot(mean(lfpresp(l1,:)),'r')
            
            if printyn
                figSize = [30 21];
                set(gcf,'paperunits','centimeters','papersize',figSize,'paperposition',[0 0 figSize])
                if cll<10, printi = ['0', int2str(cll)]; else printi = int2str(cll); end
                print([lfpprintpath ,  printi '__' files(fi).name '.pdf'],'-dpdf')
            end
            
            
%             figure
            clf
            subplot(2,1,1)
            plot(bta, squeeze(bincondresp(1,2,:)),'b') % movie 1 fc aperture ligth off
            hold on
            plot(bta, squeeze(bincondresp(2,2,:)),'r') % movie 1 fc aperture ligth on
            plot(bta, squeeze(bincondresp(1,1,:)),'c') % movie 1 fc full screen ligth off
            plot(bta, squeeze(bincondresp(2,1,:)),'m') % movie 1 fc full screen ligth off
            title([cllname{cll}, ' depth: ' num2str(depth(cll)) ' width: ' num2str(swidth(cll))])
            legend({'CRF light OFF','CRF light ON','contextual light OFF','contextual light ON'})
            ax = axis;
            axis([-300,3300,ax(3),ax(4)])
            line([0,0],[ax(3),ax(4)],'color','k')
            line([3000,3000],[ax(3),ax(4)],'color','k')
            line([1000,1000],[ax(3),ax(4)],'color','r')
            line([2000,2000],[ax(3),ax(4)],'color','r')
            
            ond(cll) = sqrt((gaussfiton(cll).shiftxcenterDeg-result.position(1)).^2+(gaussfiton(cll).shiftycenterDeg-result.position(2)).^2);
            offd(cll) = sqrt((gaussfitoff(cll).shiftxcenterDeg-result.position(1)).^2+(gaussfitoff(cll).shiftycenterDeg-result.position(2)).^2);
            meanspreadon(cll) = (gaussfiton(cll).xspreadDeg+gaussfiton(cll).yspreadDeg)/2;
            meanspreadoff(cll) = (gaussfitoff(cll).xspreadDeg+gaussfitoff(cll).yspreadDeg)/2;
            onoverlap(cll) = (meanspreadon(cll)+max(unique(result.aperture))/2-ond(cll))./(meanspreadon(cll)+max(unique(result.aperture))/2+ond(cll));
            offoverlap(cll) = (meanspreadoff(cll)+max(unique(result.aperture))/2-offd(cll))./(meanspreadoff(cll)+max(unique(result.aperture))/2+offd(cll));           
                      
            subplot(2,2,3)
            imagesc(xax(cll,:),yax(cll,:),squeeze(on(cll,:,:)))
            hold on
            plot_orrf_absdeg(gaussfiton(cll),1,'w',2)
            axis square
            axis xy
            plot_circle(result.position(1),result.position(2),max(unique(result.aperture))/2,'k',2)
            title(['on overlap: ' num2str(onoverlap(cll))])
            
            subplot(2,2,4)
            imagesc(xax(cll,:),yax(cll,:),squeeze(off(cll,:,:)))
            hold on
            plot_orrf_absdeg(gaussfitoff(cll),1,'w',2)
            axis square
            axis xy
            plot_circle(result.position(1),result.position(2),max(unique(result.aperture))/2,'k',2)
            title(['off overlap: ' num2str(offoverlap(cll))])
            
            aperturepos(cll,:) = result.position;
            
            if printyn
                figSize = [30 21];
                set(gcf,'paperunits','centimeters','papersize',figSize,'paperposition',[0 0 figSize])
                if cll<10, printi = ['0', int2str(cll)]; else printi = int2str(cll); end
                print([rfprintpath ,  printi '__' files(fi).name '.pdf'],'-dpdf')
            end
            
            
            %     subplot(2,1,2)
            %     plot(bta, squeeze(bincondresp(1,2,mci,2,:)),'b') % movie 2 fc aperture ligth off
            %     hold on
            %     plot(bta, squeeze(bincondresp(2,2,mci,2,:)),'r') % movie 2 fc aperture ligth on
            %     plot(bta, squeeze(bincondresp(1,2,mci,1,:)),'c') % movie 2 fc full screen ligth off
            %     plot(bta, squeeze(bincondresp(2,2,mci,1,:)),'m') % movie 2 fc full screen ligth off
            %     title(cllname{cll})
            %     ax = axis;
            %     axis([-300,3300,ax(3),ax(4)])
            %     line([0,0],[ax(3),ax(4)],'color','k')
            %     line([3000,3000],[ax(3),ax(4)],'color','k')
            %     line([1000,1000],[ax(3),ax(4)],'color','r')
            %     line([2000,2000],[ax(3),ax(4)],'color','r')
            
            
            disp(['done with ' files(fi).name]);
            cll = cll+1;
            
        end
    end
    save(popfile, '-v7.3'); 
else
    load(popfile);   
end

toc
%spike classification
kmeansind = kmeans([eslope',ptr',swidth',adiff'],3);
% kmeansind = kmeans([eslope',adiff'],2);
secpersamp = 1/30000;
interpf = secpersamp/10;
swidthms = swidth*interpf*1000;

% m1 = mean(swidth(kmeansind==1)); m2 = mean(swidth(kmeansind==2)); m3 = mean(swidth(kmeansind==3));
% fm = min([m1,m2,m3]); fmax = max([m1,m2,m3]);
% if fm == m1
%     pfs = find(kmeansind==1); pfsv = kmeansind==1;
% elseif fm == m2
%     pfs = find(kmeansind==2); pfsv = kmeansind == 2;
% else
%     pfs = find(kmeansind==3); pfsv = kmeansind == 3;
% end
% if fmax == m1
%     prs = find(kmeansind==1); prsv = kmeansind==1;
% elseif fmax == m2
%     prs = find(kmeansind==2); prsv = kmeansind == 2;
% else
%     prs = find(kmeansind==3); prsv = kmeansind == 3;
% end

m1 = mean(swidth(kmeansind==1)); m2 = mean(swidth(kmeansind==2)); m3 = mean(swidth(kmeansind==3));
fm = min([m1,m2,m3]);
if fm == m1
    pfs = find(kmeansind==1); prs = find(kmeansind==2|kmeansind==3); 
    pfsv = kmeansind==1; prsv = kmeansind==2|kmeansind==3;
elseif fm == m2
    pfs = find(kmeansind==2); prs = find(kmeansind==1|kmeansind==3);
    pfsv = kmeansind==2; prsv = kmeansind==1|kmeansind==3;
else
    pfs = find(kmeansind==3); prs = find(kmeansind==1|kmeansind==2);
    pfsv = kmeansind==3; prsv = kmeansind==1|kmeansind==2;
end

% if mean(swidth(find(kmeansind==1)))<mean(swidth(find(kmeansind==2)))  %1 is FS
%     pfs = find(kmeansind==1); prs = find(kmeansind==2); pfsv = kmeansind==1;
% else
%     pfs = find(kmeansind==2); prs = find(kmeansind==1); pfsv = kmeansind==2;
% end

figure
plot(swidthms(kmeansind==1),adiff(kmeansind==1),'b.')
xlabel('spike width')
ylabel('amplitude diff')
hold on
plot(swidthms(kmeansind==2),adiff(kmeansind==2),'r.')
if ~isempty(find(kmeansind==3))
    plot(swidthms(kmeansind==3),adiff(kmeansind==3),'g.')
    plot(swidthms(pfsv),adiff(pfsv),'ro');
    plot(swidthms(prsv),adiff(prsv),'o');
end

% figure
% plot(eslope(kmeansind==1),adiff(kmeansind==1),'b.')
% xlabel('end slope')
% ylabel('amplitude diff')
% hold on
% plot(eslope(kmeansind==2),adiff(kmeansind==2),'r.')
% if ~isempty(find(kmeansind==3))
%     plot(eslope(kmeansind==3),adiff(kmeansind==3),'g.')
%     plot(eslope(pfsv),adiff(pfsv),'ro');
%     plot(eslope(prsv),adiff(prsv),'o');
% end
%     % axis([5.5,20.5,-.9,.7])

prsv = swidthms>=.38; pfsv = swidthms<=.36;
prs = find(prsv); pfs = find(pfsv);

swamp = max(waveform,[],2)-min(waveform,[],2);
okwv = swamp>42;
vismod(isnan(vismod)) = 0;
% ok = vismod&okwv'&nlfr>1;
ok = okwv';

for i = 1:length(okwv)
    isodist(i) = clustqual(i).IsolationDistance;
    lratio(i) = clustqual(i).L_Ratio.Lratio;
end

% adjust depth according to penetration angle
depth = depth.*cosd(22).*cosd(pangle);

% halo expressing clls
phe = zeros(1,length(depth));

% % putative halo expressing for SOM later pop: 
phe([6,33,35,64, 75,97, 288, 326,328, 333]) = 1;

% % putative halo expressing for PV Halo pop: 
% phe([55,111,185,213,261,276,291,295, 309, 313, 315,316, 326, 329, 339, 342]) = 1;

phe = logical(phe)';

% figure out if RF mappable and overlapping with movie aperture
okrfon = [gaussfiton.rsquared]>.5 & meanspreadon>1.5; %[gaussfiton.amp]>1 & 
okrfoff = [gaussfitoff.rsquared]>.5 & meanspreadoff>1.5; %[gaussfitoff.amp]>1 & 
nonmappable = ~okrfon&~okrfoff;
apok =(okrfon&onoverlap>.33) | (okrfoff & offoverlap>.33);

% figure out which had enough spikes for reliability estimates
relyok = condrelyn>30; % need at least as many spikes as number of trials that were cross correlated

% sort out jitter that is wrong because events are only 1 spike in one trial
for i = 1:length(depth)
    for l = 1:2
        for s = 1:2
            if length(find(condallresp{i,l,s})) == nevents(i,l,s)
                meanjit(i,l,s) = NaN;
                meanjitsd(i,l,s) = NaN;
            end
        end
    end
end

[sd,si] = sort(depth);
[rsd,rsi] = sort(depth(prs));
[fsd,fsi] = sort(depth(pfs));

% fsl4 = find(depth(pfs)>=375 & depth(pfs)<=500);
% rsl4 = find(depth(prs)>=375 & depth(prs)<=500);
% fsl23 = find(depth(pfs)<375);
% rsl23 = find(depth(prs)<375);
% fsl5 = find(depth(pfs)>500 & depth(pfs)<=800);
% rsl5 = find(depth(prs)>500 & depth(prs)<=800);
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

aperturevismod(isnan(aperturevismod)) = 0;

% % SOM Halo later - for now all ~350, 1=525 - animal 4 too deep
% lfpinds = [19, 51, 62, 82, 91, 103, 112, 119, 132, 142, 153, 179, 194, 197, 216, 231, 240, 248, 253, 260, 271, 321, 338, 350, 362, 389];

% find LFPs closest to 350um
lfpinds = [];
totall = 0;
for i = 1:max(animalno)
    a = find(animalno == i);
    if ~isempty(a)
        b = depth(a)-350;
        [c,ci] = min(abs(b));
        if c<50
            lfpinds = [lfpinds,ci+totall];
        end
        totall = totall+length(a);
    end
end
   

for i = 1:length(depth)
    for l = 1:2
        for sz = 1:2            
            peakgammasfc(i,l,sz) = squeeze(C(i,l,sz,lgi(i)));
            peakgammapow(i,l,sz) = squeeze(condS(i,l,sz,lgi(i)));
            fillyspect(i,l,sz,:) = [squeeze(condSerr(i,l,sz,1,:))',fliplr(squeeze(condSerr(i,l,sz,2,:))')];
        end
    end
end
fillx = [sfx,fliplr(sfx)];

% firing rate, sparseness, reliability and precision with aperture size changes
rsp = signrank(condfr(l23rs&apok,1,1),condfr(l23rs&apok,1,2));
fsp = signrank(condfr(l23fs&apok,1,1),condfr(l23fs&apok,1,2));
figure
loglog(condfr(l23rs&apok,1,2),condfr(l23rs&apok,1,1),'ko','markerfacecolor','k','markersize',5)
hold on
loglog(condfr(l23fs&apok,1,2),condfr(l23fs&apok,1,1),'go','markerfacecolor','g','markersize',5)
% axis([0,1,0,1])
axis square
refline(1,0);
xlabel('firing rate small aperture')
ylabel('firing rate full screen')
title(['firing rate large vs small RS: ' num2str(rsp) ' FS: ' num2str(fsp)])

rsp = signrank(condsparseness(l23rs&apok,1,1),condsparseness(l23rs&apok,1,2));
fsp = signrank(condsparseness(l23fs&apok,1,1),condsparseness(l23fs&apok,1,2));
figure
plot(condsparseness(l23rs&apok,1,2),condsparseness(l23rs&apok,1,1),'ko','markerfacecolor','k','markersize',5)
hold on
plot(condsparseness(l23fs&apok,1,2),condsparseness(l23fs&apok,1,1),'go','markerfacecolor','g','markersize',5)
axis([0,1,0,1])
axis square
refline(1,0);
xlabel('sparseness small aperture')
ylabel('sparseness full screen')
title(['sparseness large vs small RS: ' num2str(rsp) ' FS: ' num2str(fsp)])

rsp = signrank(condrely(l23rs&apok&relyok(:,1,1)'&relyok(:,1,2)',1,1),condrely(l23rs&apok&relyok(:,1,1)'&relyok(:,1,2)',1,2));
fsp = signrank(condrely(l23fs&apok&relyok(:,1,1)'&relyok(:,1,2)',1,1),condrely(l23fs&apok&relyok(:,1,1)'&relyok(:,1,2)',1,2));
figure
plot(condrely(l23rs&apok&relyok(:,1,1)'&relyok(:,1,2)',1,2),condrely(l23rs&apok&relyok(:,1,1)'&relyok(:,1,2)',1,1),'ko','markerfacecolor','k','markersize',5)
hold on
plot(condrely(l23fs&apok&relyok(:,1,1)'&relyok(:,1,2)',1,2),condrely(l23fs&apok&relyok(:,1,1)'&relyok(:,1,2)',1,1),'go','markerfacecolor','g','markersize',5)
axis([0,1,0,1])
axis square
refline(1,0);
xlabel('reliability small aperture')
ylabel('reliability full screen')
title(['reliability large vs small RS: ' num2str(rsp) ' FS: ' num2str(fsp)])

rsp = signrank(meanjitsd(l23rs&apok,1,1),meanjitsd(l23rs&apok,1,2));
fsp = signrank(meanjitsd(l23fs&apok,1,1),meanjitsd(l23fs&apok,1,2));
figure
plot(meanjitsd(l23rs&apok,1,2),meanjitsd(l23rs&apok,1,1),'ko','markerfacecolor','k','markersize',5)
hold on
plot(meanjitsd(l23fs&apok,1,2),meanjitsd(l23fs&apok,1,1),'go','markerfacecolor','g','markersize',5)
% axis([0,1,0,1])
axis square
refline(1,0);
xlabel('jitter sd small aperture')
ylabel('jitter sd full screen')
title(['jitter sd large vs small RS: ' num2str(rsp) ' FS: ' num2str(fsp)])

% firing rate, sparseness, reliability and precision with light small size
rsp = signrank(condfr(l23rs&apok,1,2),condfr(l23rs&apok,2,2));
fsp = signrank(condfr(l23fs&apok,1,2),condfr(l23fs&apok,2,2));
figure
loglog(condfr(l23rs&apok,1,2),condfr(l23rs&apok,2,2),'ko','markerfacecolor','k','markersize',5)
hold on
loglog(condfr(l23fs&apok,1,2),condfr(l23fs&apok,2,2),'go','markerfacecolor','g','markersize',5)
% axis([0,1,0,1])
axis square
refline(1,0);
xlabel('firing rate small NO light')
ylabel('firing rate small light')
title(['firing rate L0 vs L1 small RS: ' num2str(rsp) ' FS: ' num2str(fsp)])

rsp = signrank(condsparsenesswin(l23rs&apok,1,2),condsparsenesswin(l23rs&apok,2,2));
fsp = signrank(condsparsenesswin(l23fs&apok,1,2),condsparsenesswin(l23fs&apok,2,2));
figure
plot(condsparsenesswin(l23rs&apok,1,2),condsparsenesswin(l23rs&apok,2,2),'ko','markerfacecolor','k','markersize',5)
hold on
plot(condsparsenesswin(l23fs&apok,1,2),condsparsenesswin(l23fs&apok,2,2),'go','markerfacecolor','g','markersize',5)
% axis([0,1,0,1])
axis square
refline(1,0);
xlabel('sparseness small No light')
ylabel('sparseness small light')
title(['sparseness L0 vs L1 small RS: ' num2str(rsp) ' FS: ' num2str(fsp)])

rsp = signrank(condrely(l23rs&apok&relyok(:,1,2)'&relyok(:,2,2)',1,2),condrely(l23rs&apok&relyok(:,1,2)'&relyok(:,2,2)',2,2));
fsp = signrank(condrely(l23fs&apok&relyok(:,1,2)'&relyok(:,2,2)',1,2),condrely(l23fs&apok&relyok(:,1,2)'&relyok(:,2,2)',2,2));
figure
plot(condrely(l23rs&apok&relyok(:,1,2)'&relyok(:,2,2)',1,2),condrely(l23rs&apok&relyok(:,1,2)'&relyok(:,2,2)',2,2),'ko','markerfacecolor','k','markersize',5)
hold on
plot(condrely(l23fs&apok&relyok(:,1,2)'&relyok(:,2,2)',1,2),condrely(l23fs&apok&relyok(:,1,2)'&relyok(:,2,2)',2,2),'go','markerfacecolor','g','markersize',5)
% axis([0,1,0,1])
axis square
refline(1,0);
xlabel('reliability small No light')
ylabel('reliability small light')
title(['reliability L0 vs. L1 small RS: ' num2str(rsp) ' FS: ' num2str(fsp)])

rsp = signrank(meanjitsd(l23rs&apok,1,2),meanjitsd(l23rs&apok,2,2));
fsp = signrank(meanjitsd(l23fs&apok,1,2),meanjitsd(l23fs&apok,2,2));
figure
plot(meanjitsd(l23rs&apok,1,2),meanjitsd(l23rs&apok,2,2),'ko','markerfacecolor','k','markersize',5)
hold on
plot(meanjitsd(l23fs&apok,1,2),meanjitsd(l23fs&apok,2,2),'go','markerfacecolor','g','markersize',5)
% axis([0,1,0,1])
axis square
refline(1,0);
xlabel('jitter sd small NO light')
ylabel('jitter sd small light')
title(['jitter sd L0 vs L1 RS: ' num2str(rsp) ' FS: ' num2str(fsp)])

% firing rate, sparseness and reliability with light large size
rsp = signrank(condfr(l23rs,1,1),condfr(l23rs,2,1));
fsp = signrank(condfr(l23fs,1,1),condfr(l23fs,2,1));
figure
loglog(condfr(l23rs,1,1),condfr(l23rs,2,1),'ko','markerfacecolor','k','markersize',5)
hold on
loglog(condfr(l23fs,1,1),condfr(l23fs,2,1),'go','markerfacecolor','g','markersize',5)
% axis([0,1,0,1])
axis square
refline(1,0);
xlabel('firing rate large NO light')
ylabel('firing rate large light')
title(['firing rate L0 vs L1 large RS: ' num2str(rsp) ' FS: ' num2str(fsp)])

rsp = signrank(condsparsenesswin(l23rs,1,1),condsparsenesswin(l23rs,2,1));
fsp = signrank(condsparsenesswin(l23fs,1,1),condsparsenesswin(l23fs,2,1));
figure
plot(condsparsenesswin(l23rs,1,1),condsparsenesswin(l23rs,2,1),'ko','markerfacecolor','k','markersize',5)
hold on
plot(condsparsenesswin(l23fs,1,1),condsparsenesswin(l23fs,2,1),'go','markerfacecolor','g','markersize',5)
% axis([0,1,0,1])
axis square
refline(1,0);
xlabel('sparseness large No light')
ylabel('sparseness large light')
title(['sparseness L0 vs L1 large RS: ' num2str(rsp) ' FS: ' num2str(fsp)])

rsp = signrank(condrely(l23rs&relyok(:,1,1)'&relyok(:,2,1)',1,1),condrely(l23rs&relyok(:,1,1)'&relyok(:,2,1)',2,1));
fsp = signrank(condrely(l23fs&relyok(:,1,1)'&relyok(:,2,1)',1,1),condrely(l23fs&relyok(:,1,1)'&relyok(:,2,1)',2,1));
figure
plot(condrely(l23rs&relyok(:,1,1)'&relyok(:,2,1)',1,1),condrely(l23rs&relyok(:,1,1)'&relyok(:,2,1)',2,1),'ko','markerfacecolor','k','markersize',5)
hold on
plot(condrely(l23fs&relyok(:,1,1)'&relyok(:,2,1)',1,1),condrely(l23fs&relyok(:,1,1)'&relyok(:,2,1)',2,1),'go','markerfacecolor','g','markersize',5)
% axis([0,1,0,1])
axis square
refline(1,0);
xlabel('reliability large No light')
ylabel('reliability large light')
title(['reliability L0 vs. L1 large RS: ' num2str(rsp) ' FS: ' num2str(fsp)])

rsp = signrank(meanjitsd(l23rs,1,1),meanjitsd(l23rs,2,1));
fsp = signrank(meanjitsd(l23fs,1,1),meanjitsd(l23fs,2,1));
figure
plot(meanjitsd(l23rs,1,1),meanjitsd(l23rs,2,1),'ko','markerfacecolor','k','markersize',5)
hold on
plot(meanjitsd(l23fs,1,1),meanjitsd(l23fs,2,1),'go','markerfacecolor','g','markersize',5)
% axis([0,1,0,1])
axis square
refline(1,0);
xlabel('jitter sd large NO light')
ylabel('jitter sd large light')
title(['jitter sd L0 vs L1 RS: ' num2str(rsp) ' FS: ' num2str(fsp)])

% LFP reliability with aperture
rsp = signrank(condlfprely(lfpinds,1,1),condlfprely(lfpinds,1,2));
figure
plot(condlfprely(lfpinds,1,2),condlfprely(lfpinds,1,1),'ko','markerfacecolor','k','markersize',5)
hold on
axis([0,1,0,1])
axis square
refline(1,0);
xlabel('LFP reliability small aperture')
ylabel('LFP reliability full screen')
title(['LFP reliability large vs small: ' num2str(rsp,3)])

% LFP reliability with light - aperture
rsp = signrank(condlfprely(lfpinds,1,2),condlfprely(lfpinds,2,2));
figure
plot(condlfprely(lfpinds,1,2),condlfprely(lfpinds,2,2),'ko','markerfacecolor','k','markersize',5)
hold on
axis([0,1,0,1])
axis square
refline(1,0);
xlabel('LFP reliability small L0')
ylabel('LFP reliability small L1')
title(['LFP reliability small L0 vs L1: ' num2str(rsp)])

rsp = signrank(condlfprely(lfpinds,1,1),condlfprely(lfpinds,2,1));
figure
plot(condlfprely(lfpinds,1,1),condlfprely(lfpinds,2,1),'ko','markerfacecolor','k','markersize',5)
hold on
axis square
refline(1,0);
xlabel('reliability large No light')
ylabel('reliability large light')
title(['reliability L0 vs. L1 large: ' num2str(rsp) ])

