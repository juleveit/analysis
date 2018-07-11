function freqpulsesana_SU

animalid = '180626';
block = 11;

basepath = ['C:\Users\Julia\work\data\' animalid '\'];
% supath = [basepath 'multiunits\'];
supath = [basepath 'singleunits\'];
basename = [animalid '_block' int2str(block) '_tet'];

files = dir([supath, basename, '*.mat']);

prestim = 300;
% poststim = 300;
poststim = 300;
respwin = 501:3500; % after stimulus onset
% respwin = 501:4500; % after stimulus onset
respwin = respwin+prestim;

cell = 1;
for fi = [10,31] %1:length(files)
    
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
    
    trialdur = result.stimduration*1000;
    msstamps = result.msstamps;
    if length(msstamps)~=length(result.light)
%         msstamps(385) = []; % for 140703 block 5
%         msstamps(62) = []; % for 140703 block 5
%         result.msstamps = msstamps;
%         save([supath, files(fi).name],'result');
        pause;
    end
    
    
    L = length(respwin); nfft = 2^nextpow2(L);
    t = (0:L-1)*(1/1000); fftx = 1000/2*linspace(0,1,nfft/2+1);
    for i = 1:length(msstamps)
        resp(i,:) = chan(msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim); 
        lfpresp(i,:) = result.lfp(msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim);
        speed(i,:) = result.runspeed(msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim); 
        [lfpspect(i,:),fax] = pmtm(lfpresp(i,respwin),3,[],1000);
        hlp = fft(lfpresp(i,respwin),nfft)/L;
        fftspect(i,:) = 2*abs(hlp(1:nfft/2+1));
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
    freqs = unique(result.frequencies);
    for f = 1:length(freqs)
        thisinds = find(result.frequencies == freqs(f));
        condn(f) = length(thisinds);
        condresp(f,:) = nanmean(resp(thisinds,:),1);
        condresperr(f,:) = nanstd(resp(thisinds,:),1,1)./sqrt(length(thisinds));
        condlfpspect(f,:) = nanmean(lfpspect(thisinds,:),1);
        condfftspect(f,:) = nanmean(fftspect(thisinds,:),1);
        condlfpresp(f,:) = nanmean(lfpresp(thisinds,:),1);
        if ~isnan(condresp(f,:))
            [bincondresp(f,:),bta] = binit(condresp(f,:),binwidth);
        else
            bincondresp(f,:) = binit(condresp(f,:),binwidth);
        end
        binconderr(f,:) = binit(condresperr(f,:),binwidth);
        cfr(f) = nanmean(frs(thisinds));
        cerr(f) =nanstd(frs(thisinds))./sqrt(length(thisinds));
        relspect(f,:) = condlfpspect(f,:)./condlfpspect(1,:);
        relfftspect(f,:) = condfftspect(f,:)./condfftspect(1,:);
        
        if freqs(f) == 0
            peakpow(f) = NaN;
            meanpow(f) = NaN;
            powerratio(f) = NaN;
            otherratio(f) = NaN;
            blratio(f) = NaN;
            fftpowerratio(f) = NaN;
            fftotherratio(f) = NaN;
            fftblratio(f) = NaN;
            relpowerratio(f) = NaN;
            relpeakpow(f) = NaN;
            relfftpowerratio(f) = NaN;
            relfftpeakpow(f) = NaN;
        else
            diffs = fax-freqs(f);
            f_ind = find(abs(diffs) == min(abs(diffs)));
            peakpow(f) = mean(condlfpspect(f,f_ind-2:f_ind+2),2);
            blpeakpow(f) = mean(condlfpspect(1,f_ind-2:f_ind+2),2);
            surroundpow(f) = mean([mean(condlfpspect(f,f_ind-7:f_ind-5),2),mean(condlfpspect(f,f_ind+5:f_ind+7),2)]);
            meanpow(f) = mean(condlfpspect(f,1:2000),2);
            powerratio(f) = peakpow(f)./meanpow(f);
            otherratio(f) = peakpow(f)./surroundpow(f);
            blratio(f) = peakpow(f)./blpeakpow(f);
            
            fftpeakpow(f) = condfftspect(f,f_ind);
            fftblpeakpow(f) = condfftspect(1,f_ind);
            fftsurroundpow(f) = mean([mean(condfftspect(f,f_ind-7:f_ind-5),2),mean(condfftspect(f,f_ind+5:f_ind+7),2)]);
            fftmeanpow(f) = mean(condfftspect(f,1:2000),2);
            fftpowerratio(f) = fftpeakpow(f)./fftmeanpow(f);
            fftotherratio(f) = fftpeakpow(f)./fftsurroundpow(f);
            fftblratio(f) = fftpeakpow(f)./fftblpeakpow(f);
            
            relpeakpow(f) = relspect(f,f_ind);
            relmeanpow(f) = mean(relspect(f,1:2000),2);
            relpowerratio(f) = relpeakpow(f)./relmeanpow(f);
            
            relfftpeakpow(f) = relfftspect(f,f_ind);
            relfftmeanpow(f) = mean(relfftspect(f,1:2000));
            relfftpowerratio(f) = relfftpeakpow(f)./relfftmeanpow(f);
        end
    end
    bincondresp = bincondresp.*(1000/binwidth);
    binconderr = binconderr.*(1000/binwidth);
    bta = bta-prestim;
    
    figure
    plot(freqs(1:10),powerratio(1:10),'o-')
    hold on
%     plot(freqs(1:10),otherratio(1:10),'go-')
    plot(freqs(1:10),blratio(1:10),'ro-')
    legend('peak/all','peak/baseline')
    title('multi-taper spectrum')
    
    figure
    plot(freqs(1:10),fftpowerratio(1:10),'o-')
    hold on
%     plot(freqs(1:10),fftotherratio(1:10),'go-')
    plot(freqs(1:10),fftblratio(1:10),'ro-')
    legend('peak/all','peak/baseline')
    title('FFT spectrum')
    
    figure
    plot(freqs(1:10),relpowerratio(1:10),'o-');
    hold on
    plot(freqs(1:10),relpeakpow(1:10),'ro-');
    legend('peak/all','just normalizezd peak')
    title('mt -div by baseline spectrum')
    
    figure
    plot(freqs(1:10),relfftpowerratio(1:10),'o-');
    hold on
    plot(freqs(1:10),relfftpeakpow(1:10),'ro-');
    legend('peak/all','just normalizezd peak')
    title('fft -div by baseline spectrum')

    
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
        plot(j,condomi(i,j,1),'co')
    end
end
for i = find(cond)
    plot(condomi(i,:,1),'c')
end
axis([.5,5.5,-1.1,1.1])
line([.5,5.5],[0,0],'color','r')

disp('')
