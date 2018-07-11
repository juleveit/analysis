function hartley_ana

animalid = '180413';
block = 11;
lcol = 'r';

supath = ['C:\Users\Julia\work\data\' animalid '\singleunits\'];
basename = [animalid '_block' int2str(block) '_tet'];

files = dir([supath, basename, '*.mat']);

respwin = 20:160; % after stimulus onset

for cell = 1:length(files)
    
    load([supath, files(cell).name]);
%     disp(['now analyzing file: ' files(cell).name]);
    cellname{cell} = files(cell).name;
    
    i = strfind(files(cell).name, 'tet');
    tetno = strread(files(cell).name(i+3));
    
    wvchan = find(var(result.waveforms) == max(var(result.waveforms)));
    spike = result.waveforms(:,wvchan);
    interpspike = spline(1:32,spike,1:.1:32);
    [adiff(cell),swidth(cell)] = spikequant(interpspike);
    
    depth(cell) = result.depth;

    msStimes = round(result.spikes);
    if ~isempty(msStimes) & msStimes(1) == 0, msStimes(1) = 1; end
    
    chan = zeros(1,length(result.lfp));
    chan(msStimes) = 1;
    
    trialdur = round(result.FramesperStim*(1000/result.frameRate));
    msstamps = result.msstamps;
    for i = 1:length(msstamps)
        resp(i,:) = chan(msstamps(i):msstamps(i) + trialdur + 60); 
        lfpresp(i,:) = result.lfp(msstamps(i):msstamps(i) + trialdur + 60);
    end
   

    %IMPORTANT: CORRECT WRONG ORIENTATIONS TO CORRECT ONES
    result.gratingInfo.orientation = rem(abs(result.gratingInfo.orientation-450),360);
    result.gratingInfo.orientation(find(result.gratingInfo.orientation > 180)) = result.gratingInfo.orientation(find(result.gratingInfo.orientation > 180)) - 180;

    %normalize image
    result.stimulus = result.stimulus - min(min(min(result.stimulus)));
    result.stimulus = result.stimulus./max(max(max(result.stimulus)));

    stims = result.stimulus;
    nstims = size(stims,3);
    stimsize = size(stims,1);
    stims(:,:,nstims+1) = zeros(stimsize);
    stimulus = ones(1,length(chan)).*(nstims+1);
    for i = 1:length(result.msstamps)
        stimulus(result.msstamps(i):length(chan)) = result.stimulusIndex(i);
    end

    % to calculate map at different time points then chose the one at
    % biggest spatial variance
    i = 1;
    for t = 0:10:150
        chan(1:t) = 0;
        map(:,:,i) = mean(stims(:,:,stimulus(find(chan)-t)),3);
        i = i + 1;
    end
    for i = 1:10
        sv(i) = var(reshape(map(:,:,i),1,numel(map(:,:,i))));
    end
    
    figure
    subplot(1,2,1)
    imagesc(map(:,:,7))
    colormap gray
    axis square
    colorbar
    
    subplot(1,2,2)
    plot(0:10:90,sv)

    nstims = size(result.stimulus,3);

    numspks = sum(resp(:,respwin),2)./((length(respwin))/1000);
    lfppow = sqrt(mean(lfpresp(:,respwin).^2,2));
    lfppeak = mean(lfpresp(:,45:55),2).*-1; % todo find exact window of first neg going peak
    
%     if mean(numspks)<3 continue; end

    mp = zeros(nstims,result.repetitions);
    lfpmp = zeros(nstims,result.repetitions);
    for rep = 1:result.repetitions
        for i = 1:nstims
            mp(result.stimulusIndex((rep-1)*nstims+i),rep) = numspks((rep-1)*nstims+i);
            lfpmp(result.stimulusIndex((rep-1)*nstims+i),rep) = lfppow((rep-1)*nstims+i);
            lfppeakmp(result.stimulusIndex((rep-1)*nstims+i),rep) = lfppeak((rep-1)*nstims+i);
        end
    end
    
    weightedimg = zeros(size(result.stimulus));
    lfpwimg = zeros(size(result.stimulus));
    lfppeakwimg = zeros(size(result.stimulus));
    for i = 1:size(result.stimulus,1)
        for j = 1:size(result.stimulus,2)
            weightedimg(i,j,:) = squeeze(weightedimg(i,j,:)) + squeeze(result.stimulus(i,j,:)).*mean(mp,2);
            lfpwimg(i,j,:) = squeeze(lfpwimg(i,j,:)) + squeeze(result.stimulus(i,j,:)).*mean(lfpmp,2);
            lfppeakwimg(i,j,:) = squeeze(lfppeakwimg(i,j,:)) + squeeze(result.stimulus(i,j,:)).*mean(lfppeakmp,2);
        end
    end
    
    spikefield = squeeze(mean(weightedimg,3));
    lfpfield = squeeze(mean(lfpwimg,3));
    lfppeakfield = squeeze(mean(lfppeakwimg,3));

    degperElem = result.stimulusSize/size(result.stimulus,1);

    figure
    clf;
    subplot(2,2,1)
    plot(mean(resp(:,1:110)),'k','LineWidth',2);
    title(['Spikes FR: ' num2str(mean(numspks)) 'Hz'])
    
    subplot(2,2,2)
    imagesc(spikefield);
    axis square
    colormap gray
    colorbar
    
    subplot(2,2,3)
    plot(mean(lfpresp(:,1:110)),'k','LineWidth',2);
    title('LFP')
    
    subplot(2,2,4)
    imagesc(lfpfield);
    axis square
    colormap gray
    colorbar
    
%     plotSFTuning(numspks, result);
%     
%     plotORTuning(numspks, result);
    

    orientations = unique(result.gratingInfo.orientation);
    contrasts = unique(result.gratingInfo.contrast);
    xshifts = unique(result.gratingInfo.xshift);
    yshifts = unique(result.gratingInfo.yshift);
    phases = unique(result.gratingInfo.phase);
    spfreqs = unique(result.gratingInfo.spFreq);
    sfmap = zeros(length(spfreqs),size(weightedimg,1),size(weightedimg,2));
    for sf = 1:length(spfreqs)
        for ph = 1:length(phases)
            for o = 1:length(orientations)
                for c = 1:length(contrasts)
                    for x = 1:length(xshifts)
                        for y = 1:length(yshifts)
                            idx = find(result.gratingInfo.orientation == orientations(o) &...
                                    result.gratingInfo.contrast == contrasts(c) &...
                                    result.gratingInfo.xshift == xshifts(x) &...
                                    result.gratingInfo.yshift == yshifts(y) &...
                                    result.gratingInfo.phase == phases(ph) &...
                                    result.gratingInfo.spFreq == spfreqs(sf));
                            condspikes(sf,ph,o,c,x,y) = mean(numspks(idx));
                            spikeerr(sf,ph,o,c,x,y) = std(numspks(idx))./sqrt(length(idx));
                            condspikeresp(sf,ph,o,c,:) = mean(resp(idx,:));
                            condspikeall(sf,ph,o,c,:,:) = resp(idx,:);
                            condlfp(sf,ph,o,c,x,y) = mean(lfppow(idx));
                            lfperr(sf,ph,o,c,x,y) = std(lfppow(idx))./sqrt(length(idx));
                            condlfpresp(sf,ph,o,c,:) = mean(lfpresp(idx,:));
                            condlfpall(sf,ph,o,c,:,:) = lfpresp(idx,:);                            
                        end
                    end
                end
            end
        end
    end

    %orientation tuning:

    %spikes mean response OTI
    [spikeori,spiketunestr] = getOrientationPref(squeeze(mean(mean(mean(mean(mean(condspikes,1),2),4),5),6))',orientations);
    %find out preferred spatial frequency, phase and contrast
    spikeprefsf = find(mean(mean(mean(condspikes,2),3),4) == max(mean(mean(mean(condspikes,2),3),4)),1);
    spikeprefph = find(mean(mean(mean(condspikes,1),3),4) == max(mean(mean(mean(condspikes,1),3),4)),1);
    spikeprefc  = find(mean(mean(mean(condspikes,1),2),3) == max(mean(mean(mean(condspikes,1),2),3)),1);
    %spikes preferred response OTI
    [spikeprefori,spikepreftunestr] = getOrientationPref(squeeze(condspikes(spikeprefsf,spikeprefph,:,spikeprefc))',orientations);
    %fit wrapped gaussians to mean and preferred curves
    spikegaussparams = fitGauss(squeeze(mean(mean(mean(mean(mean(condspikes,1),2),4),5),6)),orientations');
    spikeprefgaussparams = fitGauss(squeeze(condspikes(spikeprefsf,spikeprefph,:,spikeprefc)),orientations');
    % spikechi2 = tuning_test(squeeze(mean(mean(mean(mean(mean(condspikes,1),2),4),5),6)));
    % spikeprefchi2 = tuning_test(squeeze(condspikes(spikeprefsf,spikeprefph,:,spikeprefc)));
    %spike preferred sf and c but mean over phases (pm = phase mean)
    [spikepmori,spikepmtunestr] = getOrientationPref(squeeze(mean(condspikes(spikeprefsf,:,:,spikeprefc),2))',orientations);
    spikepmgaussparams = fitGauss(squeeze(mean(condspikes(spikeprefsf,:,:,spikeprefc))),orientations');

    figure
%     subplot(2,2,1)
    errorbar(orientations,squeeze(mean(condspikes(spikeprefsf,:,:,spikeprefc),2)),squeeze(mean(spikeerr(spikeprefsf,:,:,spikeprefc),2)))
    hold on
    plot(wrappedgauss(spikepmgaussparams,1:180),'r')
    title(['orientation tuning: OSI = ' num2str(spikepmtunestr)]);
    
    
    [lfpori, lfptunestr] = getOrientationPref(squeeze(mean(mean(mean(mean(mean(condlfp,1),2),4),5),6))',orientations);
    %find out preferred spatial frequency, phase and contrast
    lfpprefsf = find(mean(mean(mean(condlfp,2),3),4) == max(mean(mean(mean(condlfp,2),3),4)),1);
    lfpprefph = find(mean(mean(mean(condlfp,1),3),4) == max(mean(mean(mean(condlfp,1),3),4)),1);
    lfpprefc  = find(mean(mean(mean(condlfp,1),2),3) == max(mean(mean(mean(condlfp,1),2),3)),1);
    %lfps preferred response OTI
    [lfpprefori,lfppreftunestr] = getOrientationPref(squeeze(condlfp(lfpprefsf,lfpprefph,:,lfpprefc))',orientations);
    %fit wrapped gaussians to mean and preferred curves
    lfpgaussparams = fitGauss(squeeze(mean(mean(mean(mean(mean(condlfp,1),2),4),5),6)),orientations');
    lfpprefgaussparams = fitGauss(squeeze(condlfp(lfpprefsf,lfpprefph,:,lfpprefc)),orientations');
    % lfpchi2 = tuning_test(squeeze(mean(mean(mean(mean(mean(condlfp,1),2),4),5),6)));
    % lfpprefchi2 = tuning_test(squeeze(condlfp(lfpprefsf,lfpprefph,:,lfpprefc)));
    [lfppmori,lfppmtunestr] = getOrientationPref(squeeze(mean(condlfp(lfpprefsf,:,:,lfpprefc),2))',orientations);
    lfppmgaussparams = fitGauss(squeeze(mean(condlfp(lfpprefsf,:,:,lfpprefc))),orientations');

    if length(contrasts)>=3
        nrspikeparams = fit_crf_NR([0; contrasts'],[spbl; squeeze(mean(mean(mean(condspikes,1),2),3))]);
        nrlfpparams = fit_crf_NR([0; contrasts'], [lfpbl; squeeze(mean(mean(mean(condlfp,1),2),3))]);
        nrgammaparams = fit_crf_NR([0; contrasts'],[gammabl; squeeze(mean(mean(mean(condgamma,1),2),3))]);
    else
        nrspikeparams = NaN;
        nrlfpparams = NaN;
        nrgammaparams = NaN;
    end

    stimsize = result.stimulusSize;
    position = result.position;
    gaussfactor = result.gratingInfo.gf == .5;  
    
end

disp('');
    
end


function respmat = cond_overview(resp, dim1, dim2, offs, col)

    % figure
    % hold on
    ndims = length(size(resp));
    meandims = 1:ndims; meandims([dim1,dim2,ndims]) = [];

    c = size(resp,dim1);
    o = size(resp,dim2);
    triallen = 110;
    
    meanresp = squeeze(mean(mean(resp,meandims(1)),meandims(2)));

    for i = 1:c
        for j = 1:o
            axes('position',[.1+(i-1)*(.8/o),.1+(j-1)*(.8/o),.8/o,.8/o])
            plot(squeeze(meanresp(i,j,1:triallen)));
            respmat(i,j,:) = squeeze(meanresp(i,j,1:triallen));
            axis([0,triallen,min(min(min(meanresp))),max(max(max(meanresp)))])
            set(gca,'XTick',[], 'XTickLabel',[])
            set(gca,'YTick',[], 'YTickLabel',[])
        end
    end
end

function [orimean, oritunestr] = getOrientationPref(r,o)
    %orimean: vector mean angle
    %tunestrength: vector mean length / sum of all responses
    
    theta = deg2rad(o);   
    %orientation tuning
    x = sum(cos(2*theta).*r);
    y = sum(sin(2*theta).*r);
    orimean = rad2deg(atan(y/x));
    if x<0, orimean = orimean+180; end
    if orimean<0, orimean=orimean+360; end
    orimean = orimean*.5;
    oritunestr = sqrt(x^2 + y^2)/(sum(r));
    
end


function [params,rmsd] = fitGauss(response,orientations)

% [mu, sigma, maxrate, spontrate]
orientations = mod(orientations,180);
[maxor,ind] = max(response);
gaufitmargin = 0.5;
range = max(response)-min(response);
p0 = [orientations(ind) 45 range min(response)];
lb = [0,0,(1-gaufitmargin)*range,(1-gaufitmargin)*min(response)];
ub = [180,180,(1+gaufitmargin)*range,(1+gaufitmargin)*min(response)];
if lb(3) == ub(3), ub(3) = lb(3)+1; end
if lb(4) == ub(4), ub(4) = lb(4)+1; end
if min(response) == 0, lb(4) = 0; ub(4) = 5; end
%keyboard
warning off
options=optimset('TolFun',1e-8,'Display','off');
[params,resnorm,residual,exitflag] = lsqcurvefit(@(p,x) wrappedgauss(p,x),p0,...
    orientations,response,lb,ub,options);
%     [0 0 mean(response(:)) 0],...
%     [180 90 2*max(response(:)) mean(response(:))],options);
%warning on
% mu = params(1);
% sigma = params(2);
% maxrate = params(3);
% spontrate = params(4);
rmsd = sqrt(mean(residual.^2));
end

function val = wrappedgauss(params,x)

N = 5;
for ix=1:length(x)
for n=-N:N
    tmp(n+N+1) = exp(-((x(ix)-params(1)+180*n).^2)./(2*params(2)^2));
end
val(ix) = sum(tmp);
end
%val=val/sum(val);
val=params(4) + params(3)*val(:);
end


function params=fit_crf_NR(x,y)
% [rmax, n, c50, r0]

nrfitmargin = .01;
range = max(y)-min(y);
p0 = [range 2 .5 min(y)];
lb = [(1-nrfitmargin)*range, .5, 0, min(y)-.1*range];
ub = [(1+nrfitmargin)*range, 10, 1, min(y)+.1*range];
if lb(1) == ub(1), ub(1) = lb(1)+1; end
if lb(4) == ub(4), ub(4) = lb(4)+1; end
warning off
params = lsqcurvefit(@(p,x) NakaRushton(p,x),p0,x(:),y(:),lb,ub,optimset('Display','off'));
warning on
end

function val=NakaRushton(p,x)
% parameters of Naka-Rushton function as in Disney et al., Neuron, 2007
% [R_max, contrast Exponent n,  50%firing-Contrast, spontaneous rate sFR]

val = p(4)+p(1)*((x.^p(2))./(x.^p(2)+p(3).^p(2)));
end


% directed angular difference for 180 deg
function a = oridiff(x,y)
    a = x-y;
    a(find(a>90)) = -abs(a(find(a>90))-180);
    a(find(a<-90)) = abs(abs(a(find(a<-90)))-180);
end
