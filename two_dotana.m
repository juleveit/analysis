function two_dotana

animalid = '150213';
block = 11;

normalize_plot = 0;

supath = ['C:\Users\Julia\work\data\' animalid '\singleunits\'];
basename = [animalid '_block' int2str(block) '_tet'];
files = dir([supath, basename, '*.mat']);



for cell = 3:length(files)
    load([supath, files(cell).name]);
%     disp(['now analyzing file: ' files(cell).name]);
    
    spkamps = max(result.waveforms) - min(result.waveforms);
    wvchan = find(spkamps == max(spkamps));
    spike = result.waveforms(:,wvchan); 
    interpspike = spline(1:32,spike,1:.1:32);
    [adiff(cell),swidth(cell)] = spikequant(interpspike);

    sr = 1000;
    lfp = result.lfp(:,wvchan)';

    i = strfind(files(cell).name, 'tet');
    if strcmp(files(cell).name(i+4),'_')
        tetno = strread(files(cell).name(i+3)); % single character number
    else        
        tetno = strread(files(cell).name(i+3:i+4)); % number >10
    end
    if tetno>8
        v1(cell) = 0; v2(cell) = 1; cellstr = 'V2';
    else
        v1(cell) = 1; v2(cell) = 0; cellstr = 'V1';
    end

    msStimes = round(result.spikes);
    if ~isempty(msStimes)
        if msStimes(1) == 0, msStimes(1) = 1; end
    end
    
    normimg = result.stimulus;
    normimg(find(normimg == 0)) = -1;
    normimg(find(normimg == 128)) = 0;
    normimg(find(normimg == 255)) = 1;
    trialdur = (1000/result.frameRate)*result.FramesperStim;

    chan = zeros(1,length(result.lfp));
    chan(msStimes) = 1;

    sr = 1000;
    %window for calculating spikerates
    srwinoffs = 60;
    srwinsize = 79;
    trialstart = -0;
    trialend = round(trialdur);
    triallen = trialend-trialstart;
    trialnsamps = round((trialdur*sr)/1000);

    msstamps = result.msstamps;
    for trial = 1:length(msstamps)
        resp(trial,:) = chan(msstamps(trial):msstamps(trial)+triallen-1);
        lfpresp(trial,:) = lfp(msstamps(trial):msstamps(trial)+triallen-1);
    end

    
    nspikes = sum(resp(:,srwinoffs:srwinoffs+srwinsize),2).*(1000/srwinsize+1);
    tmv = find(var(lfpresp) == max(var(lfpresp(:,srwinoffs:srwinoffs+srwinsize))));
    mlfp = lfpresp(:,tmv);

    nstims = size(result.stimulus,3);
    mp = zeros(nstims,result.repetitions);
    lfpT = zeros(nstims,result.repetitions);
    for rep = 1:result.repetitions
        for i = 1:nstims
            mp(result.stimulusIndex((rep-1)*nstims+i),rep) = nspikes((rep-1)*nstims+i);
            lfpT(result.stimulusIndex((rep-1)*nstims+i),rep) = mlfp((rep-1)*nstims+i);     
        end
    end

    weightedimg = zeros(size(result.stimulus));
    weightedlfpimg = zeros(size(result.stimulus));
    for rep = 1:result.repetitions
        for i = 1:size(result.stimulus,1)
            for j = 1:size(result.stimulus,2)
                weightedimg(i,j,:) = squeeze(weightedimg(i,j,:)) + squeeze(normimg(i,j,:)).*mp(:,rep);
                weightedlfpimg(i,j,:) = squeeze(weightedlfpimg(i,j,:)) + squeeze(normimg(i,j,:)).*lfpT(:,rep);
            end
        end
    end
    weightedimg = weightedimg./result.repetitions;
    weightedlfpimg = weightedlfpimg./result.repetitions;
    pos = weightedimg; neg = weightedimg;
    pos(find(weightedimg<0)) = 0;
    neg(find(weightedimg>0)) = 0;
    lfppos = weightedlfpimg; lfpneg = weightedlfpimg;
    lfppos(find(weightedlfpimg<0)) = 0;
    lfpneg(find(weightedlfpimg>0)) = 0;

    trunc = result.sizePixel-1;
    a = trunc+1;
    b = size(pos,1)-trunc;
    degperElem = (result.shownSize*result.sizePixel)/size(result.stimulus,1);
    nelem = result.stimulusElements;
    degImg = result.stimulusSize;
    ax = linspace(1,degImg,nelem);
    degperGridPos = degImg/nelem;
    off = abs(squeeze(sum(neg(a:b,a:b,:),3)));
    on = abs(squeeze(sum(pos(a:b,a:b,:),3)));
    off = off./(result.sizePixel)^2; %not yet compensated for!!
    on = on./(result.sizePixel)^2;
    lfpoff = abs(squeeze(mean(lfpneg(a:b,a:b,:),3)));
    lfpon = abs(squeeze(mean(lfppos(a:b,a:b,:),3)));
    if normalize_plot
        off = off - min(min(off));
        off = off./max(max(off));
        on = on-min(min(on));
        on = on./max(max(on));
    end

    % fit with Gauss
    % OFF field
    %put to spikes per window
    off = off./(result.sizePixel)^2; %not yet compensated for!!
%     [rfoff, insoff, outsoff, fieldoff, snroff] = getRF_PowVepOneTet(off,-1,result,nstds,vep); % old leftover, but get the indices of pixels inside the field for background subtraction
    fit_gauss = 1;
    if fit_gauss
        %1 = center x, 2 = spatial spread, 3 = y offset, 4 = scaling (3+4 =max)
        gaussfitoff = fit_or2dgauss(off,degperGridPos,1); % fit the gaussian to the response map
    end
%     inindsoff = find(rfoff==2);
    
    %ON field
    %put to spikes per window
    on = on./(result.sizePixel)^2; %not yet compensated for!!
%     [rfon, inson, outson, fieldon, snron] = getRF_PowVepOneTet(on,1,result,nstds,vep);
    if fit_gauss
        gaussfiton = fit_or2dgauss(on,degperGridPos, 1);
    end
%     inindson = find(rfon==2);
    
%     farouts = intersect(outson, outsoff); % the pixel indices far from both ON and OFF fields
    %%%%
    
    a = ax*degperGridPos;
    w = result.stimulusSize/2;
    xax = linspace(result.position(1)-w,result.position(1)+w,result.stimulusElements);
    yax = linspace(result.position(2)-w,result.position(2)+w,result.stimulusElements);
    
    
    figure
    subplot(2,2,1)
%     imagesc(xax,yax,on);
    imagesc(on)
    hold on
    plot_orrf(gaussfiton,1,'w',2)
    axis square
    axis xy
    colorbar
    title([cellstr ' cell ' files(cell).name '  ON field'])

    subplot(2,2,2)
%     imagesc(xax,yax,off);
    imagesc(off)
    hold on
    plot_orrf(gaussfitoff,1,'w',2)
    axis('square')
    axis xy
    colorbar;
    title([cellstr ' OFF field depth: ' int2str(result.depth)])
    
    subplot(2,2,3)
    imagesc(xax,yax,on-off);
    axis square;
    colorbar;
    axis xy
    title('ON - OFF')
    
    subplot(2,2,4)
    plot(result.waveforms(:,wvchan))
    title(['swidth: ' num2str(swidth(cell))])
    
    figure    
    subplot(2,2,1)
%     imagesc(xax,yax,on);
    imagesc(lfpon)
%     hold on
%     plot_orrf(gaussfiton,1,'w',2)
    axis square
    axis xy
    colorbar
    title([cellstr ' cell ' files(cell).name '  ON field'])

    subplot(2,2,2)
%     imagesc(xax,yax,off);
    imagesc(lfpoff)
%     hold on
%     plot_orrf(gaussfitoff,1,'w',2)
    axis('square')
    axis xy
    colorbar;
    title([cellstr ' OFF field depth: ' int2str(result.depth)])
    
    subplot(2,2,3)
    imagesc(xax,yax,lfpon-lfpoff);
    axis square;
    colorbar;
    axis xy
    title('ON - OFF')
    
    
    onsnr(cell) = gaussfiton.amp/gaussfiton.offs;
    offsnr(cell) = gaussfitoff.amp/gaussfitoff.offs;
    
    oncenter(cell,:) = [gaussfiton.xcenterDeg,gaussfiton.ycenterDeg];
    offcenter(cell,:) = [gaussfitoff.xcenterDeg,gaussfitoff.ycenterDeg];
    
    onspread(cell) = mean([gaussfiton.xspreadDeg,gaussfiton.yspreadDeg]);
    offspread(cell) = mean([gaussfitoff.xspreadDeg,gaussfitoff.yspreadDeg]);
 
    disp('')
end

disp('');
