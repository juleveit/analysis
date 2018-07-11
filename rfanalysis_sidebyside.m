function rfanalysis_sidebyside

animalid = '171205';
blocks = [3,5];

normalize_plot = 0;

supath = ['C:\Users\Julia\work\data\' animalid '\singleunits\'];
basename1 = [animalid '_block' int2str(blocks(1)) '_tet'];
basename2 = [animalid '_block' int2str(blocks(2)) '_tet'];
files1 = dir([supath, basename1, '*.mat']);
files2 = dir([supath, basename2, '*.mat']);


if length(files1)~=length(files2)
    disp('not corresponding files')
    pause;
end

for cell = 1:length(files1)
    load([supath, files1(cell).name]);
    result1 = result;
    load([supath, files2(cell).name]);
    result2 = result;
    
    if strfind(files1(cell).name, 'MU')
        continue;
    end
    
    i = strfind(files1(cell).name, 'tet');
    tetno = strread(files1(cell).name(i+3));
    
    msStimes1 = round(result1.spikes);
    msStimes2 = round(result2.spikes);
    if ~isempty(msStimes1)
        if msStimes1(1) == 0, msStimes1(1) = 1; end
    end
    if ~isempty(msStimes2)
        if msStimes2(1) == 0, msStimes2(1) = 1; end
    end
    
    normimg1 = result1.stimulus;
    normimg1(find(normimg1 == 0)) = -1;
    normimg1(find(normimg1 == 128)) = 0;
    normimg1(find(normimg1 == 255)) = 1;
    result1.normimg = normimg1;
    normimg2 = result2.stimulus;
    normimg2(find(normimg2 == 0)) = -1;
    normimg2(find(normimg2 == 128)) = 0;
    normimg2(find(normimg2 == 255)) = 1;
    result2.normimg = normimg2;
    trialdur = (1000/result1.frameRate)*result2.FramesperStim;

    chan1 = zeros(1,length(result1.lfp));
    chan1(msStimes1) = 1;
    chan2 = zeros(1,length(result2.lfp));
    chan2(msStimes2) = 1;

    sr = 1000;
    %window for calculating spikerates
    srwinoffs = 20;
    srwinsize = 79;
    trialstart = 0;
    trialend = 2*round(trialdur);
    triallen = trialend-trialstart;
    trialnsamps = round((trialdur*sr)/1000);

    msstamps1 = result1.msstamps;
    msstamps2 = result2.msstamps;
    for trial = 1:length(msstamps1)
        resp1(trial,:) = chan1(msstamps1(trial):msstamps1(trial)+triallen-1);
        resp2(trial,:) = chan2(msstamps2(trial):msstamps2(trial)+triallen-1);
    end

    nspikes1 = sum(resp1(:,srwinoffs:srwinoffs+srwinsize),2).*(1000/srwinsize+1);
    nspikes2 = sum(resp2(:,srwinoffs:srwinoffs+srwinsize),2).*(1000/srwinsize+1);

    nstims = size(result1.stimulus,3);
    mp1 = zeros(nstims,result1.repetitions);
    mp2 = zeros(nstims,result1.repetitions);
    for rep = 1:result1.repetitions
        for i = 1:nstims
            mp1(result1.stimulusIndex((rep-1)*nstims+i),rep) = nspikes1((rep-1)*nstims+i);
            mp2(result2.stimulusIndex((rep-1)*nstims+i),rep) = nspikes2((rep-1)*nstims+i);
        end
    end

    weightedimg1 = zeros(size(result1.stimulus));
    weightedimg2 = zeros(size(result2.stimulus));
    for rep = 1:result1.repetitions
        for i = 1:size(result1.stimulus,1)
            for j = 1:size(result1.stimulus,2)
                weightedimg1(i,j,:) = squeeze(weightedimg1(i,j,:)) + squeeze(normimg1(i,j,:)).*mp1(:,rep);
                weightedimg2(i,j,:) = squeeze(weightedimg2(i,j,:)) + squeeze(normimg2(i,j,:)).*mp2(:,rep);
            end
        end
    end
    weightedimg1 = weightedimg1./result1.repetitions;
    pos1 = weightedimg1; neg1 = weightedimg1;
    pos1(find(weightedimg1<0)) = 0;
    neg1(find(weightedimg1>0)) = 0;
    weightedimg2 = weightedimg2./result2.repetitions;
    pos2 = weightedimg2; neg2 = weightedimg2;
    pos2(find(weightedimg2<0)) = 0;
    neg2(find(weightedimg2>0)) = 0;

    trunc = result1.sizePixel-1;
    a = trunc+1;
    b = size(pos1,1)-trunc;
    degperElem = (result1.shownSize*result1.sizePixel)/size(result1.stimulus,1);
    nelem = result1.stimulusElements;
    degImg = result1.stimulusSize;
    ax = linspace(1,degImg,nelem);
    degperGridPos = degImg/nelem;
    
    off1 = abs(squeeze(sum(neg1(a:b,a:b,:),3)));
    on1 = abs(squeeze(sum(pos1(a:b,a:b,:),3)));
    off1 = off1./(result1.sizePixel)^2; %not yet compensated for!!
    on1 = on1./(result1.sizePixel)^2;
    if normalize_plot
        off1 = off1 - min(min(off1));
        off1 = off1./max(max(off1));
        on1 = on1-min(min(on1));
        on1 = on1./max(max(on1));
    end
    off2 = abs(squeeze(sum(neg2(a:b,a:b,:),3)));
    on2 = abs(squeeze(sum(pos2(a:b,a:b,:),3)));
    off2 = off2./(result2.sizePixel)^2; %not yet compensated for!!
    on2 = on2./(result2.sizePixel)^2;
    if normalize_plot
        off2 = off2 - min(min(off2));
        off2 = off2./max(max(off2));
        on2 = on2-min(min(on2));
        on2 = on2./max(max(on2));
    end

    % fit with Gauss
    % OFF field
    %put to spikes per window
    nstds = 1;
    off1 = off1./(result1.sizePixel)^2; %not yet compensated for!!
    [rfoff1, insoff1, outsoff1, fieldoff1, snroff1] = getRF_PowVepOneTet(off1,-1,result1,nstds,0); % old leftover, but get the indices of pixels inside the field for background subtraction
    fit_gauss = 1;
    if fit_gauss
        %1 = center x, 2 = spatial spread, 3 = y offset, 4 = scaling (3+4 =max)
        gaussfitoff1 = fit_or2dgauss(off1,degperGridPos,1); % fit the gaussian to the response map
    end
%     inindsoff = find(rfoff==2);
    
    %ON field
    %put to spikes per window
    on1 = on1./(result1.sizePixel)^2; %not yet compensated for!!
    [rfon1, inson1, outson1, fieldon1, snron1] = getRF_PowVepOneTet(on1,1,result1,nstds,0);
    if fit_gauss
        gaussfiton1 = fit_or2dgauss(on1,degperGridPos, 1);
    end
%     inindson = find(rfon==2);
    
    off2 = off2./(result2.sizePixel)^2; %not yet compensated for!!
    [rfoff2, insoff2, outsoff2, fieldoff2, snroff2] = getRF_PowVepOneTet(off2,-1,result2,nstds,0); 
    fit_gauss = 1;
    if fit_gauss
        %1 = center x, 2 = spatial spread, 3 = y offset, 4 = scaling (3+4 =max)
        gaussfitoff2 = fit_or2dgauss(off2,degperGridPos,1); % fit the gaussian to the response map
    end
    on2 = on2./(result2.sizePixel)^2; %not yet compensated for!!
    [rfon2, inson2, outson2, fieldon2, snron2] = getRF_PowVepOneTet(on2,1,result2,nstds,0);
    if fit_gauss
        gaussfiton2 = fit_or2dgauss(on2,degperGridPos, 1);
    end

    a = ax*degperGridPos;
    w = result1.stimulusSize/2;
    xax = linspace(result1.position(1)-w,result1.position(1)+w,result1.stimulusElements);
    yax = linspace(result1.position(2)-w,result1.position(2)+w,result1.stimulusElements);
    
    mx = max([max(max(on1)),max(max(off1)),max(max(on2)),max(max(off2))]);
    mn = min([min(min(on1)),min(min(off1)),min(min(on2)),min(min(off2))]);
    
    figure
    subplot(2,2,1)
%     imagesc(xax,yax,on);
    imagesc(on1)
    hold on
    plot_orrf(gaussfiton1,1,'w',2)
    axis square
    axis xy
    colorbar
    colormap jet
    caxis([mn,mx])
    title(['cell ' files1(cell).name '  ON field'])

    subplot(2,2,2)
%     imagesc(xax,yax,off);
    imagesc(off1)
    hold on
    plot_orrf(gaussfitoff1,1,'w',2)
    axis('square')
    axis xy
    colorbar;
    colormap jet
    caxis([mn,mx])
    title(['OFF field depth: ' int2str(result1.depth)])    
    
    subplot(2,2,3)
%     imagesc(xax,yax,on);
    imagesc(on2)
    hold on
    plot_orrf(gaussfiton2,1,'w',2)
    axis square
    axis xy
    colorbar
    colormap jet
    caxis([mn,mx])
    title(['cell ' files2(cell).name '  ON field'])

    subplot(2,2,4)
%     imagesc(xax,yax,off);
    imagesc(off2)
    hold on
    plot_orrf(gaussfitoff2,1,'w',2)
    axis('square')
    axis xy
    colorbar;
    colormap jet
    caxis([mn,mx])
    title(['OFF field depth: ' int2str(result2.depth)])
    
end

disp('');
