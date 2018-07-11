function [onfield,offfield,onfit,offfit,rson,rsoff,xax,yax] = get_rf_MUA(filename)

load(filename);

if ~exist('result') & exist('centresult')
    result = centresult;
else
    result = surrresult;
end

sr = 1000;

normimg = result.stimulus;
normimg(find(normimg == 0)) = -1;
normimg(find(normimg == 128)) = 0;
normimg(find(normimg == 255)) = 1;
trialdur = (1000/result.frameRate)*result.FramesperStim;

%window for calculating spikerates
srwinoffs = 60;
srwinsize = 79;
trialstart = -0;
trialend = 2*round(trialdur);
triallen = trialend-trialstart;
trialnsamps = round((trialdur*sr)/1000);

msstamps = result.msstamps;

for ch = 1:size(result.lfp,1)
    lfp = result.lfp(ch,:);

    msStimes = round(result.msStimes{ch});
    if ~isempty(msStimes)
        if msStimes(1) == 0, msStimes(1) = 1; end
    end
    
    chan = zeros(1,length(lfp));
    chan(msStimes) = 1;

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

    w = result.stimulusSize/2;
    % fit with Gauss
    % OFF field
    %put to spikes per window
    off = off./(result.sizePixel)^2; %not yet compensated for!!
    [gaussfitoff,rsoff(ch)] = fit_or2dgauss(off,degperGridPos,1); % fit the gaussian to the response map
    gaussfitoff.shiftxcenterDeg = gaussfitoff.xcenterDeg + result.position(1) - w;
    gaussfitoff.shiftycenterDeg = gaussfitoff.ycenterDeg + result.position(2) - w;
    %     inindsoff = find(rfoff==2);

    %ON field
    %put to spikes per window
    on = on./(result.sizePixel)^2; %not yet compensated for!!
    [gaussfiton,rson(ch)] = fit_or2dgauss(on,degperGridPos, 1);
    gaussfiton.shiftxcenterDeg = gaussfiton.xcenterDeg + result.position(1) - w;
    gaussfiton.shiftycenterDeg = gaussfiton.ycenterDeg + result.position(2) - w;

    a = ax*degperGridPos;
    xax = linspace(result.position(1)-w,result.position(1)+w,result.stimulusElements);
    yax = linspace(result.position(2)-w,result.position(2)+w,result.stimulusElements);
    
    onfield(:,:,ch) = on; offfield(:,:,ch) = off; 
    onfit(ch) = gaussfiton; offfit(ch) = gaussfitoff;
    
end


