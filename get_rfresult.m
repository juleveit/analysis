function rfresult = get_rfresult(filename)

load(filename);

rfresult.animalid = result.animalid;
rfresult.position = result.position;
[p,rfresult.cellname] = fileparts(filename);

spkamps = max(result.waveforms) - min(result.waveforms);
wvchan = find(spkamps == max(spkamps));
spike = result.waveforms(:,wvchan);

sr = 1000;
lfp = result.lfp(:,wvchan)';

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
trialend = 2*round(trialdur);
triallen = trialend-trialstart;
trialnsamps = round((trialdur*sr)/1000);

msstamps = result.msstamps;
for trial = 1:length(msstamps)
    resp(trial,:) = chan(msstamps(trial):msstamps(trial)+triallen-1);
    lfpresp(trial,:) = lfp(msstamps(trial):msstamps(trial)+triallen-1);
end

frs = sum(resp(:,srwinoffs:srwinoffs+srwinsize),2).*(1000/srwinsize+1);
% tmv = find(var(lfpresp) == max(var(lfpresp(:,srwinoffs:srwinoffs+srwinsize))));
% hlp = mean(lfpresp,1); mn = min(hlp(srwinoffs:srwinoffs+srwinsize)); tmv = find(hlp==mn);
% mlfp = lfpresp(:,tmv);
mlfp = mean(lfpresp(:,srwinoffs+20:srwinoffs+50),2);

nstims = size(result.stimulus,3);
mp = zeros(nstims,result.repetitions);
lfpT = zeros(nstims,result.repetitions);
for rep = 1:result.repetitions
    for i = 1:nstims
        mp(result.stimulusIndex((rep-1)*nstims+i),rep) = frs((rep-1)*nstims+i);
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
lfppos(find(normimg<0)) = 0;
lfpneg(find(normimg>0)) = 0;

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
lfpoff = squeeze(mean(-lfpneg(a:b,a:b,:),3));
lfpon = squeeze(mean(lfppos(a:b,a:b,:),3));

[threshrfon, inindson, outindson, snron] = getRFinds(on, 1, result, 2, normimg);
[threshrfoff, inindsoff, outindsoff, snroff] = getRFinds(off, -1, result, 2, normimg);
[threshlfprfon, lfpinindson, lfpoutindson, lfpsnron] = getRFinds(-lfpon, 1, result, 2, normimg);
[threshlfprfoff, lfpinindsoff, lfpoutindsoff, lfpsnroff] = getRFinds(-lfpoff, -1, result, 2, normimg);

w = result.stimulusSize/2;
% fit with Gauss
% OFF field
%put to spikes per window
off = off./(result.sizePixel)^2; %not yet compensated for!!
[gaussfitoff,rsoff] = fit_or2dgauss(off,degperGridPos,1); % fit the gaussian to the response map
gaussfitoff.shiftxcenterDeg = gaussfitoff.xcenterDeg + result.position(1) - w;
gaussfitoff.shiftycenterDeg = gaussfitoff.ycenterDeg + result.position(2) - w;
%     inindsoff = find(rfoff==2);

%ON field
%put to spikes per window
on = on./(result.sizePixel)^2; %not yet compensated for!!
[gaussfiton,rson] = fit_or2dgauss(on,degperGridPos, 1);
gaussfiton.shiftxcenterDeg = gaussfiton.xcenterDeg + result.position(1) - w;
gaussfiton.shiftycenterDeg = gaussfiton.ycenterDeg + result.position(2) - w;

% LFP fields
[lfpfiton,lfprson] = fit_or2dgauss(-lfpon,degperGridPos, 1);
lfpfiton.shiftxcenterDeg = lfpfiton.xcenterDeg + result.position(1) - w;
lfpfiton.shiftycenterDeg = lfpfiton.ycenterDeg + result.position(2) - w;

[lfpfitoff,lfprsoff] = fit_or2dgauss(-lfpoff,degperGridPos, 1);
lfpfitoff.shiftxcenterDeg = lfpfitoff.xcenterDeg + result.position(1) - w;
lfpfitoff.shiftycenterDeg = lfpfitoff.ycenterDeg + result.position(2) - w;

a = ax*degperGridPos;
xax = linspace(result.position(1)-w,result.position(1)+w,result.stimulusElements);
yax = linspace(result.position(2)-w,result.position(2)+w,result.stimulusElements);

rfresult.on = on; rfresult.off = off; 
rfresult.gaussfiton = gaussfiton; rfresult.gaussfitoff = gaussfitoff;
rfresult.rson = rson; rfresult.rsoff = rsoff;
rfresult.xax = xax; rfresult.yax = yax;
rfresult.snron = snron; rfresult.snroff = snroff;
rfresult.fron = mean(frs(inindson)); rfresult.froff = mean(frs(inindsoff));

rfresult.lfpon = lfpon; rfresult.lfpoff = lfpoff;
rfresult.lfpfiton = lfpfiton; rfresult.lfpfitoff = lfpfitoff;
rfresult.lfprson = lfprson; rfresult.lfprsoff = rsoff;
rfresult.lfpsnron = lfpsnron; rfresult.lfpsnroff = snroff;


% % now make the big figures
% for r = 1:size(normimg,1)-2
%     posspikevec = []; negspikevec = [];
%     poslfpvec = []; neglfpvec = [];
%     for c = 1:size(normimg,2)-2
%         posinds = []; neginds = [];
%         posstims = find(normimg(r+1,c+1,:) == 1); negstims = find(normimg(r+1,c+1,:) == -1);
%         for i = 1:length(posstims)
%             posinds = [posinds, find(result.stimulusIndex == posstims(i))];
%             neginds = [neginds, find(result.stimulusIndex == negstims(i))];
%         end
%         posspikevec = [posspikevec, mean(resp(posinds,srwinoffs:srwinoffs+srwinsize))];
%         negspikevec = [negspikevec, mean(resp(neginds,srwinoffs:srwinoffs+srwinsize))];
%         poslfpvec = [poslfpvec, mean(lfpresp(posinds,:))];
%         neglfpvec = [neglfpvec, mean(lfpresp(neginds,:))];
%     end
%     posmat(r,:) = posspikevec; negmat(r,:) = negspikevec;
%     poslfpmat(r,:) = poslfpvec; neglfpmat(r,:) = neglfpvec;
%     poslfpmx(r) = max(poslfpmat(r,:)); poslfpmn(r) = min(poslfpmat(r,:));
%     neglfpmx(r) = max(neglfpmat(r,:)); neglfpmn(r) = min(neglfpmat(r,:));
%     poslfpamp = max(poslfpmx-poslfpmn); neglfpamp = max(neglfpmx-neglfpmn);
% end
% 
% mx = max(max(posmat));
% figure; hold on;
% for i = 1:size(posmat,1)
%     plot(posmat(i,:)+(i-1)*mx)
%     line([(i-1)*srwinsize,(i-1)*srwinsize],[0,(size(posmat,1)-1)*mx]);
% end
% mx = max(max(negmat));
% figure; hold on;
% for i = 1:size(negmat,1)
%     plot(negmat(i,:)+(i-1)*mx)
%     line([(i-1)*srwinsize,(i-1)*srwinsize],[0,(size(negmat,1)-1)*mx]);
% end
% 
% figure; hold on;
% for i = 1:size(poslfpmat,1)
%     plot(poslfpmat(i,:)+(i-1)*poslfpamp)
%     line([(i-1)*166,(i-1)*166],[0,(size(posmat,1)-1)*poslfpamp]);
%     line([0,1660],[(i-1)*poslfpamp, (i-1)*poslfpamp]);
% end
% figure; hold on;
% for i = 1:size(neglfpmat,1)
%     plot(neglfpmat(i,:)+(i-1)*neglfpamp)
%     line([(i-1)*166,(i-1)*166],[0,(size(negmat,1)-1)*neglfpamp]);
%     line([0,1660],[(i-1)*neglfpamp, (i-1)*neglfpamp]);
% end

