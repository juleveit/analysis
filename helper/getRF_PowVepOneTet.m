function [rf, ininds, faroutinds, field, snr] = getRF_PowVepOneTet(field, onoff, result, nstds, vep)

om = mean(mean(field,1),2);     %overall mean
long = reshape(field,1,numel(field));   %vector
sd = std(long);                 %std
thresh = om+nstds*sd;           %nstds sds above overall mean
rfraw = zeros(size(field));
rfraw(find(field>thresh)) = 1;
rf = bwareaopen(rfraw,3,4);

snr = mean(field(rf))/mean(field(find(rf == 0)));

field = 1;

fl = strel('square',7); %three pixels offset
expanded = imdilate(rf,fl);
farouts = find(expanded == 0);
rf = expanded+rf;

trunc = result.sizePixel-1;
allinds = [];
ins = find(rf==2);
%for every pixel in the rf
%find indices where this pixel was on/off
for i = 1:length(ins)
    [r,c] = ind2sub(size(rf),ins(i));
    r = r+trunc; c = c+trunc;
    inds = find(result.normimg(r,c,:) == onoff);
    for j = 1:length(inds)
        allinds = [allinds, find(result.stimulusIndex == inds(j))];
    end
end
ininds = unique(allinds);

allinds = [];
for i = 1:length(farouts)
    [r,c] = ind2sub(size(rf),farouts(i));
    r = r+trunc; c = c+trunc;
    inds = find(result.normimg(r,c,:) == -1 | result.normimg(r,c,:) == 1);
    for j = 1:length(inds)
        allinds = [allinds, find(result.stimulusIndex == inds(j))];
    end
end
faroutinds = unique(allinds);
