function [x,y,xerr] = runningMedian(depth,xval,uncert,avgn)

randshifts = randn(1,length(depth))*(uncert/2);
depth = depth+randshifts;
[sorteddepth,idx] = sort(depth);
sortedratio = xval(idx);
for i = 1:length(sortedratio)-avgn
    x(i) = nanmedian(sortedratio(i:i+avgn-1));
    xerr(i) = nanstd(sortedratio(i:i+avgn-1))/sqrt(avgn);
    y(i) = mean(sorteddepth(i:i+avgn-1));
end