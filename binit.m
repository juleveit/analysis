function [binned, ta] = binit(sig,binsize)

n = floor(length(sig)/binsize);

if isnan(sig); 
    binned = nan(1,n); 
    ta = nan(1,n); 
    return; 
end

for i = 1:n
    binned(i) = sum(sig(round((i-1)*binsize)+1:round(i*binsize)));
end
ta = linspace(1,n*binsize,n);