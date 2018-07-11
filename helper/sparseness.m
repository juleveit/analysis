function s = sparseness(r)

n = length(r);
s = (1 - (((sum(r)/n)^2) / (sum((r.^2)./n))) ) / (1-(1/n));