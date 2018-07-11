function [adiff,swidth,ptr,endslope,fwhh] = spikequant(spike)

t = find(spike == min(spike));
leftmax = find(spike(1:t) == max(spike(1:t)));
rightmax = find(spike(t:end) == max(spike(t:end))) + t -1;
adiff = (spike(rightmax)-spike(leftmax))/(spike(leftmax)+spike(rightmax)); %diff amplitude of positive peaks
swidth = rightmax-t;    % width trough to right max
at = spike(t);          % amplitude at through
ptr = spike(rightmax)/at;      % peak trough ratio
% endslope = diff(spike(30:31));
endslope = mean(diff(spike(250:310)));
%try width at half height
hh = (at-(spike(rightmax)+spike(leftmax))/2)/2;
hhy = at-hh;
hhx1 = find(spike(1:t)<hhy,1);
hhx2 = find(spike(t:end)>hhy,1)+t-1;
fwhh = hhx2-hhx1;