function [psd,fx] = lfppsd(signal,sr)

maxfreq = 150;
rssignal = resample(signal,1,round(sr/(maxfreq*2)));
[psd,fx] = pmtm(rssignal,3,[],maxfreq*2);

figure, semilogy(fx,psd);