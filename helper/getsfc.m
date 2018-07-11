%returns spike triggered average, average spectrum of lfp snippets,
%spectrum of spike trigered average, spike field coherence, the frequency
%axis and the number of spikes used to calculate the sta
function [sta, avgsnipspectrum, staspectrum, sfcoher, sfcfax, nspikes, snippets] = getsfc(chan, lfp, win, sr)

    if isempty(find(chan(1+win:length(chan)-win),1))    % if there is no spikes, fill everything with NaNs   
        slen = (2^nextpow2(2*win+1))/2+1;
        sta = nan(1,2*win+1);
        snippets = nan(1,2*win+1);
        avgsnipspectrum = nan(1,slen);
        staspectrum = nan(1,slen);
        sfcoher = nan(1,slen);
        nspikes = 0;
        sfcfax = nan(1,slen);
    else        
        sm = zeros(1,2*win+1);
        snippets = zeros(1,2*win+1);
        nfft = 2^nextpow2(length(sm));
        j = 1; %spike counter
        spks = find(chan(win+1:length(chan)-win))+win;
        for i = 1:length(spks)
            snippets(i,:) = lfp(spks(i)-win:spks(i)+win);
            [snippow(i,:),sfcfax] = pmtm(lfp(spks(i)-win:spks(i)+win),3,nfft,sr);
        end
        nspikes = size(snippow,1);
        sta = mean(snippets,1);
        
%         for i = win+1:length(chan)-win 
%             if chan(i) % if there is a spike add lfp snippet around to average
%                 sm = sm+chan(i).*lfp(i-win:i+win);
%                 [snippow(j,:),sfcfax] = pmtm(lfp(i-win:i+win),3,nfft,sr); % get spectrum of that snippet
%                 j = j+1;
%             end
%         end
% 
%         nspikes = j-1;
%         sta = sm./nspikes; % normalize sta by number of spikes

        avgsnipspectrum = mean(snippow,1);  % average snippet spectra
        staspectrum = pmtm(sta,3,nfft,sr);  % get spectrum of sta
        sfcoher = staspectrum'./avgsnipspectrum; % coherence
        if nspikes == 1 % if there is only one spike sfc is just a long vector of ones
            sfcoher = nan(1,length(staspectrum));
        end
    end
end
