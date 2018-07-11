function [ftspect,ftphases,ftfax,ftralp,ftppc,ftplv] = get_ft_spectstats(lfpresp,resp,trials)
% takes:   lfpresp: matrix of lfp trials (ntrials,nsamples) 
%          resp: matrix of spikes (ntrials,nsamples) all zeros and ones for spikes
%          trials: vector of indices into the rows - can be 1:size(lfpresp,1)to do all
% returns: ftspect - LFP spectrum
%          ftphases - matrix of spike phases per frequency bin
%          ftfax - frequency axis
%          ftralp - rayleigh p-value vector (per frequency)
%          ftppc - ppc spectrum
%          ftplv - phase locking value
% many parameters need to be adjusted in the code below

     if isempty(trials) || isempty(find(resp(trials,801:1800))) % adjust time window from which to take spikes
         ftspect = nan(1,30); ftphases = NaN; ftfax = nan(1,30); % adjust length of vectors for different spectrum lengths
         ftralp = nan(1,30); ftppc = nan(1,30); ftplv = nan(1,30);
     else
         clear data;
         for i = 1:length(trials)
             data.trial{i} = [lfpresp(trials(i),:);resp(trials(i),:)];
             data.time{i} = -.299:.001:2.700; % adjust time axis (in seconds)
         end
         data.fsample = 1000; % adjust if necessary - sampling rate
         cfg = [];
         data.label{1,1} = 'lfp';
         data.label{2,1} = 'spikes';
         cfg.spikechannel = 'spikes';
         cfg.channel = 'lfp';
         cfg.latency = [0.5,1.5]; % this is the time window for which spikes to take in seconds!!

         cfg.method = 'mtmfft';
         cfg.foilim = [5,100];
         cfg.timwin = [-.15, .15]; % adjust if necessary, this is for the spike triggered LFP average again in seconds and this determines how long your spectrum is.
         cfg.taper = 'hanning';
         cfg.spikechannel = 'spikes';
         cfg.channel = 'lfp';
         stsFFT           = ft_spiketriggeredspectrum(cfg, data);
         ang = angle(stsFFT.fourierspctrm{1});
         mag = abs(stsFFT.fourierspctrm{1});
         ftspect = squeeze(nanmean(mag(:,1,:)));
         ftphases = squeeze(ang);
         ftfax = stsFFT.freq;

         cfg               = [];
         cfg.method        = 'ral'; % compute the rayleigh test
         cfg.spikechannel  = stsFFT.label{1};
         cfg.channel       = stsFFT.lfplabel; % selected LFP channels
         cfg.avgoverchan   = 'unweighted'; % weight spike-LFP phases irrespective of LFP power
         cfg.timwin        = 'all'; % compute over all available spikes in the window
         cfg.latency       = [0.5 1.5]; % sustained visual stimulation period
         statSts           = ft_spiketriggeredspectrum_stat(cfg,stsFFT); % gets the Rayleigh p value
         ftralp = statSts.ral;
         cfg.method = 'ppc0';
         statSts           = ft_spiketriggeredspectrum_stat(cfg,stsFFT); % gets the ppc spectrum
         ftppc = statSts.ppc0;
         cfg.method = 'plv';
         statSts           = ft_spiketriggeredspectrum_stat(cfg,stsFFT); % gets the phase locking value
         ftplv = statSts.plv;
     end
end