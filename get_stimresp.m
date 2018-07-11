function [stimresp, stimlfpspect, stimfax] = get_stimresp(filename)

load(filename);

wvchan = find(var(result.waveforms) == max(var(result.waveforms)));
lfp = result.lfp(:,wvchan)';
sr = 1000;

msStimes = round(result.spikes);
if isempty(msStimes), msStimes(1) = 0; end
if msStimes(1) == 0, msStimes(1) = 1; end

chan = zeros(1,length(result.lfp));
chan(msStimes) = 1;

trialdur = 3000;
msstamps = result.msstamps;

if length(msstamps)~=size(result.lightStamp,2)
    disp('');
    %         msstamps([518]) = []; % for 160726 block 6
    %         result.msstamps = msstamps;
    %         save([supath, files(fi).name],'result');
    pause;
end

for i = 1:length(msstamps)
    stimresp(i,:) = chan(msstamps(i):msstamps(i) + trialdur);
    lfpresp(i,:) = result.lfp(msstamps(i):msstamps(i) + trialdur, wvchan);
    [stimlfpspect(i,:),stimfax] = pmtm(lfpresp(i,:),3,[],sr);
end


