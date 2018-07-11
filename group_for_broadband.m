function group_for_broadband

% generates the files for MClust:
% this script groups several blocks of data collected at the same cortical
% location by concatenating the spiketimes so they can be sorted together
% adjust animalid, which blocks can be grouped and the depth of the block.
% when done run MClust and sort the files. Be sure to write .t WV and CQ
% files, then run generate_SU_files and you're done

animalid = '150902';
block = 3;
channel = 3;
origpath = ['C:\Users\Julia\work\data\' animalid '\'];
recfiles = dir([origpath, '*.rec']);
contactspacing = 25;
previouselectrodechannels = 0;

for i = 1:length(recfiles)
    thisname = recfiles(i).name;
    nu(i) = strread(thisname); % gets the blocknumber
end
recfile = recfiles(nu==block);

d = dir(origpath);
for j = 1:size(d)
    if(strfind(d(j).name,[animalid '_block' int2str(block)]))
        fname = d(j).name;
        break;
    end
end
load([origpath, fname]);

filedepth = 400;

if ~exist([origpath, 'broadband/'],'dir'),mkdir([origpath, 'broadband/']),end

analogdata = readTrodesFileDigitalChannels([origpath recfile.name]);
% get the timestamps
stimons = get_timestamps(analogdata.channelData(1).data);
result.msstamps = round(stimons/30);

%get wheelspeed
[speed, ta] = get_wheelspeed(analogdata.channelData(3).data);

x_t = cumsum(speed);
[smooth_win, FWHM] = MakeGaussWindow(round(1000),23.5/2, 1000);
sw_len = length(smooth_win);

x_t(end+1:end+sw_len) = x_t(end); %pad with last value for length of smoothing kernel
d_smooth_win = [0;diff(smooth_win)]/(1/1000);
dx_dt = conv(x_t,d_smooth_win,'same');
dx_dt(end-sw_len+1:end) = []; %remove values produced by convolving kernel with padded values
mouserad = 6; %cm
scalefact = 2*pi*mouserad/360;
result.runspeed = dx_dt.*scalefact;

chans = readTrodesFileChannels([origpath recfile.name],channel);

lfp = eegfilt(chans.channelData',30000,10,100)';
result.lfp = resample(lfp,1,30);

result.filedepth = filedepth;


result.depth = result.filedepth - (channel-1)*contactspacing + previouselectrodechannels*contactspacing;

%save result for single unit
resfile = [animalid '_block' int2str(block) '_channel' int2str(channel) '_broadband.mat'];
save([origpath, 'broadband/' resfile],'result');

end
