function get_csd

animalid = '140618';

filepath = ['C:\Users\Julia\work\data\' animalid '\singleunits\'];
block = 2;
filedepth = [950];
channels = 32;

recfiles = dir([filepath, animalid '_block' int2str(block) '*.mat']);

for tet = 1:channels/4
    filebeg = [animalid '_block' int2str(block) '_tet' int2str(tet)];
    for i = 1:length(recfiles)
        if strfind(recfiles(i).name, filebeg)
            file = recfiles(i).name;
            break;
        else
            continue;
        end
    end
    load([filepath, file]);
    lfp(:,(tet-1)*4+1:tet*4) = result.lfp;
end

for i = 1:length(result.msstamps)
    lfpresp(i,:,:) = lfp(result.msstamps(i):result.msstamps(i)+3000,:);
end
mlfpresp = squeeze(mean(lfpresp));


for i = 1:channels-2
    meancsd(:,i) = mlfpresp(:,i)+mlfpresp(:,i+2)-(2.*mlfpresp(:,i+1));
end
     
for i = 1:channels-2
    csd(:,i) = lfp(:,i)+lfp(:,i+2)-(2.*lfp(:,i+1));
end

for i = 1:length(result.msstamps)
    mcsd(i,:,:) = csd(result.msstamps(i):result.msstamps(i)+3000,:);
end

save([filepath, 'sorting\', filen, 'chans_', int2str((tet-1)*4+1) '-' int2str((tet-1)*4+4) '_waveforms'],'data');
