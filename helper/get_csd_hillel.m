function csd = get_csd_hillel(smpxtrialsxchans)

load('C:\Users\Julia\work\data\others\hillels\032511b_data.mat')

for i = 1:25
    for j = 1:16
        newmat(:,i,j) = eegfilt(outputMat(:,i,j)',20000,0,200)';
        resampmat(:,i,j) = resample(newmat(:,i,j),1,20);
        % extract gamma phase (here between 30 and 80 but adjust as needed)
        gamma(:,i,j) = eegfilt(resampmat(:,i,j),1000,30,80);
        h1 = hilbert(gamma(:,i,j)); gpow(:,i,j) = abs(h1); gphas(:,i,j) = angle(h1);
    end
end
smpxtrialsxchans = resampmat;

for i = 1:size(smpxtrialsxchans,2)    
    lfpresp(i,:,:) = smpxtrialsxchans(:,i,:);
end

channels = 16;
for t = 1:size(lfpresp,1)
%     mlfpresp = squeeze(mean(lfpresp(result.randconds(1,:) == 3 & result.randconds(2,:) == 1,:,:)));
    mlfpresp = squeeze(lfpresp(t,:,:));
    % mlfpresp = squeeze(mean(lfpresp(result.randconds(2,:) == 1,:,:)));


    % %spatial smoothing either over one neighbor channel on each side
    % %first duplicate uppermost and lowermost channel
    % nmlfpresp(:,1) = mlfpresp(:,1);
    % nmlfpresp(:,2:channels+1) = mlfpresp;
    % nmlfpresp(:,channels+2) = mlfpresp(:,channels);
    % %then take weighted average of channel and neighboring channels
    % for i = 1:channels
    %     smlfpresp(:,i) = (nmlfpresp(:,i)+nmlfpresp(:,i+2)+(2.*nmlfpresp(:,i+1)))./4;
    % end

    % % or smooth over 2 adjacant
    % nmlfpresp(:,1) = mlfpresp(:,1); nmlfpresp(:,2) = mlfpresp(:,1);
    % nmlfpresp(:,3:channels+2) = mlfpresp;
    % nmlfpresp(:,channels+3) = mlfpresp(:,channels); nmlfpresp(:,channels+4) = mlfpresp(:,channels);
    % for i = 1:channels
    %     smlfpresp(:,i) = (nmlfpresp(:,i)+nmlfpresp(:,i+4)+(2.*nmlfpresp(:,i+1))+2.*nmlfpresp(:,i+3)...
    %         +4.*nmlfpresp(:,i+2))./10;
    % end

    %or smooth over 3 adjacant
    nmlfpresp(:,1) = mlfpresp(:,1); nmlfpresp(:,2) = mlfpresp(:,1); nmlfpresp(:,3) = mlfpresp(:,3);
    nmlfpresp(:,4:channels+3) = mlfpresp;
    nmlfpresp(:,channels+4) = mlfpresp(:,channels); nmlfpresp(:,channels+5) = mlfpresp(:,channels); nmlfpresp(:,channels+6) = mlfpresp(:,channels);
    for i = 1:channels
        smlfpresp(:,i) = (nmlfpresp(:,i)+nmlfpresp(:,i+6)+(2.*nmlfpresp(:,i+1))+2.*nmlfpresp(:,i+5)...
            +(3.*nmlfpresp(:,i+2))+3.*nmlfpresp(:,i+4)+4.*nmlfpresp(:,i+3))./16;
    end

    %compute actual csd
    for i = 1:channels-2
        csd(t,:,i) = (smlfpresp(:,i)-(2.*smlfpresp(:,i+1))+smlfpresp(:,i+2));
    end
    
end

figure
imagesc(squeeze(nanmean(csd,1))')
colorbar
