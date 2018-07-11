function csd = get_csd_motor(animalid, block, channels)

if nargin == 0
    animalid = '170418';
    block = 3;
    channels = 16;
    allchannels = [16];
    electrode = 1;
    depth = 780; % of deepest electrode
    penangle = 30;
end

ellength = 25*(channels-1);
highest = depth-ellength;
depthax = depth:-25:highest;
depthax = depthax.*(cosd(45-penangle));
depth = depth*(cosd(45-penangle));

printyn = 0;
singleunit = 0;
recalculate = 0;

basepath = ['C:\Users\Julia\work\data\' animalid '\'];
if printyn
    if ~exist([basepath, 'pdfs/'],'dir'),mkdir([basepath, 'pdfs/']),end
    if ~exist([basepath, 'pdfs/' 'CSD/'],'dir'),mkdir([basepath, 'pdfs/' 'CSD/']),end
end
printpath = [basepath 'pdfs\CSD\'];

if singleunit
    filepath = ['C:\Users\Julia\work\data\' animalid '\singleunits\'];
    % filepath = ['C:\Users\Julia\work\data\' animalid '\multiunits\'];

    recfiles = dir([filepath, animalid '_block' int2str(block) '*.mat']);

    for tet = 1:channels/4
        filebeg = [animalid '_block' int2str(block) '_tet' int2str(tet)];
        clear file;
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
else
    pcs = 0;
    if electrode>1        
        for i = 1:electrode-1
            pcs = allchannels(i)+pcs; % previous channels
        end
    end
    evalchans = pcs+1:pcs+allchannels(electrode);
    if ~exist([basepath 'csdresult_' int2str(block) '.mat']) || recalculate
        result = lfpdataprepare(basepath,animalid,block,evalchans);   
        save([basepath 'csdresult_' int2str(block) '.mat'], 'result')
    else
        load([basepath 'csdresult_' int2str(block) '.mat']);
    end
    lfp = result.lfp';
end

light = result.randconds(2,:); position = result.randconds(1,:);
result.stimduration = result.sweeplength;
if strcmp(animalid, '170328')
    result.randconds = result.randconds(:,1:396);
end
   
for i = 1:length(result.msstamps)
    lfpresp(i,:,:) = lfp(result.msstamps(i)+1:result.msstamps(i)+result.stimduration,:);
end
timeax = 1:result.stimduration;

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
% %or over two sites away
% for i = 1:channels-4
%     csd(:,i) = smlfpresp(:,i)+smlfpresp(:,i+4)-(2.*smlfpresp(:,i+2));
% end

cond =  result.randconds(2,:) == 0; %result.randconds(1,:) == 5 &

% get the depths right:
highestcsd = highest + 25; %because one gets lost for CSD calculation
lowestcsd = depth - 25;
depthaxcsd = lowestcsd:-25:highestcsd;

figure
imagesc(timeax,depthaxcsd,squeeze(nanmean(csd(cond,:,:),1))')
colorbar
% axis([0,200,depthaxcsd(end),depthaxcsd(1)])
if printyn
    figSize = [30 21];
    set(gcf,'paperunits','centimeters','papersize',figSize,'paperposition',[0 0 figSize])
    print(gcf,[printpath ,  '06_CSD.pdf'], '-dpdf' );
end

disp('');