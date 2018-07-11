function VEPAnaRD1_ll

animalid = '180619';
block = 3;
basepath = ['C:\Users\Julia\work\data\' animalid '\'];
penangle = 25;

supath = ['C:\Users\Julia\work\data\' animalid '\singleunits\'];
% supath = ['C:\Users\Julia\work\data\' animalid '\multiunits\'];
basename = [animalid '_block' int2str(block) '_tet'];

files = dir([supath, basename, '*.mat']);

recalculate = 0;

cll = 1;
for fi = 1:length(files)
%     
%     if strfind(files(fi).name, 'MU')
%         continue;
%     end
        
    load([supath, files(fi).name]);    

    prestim = 3000;
    poststim = 3000;
    trialdur = 5000;
    respwin = 1:5000; %501:1500; % after stimulus onset
    respwin = respwin+prestim;
    
    cllname{cll} = files(fi).name;
    
    i = strfind(files(fi).name, 'tet');
    tetno(cll) = strread(files(fi).name(i+3));
    
    wvchan = find(var(result.waveforms) == max(var(result.waveforms)));
    
    depth(cll) = result.depth;
    
    sr = 1000;
    lfp = result.lfp(:,wvchan)';
    
%     wo = 6.5/(1000/2); bw = wo/1.5;
%     [b,a] = iirnotch(wo,bw);
%     filtlfp = filtfilt(b,a,lfp);
    
    msStimes = round(result.spikes);
    if isempty(msStimes), msStimes(1) = 0; end
   if msStimes(1) == 0, msStimes(1) = 1; end  
    
    chan = zeros(1,length(result.lfp));
    chan(msStimes) = 1;
    
    msstamps = result.msstamps;

    binwidth = 100;
    for i = 1:length(msstamps)
        lfpresp(i,:) = lfp(msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim);
%         filtlfpresp(i,:) = filtlfp(msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim);
        resp(i,:) = chan(msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim);
        [binresp(i,:),bta] = binit(squeeze(resp(i,:)),binwidth);
        speed(i,:) = result.runspeed(msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim); 
    end
    ta = -prestim+1:trialdur+poststim;
    binresp = binresp.*(1000/binwidth);
    bta = bta-prestim;
    
    figure
    subplot(1,2,1)
    plot(ta,mean(lfpresp))
    subplot(1,2,2)
    plot(bta,mean(binresp))
    title(cllname{cll})
    
    stimind = result.order';
    stimind = stimind(:);
    for i = 1:length(unique(stimind))
        trials = find(stimind == i);
        condresp(i,:) = mean(resp(trials,:));
        bincondresp(i,:) = mean(binresp(trials,:));
        condlfpresp(i,:) = mean(lfpresp(trials,:));
    end
    
    figure
    plot(bta,bincondresp(1,:),'color',[.8,.8,.8],'linewidth',2)
    hold on
    plot(bta,bincondresp(2,:),'color',[.6,.6,.6],'linewidth',2)
    plot(bta,bincondresp(3,:),'color',[.4,.4,.4],'linewidth',2)
    plot(bta,bincondresp(4,:),'color',[0,0,0],'linewidth',2)
    ax = axis;
    axis([-500,4000,ax(3),ax(4)])
    xlabel('time (ms)')
    ylabel('firing rate (Hz)')
    title(cllname{cll})
    legend('level 1','level 2','level 3','level 4')

    cll = cll+1;
    
    disp('');
    
end

figure
plot(lfp);
hold on
for i = 1:length(msStimes)
    line([msStimes(i),msStimes(i)],[0,100],'color','r')
end
for i = 1:length(msstamps)
    line([msstamps(i),msstamps(i)],[300,400],'color','g')
end
