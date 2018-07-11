function hartley_mauro_ana

animalid = '160122';  % the ones that work are 160122 - block 5 and 160129 - block 4
block = 5;
lcol = 'r';

supath = ['C:\Users\Julia\work\data\' animalid '\singleunits\'];  % adjust the path to point to the data
basename = [animalid '_block' int2str(block) '_tet'];
printpath = ['C:\Users\Julia\work\data\populations\mauro\mauro_disc\' animalid '\'];
printyn = 0;  % put to 1 if you want to print pdfs for each tuning curve

files = dir([supath, basename, '*.mat']);

% response window - over what time do you want to count spikes. You can adjust this here
respwin = 41:440; % after stimulus onset 
trialdur = 2499; % this is stimulus plus mask plus isi
prestim = 300;  % baseline period before stimulus on
respwin = respwin + prestim; 

for cell = 1:length(files)
    
   
    load([supath, files(cell).name]);  % this gives you the structure 'result' which has all the information on the stimulus as well as the neural data
%     disp(['now analyzing file: ' files(cell).name]);
    cellname{cell} = files(cell).name;
    printname = files(cell).name;
    printname(find(printname=='_')) = ' ';
    
    i = strfind(files(cell).name, 'tet');
    tetno = strread(files(cell).name(i+3));
    
    wvchan = find(var(result.waveforms) == max(var(result.waveforms))); % because we sort tetrodes, this gives you the channel on which the spike was biggest to determine cortical depth - don't worry about it
    spike = result.waveforms(:,wvchan);  % stores the average waveform of the spike of that cell
    interpspike = spline(1:32,spike,1:.1:32);  % interpolates it
%     [adiff(cell),swidth(cell)] = spikequant(interpspike);  % gets some waveform features to determine whether it is an internueron or a pyramidal cell - you don't need this I guess
    
    depth(cell) = result.depth;

    msStimes = round(result.spikes);  % round all the spike-times to the nearest millisecond
    if ~isempty(msStimes) & msStimes(1) == 0, msStimes(1) = 1; end
    
    chan = zeros(1,length(result.lfp));   % make a long vector of zeros
    chan(msStimes) = 1;                   % and put ones in everyt time a spikes happened
    
    msstamps = result.msstamps;  % these are the trial onset timestamps
    for i = 1:length(msstamps)  % for each timestamp, cut the appropriate snippet of 
        resp(i,:) = chan(msstamps(i)-prestim:msstamps(i) + trialdur); % spiking vector
        lfpresp(i,:) = result.lfp(msstamps(i)-prestim:msstamps(i) + trialdur);  % and LFP-vector - you don't need this
    end
    ta = -prestim:trialdur;  % time axis taking baseline period into account
    
    fr = sum(resp(:,respwin),2).*(1000/length(respwin));  % get the firing rate in Hz for each trial
    oris = unique(result.gratingInfo.orientation);        % get all the orientations that were presented
    for i = 1:length(result.gratingInfo.orientation)      % and for each trial find the angular difference to preceding trial
        if i == 1
            difftoprev(i) = NaN;
        else
            difftoprev(i) = oridiff(result.gratingInfo.orientation(i-1),result.gratingInfo.orientation(i));
        end
    end
    difftoprev(difftoprev == -90) = 90;  % 90 == -90 so keep it consitant
    diffs = unique(difftoprev); diffs(isnan(diffs)) = [];  % all occuring angular differences between trials
    for i = 1:length(oris)   % for each orienation
        norm(i) = mean(fr(result.gratingInfo.orientation == oris(i)));  % get the normal tuning curve point by just averaging all trials with that orientation
        normerr(i) = std(fr(result.gratingInfo.orientation == oris(i)))./sqrt(length(find(result.gratingInfo.orientation == oris(i))));  % and get error bars on that
        for j = 1:length(diffs)  % and then for each possible preceding trial difference
            oneback(i,j) = mean(fr(result.gratingInfo.orientation == oris(i) & difftoprev == diffs(j))); % fill up a matrix (i,j) where i is ori of this trial and j is difference to preceding trial
            onebackerr(i,j) = std(fr(result.gratingInfo.orientation == oris(i) & difftoprev == diffs(j)))./sqrt(length(find(result.gratingInfo.orientation == oris(i) & difftoprev == diffs(j))));  % and the error bars on that
        end
    end
    
    % plot the figure
    figure
    errorbar(oris,norm,normerr)
    hold on
    errorbar(oris,oneback(:,4),onebackerr(:,4),'r')
    errorbar(oris,oneback(:,6),onebackerr(:,6),'g')
    legend('original','previous 18 less', 'previous 18 more')
    set(gca,'xtick',oris)
    title(printname)
    
    if printyn   % and save as a pdf if printyn is 1 - if so adjust path
        figSize = [30 21];
        set(gcf,'paperunits','centimeters','papersize',figSize,'paperposition',[0 0 figSize])
        if cell<10, printi = ['0', int2str(cell)]; else printi = int2str(cell); end
        print([printpath ,  printi '__' printname '.pdf'],'-dpdf')
    end
    
end
end

% directed angular difference for 180 deg
function a = oridiff(x,y)
    a = x-y;
    a(find(a>90)) = -abs(a(find(a>90))-180);
    a(find(a<-90)) = abs(abs(a(find(a<-90)))-180);
end