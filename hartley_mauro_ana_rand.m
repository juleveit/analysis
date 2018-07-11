function hartley_mauro_ana_rand   % use this for the 180 orientations

animalid = '160125';       % the ones that work are 160125 - block 12, 160128 - block 13 and 160129 - block 5
block = 12;
lcol = 'r';

% adjust path to wherever the data is
supath = ['C:\Users\Julia\work\data\' animalid '\singleunits\'];
basename = [animalid '_block' int2str(block) '_tet'];
printpath = ['C:\Users\Julia\work\data\populations\mauro\mauro_rand\' animalid '\']; 
printyn = 0; % put to 1 if you want to generate pdf output. I f yes, adjust printpath, too

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
    shownori = result.gratingInfo.orientation;            % get the orientation for each trial
    shownori(shownori == 0) = 180;                        % 0 == 180 so make it consitent
    oris = unique(shownori);      % all orientations that were shown
    for i = 1:length(shownori)      % for each of them get the angular difference to previous trial
        if i == 1
            difftoprev(i) = NaN;
        else
            difftoprev(i) = oridiff(shownori(i-1),shownori(i));   % difftoprev now has the difference to previous shown trial
        end
    end
    difftoprev(difftoprev == -90) = 90;   % 90 == -90 so be consistent
    diffs = unique(difftoprev); diffs(isnan(diffs)) = [];  % all occuring angular differences
    
    smoothwinpm = 10;   % this is how broad in degrees your smoothing kernel for the normal curve will be (shown ori == 50 then with smoothwinpm 10 you would average all trials showing 40-60 degrees)
    distwinprev = 20;   % this is the tolerance for the angular differences to previous trial
    for i = 1:180  % for each orientation
        difftothis = abs(oridiff(i,shownori));  % calculate the vector of angular differences to that angle
        thiswin = find(difftothis<smoothwinpm); % find all trials where angle was sufficiently similar
        curve(i) = mean(fr(thiswin));   % and then average those together to get that point in the normal tuning curve 'curve'
        curveerr(i) = std(fr(thiswin))./sqrt(length(thiswin)); % also get the error bars on that
        
        thisdistwinmin = find(difftothis<smoothwinpm & difftoprev<0 & difftoprev>-distwinprev);  % this is for the minus curve: find all trials where orientation on this trial was sufficiently close (as before) but also preceded by a trial with an angle smaller than this within the allowed 'distwinprev'
        mincurve(i) = mean(fr(thisdistwinmin)); % then average those
        mincurveerr(i) = std(fr(thisdistwinmin))./sqrt(length(thisdistwinmin)); % and get the error bars
        thisdistwinplus = find(difftothis<smoothwinpm & difftoprev>0 & difftoprev<distwinprev); % same for the larger angle preceded trials
        pluscurve(i) = mean(fr(thisdistwinplus));
        pluscurveerr(i) = std(fr(thisdistwinplus))./sqrt(length(thisdistwinplus));
    end
    
    % plot the tuning curve and the minus and plus curves for trials
    % preceded by smaller or larger angles
    figure
    errorbar(1:180,curve,curveerr);
    hold on
    errorbar(1:180,mincurve,mincurveerr,'r')
    errorbar(1:180,pluscurve,pluscurveerr,'g')
    errorbar(1:180,curve,curveerr);
    legend('original','previous <20 less', 'previous <20 more')
    title(printname)
    
    if printyn  % this will make a pdf if you put printyn to 1
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