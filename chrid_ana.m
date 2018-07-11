function chrid_ana

animalid = '170414';
block = 3;
lcol = 'b'; %lasercolor

onlymod = 0;
printyn = 1;
sfc = 0;

supath = ['C:\Users\Julia\work\data\' animalid '\singleunits\'];
basename = [animalid '_block' int2str(block) '_tet'];

files = dir([supath, basename, '*.mat']);

cll = 1;
for fi = 1:length(files)
%     
%     if strfind(files(fi).name, 'MU')
%         continue;
%     end
        
    load([supath, files(fi).name]);    
    
    prestim = 0;
    poststim = 0;    
    respwin = 1500:2250;
    respwin = respwin-prestim;
    blwin = 1:1000;
    
%     disp(['now analyzing file: ' files(cll).name]);
    cllname{cll} = files(fi).name;
    
    i = strfind(files(fi).name, 'tet');
    tetno(cll) = strread(files(fi).name(i+3));
    
    wvchan = find(var(result.waveforms) == max(var(result.waveforms)));
    
    sr = 1000;
    lfp = result.lfp(:,wvchan)';
    
    msStimes = round(result.spikes);
    if isempty(msStimes), msStimes(1) = 0; end
    if msStimes(1) == 0, msStimes(1) = 1; end  
    
    chan = zeros(1,length(result.lfp));
    chan(msStimes) = 1;
    
    trialdur = 3000;
    msstamps = result.msstamps;
    
%     if strcmp(animalid, '170328')
%         result.randconds = result.randconds(:,1:396);
%         msstamps = msstamps(1:396);
%     end
    if length(msstamps)~=size(result.lightStamp,2)
        disp('');
        %         msstamps([518]) = []; % for 160726 block 6
        %         result.msstamps = msstamps;
        %         save([supath, files(fi).name],'result');
        pause;
    end
    
    
    for i = 1:length(msstamps)
        resp(i,:) = chan(msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim); 
        lfpresp(i,:) = result.lfp(msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim, wvchan);
        [lfpspect(i,:),trialfax] = pmtm(lfpresp(i,respwin(1)+200:respwin(end)),3,[],sr);
    end
    
    
    spike = result.waveforms(:,wvchan);
    interpspike = spline(1:32,spike,1:.1:32);
    [adiff(cll),swidth(cll)] = spikequant(interpspike);
       
    msta = linspace(-blwin(end),trialdur-blwin(end),size(resp,2));
 
       
    figure
    rasterplot(resp,msta);
    title(cllname{cll})
end


disp('')