function MUAana

animalid = '180522';
block = 1;
basepath = ['C:\Users\Julia\work\data\' animalid '\'];
channels = 1:16;
recalculate = 0;
recdepth = 600;
spacing = 25;
penangle = 65;

prestim = 300;
poststim = 300;
respwin = 1:500; %501:1500; % after stimulus onset
respwin = respwin+prestim;

if ~exist([basepath 'muaresult_' int2str(block) '.mat']) || recalculate
    result = MUAdataprepare(basepath,animalid,block,channels);
    save([basepath 'muaresult_' int2str(block) '.mat'],'result')
else
    load([basepath 'muaresult_' int2str(block) '.mat']);
end

for ch = channels
    
    depth(ch) = recdepth-(ch-1)*spacing;
    
    msStimes = round(result.msStimes{ch});
    if ~isempty(msStimes) & msStimes(1) == 0, msStimes(1) = 1; end
    
    chan = zeros(1,size(result.lfp,2));
    chan(msStimes) = 1;
    
    trialdur = result.stimduration*1000;
    msstamps = result.msstamps;
    if length(msstamps)~=length(result.light)
%         msstamps(385) = []; % for 140703 block 5
%         result.msstamps = msstamps;
%         save([supath, files(fi).name],'result');
        pause;
    end
    
    for i = 1:length(msstamps)
        resp(i,:) = chan(msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim); 
        lfpresp(i,:) = result.lfp(ch,msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim);
        muacresp(i,:) =  result.muac(ch,msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim);
        speed(i,:) = result.runspeed(msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim);     
    end
    
    muacresp = muacresp.^2;
    
    frs = sum(resp(:,respwin),2)./(length(respwin)/1000);
    bl = sum(resp(:,1:prestim),2)./(prestim/1000);
    muacmean = mean(muacresp(:,respwin),2);
    
    %determine if cell is visually modulated
    blfr = sum(resp(:,1:prestim),2);
    
    binwidth = 50;
    oris = unique(result.gratingInfo.Orientation); oris(find(oris == 180)) = [];
    lightlevels = unique(result.light);
    for l = 1:length(lightlevels)
        for ori = 1:length(oris)            
            thisinds = find((result.gratingInfo.Orientation == oris(ori) | result.gratingInfo.Orientation ==oris(ori)+180) &...
                result.light == lightlevels(l));
            condn(l,ori) = length(thisinds);
            condresp(l,ori,:) = nanmean(resp(thisinds,:),1);
            condmuacresp(l,ori,:) = nanmean(muacresp(thisinds,:));
            condresperr(l,ori,:) = nanstd(resp(thisinds,:),1,1)./sqrt(length(thisinds));
            if ~isnan(condresp(l,ori,:))
                [bincondresp(l,ori,:),bta] = binit(condresp(l,ori,:),binwidth);
            else
                bincondresp(l,ori,:) = binit(condresp(l,ori,:),binwidth);
            end
            binconderr(l,ori,:) = binit(condresperr(l,ori,:),binwidth);
            cfr(l,ori) = nanmean(frs(thisinds));
            cmua(l,ori) = nanmean(muacmean(thisinds));
            cerr(l,ori) =nanstd(frs(thisinds))./sqrt(length(thisinds));
            
        end
    end
    bincondresp = bincondresp.*(1000/binwidth);
    binconderr = binconderr.*(1000/binwidth);
    bta = bta-prestim;
   
    condfr(ch,:,:) = cfr; conderr(ch,:,:) = cerr;
    cresp(ch,:,:,:) = condresp;
    cmuac(ch,:,:,:) = condmuacresp;
    conmua(ch,:,:) = cmua;
    
end

depth = depth.*(cosd(90-penangle)*cosd(22));

figure
imagesc(0:5,depth,condfr(:,:,1))

figure
plot(condfr(:,6,1)-condfr(:,1,1),depth,'ko','markerfacecolor','k')
axis ij
ax = axis;
line([ax(1),ax(2)],[800,800],'color','k','linestyle',':');
line([ax(1),ax(2)],[500,500],'color','k','linestyle',':');
line([ax(1),ax(2)],[375,375],'color','k','linestyle',':');
line([0,0],[0,1000],'color','k','linewidth',2);

 disp('');