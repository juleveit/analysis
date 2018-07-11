function position_populationana

% X98 ChR2
animalids = {'170404','170501','170501','170501','170502','170502'};
blocks    = [3,        1,       3,       5,       1,       3];
animal    = [1,        2,       2,       2,       3,       3];
electrodes =[[1,32];  [1,16];  [1,16];  [1,16];  [1,16];  [1,16]];
penangle =  [30,       30,      30,      30,      30,      30];
popfile = 'somefile.mat'; % adjust to where you want it to go

% % X94 ChR2
% animalids = {'170322', '170322','170328','170330','170414','170418','170418'};
% blocks    = [1,         2,       1,       1,       2,       1,       3];
% animal    = [1,         1,       2,       3,       4,       5,       5];
% electrodes =[[1,16];   [1,32];  [1,32];  [1,32];  [1,32];  [1,16];  [1,16]];
% penangle =  [25,        25,      30,      30,      30,      30,      30];

recalculate = 0;

sr = 1000;
prestim = 0;
poststim = 0;
respwin = 1500:2250;
respwin = respwin+prestim;
blwin = 1:1000;

if ~exist(popfile) || recalculate

    cll = 1;
    for blck = 1:length(blocks)

        supath = ['C:\Users\Julia\work\data\' animalids{blck} '\singleunits\']; %adjust
        basename = [animalids{blck} '_block' int2str(blocks(blck)) '_tet'];
        files = dir([supath, basename, '*.mat']);
        
        for fi = 1:length(files)

            if strfind(files(fi).name, 'MU')
                mu(cll) = 1; su(cll) = 0;
            else
                mu(cll) = 0; su(cll) = 1;
            end
            
            load([supath, files(fi).name]);            
            
            i = strfind(files(fi).name, 'tet');
            if strcmp(files(fi).name(i+4),'_')
                tetno = strread(files(fi).name(i+3)); % single character number
            else
                tetno = strread(files(fi).name(i+3:i+4)); % number >10
            end
            
            cllname{cll} = files(fi).name;
            printname = files(fi).name;
            printname(find(printname=='_')) = ' ';
                        
            % important stuff
            depth(cll) = result.depth;
            pangle(cll) = penangle(blck);
            recording(cll) = blck;
            animalno(cll) = animal(blck);
            
            light = result.randconds(2,:); position = result.randconds(1,:);
            
            wvchan = find(var(result.waveforms) == max(var(result.waveforms)));
            lfp = result.lfp(:,wvchan)';
            spike = result.waveforms(:,wvchan);
            interpspike = spline(1:32,spike,1:.1:32);
            [adiff(cll),swidth(cll)] = spikequant(interpspike);            
            secpersamp = 1/30000;
            interpf = secpersamp/10;
            swidthms(cll) = swidth(cll)*interpf*1000;
            waveform(cll,:) = spike;
            clustqual(cll) = result.clusterquality;

            % get spiketimes
            msStimes = round(result.spikes);
            chan = zeros(1,length(result.lfp));
            chan(msStimes) = 1;    
            
            trialdur = result.sweeplength;
            msstamps = result.msstamps;
            
            % running
            for i = 1:length(msstamps)
                speed(i,:) = result.runspeed(msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim);
            end
       
            % figure out sufficiently high and nonvariable runspeed trials
            meanspeed = mean(speed(:,respwin),2);
            stdspeed = std(speed(:,respwin),1,2);
            notstill = find(meanspeed>1);
            okspeed = find(meanspeed>( mean(meanspeed(notstill))-(1.5*std(meanspeed(notstill))) ) & meanspeed>1);
            okvar = find(stdspeed<( mean(stdspeed(notstill))+(1.5*std(stdspeed(notstill)))) & stdspeed>.5);
            oktrials = intersect(okspeed,okvar);
            nonoktrials = 1:size(speed,1); nonoktrials(oktrials) = [];
            stilltrials = 1:size(speed,1); stilltrials(notstill) = [];
            
            nrun(cll) = length(oktrials);
            nstill(cll) = length(stilltrials);
            if ~isempty(oktrials)
                cellmeanspeed(cll) = mean(meanspeed(oktrials));
                cellminspeed(cll) = min(meanspeed(oktrials));
            else
                cellmeanspeed(cll) = NaN;
                cellminspeed(cll) = NaN;
            end
            
            for i = 1:length(msstamps)
                resp(i,:) = chan(msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim);
                lfpresp(i,:) = lfp(msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim);               
                msl(i) = mean(find(resp(i,respwin))); % mean latency
            end                        
            msta = linspace(-blwin(end),trialdur-blwin(end),size(resp,2));
            
            frs = sum(resp(:,respwin),2)./(length(respwin)/1000);
            bl = sum(resp(:,blwin),2)./(length(blwin)/1000);
            sc = sum(resp(:,respwin),2); % spike count            
            
            %determine if cll is touch modulated
            touchfr = sum(resp(:,blwin(end)+40:blwin(end)+40+200),2);
            shortbl = sum(resp(:,blwin(end)-40-200:blwin(end)-40),2);
            touchmod(cll) = ttest2(shortbl,touchfr);
            
            %determine if cll is modulated by light
            lightmod(cll) = ttest2(frs(find(light)),frs(find(light == 0)));
            
            lfr(cll) = mean(frs(find(light)));
            nlfr(cll) = mean(frs(find(light == 0)));
            
            binwidth = 20;
            [binnedlight,bta] = binit(mean(resp(find(light),:)),binwidth);
            
            for l = 1:length(unique(result.light))
                for pos = 1:length(result.positions)
                    trials = find(position == result.positions(pos) & light == l-1);
                    
                    condresp(cll,l,pos,:) = mean(resp(trials,:),1);
                    condresperr(cll,l,pos,:) = nanstd(resp(trials,:),1,1)./sqrt(length(trials));
                    condlfpresp(cll,l,pos,:) = mean(lfpresp(trials,:),1);
                    
                    condfr(cll,l,pos) = mean(frs(trials));%-mean(bl);
                    conderr(cll,l,pos) =std(frs(trials))./sqrt(length(trials));
                    
                    hlp = ones(size(condresp,4),1);
                    [xx,bta] = binit(hlp,binwidth);
                    bincondresp(cll,l,pos,:) = binit(condresp(cll,l,pos,:),binwidth).*(1000/binwidth);
                    binconderr(cll,l,pos,:) = binit(condresperr(cll,l,pos,:),binwidth).*(1000/binwidth);
                    
                    thisruninds = intersect(trials,oktrials);
                    if ~isempty(thisruninds)
                        thisrunn(cll,l,pos) = length(thisruninds);
                        runcondresp(cll,l,pos,:) = mean(resp(thisruninds,:),1);
                        runcondresperr(cll,l,pos,:) = nanstd(resp(thisruninds,:),1,1)./sqrt(length(thisruninds));
                        runcondfr(cll,l,pos) = mean(frs(thisruninds));
                        runconderr(cll,l,pos) = std(frs(thisruninds))./sqrt(length(thisruninds));
                        runbincondresp(cll,l,pos,:) = binit(runcondresp(cll,l,pos,:),binwidth).*(1000/binwidth);
                        runbinconderr(cll,l,pos,:) = binit(runcondresperr(cll,l,pos,:),binwidth).*(1000/binwidth);
                    else
                        thisrunn(cll,l,pos) = 0;
                        runcondresp(cll,l,pos,:) = nan(1,size(resp,2));
                        runcondresperr(cll,l,pos,:) = nan(1,size(resp,2));
                        runcondfr(cll,l,pos) = NaN;
                        runconderr(cll,l,pos) = NaN;
                        runbincondresp(cll,l,pos,:) = nan(1,length(bta));
                        runbinconderr(cll,l,pos,:) = nan(1,length(bta));
                    end
                    
                    thisstillinds = intersect(trials,stilltrials);
                    if ~isempty(thisstillinds)
                        thisstilln(cll,l,pos) = length(thisstillinds);
                        stillcondresp(cll,l,pos,:) = mean(resp(thisstillinds,:),1);
                        stillcondresperr(cll,l,pos,:) = nanstd(resp(thisstillinds,:),1,1)./sqrt(length(thisstillinds));
                        stillcondfr(cll,l,pos) = mean(frs(thisstillinds));
                        stillconderr(cll,l,pos) = std(frs(thisstillinds))./sqrt(length(thisstillinds));
                        stillbincondresp(cll,l,pos,:) = binit(stillcondresp(cll,l,pos,:),binwidth).*(1000/binwidth);
                        stillbinconderr(cll,l,pos,:) = binit(stillcondresperr(cll,l,pos,:),binwidth).*(1000/binwidth);
                    else
                        thisstilln(cll,l,pos) = 0;
                        stillcondresp(cll,l,pos,:) = nan(1,size(resp,2));
                        stillcondresperr(cll,l,pos,:) = nan(1,size(resp,2));
                        stillcondfr(cll,l,pos) = NaN;
                        stillconderr(cll,l,pos) = NaN;
                        stillbincondresp(cll,l,pos,:) = nan(1,length(bta));
                        stillbinconderr(cll,l,pos,:) = nan(1,length(bta));
                    end
                end
            end            
            ta = bta-blwin(end);
            
            %  figure
            clf;
            subplot(2,2,1)
            plot(ta,squeeze(mean(runbincondresp(cll,1,:,:),3)),'k','linewidth',2);
            hold on
            plot(ta,squeeze(mean(runbincondresp(cll,2,:,:),3)),'b','linewidth',2);
            mx = max([max(squeeze(mean(runbincondresp(cll,1,:,:),3))),max(squeeze(mean(runbincondresp(cll,2,:,:),3))),.01]);
            axis([-blwin(end),trialdur-blwin(end),0,mx]);
            line([0,0],[0,mx],'color','k','linewidth',2);
            line([1500,1500],[0,mx],'color','k','linewidth',2);
            line([500,500],[0,mx],'color','b','linewidth',2)
            line([1250,1250],[0,mx],'color','b','linewidth',2);
            legend({'Light OFF','Light ON'})
            xlabel('time [ms]')
            ylabel('firing rate [Hz]')
            title(['cll ' int2str(cll) ' depth: ' int2str(result.depth), 'cll ' cllname{cll} ])
            
            subplot(2,2,2)
            errorbar(result.positions,squeeze(runcondfr(cll,1,:)),squeeze(runconderr(cll,1,:)),'ko-','markerfacecolor','k','linewidth',2)
            hold on
            errorbar(result.positions,squeeze(runcondfr(cll,2,:)),squeeze(runconderr(cll,2,:)),'ro-','markerfacecolor','r','linewidth',2)
            ax = axis;
            axis([-1,9,ax(3),ax(4)]);
            legend('light off','light on')
            xlabel('bar position')
            ylabel('firing rate [Hz]')
            
            subplot(2,2,3)
            plot(ta,squeeze(mean(stillbincondresp(cll,1,:,:),3)),'k','linewidth',2);
            hold on
            plot(ta,squeeze(mean(stillbincondresp(cll,2,:,:),3)),'b','linewidth',2);
            mx = max([max(squeeze(mean(stillbincondresp(cll,1,:,:),3))),max(squeeze(mean(stillbincondresp(cll,2,:,:),3))),.01]);
            axis([-blwin(end),trialdur-blwin(end),0,mx]);
            line([0,0],[0,mx],'color','k','linewidth',2);
            line([1500,1500],[0,mx],'color','k','linewidth',2);
            line([500,500],[0,mx],'color','b','linewidth',2)
            line([1250,1250],[0,mx],'color','b','linewidth',2);
            legend({'Light OFF','Light ON'})
            xlabel('time [ms]')
            ylabel('firing rate [Hz]')
            title('non-running')
            
            subplot(4,4,11)
            plot(waveform(cll,:))
            axis([0,40,-100,100])
            legend(['width: ' num2str(swidthms(cll))])
            
            disp([files(fi).name '   done'])
            cll = cll + 1;
            
        end
    end
    save(popfile, '-v7.3');
else
    load(popfile);   
end

prsv = swidthms>.36; pfsv = swidthms<=.36;
prs = find(prsv); pfs = find(pfsv);

% adjust depth according to penetration angle
depth = depth.*cosd(45-pangle);

phe = zeros(1,length(depth));

% putative halo expressing for X94 ChR2
phe([10, 20, 102,115]) = 1; 

phe = logical(phe);

l23 = depth<375;
l4 = depth>=375&depth<=550;
l5 = depth>550&depth<=800;
l5a = depth>550&depth<=650;
l5b = depth>650&depth<=800;
l6 = depth>800;
l23rs = l23&prsv&~phe;
l23fs = l23&pfsv&~phe;
l4rs = l4&prsv&~phe;
l4fs = l4&pfsv&~phe;
l5rs = l5&prsv&~phe;
l5fs = l5&pfsv&~phe;
l6rs = l6&prsv&~phe;
l6fs = l6&pfsv&~phe;
okrs = prsv&~phe; 

for i = 1:length(depth)
    mfr(i) = mean(runcondfr(i,1,2:9),3);
    omi(i) = (nanmean(runcondfr(i,2,2:9),3)-nanmean(runcondfr(i,1,2:9),3))./(nanmean(runcondfr(i,2,2:9),3)+nanmean(runcondfr(i,1,2:9),3));
    ominr(i) = (nanmean(stillcondfr(i,2,2:9),3)-nanmean(stillcondfr(i,1,2:9),3))./(nanmean(stillcondfr(i,2,2:9),3)+nanmean(stillcondfr(i,1,2:9),3));
    deltaspikes(i) = nanmean(runcondfr(i,2,:),3)-nanmean(runcondfr(i,1,:),3);
    for l = 1:2
        ssi(i,l) = get_ssi(squeeze(runcondfr(i,l,2:8))); % spatial selectivity
    end
    for pos = 1:9
        omicurve(i,pos) = (runcondfr(i,2,pos)-runcondfr(i,1,pos))./(runcondfr(i,2,pos)+runcondfr(i,1,pos));
        deltacurve(i,pos) = runcondfr(i,2,pos)-runcondfr(i,1,pos);
    end
    [ranked(i,1,:),si] = sort(squeeze(runcondfr(i,1,2:9)),'descend');
    ranked(i,2,:) = squeeze(runcondfr(i,2,si));
    normranked(i,:,:) = ranked(i,:,:)./ranked(i,1,1);
    rankedomi(i,:) = (ranked(i,2,:)-ranked(i,1,:))./(ranked(i,2,:)+ranked(i,1,:));
end

%running OMI no touch - scatter
figure
plot(omicurve(prsv&touchmod,1),depth(prsv&touchmod),'ko','markerfacecolor','k')
hold on
plot(omicurve(pfsv&touchmod,1),depth(pfsv&touchmod),'go','markerfacecolor','g')
plot(omicurve(phe&touchmod,1),depth(phe&touchmod),'ro','markerfacecolor','r')
plot(omicurve(find(mu),1),depth(find(mu)),'co','markerfacecolor','c')
axis ij
line([0,0],[100,900],'color','k')
line([-1,1],[575,575],'color','k','linestyle',':')
line([-1,1],[375,375],'color','k','linestyle',':')
legend('RS','FS','MU')
xlabel('OMI')
ylabel('cortical depth')
title('OMI running - no touch')

% running OMI no touch - running average
cond = ~phe & mu == 0 & mfr>1.5;
figure
[x,y,xerr] = runningAverage(depth(cond),omicurve(cond,1),20,15);
plot(x,y,'k','linewidth',2)
hold on
plot(x+xerr,y,'k')
plot(x-xerr,y,'k')
axis ij
line([0,0],[100,900],'color','k')
line([-.8,.1],[575,575],'color','k','linestyle',':')
line([-.8,.1],[375,375],'color','k','linestyle',':')
axis([-.8,.1,200,950])
xlabel('average OMI')
ylabel('cortical depth um')

% running OMI mean touch scatter
figure
plot(omi(prsv&touchmod&lightmod),depth(prsv&touchmod&lightmod),'ko','markerfacecolor','k')
hold on
plot(omi(pfsv&touchmod&lightmod),depth(pfsv&touchmod&lightmod),'go','markerfacecolor','g')
plot(omi(phe&touchmod&lightmod),depth(phe&touchmod&lightmod),'ro','markerfacecolor','r')
plot(omi(find(mu)),depth(find(mu)),'co','markerfacecolor','c')
axis ij
line([0,0],[100,900],'color','k')
line([-1,1],[575,575],'color','k','linestyle',':')
line([-1,1],[375,375],'color','k','linestyle',':')
legend('RS','FS','MU')
xlabel('OMI')
ylabel('cortical depth')
title('OMI running')

figure
cond = ~phe & mu == 0;
[x,y,xerr] = runningAverage(depth(cond),omi(cond),20,15);
plot(x,y,'k','linewidth',2)
hold on
plot(x+xerr,y,'k')
plot(x-xerr,y,'k')
axis ij
line([0,0],[100,900],'color','k')
line([-.5,.1],[575,575],'color','k','linestyle',':')
line([-.5,.1],[375,375],'color','k','linestyle',':')
axis([-.5,.1,200,950])
xlabel('average OMI')
ylabel('cortical depth um')

cond = prsv&~phe&mu == 0;
figure
[x,y,xerr] = runningAverage(depth(cond),omi(cond),20,15);
plot(x,y,'k','linewidth',2)
hold on
plot(x+xerr,y,'k')
plot(x-xerr,y,'k')
axis ij
line([0,0],[100,900],'color','k')
line([-.5,.1],[575,575],'color','k','linestyle',':')
line([-.5,.1],[375,375],'color','k','linestyle',':')
axis([-.5,.1,200,950])
xlabel('average OMI')
ylabel('cortical depth um')

% non-running OMI in depth
figure
plot(ominr(prs),depth(prs),'ko','markerfacecolor','k')
hold on
plot(ominr(pfs),depth(pfs),'go','markerfacecolor','g')
plot(omi(phe),depth(phe),'ro','markerfacecolor','r')
axis ij
line([0,0],[100,900],'color','k')
line([-1,1],[550,550],'color','k','linestyle',':')
line([-1,1],[375,375],'color','k','linestyle',':')
legend('RS','FS','ChR2 cell')
xlabel('OMI')
ylabel('cortical depth')
title('OMI non-running')

figure
[x,y,xerr] = runningAverage(depth(~phe),ominr(~phe),20,15);
plot(x,y,'k','linewidth',2)
hold on
plot(x+xerr,y,'k')
plot(x-xerr,y,'k')
axis ij
line([0,0],[100,900],'color','k')
line([-.5,.1],[550,550],'color','k','linestyle',':')
line([-.5,.1],[375,375],'color','k','linestyle',':')
axis([-.6,.1,200,950])
xlabel('average OMI (non-runnning')
ylabel('cortical depth um')


stepsize = 50;
a = 100:stepsize:1000;
for i = 1:length(a)
    dinds = find(depth>=a(i)-stepsize/2 & depth<a(i)+stepsize/2);
    momi(i) = mean(omi(dinds)); % meanomi for bin
    momierr(i) = std(omi(dinds))./sqrt(length(dinds));
end

figure
barh(a,momi,'k')
hold on
herrorbar(momi,a,momierr,'k.')
axis ij
axis([-.5,.2,0,1050])
line([-.5,.2],[375,375],'linestyle',':','color','k')
line([-.5,.2],[575,575],'linestyle',':','color','k')
xlabel('OMI')
ylabel('depth (um)')


% spatial selectivity
figure
plot(ssi(prs,1),depth(prs),'ko','markerfacecolor','k')
hold on
plot(ssi(pfs,1),depth(pfs),'go','markerfacecolor','g')
axis ij
line([0,1],[550,550],'color','k','linestyle',':')
line([0,1],[375,375],'color','k','linestyle',':')
legend('RS','FS')
xlabel('SSI')
ylabel('cortical depth')
title('SSI running')

% spatial selectivity change
figure
plot(ssi(prs,2)-ssi(prs,1),depth(prs),'ko','markerfacecolor','k')
hold on
plot(ssi(pfs,2)-ssi(pfs,1),depth(pfs),'go','markerfacecolor','g')
axis ij
line([0,0],[100,900],'color','k')
line([-.6,.6],[550,550],'color','k','linestyle',':')
line([-.6,.6],[375,375],'color','k','linestyle',':')
axis([-.6,.6,100,900])
legend('RS','FS')
xlabel('SSI')
ylabel('cortical depth')
title('change in SSI with light')



function ssi = get_ssi(curve)
ssi  = 1 - (((norm(curve)/max(curve)) - 1)./((sqrt(length(curve)))-1));