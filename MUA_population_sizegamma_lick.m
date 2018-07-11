function MUA_population_sizegamma_lick

% %lick animals lickport task
% animalids = {'161011_2'};
% blocks    = [1];
% animal    = [1];
% electrodes =[[1,16]];
% bldepth    =[500];
% penangle =  [25];
% spacing =   [25];
% popfile = 'C:\Users\Julia\work\data\populations\tmp.mat';

% %lick animals lickport detection task
% animalids = {'161013'};
% blocks    = [2];
% animal    = [1];
% electrodes =[[1,16]];
% bldepth    =[500];
% penangle =  [25];
% spacing =   [25];
% popfile = 'C:\Users\Julia\work\data\populations\tmp.mat';

%lick animals lickport detection task
animalids = {'161019'};
blocks    = [3];
animal    = [1];
electrodes =[[1,16]];
bldepth    =[500];
penangle =  [25];
spacing =   [25];
popfile = 'C:\Users\Julia\work\data\populations\tmp.mat';


depth = bldepth.*(cosd(penangle)*cosd(22));
spacing = spacing.*(cosd(penangle)*cosd(22));
evaldepth = 300;
for i = 1:length(depth)
    for j = 1:length(electrodes(i,1):electrodes(i,2))
        dm(i,j) = depth(i)-((j-1)*spacing(i)); % depth matrix
    end
    [c,di(i)] = min(abs(dm(i,:)-evaldepth)); % get depth index, least distance to evaldepth
end

recalculate_pop = 0;
recalculate_muafile = 0;

%chronux parameters
params.tapers = [5,9]; params.Fs = 1000; params.err = [2, 0.05]; params.trialave = 1;


if ~exist(popfile) || recalculate_pop

    cll = 1;
    for blck = 1:length(blocks)
        
        basepath = strcat('C:\Users\Julia\work\data\', animalids{blck}, '\');
        file = strcat(basepath, 'muaresult_', int2str(blocks(blck)), '_', int2str(electrodes(blck,1)), '_', int2str(electrodes(blck,2)), '.mat');
        
        clear result;
        if ~exist(file) || recalculate_muafile
            result = MUAdataprepare(basepath,animalids{blck},blocks(blck),electrodes(blck,1):electrodes(blck,2));
            save(file,'result')
        else
            clear centresult; clear surrresult;
            load(file);
            if exist('centresult')
                result = centresult;
            elseif exist('surrresult')
                result = surrresult;
            end            
        end        
        
        prestim = 300;
        poststim = 1700;
        respwin = 1:1500; % after stimulus onset
        respwin = respwin+prestim;
        
        ch = di(blck);
        
        disp(['Block ' int2str(blck) '/' int2str(length(blocks)) '   channel ' int2str(ch) '/' int2str(length(electrodes(blck,1):electrodes(blck,2)))]);
              
        trialdur = result.stimduration*1000;
        msstamps = result.msstamps;
        if length(msstamps)~=length(result.light)
%             msstamps([193,291]) = []; % for 151023 block 5
%             result.msstamps = msstamps;
%             save(file,'result');
            pause;
        end
        
        if ~isempty(result.licks)
            if result.licks(1) == 0; result.licks(1) = 1; end
        end
        rewvec = zeros(1,size(result.lfp,2)); rewvec(result.rewards) = 1;   
        lickvec = zeros(1,size(result.lfp,2)); lickvec(result.licks) = 1;
        for i = 1:length(msstamps)
            speed(i,:) = result.runspeed(msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim);
            lickmat(i,:) = lickvec(msstamps(i) - prestim+1:msstamps(i) + trialdur + poststim);
            if find(rewvec(msstamps(i):msstamps(i)+trialdur))
                autotrial(i) = 1;
            else
                autotrial(i) = 0;
            end
            if find(rewvec(msstamps(i)+trialdur:msstamps(i)+trialdur+1000))
                rewarded(i) = 1;
            else
                rewarded(i) = 0;
            end
        end
        
        if find(result.light)
            a0 = find(rewarded&result.light == 0);
            a1 = find(rewarded&result.light == 1);
            b0 = find(~rewarded&~autotrial&result.light == 0);
            b1 = find(~rewarded&~autotrial&result.light == 1);
            for i = 1:length(a0) % power rewarded
                [pr0(i,:),f] = pmtm(result.lfp(ch,result.msstamps(a0(i))+300:result.msstamps(a0(i))+1300),3,[],1000);
            end
            for i = 1:length(a1) % power rewarded
                [pr1(i,:),f] = pmtm(result.lfp(ch,result.msstamps(a1(i))+300:result.msstamps(a1(i))+1300),3,[],1000);
            end
            for i = 1:length(b0) % power non-rewarded
                [pn0(i,:),f] = pmtm(result.lfp(ch,result.msstamps(b0(i))+300:result.msstamps(b0(i))+1300),3,[],1000);
            end
            for i = 1:length(b1) % power non-rewarded
                [pn1(i,:),f] = pmtm(result.lfp(ch,result.msstamps(b1(i))+300:result.msstamps(b1(i))+1300),3,[],1000);
            end            

            figure
            semilogy(f,nanmean(pr0,1))
            hold on
            semilogy(f,nanmean(pr1,1),'r')
            semilogy(f,nanmean(pn0,1),'c')
            semilogy(f,nanmean(pn1,1),'m')
%             mx = max([max(nanmean(pn(:,1:104),1)),max(nanmean(pr(:,1:104),1))]);
%             mn = min([min(nanmean(pn(:,1:104),1)),min(nanmean(pr(:,1:104),1))]);
%             axis([0,100,mn,mx]);
%             legend('reward no light','reward light','no reward no light','no reward light')

            fillx = [f(1:104)',fliplr(f(1:104)')];
            fillyr0 = [nanmean(pr0(:,1:104),1)+nanstd(pr0(:,1:104),1)./sqrt(size(pr0,1)),fliplr(nanmean(pr0(:,1:104),1)-nanstd(pr0(:,1:104),1)./sqrt(size(pr0,1)))];
            fillyr1 = [nanmean(pr1(:,1:104),1)+nanstd(pr1(:,1:104),1)./sqrt(size(pr1,1)),fliplr(nanmean(pr1(:,1:104),1)-nanstd(pr1(:,1:104),1)./sqrt(size(pr1,1)))];
            fillyn0 = [nanmean(pn0(:,1:104),1)+nanstd(pn0(:,1:104),1)./sqrt(size(pn0,1)),fliplr(nanmean(pn0(:,1:104),1)-nanstd(pn0(:,1:104),1)./sqrt(size(pn0,1)))];
            fillyn1 = [nanmean(pn1(:,1:104),1)+nanstd(pn1(:,1:104),1)./sqrt(size(pn1,1)),fliplr(nanmean(pn1(:,1:104),1)-nanstd(pn1(:,1:104),1)./sqrt(size(pn1,1)))];

            figure
            fill(fillx,fillyr0,'b')
            hold on
            fill(fillx,fillyr1,'r')
            fill(fillx,fillyn0,'c')
            fill(fillx,fillyn1,'m')
            set(gca,'yscale','log')
            legend('reward no light','reward light','no reward no light','no reward light')
        else
            
            a = find(rewarded);
            b = find(~rewarded&~autotrial);
            for i = 1:length(a) % power rewarded
                [pr(i,:),f] = pmtm(result.lfp(ch,result.msstamps(a(i))+300:result.msstamps(a(i))+1300),3,[],1000);
            end
            for i = 1:length(b) % power non-rewarded
                [pn(i,:),f] = pmtm(result.lfp(ch,result.msstamps(b(i))+300:result.msstamps(b(i))+1300),3,[],1000);
            end

            % throw out dryed out headplate trials for 161018-1
            if strcmp(animalids{blck},'161018') & blocks(blck) == 1
                pr(88:91,:) = nan(4,513);
                pn(29:30,:) = nan(2,513);
            end

            figure
            semilogy(f,nanmean(pr,1))
            hold on
            semilogy(f,nanmean(pn,1),'r')
            mx = max([max(nanmean(pn(:,1:104),1)),max(nanmean(pr(:,1:104),1))]);
            mn = min([min(nanmean(pn(:,1:104),1)),min(nanmean(pr(:,1:104),1))]);
            axis([0,100,mn,mx]);

            fillx = [f(1:104)',fliplr(f(1:104)')];
            fillyr = [nanmean(pr(:,1:104),1)+nanstd(pr(:,1:104),1)./sqrt(size(pr,1)),fliplr(nanmean(pr(:,1:104),1)-nanstd(pr(:,1:104),1)./sqrt(size(pr,1)))];
            fillyn = [nanmean(pn(:,1:104),1)+nanstd(pn(:,1:104),1)./sqrt(size(pn,1)),fliplr(nanmean(pn(:,1:104),1)-nanstd(pn(:,1:104),1)./sqrt(size(pn,1)))];

            figure
            fill(fillx,fillyr,'b')
            hold on
            fill(fillx,fillyn,'r')
            set(gca,'yscale','log')
            legend('rewarded','non-rewarded')
        end
        
    end
end
