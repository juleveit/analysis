function S1LFPana
clear all
% scotts
% load('C:\Users\Julia\work\data\others\Yoolizzle3.mat')
% 
%  for i = 1:size(LFP200,1)   
%      [lfpspect(i,1,:),fax] = pmtm(LFP200(i,:),3,[],1000); 
%      [lfpspect(i,2,:),fax] = pmtm(LFP300(i,:),3,[],1000); 
%      [lfpspect(i,3,:),fax] = pmtm(LFP400(i,:),3,[],1000); 
%      [lfpspect(i,4,:),fax] = pmtm(LFP500(i,:),3,[],1000); 
%      [lfpspect(i,5,:),fax] = pmtm(LFP600(i,:),3,[],1000); 
%      [lfpspect(i,6,:),fax] = pmtm(LFP700(i,:),3,[],1000); 
%      [lfpspect(i,7,:),fax] = pmtm(LFP800(i,:),3,[],1000); 
%      [lfpspect(i,8,:),fax] = pmtm(LFP900(i,:),3,[],1000); 
%      speed(i,:) = runspeed{i}(:);
%  end
%  
%  meanspeed = mean(speed,2);
%  stdspeed = std(speed,1,2);
%  notstill = find(meanspeed>1);
%  okspeed = find(meanspeed>( mean(meanspeed(notstill))-(1.5*std(meanspeed(notstill))) ) );
%  okvar = find(stdspeed<( mean(stdspeed(notstill))+(1.5*std(stdspeed(notstill)))) & stdspeed>.5);
%  oktrials = intersect(okspeed,okvar);
%  nonoktrials = 1:size(LFP200,1); nonoktrials(oktrials) = [];
%  stilltrials = 1:size(LFP200,1); stilltrials(notstill) = [];
 
 
 %S1M1 Gregs
 if exist('C:\Users\Julia\work\data\others\gregs\s1m1lfp\m1lfp.mat') == 0
     load('C:\Users\Julia\work\data\others\gregs\s1m1lfp\0986-0987-29-May-2015_e1.phy.mat')
     for i = 1:length(MCdata)
         m1suplfp(i,:) = resample(MCdata{i}(:,13),1,30);
         m1inflfp(i,:) = resample(MCdata{i}(:,14),1,30);
     end
     save('C:\Users\Julia\work\data\others\gregs\s1m1lfp\m1lfp.mat','m1suplfp','m1inflfp')
 else
     load('C:\Users\Julia\work\data\others\gregs\s1m1lfp\m1lfp.mat')
 end
 
 if exist('C:\Users\Julia\work\data\others\gregs\s1m1lfp\s1lfp.mat') == 0
     load('C:\Users\Julia\work\data\others\gregs\s1m1lfp\0986-0987-29-May-2015_e2.phy.mat')
     for i = 1:length(MCdata)
         s1suplfp(i,:) = resample(MCdata{i}(:,9),1,30);
     end
     save('C:\Users\Julia\work\data\others\gregs\s1m1lfp\s1lfp.mat','s1suplfp')
 else
     load('C:\Users\Julia\work\data\others\gregs\s1m1lfp\s1lfp.mat')
 end
 
 load('C:\Users\Julia\work\data\others\gregs\s1m1lfp\0986-0987-29-May-2015.dat.mat')
 for i = 1:size(run_data,2)
     runspeed(i,:) = resample(run_data(:,1),1,30);
 end
 
 for i = 1:size(m1suplfp,1)   
     [m1suplfpspect(i,:),fax] = pmtm(m1suplfp(i,1000:2500),3,[],1000); 
     [m1inflfpspect(i,:),fax] = pmtm(m1inflfp(i,1000:2500),3,[],1000); 
     [s1suplfpspect(i,:),fax] = pmtm(s1suplfp(i,1000:2500),3,[],1000); 
     speed(i,:) = runspeed(i,1001:2500);
 end
 
 
 meanspeed = mean(speed,2);
 stdspeed = std(speed,1,2);
 notstill = find(meanspeed>1);
 okspeed = find(meanspeed>( mean(meanspeed(notstill))-(1.5*std(meanspeed(notstill))) ) );
 okvar = find(stdspeed<( mean(stdspeed(notstill))+(1.5*std(stdspeed(notstill)))) & stdspeed>.5);
 oktrials = intersect(okspeed,okvar);
 nonoktrials = 1:size(m1suplfp,1); nonoktrials(oktrials) = [];
 stilltrials = 1:size(m1suplfp,1); stilltrials(notstill) = [];
 