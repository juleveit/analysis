function psthcondplot(condresp,condresperr,bta)

figure
ncond1 = size(condresp,2);
ncond2 = size(condresp,3);

for i = 1:ncond1
    for j = 1:ncond2
        subplot(ncond1,ncond2,(i-1)*ncond2+j)
%         [binl0,bta] = binit(condresp(1,i,j,:),binwidth); binl0 = binl0.*(1000/binwidth);
%         [binl1,bta] = binit(condresp(2,i,j,:),binwidth); binl1 = binl1.*(1000/binwidth);
%         [binerrl0,bta] = binit(condresperr(1,i,j,:),binwidth);
%         [binerrl1,bta] = binit(condresperr(2,i,j,:),binwidth); 
%         bta = bta-300;
        boundedline(bta,squeeze(condresp(1,i,j,:)),squeeze(condresperr(1,i,j,:)),'k');
        if size(condresp,1)>1
            boundedline(bta,squeeze(condresp(2,i,j,:)),squeeze(condresperr(2,i,j,:)),'r');
        end
        mx = max([max(max(max(max(condresp,[],1),[],2),[],3),[],4),0.1]);
        axis([bta(1),bta(end),0,mx]);
        line([0,0],[0,mx],'color','k')
        line([2000,2000],[0,mx],'color','k')
        line([500,500],[0,mx],'color','r')
        line([1500,1500],[0,mx],'color','r') 
    end
end
