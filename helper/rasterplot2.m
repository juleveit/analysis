function rasterplot2(resp,condl0,condl1,ta)

binwidth = 30;

indsl0 = find(condl0);
subplot(3,1,1)
hold on
for i = 1:length(indsl0)
    if ~isempty(find(resp(indsl0(i),:)))
        plot(find(resp(indsl0(i),:)),i,'ko','MarkerSize',1.5,'MarkerFaceColor','k')
    end
end
% axis([0, size(resp,2), 0, length(indsl0)+1])
xlim([0 size(resp,2)]),ylim([0 length(indsl0)+1]),set(gca,'visible','off','Xtick',1:1:length(indsl0));

indsl1 = find(condl1);
subplot(3,1,2)
hold on
for i = 1:length(indsl1)
    if ~isempty(find(resp(indsl1(i),:)))
        plot(find(resp(indsl1(i),:)),i,'bo','MarkerSize',1.5,'MarkerFaceColor','b');
    end
end
xlim([0 size(resp,2)]),ylim([1 length(indsl1)+1]),set(gca,'visible','off','Xtick',1:1:length(indsl1));

subplot(3,1,3)
[binl0,bta] = binit(mean(resp(indsl0,:),1),binwidth); binl0 = binl0.*(1000/binwidth);
[binl1,bta] = binit(mean(resp(indsl1,:),1),binwidth); binl1 = binl1.*(1000/binwidth);
% bta = bta-300;
plot(bta,binl0,'k','LineWidth',1.5),...
    xlim([0 size(resp,2)])
hold on
plot(bta,binl1,'b','LineWidth',1.5),...
    xlim([0 size(resp,2)]),set(gca,'box','off');%set(gcf,'Position',[680 558 560/2 418/1.6]);
line([1000,2250],[14,14],'color','k','linewidth',2)
line([1500,2250],[12,12],'color','b','linewidth',2)


