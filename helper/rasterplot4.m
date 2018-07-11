function rasterplot4(resp,cond1l0,cond2l0,cond1l1,cond2l1,ta)

binwidth = 30;

inds1l0 = find(cond1l0);
subplot(5,1,1)
hold on
for i = 1:length(inds1l0)
    if ~isempty(find(resp(inds1l0(i),:)))
        plot(find(resp(inds1l0(i),:)),i,'bo','MarkerSize',1.5,'MarkerFaceColor','b')
    end
end
% axis([0, size(resp,2), 0, length(inds1l0)+1])
xlim([0 size(resp,2)]),ylim([0 length(inds1l0)+1]),set(gca,'visible','off','Xtick',1:1:length(inds1l0));

inds2l0 = find(cond2l0);
subplot(5,1,2)
hold on
for i = 1:length(inds2l0)
    if ~isempty(find(resp(inds2l0(i),:)))
        plot(find(resp(inds2l0(i),:)),i,'ro','MarkerSize',1.5,'MarkerFaceColor','r')
    end
end
% axis([0, size(resp,2), 0, length(inds2l0)+1])
xlim([0 size(resp,2)]),ylim([0 length(inds2l0)+1]),set(gca,'visible','off','Xtick',1:1:length(inds2l0));

inds1l1 = find(cond1l1);
subplot(5,1,3)
hold on
for i = 1:length(inds1l1)
    if ~isempty(find(resp(inds1l1(i),:)))
        plot(find(resp(inds1l1(i),:)),i,'co','MarkerSize',1.5,'MarkerFaceColor','c')
    end
end
% axis([0, size(resp,2), 0, length(inds1l1)+1])
xlim([0 size(resp,2)]),ylim([0 length(inds1l1)+1]),set(gca,'visible','off','Xtick',1:1:length(inds1l1));

inds2l1 = find(cond2l1);
subplot(5,1,4)
hold on
for i = 1:length(inds2l1)
    if ~isempty(find(resp(inds2l1(i),:)))
        plot(find(resp(inds2l1(i),:)),i,'mo','MarkerSize',1.5,'MarkerFaceColor','m')
    end
end
% axis([0, size(resp,2), 0, length(inds2l1)+1])
xlim([0 size(resp,2)]),ylim([0 length(inds2l1)+1]),set(gca,'visible','off','Xtick',1:1:length(inds2l1));


subplot(5,1,5)
[bin1l0,bta] = binit(mean(resp(inds1l0,:),1),binwidth); bin1l0 = bin1l0.*(1000/binwidth);
[bin2l0,bta] = binit(mean(resp(inds2l0,:),1),binwidth); bin2l0 = bin2l0.*(1000/binwidth);
[bin1l1,bta] = binit(mean(resp(inds1l1,:),1),binwidth); bin1l1 = bin1l1.*(1000/binwidth);
[bin2l1,bta] = binit(mean(resp(inds2l1,:),1),binwidth); bin2l1 = bin2l1.*(1000/binwidth);
bta = bta-300;
plot(bta,bin1l0,'b','LineWidth',1.5),...
    xlim([-300 size(resp,2)-300])
hold on
plot(bta,bin2l0,'r','LineWidth',1.5),...
    xlim([-300 size(resp,2)-300])
plot(bta,bin1l1,'c','LineWidth',1.5),...
    xlim([-300 size(resp,2)-300])
plot(bta,bin2l1,'m','LineWidth',1.5),...
    xlim([-300 size(resp,2)-300]),set(gca,'box','off');%set(gcf,'Position',[680 558 560/2 418/1.6]);

