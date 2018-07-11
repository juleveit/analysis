function rasterplot(resp,ta)

subplot(2,1,1)
hold on;
for i = 1:size(resp,1)
    if ~isempty(find(resp(i,:)))
        plot(find(resp(i,:)),i,'ko','MarkerSize',1.5,'MarkerFaceColor','k')
    end
end
% xlim([0 size(resp,2)]),ylim([0 size(resp,1)+1]),set(gca,'visible','off','Xtick',1:1:size(resp,1));

binwidth = 10;

subplot(2,1,2)
[binr,bta] = binit(mean(resp,1),binwidth); binr = binr.*(1000/binwidth);
plot(bta,binr,'k','LineWidth',1.5)

