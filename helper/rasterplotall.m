function rasterplotall(allresp,ta)

hold on;
ii = 1;
for l = 1:size(allresp,1)
    for s = 1:size(allresp,3)
        for o = 1:size(allresp,2)
            for i = 1:size(allresp{l,o,s},1)
                if ~isempty(find(allresp{l,o,s}(i,:)))
                    if l == 1 lcol = 'k'; else lcol = 'r'; end
                    found = find(allresp{l,o,s}(i,:));
                    plot(found-300,ii,'o','color',lcol,'MarkerSize',1,'MarkerFaceColor',lcol)
                    ii = ii+1;
                end
            end
%             line([-299,2700],[ii+.5,ii+.5],'color','k','linestyle',':')
        end        
        line([-299,2700],[ii+.5,ii+.5],'color','r')
    end
end
axis ij 
axis([-300,2300,0,ii])

% xlim([0 size(resp,2)]),ylim([0 size(resp,1)+1]),set(gca,'visible','off','Xtick',1:1:size(resp,1));

