function rasterplotallmv(allresp,ta)

hold on;
ii = 1;
for l = 1:size(allresp,1)
    for s = 1:size(allresp,2)
        for i = 1:size(allresp{l,s},1)
            if ~isempty(find(allresp{l,s}(i,:)))
                if l == 1 lcol = 'k'; else lcol = 'r'; end
                found = find(allresp{l,s}(i,:));
                plot(found-300,ii,'o','color',lcol,'MarkerSize',1.5,'MarkerFaceColor',lcol)
                ii = ii+1;
            end
        end
        line([-299,2700],[ii+.5,ii+.5],'color','r')
    end
end
axis ij
axis([-300,3300,0,ii])
line([1000,1000],[0,ii],'color','k')
line([2000,2000],[0,ii],'color','k')