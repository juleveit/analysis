function [r,p] = nancorr(a,b)
    newa = []; newb = [];
    for i = 1:length(a)
        if ~isnan(a(i))&&~isnan(b(i))&&~isinf(a(i))&&~isinf(b(i))
            newa = [newa, a(i)];
            newb = [newb, b(i)];
        end
    end
    if ~isempty(newa) & ~isempty(newb)
        [r,p] = corr(newa',newb');
    else
        r = NaN; p = NaN;
    end