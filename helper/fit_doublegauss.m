

function [p, resnorm, rsquared] = fit_doublegauss(x,y)

    %[C Rp Rn theta, sigma]
    [m,ind] = max(y); % find index of max response
    half = length(x)/2;
    orth = mod(ind+half-1,length(x))+1;  % find orthogonal orientation
    range = max(y)-min(y);
    if ~isempty(find(isnan(y))) || range==0 || length(x)<8
        p = [NaN,NaN,NaN,NaN,NaN];
        resnorm = NaN;
        rsquared = NaN;
    else
        alph = mean(diff(x));
        sigmas = [alph/2, alph, 40, 60, 90];
        for i = 1:length(sigmas)
            p0 = [range, m, y(orth), x(ind), sigmas(i)]; 
            lb = [-m,0,0,0,alph/2];
            ub = [m,1.5*range, 1.5*range,360,180];
            [p,resnorm,residual,exitflag] = lsqcurvefit(@(p,x) doublegauss(p,x), p0,x,y,lb,ub,optimset('Display','off'));
            restot = sum((y-mean(y)).^2);
            rs(i) = 1 - (resnorm/restot);
        end
        [minv,minind] = min(rs);
        p0 = [range, m, y(orth), x(ind), sigmas(minind)];
        lb = [0,0,0,0,alph/2];
        ub = [m,1.5*range, 1.5*range,360,180];
        [p,resnorm,residual,exitflag] = lsqcurvefit(@(p,x) doublegauss(p,x), p0,x,y,lb,ub,optimset('Display','off'));
        restot = sum((y-mean(y)).^2);
        rsquared = 1 - (resnorm/restot);
    end
end

function val = doublegauss(p,x)
    %[C Rp Rn theta, sigma]
%     val = C + Rp* exp(-(dirdist(x-theta)^2)/(2*(sigma^2)) )+ Rn* exp(-(dirdist(x+pi-theta)^2)/(2*(sigma^2)) )
    val = p(1) + p(2)* exp(-(dirdiff(x-p(4)).^2)/(2*(p(5).^2)) ) + p(3)* exp(-(dirdiff(x+180-p(4)).^2)/(2*(p(5).^2)) );
end

function a = dirdiff(a)
% a(a<0) = a(a<0)+2*pi;
    a(a<0) = abs(a(a<0));
    a(a>180) = abs(a(a>180)-360);
    
end
