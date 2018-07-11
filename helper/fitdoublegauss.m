function l = fitdoublegauss(x,y)
testx = 1:360;
p = [2,10,5,90,30];
a = [0    30    60    90   120   150   180   210   240   270   300   330];
b = [0.1667    1.0000    1.3333    0.4000   0   0  2.7500   12.6000 4.7500    1.4000    0.2500   0];

end

function [p, resnorm, rsquared] = fit_doublegauss(x,y)

    %[C Rp Rn theta, sigma]
    [m,i] = max(y); % find index of max response
    half = length(x)/2;
    orth = mod(i+half-1,length(x))+1;  % find orthogonal orientation
    range = max(y)-min(y);
    if range~=0
        p0 = [range, m, y(orth), x(i), 30]; 
        lb = [0,0,0,0,10];
        ub = [1.5*(range/2),1.5*range, 1.5*range,360,180];
        [p,resnorm,residual,exitflag] = lsqcurvefit(@(p,x) doublegauss(p,x), p0,x,y,lb,ub,optimset('Display','off'));
        restot = sum((y-mean(y)).^2);
        rsquared = 1 - (resnorm/restot);
    else
        p = [NaN,NaN,NaN];
        resnorm = NaN;
        rsquared = NaN;
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