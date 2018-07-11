function l = fitvanmises(x,y)
testx = 1:180;
p = [10,2,pi/2,2];

end

function [p, resnorm, rsquared] = fit_vanMises(x,y)

    %[A k phi offs]
    range = max(y)-min(y);
    if range~=0
        p0 = [range/2, pi, (max(y)+min(y))/2]; 
        lb = [.5*(range/2),0,min(y)];
        ub = [1.5*(range/2),2*pi,max(y)];
        [p,resnorm,residual,exitflag] = lsqcurvefit(@(p,x) vanMises(p,x), p0,x,y,lb,ub,optimset('Display','off'));
        restot = sum((y-mean(y)).^2);
        rsquared = 1 - (resnorm/restot);
    else
        p = [NaN,NaN,NaN];
        resnorm = NaN;
        rsquared = NaN;
    end
end

function val = vanmises(p,x)
    %[A k phi offs]
%     val = A* exp( k * cos( 2* (x-phi))-1) + offs; 
    val = p(1)* exp( p(2) * cos( 2* (x-p(3)))-1) + p(4); 
end