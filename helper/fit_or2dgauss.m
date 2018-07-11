function [f, rsquared] = fit_or2dgauss(rf, degperGridPos, rfok)
    %function [f] = autoGaussianSurfML(xi,yi,zi)
    %
    %Fit a surface zi = a*exp(-((xi-x0).^2/2/sigmax^2 + ...
    %                           (yi-y0).^2/2/sigmay^2)) + b
    %
    %On the assumption of Gaussian noise through maximum likelihood
    %(least-squares).
    %
    %The procedure is "robust" in the sense that it is unlikely to get 
    %stuck in local minima which makes it appropriate in fully automated 
    %scenarios. This is accomplished through an initial exhaustive search 
    %for the parameters, followed by refinement with lsqcurvefit
    %
    %Currently only regular grids (generated through meshgrid) are accepted
    %for xi and yi. 
    %Example use:
    %
    %[xi,yi] = meshgrid(-10:10,-20:20);
    %zi = exp(-(xi-3).^2-(yi+4).^2) + randn(size(xi));
    %f = autoGaussianSurfML(xi,yi,zi)
    
    if rfok == 0
        f.amp = NaN;
        f.offs = NaN;
        f.xcenterPix = NaN;
        f.ycenterPix = NaN;
        f.xspreadPix = NaN;
        f.yspreadPix = NaN;
        f.theta = NaN;
        f.thetadeg = NaN;
        f.rmse = NaN;
        f.sse = NaN;

        f.rf = NaN;
        f.xcenterDeg = NaN;
        f.xspreadDeg = NaN;
        f.ycenterDeg = NaN;
        f.yspreadDeg = NaN;
        return;
    end
    
    [a,b] = size(rf);
    [xi,yi] = meshgrid(1:a,1:b);
    zi = reshape(rf,numel(rf),1);
    
    sz = size(zi);
    
    %Verify that the grid is regular
    if any(any(diff(xi,1,1) ~= 0)) || any(any(diff(yi,1,2) ~= 0))
        error('xi or yi is not a regular grid');
    end
    
    xi = xi(:);
    yi = yi(:);
    boundx = [min(xi),max(xi)];
    boundy = [min(yi),max(yi)];
    
    %Find a minimum sigma based on number of elements, range of x and y
    rg = (range(boundx)+range(boundy))/2;
    minsigma = rg/sqrt(length(xi))/5;
    maxsigma = rg;
    sigmas = exp(log(minsigma):.3:log(maxsigma));
    
    rgx = [0:sz(2)/2,-ceil(sz(2)/2)+1:-1]';
    rgy = [0:sz(1)/2,-ceil(sz(1)/2)+1:-1]';
    
    res = zeros(length(sigmas),7);
    
    %Run through all the different values for sigma
    for ii = 1:length(sigmas)
        thefiltx = exp(-rgx.^2/2/sigmas(ii));
        thefilty = exp(-rgy.^2/2/sigmas(ii));
        %convolve zi with gaussian filters and find the maximum response
        %(matched filtering)
        zi2 = reflectconv(reflectconv(zi,thefilty)',thefiltx)';
        [aaa,pos] = max(zi2(:));
        x0 = xi(pos);
        y0 = yi(pos);
        %[y0,x0] = ind2sub(sz,pos);
        
        %Determine the residual error for the optimal x, y for this sigma
        G = exp(-((xi-x0).^2+(yi-y0).^2)/2/sigmas(ii)^2);
        X = [G,ones(length(G),1)];
        ps = X\zi(:);
        res(ii,:) = [sum((zi(:) - X*ps).^2),ps(:)',x0,y0,sigmas(ii),sigmas(ii)];
    end
    
    %Find sigma with corresponding least error
    [aaa,optsigma] = min(res(:,1));
    
    %Fit the parameters again through lsqcurvefit
    opts = optimset('Display','off');
    lb = [-Inf,-Inf,boundx(1),boundy(1),minsigma/100,minsigma/100,0]';
    ub = [ Inf, Inf,boundx(2),boundy(2),maxsigma + 1,maxsigma + 1,pi]';
    [params,resnorm,residual] = lsqcurvefit(@(x,xdata) pointgaussian(x,[xi,yi]),[res(optsigma,2:end),1.4]',xi,zi(:),lb,ub,opts);
    
    restot = sum(sum((rf-mean(mean(rf))).^2));
    rsquared = 1 - (resnorm/restot);
    
    %Collect the f
    f.amp = params(1);
    f.offs = params(2);
    f.xcenterPix = params(3);
    f.ycenterPix = params(4);
    f.xspreadPix = params(5);
    f.yspreadPix = params(6);
    f.theta = pi-params(7);
    f.thetadeg = rad2deg(f.theta);
    f.rmse = sqrt(mean(mean((rf-reshape(pointgaussian(params,[xi,yi]),size(rf,1),size(rf,2))).^2)));
    f.sse = residual;
    f.rsquared = rsquared;
    
    f.rf = reshape(pointgaussian(params,[xi,yi]),size(rf,1),size(rf,2));
    f.xcenterDeg = f.xcenterPix.*degperGridPos;
    f.xspreadDeg = f.xspreadPix.*degperGridPos;
    f.ycenterDeg = f.ycenterPix.*degperGridPos;
    f.yspreadDeg = f.yspreadPix.*degperGridPos;
end

function [thef] = pointgaussian(x,xdata)
    u = cos(x(7))*xdata(:,1) + sin(x(7))*xdata(:,2);
    v = -sin(x(7))*xdata(:,1) + cos(x(7))*xdata(:,2);
    u0 = cos(x(7))*x(3) + sin(x(7))*x(4);
    v0 = -sin(x(7))*x(3) + cos(x(7))*x(4);
    thef = x(1)*exp(-( ((u-u0).^2)/(2*(x(5).^2)) + ((v-v0).^2)/(2*(x(6).^2)) )) + x(2);
end

%Convolution with reflected edge handling
function A = reflectconv(A,f)
    A = bsxfun(@times,fft([A(end:-1:1,:);A;A(end:-1:1,:)]),fft([f(1:floor(end/2));zeros(length(f)*2,1);f(floor(end/2)+1:end)]));
    A = ifft(A);
    A = A(length(f)+1:2*length(f),:);
end