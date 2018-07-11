
function a = dirdiff(a)
% a(a<0) = a(a<0)+2*pi;
    a(a<0) = abs(a(a<0));
    a(a>180) = abs(a(a>180)-360);
    
end
