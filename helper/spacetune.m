function sti = spacetune(r)

n = length(r);
% a = mod(0:n-1,n/2)+1;
anglebin = 360/n;
theta = deg2rad(0:anglebin:360-anglebin);

x = sum(cos(theta).*r);
y = sum(sin(theta).*r);
meandir = rad2deg(atan(y/x));
if x<0, meandir = meandir+180; end
if meandir<0, meandir=meandir+360; end
dsi = sqrt(x^2 + y^2)/(sum(r));



