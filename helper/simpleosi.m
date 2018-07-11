function [mo,osi] = simpleosi(s,o)

%orinetation tuning
th = deg2rad(o);
x = sum(cos(2*th).*s);
y = sum(sin(2*th).*s);
meanori = rad2deg(atan(y/x));
if x<0, meanori = meanori+180; end
if meanori<0, meanori = meanori + 360; end
mo = meanori*.5;
osi = sqrt(x^2 + y^2)/sum(s);