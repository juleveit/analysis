function [oriprefratio, dirprefratio, prefori, meanori, osi, meandir, dsi] = getOSI(r,o)
    %prefratio: ratio of mean(max response, 180)/mean(orthogonals)
    %orientation: max response orientation mod 180
    %meanori: vector mean angle
    %tunestrength: vector mean length / sum of all responses
    if isnan(r)
        oriprefratio = NaN; dirprefratio = NaN; prefori = NaN; meanori = NaN; osi = NaN; meandir = NaN; dsi = NaN;
        return;
    end
    n = length(o);
    a = mod(0:n-1,n/2)+1;
    theta = deg2rad(o);
    
    prefs = find(a == a(find(r == max(r),1)));
    ortho = mod((prefs-1)-n/4,n)+1 ;
    pref = find(r == max(r),1);
    opposite = mod((pref-1)+n/2,n)+1;
    oriprefratio = (r(pref)-mean(r(ortho)))/(r(pref)+mean(r(ortho)));
    dirprefratio = (r(pref)-r(opposite))/(r(pref)+r(opposite));
    prefori = o(prefs(1));
    
    %direction tuning
    x = sum(cos(theta).*r);
    y = sum(sin(theta).*r);
    meandir = rad2deg(atan(y/x));
    if x<0, meandir = meandir+180; end
    if meandir<0, meandir=meandir+360; end
    dsi = sqrt(x^2 + y^2)/(sum(r));
    
    %orinetation tuning
    s = reshape(r,n/2,2);
    s = mean(s,2);
    th = theta(1:n/2);
    x = sum(cos(2*th).*s');
    y = sum(sin(2*th).*s');
    meanori = rad2deg(atan(y/x));
    if x<0, meanori = meanori+180; end
    if meanori<0, meanori = meanori + 360; end
    meanori = meanori*.5;
    osi = sqrt(x^2 + y^2)/sum(s);
end