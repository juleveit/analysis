
function val = doublegauss(p,x)
    %[C Rp Rn theta, sigma]
%     val = C + Rp* exp(-(dirdist(x-theta)^2)/(2*(sigma^2)) )+ Rn* exp(-(dirdist(x+pi-theta)^2)/(2*(sigma^2)) )
    val = p(1) + p(2)* exp(-(dirdiff(x-p(4)).^2)/(2*(p(5).^2)) ) + p(3)* exp(-(dirdiff(x+180-p(4)).^2)/(2*(p(5).^2)) );
end
