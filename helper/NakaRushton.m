
function val=NakaRushton(p,x)
% parameters of Naka-Rushton function as in Disney et al., Neuron, 2007
% [R_max, contrast Exponent n,  50%firing-Contrast, spontaneous rate sFR]

val = p(4)+p(1)*((x.^p(2))./(x.^p(2)+p(3).^p(2)));
end