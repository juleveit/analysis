function val = ppc(phasevec)

if length(phasevec)>1
    k = 1;
    for i = 1:length(phasevec)-1
        for j = i+1:length(phasevec)
            cf(k) = cos(phasevec(i))*cos(phasevec(j)) + sin(phasevec(i))*sin(phasevec(j));
            k = k+1;
        end
    end
    val = (2*sum(cf))/(length(phasevec)*(length(phasevec)-1));
else
    val = NaN;
end