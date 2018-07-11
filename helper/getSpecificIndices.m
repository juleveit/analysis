function allinds = getSpecificIndices(result, fieldname)

design = getfield(result,fieldname);
fact = fullfact([length(design.orientations),length(design.sizes),length(design.contrasts),...
    length(design.spfreqs),length(design.light),length(design.tfreqs)]);


for i = 1:size(fact,1)
    conds(i,:) = [design.orientations(fact(i,1)),design.sizes(fact(i,2)),design.contrasts(fact(i,3)),...
        design.spfreqs(fact(i,4)),design.light(fact(i,5)),design.tfreqs(fact(i,6))];
    
    if strcmp(fieldname,'sfconds')
        factors = .04./design.spfreqs;
        alltfreqs = design.tfreqs./factors;
        conds(i,6) = alltfreqs(fact(i,4)); % take the one to the accroding spatial frequency
    end
    
    indices{i} = find((result.gratingInfo.Orientation == conds(i,1)  &...
        result.gratingInfo.size == conds(i,2) & result.gratingInfo.Contrast == conds(i,3) &...
        result.gratingInfo.spFreq == conds(i,4) & result.light == conds(i,5) & result.gratingInfo.tFreq == conds(i,6)) |...
        (result.gratingInfo.Orientation == conds(i,1)+180  &...
        result.gratingInfo.size == conds(i,2) & result.gratingInfo.Contrast == conds(i,3) &...
        result.gratingInfo.spFreq == conds(i,4) & result.light == conds(i,5) & result.gratingInfo.tFreq == conds(i,6)));
    
    nind(i) = length(indices{i});    
end
% add control condition

for j = 1:length(unique(fact(:,5)))
    indices{i+1} = find((result.gratingInfo.Orientation == -1) & result.light == j-1);
    i = i+1;
end

minindn = min(nind);
if minindn ~= result.repetitions
    disp('careful somethings foul here - not the same number of conditions as repetitions')
end

allinds = [];
for i = 1:length(indices)
    rp = randperm(length(indices{i}));
    allinds = [allinds, indices{i}(rp(1:minindn))];
end

disp('');