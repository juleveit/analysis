function ph_result = makeNaNph(number, id)

for i = 1:number
    ph_result(i).animalid = id;
    ph_result(i).position = [NaN NaN];
    ph_result(i).cellname = [];
    ph_result(i).tetno = NaN;
    ph_result(i).depth = NaN;
    ph_result(i).interpspike = single(nan(1,311));
    ph_result(i).adiff = NaN;
    ph_result(i).swidth = NaN;
    ph_result(i).vismod = NaN;
    ph_result(i).lightmod = NaN;
    ph_result(i).lfr = NaN;
    ph_result(i).nlfr = NaN;
    ph_result(i).printname = [];
    ph_result(i).condn = nan(2,8,4);
    ph_result(i).condresp = nan(2,8,4,3000);
    ph_result(i).condresperr = nan(2,8,4,3000);
    ph_result(i).allresp = cell(2,8,4);
    ph_result(i).condlfpresp = nan(2,8,4,3000);
    ph_result(i).alllfpresp = cell(2,8,4);
    ph_result(i).condlfpspect = nan(2,8,4,513);
    ph_result(i).condlfpspecterr = nan(2,8,4,513);
    ph_result(i).condfr = nan(2,8,4);
    ph_result(i).conderr = nan(2,8,4);
    ph_result(i).condz = cell(2,8,4);
    ph_result(i).condsc = cell(2,8,4);
    ph_result(i).ff = nan(2,8,4);
    ph_result(i).bincondresp = nan(2,8,4,90);
    ph_result(i).bta = nan(1,90);
    ph_result(i).fax = nan(1,513);
    ph_result(i).binconderr = nan(2,8,4,90);
    ph_result(i).condsparseness = nan(2,8,4);
    ph_result(i).condsparsenesswin = nan(2,8,4);
    ph_result(i).runcondresp = nan(2,8,4,3000);
    ph_result(i).runcondfr = nan(2,8,4);
    ph_result(i).runconderr = nan(2,8,4);
    ph_result(i).runn = nan(2,8,4);
    ph_result(i).stillcondresp = nan(2,8,4,3000);
    ph_result(i).stillcondfr = nan(2,8,4);
    ph_result(i).stillconderr = nan(2,8,4);
    ph_result(i).stilln = nan(2,8,4);
    ph_result(i).condrely = nan(2,8,4);
    ph_result(i).condrelyn = nan(2,8,4);
    ph_result(i).condfiltrely = nan(2,8,4);
    ph_result(i).condlfprely = nan(2,8,4);
    ph_result(i).ordf = nan(1,4); 
    ph_result(i).oris = nan(1,8);
    ph_result(i).anova_dp = NaN;
    ph_result(i).anova_op = NaN; 
    ph_result(i).anova_doip = NaN;
    ph_result(i).omlfpspect = nan(2,4,513);
    ph_result(i).omlfpspecterr = nan(2,4,513);
    ph_result(i).controlresp = nan(2,3000);
    ph_result(i).allcontrolresp = nan(10,3000);
    ph_result(i).controlfr = [NaN, NaN];
    ph_result(i).controlerr = [NaN, NaN];
    ph_result(i).controllfpspect = nan(2,513);
    ph_result(i).binnedctr = nan(2,90);
    ph_result(i).sizecenter = NaN;
    ph_result(i).sizesurround = NaN;
    ph_result(i).anova_lp = NaN;
    ph_result(i).anova_rp = NaN;
    ph_result(i).anova_rlip = NaN;
end

