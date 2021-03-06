function sz_result = makeNaNsz(number, id)

for i = 1:number
    sz_result(i).animalid = id;
    sz_result(i).position = [NaN NaN];
    sz_result(i).cellname = [];
    sz_result(i).tetno = NaN;
    sz_result(i).depth = NaN;
    sz_result(i).interpspike = single(nan(1,311));
    sz_result(i).adiff = NaN;
    sz_result(i).swidth = NaN;
    sz_result(i).vismod = NaN;
    sz_result(i).lightmod = NaN;
    sz_result(i).lfr = NaN;
    sz_result(i).nlfr = NaN;
    sz_result(i).printname = [];
    sz_result(i).condn = nan(2,8,5);
    sz_result(i).condresp = nan(2,8,5,3000);
    sz_result(i).condresperr = nan(2,8,5,3000);
    sz_result(i).allresp = cell(2,8,5);
    sz_result(i).condlfpresp = nan(2,8,5,3000);
    sz_result(i).alllfpresp = cell(2,8,5);
    sz_result(i).condlfpspect = nan(2,8,5,513);
    sz_result(i).condfr = nan(2,8,5);
    sz_result(i).conderr = nan(2,8,5);
    sz_result(i).condz = cell(2,8,5);
    sz_result(i).condsc = cell(2,8,5);
    sz_result(i).ff = nan(2,8,5);
    sz_result(i).bincondresp = nan(2,8,5,90);
    sz_result(i).bta = nan(1,90);
    sz_result(i).binconderr = nan(2,8,5,90);
    sz_result(i).condsparseness = nan(2,8,5);
    sz_result(i).condsparsenesswin = nan(2,8,5);
    sz_result(i).runcondresp = nan(2,8,5,3000);
    sz_result(i).runcondfr = nan(2,8,5);
    sz_result(i).runconderr = nan(2,8,5);
    sz_result(i).runn = nan(2,8,5);
    sz_result(i).stillcondresp = nan(2,8,5,3000);
    sz_result(i).stillcondfr = nan(2,8,5);
    sz_result(i).stillconderr = nan(2,8,5);
    sz_result(i).stilln = nan(2,8,5);
    sz_result(i).condrely = nan(2,8,5);
    sz_result(i).condrelyn = nan(2,8,5);
    sz_result(i).condfiltrely = nan(2,8,5);
    sz_result(i).condlfprely = nan(2,8,5);
    sz_result(i).condoriprefratio = nan(2,5);
    sz_result(i).conddirprefratio = nan(2,5);
    sz_result(i).condmeanori = nan(2,5);
    sz_result(i).condosi = nan(2,5);
    sz_result(i).condmeandir = nan(2,5);
    sz_result(i).conddsi = nan(2,5);
    sz_result(i).gaussparams = nan(2,5,5);
    sz_result(i).gaussrsquared = nan(2,5);
    sz_result(i).sizes = nan(1,5); 
    sz_result(i).oris = nan(1,8);
    sz_result(i).omlfpspect = nan(2,5,513);
    sz_result(i).omlfpspecterr = nan(2,5,513);
    sz_result(i).anova_sp = NaN;
    sz_result(i).anova_op = NaN; 
    sz_result(i).anova_soip = NaN;
    sz_result(i).controlresp = nan(2,3000);
    sz_result(i).allcontrolresp = nan(10,3000);
    sz_result(i).controlfr = [NaN, NaN];
    sz_result(i).controlerr = [NaN, NaN];
    sz_result(i).controllfpspect = nan(2,513);
    sz_result(i).binnedctr = nan(2,90);
    sz_result(i).r0contolresp = nan(2,3000);
    sz_result(i).r0controlfr = [NaN, NaN];
    sz_result(i).r0controlerr = [NaN, NaN];
    sz_result(i).r0controln = [NaN, NaN];
    sz_result(i).r1contolresp = nan(2,3000);
    sz_result(i).r1controlfr = [NaN, NaN];
    sz_result(i).r1controlerr = [NaN, NaN];
    sz_result(i).r1controln = [NaN, NaN];    
    sz_result(i).prefsize = NaN;
    sz_result(i).anova_lp = NaN;
    sz_result(i).anova_rp = NaN;
    sz_result(i).anova_rlip = NaN;
end