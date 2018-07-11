
function ssi = get_ssi(curve)
ssi  = 1 - (((norm(curve)/max(curve)) - 1)./((sqrt(length(curve)))-1));