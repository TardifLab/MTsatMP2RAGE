function uni = calculate_UNI_from_sims(inv1, inv2)
% From Marques et al., 2010

uni = (inv1.*inv2)./(inv1.^2 + inv2.^2);