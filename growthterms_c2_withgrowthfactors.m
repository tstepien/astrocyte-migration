function g = growthterms_c2_withgrowthfactors(c1,c2,q1,q2,PO2,Pm,alpha2,beta,gamma2,cmax,hyaloid)
% g = growthterms_c2(c1,c2,q1,q2,PO2,Pm,alpha2,beta,gamma2,cmax,hyaloid)
%
% growth terms

hyaloid = hyaloid(1:length(c1));

g = alpha2*(PO2./(Pm+PO2)+q1).*c2.*(1 - (c1+c2)/cmax) ...
    + beta*(hyaloid + PO2./(Pm+PO2)+q2).*c1 ...
    - gamma2.*c2;