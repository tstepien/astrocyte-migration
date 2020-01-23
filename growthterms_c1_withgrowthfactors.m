function g = growthterms_c1_withgrowthfactors(c1,c2,q1,q2,PO2,Pm,alpha1,beta,gamma1,cmax,hyaloid)
% g = growthterms_c1(c1,c2,q1,q2,PO2,Pm,alpha1,beta,gamma1,cmax,hyaloid)
%
% growth terms

hyaloid = hyaloid(1:length(c1));

g = alpha1*(PO2./(Pm+PO2)+q1).*c1.*(1 - (c1+c2)/cmax) ...
    - beta*(hyaloid + PO2./(Pm+PO2)+q2).*c1 ...
    - gamma1*c1;