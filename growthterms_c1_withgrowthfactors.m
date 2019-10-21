function g = growthterms_c1_withgrowthfactors(c1,c2,q1,q2,PO2,Pm,alpha1,beta,gamma1,cmax)
% g = growthterms_c1(c1,c2,q1,q2,PO2,Pm,alpha1,beta,gamma1,cmax)
%
% growth terms

g = alpha1*(PO2./(Pm+PO2)+q1).*c1.*(1 - (c1+c2)/cmax) ...
    - beta*(PO2./(Pm+PO2)+q2).*c1 ...
    - gamma1*c1;