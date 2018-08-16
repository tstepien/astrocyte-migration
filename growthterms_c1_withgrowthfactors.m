function g = growthterms_c1_withgrowthfactors(c1,c2,p1,p2,PO2,Pm,alpha1,beta,gamma1,rbar)
% g = growthterms_c1(c1,c2,p1,p2,PO2,Pm,alpha1,beta,gamma1,rbar)
%
% growth terms

cmax = 1/(pi*rbar^2);

g = alpha1*(PO2./(Pm+PO2)+p1).*c1.*(1 - (c1+c2)/cmax) ...
    - beta*(PO2./(Pm+PO2)+p2).*c1 ...
    - gamma1*c1;