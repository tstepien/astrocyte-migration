function g = growthterms_c2_withgrowthfactors(c1,c2,p1,p2,PO2,Pm,alpha2,beta,gamma2,rbar)
% g = growthterms_c2(c1,c2,p1,p2,PO2,Pm,alpha2,beta,gamma2,rbar)
%
% growth terms

cmax = 1/(pi*rbar^2);

g = alpha2*(PO2./(Pm+PO2)+p1).*c2.*(1 - (c1+c2)/cmax) ...
    + (beta*PO2./(Pm+PO2)+p2).*c1 ...
    - gamma2.*c2;