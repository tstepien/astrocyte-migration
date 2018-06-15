function g = growthterms_sum(c1,c2,PO2,Pm,alpha1,alpha2,gamma1,gamma2,rbar)
% g = growthterms_sum(c1,c2,PO2,Pm,alpha1,alpha2,beta,gamma1,gamma2,rbar)
%
% growth terms

cmax = 1/(pi*rbar^2);

g = PO2./(Pm+PO2).*(alpha1*c1 + alpha2*c2).*(1 - (c1+c2)/cmax) ...
    - (gamma1*c1 + gamma2*c2);