function g = growthterms_c2(c1,c2,PO2,Pm,alpha2,beta,gamma2,rbar)
% g = growthterms_c2(c1,c2,PO2,Pm,alpha2,beta,gamma2,rbar)
%
% growth terms

cmax = 1/(pi*rbar^2);

g = alpha2*PO2./(Pm+PO2).*c2.*(1 - (c1+c2)/cmax) + beta*PO2./(Pm+PO2).*c1 ...
    - gamma2.*c2;