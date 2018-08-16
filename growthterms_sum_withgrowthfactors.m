function g = growthterms_sum_withgrowthfactors(c1,c2,p1,p2,PO2,Pm,alpha1,alpha2,gamma1,gamma2,rbar)
% g = growthterms_sum(c1,c2,p1,p2,PO2,Pm,alpha1,alpha2,beta,gamma1,gamma2,rbar)
%
% growth terms

global tcurr;

cmax = 1/(pi*rbar^2);

g = (PO2./(Pm+PO2)+p1).*(alpha1*c1 + alpha2*c2).*(1 - (c1+c2)/cmax) ...
    - (gamma1*c1 + gamma2*c2);

% [tcurr/24 max(g) max(PO2./(Pm+PO2)) max(p1)]
% plot(g)
% ylim([0,200])
% pause(0.1)