function g = growthterms_sum_withgrowthfactors(c1,c2,q1,PO2,Pm,alpha11,...
    alpha12,alpha21,alpha22,gamma1,gamma2,cmax,hyaloid)
% g = growthterms_sum_withgrowthfactors(c1,c2,q1,PO2,Pm,alpha11,...
%     alpha12,alpha21.alpha22,gamma1,gamma2,cmax,hyaloid)
%
% growth terms summed: g=g1+g2

global tcurr;

hyaloid = hyaloid(1:length(c1));
choroid = PO2./(Pm+PO2);

g = ( (alpha11*(hyaloid + choroid) + alpha12*q1).*c1 ...
    + (alpha21*PO2./(Pm+PO2) + alpha22*q1).*c2 ).*(1 - (c1+c2)/cmax) ...
    - (gamma1*c1 + gamma2*c2);

% [tcurr/24 max(g) max(PO2./(Pm+PO2)) max(q1)]
% plot(g)
% ylim([0,200])
% pause(0.1)