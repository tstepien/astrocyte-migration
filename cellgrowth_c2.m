function g = cellgrowth_c2(c1,c2,PO2,Pm,alpha2,beta,gamma2)
% g = cellgrowth_c2(c1,c2,PO2,Pm,alpha2,beta,gamma2)
%
% growth function for IPAs

% g = alpha2*p1.*c2 + beta*p2.*c1 - gamma2*c2;

g = alpha2 * PO2/(Pm+PO2) * c2 + beta * Pm/(Pm+PO2) * c1 - gamma2 * c2;