function g = cellgrowth_c1(c1,PO2,Pm,alpha1,beta,gamma1)
% g = cellgrowth_c1(c1,PO2,Pm,alpha1,beta,gamma1)
%
% growth function for APCs

% g = alpha1*p1.*c1 - beta*p2.*c1 - gamma1*c1;

g = alpha1 * PO2/(Pm+PO2) * c1 - beta * Pm/(Pm+PO2) * c1 - gamma1 * c1;