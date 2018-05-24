function ve = ve_calc(j,tcurr,r,c1,c2,Pm,alpha1,alpha2,gamma1,gamma2,ce)
% ve = ve_calc(j,tcurr,r,c1,c2,Pm,alpha1,alpha2,gamma1,gamma2,ce)
%
% inputs:
%   j     = node that moving boundary is located at
%   tcurr = current time
%   r     = spatial mesh
%   c1    = APC cell concentration at current time
%   c2    = IPA cell concentration at current time
%   {others} = parameters
%
% outputs:
%   ve = sum of circumferential and radial spreading velocity at edge at
%        current time

PO2 = oxygen(tcurr,r);

g = PO2(j)/(Pm+PO2(j))*(alpha1*c1(j) + alpha2*c2(j)) ...
    - (gamma1*c1(j)+gamma2*c2(j));

ve = 1/ce*g;