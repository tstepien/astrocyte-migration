function PO2 = oxygen(r,thickness_ret,P0,Dalpha,M0)
% PO2 = oxygen(r,thickness_ret,P0,Dalpha,M0)
% choroid PO2
%
% inputs:
%   r             = spatial mesh
%   thickness_ret = total retinal thickness (mm)
%   {others}      = parameters
%
% outputs:
%   PO2 = partial pressure of oxygen (O2)

%%% partial pressure of O2
ind = (0 < thickness_ret) & (thickness_ret <= sqrt(2*P0*Dalpha/M0));
if length(thickness_ret)==1
    PO2 = ( P0 - M0/(2*Dalpha)*thickness_ret.^2 ) .*ind .* ones(size(r));
else
    PO2 = ( P0 - M0/(2*Dalpha)*thickness_ret.^2 ) .*ind;
end