function hy = hyaloid(r,P_hy,r_hy)
% function hy = hyaloid(r,P_hy,r_hy)
%
% function for the hyaloid artery partial pressure due to oxygen
%
% inputs:
%   r    = spatial grid
%   P_hy = partial pressure due to oxygen at optic nerve head (hyaloid artery)
%   r_hy = radius at half-maximum of Hill function
%
% outputs:
%   hy = function

hy = P_hy*(1- r.^2./ (r_hy^2+r.^2));