function hy = hyaloid(r,Ph)
% function hy = hyaloid(r,Ph)
%
% function for the hyaloid artery partial pressure due to oxygen
%
% inputs:
%   r  = spatial grid
%   Ph = partial pressure due to oxygen at optic nerve head (hyaloid artery)
%
% outputs:
%   hy = function

hy = Ph*(1- r.^2./ ((1/100)^2+r.^2));