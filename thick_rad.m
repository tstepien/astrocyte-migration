function [thickness_ret,thickness_RGC,radius_endo,radius_ret] = thick_rad(t,r)
% Calculates the thickness of the retina and the retinal ganglion cell (RGC)
% layer, and the radius of the endothelial cell spread
% inputs:
%   t = current time
%   r = radial position
%
% outputs:
%   thickness_ret = thickness of the retina
%   thickness_RGC = thickness of the RGC layer
%   radius_endo   = radius of endothelial cell spread
%   radius_ret    = radius of retina spread

global maxRGCthick;

tday = t/24;    % time in days

% thickness at edge of retina (mm)
thickness_peripheral = (13.77 * tday + 72.8) * 0.001;
% thickness at center of retina (optic nerve head) (mm)
thickness_origin = (14.33 * tday + 98.78) * 0.001;
% radius of retina (mm)
radius_ret = (414.17 * tday + 1029.17) * 0.001;
% parabolic thickness throughout retina
thickness_ret = ((thickness_peripheral - thickness_origin)./radius_ret.^2 .* r.^2 ...
    + thickness_origin ) .* (r<=radius_ret);

% thickness of retinal ganglion cell layer (Braekevelt and Hollenburg) (mm)
thickness_RGC_origin = max(-3.79*tday.^2 + 31.02*tday - 23.16, 0) * 0.001;
thickness_RGC_peripheral = max(-2.49*tday.^2 + 23.81*tday - 24.12, 0) * 0.001;
thickness_RGC = max((thickness_RGC_peripheral-thickness_RGC_origin)./radius_ret.^2 ...
    .* r.^2 + thickness_RGC_origin, 0) .* (r<=radius_ret);
thickness_RGC = thickness_RGC/maxRGCthick; % normalized

% radius of endothelial cells (mm)
radius_endo = 0.185 * tday;