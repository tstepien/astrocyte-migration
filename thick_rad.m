function [thickness_ret,thickness_RGC,radius_endo,radius_ret] = thick_rad(time,r)
% [thickness_ret,thickness_RGC,radius_endo,radius_ret] = thick_rad(time,r)
%
% Calculates the thickness of the retina and the retinal ganglion cell (RGC)
% layer, and the radius of the endothelial cell spread
%
% inputs:
%   time = (tcurr+dt) = current time
%   r    = spatial mesh
%
% outputs:
%   thickness_ret = thickness of the retina
%   thickness_RGC = thickness of the RGC layer
%   radius_endo   = radius of endothelial cell spread
%   radius_ret    = radius of retina spread

%%% convert current time from hours to days
tday = time/24;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% retina %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% thickness at edge of retina (convert from micron to mm)
thickness_peripheral = (13.77 * tday + 72.8) * 0.001;
%%% thickness at center of retina (optic nerve head) (convert from
%%% micron to mm)
thickness_origin = (14.33 * tday + 98.78) * 0.001;
%%% radius of retina (convert from micron to mm)
radius_ret = (414.17 * tday + 1029.17) * 0.001;
%%% parabolic nonuniform thickness throughout retina
thickness_ret = (( thickness_peripheral - thickness_origin )./radius_ret.^2 .*r.^2 ...
    + thickness_origin ) .* (r<=radius_ret);

%%%%%%%%%%%%%%%%%%%%%%%%% retinal ganglion cells %%%%%%%%%%%%%%%%%%%%%%%%%%
%%% thickness of retinal ganglion cell layer (from Braekevelt and Hollenburg)
%%% in microns, converted to mm
thickness_RGC_origin = max(-3.79*tday.^2 + 31.02*tday - 23.16 , 0) * 0.001;
thickness_RGC_peripheral = max(-2.49*tday.^2 + 23.81*tday - 24.12 , 0) * 0.001;

thickness_RGC = max( (thickness_RGC_peripheral-thickness_RGC_origin)./radius_ret.^2 ...
    .* r.^2 + thickness_RGC_origin, 0) .* (r<=radius_ret);


%%%%%%%%%%%%%%%%%%%%%%%%%%%% endothelial cells %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% radius of endothelial cells (in microns, converted to mm)
radius_endo = 185*tday * 0.001;