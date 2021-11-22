function [c,f,s] = GF_PDE(r,t,u,DuDr)
% [c,f,s] = GF_PDE(r,t,u,DuDr)
%
% growth factors: partial differential equations

global xi1p xi2p gamma3 gamma4 D1 D2

timeramp = max((t/24 - 3) / 4, 0);
[thickness_ret,thickness_RGC,radius_endo,~] = thick_rad(t,r);
xi1 = xi1p * thickness_RGC * timeramp * (thickness_ret > 0);
% use smoothed representation of step function to avoid jumpy output
xi2 = xi2p/(1 + exp(10*(r - radius_endo)));

c = [1; 1];
f = [D1 * DuDr(1); D2 * DuDr(2)];
s = [xi1 - gamma3 * u(1); xi2 - gamma4 * u(2)];

end