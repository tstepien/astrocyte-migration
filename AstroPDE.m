function [c,f,s] = AstroPDE(x,t,u,DuDx)
global kTprime1 kTprime2 cmax Lvec Pvec Pm LIF PDGFA nxpts ntpts tmax rmax
global alpha10 alpha11 alpha12 alpha20 alpha21 alpha22 beta1 beta2 beta3 gamma1 gamma2

if ~isreal(u) || ~isreal(DuDx)
    c = zeros(size(u));
    f = zeros(size(u));
    s = zeros(size(u));
    return
end

usum = u(1) + u(2);
if usum <= 0
    fprintf('*** Error: negative usum = %f\n',usum);
end
Dusum = DuDx(1) + DuDx(2);

% oxygen, PDGFA, LIF - use previously calculated interpolation functions
[thickness_ret,~,~,~] = thick_rad(t,x * u(4));

nthpts = size(Pvec,1);  % do interpolation by hand for speed
thmax = Lvec(nthpts);
thgrid = 1 + (nthpts - 1) * thickness_ret / thmax;
ith = floor(thgrid);
dth = thgrid - ith;
ith1 = min(ith+1,nthpts);
PO2 = (1-dth) * Pvec(ith) + dth * Pvec(ith1);
% PO2a = interp1(Lvec,Pvec,thickness_ret);
choroid = PO2 / (Pm + PO2);
tgrid = 1 + (ntpts - 1) * t / tmax;
xgrid = 1 + (nxpts - 1) * x * u(4) / rmax;
it = floor(tgrid);
ix = floor(xgrid);
dt = tgrid - it;
dx = xgrid - ix;
% 2D interpolation (Matlab routine is slow)
ix1 = min(ix+1,nxpts);
it1 = min(it+1,ntpts);
PDGFA1 = (1-dt) * ((1-dx) * PDGFA(it,ix) + dx * PDGFA(it,ix1))...
    + dt * ((1-dx) * PDGFA(it1,ix) + dx * PDGFA(it1,ix1));
LIF1 = (1-dt) * ((1-dx) * LIF(it,ix) + dx * LIF(it,ix1))...
    + dt * ((1-dx) * LIF(it1,ix) + dx * LIF(it1,ix1));

g1 = u(1) * ((alpha10 + alpha11 * choroid + alpha12 * PDGFA1) * (1 - usum/cmax) ...
     - gamma1 - (beta1 + beta2 * choroid + beta3 * LIF1));
g2 = u(2) * ((alpha20 + alpha21 * choroid + alpha22 * PDGFA1) * (1 - usum/cmax)...
    - gamma2) + u(1) * (beta1 + beta2 * choroid + beta3 * LIF1);

c = [u(4)^2; u(4)^2; 0; 1];
f = [kTprime1 * u(1) * usum^(-3/2) * Dusum;...
    kTprime2 * u(2) * usum^(-3/2) * Dusum;...
    DuDx(3); DuDx(4)];
s = [x * u(4) * u(3) * DuDx(1) + u(4)^2 * g1;...
    x * u(4) * u(3) * DuDx(2) + u(4)^2 * g2;...
    0; u(3)];
end
