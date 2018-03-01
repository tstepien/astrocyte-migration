function Tp = Tderivative(c,kappa,cmin,rbar)
% Tp = Tderivative(c,kappa,cmin,rbar)
%
% Derivative of tension function evaluated at c: T'(c)

Tp = -kappa/2 * sqrt(c/pi) ...
    .* ( c.^2 - 3*cmin^2 + 4*rbar*cmin^2.*sqrt(pi*c) )./(cmin^2+c.^2).^2;