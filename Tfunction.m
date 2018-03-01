function T = Tfunction(c,kappa,cmin,rbar)
% T = Tfunction(c,kappa,cmin,rbar)
%
% Tension function evaluated at c: T(c)

T = kappa * c^2/(cmin^2+c^2) * ( 1/sqrt(pi*c) - rbar);