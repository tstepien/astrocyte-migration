function ce = densityatbdy(Te,kappa,cmin,rbar)
% ce = densityatbdy(Te,kappa,cmin,rbar)
%
% Numerically calculate the solution of the inverse problem T(ce)=Te for ce

cmax = 1/(pi*rbar^2);

syms cc

S = vpasolve(kappa*cc^2/(cmin^2 + cc^2)*(1/sqrt(pi*cc) - rbar)==Te, ...
    cc, [cmin,cmax]);

ce = double(S);