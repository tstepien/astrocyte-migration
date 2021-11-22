function u0 = AstroIC(x)
% u0 = AstroIC(x)
%
% astrocytes: initial conditions

global ce kTprime1 s0 ce1

u0 = [ ce + (ce1 - ce) * (1 - x^2); 0;...
   2 * (ce1/ce - 1) * kTprime1 / s0 / sqrt(ce); s0];

end
