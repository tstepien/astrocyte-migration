function u0 = AstroIC(x)
global ce kTprime1 s0 ce1
% u0 = [ ce * (1 + ce1 * (1.3 - sqrt(0.69 + x^2))); 0;...
%    0.5 * ce1 * kTprime1 / s0 / sqrt(ce); s0];
u0 = [ ce * (1 + ce1 * (1 - x^2)); 0;...
    2 * ce1 * kTprime1 / s0 / sqrt(ce); s0];
end
