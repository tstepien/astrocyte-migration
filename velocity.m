function [vel_cir,vel_rad] = velocity(j,c1,c2,r,kappa,cmin,rbar,mu)
% vel_cir : circumferential spreading (v/r)
% vel_rad : radial spreading (partial v/partial r)

% at a single time - time step j

%%% spatial mesh
R = length(r);
dr = r(2)-r(1);

%%% initialize
k = c1 + c2;
khalf = (k(1:R-1)+k(2:R))/2;
Tp = Tderivative(k,kappa,cmin,rbar);
Tphalf = Tderivative(khalf,kappa,cmin,rbar);

vel_cir = zeros(1,R);
vel_rad = zeros(1,R);

%%% interior
for i = 2:j-1
    vel_cir(i) = 1/(2*mu*r(i)*dr) * Tp(i) * (k(i+1)-k(i-1));
    vel_rad(i) = 1/(mu*dr^2) * ( Tphalf(i)*(k(i+1)-k(i)) - Tphalf(i-1)*(k(i)-k(i-1)) );
end

%%% origin
vel_cir(1) = 2/(mu*dr^2) * Tphalf(1)*(k(2)-k(1));
vel_rad(1) = vel_cir(1);

%%% moving boundary
vel_cir(j) = 1/(mu*r(j)*dr) * Tp(j)*(k(j)-k(j-1));
vel_rad(j) = 1/(mu*dr^2) * ( Tp(j)*(k(j)-k(j-1)) - Tp(j-1)*(k(j-1)-k(j-2)) );