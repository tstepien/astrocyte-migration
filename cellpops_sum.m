function k_new = cellpops_sum(j,c1_old,c2_old,PO2,dt,r,Pm,kappa,mu,...
    alpha1,alpha2,gamma1,gamma2,cmin,rbar,ce)
% c_new = cellpops_sum(j,c1_old,c2_old,PO2,dt,r,Pm,kappa,mu,...
%     alpha1,alpha2,gamma1,gamma2,cmin,rbar,ce)
%
% inputs:
%   j      = node that moving boundary is located at
%   c1_old = APC cell concentration at previous time
%   c2_old = IPA cell concentration at previous time
%   PO2    = partial pressure of oxygen at next time
%   dt     = time step size
%   r      = spatial mesh
%   {others} = parameters
%
% outputs:
%   k_new = APC+IPA cell concentration at next time


%%% spatial mesh
R = length(r);
dr = r(2)-r(1);
rhalf = (r(1:R-1)+r(2:R))/2;

%%% initialize
k_old = c1_old + c2_old;
khalf = (k_old(1:R-1)+k_old(2:R))/2;
% keyboard
% Tp = Tderivative(khalf,kappa,cmin,rbar);
% Psi = rhalf.*Tp.*khalf;

%%% extrapolate Psi to get node j+1/2
% Psi(j) = interp1(rhalf(1:j-1),Psi(1:j-1),rhalf(j),'pchip','extrap');
% Psi(j) = interp1(rhalf(j-2:j-1),Psi(j-2:j-1),rhalf(j),'pchip','extrap');

rinterp = linspace(0,r(j),j+1);
rinterp_half = (rinterp(1:end-1)+rinterp(2:end))/2;
kinterp = [interp1(r(1:j),k_old(1:j),rinterp_half,'pchip') , zeros(1,R-j-1)];
Tp = Tderivative(kinterp,kappa,cmin,rbar);
Psi = rhalf.*Tp.*kinterp;

%%% cells at time step j are at mesh points 1:j
%%% cells at time step j are at mesh points 1:j+1

omega = 1./(mu*r(2:j)) * dt/dr^2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% growth function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
g = PO2(1:j)./(Pm+PO2(1:j)).*(alpha1*c1_old(1:j) + alpha2*c2_old(1:j)) ...
    - (gamma1*c1_old(1:j)+gamma2*c2_old(1:j));


%%%%%%%%%%%%%%%%%%%%%%%%%%%% construct matrix %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% upperdiag = [2/mu * dt/dr^2 * khalf(1)*Tp(1) , ...
%     omega .* Psi(2:j)];
% maindiag = [1 - 2/mu * dt/dr^2 * khalf(1)*Tp(1) , ...
%     1 - omega .* ( Psi(2:j) + Psi(1:j-1) ) , ...
%     1];
upperdiag = [2/mu * dt/dr^2 * kinterp(1)*Tp(1) , ...
    omega .* Psi(2:j)];
maindiag = [1 - 2/mu * dt/dr^2 * kinterp(1)*Tp(1) , ...
    1 - omega .* ( Psi(2:j) + Psi(1:j-1) ) , ...
    1];
lowerdiag = [omega .* Psi(1:j-1) , ...
    0];

thetamatrix = diag(maindiag) + diag(upperdiag,1) + diag(lowerdiag,-1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%% right hand side %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bvector = [k_old(1:j)' + dt*g' ; ce];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% solve system %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k_new = ( thetamatrix \ bvector )';

k_new = [k_new(1:j+1) , zeros(1,R-(j+1))];

% keyboard