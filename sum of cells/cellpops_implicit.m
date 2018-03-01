function c_new = cellpops_implicit(j,c_old,p1,dt,r,kappa,mu,alpha,gamma,...
    cmin,rbar,ce)
% c_new = cellpops_implicit(j,c_old,p1,dt,r,kappa,mu,alpha,gamma,...
%   cmin,rbar,ce)
%
% inputs:
%   j     = node that moving boundary is located at
%   c_old = APC+IPA cell concentration at previous time
%   p1    = PDGFA growth factor concentration at previous time
%   p2    = LIF growth factor concentration at previous time
%   dt    = time step size
%   r     = spatial mesh
%   {others} = parameters
%
% outputs:
%   c_new = APC+IPA cell concentration at next time

%%% spatial mesh
R = length(r);
dr = r(2)-r(1);
rhalf = (r(1:R-1)+r(2:R))/2;

%%% initialize
chalf = (c_old(1:R-1)+c_old(2:R))/2;
Tp = Tderivative(chalf,kappa,cmin,rbar);
Psi = rhalf.*chalf.*Tp;

%%% extrapolate Psi to get node j+1/2
Psi(j) = interp1(rhalf(1:j-1),Psi(1:j-1),rhalf(j),'pchip','extrap');

%%% cells at time step j are at mesh points 1:j
%%% cells at time step j are at mesh points 1:j+1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
theta1 = 1./(mu*r(2:j)) * dt/dr^2 .* Psi(2:j);
theta2 = 1 - 1./(mu*r(2:j)) * dt/dr^2 .* ( Psi(2:j) + Psi(1:j-1) ) ...
    - dt*( alpha*p1(2:j) - gamma );
theta3 = 1./(mu*r(2:j)) * dt/dr^2 .* Psi(1:j-1);

theta8 = 2/mu * dt/dr^2 * chalf(1)*Tp(1);
theta9 = 1 - 2/mu * dt/dr^2 * chalf(1)*Tp(1) - dt*( alpha*p1(1) - gamma );

maindiag = [theta9 , theta2 , 1];
upperdiag = [theta8 , theta1];
lowerdiag = [theta3 , 0];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
thetamatrix = diag(maindiag) + diag(upperdiag,1) + diag(lowerdiag,-1);

bc = zeros(j+1,1);
bc(j+1) = ce;

bvector = [c_old(1:j)' ; 0] + bc;

c_new = ( thetamatrix \ bvector )';

c_new = [c_new(1:j+1) , zeros(1,R-(j+1))];