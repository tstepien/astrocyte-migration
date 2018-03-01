function [c1_new,c2_new] = cellpops_implicit_fulleqntry_bc(j,c1_old,...
    c2_old,p1,p2,...
    dt,r,kappa,mu,alpha1,alpha2,beta,gamma1,gamma2,cmin,rbar,ce)
% [c1_new,c2_new] = cellpops_implicit(j,c1_old,c2_old,p1,p2,...
%     dt,r,kappa,mu,alpha1,alpha2,beta,gamma1,gamma2,cmin,rbar,ce)
%
% inputs:
%   j      = node that moving boundary is located at
%   c1_old = APC cell concentration at previous time
%   c2_old = IPA cell concentration at previous time
%   p1     = PDGFA growth factor concentration at previous time
%   p2     = LIF growth factor concentration at previous time
%   dt     = time step size
%   r      = spatial mesh
%   {others} = parameters
%
% outputs:
%   c1_new   = APC cell concentration at next time
%   c2_new   = IPA cell concentration at next time

%%% spatial mesh
R = length(r);
dr = r(2)-r(1);
rhalf = (r(1:R-1)+r(2:R))/2;

%%% initialize
k = c1_old + c2_old;
khalf = (k(1:R-1)+k(2:R))/2;
c1half = (c1_old(1:R-1)+c1_old(2:R))/2;
c2half = (c2_old(1:R-1)+c2_old(2:R))/2;
Tp = Tderivative(khalf,kappa,cmin,rbar);
Psi1 = rhalf.*c1half.*Tp;
Psi2 = rhalf.*c2half.*Tp;

%%% extrapolate Psi to get node j+1/2
Psi1(j) = interp1(rhalf(1:j-1),Psi1(1:j-1),rhalf(j),'pchip','extrap');
Psi2(j) = interp1(rhalf(1:j-1),Psi2(1:j-1),rhalf(j),'pchip','extrap');

%%% cells at time step j are at mesh points 1:j
%%% cells at time step j are at mesh points 1:j+1

Tprimeatce = Tderivative(ce,kappa,cmin,rbar); % T'(ce)
Tprime2atce = Tsecderivative(ce,kappa,cmin,rbar); % T''(ce)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% c1 - APC %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
theta1_1 = dt./(mu*r(2:j)*dr^2) .* Psi1(2:j);
theta1_2 = dt./(mu*r(2:j)*dr^2) .* ( Psi1(2:j) + Psi1(1:j-1) );
theta1_3 = dt./(mu*r(2:j)*dr^2) .* Psi1(1:j-1);

theta1_4 = 2*dt/(mu*dr^2) * c1half(1)*Tp(1);

upperdiag11 = [theta1_4 , theta1_1];
maindiag11 = [1-theta1_4 , 1-theta1_2 , 1];
lowerdiag11 = [theta1_3 , 0];

upperdiag12 = [theta1_4 , theta1_1];
maindiag12 = [-theta1_4 , -theta1_2 , 1];
lowerdiag12 = [theta1_3, 0];

block11 = diag(maindiag11) + diag(upperdiag11,1) + diag(lowerdiag11,-1);
block12 = diag(maindiag12) + diag(upperdiag12,1) + diag(lowerdiag12,-1);

g1 = [cellgrowth_c1(c1_old(1:j),p1(1:j),p2(1:j),alpha1,beta,gamma1) , 0];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% c2 - IPA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
theta2_1 = dt./(mu*r(2:j)*dr^2) .* Psi2(2:j);
theta2_2 = dt./(mu*r(2:j)*dr^2) .* ( Psi2(2:j) + Psi2(1:j-1) );
theta2_3 = dt./(mu*r(2:j)*dr^2) .* Psi2(1:j-1);

theta2_4 = 2*dt/(mu*dr^2) * c2half(1)*Tp(1);

% Tpinterp = interp1([rhalf(1:j-1),r(j+1)],[Tp(1:j-1),Tprimeatce],rhalf(j),...
%     'pchip','extrap');
% theta5 = -Tprimeatce/(mu*r(j+1)*dr) - Tprimeatce/(mu*dr^2) - Tp(j-1)/(mu*dr^2);
% theta6 = Tp(j-1)/(mu*dr^2);
% theta7 = -ce*( Tprimeatce/(mu*r(j+1)*dr) + Tprimeatce/(mu*dr^2) ) ...
%     + 1/ce*( cellgrowth_c1(c1_old(j),p1(j+1),...
%     p2(j+1),alpha1,beta,gamma1) + cellgrowth_c2(c1_old(j),c2_old(j),...
%     p1(j+1),p2(j+1),alpha2,beta,gamma2) );

% theta5 = -g1plusg2_ord2(p1(j+1),alpha1,alpha2,gamma1,gamma2);
% theta6 = ce/(mu*r(j+1)*dr^2) * ( r(j+1)*Tderivative(k(j),kappa,cmin,rbar) ...
%     - r(j)*Tderivative(k(j-1),kappa,cmin,rbar) );
% theta7 = -ce/(mu*r(j+1)*dr^2) * r(j)*Tderivative(k(j-1),kappa,cmin,rbar);

% theta5 = (1/r(j+1)+1/dr)*ce/(mu*dr)*Tprimeatce;
% theta6 = -ce/(mu*dr)*( Tprimeatce/r(j+1) + Tprimeatce/dr + Tp(j-1)/dr );
% theta7 = ce/(mu*dr^2)*Tp(j-1);

theta5=1;
theta6=0;
theta7=0;

upperdiag21 = [theta2_4 , theta2_1];
maindiag21 = [-theta2_4 , -theta2_2 , 0];
lowerdiag21 = [theta2_3 , theta6];

upperdiag22 = [theta2_4 , theta2_1];
maindiag22 = [1-theta2_4 , 1-theta2_2 , theta5];
lowerdiag22 = [theta2_3 , theta6];

block21 = diag(maindiag21) + diag(upperdiag21,1) + diag(lowerdiag21,-1);
block21(end,end-2) = theta7;
block22 = diag(maindiag22) + diag(upperdiag22,1) + diag(lowerdiag22,-1);
block22(end,end-2) = theta7;

% c2endinterp = interp1(r(j-1:j),c2_old(j-1:j),r(j+1),'pchip','extrap');

g2 = [cellgrowth_c2(c1_old(1:j),c2_old(1:j),p1(1:j),p2(1:j),alpha2,...
    beta,gamma2) , 0];
    
%     cellgrowth_c2(ce-c2endinterp,c2endinterp,p1(j+1),p2(j+1),...
%     alpha2,beta,gamma2)];
% cellgrowth_c2(c1_old(j),c2_old(j),p1(j),p2(j),alpha2,beta,gamma2)];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
thetamatrix = [block11 , block12 ; ...
    block21 , block22];

bc = zeros(2*(j+1),1);
bc(j+1) = ce;
% bc(end) = c2_old(j);
% bc(end) = c2endinterp;
% bc(end) = -ce^2/(mu*dr^2)*Tderivative(k(j),kappa,cmin,rbar) ...
%     + g1plusg2_ord1(ce,p1(j+1),alpha1,gamma1);
% bc(end) = cellgrowth_c1(c1_old(j),p1(j+1),p2(j+1),alpha1,beta,gamma1) ...
%     + cellgrowth_c2(c1_old(j),c2_old(j),p1(j+1),p2(j+1),alpha2,beta,gamma2);
bc(end)=0;

bvector = [c1_old(1:j)' ; 0 ; c2_old(1:j)' ; 0] + bc + dt*[g1,g2]';

c_new = ( thetamatrix \ bvector )';

c1_new = [c_new(1:j+1) , zeros(1,R-(j+1))];
c2_new = [c_new(j+2:end) , zeros(1,R-(j+1))];

c_newT = c_new';

keyboard