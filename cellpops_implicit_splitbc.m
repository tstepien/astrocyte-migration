function [c1_new,c2_new] = cellpops_implicit_splitbc(j,c1_old,...
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

%%%%%%%%%%%%%%%% adding v_e as unknown

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
%%% cells at time step j+1 are at mesh points 1:j+1

Tprimeatce = Tderivative(ce,kappa,cmin,rbar); % T'(ce)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% c1 - APC %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
theta1_1 = dt./(mu*r(2:j)*dr^2) .* Psi1(2:j);
theta1_2 = dt./(mu*r(2:j)*dr^2) .* ( Psi1(2:j) + Psi1(1:j-1) );
theta1_3 = dt./(mu*r(2:j)*dr^2) .* Psi1(1:j-1);

theta1_4 = 2*dt/(mu*dr^2) * c1half(1)*Tp(1);

theta1_5 = dt * c1_old(j);

upperdiag11 = [theta1_4 , theta1_1];
maindiag11 = [1-theta1_4 , 1-theta1_2 , 1];
lowerdiag11 = [theta1_3 , 0];

block11 = diag(maindiag11) + diag(upperdiag11,1) + diag(lowerdiag11,-1);

upperdiag12 = [theta1_4 , theta1_1];
maindiag12 = [-theta1_4 , -theta1_2 , 0];
lowerdiag12 = [theta1_3, 0];

block12 = diag(maindiag12) + diag(upperdiag12,1) + diag(lowerdiag12,-1);
block12 = [block12 , [zeros(j,1) ; theta1_5] ];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% c2 - IPA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
theta2_1 = dt./(mu*r(2:j)*dr^2) .* Psi2(2:j);
theta2_2 = dt./(mu*r(2:j)*dr^2) .* ( Psi2(2:j) + Psi2(1:j-1) );
theta2_3 = dt./(mu*r(2:j)*dr^2) .* Psi2(1:j-1);

theta2_4 = 2*dt/(mu*dr^2) * c2half(1)*Tp(1);

theta2_5 = dt * c2_old(j);
    
upperdiag21 = [theta2_4 , theta2_1];
maindiag21 = [-theta2_4 , -theta2_2 , 0];
lowerdiag21 = [theta2_3 , 0];

block21 = diag(maindiag21) + diag(upperdiag21,1) + diag(lowerdiag21,-1);
block21 = [block21 ; [zeros(1,j) , 1] ];

upperdiag22 = [theta2_4 , theta2_1];
maindiag22 = [1-theta2_4 , 1-theta2_2 , 1];
lowerdiag22 = [theta2_3 , 0];

block22 = diag(maindiag22) + diag(upperdiag22,1) + diag(lowerdiag22,-1);
block22 = [block22 ; [zeros(1,j) , 1] ];
block22 = [block22 , [zeros(j,1) ; theta2_5 ; 0] ];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% right hand side %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
g1 = [cellgrowth_c1(c1_old(1:j),p1(1:j),p2(1:j),alpha1,beta,gamma1) , ...
    cellgrowth_c1(c1_old(j),p1(j),p2(j),alpha1,beta,gamma1)];

g2 = [cellgrowth_c2(c1_old(1:j),c2_old(1:j),p1(1:j),p2(1:j),alpha2,beta,gamma2) , ...
    cellgrowth_c2(c1_old(j),c2_old(j),p1(j),p2(j),alpha2,beta,gamma2)];
    
bvector = [c1_old(1:j)' ; c1_old(j) ; ...
    c2_old(1:j)' ; c2_old(j) ; ce] ...
    + dt*[g1,g2,0]';

thetamatrix = [block11 , block12 ; ...
    block21 , block22];
    
c_new = ( thetamatrix \ bvector )';
    
c1_new = [c_new(1:j+1) , zeros(1,R-(j+1))];
c2_new = [c_new(j+2:end-1) , zeros(1,R-(j+1))];

c_newT = c_new(1:end-1)';


% keyboard