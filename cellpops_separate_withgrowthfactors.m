function [c1_new,c2_new] = cellpops_separate_withgrowthfactors(j,c1_old,...
    c2_old,k_new,q1,q2,PO2,dt,r,Pm,kappa,mu,alpha11,alpha12,alpha21,...
    alpha22,beta1,beta2,beta3,gamma1,gamma2,cmin,rbar,ce,cmax,hyaloid)
% [c1_new,c2_new] = cellpops_separate_withgrowthfactors(j,c1_old,...
%     c2_old,k_new,q1,q2,PO2,dt,r,Pm,kappa,mu,alpha11,alpha12,alpha21,...
%     alpha22,beta1,beta2,beta3,gamma1,gamma2,cmin,rbar,ce,cmax,hyaloid)
%
% EXPLICIT method for growth terms
%
% inputs:
%   j      = node that moving boundary is located at
%   c1_old = APC cell concentration at previous time
%   c2_old = IPA cell concentration at previous time
%   k_new  = APC+IPA concentration at next time
%   q1     = PDGFA at next time
%   q2     = LIF at next time
%   PO2    = partial pressure of oxygen at next time
%   dt     = time step size
%   r      = spatial mesh
%   {others} = parameters
%
% outputs:
%   c1_new   = APC cell concentration at next time
%   c2_new   = IPA cell concentration at next time
%
% note: cells at time step j are at mesh points 1:j
%       cells at time step j+1 are at mesh points 1:j+1
%
% %%%%%%%%%%%%%%% adding v_e as unknown ------------------------

global whatstep tcurr;

%%% spatial mesh
R = length(r);
dr = r(2)-r(1);
rhalf = (r(1:R-1)+r(2:R))/2;

%%% initialize
khalf = (k_new(1:R-1)+k_new(2:R))/2;
% k_old = c1_old+c2_old;
% khalf_old = (k_old(1:R-1)+k_old(2:R))/2;
% kinterp = interp1(r,k_new,rhalf,'pchip');
Tp = Tderivative(khalf,kappa,cmin,rbar);
%     Tp(j) = interp1(rhalf(1:j-1),Tp(1:j-1),rhalf(j),'pchip','extrap');
Psi = rhalf.*Tp.*diff(k_new);

% Psi(j) = interp1(rhalf(1:j-1),Psi(1:j-1),rhalf(j),'pchip','extrap');
% Psi(j) = Psi(j-1)/2;

% Psi(j) = 1.81*10^-5;
% Psi(j-1) = 1.79*10^-5;

omega = 1./(2*mu*r(2:j)) * dt/dr^2;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%% growth function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
g1 = growthterms_c1_withgrowthfactors(c1_old(1:j),c2_old(1:j),q1(1:j),...
    q2(1:j),PO2(1:j),Pm,alpha11,alpha12,beta1,beta2,beta3,gamma1,cmax,hyaloid);
g2 = growthterms_c2_withgrowthfactors(c1_old(1:j),c2_old(1:j),q1(1:j),...
    q2(1:j),PO2(1:j),Pm,alpha21,alpha22,beta1,beta2,beta3,gamma2,cmax,hyaloid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%% construct matrix %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% c1 - APC %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% block 11
upperdiag11 = [dt/(mu*dr^2)*Tp(1)*(k_new(2)-k_new(1)) , ...
    omega.*Psi(2:j)];
maindiag11 = [1 + dt/(mu*dr^2)*Tp(1)*(k_new(2)-k_new(1)), ...
    1 + omega.*(Psi(2:j) - Psi(1:j-1)),...
    1];
lowerdiag11 = [-omega.*Psi(1:j-1) , ...
    0];
block11 = diag(maindiag11) + diag(upperdiag11,1) + diag(lowerdiag11,-1);

%%% block 12
block12 = zeros(size(block11));
block12 = [block12 , [zeros(j,1) ; dt*c1_old(j)] ]; % add on right-most column


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% c2 - IPA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% block 21
block21 = zeros(size(block11));
block21 = [block21 ; [zeros(1,j) , 1] ]; % add on bottom row

%%% block 22
upperdiag22 = [dt/(mu*dr^2)*Tp(1)*(k_new(2)-k_new(1)) , ...
    omega.*Psi(2:j)];
maindiag22 = [1 + dt/(mu*dr^2)*Tp(1)*(k_new(2)-k_new(1)),...
    1 + omega .* (Psi(2:j) - Psi(1:j-1)), ...
    1];
lowerdiag22 = [-omega.*Psi(1:j-1) , ...
    0];
block22 = diag(maindiag22) + diag(upperdiag22,1) + diag(lowerdiag22,-1);
block22 = [block22 ; [zeros(1,j) , 1] ]; % add on bottom row
block22 = [block22 , [zeros(j,1) ; dt*c2_old(j) ; 0] ]; % add on right-most column

%%%%%%%%%%%%%%%%%%%%%%%%%%% c1 and c2 together %%%%%%%%%%%%%%%%%%%%%%%%%%%%
thetamatrix = [block11 , block12 ; ...
    block21 , block22];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% right hand side %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bvector = [c1_old(1:j)' ; ...
    c1_old(j) ; ...
    c2_old(1:j)'; ...
    c2_old(j) ; ...
    ce] ...
    + [dt*g1' ; ...
    dt*g1(end); ...
    dt*g2';...
    dt*g2(end); ...
    0];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% solve system %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c_new = ( thetamatrix \ bvector )';
    
c1_new = [c_new(1:j+1) , zeros(1,R-(j+1))];
c2_new = [c_new(j+2:end-1) , zeros(1,R-(j+1))];

c_newT = c_new(1:end-1)';

% if tcurr>3*24
%     disp(whatstep)
%     keyboard
% end

% if min(c1_new)<0
%     keyboard
% end

% if g11(1)<0
%     keyboard
% end

%%% check value of ve obtained
% [c_new(end) , (c1_old(j)-(1-dt*g11(j+1))*c1_new(j+1))/(dt*c1_old(j)) , g11(j+1)]


% if strcmp(whatstep,'corrector')==1
% s=c1_new+c2_new;
% % figure
% if dr>0.01
%     plot(r,s,'-o',r,k_new,'-o')
% else
%     hold on
%     plot(r,s,'Linewidth',1.5)
%     plot(r,k_new,'--','LineWidth',1.5)
%     hold off
% end
% legend('c1+c1','k')
% 
% keyboard
% end


% if strcmp(whatstep,'predictor')==1
%     1
% elseif strcmp(whatstep,'corrector')==1
%     keyboard
% end

% if j==11+5 && c1_new(1)<1000
%     keyboard
% end





%%% debugging - check LHS and RHS
% -1/(mu*r(j)*dr^2) * ...
% ( rhalf(j)*(k_new(j)+k_new(j+1))/2*Tderivative((k_new(j)+k_new(j+1))/2,kappa,cmin,rbar)*(k_new(j+1)-k_new(j)) ...
% - rhalf(j-1)*(k_new(j-1)+k_new(j))/2*Tderivative((k_new(j-1)+k_new(j))/2,kappa,cmin,rbar)*(k_new(j)-k_new(j-1)) ) ...
% +g11(j)*k_new(j)
% 
% (k_new(j)-k_old(j))/dt
% 
% 
% -1/(mu*r(j)*dr^2) * ...
% ( rhalf(j)*(c1_new(j)+c1_new(j+1))/2*Tderivative((k_new(j)+k_new(j+1))/2,kappa,cmin,rbar)*(k_new(j+1)-k_new(j)) ...
% - rhalf(j-1)*(c1_new(j-1)+c1_new(j))/2*Tderivative((k_new(j-1)+k_new(j))/2,kappa,cmin,rbar)*(k_new(j)-k_new(j-1)) ) ...
% +g11(j)*c1_new(j)
% 
% (c1_new(j)-c1_old(j))/dt