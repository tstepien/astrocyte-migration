function [c1_new,c2_new] = cellpopulations_vecalc(j,c1_old,c2_old,k_new,PO2,dt,...
    r,ve_old,Pm,kappa,mu,alpha1,alpha2,beta,gamma1,gamma2,cmin,rbar,ce)
% [c1_new,c2_new] = cellpopulations(j,c1_old,c2_old,k_hat,PO2,dt,...
%     r,Pm,kappa,mu,alpha1,alpha2,beta,gamma1,gamma2,cmin,rbar,ce)
%
% inputs:
%   j      = node that moving boundary is located at
%   c1_old = APC cell concentration at previous time
%   c2_old = IPA cell concentration at previous time
%   k_new  = APC+IPA concentration at next time
%   PO2    = partial pressure of oxygen at next time
%   dt     = time step size
%   r      = spatial mesh
%   {others} = parameters
%
% outputs:
%   c1_new   = APC cell concentration at next time
%   c2_new   = IPA cell concentration at next time
%
%
% %%%%%%%%%%%%%%% adding v_e as unknown ------------------------

global whatstep;

%%% spatial mesh
R = length(r);
dr = r(2)-r(1);
rhalf = (r(1:R-1)+r(2:R))/2;

%%% initialize
khalf = (k_new(1:R-1)+k_new(2:R))/2;
Tp = Tderivative(khalf,kappa,cmin,rbar);
Psi = rhalf.*Tp.*diff(k_new);

omega = dt./(2*mu*r(2:j)*dr^2);

%%% cells at time step j are at mesh points 1:j
%%% cells at time step j+1 are at mesh points 1:j+1


%%%%%%%%%%%%%%%%%%%%%%%%%%%%% growth function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
g11 = alpha1 * PO2./(Pm+PO2) - beta * Pm./(Pm+PO2) - gamma1; %mult by c1
g22 = alpha2 * PO2./(Pm+PO2) - gamma2; %mult by c2
g21 = beta * Pm./(Pm+PO2); %mult by c1

% if strcmp(whatstep,'corrector')==1
%     [PO2,g11,g21,g22]
% end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% c1 - APC %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% block 11
upperdiag11 = [dt/(mu*dr^2)*Tp(1)*(k_new(2)-k_new(1)) , ...
    omega.*Psi(2:j)];
maindiag11 = [1 + dt/(mu*dr^2)*Tp(1)*(k_new(2)-k_new(1)) - dt*g11(1) , ...
    1 + omega.*(Psi(2:j) - Psi(1:j-1)) - dt*g11(2:j) , ...
    1];
lowerdiag11 = [-omega.*Psi(1:j-1) , ...
    0];

block11 = diag(maindiag11) + diag(upperdiag11,1) + diag(lowerdiag11,-1);

%%% block 12
block12 = zeros(size(block11));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% c2 - IPA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% block 21
block21 = diag([-dt*g21(1:j) , 0]);


%%% block 22
upperdiag22 = [dt/(mu*dr^2)*Tp(1)*(k_new(2)-k_new(1)) , ...
    omega.*Psi(2:j)];
maindiag22 = [1 + dt/(mu*dr^2)*Tp(1)*(k_new(2)-k_new(1)) - dt*g22(1) , ...
    1 + omega .* (Psi(2:j) - Psi(1:j-1)) - dt*g22(2:j) , ...
    1];
lowerdiag22 = [-omega.*Psi(1:j-1) , ...
    0];

block22 = diag(maindiag22) + diag(upperdiag22,1) + diag(lowerdiag22,-1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%% construct matrix %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% right hand side %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
thetamatrix = [block11 , block12 ; ...
    block21 , block22];
    
bvector = [c1_old(1:j)' ; ...
    c1_old(j) + dt*( -ve_old*c1_old(j) + g11(j)*c1_old(j) ); ...
    c2_old(1:j)' ; ...
    c2_old(j) + dt*( -ve_old*c2_old(j) + g21(j)*c1_old(j) + g22(j)*c2_old(j) )] ;
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% solve system %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c_new = ( thetamatrix \ bvector )';
    
c1_new = [c_new(1:j+1) , zeros(1,R-(j+1))];
c2_new = [c_new(j+2:end) , zeros(1,R-(j+1))];

c_newT = c_new(1:end-1)';

% if g11(1)<0
%     keyboard
% end

keyboard


% if strcmp(whatstep,'predictor')==1
%     1
% elseif strcmp(whatstep,'corrector')==1
%     keyboard
% end

% if j==11+5 && c1_new(1)<1000
%     keyboard
% end