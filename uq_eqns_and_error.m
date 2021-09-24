function Y = uq_eqns_and_error(X)
% Y = uq_eqns_and_error(X)
%
% This function returns the computed values of the moving boundary location
% astrocyte migration for a single population of cells and the final
% simulated time (days)
%
% inputs:
%   X = [mu, alpha11, alpha12, alpha21, alpha22, beta1, beta2, beta3, ...
%           gamma1, gamma2, Te, P_hy, r_hy]
%
% output:
%   Y = total error

%%%%%%%%%%%%%%%%%%%%%%%%%% parameters to examine %%%%%%%%%%%%%%%%%%%%%%%%%%
%%% APCs and IPAs
alpha10 = X(:,1);
alpha11 = X(:,2);
alpha12 = X(:,3);
alpha20 = X(:,4);
alpha21 = X(:,5);
alpha22 = X(:,6);
beta1 = X(:,7);
beta2 = X(:,8);
beta3 = X(:,9);

%%%%%%%%%%%%%%%%%%%%%%%%%%%% fixed parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
parameters_fixed

%%% mesh parameters
m.rmax = 5; % max radius (mm) (estimate rat retinal radius = 4.1 mm)
m.tmax = 7 * 24; % max time (hr) (7 days)

%%% tension parameters 
mu1 = 1; %1.75 % adhesion constant (mm/hr/(mN/mm^2))
mu2 = 2; %1.75
kappa = 30; % tension function scaling (Pa, microN/mm^2)
% tension on boundary (mN/mm)
Te = kappa * (1/sqrt(pi*ce) - rbar) ; % simpler form, June 2021
m.kTprime1 = kappa / (2 * mu1 * sqrt(pi));
m.kTprime2 = kappa / (2 * mu2 * sqrt(pi));
% parameter in initial condition, chosen to match subsequent solution
m.ce1 = ce*(1 + 0.5 * alpha11 * s0^2 * sqrt(ce) * (1 - ce/cmax) / m.kTprime1);

%%%%%%%%%%%%%%%%%%% solve equation and calculate error %%%%%%%%%%%%%%%%%%%%
%%% initialize
m.nxpts = 26; %751
m.ntpts = 29;
x = linspace(0,1,m.nxpts);    % 0 <= x <= 1 by definition
t = linspace(0,m.tmax,m.ntpts);

% Set up vectors for interpolation to obtain PO2 
[m.Lvec, m.Pvec] = oxygen_setup(M0, Dalpha, Pm, P0);

% Calculate growth factors for interpolation using fixed domain
r = linspace(0,m.rmax,m.nxpts);
sol = pdepe(1,@GF_PDE,@GF_IC,@GF_BC,r,t);
m.PDGFA = sol(:,:,1);
m.LIF = sol(:,:,2);

N = size(X,1);
Y = zeros(N,1);

for i=1:N
    % Calculate cell densities using time-dependent domain stretching
    sol = pdepe(1,@(x,t,u,DuDx) AstroPDE_uq(x,t,u,DuDx,alpha10(i),...
        alpha11(i),alpha12(i),alpha20(i),alpha21(i),alpha22(i),beta1(i),...
        beta2(i),beta3(i),m),@AstroIC,@AstroBC,x,t);
    
    [Y(i),~,~] = errorfunction(t,x,sol);
end