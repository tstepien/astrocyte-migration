%   AstrocyteModel_TWS
%   Model for spread of astrocytes of surface of retina
%   Based on work of Tracy Stepien
%   TWS, June 2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Use the MatLab PDE solver PDEPE. The PDE is
%   [A] .*  D_ [u] = D_ [ A * D * Du/Dx ] +      [R]
%           Dt       Dx
%    ----     ---          -------------    -------------
%     c        u          f(x,t,u,Du/Dx)   s(x,t,u,Du/Dx)
%   where R is the reaction rate. The left bc is:
%     [0]    +   [1] .* [ A * D * Du/Dx ] = [0] 
%     ---        ---    -----------------   ---
%   p(0,t,u)     q(0,t)    f(0,t,u,Du/Dx)    0
%   and similarly at x = L

%  Field variables:
%  1 density of APCs
%  2 density of IPAs
%  3 speed of cell boundary (forced to be independent of x)
%  4 position of cell boundary (forced to be independent of x)
%  5 concentration of PDGFA
%  6 concentration of LIF

% Units are mm and hour

clc
clear variables global;
addpath plot_simulations
global kTprime1 kTprime2 ce ce1 cmax Lvec Pvec Pm rmax
global alpha10 alpha11 alpha12 alpha20 alpha21 alpha22 beta1 beta2 beta3 gamma1 gamma2
global LIF PDGFA nxpts ntpts tmax

parameters_fixed;

%% configuration parameters
rmax = 5; % max radius (mm) (estimate rat retinal radius = 4.1 mm)
tmax = 7 * 24; % max time (hr) (7 days)

%% cell growth parameters
alpha10 = 0.08; % (/hr) basal proliferation rate (0.1)
alpha11 = 0.08; %  (/hr) proliferation rate APC wrt oxygen (0.08)
alpha12 = 0.12; %  (/hr) proliferation rate APC wrt PDGFA (0.1)
alpha20 = 0.0; % (/hr) basal proliferation rate (zero) 
alpha21 = 0.0; % (/hr) proliferation rate IPA wrt oxygen
alpha22 = 0.0; % (/hr) proliferation rate IPA wrt PDGFA
beta1 = 0.07; % (/hr) basal differentiation rate (0.08)
beta2 = 0.03; % (/hr) differentiation rate wrt oxygen (0.03)
beta3 = 0.02; % (/hr) differentiation rate wrt LIF (0.02)
gamma1 = 0.0; % (/hr) apoptosis rate APC
gamma2 = 0.0; % (/hr) apoptosis rate IPA

%% tension parameters 
mu1 = 1; %1.75 % adhesion constant (mm/hr/(mN/mm^2))
mu2 = 2; %1.75
kappa = 30; % tension function scaling (Pa, microN/mm^2)
% tension on boundary (mN/mm)
Te = kappa * (1/sqrt(pi*ce) - rbar) ; % simpler form, June 2021
kTprime1 = kappa / (2 * mu1 * sqrt(pi));
kTprime2 = kappa / (2 * mu2 * sqrt(pi));
% parameter in initial condition, chosen to match subsequent solution
ce1 = ce*(1 + 0.5 * alpha11 * s0^2 * sqrt(ce) * (1 - ce/cmax) / kTprime1);

%% set up and solve
nxpts = 26;
ntpts = 29;
x = linspace(0,1,nxpts);    % 0 <= x <= 1 by definition
t = linspace(0,tmax,ntpts);

% Set up vectors for interpolation to obtain PO2 
[Lvec, Pvec] = oxygen_setup(M0, Dalpha, Pm, P0);

% Calculate growth factors for interpolation using fixed domain
r = linspace(0,rmax,nxpts);
sol = pdepe(1,@GF_PDE,@GF_IC,@GF_BC,r,t);
PDGFA = sol(:,:,1);
LIF = sol(:,:,2);

% Calculate cell densities using time-dependent domain stretching
sol = pdepe(1,@AstroPDE,@AstroIC,@AstroBC,x,t);

%% plot results

% plot_the_plots
% plot_the_plots_5panels
plot_the_plots_6panels
% plot_the_plots_APCIPA
% plot_growthfactors_O2_thickness

% plot_timversion

%% calculate error

[err_tot,err_rad,err_dens] = errorfunction(t,x,sol)