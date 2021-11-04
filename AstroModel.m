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
alpha12 = 0.09; %  (/hr) proliferation rate APC wrt PDGFA (0.1)
alpha20 = 0.0; % (/hr) basal proliferation rate (zero) 
alpha21 = 0.005; % (/hr) proliferation rate IPA wrt oxygen
alpha22 = 0.01; % (/hr) proliferation rate IPA wrt PDGFA
beta1 = 0.07; % (/hr) basal differentiation rate (0.08)
beta2 = 0.03; % (/hr) differentiation rate wrt oxygen (0.03)
beta3 = 0.02; % (/hr) differentiation rate wrt LIF (0.02)
gamma1 = 0.0; % (/hr) apoptosis rate APC
gamma2 = 0.0; % (/hr) apoptosis rate IPA

%%% from september 27 afternoon
% alpha10 = 0.0746; % (/hr) basal proliferation rate (0.1)
% alpha11 = 0.0667; %  (/hr) proliferation rate APC wrt oxygen (0.08)
% alpha12 = 0.1005; %  (/hr) proliferation rate APC wrt PDGFA (0.1)
% alpha20 = 0.0006; % (/hr) basal proliferation rate (zero) 
% alpha21 = 0.0049; % (/hr) proliferation rate IPA wrt oxygen
% alpha22 = 0.0109; % (/hr) proliferation rate IPA wrt PDGFA
% beta1 = 0.0589; % (/hr) basal differentiation rate (0.08)
% beta2 = 0.0337; % (/hr) differentiation rate wrt oxygen (0.03)
% beta3 = 0.0184; % (/hr) differentiation rate wrt LIF (0.02)
% gamma1 = 0.0; % (/hr) apoptosis rate APC
% gamma2 = 0.0; % (/hr) apoptosis rate IPA

%%% from oct 4 - starting from tim's values
% alpha10 = 0.0584; % (/hr) basal proliferation rate (0.1)
% alpha11 = 0.0673; %  (/hr) proliferation rate APC wrt oxygen (0.08)
% alpha12 = 0.1761; %  (/hr) proliferation rate APC wrt PDGFA (0.1)
% alpha20 = 0.0019; % (/hr) basal proliferation rate (zero) 
% alpha21 = 0.0019; % (/hr) proliferation rate IPA wrt oxygen
% alpha22 = 0.0030; % (/hr) proliferation rate IPA wrt PDGFA
% beta1 = 0.0712; % (/hr) basal differentiation rate (0.08)
% beta2 = 0.0049; % (/hr) differentiation rate wrt oxygen (0.03)
% beta3 = 0.0078; % (/hr) differentiation rate wrt LIF (0.02)
% gamma1 = 0.0; % (/hr) apoptosis rate APC
% gamma2 = 0.0; % (/hr) apoptosis rate IPA

%%% from oct 4 - starting from my values
% alpha10 = 0.0774; % (/hr) basal proliferation rate (0.1)
% alpha11 = 0.0262; %  (/hr) proliferation rate APC wrt oxygen (0.08)
% alpha12 = 0.1286; %  (/hr) proliferation rate APC wrt PDGFA (0.1)
% alpha20 = 0.00039; % (/hr) basal proliferation rate (zero) 
% alpha21 = 0.0053; % (/hr) proliferation rate IPA wrt oxygen
% alpha22 = 0.0220; % (/hr) proliferation rate IPA wrt PDGFA
% beta1 = 0.0547; % (/hr) basal differentiation rate (0.08)
% beta2 = 0.0171; % (/hr) differentiation rate wrt oxygen (0.03)
% beta3 = 0.0060; % (/hr) differentiation rate wrt LIF (0.02)
% gamma1 = 0.0; % (/hr) apoptosis rate APC
% gamma2 = 0.0; % (/hr) apoptosis rate IPA

%%% from oct 5 - starting from tim's values
% alpha10 = 0.0820; % (/hr) basal proliferation rate (0.1)
% alpha11 = 0.0681; %  (/hr) proliferation rate APC wrt oxygen (0.08)
% alpha12 = 0.1976; %  (/hr) proliferation rate APC wrt PDGFA (0.1)
% alpha20 = 0.0002; % (/hr) basal proliferation rate (zero) 
% alpha21 = 0.0033; % (/hr) proliferation rate IPA wrt oxygen
% alpha22 = 0.0088; % (/hr) proliferation rate IPA wrt PDGFA
% beta1 = 0.0771; % (/hr) basal differentiation rate (0.08)
% beta2 = 0.0338; % (/hr) differentiation rate wrt oxygen (0.03)
% beta3 = 0.0097; % (/hr) differentiation rate wrt LIF (0.02)
% gamma1 = 0.0; % (/hr) apoptosis rate APC
% gamma2 = 0.0; % (/hr) apoptosis rate IPA

%%% from oct 5 - starting from my values
% alpha10 = 0.0734; % (/hr) basal proliferation rate (0.1)
% alpha11 = 0.0701; %  (/hr) proliferation rate APC wrt oxygen (0.08)
% alpha12 = 0.1835; %  (/hr) proliferation rate APC wrt PDGFA (0.1)
% alpha20 = 0.0019; % (/hr) basal proliferation rate (zero) 
% alpha21 = 0.0018; % (/hr) proliferation rate IPA wrt oxygen
% alpha22 = 0.0016; % (/hr) proliferation rate IPA wrt PDGFA
% beta1 = 0.0727; % (/hr) basal differentiation rate (0.08)
% beta2 = 0.0336; % (/hr) differentiation rate wrt oxygen (0.03)
% beta3 = 0.00036; % (/hr) differentiation rate wrt LIF (0.02)
% gamma1 = 0.0; % (/hr) apoptosis rate APC
% gamma2 = 0.0; % (/hr) apoptosis rate IPA

%%% from oct 5 - starting from tim's values - radius weighted points 5-6
% alpha10 = 0.0004; % (/hr) basal proliferation rate (0.1)
% alpha11 = 0.0859; %  (/hr) proliferation rate APC wrt oxygen (0.08)
% alpha12 = 0.1699; %  (/hr) proliferation rate APC wrt PDGFA (0.1)
% alpha20 = 0.0024; % (/hr) basal proliferation rate (zero) 
% alpha21 = 0.0028; % (/hr) proliferation rate IPA wrt oxygen
% alpha22 = 0.0035; % (/hr) proliferation rate IPA wrt PDGFA
% beta1 = 0.0105; % (/hr) basal differentiation rate (0.08)
% beta2 = 0.0308; % (/hr) differentiation rate wrt oxygen (0.03)
% beta3 = 0.00004; % (/hr) differentiation rate wrt LIF (0.02)
% gamma1 = 0.0; % (/hr) apoptosis rate APC
% gamma2 = 0.0; % (/hr) apoptosis rate IPA

%%% from oct 5 - starting from my values - radius weighted points 5-6
% alpha10 = 0.0708; % (/hr) basal proliferation rate (0.1)
% alpha11 = 0.0398; %  (/hr) proliferation rate APC wrt oxygen (0.08)
% alpha12 = 0.0949; %  (/hr) proliferation rate APC wrt PDGFA (0.1)
% alpha20 = 0.0018; % (/hr) basal proliferation rate (zero) 
% alpha21 = 0.0015; % (/hr) proliferation rate IPA wrt oxygen
% alpha22 = 0.0084; % (/hr) proliferation rate IPA wrt PDGFA
% beta1 = 0.0508; % (/hr) basal differentiation rate (0.08)
% beta2 = 0.0139; % (/hr) differentiation rate wrt oxygen (0.03)
% beta3 = 0.0007; % (/hr) differentiation rate wrt LIF (0.02)
% gamma1 = 0.0; % (/hr) apoptosis rate APC
% gamma2 = 0.0; % (/hr) apoptosis rate IPA

%%% from nov 2 - starting from tim's values - radius weighted points 4-6
% alpha10 = 0.0794; % (/hr) basal proliferation rate (0.1)
% alpha11 = 0.0343; %  (/hr) proliferation rate APC wrt oxygen (0.08)
% alpha12 = 0.2187; %  (/hr) proliferation rate APC wrt PDGFA (0.1)
% alpha20 = 0.001; % (/hr) basal proliferation rate (zero) 
% alpha21 = 0.0005; % (/hr) proliferation rate IPA wrt oxygen
% alpha22 = 0.0007; % (/hr) proliferation rate IPA wrt PDGFA
% beta1 = 0.0483; % (/hr) basal differentiation rate (0.08)
% beta2 = 0.0331; % (/hr) differentiation rate wrt oxygen (0.03)
% beta3 = 0.0184; % (/hr) differentiation rate wrt LIF (0.02)
% gamma1 = 0.0; % (/hr) apoptosis rate APC
% gamma2 = 0.0; % (/hr) apoptosis rate IPA

%%% from nov 2 - starting from my values - radius weighted points 4-6
% alpha10 = 0.0973; % (/hr) basal proliferation rate (0.1)
% alpha11 = 0.00002; %  (/hr) proliferation rate APC wrt oxygen (0.08)
% alpha12 = 0.1288; %  (/hr) proliferation rate APC wrt PDGFA (0.1)
% alpha20 = 0.0013; % (/hr) basal proliferation rate (zero) 
% alpha21 = 0.0058; % (/hr) proliferation rate IPA wrt oxygen
% alpha22 = 0.0141; % (/hr) proliferation rate IPA wrt PDGFA
% beta1 = 0.0434; % (/hr) basal differentiation rate (0.08)
% beta2 = 0.0210; % (/hr) differentiation rate wrt oxygen (0.03)
% beta3 = 0.0250; % (/hr) differentiation rate wrt LIF (0.02)
% gamma1 = 0.0; % (/hr) apoptosis rate APC
% gamma2 = 0.0; % (/hr) apoptosis rate IPA

%%% from nov 3 - starting from tim's values - radius weighted points 4-7
% alpha10 = 0.0751; % (/hr) basal proliferation rate (0.1)
% alpha11 = 0.0581; %  (/hr) proliferation rate APC wrt oxygen (0.08)
% alpha12 = 0.1550; %  (/hr) proliferation rate APC wrt PDGFA (0.1)
% alpha20 = 0.0003; % (/hr) basal proliferation rate (zero) 
% alpha21 = 0.0006; % (/hr) proliferation rate IPA wrt oxygen
% alpha22 = 0.0007; % (/hr) proliferation rate IPA wrt PDGFA
% beta1 = 0.0481; % (/hr) basal differentiation rate (0.08)
% beta2 = 0.0518; % (/hr) differentiation rate wrt oxygen (0.03)
% beta3 = 0.0118; % (/hr) differentiation rate wrt LIF (0.02)
% gamma1 = 0.0; % (/hr) apoptosis rate APC
% gamma2 = 0.0; % (/hr) apoptosis rate IPA

%%% from nov 3 - starting from my values - radius weighted points 4-7
% alpha10 = 0.0847; % (/hr) basal proliferation rate (0.1)
% alpha11 = 0.0562; %  (/hr) proliferation rate APC wrt oxygen (0.08)
% alpha12 = 0.1987; %  (/hr) proliferation rate APC wrt PDGFA (0.1)
% alpha20 = 0.0008; % (/hr) basal proliferation rate (zero) 
% alpha21 = 0.0048; % (/hr) proliferation rate IPA wrt oxygen
% alpha22 = 0.0006; % (/hr) proliferation rate IPA wrt PDGFA
% beta1 = 0.0680; % (/hr) basal differentiation rate (0.08)
% beta2 = 0.0364; % (/hr) differentiation rate wrt oxygen (0.03)
% beta3 = 0.0189; % (/hr) differentiation rate wrt LIF (0.02)
% gamma1 = 0.0; % (/hr) apoptosis rate APC
% gamma2 = 0.0; % (/hr) apoptosis rate IPA

%%% from nov 3 - starting from tim's values - radius weighted points 3,5-7
% alpha10 = 0.0625; % (/hr) basal proliferation rate (0.1)
% alpha11 = 0.0428; %  (/hr) proliferation rate APC wrt oxygen (0.08)
% alpha12 = 0.1411; %  (/hr) proliferation rate APC wrt PDGFA (0.1)
% alpha20 = 0.0013; % (/hr) basal proliferation rate (zero) 
% alpha21 = 0.0011; % (/hr) proliferation rate IPA wrt oxygen
% alpha22 = 0.0014; % (/hr) proliferation rate IPA wrt PDGFA
% beta1 = 0.0423; % (/hr) basal differentiation rate (0.08)
% beta2 = 0.0294; % (/hr) differentiation rate wrt oxygen (0.03)
% beta3 = 0.0087; % (/hr) differentiation rate wrt LIF (0.02)
% gamma1 = 0.0; % (/hr) apoptosis rate APC
% gamma2 = 0.0; % (/hr) apoptosis rate IPA

%%% from nov 3 - starting from my values - radius weighted points 3,5-7
% alpha10 = 0.0798; % (/hr) basal proliferation rate (0.1)
% alpha11 = 0.0795; %  (/hr) proliferation rate APC wrt oxygen (0.08)
% alpha12 = 0.1258; %  (/hr) proliferation rate APC wrt PDGFA (0.1)
% alpha20 = 0.0000009; % (/hr) basal proliferation rate (zero) 
% alpha21 = 0.000005; % (/hr) proliferation rate IPA wrt oxygen
% alpha22 = 0.0155; % (/hr) proliferation rate IPA wrt PDGFA
% beta1 = 0.0698; % (/hr) basal differentiation rate (0.08)
% beta2 = 0.0339; % (/hr) differentiation rate wrt oxygen (0.03)
% beta3 = 0.0191; % (/hr) differentiation rate wrt LIF (0.02)
% gamma1 = 0.0; % (/hr) apoptosis rate APC
% gamma2 = 0.0; % (/hr) apoptosis rate IPA

%% tension parameters
mu1 = 15; %15 % adhesion constant (mm/hr/(mN/mm^2))
mu2 = 15; %15
kTprime1 = kappa / (2 * mu1 * sqrt(pi));
kTprime2 = kappa / (2 * mu2 * sqrt(pi));

% parameter in initial condition, chosen to match subsequent solution
ce1 = ce*(1 + 0.5 * alpha11 * s0^2 * sqrt(ce) * (1 - ce/cmax) / kTprime1);

%% set up and solve
% mesh
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
%options = odeset('RelTol',1e-10,'AbsTol',1e-10,'NormControl',1);
sol = pdepe(1,@AstroPDE,@AstroIC,@AstroBC,x,t);%,options);

%% plot results

% plot_the_plots
% plot_the_plots_5panels
plot_the_plots_6panels
% plot_the_plots_APCIPA
% plot_growthfactors_O2_thickness

% plot_timversion

%% calculate error
[err_tot,err_rad,err_dens] = errorfunction(t,x,sol)