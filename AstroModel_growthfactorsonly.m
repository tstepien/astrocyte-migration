%   AstrocyteModel_TWS
%   Model for spread of astrocytes of surface of retina
%   TWS and TLS, 2021
%
%   ***Simulation of growth factors only***
%
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

clear variables global;
addpath plot_simulations
global rmax LIF PDGFA nxpts ntpts tmax

parameters_fixed;

%% configuration parameters
rmax = 5; % max radius (mm) (estimate rat retinal radius = 4.1 mm)
tmax = 7 * 24; % max time (hr) (7 days)

%% set up and solve
nxpts = 100;
ntpts = 29;
x = linspace(0,1,nxpts);    % 0 <= x <= 1 by definition
t = linspace(0,tmax,ntpts);

% Calculate growth factors for interpolation using fixed domain
r = linspace(0,rmax,nxpts);
sol = pdepe(1,@GF_PDE,@GF_IC,@GF_BC,r,t);
PDGFA = sol(:,:,1);
LIF = sol(:,:,2);

%% plot results

plot_growthfactorsonly