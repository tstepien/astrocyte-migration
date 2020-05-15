clear variables global;
clc;
close all;

%%% time unit: hr
%%% space unit: mm

%%%%%%%%%%%%%%%%%%%%%%%%%% astrocyte parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%
p.mu = 0.1; %%% adhesion constant
p.alpha11 = 0.111; %%% (/hr) proliferation rate APC wrt oxygen
p.alpha12 = 0.22; %%% (/hr) proliferation rate APC wrt PDGFA
p.alpha21 = 0.15; %%% (/hr) proliferation rate IPA wrt oxygen
p.alpha22 = 0.15; %%% (/hr) proliferation rate IPA wrt PDGFA
p.beta = 0.01; %%% (/hr) differentiation rate
p.beta_hat = 0.005; %%% (/hr) mass action rate
p.gamma1 = 0.0001; %%% (/hr) apoptosis rate APC
p.gamma2 = 0.0001; %%% (/hr) apoptosis rate IPA

%%%%%%%%%%%%%%%%%% moving boundary condition parameters %%%%%%%%%%%%%%%%%%%
p.Te = 0.0035; %%% tension on boundary

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% hyaloid artery %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p.P_hy = 0.01; %%% partial pressure of oxygen due to hyaloid artery
p.r_hy = 0.1; %%% radius at half-maximum of Hill function for hyaloid

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% mesh parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
m.dr = 0.01;
m.rmax = 5; %%% max radius (mm) (estimate rat retinal radius = 4.1 mm)
m.tmax = 7*24; %%% max time (hr) (7 days = 168 hr)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% plots on or off %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plotsonoff = 'on';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% solve equation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
[t,r,c1,c2,q1,q2,mvgbdy,vel_cir,vel_rad] = eqnsolver(p,m,plotsonoff);
toc


%% error
%%%%%%%%%%%%%%%%%%%%%%%%%%%% error calculation %%%%%%%%%%%%%%%%%%%%%%%%%%%%
[err_rad,err_dens,err_tot] = errorfunction(t,r,mvgbdy,c1,c2);


%% area under curve - trapezoidal rule
%%%%%%%%%%%%%%%%%%%%%% testing: conservation of mass %%%%%%%%%%%%%%%%%%%%%%
areaundercurve = zeros(length(t),1);
dr = m.dr;
k = c1+c2;
for i = 1:length(t)
    areaundercurve(i) = dr*sum( k(1:mvgbdy(1)/dr+(i-1)) ...
        + k(2:mvgbdy(1)/dr+1+(i-1)) ) / 2;
end