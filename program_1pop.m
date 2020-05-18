clear variables global;
clc;
close all;

%%% time unit: hr
%%% space unit: mm

%%%%%%%%%%%%%%%%%%%%%%%%%% astrocyte parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% only APCs
p.mu = 0.1; %%% adhesion constant
p.alpha11 = 0.03; %%% (/hr) proliferation rate APC wrt oxygen
p.alpha12 = 0.22; %%% (/hr) proliferation rate APC wrt PDGFA
p.gamma1 = 0.01; %%% (/hr) apoptosis rate APC

%%% set all IPA stuff = 0
p.alpha21 = 0; %%% (/hr) proliferation rate IPA wrt oxygen
p.alpha22 = 0; %%% (/hr) proliferation rate IPA wrt PDGFA
p.beta = 0; %%% (/hr) differentiation rate
p.beta_hat = 0; %%% (/hr) mass action rate
p.gamma2 = 0; %%% (/hr) apoptosis rate IPA

%%%%%%%%%%%%%%%%%% moving boundary condition parameters %%%%%%%%%%%%%%%%%%%
p.Te = 0.0035; %%% tension on boundary

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% hyaloid artery %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p.P_hy = 0; %%% partial pressure of oxygen due to hyaloid artery
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%% error calculation %%%%%%%%%%%%%%%%%%%%%%%%%%%%
err_rad = errorfunction_1pop(t,mvgbdy);