clear variables global;
clc;


%%% time unit: hr
%%% space unit: mm

m.dr = 0.01;
m.rmax = 5; %%% max radius (mm) (estimate rat retinal radius = 4.1 mm)
m.tmax = 7*24; %%% max time (hr) (7 days = 168 hr)
plotsonoff = 'on';

%%%%%%%%%%%%%%%%%%%%%%%%%%% baseline parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%
p.mu = 0.1; %%% adhesion constant
p.alpha1 = 0.08; %%% (/hr) proliferation rate APC
p.alpha2 = 0.08; %%% (/hr) proliferation rate IPA
p.beta = 0.003; %%% (/hr) differentiation rate
p.gamma1 = 0.0001; %%% (/hr) apoptosis rate APC
p.gamma2 = 0.0001; %%% (/hr) apoptosis rate IPA
p.xibar_PDGFA = 1; %%% PDGFA production
p.xibar_LIF = 1; %%% LIF production
p.Te = 0.0035; %%% tension on boundary


%%% solve equation
tic
[t,r,c1,c2,q1,q2,mvgbdy,vel_cir,vel_rad] = eqnsolver(p,m,plotsonoff);
toc

%%% error calculation
[err_rad,err_dens,err_time,err_tot] = errorfunction(t,r,mvgbdy,c1,c2,q1,q2);