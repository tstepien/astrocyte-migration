clear variables global;
clc;
close all;
addpath plot_simulations

%%% time unit: hr
%%% space unit: mm

%%%%%%%%%%%%%%%%%%%%%%%%%% astrocyte parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% only APCs
mu = 2.5017;%0.1; %%% adhesion constant
alpha11 = 0.4456;%0.03; %%% (/hr) proliferation rate APC wrt oxygen
alpha12 = 0.0299;%0.22; %%% (/hr) proliferation rate APC wrt PDGFA
gamma1 = 0.8259;%0.01; %%% (/hr) apoptosis rate APC

%%% set all IPA stuff = 0
alpha21 = 0; %%% (/hr) proliferation rate IPA wrt oxygen
alpha22 = 0; %%% (/hr) proliferation rate IPA wrt PDGFA
beta = 0; %%% (/hr) differentiation rate
beta_hat = 0; %%% (/hr) mass action rate
gamma2 = 0; %%% (/hr) apoptosis rate IPA

%%%%%%%%%%%%%%%%%% moving boundary condition parameters %%%%%%%%%%%%%%%%%%%
Te = 0.0035; %%% tension on boundary

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% hyaloid artery %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
P_hy = 19.2786;%0; %%% partial pressure of oxygen due to hyaloid artery
r_hy = 0.5649;%0.1; %%% radius at half-maximum of Hill function for hyaloid

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% mesh parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
m.dr = 0.01;
m.rmax = 5; %%% max radius (mm) (estimate rat retinal radius = 4.1 mm)
m.tmax = 7*24; %%% max time (hr) (7 days = 168 hr)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% solve equation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
[t,r,c1,c2,q1,q2,mvgbdy,vel_cir,vel_rad] = eqnsolver(mu,alpha11,alpha12,...
    alpha21,alpha22,beta,beta_hat,gamma1,gamma2,Te,P_hy,r_hy,m);
toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%% error calculation %%%%%%%%%%%%%%%%%%%%%%%%%%%%
[err_tot,err_time,err_rad] = errorfunction_1pop(t,mvgbdy);

%% plots
plot_the_plots
plot_the_plots_APCIPA
plot_movingbdy