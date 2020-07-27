clear variables global;
clc;
close all;

%%% time unit: hr
%%% space unit: mm

%%%%%%%%%%%%%%%%%%%%%%%%%% astrocyte parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%
mu = [1.26090531835398];%0.1; %%% adhesion constant
alpha11 = [0.449078134936709];%0.111; %%% (/hr) proliferation rate APC wrt oxygen
alpha12 = [0.900407768950013];%0.22; %%% (/hr) proliferation rate APC wrt PDGFA
alpha21 = [0.319859742866316];%0.15; %%% (/hr) proliferation rate IPA wrt oxygen
alpha22 = [0.511353354168965];%0.15; %%% (/hr) proliferation rate IPA wrt PDGFA
beta1 = [0.00940651778996198];%0.005; %%% (/hr) mass action rate
beta2 = [0.0462454958915238];%0.01; %%% (/hr) differentiation rate wrt oxygen
beta3 = [0.590785448840143];%0.01; %%% (/hr) differentiation rate wrt LIF
gamma1 = [0.00908065628958660];%0.0001; %%% (/hr) apoptosis rate APC
gamma2 = [0.467004259082613];%0.0001; %%% (/hr) apoptosis rate IPA

%%%%%%%%%%%%%%%%%% moving boundary condition parameters %%%%%%%%%%%%%%%%%%%
Te = [0.00165058638547420];%0.0035; %%% tension on boundary

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% hyaloid artery %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
P_hy = [4.63425332322229];%0.01; %%% partial pressure of oxygen due to hyaloid artery
r_hy = [0.696259214684239];%0.1; %%% radius at half-maximum of Hill function for hyaloid

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% mesh parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
m.dr = 0.01;
m.rmax = 5; %%% max radius (mm) (estimate rat retinal radius = 4.1 mm)
m.tmax = 7*24; %%% max time (hr) (7 days = 168 hr)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% solve equation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
[t,r,c1,c2,q1,q2,mvgbdy,vel_cir,vel_rad] = eqnsolver(mu,alpha11,alpha12,...
    alpha21,alpha22,beta1,beta2,beta3,gamma1,gamma2,Te,P_hy,r_hy,m);
toc


%% error
%%%%%%%%%%%%%%%%%%%%%%%%%%%% error calculation %%%%%%%%%%%%%%%%%%%%%%%%%%%%
[err_tot,err_time,err_rad,err_dens] = errorfunction(t,r,mvgbdy,c1,c2);


%% area under curve - trapezoidal rule
%%%%%%%%%%%%%%%%%%%%%% testing: conservation of mass %%%%%%%%%%%%%%%%%%%%%%
areaundercurve = zeros(length(t),1);
dr = m.dr;
k = c1+c2;
for i = 1:length(t)
    areaundercurve(i) = dr*sum( k(1:mvgbdy(1)/dr+(i-1)) ...
        + k(2:mvgbdy(1)/dr+1+(i-1)) ) / 2;
end

%% plots
plot_the_plots
plot_the_plots_APCIPA
plot_movingbdy