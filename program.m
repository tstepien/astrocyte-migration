clear variables global;
clc;
% close all;

%%% time unit: hr
%%% space unit: mm

%%%%%%%%%%%%%%%%%%%%%%%%%% astrocyte parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%
p.kappa = 1; %%% tension function scaling
p.mu = 0.1; %%% adhesion constant
p.alpha1 = 0.1975; %%% (/hr) proliferation rate APC
p.alpha2 = p.alpha1 *(7/12); %%% (/hr) proliferation rate IPA
p.beta = 0.03; %%% (/hr) differentiation rate
p.gamma1 = 0;%0.0001;%0; apoptosis rate APC
p.gamma2 = 0;%0.0001;%0.5; apoptosis rate IPA

%%%%%%%%%%%%%%%%%%%%%%%% growth factor parameters %%%%%%%%%%%%%%%%%%%%%%%%%
p.xibar_PDGFA = 0.015/15; %%% PDGFA production
p.xibar_LIF = 0.015/7; %%% LIF production

%%%%%%%%%%%%%%%%%% moving boundary condition parameters %%%%%%%%%%%%%%%%%%%
p.Te = 0.0035; %%% tension on boundary

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% mesh parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
m.dr = 0.01;
m.rmax = 5; %%% max radius (mm) (estimate rat retinal radius = 4.1 mm)
m.tmax = 7*24; %%% max time (hr) (7 days = 168 hr)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% plots on or off %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plotsonoff = 'on';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% solve equation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[t,r,c1,c2,q1,q2,mvgbdy,vel_cir,vel_rad] = eqnsolver(p,m,plotsonoff);


%% error
%%%%%%%%%%%%%%%%%%%%%%%%%%%% error calculation %%%%%%%%%%%%%%%%%%%%%%%%%%%%
[err_rad,err_dens,err_time,err_tot] = errorfunction(t,r,mvgbdy,c1,c2);


%% area under curve - trapezoidal rule
%%%%%%%%%%%%%%%%%%%%%% testing: conservation of mass %%%%%%%%%%%%%%%%%%%%%%
areaundercurve = zeros(length(t),1);
dr = m.dr;
k = c1+c2;
for i = 1:length(t)
    areaundercurve(i) = dr*sum( k(1:mvgbdy(1)/dr+(i-1)) ...
        + k(2:mvgbdy(1)/dr+1+(i-1)) ) / 2;
end