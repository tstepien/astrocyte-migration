clear variables global;
clc;
addpath plot_simulations

%%%%%%%%%%%%%%%%%%% input all fixed parameters that are %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% known/derived from literature %%%%%%%%%%%%%%%%%%%%%%
parameters_fixed

%%%%%%%%%%%%%%%%%%%%% parameter scalings/calculations %%%%%%%%%%%%%%%%%%%%%
xi1 = xibar_PDGFA / phi; %%% production/release rate of PDGFA
xi2 = xibar_LIF / phi; %%% production/release rate of LIF

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% mesh %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dr = 0.1;
dt = 0.01;

%check CFL
if dt >=dr^2/(2*max(D1,D2))
    disp('check CFL')
end

rmax = 5; %%% max radius (mm) (estimate rat retinal radius = 4.1 mm)
tmax = 7*24; %%% max time (hr) (7 days = 168 hr)

r = 0:dr:rmax;
R = length(r);


%%%%%%%%%%%%%%%%%%%%%%%%%%% initial conditions %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% growth factors initial condition
q1_old = zeros(1,R);
q2_old = zeros(1,R);

%%% initialize variables
t = 0:dt:tmax;

q1_exp = zeros(length(t),length(r));
q2_exp = zeros(length(t),length(r));
q1_imp = zeros(length(t),length(r));
q2_imp = zeros(length(t),length(r));

radius_ret = zeros(length(t),1);
radius_endo = zeros(length(t),1);
[~,~,radius_endo(1),radius_ret(1)] = thick_rad(0,r);

q1_exp(1,:) = q1_old;
q2_exp(1,:) = q2_old;
q1_imp(1,:) = q1_old;
q2_imp(1,:) = q2_old;

tcurr = 0;

for i = 2:length(t)
    %%% cell layer thickness and radius
    [thickness_ret,thickness_RGC,radius_endo(i),radius_ret(i)] = thick_rad(tcurr+dt,r);
    
    %%% growth factors
%     [q1_exp(i,:),q2_exp(i,:)] = growthfactors_explicit(q1_exp(i-1,:),...
%         q2_exp(i-1,:),dt,r,D1,D2,xi1,xi2);
    [q1_imp(i,:),q2_imp(i,:)] = growthfactors_implicit(q1_imp(i-1,:),...
        q2_imp(i-1,:),dt,tcurr,r,dr,R,thickness_RGC,radius_endo(i),...
        maxRGCthick,thickness_ret,D1,D2,xi1,xi2,gamma3,gamma4);
    
    tcurr = tcurr + dt;
end

diffq1 = abs(q1_exp-q1_imp);
diffq2 = abs(q2_exp-q2_imp);



%% plots
plot_growthfactorsonly