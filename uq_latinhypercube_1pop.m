clear variables global;
clc;

%%% this script file is to run a sensitivity analysis of the parameters
diary results_latinhypercube.txt

%%% time unit: hr
%%% space unit: mm

m.dr = 0.01;
m.rmax = 5; %%% max radius (mm) (estimate rat retinal radius = 4.1 mm)
m.tmax = 7*24; %%% max time (hr) (7 days = 168 hr)

%%% set all IPA stuff = 0
alpha21 = 0; %%% (/hr) proliferation rate IPA wrt oxygen
alpha22 = 0; %%% (/hr) proliferation rate IPA wrt PDGFA
beta = 0; %%% (/hr) differentiation rate
beta_hat = 0; %%% (/hr) mass action rate
gamma2 = 0; %%% (/hr) apoptosis rate IPA

%%%%%%%%%%%%%%%%%%%%%%%%%%%% parameters bounds %%%%%%%%%%%%%%%%%%%%%%%%%%%%
bound = [0.01 5; %mu - adhesion constant
    0.01 1; %alpha11 - (/hr) proliferation rate APC wrt oxygen
    0 1; %alpha12 - (/hr) proliferation rate APC wrt PDGFA
    0 1; %gamma1 - (/hr) apoptosis rate APC
    0.0001 0.0038; %Te - tension on boundary
    0 20; %P_hy - partial pressure of oxygen due to hyaloid artery
    0.001 1]; %r_hy - radius at half-maximum of Hill function for hyaloid
numpar = length(bound);

N = 5000;
LHpts = lhsdesign(N,numpar);

mu = (bound(1,2) - bound(1,1))*LHpts(:,1) + bound(1,1);
alpha11 = (bound(2,2) - bound(2,1))*LHpts(:,2) + bound(2,1);
alpha12 = (bound(3,2) - bound(3,1))*LHpts(:,3) + bound(3,1);
gamma1 = (bound(4,2) - bound(4,1))*LHpts(:,4) + bound(4,1);
Te = (bound(5,2) - bound(5,1))*LHpts(:,5) + bound(5,1);
P_hy = (bound(6,2) - bound(6,1))*LHpts(:,6) + bound(6,1);
r_hy = (bound(7,2) - bound(7,1))*LHpts(:,7) + bound(7,1);

%%% preallocate
err_tot = zeros(N,1);
err_time = zeros(N,1);
err_rad = zeros(N,1);

parfor i=1:N    
    %%% solve equation
    tic
    [t,~,~,~,~,~,mvgbdy,~,~] = eqnsolver(mu(i),alpha11(i),alpha12(i),...
        alpha21,alpha22,beta,beta_hat,gamma1(i),gamma2,Te(i),P_hy(i),...
        r_hy(i),m);
    toc
    
    %%% error calculation
    [err1,err2,err3]  = errorfunction_1pop(t,mvgbdy);
    
    err_tot(i) = err1;
    err_time(i) = err2;
    err_rad(i) = err3;
end

% correlation_param = zeros(numpar,1);
% correlation_param(1) = corr(mu,err_time);
% correlation_param(2) = corr(alpha11,err_time);
% correlation_param(3) = corr(alpha12,err_time);
% correlation_param(4) = corr(gamma1,err_time);
% correlation_param(5) = corr(Te,err_time);
% correlation_param(6) = corr(P_hy,err_time);
% correlation_param(7) = corr(r_hy,err_time);

param_bad = [];
param_good = [];

for i=1:N
    if err_time(i)>=0.7
        param_bad = [param_bad ; 
            mu(i) alpha11(i) alpha12(i) gamma1(i) Te(i) P_hy(i) r_hy(i)];
    else
        param_good = [param_good ; 
            mu(i) alpha11(i) alpha12(i) gamma1(i) Te(i) P_hy(i) r_hy(i)];
    end
end

new_range = zeros(numpar,2);
for i = 1:numpar
    new_range(i,:) = [min(param_good(:,i)) max(param_good(:,i))];
end

save('results_latinhypercube.mat');

diary off