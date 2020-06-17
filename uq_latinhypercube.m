clear variables global;
clc;

savefiles = 'yes';

N = 1000000;

LHpts = lhsdesign(N,13);
save(strcat('parameter_analysis/LHpts_',num2str(N),'.mat'),'LHpts')
stop

load(strcat('parameter_analysis/LHpts_',num2str(N),'.mat'))

if strcmp(savefiles,'yes')==1
    doublecheck = input('Are you sure you would like to save the output files? (it may overwrite): ');
    if strcmp(doublecheck,'y')==1
        diary(strcat('parameter_analysis/latinhypercube_',num2str(N),'.txt'));
    else
        return;
    end
end

%%% mesh
m.dr = 0.01;
m.rmax = 5; %%% max radius (mm) (estimate rat retinal radius = 4.1 mm)
m.tmax = 7*24; %%% max time (hr) (7 days = 168 hr)

%%% parameters to investigate
bound = [0.01 5; %mu - adhesion constant
    0.01 1; %alpha11 - (/hr) proliferation rate APC wrt oxygen
    0 1; %alpha12 - (/hr) proliferation rate APC wrt PDGFA
    0 1; %alpha21 - (/hr) proliferation rate IPA wrt oxygen
    0 1; %alpha22 - (/hr) proliferation rate IPA wrt PDGFA
    0 1; %beta1 - (/hr) mass action rate
    0 1; %beta2 -  (/hr) differentiation rate wrt oxygen
    0 1; %beta3 -  (/hr) differentiation rate wrt LIF
    0 1; %gamma1 - (/hr) apoptosis rate APC
    0 1; %gamma2 - (/hr) apoptosis rate IPA
    0.0001 0.0038; %Te - tension on boundary
    0 20; %P_hy - partial pressure of oxygen due to hyaloid artery
    0.001 1]; %r_hy - radius at half-maximum of Hill function for hyaloid
numpar = length(bound);

mu      = (bound(1,2) - bound(1,1))*LHpts(:,1) + bound(1,1);
alpha11 = (bound(2,2) - bound(2,1))*LHpts(:,2) + bound(2,1);
alpha12 = (bound(3,2) - bound(3,1))*LHpts(:,3) + bound(3,1);
alpha21 = (bound(4,2) - bound(4,1))*LHpts(:,4) + bound(4,1);
alpha22 = (bound(5,2) - bound(5,1))*LHpts(:,5) + bound(5,1);
beta1   = (bound(6,2) - bound(6,1))*LHpts(:,6) + bound(6,1);
beta2   = (bound(7,2) - bound(7,1))*LHpts(:,7) + bound(7,1);
beta3   = (bound(8,2) - bound(8,1))*LHpts(:,8) + bound(8,1);
gamma1  = (bound(9,2) - bound(9,1))*LHpts(:,9) + bound(9,1);
gamma2  = (bound(10,2) - bound(10,1))*LHpts(:,10) + bound(10,1);
Te      = (bound(11,2) - bound(11,1))*LHpts(:,11) + bound(11,1);
P_hy    = (bound(12,2) - bound(12,1))*LHpts(:,12) + bound(12,1);
r_hy    = (bound(13,2) - bound(13,1))*LHpts(:,13) + bound(13,1);

%%% preallocate
err_tot = zeros(N,1);
err_time = zeros(N,1);
err_rad = zeros(N,1);
err_dens = zeros(N,1);

parfor i=1:N
    %%% solve equation    
    [t,r,c1,c2,~,~,mvgbdy,~,~] = eqnsolver(mu(i),alpha11(i),alpha12(i),...
        alpha21(i),alpha22(i),beta1(i),beta2(i),beta3(i),gamma1(i),...
        gamma2(i),Te(i),P_hy(i),r_hy(i),m);
    
    %%% error calculation
    [err_tot(i),err_time(i),err_rad(i),err_dens(i)]  = errorfunction(t,r,mvgbdy,c1,c2);
end

if strcmp(savefiles,'yes')==1
    save(strcat('parameter_analysis/latinhypercube_',num2str(N),'pts.mat'));
    diary off
end