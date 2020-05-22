clear variables global;
clc;

savefiles = 'yes';

N = 10000;

% LHpts = lhsdesign(N,7);
% save(strcat('LHpts_',num2str(N),'.mat'),'LHpts')
% stop

load(strcat('parameter analysis/LHpts_',num2str(N),'.mat'))

if strcmp(savefiles,'yes')==1
    doublecheck = input('Are you sure you would like to save the output files? (it may overwrite): ');
    if strcmp(doublecheck,'y')==1
        diary(strcat('parameter analysis/latinhypercube_',num2str(N),'.txt'));
    else
        return;
    end
end

%%% mesh
m.dr = 0.01;
m.rmax = 5; %%% max radius (mm) (estimate rat retinal radius = 4.1 mm)
m.tmax = 7*24; %%% max time (hr) (7 days = 168 hr)

%%% set all IPA stuff = 0
alpha21 = 0; %%% (/hr) proliferation rate IPA wrt oxygen
alpha22 = 0; %%% (/hr) proliferation rate IPA wrt PDGFA
beta = 0; %%% (/hr) differentiation rate
beta_hat = 0; %%% (/hr) mass action rate
gamma2 = 0; %%% (/hr) apoptosis rate IPA

%%% parameters to investigate
bound = [0.01 5; %mu - adhesion constant
    0.01 1; %alpha11 - (/hr) proliferation rate APC wrt oxygen
    0 1; %alpha12 - (/hr) proliferation rate APC wrt PDGFA
    0 1; %gamma1 - (/hr) apoptosis rate APC
    0.0001 0.0038; %Te - tension on boundary
    0 20; %P_hy - partial pressure of oxygen due to hyaloid artery
    0.001 1]; %r_hy - radius at half-maximum of Hill function for hyaloid
numpar = length(bound);

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

parfor i=1745:N
    %%% solve equation
    [t,~,~,~,~,~,mvgbdy,~,~] = eqnsolver(mu(i),alpha11(i),alpha12(i),...
        alpha21,alpha22,beta,beta_hat,gamma1(i),gamma2,Te(i),P_hy(i),...
        r_hy(i),m);
    
    %%% error calculation
    [err_tot(i),err_time(i),err_rad(i)]  = errorfunction_1pop(t,mvgbdy);
end

if strcmp(savefiles,'yes')==1
    save(strcat('parameter analysis/latinhypercube_',num2str(N),'pts.mat'));
    diary off
end