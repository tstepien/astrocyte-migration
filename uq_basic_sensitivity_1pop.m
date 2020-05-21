clear variables global;
clc;

%%% this script file is to run a sensitivity analysis of the parameters
diary results_sensitivity_analysis.txt

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

%%%%%%%%%%%%%%%%%%%%%%%%%%% baseline parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%
base_mu = 0.1; %%% adhesion constant
base_alpha11 = 0.03; %%% (/hr) proliferation rate APC wrt oxygen
base_alpha12 = 0.22; %%% (/hr) proliferation rate APC wrt PDGFA
base_gamma1 = 0.01; %%% (/hr) apoptosis rate APC
base_Te = 0.0035; %%% tension on boundary
base_P_hy = 10; %%% partial pressure of oxygen due to hyaloid artery
base_r_hy = 0.1; %%% radius at half-maximum of Hill function for hyaloid

%%%%%%%%%%%%%%%%%%%%%%%%%%%% parameters bounds %%%%%%%%%%%%%%%%%%%%%%%%%%%%
bound = [0.01 5; %mu
    0.01 1; %alpha11
    0 1; %alpha12
    0 1; %gamma1
    0.0001 0.0038; %Te
    0 20; %P_hy
    0.001 1]; %r_hy
numpar = length(bound);

N = 21;

intrange = zeros(numpar,N);
for i=1:numpar
    intrange(i,:) = linspace(bound(i,1),bound(i,2),N);
end

%%% preallocate
err = zeros(size(intrange));

for j = 1:numpar %%% parameter
    for i = 1:N %%% split up interval range
        [j i]
        
        mu = base_mu;
        alpha11 = base_alpha11;
        alpha12 = base_alpha12;
        gamma1 = base_gamma1;
        Te = base_Te;
        P_hy = base_P_hy;
        r_hy = base_r_hy;
            
        if j==1
            mu = intrange(j,i);
        elseif j==2
            alpha11 = intrange(j,i);
        elseif j==3
            alpha12 = intrange(j,i);
        elseif j==4
            gamma1 = intrange(j,i);
        elseif j==5
            Te = intrange(j,i);
        elseif j==6
            P_hy = intrange(j,i);
        elseif j==7
            r_hy = intrange(j,i);
        end
        
        %%% solve equation
        tic
        [t,~,~,~,~,~,mvgbdy,~,~] = eqnsolver(mu,alpha11,alpha12,...
            alpha21,alpha22,beta,beta_hat,gamma1,gamma2,Te,P_hy,r_hy,m);
        toc

        %%% error calculation
        err_rad = errorfunction_1pop(t,mvgbdy);
        
        err(j,i) = err_rad;

        save('results_sensitivity_analysis.mat');
    end
end

save('results_sensitivity_analysis.mat');

diary off