clear variables global;
clc;

%%% this script file is to run a sensitivity analysis of the parameters
diary results_sensitivity_analysis.txt

%%% time unit: hr
%%% space unit: mm

m.dr = 0.01;
m.rmax = 5; %%% max radius (mm) (estimate rat retinal radius = 4.1 mm)
m.tmax = 7*24; %%% max time (hr) (7 days = 168 hr)

%%%%%%%%%%%%%%%%%%%%%%%%%%% baseline parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%
base_mu = 0.1; %%% adhesion constant
base_alpha11 = 0.15; %%% (/hr) proliferation rate APC wrt oxygen
base_alpha12 = 0; %%% (/hr) proliferation rate APC wrt PDGFA
base_alpha21 = 0.15; %%% (/hr) proliferation rate IPA wrt oxygen
base_alpha22 = 0; %%% (/hr) proliferation rate IPA wrt PDGFA
base_beta = 0.01;%0.003; %%% (/hr) differentiation rate
base_beta_hat = 0.005; %%% (/hr) mass action rate
base_gamma1 = 0.0001; %%% (/hr) apoptosis rate APC
base_gamma2 = 0.0001; %%% (/hr) apoptosis rate IPA
base_Te = 0.0035; %%% tension on boundary
base_P_hy = 0.01; %%% partial pressure of oxygen due to hyaloid artery
base_r_hy = 0.1; %%% radius at half-maximum of Hill function for hyaloid

%%%%%%%%%%%%%%%%%%%%%%%%%%%% parameters bounds %%%%%%%%%%%%%%%%%%%%%%%%%%%%
bound = [0.01 5; %mu
    0 1; %alpha11
    0 1; %alpha12
    0 1; %alpha21
    0 1; %alpha22
    0 0.15; %beta
    0 0.5; %beta_hat
    0 0.005; %gamma1
    0 0.005; %gamma2
    0 0.0038; %Te
    0 1; %P_hy
    1/1000 1]; %r_hy
numpar = length(bound);

N = 3;

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
        alpha21 = base_alpha21;
        alpha22 = base_alpha22;
        beta = base_beta;
        beta_hat = base_beta_hat;
        gamma1 = base_gamma1;
        gamma2 = base_gamma2;
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
            alpha21 = intrange(j,i);
        elseif j==5
            alpha22 = intrange(j,i);
        elseif j==6
            beta = intrange(j,i);
        elseif j==7
            beta_hat = intrange(j,i);
        elseif j==8
            gamma1 = intrange(j,i);
        elseif j==9
            gamma2 = intrange(j,i);
        elseif j==10
            Te = intrange(j,i);
        elseif j==11
            P_hy = intrange(j,i);
        elseif j==12
            r_hy = intrange(j,i);
        end
        
        %%% solve equation
        tic
        [t,r,c1,c2,q1,q2,mvgbdy,vel_cir,vel_rad] = eqnsolver(mu,alpha11,...
            alpha12,alpha21,alpha22,beta,beta_hat,gamma1,gamma2,Te,P_hy,r_hy,m);
        toc

        %%% error calculation
        [err_rad,err_dens,err_tot] = errorfunction(t,r,mvgbdy,c1,c2);
        
        err(j,i) = err_tot;

        save('results_sensitivity_analysis.mat');
    end
end

save('results_sensitivity_analysis.mat');

diary off