clear variables global;
clc;

%%% this script file is to run a sensitivity analysis of the parameters
diary sensitivity_analysis.txt

%%% time unit: hr
%%% space unit: mm

m.dr = 0.01;
m.rmax = 5; %%% max radius (mm) (estimate rat retinal radius = 4.1 mm)
m.tmax = 7*24; %%% max time (hr) (7 days = 168 hr)
plotsonoff = 'off';

%%%%%%%%%%%%%%%%%%%%%%%%%%% baseline parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%
p.mu = 0.5; %%% adhesion constant
p.alpha1 = 0.15; %%% (/hr) proliferation rate APC
p.alpha2 = 0.15; %%% (/hr) proliferation rate IPA
p.beta = 0.02;%0.003; %%% (/hr) differentiation rate
p.beta_hat = 0.05; %%% (/hr) mass action rate
p.gamma1 = 0.0001; %%% (/hr) apoptosis rate APC
p.gamma2 = 0.0001; %%% (/hr) apoptosis rate IPA
p.Te = 0.0035; %%% tension on boundary
p.P_hy = 0.2; %%% partial pressure of oxygen due to hyaloid artery
p.r_hy = 1/10; %%% radius at half-maximum of Hill function for hyaloid

%%%%%%%%%%%%%%%%%%%%%%%%%%%% parameters bounds %%%%%%%%%%%%%%%%%%%%%%%%%%%%
bound = [0.01 5; %mu
    0 1; %alpha1
    0 1; %alpha2
    0 0.15; %beta
    0 0.5; %beta_hat
    0 0.005; %gamma1
    0 0.005; %gamma2
    0 0.0038; %Te
    0 1; %P_hy
    1/1000 1]; %r_hy
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
        paramval = p;
        if j==1
            paramval.mu = intrange(j,i);
        elseif j==2
            paramval.alpha1 = intrange(j,i);
        elseif j==3
            paramval.alpha2 = intrange(j,i);
        elseif j==4
            paramval.beta = intrange(j,i);
        elseif j==5
            paramval.beta_hat = intrange(j,i);
        elseif j==6
            paramval.gamma1 = intrange(j,i);
        elseif j==7
            paramval.gamma2 = intrange(j,i);
        elseif j==8
            paramval.Te = intrange(j,i);
        elseif j==9
            paramval.P_hy = intrange(j,i);
        elseif j==10
            paramval.r_hy = intrange(j,i);
        end
        
        %%% solve equation
        tic
        [t,r,c1,c2,q1,q2,mvgbdy,vel_cir,vel_rad] = eqnsolver(paramval,m,plotsonoff);
        toc

        %%% error calculation
        [err_rad,err_dens,err_tot] = errorfunction(t,r,mvgbdy,c1,c2);
        
        err(j,i) = err_tot;

        save('sensitivity_analysis.mat','p','intrange','err','bound');
    end
end

save('sensitivity_analysis.mat','p','intrange','err','bound')

diary off