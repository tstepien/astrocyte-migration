clear variables global;
clc;

%%% this script file is to run a sensitivity analysis of the parameters

%%% time unit: hr
%%% space unit: mm

m.dr = 0.01;
m.rmax = 5; %%% max radius (mm) (estimate rat retinal radius = 4.1 mm)
m.tmax = 7*24; %%% max time (hr) (7 days = 168 hr)
plotsonoff = 'off';

%%%%%%%%%%%%%%%%%%%%%%%%%%% baseline parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%
p.kappa = 1; %%% tension function scaling
p.mu = 0.1; %%% adhesion constant
p.alpha1 = 0.23; %%% (/hr) proliferation rate APC
p.alpha2 = 0.13; %%% (/hr) proliferation rate IPA
p.beta = 0.03; %%% (/hr) differentiation rate
p.gamma1 = 0.0001; %%% (/hr) apoptosis rate APC
p.gamma2 = 0.0001; %%% (/hr) apoptosis rate IPA
p.xibar_PDGFA = 0.001; %%% PDGFA production
p.xibar_LIF = 0.002; %%% LIF production
p.Te = 0.0035; %%% tension on boundary

%%%%%%%%%%%%%%%%%%%%%%%%%%%% parameters bounds %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% note that parameters go in the following order:
%%%     kappa, mu, alpha1, alpha2, beta, gamma1, gamma2, xibar_PDGFA,
%%%     xibar_LIF, Te
bound = [1 10;
    0.01 2;
    0 1;
    0 1;
    0 1;
    0 0.0001;
    0 0.0001;
    0 0.1;
    0 0.1;
    0 0.0035];
numpar = length(bound);

N = 2;

intrange = zeros(numpar,N);
for i=1:numpar
    intrange(i,:) = linspace(bound(i,1),bound(i,2),N);
end

%%% preallocate
err = zeros(size(intrange));

for j = 1:numpar
    for i = 1:N
        paramval = p;
        if j==1
            paramval.kappa = intrange(j,i);
        elseif j==2
            paramval.mu = intrange(j,i);
        elseif j==3
            paramval.alpha1 = intrange(j,i);
        elseif j==4
            paramval.alpha2 = intrange(j,i);
        elseif j==5
            paramval.beta = intrange(j,i);
        elseif j==6
            paramval.gamma1 = intrange(j,i);
        elseif j==7
            paramval.gamma2 = intrange(j,i);
        elseif j==8
            paramval.xibar_PDGFA = intrange(j,i);
        elseif j==9
            paramval.xibar_LIF = intrange(j,i);
        elseif j==10
            paramval.Te = intrange(j,i);
        end
        
        %%% solve equation
        [t,r,c1,c2,q1,q2,mvgbdy,vel_cir,vel_rad] = eqnsolver(paramval,m,plotsonoff);

        %%% error calculation
        [err_rad,err_dens,err_tot] = errorfunction(t,r,mvgbdy,c1,c2);
        
        err(j,i) = err_tot;

        save('sensitivity_analysis.mat','p','intrange','err');
    end
end