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

%%% set all IPA stuff = 0
p.alpha21 = 0; %%% (/hr) proliferation rate IPA wrt oxygen
p.alpha22 = 0; %%% (/hr) proliferation rate IPA wrt PDGFA
p.beta = 0; %%% (/hr) differentiation rate
p.beta_hat = 0; %%% (/hr) mass action rate
p.gamma2 = 0; %%% (/hr) apoptosis rate IPA

%%%%%%%%%%%%%%%%%%%%%%%%%%% baseline parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%
p.mu = 0.1; %%% adhesion constant
p.alpha11 = 0.03; %%% (/hr) proliferation rate APC wrt oxygen
p.alpha12 = 0.22; %%% (/hr) proliferation rate APC wrt PDGFA
p.gamma1 = 0.01; %%% (/hr) apoptosis rate APC
p.Te = 0.0035; %%% tension on boundary
p.P_hy = 10; %%% partial pressure of oxygen due to hyaloid artery
p.r_hy = 0.1; %%% radius at half-maximum of Hill function for hyaloid

%%%%%%%%%%%%%%%%%%%%%%%%%%%% parameters bounds %%%%%%%%%%%%%%%%%%%%%%%%%%%%
bound = [0.01 5; %mu
    0 1; %alpha11
    0 1; %alpha12
    0 0.1; %gamma1
    0 0.0038; %Te
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
        paramval = p;
        if j==1
            paramval.mu = intrange(j,i);
        elseif j==2
            paramval.alpha11 = intrange(j,i);
        elseif j==3
            paramval.alpha12 = intrange(j,i);
        elseif j==4
            paramval.gamma1 = intrange(j,i);
        elseif j==5
            paramval.Te = intrange(j,i);
        elseif j==6
            paramval.P_hy = intrange(j,i);
        elseif j==7
            paramval.r_hy = intrange(j,i);
        end
        
        %%% solve equation
        tic
        [t,~,~,~,~,~,mvgbdy,~,~] = eqnsolver(paramval,m,plotsonoff);
        toc

        %%% error calculation
        err_rad = errorfunction_1pop(t,mvgbdy);
        
        err(j,i) = err_rad;

        save('sensitivity_analysis.mat','p','intrange','err','bound');
    end
end

save('sensitivity_analysis.mat','p','intrange','err','bound')

diary off