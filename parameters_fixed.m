%%% this file contains all of the fixed parameters that are known/derived
%%% from literature

%%%%%%%%%%%%%%%%%%%%%%%%%%%% oxygen parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%
Pm = 10.5; % (mmHg)


%%%%%%%%%%%%%%%%%%%%%%%% growth factor parameters %%%%%%%%%%%%%%%%%%%%%%%%%
lambda = 1.6; %%% tortuosity of medium (dimensionless)
phi = 0.2; %%% porosity/volume fraction in extracellular space (%)

%%% diffusion of PDGFA in water at 37C (cm^2/s), converted to (mm^2/hr)
Dwater_PDGFA = 1.32*10^(-6) *(60^2*10^2);
%%% diffusion of LIF in water at 37C (cm^2/s), converted to (mm^2/hr)
Dwater_LIF = 1.33*10^(-6) *(60^2*10^2);

D1 = Dwater_PDGFA / lambda^2; %%% effective diffusivity of PDGFA
D2 = Dwater_LIF / lambda^2; %%% effective diffusivity of LIF

%%% degradation rates
quasilength = 0.2;
gamma3 = D1/quasilength^2;
gamma4 = D2/quasilength^2;


%%%%%%%%%%%%%%%%%%%%%%%%%%% tension parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%
rbar = 7.75*10^(-3); %%% reference radius (mm)
rproc = 15.5*10^(-3); %%% reference radius with processes included (mm)
cmin = 1/(pi*rproc^2); %%% reference cell density that includes processes
                       %%% (cells/mm^2)


%%%%%%%%%%%%%%%%%%%% moving boundary initial location %%%%%%%%%%%%%%%%%%%%%
s0 = 0.2; %%% based on Chan-Ling et al. (2009) Fig 2G: E15