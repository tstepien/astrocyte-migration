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
gamma3 = 4.6; %D1/quasilength^2;
gamma4 = 4.7; %D2/quasilength^2;

%%% production rates
xibar_PDGFA = 0.92; %phi*gamma3;
xibar_LIF = 0.94; %phi*gamma4;


%%%%%%%%%%%%%%%%%%%%%%%%%%% tension parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%
rbar = 7.75*10^(-3); %%% reference radius (mm)
rproc = 15.5*10^(-3); %%% reference radius with processes included (mm)
cmin = 1/(pi*rproc^2); %%% reference cell density that includes processes
                       %%% (cells/mm^2)
cmax = 1/(pi*rbar^2); %%% reference cell density that includes only the 
                      %%% cell body (cells/mm^2)
kappa = 1; %%% tension function scaling (mN/mm^2)


%%%%%%%%%%%%%%%%%%%% moving boundary initial location %%%%%%%%%%%%%%%%%%%%%
s0 = 0.17; %%% based on Chan-Ling et al. (2009) Fig 2G: E15


%%%%%%%%%%%%%%%%%%%%%%%%% retinal ganglion cells %%%%%%%%%%%%%%%%%%%%%%%%%%
maxRGCthick = 46*0.001; %%% maximum thickness of RGC layer: 46 micon converted to mm


%%%%%%%%%%%%%%%%%%%%%%% partial pressure of oxygen %%%%%%%%%%%%%%%%%%%%%%%%
P0 = 60;

Dalpha = 4.73*10^(-10); % cm^3 O2/cm/s/mmHg
M0 = 1.8; % cm^3 O2/100g/min

%%% convert parameters to time units of hours and space units of mm
Dalpha = Dalpha * (60*60*0.1); % cm^3 O2/mm/hr/mmHg
M0 = M0 * (60/100*0.1^3); % cm^3 O2/hr

%%% check calculations
%     M0to_s = 1.8/(60*100);
%     Dalpha/M0to_s
%     check1 = 2*P0*Dalpha/M0to_s
%     sq_check1 = sqrt(check1)

%%% check calculations
%     Dalpha/M0
%     check2 = 2*P0*Dalpha/M0
%     sq_check2 = sqrt(check2)