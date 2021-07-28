%   AstrocyteModel_TWS
%   Model for spread of astrocytes of surface of retina
%   Based on work of Tracy Stepien
%   TWS, June 2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Use the MatLab PDE solver PDEPE. The PDE is
%   [A] .*  D_ [u] = D_ [ A * D * Du/Dx ] +      [R]
%           Dt       Dx
%    ----     ---          -------------    -------------
%     c        u          f(x,t,u,Du/Dx)   s(x,t,u,Du/Dx)
%   where R is the reaction rate. The left bc is:
%     [0]    +   [1] .* [ A * D * Du/Dx ] = [0] 
%     ---        ---    -----------------   ---
%   p(0,t,u)     q(0,t)    f(0,t,u,Du/Dx)    0
%   and similarly at x = L

%  Field variables:
%  1 density of APCs
%  2 density of IPAs
%  3 speed of cell boundary (forced to be independent of x)
%  4 position of cell boundary (forced to be independent of x)
%  5 concentration of PDGFA
%  6 concentration of LIF

% Units are mm and hour

clear variables global;
global kTprime1 kTprime2 ce s0 ce1 cmax Lvec Pvec Pm D1 D2 rmax
global alpha10 alpha11 alpha12 alpha20 alpha21 alpha22 beta1 beta2 beta3 gamma1 gamma2
global xi1p xi2p gamma3 gamma4 LIF PDGFA nxpts ntpts tmax

% kTprime = 0.01;       % kappa / (2 mu sqrt(pi))
% ce = 1;             % cell concentration at edge
% cmax = 5;       % maximal concentration 
% s0  = 0.1;          % initial radius of cell disk
% alpha11 = 5;      % cell growth rate
% tmax = 10;

%% configuration parameters
s0 = 0.17; % initial boundary location, Chan-Ling et al. (2009) Fig 2G: E15
rmax = 5; % max radius (mm) (estimate rat retinal radius = 4.1 mm)
tmax = 7 * 24; % max time (hr) (7 days)

%% cell growth parameters
alpha10 = 0.1; % (/hr) basal proliferation rate 
alpha11 = 0.08; % 0.2 % (/hr) proliferation rate APC wrt oxygen
alpha12 = 0.1; % 0.5 % (/hr) proliferation rate APC wrt PDGFA
alpha20 = 0.0; % (/hr) basal proliferation rate 
alpha21 = 0.0; % 0.2 % (/hr) proliferation rate IPA wrt oxygen
alpha22 = 0.0; % 0.5 % (/hr) proliferation rate IPA wrt PDGFA
beta1 = 0.08; % (/hr) basal differentiation rate 
beta2 = 0.03; % (/hr) differentiation rate wrt oxygen
beta3 = 0.02; % (/hr) differentiation rate wrt LIF
gamma1 = 0.0; % (/hr) apoptosis rate APC
gamma2 = 0.0; % (/hr) apoptosis rate IPA

%% tension parameters 
rbar = 7.75*10^(-3); % reference radius (mm)
rproc = 15.5*10^(-3); % reference radius with processes included (mm)
cmin = 1/(pi*rproc^2); % reference cell density with processes (cells/mm^2)
cmax = 1/(pi*rbar^2); % reference cell density - cell body (cells/mm^2)
kappa = 1; % tension function scaling (mN/mm^2)
mu1 = 0.2; % adhesion constant (mm/hr/(mN/mm^2))
mu2 = 0.5;
% fun = @(x) kappa*x^2/(cmin^2 + x^2)*(1/sqrt(pi*x) - rbar) - Te ;
ce = 1000; % /mm^2, estimated from Figure 2B of Fruttiger 2007
% tension on boundary (mN/mm)
Te = kappa * (1/sqrt(pi*ce) - rbar) ; % simpler form, June 2021
kTprime1 = kappa / (2 * mu1 * sqrt(pi));
kTprime2 = kappa / (2 * mu2 * sqrt(pi));
% parameter in initial condition, chosen to match subsequent solution
ce1 = 0.5 * alpha11 * s0^2 * sqrt(ce) * (1 - ce/cmax) / kTprime1;

%% Oxygen parameters
P0 = 60; % mmHg
Pm = 10.5; % mmHg
Dalpha = 4.73*10^(-10); % cm^3 O2/cm/s/mmHg
Dalpha = Dalpha * (60*60*0.1); % cm^3 O2/mm/hr/mmHg
M0 = 1.8; % cm^3 O2/100g/min
M0 = M0 * (60/100*0.1^3); % cm^3 O2/mm^3/hr

%% Growth factor production and diffusion
lambda = 1.6; % tortuosity of medium (dimensionless)
phi = 0.2; % porosity/volume fraction in extracellular space
% diffusion of PDGFA in water at 37C (cm^2/s), converted to (mm^2/hr)
Dwater_PDGFA = 1.32*10^(-6) *(60^2*10^2);
D1 = Dwater_PDGFA / lambda^2; % effective diffusivity of PDGFA
% diffusion of LIF in water at 37C (cm^2/s), converted to (mm^2/hr)
Dwater_LIF = 1.33*10^(-6) *(60^2*10^2);
D2 = Dwater_LIF / lambda^2; % effective diffusivity of LIF

% degradation rates
quasilength = 0.2;
gamma3 = 4.6; %D1/quasilength^2;
gamma4 = 4.7; %D2/quasilength^2;

% production rates
xi1p = 4.6; %gamma3;
xi2p = 4.7; %gamma4;

%% set up and solve
nxpts = 26;
ntpts = 29;
x = linspace(0,1,nxpts);    % 0 <= x <= 1 by definition
t = linspace(0,tmax,ntpts);

% Set up vectors for interpolation to obtain PO2 
[Lvec, Pvec] = oxygen_setup(M0, Dalpha, Pm, P0);

% Calculate growth factors for interpolation using fixed domain
r = linspace(0,rmax,26);
sol = pdepe(1,@GF_PDE,@GF_IC,@GF_BC,r,t);
PDGFA = sol(:,:,1);
LIF = sol(:,:,2);

% Calculate cell densities using time-dependent domain stretching
sol = pdepe(1,@AstroPDE,@AstroIC,@AstroBC,x,t);

%% plot results

figure(1);  % check boundary condition at outer edge satisfied
plot(t,sol(:,nxpts,1)+sol(:,nxpts,2))
figure(2);
surf(x,t,sol(:,:,1));
figure(3);
surf(x,t,sol(:,:,2));
figure(4);
surf(x,t,sol(:,:,3));
figure(5);
surf(x,t,sol(:,:,4));

% color order
co = [0    0.4470    0.7410
    0.8500    0.3250    0.0980
    0.9290    0.6940    0.1250
    0.4940    0.1840    0.5560
    0.4660    0.6740    0.1880
    0.3010    0.7450    0.9330
    0.6350    0.0780    0.1840
    0         0.3725    0];

fslabel = 12;
numcurvesplot = floor(ntpts/4)+1;

figure(6);
tiledlayout(3,3,'TileSpacing','Compact')
% extra points to make plots go to axis at outer edge
r = zeros(numcurvesplot,nxpts+2);
y1 = zeros(numcurvesplot,nxpts+2);
y2 = zeros(numcurvesplot,nxpts+2);

% set up retinal thickness and choroid oxygen
for i=1:numcurvesplot
    [~,~,~,radius_ret] = thick_rad(t(i*4-3),0);
    for j = 1:nxpts
        r(i,j) = radius_ret * x(j);
        [thickness_ret,~,~,~] = thick_rad(t(i*4-3),r(i,j));
        y1(i,j) = thickness_ret;
        y2(i,j) = interp1(Lvec,Pvec,thickness_ret);        
    end
    r(i,nxpts+1) = r(i,nxpts);   %to make plots go to axis at outer edge 
    r(i,nxpts+2) = rmax;  
end

% retinal thickness
ax1 = nexttile;
hold on
for i=1:numcurvesplot
    plot(r(i,:),y1(i,:),'LineWidth',1.5,'Color',co(i,:)) 
end
xlabel('Radius (mm)','FontSize',fslabel)
ylabel('Retinal thickness (mm)','FontSize',fslabel)
hold off

% choroid oxygen
ax2 = nexttile;
hold on
for i=1:numcurvesplot
    plot(r(i,:),y2(i,:),'LineWidth',1.5,'Color',co(i,:)) 
end
xlabel('Radius (mm)','FontSize',fslabel)
ylabel('Choroid P_{O2} (mmHg)','FontSize',fslabel)
hold off

% set up PDGFA and LIF
for i=1:numcurvesplot
    for j = 1:nxpts
        r(i,j) = rmax * x(j);
        y1(i,j) = PDGFA(i*4-3,j);
        y2(i,j) = LIF(i*4-3,j);        
    end
    r(i,nxpts+1) = r(i,nxpts);
end

% PDGFA
ax3 = nexttile;
hold on
for i=1:numcurvesplot
   plot(r(i,:),y1(i,:),'LineWidth',1.5,'Color',co(i,:)) 
end
xlabel('Radius (mm)','FontSize',fslabel)
ylabel('PDGFA (ng/mL)','FontSize',fslabel)
hold off

% LIF
ax4 = nexttile;
hold on
for i=1:numcurvesplot
   plot(r(i,:),y2(i,:),'LineWidth',1.5,'Color',co(i,:)) 
end
xlabel('Radius (mm)','FontSize',fslabel)
ylabel('LIF (ng/mL)','FontSize',fslabel)
hold off

% set up cell densities
rm = sol(:,nxpts,4);
for i=1:numcurvesplot
    for j = 1:nxpts
        r(i,j) = rm(i*4-3) * x(j);
        y1(i,j) = sol(i*4-3,j,1);
        y2(i,j) = sol(i*4-3,j,2);        
    end
    r(i,nxpts+1) = r(i,nxpts);
end

% c1 - APCs
ax5 = nexttile;
hold on
for i=1:numcurvesplot
    plot(r(i,:),y1(i,:),'LineWidth',1.5,'Color',co(i,:))
end
xlabel('Radius (mm)','FontSize',fslabel)
ylabel('APCs (cells/mm^2)','FontSize',fslabel)
hold off

% c2 - IPAs
ax6 = nexttile;
hold on
for i=1:numcurvesplot
    plot(r(i,:),y2(i,:),'LineWidth',1.5,'Color',co(i,:))
end
xlabel('Radius (mm)','FontSize',fslabel)
ylabel('IPAs (cells/mm^2)','FontSize',fslabel)
hold off

% c1 + c2 - APCs + IPAs
ax7 = nexttile;
hold on
for i=1:numcurvesplot
    plot(r(i,:),y1(i,:)+y2(i,:),'LineWidth',1.5,'Color',co(i,:))
end
xlabel('Radius (mm)','FontSize',fslabel)
ylabel('APCs + IPAs (cells/mm^2)','FontSize',fslabel)
hold off

%set up radial velocity in interior of domain
for i=1:numcurvesplot
    usum = sol(i*4-3,:,1) + sol(i*4-3,:,2);
    for j = 1:nxpts
    if j == 1 
        DuDx = (nxpts - 1) * (usum(j+1) - usum(j));
    elseif  j == nxpts
        DuDx = (nxpts - 1) * (usum(j) - usum(j-1));
    else
        DuDx = (nxpts - 1)/2 * (usum(j+1) - usum(j-1));
    end
    y1(i,j) = -kTprime1 * usum(j)^(-3/2) / sol(i*4-3,j,4) / mu1 * DuDx;
    end
end

% radial velocity
ax8 = nexttile;
hold on
for i=1:numcurvesplot
    plot(r(i,:),y1(i,:),'LineWidth',1.5,'Color',co(i,:))
end
xlabel('Radius (mm)','FontSize',fslabel)
ylabel('APC velocity (mm/hr)','FontSize',fslabel)
hold off

axis([ax1,ax2,ax3,ax4,ax5,ax6,ax7,ax8],[0,5,-inf,inf]);

% Radius vs. time
% distance astrocytes spread
rad_APC = [0.17; 0.33; 0.5; 0.67; 1.67; 2.17; 2.67]; %not including E17
rad_days = [0; 1; 3; 4; 5; 6; 7];
nexttile
hold on
plot(t/24,rm,'k','LineWidth',1.5)
scatter(rad_days,rad_APC,150,[0.5 0.5 0.5],'x','LineWidth',1.5)
hold off
xlabel('t (days)','FontSize',fslabel)
ylabel('Cell boundary (mm)','FontSize',fslabel)

set(gcf,'Units','inches','Position',[1,1,10,7],'PaperPositionMode','auto');
