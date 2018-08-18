%%%%%%%%%%%%% assume BC is on rmax - no-flux condition/0 Dirichlet condition

clear variables global;
clc;
% close all;

%%% time unit: hr
%%% space unit: mm

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% oxygen parameters
Pm = 1; % (mmHg)

%%% astrocyte parameters
kappa = 1;
mu = 0.000001;
alpha1 = 0.028; %%% (/hr)
alpha2 = 0.010; %%% (/hr)
beta = 1000.0525; %%% (/hr)
gamma1 = 0;%0;
gamma2 = 0;%0.5;

%%% tension parameters
rbar = 7.75*10^(-3); %%% reference radius (mm)
rproc = 15.5*10^(-3); %%% reference radius with processes included (mm)
cmin = 1/(pi*rproc^2); %%% reference cell density that includes processes (cells/mm^2)

%%% moving boundary condition parameters
Te = 0.0035; %%% tension on boundary
ce = densityatbdy(Te,kappa,cmin,rbar); % c1+c2 on boundary

Tprimeatce = Tderivative(ce,kappa,cmin,rbar); % T'(ce)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% mesh %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dr = 0.01;
dt = 0.001;

rmax = 5; %%% max radius (mm) (estimate rat retinal radius = 4.1 mm)
tmax = 50*dt; %%% max time (hr) (7 days = 168 hr)

r = 0:dr:rmax;
R = length(r);
rhalf = (r(1:R-1)+r(2:R))/2;

tol = 10^(-6); % tolerance for predictor-corrector scheme

%%%%%%%%%%%%%%%%%%%%%%%%%%% initial conditions %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% moving boundary location
s = 1;%10*dr;

if abs(s/dr - round(s/dr))>2*eps
    error('error specifying s(0): moving boundary must be located on grid node');
end

%%% time
tcurr = 0;

%%% astrocytes
c1_init = zeros(1,R);
for i=1:R
    xval = s;
    yval = 1.5*ce;
    if r(i)<=xval %fitted parabola
        c1_init(i) = (ce-yval)/s^2*r(i)^2 + yval;
    else
        c1_init(i) = 0;
    end
end
c2_init = zeros(1,R);

c1_old = c1_init;
c2_old = c2_init;

%%%%%%%%%%%%%%%%%%%%%%% initialize final variables %%%%%%%%%%%%%%%%%%%%%%%%
c1 = c1_init;
c2 = c2_init;
t = tcurr;

%%% subscript i is for space, j is for time (write in the order (j,i)) 


while tcurr < tmax
    
    %%%%%%%%%%%%%%%%%%%%%%%% solve next time step %%%%%%%%%%%%%%%%%%%%%%%%%
    PO2 = oxygen(tcurr + dt);
    
    [c1_new,c2_new] = cellpop_nomvgbdy(c1_old,c2_old,PO2,dt,r,Pm,...
        kappa,mu,alpha1,alpha2,beta,gamma1,gamma2,cmin,rbar,ce);
    
    %%%%%%%%%%%%%%%%%%%%%% reset for next time step %%%%%%%%%%%%%%%%%%%%%%%
    tcurr = tcurr + dt;
    c1_old = c1_new;
    c2_old = c2_new;
    
    %%% save variables
    c1 = [c1 ; c1_new];
    c2 = [c2 ; c2_new];
    t = [t ; tcurr];
end


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% plotting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[PO2,thickness] = oxygen(t);

T = length(t);
numcurvesplot = 10;
if T<=numcurvesplot
    ind = 1:T;
else
    ind = 1:floor(T/numcurvesplot):T;
%     ind = 1:10;
end

%%% color order
co = [0    0.4470    0.7410
    0.8500    0.3250    0.0980
    0.9290    0.6940    0.1250
    0.4940    0.1840    0.5560
    0.4660    0.6740    0.1880
    0.3010    0.7450    0.9330
    0.6350    0.0780    0.1840];

figure
subaxis(3,2,1,'MarginLeft',0.05,'MarginRight',0.01,'MarginTop',0.03,'MarginBottom',0.05)
hold on
for i=1:T
    plot(r,c1(i,:))
end
hold off
xlabel('radius r (mm)')
ylabel('c1 (APCs)')
ylims_c1 = get(gca,'YLim');

subaxis(3,2,2,'MarginLeft',0.05,'MarginRight',0.01)
hold on
for i=1:T
    plot(r,c2(i,:))
end
hold off
xlabel('radius r (mm)')
ylabel('c2 (IPAs)')
set(gca,'YLim',ylims_c1)

subaxis(3,2,3,'MarginLeft',0.05,'MarginRight',0.01)
hold on
for i=1:T
    plot(r,c1(i,:))%,'Color',co(i,:))
end
for i=1:T
    plot(r,c2(i,:),'--')%,'Color',co(i,:))
end
hold off
xlabel('radius r (mm)')
ylabel('Solid: c1 (APCs), Dashed: c2 (IPAs)')

subaxis(3,2,5,'MarginTop',0.03,'MarginBottom',0.05)
plot(t,thickness,'-o')
xlabel('t (hr)')
ylabel('total retinal thickness (mm)')

subaxis(3,2,6)
plot(thickness,PO2,'-o')
xlabel('total retinal thickness (mm)')
ylabel('PO2 (mmHg)')

set(gcf,'Units','inches','Position',[2,2,12,8],'PaperPositionMode','auto')
% legend(['t=',num2str(tplot(1))],['t=',num2str(tplot(2))],...
%     ['t=',num2str(tplot(3))],['t=',num2str(tplot(4))],...
%     ['t=',num2str(tplot(5))])
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure
% subaxis(2,2,1,'MarginLeft',0.05,'MarginRight',0.01,'MarginTop',0.03,'MarginBottom',0.05)
% plot(r,vel_cir)
% xlabel('r')
% ylabel('circumferential velocity')
% set(gca,'XLim',[0,mvgbdy(end)+5*dr])
% 
% subaxis(2,2,2,'MarginLeft',0.05,'MarginRight',0.01)
% plot(r,vel_rad)
% xlabel('r')
% ylabel('radial velocity')
% set(gca,'XLim',[0,mvgbdy(end)+5*dr])
% 
% subaxis(2,2,3,'MarginTop',0.03,'MarginBottom',0.05)
% plot(r,vel_cir+vel_rad)
% xlabel('r')
% ylabel('circumferential + radial velocity')
% set(gca,'XLim',[0,mvgbdy(end)+5*dr])
% 
% subaxis(2,2,4,'MarginTop',0.03,'MarginBottom',0.05)
% plot(r,vel_cir.*r)
% xlabel('r')
% ylabel('velocity')
% set(gca,'XLim',[0,mvgbdy(end)+5*dr])
% 
% 
% set(gcf,'Units','inches','Position',[2,2,12,8],'PaperPositionMode','auto')

