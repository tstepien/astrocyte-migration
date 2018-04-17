clear variables global;
clc;

%%% time unit: hr
%%% space unit: mm

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% astrocyte parameters
kappa = 1;
mu = 0.01;
alpha1 = 0.028; %%% (/hr)
alpha2 = 0.010; %%% (/hr)
beta = 0.0525; %%% (/hr)
gamma1 = 0;%0;
gamma2 = 0;%0.5;

%%% growth factor parameters
lambda = 1.6; %%% tortuosity of medium (unitless)
phi = 0.2; %%% porosity/volume fraction in extracellular space (%)
Dwater_PDGFA = 1.2*10^(-6) *(60*60*10^2); %%% diffusion of PDGFA in water at 37C (cm^2/s)
                                     %%% but then converted to (mm^2/hr)
Dwater_LIF = 1.38*10^(-6) *(60*60*10^2); %%% diffusion of LIF in water at 37C (cm^2/s)
                                    %%% but then converted to (mm^2/hr)
eta_PDGFA = 0;
eta_LIF = 0;

D1 = Dwater_PDGFA / lambda^2; %%% effective diffusivity of PDGFA
D2 = Dwater_LIF / lambda^2; %%% effective diffusivity of LIF
eta1 = eta_PDGFA / phi; %%% production/release rate of PDGFA
eta2 = eta_LIF / phi; %%% production/release rate of LIF

p1max = 1;
p1beta = 0.45;
p1epsilon = 0.45;

%%% tension parameters
rbar = 7.75*10^(-3); %%% reference radius (mm)
rproc = 15.5*10^(-3); %%% reference radius with processes included (mm)
cmin = 1/(pi*rproc^2); %%% reference cell density that includes processes (cells/mm^2)

%%% moving boundary condition parameters
Te = 0.0035; %%% tension on boundary
ce = densityatbdy(Te,kappa,cmin,rbar); % c1+c2 on boundary

Tprimeatce = Tderivative(ce,kappa,cmin,rbar); % T'(ce)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% mesh %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dr = 0.001;

rmax = 5; %%% max radius (mm) (estimate rat retinal radius = 4.1 mm)
tmax = (1/24)*24; %%% max time (hr) (7 days = 168 hr)

r = 0:dr:rmax;
R = length(r);
rhalf = (r(1:R-1)+r(2:R))/2;

tol = 10^(-6); % tolerance for predictor-corrector scheme

%%%%%%%%%%%%%%%%%%%%%%%%%%% initial conditions %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% moving boundary location
s = 10*dr;

if abs(s/dr - round(s/dr))>2*eps
    error('error specifying s(0): moving boundary must be located on grid node');
end

%%% time
tcurr = 0;

%%% astrocytes
c1_init = zeros(1,R);
for i=1:R
%     xval1 = 5*dr;
%     xval2 = s;
%     yval1 = 1.5*ce;
%     yval2 = ce;
%     if r(i)<=xval1
% %         c1_init(i) = ce;
%         c1_init(i) = 1.5*ce;
%     elseif r(i)<=xval2 && r(i)>xval1
% %         c1_init(i) = ce;
%         sl = (yval2-yval1)/(xval2-xval1);
%         c1_init(i) = sl*(r(i)-xval1) + yval1;
    xval2 = s;
    yval = 1.5*ce;
    if r(i)<=xval2 %fitted parabola
        c1_init(i) = (ce-yval)/s^2*r(i)^2 + yval;
    else
        c1_init(i) = 0;
    end
end
c2_init = zeros(1,R);

c1_old = c1_init;
c2_old = c2_init;

%%% growth factors
p1_old = p1max * (1 - p1beta*exp(-r.^2/p1epsilon));%100*ones(1,R);
p1_old(1) = 0;
p2_old = 10*ones(1,R);%smooth(100*heaviside(1-r))';

%%%%%%%%%%%%%%%%%%%%%%% initialize final variables %%%%%%%%%%%%%%%%%%%%%%%%
mvgbdy = s;
c1 = c1_init;
c2 = c2_init;
p1 = p1_old;
p2 = p2_old;
t = tcurr;
% Tension = [];
% vhalf = [];

%%% subscript i is for space, j is for time (write in the order (j,i)) 

j_init = s/dr+1;
j = j_init;

%velocity
[vel_cir,vel_rad] = velocity(j,c1,c2,r,kappa,cmin,rbar,mu);

etaf = [];
c1mb = c1(j);
c2mb = c2(j);

while tcurr < tmax && j<R-1
    %%%%%%%%%%%%%%%%%%%%%%%%%%% predictor step %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    aa = 1;
    k_old = c1_old + c2_old;
    if s==0
        dt_p = dr;
    elseif s==dr
        dt_p = 1/aa * mu/Tprimeatce * dr^2 / ( k_old(j) - k_old(j-1) );
    else
%         disp('new dt_p')
        dt_p = 1/aa * mu/Tprimeatce * 2*dr^2 / ...
            ( 3*k_old(j) - 4*k_old(j-1) + k_old(j-2) );
    end
%     keyboard
    
    dt_c = 0;
    while abs(dt_p-dt_c)>=tol        
        [p1_hat,p2_hat] = growthfactors_implicit(p1_old,p2_old,dt_p,r,D1,...
            D2,eta1,eta2);
    
        [c1_hat,c2_hat] = cellpops_implicit_splitbc(j,c1_old,c2_old,p1_hat,...
            p2_hat,dt_p,r,kappa,mu,alpha1,alpha2,beta,gamma1,gamma2,...
            cmin,rbar,ce);%,tcurr,a,n);
        
%         if c1_hat==c1_old & c2_hat==c2_old
%             continue;
%         end
        
%         if abs(vhalf_hat(j-1)-vhalf_hat(j))<eps
%             break;
%         end
        
        %%%%%%%%%%%%%%%%%%%%%%%%% corrector step %%%%%%%%%%%%%%%%%%%%%%%%%%
        bb = 1/2;
        k_hat = c1_hat + c2_hat;
        if s==0
            dt_c = mu/Tprimeatce * dr / ( ...
                bb*( k_hat(j+1) - k_hat(j) )/dr ...
                + (1-bb)*( 1/aa*mu/Tprimeatce ) );
        elseif s==dr
            dt_c = mu/Tprimeatce * dr^2 / ( ...
                bb*( 3*k_hat(j+1) - 4*k_hat(j) + k_hat(j-1) )/2 ...
                + (1-bb)*( k_old(j) - k_old(j-1) ) );
        else
            dt_c = mu/Tprimeatce * 2*dr^2 / ( ...
                bb*( 3*k_hat(j+1) - 4*k_hat(j) + k_hat(j-1) ) ...
                + (1-bb)*( 3*k_old(j) - 4*k_old(j-1) + k_old(j-2) ) );
        end
        
%         [dt_p, dt_c]
%         keyboard
        
        if abs(dt_p-dt_c)<tol
            break;
        else
            dt_hold = dt_c;
            dt_c = 0;
            dt_p = dt_hold;
        end
    end
    
%     if abs(vhalf_hat(j-1)-vhalf_hat(j))<eps
%         disp('moving boundary has speed=0');
%         break;
%     end
    
    %%%%%%%%%%%%%%%%%%%%%%%% solve next time step %%%%%%%%%%%%%%%%%%%%%%%%%    
    [p1_new,p2_new] = growthfactors_implicit(p1_old,p2_old,dt_c,r,D1,...
        D2,eta1,eta2);
    
    [c1_new,c2_new] = cellpops_implicit_splitbc(j,c1_old,c2_old,p1_new,...
        p2_new,dt_c,r,kappa,mu,alpha1,alpha2,beta,gamma1,gamma2,...
        cmin,rbar,ce);%,tcurr,a,n);
%     etafrac=0;
        
%     if c1_new==c1_old & c2_new==c2_old
%         continue;
%     end

    
    %%% reset for next time step
    j = j+1;
    s = s + dr;
    tcurr = tcurr + dt_c;
    c1_old = c1_new;
    c2_old = c2_new;
    p1_old = p1_new;
    p2_old = p2_new;
    
    %%% save variables
    mvgbdy = [mvgbdy ; s];
    c1 = [c1 ; c1_new];
    c2 = [c2 ; c2_new];
    p1 = [p1 ; p1_new];
    p2 = [p2 ; p2_new];
    t = [t ; tcurr];
%     Tension = [Tension ; Tension_old];
%     vhalf = [vhalf ; vhalf_old];
%     etaf = [etaf ; etafrac];
    c1mb = [c1mb ; c1_new(j)];
    c2mb = [c2mb ; c2_new(j)];
    
    %%% velocity calculation
    [vel_cir_new,vel_rad_new] = velocity(j,c1_new,c2_new,r,kappa,cmin,rbar,mu);
    vel_cir = [vel_cir ; vel_cir_new];
    vel_rad = [vel_rad ; vel_rad_new];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% plotting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plots
T = length(t);		
numcurvesplot = 10;		
if T<=numcurvesplot		
    ind = 1:T;		
else		
    ind = 1:floor(T/numcurvesplot):T;		
%     ind = 1:10;		
end

figure
subaxis(3,2,1,'MarginLeft',0.05,'MarginRight',0.01,'MarginTop',0.03,'MarginBottom',0.05)
plot(r,c1)
xlabel('r')
ylabel('c1 (APCs)')
set(gca,'XLim',[0,mvgbdy(end)+5*dr])

subaxis(3,2,2,'MarginLeft',0.05,'MarginRight',0.01)
plot(r,c2)
xlabel('r')
ylabel('c2 (IPAs)')
set(gca,'XLim',[0,mvgbdy(end)+5*dr])

subaxis(3,2,3,'MarginTop',0.03,'MarginBottom',0.05)
plot(r,p1)
xlabel('r')
ylabel('p1 (PDGFA)')

subaxis(3,2,4)
plot(r,p2)
xlabel('r')
ylabel('p2 (LIF)')

subaxis(3,2,5,'MarginTop',0.03,'MarginBottom',0.05)
plot(t,mvgbdy)
xlabel('t')
ylabel('moving boundary s(t)')

set(gcf,'Units','inches','Position',[2,2,12,8],'PaperPositionMode','auto')
% legend(['t=',num2str(tplot(1))],['t=',num2str(tplot(2))],...
%     ['t=',num2str(tplot(3))],['t=',num2str(tplot(4))],...
%     ['t=',num2str(tplot(5))])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
subaxis(2,2,1,'MarginLeft',0.05,'MarginRight',0.01,'MarginTop',0.03,'MarginBottom',0.05)
plot(r,vel_cir)
xlabel('r')
ylabel('circumferential velocity')
set(gca,'XLim',[0,mvgbdy(end)+5*dr])

subaxis(2,2,2,'MarginLeft',0.05,'MarginRight',0.01)
plot(r,vel_rad)
xlabel('r')
ylabel('radial velocity')
set(gca,'XLim',[0,mvgbdy(end)+5*dr])

subaxis(2,2,3,'MarginTop',0.03,'MarginBottom',0.05)
plot(r,vel_cir+vel_rad)
xlabel('r')
ylabel('circumferential + radial velocity')
set(gca,'XLim',[0,mvgbdy(end)+5*dr])

subaxis(2,2,4,'MarginTop',0.03,'MarginBottom',0.05)
plot(r,vel_cir.*r)
xlabel('r')
ylabel('velocity')
set(gca,'XLim',[0,mvgbdy(end)+5*dr])


set(gcf,'Units','inches','Position',[2,2,12,8],'PaperPositionMode','auto')


%%%%%%%%%%%%%%%%%%%%%% testing: conservation of mass %%%%%%%%%%%%%%%%%%%%%%
% area under curve - trapezoidal rule
areaundercurve = zeros(length(t),1);
for i = 1:length(t)
    areaundercurve(i) = dr*sum( c1(1:mvgbdy(1)/dr+(i-1)) ...
        + c1(2:mvgbdy(1)/dr+1+(i-1)) ) / 2;
end