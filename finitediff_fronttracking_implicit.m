clear variables global;
clc;
close all;

global whatstep tcurr;

%%% time unit: hr
%%% space unit: mm

p1=[]; %set this empty so that plots will work

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% oxygen parameters
Pm = 10; % (mmHg)

%%% astrocyte parameters
kappa = 1;
mu = 0.1;
alpha1 = 0.3; %%% (/hr)
alpha2 = alpha1/2; %%% (/hr)
beta = 0.01; %%% (/hr)
gamma1 = 0;%0.0001;%0;
gamma2 = 0;%0.0001;%0.5;

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

rmax = 5; %%% max radius (mm) (estimate rat retinal radius = 4.1 mm)
tmax = 7*24; %%% max time (hr) (7 days = 168 hr)

r = 0:dr:rmax;
R = length(r);
rhalf = (r(1:R-1)+r(2:R))/2;

tol = 10^(-6); % tolerance for predictor-corrector scheme

%%%%%%%%%%%%%%%%%%%%%%%%%%% initial conditions %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% moving boundary location
s = 0.2; %based on Chan-Ling et al. (2009) Fig 2G: E15

if abs(s/dr - round(s/dr))>2*eps
    error('error specifying s(0): moving boundary must be located on grid node');
end

%%% time
tcurr = 0;

%%% astrocytes
c1_init = zeros(1,R);
for i=1:R
    xval = s;
    yval = 1.05*ce;
    if r(i)<=xval %fitted parabola
        c1_init(i) = (ce-yval)/s^2*r(i)^2 + yval;
    else
        c1_init(i) = 0;
    end
end
c2_init = zeros(1,R);

c1_old = c1_init;
c2_old = c2_init;
k_old = c1_old + c2_old;

%%%%%%%%%%%%%%%%%%%%%%% initialize final variables %%%%%%%%%%%%%%%%%%%%%%%%
mvgbdy = s;
c1 = c1_init;
c2 = c2_init;
k = k_old;
t = tcurr;

%%% subscript i is for space, j is for time (write in the order (j,i)) 

j_init = s/dr+1;
j = j_init;

%%% velocity
[vel_cir,vel_rad] = velocity(j,c1,c2,r,kappa,cmin,rbar,mu);

mvgbdy_vel = [];

%%% concentration on the moving boundary (mb)
c1mb = c1(j);
c2mb = c2(j);

while tcurr < tmax && j<R-1
    %%%%%%%%%%%%%%%%%%%%%%%%%%% predictor step %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    whatstep = 'predictor';
    
    aa = 1;
    if s==0
        dt_p = dr;
    elseif s==dr
        dt_p = 1/aa * mu/Tprimeatce * dr^2 / ( k_old(j) - k_old(j-1) );
    else
        dt_p = 1/aa * mu/Tprimeatce * 2*dr^2 / ...
            ( 3*k_old(j) - 4*k_old(j-1) + k_old(j-2) );
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%% solve eqn's with dt_p %%%%%%%%%%%%%%%%%%%%%%%%
    dt_c = 0;
    while abs(dt_p-dt_c)>=tol
        PO2 = oxygen(tcurr + dt_p,r);
        
        ve_old = ve_calc(j,tcurr,r,c1_old,c2_old,Pm,alpha1,alpha2,gamma1,gamma2,ce);
        
        k_hat = cellpops_sum(j,c1_old,c2_old,PO2,dt_p,r,Pm,kappa,mu,...
            alpha1,alpha2,gamma1,gamma2,cmin,rbar,ce);
    
        [c1_hat,c2_hat] = cellpops_separate(j,c1_old,c2_old,k_hat,PO2,dt_p,...
            r,Pm,kappa,mu,alpha1,alpha2,beta,gamma1,gamma2,cmin,rbar,ce);
        
        %%%%%%%%%%%%%%%%%%%%%%%%% corrector step %%%%%%%%%%%%%%%%%%%%%%%%%%
        bb = 1/2;
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
        
%         [3*k_old(j) - 4*k_old(j-1) + k_old(j-2) , 3*k_hat(j+1) - 4*k_hat(j) + k_hat(j-1)]
        
%         [dt_p, dt_c]
%         keyboard
        
        if abs(dt_p-dt_c)<tol
            break;
        else
            dt_p = dt_c;
            dt_c = 0;
%             dt_hold = dt_c;
%             dt_c = 0;
%             dt_p = dt_hold;
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%% solve next time step %%%%%%%%%%%%%%%%%%%%%%%%%
    whatstep = 'corrector';
    
    PO2 = oxygen(tcurr + dt_c,r);
    
    k_new = cellpops_sum(j,c1_old,c2_old,PO2,dt_c,r,Pm,kappa,mu,...
        alpha1,alpha2,gamma1,gamma2,cmin,rbar,ce);
    
    [c1_new,c2_new] = cellpops_separate(j,c1_old,c2_old,k_new,PO2,dt_c,...
        r,Pm,kappa,mu,alpha1,alpha2,beta,gamma1,gamma2,cmin,rbar,ce);
    
    %%%%%%%%%%%%%%%%%%%%%% reset for next time step %%%%%%%%%%%%%%%%%%%%%%%
    j = j+1;
    s = s + dr;
    tcurr = tcurr + dt_c;
    c1_old = c1_new;
    c2_old = c2_new;
    k_old = k_new;
    
    %%% save variables
    mvgbdy = [mvgbdy ; s];
    c1 = [c1 ; c1_new];
    c2 = [c2 ; c2_new];
    k = [k ; k_new];
    t = [t ; tcurr];
    c1mb = [c1mb ; c1_new(j)];
    c2mb = [c2mb ; c2_new(j)];
    
    %%% velocity calculation
    [vel_cir_new,vel_rad_new] = velocity(j,c1_new,c2_new,r,kappa,cmin,rbar,mu);
    vel_cir = [vel_cir ; vel_cir_new];
    vel_rad = [vel_rad ; vel_rad_new];
    
    mvgbdy_vel = [mvgbdy_vel ; (mvgbdy(end)-mvgbdy(end-1))/(t(end)-t(end-1))];
    
    if mvgbdy_vel(end)<10^(-4) || isreal(mvgbdy_vel(end))==0
        disp('***stopping simulation since moving boundary velocity is <10^(-4)***')
        break;
    elseif mvgbdy_vel(end)<10^(-3)
        disp('**fyi moving boundary is moving slow**')
    end
end

disp(['location of moving boundary at last time step: ',num2str(mvgbdy(end))])
% disp('PO2 at last time step')
% PO2(end)
disp(['ending time in hours: ',num2str(t(end)/24)])
disp(['velocity of moving boundary at last time step: ',num2str(mvgbdy_vel(end))])

%% plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% plotting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot_the_plots


%% area under curve - trapezoidal rule
%%%%%%%%%%%%%%%%%%%%%% testing: conservation of mass %%%%%%%%%%%%%%%%%%%%%%
areaundercurve = zeros(length(t),1);
for i = 1:length(t)
    areaundercurve(i) = dr*sum( k(1:mvgbdy(1)/dr+(i-1)) ...
        + k(2:mvgbdy(1)/dr+1+(i-1)) ) / 2;
end