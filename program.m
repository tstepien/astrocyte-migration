clear variables global;
clc;
% close all;

global whatstep tcurr;

%%% time unit: hr
%%% space unit: mm

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% VARIABLE CODE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%% astrocyte parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%
kappa = 1; %%% tension function scaling
mu = 0.1; %%% adhesion constant
alpha1 = 0.1975; %%% (/hr) proliferation rate APC
alpha2 = alpha1 *(7/12); %%% (/hr) proliferation rate IPA
beta = 0.03; %%% (/hr) differentiation rate
gamma1 = 0;%0.0001;%0; apoptosis rate APC
gamma2 = 0;%0.0001;%0.5; apoptosis rate IPA


%%%%%%%%%%%%%%%%%%%%%%%% growth factor parameters %%%%%%%%%%%%%%%%%%%%%%%%%
xibar_PDGFA = 0.015/15; %%% PDGFA production
xibar_LIF = 0.015/7; %%% LIF production


%%%%%%%%%%%%%%%%%% moving boundary condition parameters %%%%%%%%%%%%%%%%%%%
Te = 0.0035; %%% tension on boundary


%%%%%%%%%%%%%%%%%%%%%%%%%%%%% mesh parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dr = 0.01;

rmax = 5; %%% max radius (mm) (estimate rat retinal radius = 4.1 mm)
tmax = 7*24; %%% max time (hr) (7 days = 168 hr)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% FIXED CODE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%% input all fixed parameters that are %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% known/derived from literature %%%%%%%%%%%%%%%%%%%%%%
parameters_fixed


%%%%%%%%%%%%%%%%%%%%% parameter scalings/calculations %%%%%%%%%%%%%%%%%%%%%
xi1 = xibar_PDGFA / phi;
xi2 = xibar_LIF / phi;

ce = densityatbdy(Te,kappa,cmin,rbar); % c1+c2 on boundary
Tprimeatce = Tderivative(ce,kappa,cmin,rbar); % T'(ce)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% mesh set up %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
r = 0:dr:rmax;
R = length(r);
rhalf = (r(1:R-1)+r(2:R))/2;

tol = 10^(-6); % tolerance for predictor-corrector scheme

if abs(s0/dr - round(s0/dr))>2*eps
    error('error specifying s(0): moving boundary must be located on grid node');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%% initial conditions %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% time
tcurr = 0;

%%% astrocytes
c1_init = zeros(1,R);
for i=1:R
    xval = s0;
    yval = 1.05*ce;
    if r(i)<=xval %fitted parabola
        c1_init(i) = (ce-yval)/s0^2*r(i)^2 + yval;
    else
        c1_init(i) = 0;
    end
end
c2_init = zeros(1,R);

c1_old = c1_init;
c2_old = c2_init;
k_old = c1_old + c2_old;

%%% growth factors
q1_old = zeros(1,R);
q2_old = zeros(1,R);


%%%%%%%%%%%%%%%%%%%%%%% initialize final variables %%%%%%%%%%%%%%%%%%%%%%%%
mvgbdy = s0;
s = s0;
c1 = c1_init;
c2 = c2_init;
k = k_old;
q1 = q1_old;
q2 = q2_old;
t = tcurr;

%%% subscript i is for space, j is for time (write in the order (j,i))
j_init = s0/dr+1;
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
    %%% should have s0 start further than 2 nodes in
    %%% lines commented so code runs faster, but left in case they're
    %%% desired
%     if s0==0
%         dt_p = dr;
%     elseif s0==dr
%         dt_p = 1/aa * mu/Tprimeatce * dr^2 / ( k_old(j) - k_old(j-1) );
%     else
        dt_p = 1/aa * mu/Tprimeatce * 2*dr^2 / ...
            ( 3*k_old(j) - 4*k_old(j-1) + k_old(j-2) );
%     end
    
    %%%%%%%%%%%%%%%%%%%%%%%% solve eqn's with dt_p %%%%%%%%%%%%%%%%%%%%%%%%
    dt_c = 0;
    while abs(dt_p-dt_c)>=tol
        [PO2,thickness,width_retina] = oxygen(tcurr + dt_p,r);
        
        [q1_hat,q2_hat] = growthfactors_implicit(q1_old,q2_old,dt_p,tcurr,...
            r,D1,D2,xi1,xi2,gamma3,gamma4,thickness,width_retina);
        
%         ve_old = ve_calc(j,tcurr,r,c1_old,c2_old,Pm,alpha1,alpha2,gamma1,gamma2,ce);
        
        k_hat = cellpops_sum_withgrowthfactors(j,c1_old,c2_old,q1_hat,q2_hat,...
            PO2,dt_p,r,Pm,kappa,mu,alpha1,alpha2,gamma1,gamma2,cmin,rbar,ce);
    
        [c1_hat,c2_hat] = cellpops_separate_withgrowthfactors(j,c1_old,...
            c2_old,k_hat,q1_hat,q2_hat,PO2,dt_p,r,Pm,kappa,mu,alpha1,alpha2,...
            beta,gamma1,gamma2,cmin,rbar,ce);
        
        %%%%%%%%%%%%%%%%%%%%%%%%% corrector step %%%%%%%%%%%%%%%%%%%%%%%%%%
        bb = 1/2;
        
        %%% should have s0 start further than 2 nodes in
        %%% lines commented so code runs faster, but left in case they're
        %%% desired
%         if s0==0
%             dt_c = mu/Tprimeatce * dr / ( ...
%                 bb*( k_hat(j+1) - k_hat(j) )/dr ...
%                 + (1-bb)*( 1/aa*mu/Tprimeatce ) );
%         elseif s0==dr
%             dt_c = mu/Tprimeatce * dr^2 / ( ...
%                 bb*( 3*k_hat(j+1) - 4*k_hat(j) + k_hat(j-1) )/2 ...
%                 + (1-bb)*( k_old(j) - k_old(j-1) ) );
%         else
            dt_c = mu/Tprimeatce * 2*dr^2 / ( ...
                bb*( 3*k_hat(j+1) - 4*k_hat(j) + k_hat(j-1) ) ...
                + (1-bb)*( 3*k_old(j) - 4*k_old(j-1) + k_old(j-2) ) );
%         end
        
%         [dt_p, dt_c]
%         keyboard
        
        if abs(dt_p-dt_c)<tol
            break;
        else
            dt_p = dt_c;
            dt_c = 0;
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%% solve next time step %%%%%%%%%%%%%%%%%%%%%%%%%
    whatstep = 'corrector';
    
    [PO2,thickness] = oxygen(tcurr + dt_c,r);
    
    [q1_new,q2_new] = growthfactors_implicit(q1_old,q2_old,dt_c,tcurr,r,...
        D1,D2,xi1,xi2,gamma3,gamma4,thickness,width_retina);
    
    k_new = cellpops_sum_withgrowthfactors(j,c1_old,c2_old,q1_new,q2_new,...
        PO2,dt_c,r,Pm,kappa,mu,alpha1,alpha2,gamma1,gamma2,cmin,rbar,ce);
    
    [c1_new,c2_new] = cellpops_separate_withgrowthfactors(j,c1_old,...
        c2_old,k_new,q1_new,q2_new,PO2,dt_c,r,Pm,kappa,mu,alpha1,alpha2,beta,...
        gamma1,gamma2,cmin,rbar,ce);

    %%%%%%%%%%%%%%%%%%%%%% reset for next time step %%%%%%%%%%%%%%%%%%%%%%%
    j = j+1;
    s = s + dr;
    tcurr = tcurr + dt_c;
    c1_old = c1_new;
    c2_old = c2_new;
    k_old = k_new;
    q1_old = q1_new;
    q2_old = q2_new;
    
    %%% save variables
    mvgbdy = [mvgbdy ; s];
    c1 = [c1 ; c1_new];
    c2 = [c2 ; c2_new];
    k = [k ; k_new];
    q1 = [q1 ; q1_new];
    q2 = [q2 ; q2_new];
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