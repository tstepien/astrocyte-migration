clear variables global;
clc;
% close all;

global whatstep;

%%% time unit: hr
%%% space unit: mm

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% oxygen parameters
Pm = 10; % (mmHg)

%%% astrocyte parameters
kappa = 1;
mu = 0.1;
alpha1 = 0.08; %%% (/hr)
alpha2 = 0.06; %%% (/hr)
beta = 0.2; %%% (/hr)
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

rmax = 5; %%% max radius (mm) (estimate rat retinal radius = 4.1 mm)
tmax = 0.1*24; %%% max time (hr) (7 days = 168 hr)

r = 0:dr:rmax;
R = length(r);
rhalf = (r(1:R-1)+r(2:R))/2;

tol = 10^(-6); % tolerance for predictor-corrector scheme

%%%%%%%%%%%%%%%%%%%%%%%%%%% initial conditions %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% moving boundary location
s = 0.5;%10*dr;

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
mvgbdy = s;
c1 = c1_init;
c2 = c2_init;
t = tcurr;

%%% subscript i is for space, j is for time (write in the order (j,i)) 

j_init = s/dr+1;
j = j_init;

%%% velocity
[vel_cir,vel_rad] = velocity(j,c1,c2,r,kappa,cmin,rbar,mu);

%%% concentration on the moving boundary (mb)
c1mb = c1(j);
c2mb = c2(j);

while tcurr < tmax && j<R-1
    %%%%%%%%%%%%%%%%%%%%%%%%%%% predictor step %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    whatstep = 'predictor';
    
    aa = 1;
    k_old = c1_old + c2_old;
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
    
%         [c1_hat,c2_hat] = cellpopulations_may14(j,c1_old,c2_old,PO2,dt_p,r,Pm,...
%             kappa,mu,alpha1,alpha2,beta,gamma1,gamma2,cmin,rbar,ce);
        [c1_hat,c2_hat] = cellpopulations(j,c1_old,c2_old,k_hat,PO2,dt_p,...
            r,Pm,kappa,mu,alpha1,alpha2,beta,gamma1,gamma2,cmin,rbar,ce);
%         [c1_hat,c2_hat] = cellpopulations_vecalc_trapez(j,c1_old,c2_old,k_hat,PO2,dt_p,...
%             r,ve_old,Pm,kappa,mu,alpha1,alpha2,beta,gamma1,gamma2,cmin,rbar,ce);
%         stop
        %%%%%%%%%%%%%%%%%%%%%%%%% corrector step %%%%%%%%%%%%%%%%%%%%%%%%%%
        
        bb = 1;%1/2;
%         k_hat = c1_hat + c2_hat;
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
    
    [c1_new,c2_new] = cellpopulations(j,c1_old,c2_old,k_new,PO2,dt_c,...
        r,Pm,kappa,mu,alpha1,alpha2,beta,gamma1,gamma2,cmin,rbar,ce);
    
    %%%%%%%%%%%%%%%%%%%%%% reset for next time step %%%%%%%%%%%%%%%%%%%%%%%
    j = j+1;
    s = s + dr;
    tcurr = tcurr + dt_c;
    c1_old = c1_new;
    c2_old = c2_new;
    
    %%% save variables
    mvgbdy = [mvgbdy ; s];
    c1 = [c1 ; c1_new];
    c2 = [c2 ; c2_new];
    t = [t ; tcurr];
    c1mb = [c1mb ; c1_new(j)];
    c2mb = [c2mb ; c2_new(j)];
    
    %%% velocity calculation
    [vel_cir_new,vel_rad_new] = velocity(j,c1_new,c2_new,r,kappa,cmin,rbar,mu);
    vel_cir = [vel_cir ; vel_cir_new];
    vel_rad = [vel_rad ; vel_rad_new];
end

mvgbdy(end)
PO2(end)
t(end)/24

mvgbdy_vel = (mvgbdy(2:end) - mvgbdy(1:end-1)) ./ (t(2:end)-t(1:end-1));
mvgbdy_vel(end)

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% plotting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[PO2,thickness] = oxygen(t,r);

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

%%% make moving boundary sharp
c1plot = zeros(T,R+1);
c2plot = zeros(T,R+1);
rplot = zeros(T,R+1);
for i = 1:T
    c1plot(i,:) = [c1(i,1:j_init+(i-1)) , zeros(1,R+1-(j_init+(i-1)))];
    c2plot(i,:) = [c2(i,1:j_init+(i-1)) , zeros(1,R+1-(j_init+(i-1)))];
    
    rplot(i,:) = [r(1:j_init+(i-1)) , r(j_init+(i-1)) , r(j_init+i:end)];
end

fslabel = 16;
fsticks = 14;

figure
subaxis(2,3,1,'MarginLeft',0.055,'MarginRight',0.01,'MarginTop',0.03,'MarginBottom',0.15)
hold on
for i=1:T
    plot(rplot(i,ind),c1plot(i,ind))
end
hold off
xlabel('radius r (mm)','FontSize',fslabel)
ylabel('APCs (cells/mm^2)','FontSize',fslabel)
set(gca,'XLim',[0,mvgbdy(end)+5*dr],'FontSize',fsticks)
ylims_c1 = get(gca,'YLim');

subaxis(2,3,2,'MarginLeft',0.06,'MarginRight',0.01)
hold on
for i=1:T
    plot(rplot(i,ind),c2plot(i,ind))
end
hold off
xlabel('radius r (mm)','FontSize',fslabel)
ylabel('IPAs (cells/mm^2)','FontSize',fslabel)
set(gca,'XLim',[0,mvgbdy(end)+5*dr],'FontSize',fsticks)%,'YLim',ylims_c1)

subaxis(2,3,3,'MarginLeft',0.06,'MarginRight',0.01)
hold on
for i=1:T
    plot(rplot(i,ind),c1plot(i,ind)+c2plot(i,ind))%,'Color',co(i,:))
end
% for i=1:T
%     plot(rplot(i,:),c2plot(i,:),'--')%,'Color',co(i,:))
% end
hold off
xlabel('radius r (mm)','FontSize',fslabel)
% ylabel('Solid: c1 (APCs), Dashed: c2 (IPAs)')
ylabel('APCs + IPAs','FontSize',fslabel)
set(gca,'XLim',[0,mvgbdy(end)+5*dr],'FontSize',fsticks)

subaxis(2,3,4,'MarginTop',0.03,'MarginBottom',0.08)
plot(t/24,mvgbdy,'-o')
xlabel('t (days)','FontSize',fslabel)
ylabel('moving boundary (mm)','FontSize',fslabel)
set(gca,'FontSize',fsticks)

subaxis(2,3,5,'MarginTop',0.03,'MarginBottom',0.08)
if size(thickness,1)==1 || size(thickness,2)==1
    plot(t/24,thickness,'-o')
    xlabel('t (days)','FontSize',fslabel)
else
    plot(r,thickness)
    xlabel('radius r (mm)','FontSize',fslabel)
end
ylabel('average total retinal thickness (mm)','FontSize',fslabel)
set(gca,'FontSize',fsticks)

subaxis(2,3,6)
if size(thickness,1)==1 || size(thickness,2)==1
    plot(thickness,PO2,'-o')
    xlabel('average total retinal thickness (mm)','FontSize',fslabel)
else
    plot(r,PO2)
    xlabel('radius r (mm)','FontSize',fslabel)
end
ylabel('PO_2 (mmHg)','FontSize',fslabel)
set(gca,'FontSize',fsticks)

set(gcf,'Units','inches','Position',[2,2,16,8],'PaperPositionMode','auto')
% legend(['t=',num2str(tplot(1))],['t=',num2str(tplot(2))],...
%     ['t=',num2str(tplot(3))],['t=',num2str(tplot(4))],...
%     ['t=',num2str(tplot(5))])
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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