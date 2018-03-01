clear variables global;
clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% astrocyte parameters
kappa = 1;
mu = 0.1;
alpha = 0.001;%0.01;%0.028;
beta = 0.01;%0.0525;  %%% bigger than alpha
gamma = 0;%0;

%%% growth factor parameters
D1 = 0;
D2 = 0;
eta1 = 0;%0.1;
eta2 = 0;%0.1;

%%% tension parameters
rbar = 7.75*10^(-3);
rproc = 15.5*10^(-3);
cmin = 1/(pi*rproc^2);

%%% moving boundary condition parameters
Te = 0.0035; % tension on boundary
ce = densityatbdy(Te,kappa,cmin,rbar); % c1+c2 on boundary

Tprimeatce = Tderivative(ce,kappa,cmin,rbar); % T'(ce)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% mesh %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dr = 0.01;

rmax = 5;
tmax = 15;

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

c_old = c1_init + c2_init;

%%% growth factors
p1_old = 100*ones(1,R);%smooth(100*( heaviside(0.25-r)+(exp(-r+0.25)).*heaviside(r-0.25) ))';
p2_old = 100*ones(1,R);%smooth(100*heaviside(0.3-r))';

%%%%%%%%%%%%%%%%%%%%%%% initialize final variables %%%%%%%%%%%%%%%%%%%%%%%%
mvgbdy = s;
c = c_old;
p1 = p1_old;
p2 = p2_old;
t = tcurr;
% Tension = [];
% vhalf = [];

%%% subscript i is for space, j is for time (write in the order (j,i)) 

j_init = s/dr+1;
j = j_init;

while tcurr < tmax && j<R-1
    %%%%%%%%%%%%%%%%%%%%%%%%%%% predictor step %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    aa = 1;
    if s==0
        dt_p = dr;
    elseif s==dr
        dt_p = 1/aa * mu/Tprimeatce * dr^2 / ( c_old(j) - c_old(j-1) );
    else
%         disp('new dt_p')
        dt_p = 1/aa * mu/Tprimeatce * 2*dr^2 / ...
            ( 3*c_old(j) - 4*c_old(j-1) + c_old(j-2) );
    end
%     keyboard
    
    dt_c = 0;
    while abs(dt_p-dt_c)>=tol        
        [p1_hat,p2_hat] = growthfactors_implicit(p1_old,p2_old,dt_p,r,D1,...
            D2,eta1,eta2);
    
        c_hat = cellpops_implicit(j,c_old,p1_hat,dt_p,r,kappa,mu,alpha,...
            gamma,cmin,rbar,ce);
        
%         if c1_hat==c1_old & c2_hat==c2_old
%             continue;
%         end
        
%         if abs(vhalf_hat(j-1)-vhalf_hat(j))<eps
%             break;
%         end
        
        %%%%%%%%%%%%%%%%%%%%%%%%% corrector step %%%%%%%%%%%%%%%%%%%%%%%%%%
        bb = 1/2;
        if s==0
            dt_c = mu/Tprimeatce * dr / ( ...
                bb*( c_hat(j+1) - c_hat(j) )/dr ...
                + (1-bb)*( 1/aa*mu/Tprimeatce ) );
        elseif s==dr
            dt_c = mu/Tprimeatce * dr^2 / ( ...
                bb*( 3*c_hat(j+1) - 4*c_hat(j) + c_hat(j-1) )/2 ...
                + (1-bb)*( c_old(j) - c_old(j-1) ) );
        else
            dt_c = mu/Tprimeatce * 2*dr^2 / ( ...
                bb*( 3*c_hat(j+1) - 4*c_hat(j) + c_hat(j-1) ) ...
                + (1-bb)*( 3*c_old(j) - 4*c_old(j-1) + c_old(j-2) ) );
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
    
    c_new = cellpops_implicit(j,c_old,p1_new,dt_c,r,kappa,mu,alpha,...
        gamma,cmin,rbar,ce);
        
    %%% reset for next time step
    j = j+1;
    s = s + dr;
    tcurr = tcurr + dt_c;
    c_old = c_new;
    p1_old = p1_new;
    p2_old = p2_new;
    
    %%% save variables
    mvgbdy = [mvgbdy ; s];
    c = [c ; c_new];
    p1 = [p1 ; p1_new];
    p2 = [p2 ; p2_new];
    t = [t ; tcurr];
%     Tension = [Tension ; Tension_old];
%     vhalf = [vhalf ; vhalf_old];
end
    
%%%%%%%%%%%%%%%%%%%% solve ODE for c2 - c1 populations %%%%%%%%%%%%%%%%%%%%
[k,velplot] = cellpopsdiff_upwind(j_init,c1_init,c2_init,c,r,t,p1,p2,kappa,mu,alpha,beta,gamma,cmin,rbar);

c1 = (c-k)/2;
c2 = (c+k)/2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% plotting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T = length(t);
numcurvesplot = 10;
if T<=numcurvesplot
    ind = 1:T;
else
    ind = 1:floor(T/numcurvesplot):T;
%     ind = 1:10;
end
% tplot = t;
% [tf,ind] = ismember(tplot,t);
% c1plot = c1(ind,:);
% c2plot = c2(ind,:);
% p1plot = p1(ind,:);
% p2plot = p2(ind,:);
% Tensionplot = Tension(ind,:);
% vhalfplot = vhalf(ind,:);

xL = min(mvgbdy(end) + 5*dr , r(end));

figure
subaxis(4,2,1,'MarginLeft',0.05,'MarginRight',0.01,'MarginTop',0.03,'MarginBottom',0.05)
plot(r,c(ind,:))
xlabel('r')
ylabel('c1 (APCs) + c2 (IPAs)')
set(gca,'Xlim',[0,xL])

subaxis(4,2,2,'MarginLeft',0.05,'MarginRight',0.01,'MarginTop',0.03,'MarginBottom',0.05)
plot(r,k(ind,:))
xlabel('r')
ylabel('c2 (IPAs) - c1 (APCs)')
set(gca,'Xlim',[0,xL])

subaxis(4,2,3,'MarginLeft',0.05,'MarginRight',0.01,'MarginTop',0.03,'MarginBottom',0.05)
plot(r,c1(ind,:))
xlabel('r')
ylabel('c1 (APCs)')
set(gca,'Xlim',[0,xL])
% ylim(yl);

subaxis(4,2,4,'MarginLeft',0.05,'MarginRight',0.01,'MarginTop',0.03,'MarginBottom',0.05)
plot(r,c2(ind,:))
xlabel('r')
ylabel('c2 (IPAs)')
set(gca,'Xlim',[0,xL])
% ylim(yl);
% yl = ylim;

subaxis(4,2,5,'MarginTop',0.03,'MarginBottom',0.05)
plot(t,mvgbdy)
xlabel('t')
ylabel('moving boundary s(t)')


subaxis(4,2,6,'MarginTop',0.03,'MarginBottom',0.05)
plot(r,velplot(ind,:))
xlabel('r')
ylabel('velocity')
set(gca,'Xlim',[0,xL])

subaxis(4,2,7,'MarginTop',0.03,'MarginBottom',0.05)
plot(r,p1(ind,:))
xlabel('r')
ylabel('p1 (PDGFA)')

subaxis(4,2,8)
plot(r,p2(ind,:))
xlabel('r')
ylabel('p2 (LIF)')

set(gcf,'Units','inches','Position',[2,2,12,8],'PaperPositionMode','auto')
% legend(['t=',num2str(tplot(1))],['t=',num2str(tplot(2))],...
%     ['t=',num2str(tplot(3))],['t=',num2str(tplot(4))],...
%     ['t=',num2str(tplot(5))])


%%%%%%%%%%%%%%%%%%%%%% testing: conservation of mass %%%%%%%%%%%%%%%%%%%%%%
% area under curve - trapezoidal rule
areaundercurve = zeros(length(t),1);
for i = 1:length(t)
    areaundercurve(i) = dr*sum( c(1:mvgbdy(1)/dr+(i-1)) ) / 2;
end