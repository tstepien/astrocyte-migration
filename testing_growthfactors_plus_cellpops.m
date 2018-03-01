clear variables global;
clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% astrocyte parameters
kappa = 1;
mu = 0.01;
alpha1 = 0;%0.028;
alpha2 = 0;%0.010;
beta = 0;%0.0525;
gamma1 = 0;%0;
gamma2 = 0;%0.5;

%%% growth factor parameters
D1 = 1;
D2 = 1;
eta1 = 0.1;
eta2 = 0.1;

%%% tension parameters
rbar = 7.75*10^(-3);
rproc = 15.5*10^(-3);
cmin = 1/(pi*rproc^2);

%%% moving boundary condition parameters
Te = 0.0035; % tension on boundary
ce = densityatbdy(Te,kappa,cmin,rbar); % c1+c2 on boundary

Tprimeatce = Tderivative(ce,kappa,cmin,rbar); % T'(ce)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% mesh %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dr = 0.05;
dt = 0.001;

%check CFL
if dt >=dr^2/(2*max(D1,D2))
    disp('check CFL')
end

rmax = 5;
tmax = 1;

r = 0:dr:rmax;
R = length(r);
rhalf = (r(1:R-1)+r(2:R))/2;

%%% time
t = 0:dt:tmax;

%%%%%%%%%%%%%%%%%%%%%%%%%%% initial conditions %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% moving boundary location
j = 11;
s = (j-1)*dr;

if abs(s/dr - round(s/dr))>2*eps
    error('error specifying s(0): moving boundary must be located on grid node');
end

%%% astrocytes
c1_old = zeros(1,R);
for i=1:R
    xval2 = s;
    yval = 1.5*ce;
    if r(i)<=xval2 %fitted parabola
        c1_old(i) = (ce-yval)/s^2*r(i)^2 + yval;
    else
        c1_old(i) = 0;
    end
end
c2_old = zeros(1,R);

c1 = zeros(length(t),length(r));
c2 = zeros(length(t),length(r));

c1(1,:) = c1_old;
c2(1,:) = c2_old;

%%% growth factors
% p1_old = 100*ones(size(r));
p1_old = smooth(100*( heaviside(1-r)+(exp(-r+1)).*heaviside(r-1) ))';
% p2_old = 100*ones(size(r));
p2_old = smooth(100*heaviside(1-r))';

p1 = zeros(length(t),length(r));
p2 = zeros(length(t),length(r));

p1(1,:) = p1_old;
p2(1,:) = p2_old;

for i = 2:length(t)
    [p1(i,:),p2(i,:)] = growthfactors_implicit(p1(i-1,:),...
        p2(i-1,:),dt,r,D1,D2,eta1,eta2);
    
    [c1(i,:),c2(i,:)] = cellpops_implicit(j,c1(i-1,:),c2(i-1,:),...
        p1(i,:),p2(i,:),dt,r,kappa,mu,alpha1,alpha2,beta,gamma1,gamma2,...
        cmin,rbar,ce);
    
    j = j+1;
    if j==101
        break;
    end
end

if R<length(t)
    ind_end = find(c1(:,1)==0,1);
else
    ind_end = length(t);
end

tplot = 1:floor(ind_end/10):ind_end;
figure
subplot(2,2,1)
hold on
for i=tplot
    plot(r,c1(i,:))
end
hold off
title('c1 - Implicit Method')

subplot(2,2,2)
hold on
for i=tplot
    plot(r,p1(i,:))
end
hold off
title('p1 - Implicit Method')

subplot(2,2,3)
hold on
for i=tplot
    plot(r,c2(i,:))
end
hold off
title('c2 - Implicit Method')

subplot(2,2,4)
hold on
for i=tplot
    plot(r,p2(i,:))
end
hold off
title('p2 - Implicit Method')