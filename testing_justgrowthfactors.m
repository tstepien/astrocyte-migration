clear variables global;
clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% growth factor parameters
D1 = 1;
D2 = 1;
eta1 = 0.1;
eta2 = 0.1;

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

%%% growth factors
% p1_old = 100*ones(size(r));
p1_old = smooth(100*( heaviside(1-r)+(exp(-r+1)).*heaviside(r-1) ))';
% p2_old = 100*ones(size(r));
p2_old = smooth(100*heaviside(1-r))';

t = 0:dt:tmax;

p1_exp = zeros(length(t),length(r));
p2_exp = zeros(length(t),length(r));
p1_imp = zeros(length(t),length(r));
p2_imp = zeros(length(t),length(r));

p1_exp(1,:) = p1_old;
p2_exp(1,:) = p2_old;
p1_imp(1,:) = p1_old;
p2_imp(1,:) = p2_old;

for i = 2:length(t)
    [p1_exp(i,:),p2_exp(i,:)] = growthfactors_explicit(p1_exp(i-1,:),...
        p2_exp(i-1,:),dt,r,D1,D2,eta1,eta2);
    [p1_imp(i,:),p2_imp(i,:)] = growthfactors_implicit(p1_imp(i-1,:),...
        p2_imp(i-1,:),dt,r,D1,D2,eta1,eta2);
end

diffp1 = abs(p1_exp-p1_imp);
diffp2 = abs(p2_exp-p2_imp);


tplot = 1:floor(length(t)/10):length(t);
figure
subplot(2,3,1)
hold on
for i=tplot
    plot(r,p1_exp(i,:))
end
hold off
title('p1 - Explicit Method')

subplot(2,3,2)
hold on
for i=tplot
    plot(r,p1_imp(i,:))
end
hold off
title('p1 - Implicit Method')

subplot(2,3,3)
hold on
for i=tplot
    plot(r,diffp1(i,:))
end
hold off
title('p1 - Difference between methods')

subplot(2,3,4)
hold on
for i=tplot
    plot(r,p2_exp(i,:))
end
hold off
title('p2 - Explicit Method')

subplot(2,3,5)
hold on
for i=tplot
    plot(r,p2_imp(i,:))
end
hold off
title('p2 - Implicit Method')

subplot(2,3,6)
hold on
for i=tplot
    plot(r,diffp2(i,:))
end
hold off
title('p2 - Difference between methods')