clear variables global;
% clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% growth factor parameters
lambda = 1.6; %%% tortuosity of medium (unitless)
phi = 0.2; %%% porosity/volume fraction in extracellular space (%)
Dwater_PDGFA = 1.2*10^(-6) *(60*60*10^2); %%% diffusion of PDGFA in water at 37C (cm^2/s)
                                     %%% but then converted to (mm^2/hr)
eta_PDGFA = 0.02;

D1 = Dwater_PDGFA / lambda^2; %%% effective diffusivity of PDGFA
D2 = 1;
eta1 = eta_PDGFA / phi; %%% production/release rate of PDGFA
eta2 = 0.1;

p1max = 0.0001;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% mesh %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dr = 0.1;
dt = 0.1;

%check CFL
if dt >=dr^2/(2*D1)%dr^2/(2*max(D1,D2))
    disp('check CFL')
end

rmax = 5;
tmax = 7*24;

r = 0:dr:rmax;
R = length(r);
rhalf = (r(1:R-1)+r(2:R))/2;

%%% growth factors
% p1_old = 100*ones(size(r));
% p1_old = smooth(100*( heaviside(1-r)+(exp(-r+1)).*heaviside(r-1) ))';
sympref('HeavisideAtOrigin', 1); %%% set Heaviside at origin to be 1 
                                 %%% instead of MATLAB's default 1/2
initwidth_retina = 1029.17*0.001;
p1_old = p1max * (-heaviside(r-initwidth_retina)+1);

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
%     [p1_exp(i,:),p2_exp(i,:)] = growthfactors_explicit(p1_exp(i-1,:),...
%         p2_exp(i-1,:),dt,r,D1,D2,eta1,eta2);
    [p1_imp(i,:),p2_imp(i,:)] = growthfactors_implicit(p1_imp(i-1,:),...
        p2_imp(i-1,:),dt,r,D1,D2,eta1,eta2);
end

diffp1 = abs(p1_exp-p1_imp);
diffp2 = abs(p2_exp-p2_imp);


tplot = length(t);%1:floor(length(t)/10):length(t);

max(p1_imp(tplot,:))


figure
% subplot(2,3,1)
% hold on
% for i=tplot
%     plot(r,p1_exp(i,:))
% end
% hold off
% title('p1 - Explicit Method')

% subplot(2,3,2)
hold on
for i=tplot
    plot(r,p1_imp(i,:))
end
hold off
title('p1 - Implicit Method')

% subplot(2,3,3)
% hold on
% for i=tplot
%     plot(r,diffp1(i,:))
% end
% hold off
% title('p1 - Difference between methods')

% subplot(2,3,4)
% hold on
% for i=tplot
%     plot(r,p2_exp(i,:))
% end
% hold off
% title('p2 - Explicit Method')
% 
% subplot(2,3,5)
% hold on
% for i=tplot
%     plot(r,p2_imp(i,:))
% end
% hold off
% title('p2 - Implicit Method')
% 
% subplot(2,3,6)
% hold on
% for i=tplot
%     plot(r,diffp2(i,:))
% end
% hold off
% title('p2 - Difference between methods')