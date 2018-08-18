clear variables global;
% clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% growth factor parameters
lambda = 1.6; %%% tortuosity of medium (unitless)
phi = 0.2; %%% porosity/volume fraction in extracellular space (%)
Dwater_PDGFA = 1.2*10^(-6) *(60^2*10^2); %%% diffusion of PDGFA in water at 37C (cm^2/s)
                                     %%% but then converted to (mm^2/hr)
Dwater_LIF = 1.38*10^(-6) *(60^2*10^2); %%% diffusion of LIF in water at 37C (cm^2/s)
                                    %%% but then converted to (mm^2/hr)
xibar_PDGFA = 0.02;
xibar_LIF = 0.015;

D1 = Dwater_PDGFA / lambda^2; %%% effective diffusivity of PDGFA
D2 = Dwater_LIF / lambda^2; %%% effective diffusivity of LIF
xi1 = xibar_PDGFA / phi; %%% production/release rate of PDGFA
xi2 = xibar_LIF / phi; %%% production/release rate of LIF

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% mesh %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dr = 0.1;
dt = 0.1;

%check CFL
if dt >=dr^2/(2*D1)%dr^2/(2*max(D1,D2))
    disp('check CFL')
end

rmax = 5;
tmax = 1*24;

r = 0:dr:rmax;
R = length(r);
rhalf = (r(1:R-1)+r(2:R))/2;

%%% growth factors initial condition
% q1_old = 100*ones(size(r));
% q1_old = smooth(100*( heaviside(1-r)+(exp(-r+1)).*heaviside(r-1) ))';

% sympref('HeavisideAtOrigin', 1); %%% set Heaviside at origin to be 1 
%                                  %%% instead of MATLAB's default 1/2
% initwidth_retina = 1029.17*0.001;
% q1_old = q1max * (-heaviside(r-initwidth_retina)+1);
q1_old = zeros(1,R);

% q2_old = 100*ones(size(r));
% q2_old = smooth(100*heaviside(1-r))';
q2_old = zeros(1,R);

t = 0:dt:tmax;

q1_exp = zeros(length(t),length(r));
q2_exp = zeros(length(t),length(r));
q1_imp = zeros(length(t),length(r));
q2_imp = zeros(length(t),length(r));

q1_exp(1,:) = q1_old;
q2_exp(1,:) = q2_old;
q1_imp(1,:) = q1_old;
q2_imp(1,:) = q2_old;

for i = 2:length(t)
    tday = t(i)/24;
    %%% thickness at edge of retina (convert from micron to mm)
    thickness_peripheral = (13.77 * tday + 72.8) * 0.001;
    %%% thickness at center of retina (optic nerve head) (convert from
    %%% micron to mm)
    thickness_origin = (14.33 * tday + 98.78) * 0.001;
    %%% width/radius of retina (convert from micron to mm)
    width_retina = (414.17 * tday + 1029.17) * 0.001;
    %%% parabolic nonuniform thickness throughout retina
    thickness = (( thickness_peripheral - thickness_origin )./width_retina.^2 .*r.^2 ...
        + thickness_origin ) .* (r<=width_retina);
    
%     [q1_exp(i,:),q2_exp(i,:)] = growthfactors_explicit(q1_exp(i-1,:),...
%         q2_exp(i-1,:),dt,r,D1,D2,xi1,xi2);
    [q1_imp(i,:),q2_imp(i,:)] = growthfactors_implicit(q1_imp(i-1,:),...
        q2_imp(i-1,:),dt,r,D1,D2,xi1,xi2,thickness);
end

diffq1 = abs(q1_exp-q1_imp);
diffq2 = abs(q2_exp-q2_imp);


tplot = 1:floor(length(t)/10):length(t); %length(t);

% max(q1_imp(tplot,:))


figure
% subplot(2,3,1)
% hold on
% for i=tplot
%     plot(r,q1_exp(i,:))
% end
% hold off
% title('q1 - Explicit Method')

% subplot(2,3,2)
subplot(1,2,1)
hold on
for i=tplot
    plot(r,q1_imp(i,:))
end
hold off
title('q1 - Implicit Method')

% subplot(2,3,3)
% hold on
% for i=tplot
%     plot(r,diffq1(i,:))
% end
% hold off
% title('q1 - Difference between methods')

% subplot(2,3,4)
% hold on
% for i=tplot
%     plot(r,q2_exp(i,:))
% end
% hold off
% title('q2 - Explicit Method')
% 
% subplot(2,3,5)
subplot(1,2,2)
hold on
for i=tplot
    plot(r,q2_imp(i,:))
end
hold off
title('q2 - Implicit Method')
% 
% subplot(2,3,6)
% hold on
% for i=tplot
%     plot(r,diffq2(i,:))
% end
% hold off
% title('q2 - Difference between methods')