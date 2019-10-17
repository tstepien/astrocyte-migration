clear variables global;
% clc;

%%%%%%%%%%%%%%%%%%% input all fixed parameters that are %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% known/derived from literature %%%%%%%%%%%%%%%%%%%%%%
parameters_fixed

%%%%%%%%%%%%%%%%%%%%%%%% growth factor parameters %%%%%%%%%%%%%%%%%%%%%%%%%
xibar_PDGFA = 0.000000001;
xibar_LIF = 0.001;

xi1 = xibar_PDGFA / phi; %%% production/release rate of PDGFA
xi2 = xibar_LIF / phi; %%% production/release rate of LIF

%%% degradation rates
quasilength = 0.1;
gamma3 = D1/quasilength^2;
gamma4 = D2/quasilength^2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% mesh %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dr = 0.1;
dt = 0.01;

%check CFL
if dt >=dr^2/(2*max(D1,D2))
    disp('check CFL')
end

rmax = 5; %%% max radius (mm) (estimate rat retinal radius = 4.1 mm)
tmax = 7*24; %%% max time (hr) (7 days = 168 hr)

r = 0:dr:rmax;
R = length(r);
rhalf = (r(1:R-1)+r(2:R))/2;


%%%%%%%%%%%%%%%%%%%%%%%%%%% initial conditions %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% growth factors initial condition
q1_old = zeros(1,R);
q2_old = zeros(1,R);

%%% initialize variables
t = 0:dt:tmax;

q1_exp = zeros(length(t),length(r));
q2_exp = zeros(length(t),length(r));
q1_imp = zeros(length(t),length(r));
q2_imp = zeros(length(t),length(r));

q1_exp(1,:) = q1_old;
q2_exp(1,:) = q2_old;
q1_imp(1,:) = q1_old;
q2_imp(1,:) = q2_old;

tcurr = 0;

for i = 2:length(t)
    [PO2,thickness,width_retina] = oxygen(tcurr + dt,r);
    [thickness_RGC,radius_endo] = thick_rad(dt,tcurr,r,width_retina);
    
%     [q1_exp(i,:),q2_exp(i,:)] = growthfactors_explicit(q1_exp(i-1,:),...
%         q2_exp(i-1,:),dt,r,D1,D2,xi1,xi2);
    [q1_imp(i,:),q2_imp(i,:)] = growthfactors_implicit(q1_imp(i-1,:),...
        q2_imp(i-1,:),dt,tcurr,r,dr,R,thickness_RGC,radius_endo,...
        maxRGCthick,thickness,D1,D2,xi1,xi2,gamma3,gamma4);
    
    
    tcurr = tcurr + dt;
end

diffq1 = abs(q1_exp-q1_imp);
diffq2 = abs(q2_exp-q2_imp);



%% plots
%%% plot times: day 0, 1, 2, 3, 4, 5, 6, 7
numcurvesplot = 8;
tplot = zeros(1,numcurvesplot);
for i=0:numcurvesplot-1
    tplot(i+1) = find(abs((t/24-i))==min(abs(t/24-i)));
end


figure
% subplot(2,3,1)
% hold on
% for i=tplot
%     plot(r,q1_exp(i,:))
% end
% hold off
% title('PDGFA - Explicit Method')

% subplot(2,3,2)
subplot(1,2,1)
hold on
for i=tplot
    plot(r,q1_imp(i,:))
end
hold off
title('PDGFA - Implicit Method')

% subplot(2,3,3)
% hold on
% for i=tplot
%     plot(r,diffq1(i,:))
% end
% hold off
% title('PDGFA - Difference between methods')

% subplot(2,3,4)
% hold on
% for i=tplot
%     plot(r,q2_exp(i,:))
% end
% hold off
% title('LIF - Explicit Method')
% 
% subplot(2,3,5)
subplot(1,2,2)
hold on
for i=tplot
    plot(r,q2_imp(i,:))
end
hold off
title('LIF - Implicit Method')
% 
% subplot(2,3,6)
% hold on
% for i=tplot
%     plot(r,diffq2(i,:))
% end
% hold off
% title('LIF - Difference between methods')

fsticks = 14;
h = legend('0 days (E15)','1 day (E16)','2 days (E17)','3 days (E18)',...
    '4 days (E19)','5 days (E20)','6 days (E21)','7 days (E22/P0)');
set(h,'FontSize',fsticks,'Position',[0.85 0.57827 0.0833 0.372]);