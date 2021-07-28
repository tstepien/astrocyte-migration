clear variables;
clc;
close all;
addpath ../

parameters_fixed;

rmax = 5; %%% max radius (mm) (estimate rat retinal radius = 4.1 mm)
tmax = 7*24; %%% max time (hr) (7 days = 168 hr)

nxpts = 26;
ntpts = 29;

x = linspace(0,1,nxpts);    % 0 <= x <= 1 by definition
t = linspace(0,tmax,ntpts);

% Set up vectors for interpolation to obtain PO2 
[Lvec, Pvec] = oxygen_setup(M0, Dalpha, Pm, P0);

%%% plot times: day 0, 1, 2, 3, 4, 5, 6, 7
numcurvesplot = 8;

% set up retinal thickness and choroid oxygen
r = zeros(numcurvesplot,nxpts+2);
PO2 = zeros(numcurvesplot,nxpts+2);

for i=1:numcurvesplot
    [~,~,~,radius_ret] = thick_rad(t(i*4-3),0);
    for j = 1:nxpts
        r(i,j) = radius_ret * x(j);
        [thickness_ret,~,~,~] = thick_rad(t(i*4-3),r(i,j));
        PO2(i,j) = interp1(Lvec,Pvec,thickness_ret);
    end
    r(i,nxpts+1) = r(i,nxpts);   %to make plots go to axis at outer edge 
    r(i,nxpts+2) = rmax;  
end

%%% color order
co = [0    0.4470    0.7410
    0.8500    0.3250    0.0980
    0.9290    0.6940    0.1250
    0.4940    0.1840    0.5560
    0.4660    0.6740    0.1880
    0.3010    0.7450    0.9330
    0.6350    0.0780    0.1840
    0         0.3725    0];

fslabel = 22;
fsticks = 18;
fslegend = 14;

figure
hold on
for i=1:numcurvesplot
    plot(r(i,:),PO2(i,:),'LineWidth',1.5,'Color',co(i,:)) 
end
hold off
xlabel('Radius (mm)','Interpreter','latex','FontSize',fslabel)
ylabel('Choroid\, $\mathrm{P}_{\mathrm{O}_2}$ (mmHg)','Interpreter','latex','FontSize',fslabel)
box on
set(gca,'XLim',[0,rmax],'FontSize',fsticks,'Position',[0.14 0.14 0.82 0.76])

h = legend('E15','E16','E17','E18','E19','E20','E21','E22/P0');
set(h,'FontSize',fslegend);

title('B                                                ',...
    'FontSize',26)

set(gcf,'Units','inches','Position',[2 2 7.75 5.75],'PaperPositionMode','auto')