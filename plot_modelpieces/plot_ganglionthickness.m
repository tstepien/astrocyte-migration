clear variables;
clc;
close all;
addpath ../

parameters_fixed;

dr = 0.01;
rmax = 5; %%% max radius (mm) (estimate rat retinal radius = 4.1 mm)
dt = 1;
tmax = 7*24; %%% max time (hr) (7 days = 168 hr)

r = 0:dr:rmax;
t = (0:dt:tmax)';

%%% cell layer thickness and radius
[~,thickness_RGC,~,~] = thick_rad(t,r);

%%% dimensionalize
thickness_RGC = thickness_RGC * maxRGCthick; % normalized

%%% plot times: day 0, 1, 2, 3, 4, 5, 6, 7
numcurvesplot = 8;
plotind = zeros(1,numcurvesplot);
for i=0:numcurvesplot-1
    plotind(i+1) = find(abs((t/24-i))==min(abs(t/24-i)));
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

fstitle = 26;
fslabel = 26;
fsticks = 22;
fslegend = 18;

figure
hold on
for i=1:numcurvesplot
    plot(r,thickness_RGC(plotind(i),:),'LineWidth',2.5,'Color',co(i,:))
end
hold off
xlabel('$r$ (mm)','Interpreter','latex','FontSize',fslabel)
ylabel('Ganglion cell layer thickness (mm)','Interpreter','latex','FontSize',fslabel)
box off
set(gca,'XLim',[0,rmax],'FontSize',fsticks,'Position',[0.175 0.14 0.79 0.76])

h = legend('E15','E16','E17','E18','E19','E20','E21','E22/P0');
set(h,'FontSize',fslegend,'Position',[0.769 0.467 0.181 0.434]);

title('A                                                ',...
    'FontSize',fstitle)

set(gcf,'Units','inches','Position',[2 2 7.75 5.75],'PaperPositionMode','auto')