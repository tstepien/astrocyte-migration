clear variables;
clc;
close all;
addpath ../

dr = 0.01;
rmax = 5; %%% max radius (mm) (estimate rat retinal radius = 4.1 mm)
dt = 1;
tmax = 7*24; %%% max time (hr) (7 days = 168 hr)

r = 0:dr:rmax;
t = (0:dt:tmax)';

%%% cell layer thickness and radius
[thickness_ret,~,~,~] = thick_rad(t,r);

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

fslabel = 22;
fsticks = 18;
fslegend = 14;

figure
hold on
for i=1:numcurvesplot
    plot(r,thickness_ret(plotind(i),:),'LineWidth',1.5,'Color',co(i,:))
end
hold off
xlabel('Radius (mm)','Interpreter','latex','FontSize',fslabel)
ylabel('Retinal Thickness (mm)','Interpreter','latex','FontSize',fslabel)
box on
set(gca,'XLim',[0,rmax],'FontSize',fsticks,'Position',[0.14 0.14 0.82 0.76])

h = legend('E15','E16','E17','E18','E19','E20','E21','E22/P0');
set(h,'FontSize',fslegend);

title('C                                                ',...
    'FontSize',26)

set(gcf,'Units','inches','Position',[2 2 7.75 5.75],'PaperPositionMode','auto')