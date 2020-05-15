clear variables;
clc;
close all;

kappa = 1;
rbar = 7.75*10^(-3);
rproc = 15.5*10^(-3);
cmin = 1/(pi*rproc^2);
cmax = 1/(pi*rbar^2);

cc = 0:1:5500;

tension = kappa*cc.^2./(cmin^2 + cc.^2).*(1./sqrt(pi*cc) - rbar);
tmax = 0.004;

fslabel = 22;
fsticks = 18;
graycolor = [0.4 0.4 0.4];

figure
hold on
plot(cc,tension,'LineWidth',1.5,'Color','k')
line([cmin,cmin],[0,tmax],'LineStyle','--','LineWidth',1.5,'Color',graycolor)
line([cmax,cmax],[0,tmax],'LineStyle','--','LineWidth',1.5,'Color',graycolor)
hold off
text(cmin-550,0.002,'$c_\mathrm{min}$','Color',graycolor,...
    'FontSize',fslabel,'Interpreter','latex')
text(cmax-550,0.002,'$c_\mathrm{max}$','Color',graycolor,...
    'FontSize',fslabel,'Interpreter','latex')
xlabel('$c_1+c_2$ (cells/mm$^2$)','Interpreter','latex','FontSize',fslabel)
ylabel('$T(c_1+c_2)$ (mN/mm)','Interpreter','latex','FontSize',fslabel)
box on
set(gca,'YLim',[0,tmax],'FontSize',fsticks,'Position',[0.14 0.14 0.82 0.76])

set(gcf,'Units','inches','Position',[2 2 7.75 5.75],'PaperPositionMode','auto')