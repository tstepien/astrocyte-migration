%%% plot times: day 0, 1, 2, 3, 4, 5, 6, 7
numcurvesplot = 8;

% color order
co = [0    0.4470    0.7410
    0.8500    0.3250    0.0980
    0.9290    0.6940    0.1250
    0.4940    0.1840    0.5560
    0.4660    0.6740    0.1880
    0.3010    0.7450    0.9330
    0.6350    0.0780    0.1840
    0         0.3725    0];

fstitle = 18;
fslabel = 18;
fsticks = 16;
fslegend = 14;

figure
tiledlayout(2,3,'TileSpacing','Compact','Padding','none')
% extra points to make plots go to axis at outer edge
r = zeros(numcurvesplot,nxpts+2);
APC_plot = zeros(numcurvesplot,nxpts+2);
IPA_plot = zeros(numcurvesplot,nxpts+2);
tension_plot = zeros(numcurvesplot,nxpts+2);
velocity_plot = zeros(numcurvesplot,nxpts+2);

%%%%%%%%%%%%%%%%%%%%% top row: set up cell densities %%%%%%%%%%%%%%%%%%%%%%
rm = sol(:,nxpts,4);
for i=1:numcurvesplot
    for j = 1:nxpts
        r(i,j) = rm(i*4-3) * x(j);
        APC_plot(i,j) = sol(i*4-3,j,1);
        IPA_plot(i,j) = sol(i*4-3,j,2);        
    end
    r(i,nxpts+1) = r(i,nxpts);   %to make plots go to axis at outer edge 
    r(i,nxpts+2) = rmax;
end

% find max value of APC+IPA and round up to nearest 500
ylim_sum = ceil(max(max(APC_plot+IPA_plot))/500)*500;

%%%%%%%%%%%%%%%% PANEL: c1 + c2 - APCs + IPAs
nexttile
hold on
for i=1:numcurvesplot
    plot(r(i,:),APC_plot(i,:)+IPA_plot(i,:),'LineWidth',2.5,'Color',co(i,:))
end
hold off
box off
xlabel('$r$ (mm)','FontSize',fslabel,'Interpreter','latex')
ylabel('APCs + IPAs (cells/mm$^2$)','FontSize',fslabel,'Interpreter','latex')
set(gca,'XLim',[0,rmax],'YLim',[0,ylim_sum],'FontSize',fsticks)
xticks(0:rmax)

h = legend('E15','E16','E17','E18','E19','E20','E21','E22/P0');
set(h,'FontSize',fslegend,'Position',[0.232 0.658 0.1 0.288]);
title('A                              ',...
    'FontSize',fstitle)

%%%%%%%%%%%%%%%% PANEL: c1 - APCs
nexttile
hold on
for i=1:numcurvesplot
    plot(r(i,:),APC_plot(i,:),'LineWidth',2.5,'Color',co(i,:))
end
hold off
box off
xlabel('$r$ (mm)','FontSize',fslabel,'Interpreter','latex')
ylabel('APCs (cells/mm$^2$)','FontSize',fslabel,'Interpreter','latex')
set(gca,'XLim',[0,rmax],'YLim',[0,ylim_sum],'FontSize',fsticks)
xticks(0:rmax)
title('B                              ',...
    'FontSize',fstitle)

%%%%%%%%%%%%%%%% PANEL: c2 - IPAs
nexttile
hold on
for i=1:numcurvesplot
    plot(r(i,:),IPA_plot(i,:),'LineWidth',2.5,'Color',co(i,:))
end
hold off
box off
xlabel('$r$ (mm)','FontSize',fslabel,'Interpreter','latex')
ylabel('IPAs (cells/mm$^2$)','FontSize',fslabel,'Interpreter','latex')
set(gca,'XLim',[0,rmax],'YLim',[0,ylim_sum],'FontSize',fsticks)
xticks(0:rmax)
title('C                              ',...
    'FontSize',fstitle)

%%%%%%%%%%%%%%%%%%%%%%%% bottom row: set up tension %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% and radial velocity in interior of domain %%%%%%%%%%%%%%%%
%%% tension
for i=1:numcurvesplot
    for j = 1:nxpts
        usumT = sol(i*4-3,j,1) + sol(i*4-3,j,2);
        tension_plot(i,j) = kappa * (1./sqrt(pi*usumT) - rbar/1000) ;
    end
end

%%% radial velocity
for i=1:numcurvesplot
    usumV = sol(i*4-3,:,1) + sol(i*4-3,:,2);
    for j = 1:nxpts
    if j == 1 
        DuDx = (nxpts - 1) * (usumV(j+1) - usumV(j));
    elseif  j == nxpts
        DuDx = (nxpts - 1) * (usumV(j) - usumV(j-1));
    else
        DuDx = (nxpts - 1)/2 * (usumV(j+1) - usumV(j-1));
    end
    velocity_plot(i,j) = -kTprime1 * usumV(j)^(-3/2) / sol(i*4-3,j,4) * DuDx;
    end
end


%%%%%%%%%%%%%%%% PANEL0.215 0.614 0.101 0.332: tension
nexttile

hold on
for i=1:numcurvesplot
    plot(r(i,:),tension_plot(i,:),'LineWidth',2.5,'Color',co(i,:))
end
hold off
xlabel('$r$ (mm)','FontSize',fslabel,'Interpreter','latex')
ylabel('Tension (Pa mm)','FontSize',fslabel,'Interpreter','latex')
set(gca,'XLim',[0,rmax],'FontSize',fsticks)
xticks(0:rmax)
title('D                              ',...
    'FontSize',fstitle)

%%%%%%%%%%%%%%%% PANEL: radial velocity
nexttile
hold on
for i=1:numcurvesplot
    plot(r(i,:),velocity_plot(i,:),'LineWidth',2.5,'Color',co(i,:))
end
hold off
box off
xlabel('$r$ (mm)','FontSize',fslabel,'Interpreter','latex')
ylabel('APC velocity (mm/hr)','FontSize',fslabel,'Interpreter','latex')
set(gca,'XLim',[0,rmax],'YLim',[0 0.06],'FontSize',fsticks)
xticks(0:rmax)
title('E                              ',...
    'FontSize',fstitle)

%%%%%%%%%%%%%%%% PANEL: Radius vs. time / distance astrocytes spread
rad_APC = [0.17; 0.33; 0.5; 0.67; 1.67; 2.17; 2.67]; %not including E17
rad_days = [0; 1; 3; 4; 5; 6; 7];
nexttile(6)
hold on
plot(t/24,rm,'k','LineWidth',2.5)
scatter(rad_days,rad_APC,150,[0.5 0.5 0.5],'x','LineWidth',2.5)
hold off
box off
xlabel('$t$ (days)','FontSize',fslabel,'Interpreter','latex')
ylabel('Cell boundary (mm)','FontSize',fslabel,'Interpreter','latex')
set(gca,'XLim',[0,8],'FontSize',fsticks)
xticks(0:2:8)
title('F                              ',...
    'FontSize',fstitle)

set(gcf,'Units','inches','Position',[1,1,12,7],'PaperPositionMode','auto');
