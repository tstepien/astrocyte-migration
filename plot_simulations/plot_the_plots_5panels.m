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

fslabel = 16;
fsticks = 14;
fslegend = 14;

figure(7);
tiledlayout(2,3,'TileSpacing','Compact','Padding','compact')
% extra points to make plots go to axis at outer edge
r = zeros(numcurvesplot,nxpts+1);
y1 = zeros(numcurvesplot,nxpts+1);
y2 = zeros(numcurvesplot,nxpts+1);

%%%%%%%%%%%%%%%%%%%%% top row: set up cell densities %%%%%%%%%%%%%%%%%%%%%%
rm = sol(:,nxpts,4);
for i=1:numcurvesplot
    for j = 1:nxpts
        r(i,j) = rm(i*4-3) * x(j);
        y1(i,j) = sol(i*4-3,j,1);
        y2(i,j) = sol(i*4-3,j,2);        
    end
    r(i,nxpts+1) = r(i,nxpts);
end

% find max value of APC+IPA and round up to nearest thousand
ylim_sum = ceil(max(max(y1+y2))/1000)*1000;

%%%%%%%%%%%%%%%% PANEL: c1 - APCs
nexttile
hold on
for i=1:numcurvesplot
    plot(r(i,:),y1(i,:),'LineWidth',1.5,'Color',co(i,:))
end
hold off
xlabel('Radius (mm)','FontSize',fslabel,'Interpreter','latex')
ylabel('APCs (cells/mm$^2$)','FontSize',fslabel,'Interpreter','latex')
set(gca,'XLim',[0,rmax],'YLim',[0,ylim_sum],'FontSize',fsticks)
xticks(0:rmax)

%%%%%%%%%%%%%%%% PANEL: c2 - IPAs
nexttile
hold on
for i=1:numcurvesplot
    plot(r(i,:),y2(i,:),'LineWidth',1.5,'Color',co(i,:))
end
hold off
xlabel('Radius (mm)','FontSize',fslabel,'Interpreter','latex')
ylabel('IPAs (cells/mm$^2$)','FontSize',fslabel,'Interpreter','latex')
set(gca,'XLim',[0,rmax],'YLim',[0,ylim_sum],'FontSize',fsticks)
xticks(0:rmax)

%%%%%%%%%%%%%%%% PANEL: c1 + c2 - APCs + IPAs
nexttile
hold on
for i=1:numcurvesplot
    plot(r(i,:),y1(i,:)+y2(i,:),'LineWidth',1.5,'Color',co(i,:))
end
hold off
xlabel('Radius (mm)','FontSize',fslabel,'Interpreter','latex')
ylabel('APCs + IPAs (cells/mm$^2$)','FontSize',fslabel,'Interpreter','latex')
set(gca,'XLim',[0,rmax],'YLim',[0,ylim_sum],'FontSize',fsticks)
xticks(0:rmax)

%%%%%%%% bottom row: set up radial velocity in interior of domain %%%%%%%%%
for i=1:numcurvesplot
    usum = sol(i*4-3,:,1) + sol(i*4-3,:,2);
    for j = 1:nxpts
    if j == 1 
        DuDx = (nxpts - 1) * (usum(j+1) - usum(j));
    elseif  j == nxpts
        DuDx = (nxpts - 1) * (usum(j) - usum(j-1));
    else
        DuDx = (nxpts - 1)/2 * (usum(j+1) - usum(j-1));
    end
    y1(i,j) = -kTprime1 * usum(j)^(-3/2) / sol(i*4-3,j,4) / mu1 * DuDx;
    end
end

%%%%%%%%%%%%%%%% PANEL: radial velocity
nexttile
hold on
for i=1:numcurvesplot
    plot(r(i,:),y1(i,:),'LineWidth',1.5,'Color',co(i,:))
end
hold off
xlabel('Radius (mm)','FontSize',fslabel,'Interpreter','latex')
ylabel('APC Velocity (mm/hr)','FontSize',fslabel,'Interpreter','latex')
set(gca,'XLim',[0,rmax],'FontSize',fsticks)
xticks(0:rmax)

h = legend('E15','E16','E17','E18','E19','E20','E21','E22/P0');
set(h,'FontSize',fslegend,'Position',[0.355 0.12 0.1 0.34]);

%%%%%%%%%%%%%%%% PANEL: Radius vs. time / distance astrocytes spread
rad_APC = [0.17; 0.33; 0.5; 0.67; 1.67; 2.17; 2.67]; %not including E17
rad_days = [0; 1; 3; 4; 5; 6; 7];
nexttile(6)
hold on
plot(t/24,rm,'k','LineWidth',1.5)
scatter(rad_days,rad_APC,150,[0.5 0.5 0.5],'x','LineWidth',1.5)
hold off
xlabel('Time (days)','FontSize',fslabel,'Interpreter','latex')
ylabel('Cell Boundary (mm)','FontSize',fslabel,'Interpreter','latex')
set(gca,'FontSize',fsticks)

set(gcf,'Units','inches','Position',[1,1,12,6],'PaperPositionMode','auto');
