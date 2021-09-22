% figure  % check boundary condition at outer edge satisfied
% plot(t,sol(:,nxpts,1)+sol(:,nxpts,2))
% figure
% surf(x,t,sol(:,:,1));
% figure
% surf(x,t,sol(:,:,2));
% figure
% surf(x,t,sol(:,:,3));
% figure
% surf(x,t,sol(:,:,4));

% color order
co = [0    0.4470    0.7410
    0.8500    0.3250    0.0980
    0.9290    0.6940    0.1250
    0.4940    0.1840    0.5560
    0.4660    0.6740    0.1880
    0.3010    0.7450    0.9330
    0.6350    0.0780    0.1840
    0         0.3725    0];

fslabel = 12;
numcurvesplot = floor(ntpts/4)+1;

figure
tiledlayout(3,3,'TileSpacing','Compact')
% extra points to make plots go to axis at outer edge
rplot = zeros(numcurvesplot,nxpts+2);
y1 = zeros(numcurvesplot,nxpts+2);
y2 = zeros(numcurvesplot,nxpts+2);

% set up retinal thickness and choroid oxygen
for i=1:numcurvesplot
    [~,~,~,radius_ret] = thick_rad(t(i*4-3),0);
    for j = 1:nxpts
        rplot(i,j) = radius_ret * x(j);
        [thickness_ret,~,~,~] = thick_rad(t(i*4-3),rplot(i,j));
        y1(i,j) = thickness_ret;
        y2(i,j) = interp1(Lvec,Pvec,thickness_ret);        
    end
    rplot(i,nxpts+1) = rplot(i,nxpts);   %to make plots go to axis at outer edge 
    rplot(i,nxpts+2) = rmax;  
end

% retinal thickness
ax1 = nexttile;
hold on
for i=1:numcurvesplot
    plot(rplot(i,:),y1(i,:),'LineWidth',1.5,'Color',co(i,:)) 
end
xlabel('Radius (mm)','FontSize',fslabel)
ylabel('Retinal thickness (mm)','FontSize',fslabel)
hold off

% choroid oxygen
ax2 = nexttile;
hold on
for i=1:numcurvesplot
    plot(rplot(i,:),y2(i,:),'LineWidth',1.5,'Color',co(i,:)) 
end
xlabel('Radius (mm)','FontSize',fslabel)
ylabel('Choroid P_{O2} (mmHg)','FontSize',fslabel)
hold off

% set up PDGFA and LIF
for i=1:numcurvesplot
    for j = 1:nxpts
        rplot(i,j) = rmax * x(j);
        y1(i,j) = PDGFA(i*4-3,j);
        y2(i,j) = LIF(i*4-3,j);        
    end
    rplot(i,nxpts+1) = rplot(i,nxpts);
end

% PDGFA
ax3 = nexttile;
hold on
for i=1:numcurvesplot
   plot(rplot(i,:),y1(i,:),'LineWidth',1.5,'Color',co(i,:)) 
end
xlabel('Radius (mm)','FontSize',fslabel)
ylabel('PDGFA (ng/mL)','FontSize',fslabel)
hold off

% LIF
ax4 = nexttile;
hold on
for i=1:numcurvesplot
   plot(rplot(i,:),y2(i,:),'LineWidth',1.5,'Color',co(i,:)) 
end
xlabel('Radius (mm)','FontSize',fslabel)
ylabel('LIF (ng/mL)','FontSize',fslabel)
hold off

% set up cell densities
rm = sol(:,nxpts,4);
for i=1:numcurvesplot
    for j = 1:nxpts
        rplot(i,j) = rm(i*4-3) * x(j);
        y1(i,j) = sol(i*4-3,j,1);
        y2(i,j) = sol(i*4-3,j,2);        
    end
    rplot(i,nxpts+1) = rplot(i,nxpts);
end

% c1 - APCs
ax5 = nexttile;
hold on
for i=1:numcurvesplot
    plot(rplot(i,:),y1(i,:),'LineWidth',1.5,'Color',co(i,:))
end
xlabel('Radius (mm)','FontSize',fslabel)
ylabel('APCs (cells/mm^2)','FontSize',fslabel)
hold off

% c2 - IPAs
ax6 = nexttile;
hold on
for i=1:numcurvesplot
    plot(rplot(i,:),y2(i,:),'LineWidth',1.5,'Color',co(i,:))
end
xlabel('Radius (mm)','FontSize',fslabel)
ylabel('IPAs (cells/mm^2)','FontSize',fslabel)
hold off

% c1 + c2 - APCs + IPAs
ax7 = nexttile;
hold on
for i=1:numcurvesplot
    plot(rplot(i,:),y1(i,:)+y2(i,:),'LineWidth',1.5,'Color',co(i,:))
end
xlabel('Radius (mm)','FontSize',fslabel)
ylabel('APCs + IPAs (cells/mm^2)','FontSize',fslabel)
hold off

%set up radial velocity in interior of domain
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

% radial velocity
ax8 = nexttile;
hold on
for i=1:numcurvesplot
    plot(rplot(i,:),y1(i,:),'LineWidth',1.5,'Color',co(i,:))
end
xlabel('Radius (mm)','FontSize',fslabel)
ylabel('APC velocity (mm/hr)','FontSize',fslabel)
hold off

axis([ax1,ax2,ax3,ax4,ax5,ax6,ax7,ax8],[0,5,-inf,inf]);

% Radius vs. time
% distance astrocytes spread
rad_APC = [0.17; 0.33; 0.5; 0.67; 1.67; 2.17; 2.67]; %not including E17
rad_days = [0; 1; 3; 4; 5; 6; 7];
nexttile
hold on
plot(t/24,rm,'k','LineWidth',1.5)
scatter(rad_days,rad_APC,150,[0.5 0.5 0.5],'x','LineWidth',1.5)
hold off
xlabel('t (days)','FontSize',fslabel)
ylabel('Cell boundary (mm)','FontSize',fslabel)

set(gcf,'Units','inches','Position',[1,1,10,7],'PaperPositionMode','auto');
