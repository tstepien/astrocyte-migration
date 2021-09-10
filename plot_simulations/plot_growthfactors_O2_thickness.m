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
tiledlayout(2,2,'TileSpacing','Compact')
% extra points to make plots go to axis at outer edge
r = zeros(numcurvesplot,nxpts+2);
y1 = zeros(numcurvesplot,nxpts+2);
y2 = zeros(numcurvesplot,nxpts+2);

% set up retinal thickness and choroid oxygen
for i=1:numcurvesplot
    [~,~,~,radius_ret] = thick_rad(t(i*4-3),0);
    for j = 1:nxpts
        r(i,j) = radius_ret * x(j);
        [thickness_ret,~,~,~] = thick_rad(t(i*4-3),r(i,j));
        y1(i,j) = thickness_ret;
        y2(i,j) = interp1(Lvec,Pvec,thickness_ret);        
    end
    r(i,nxpts+1) = r(i,nxpts);   %to make plots go to axis at outer edge 
    r(i,nxpts+2) = rmax;  
end

% retinal thickness
ax1 = nexttile;
hold on
for i=1:numcurvesplot
    plot(r(i,:),y1(i,:),'LineWidth',1.5,'Color',co(i,:)) 
end
xlabel('Radius (mm)','FontSize',fslabel)
ylabel('Retinal thickness (mm)','FontSize',fslabel)
hold off

% choroid oxygen
ax2 = nexttile;
hold on
for i=1:numcurvesplot
    plot(r(i,:),y2(i,:),'LineWidth',1.5,'Color',co(i,:)) 
end
xlabel('Radius (mm)','FontSize',fslabel)
ylabel('Choroid P_{O2} (mmHg)','FontSize',fslabel)
hold off

% set up PDGFA and LIF
for i=1:numcurvesplot
    for j = 1:nxpts
        r(i,j) = rmax * x(j);
        y1(i,j) = PDGFA(i*4-3,j);
        y2(i,j) = LIF(i*4-3,j);        
    end
    r(i,nxpts+1) = r(i,nxpts);
end

% PDGFA
ax3 = nexttile;
hold on
for i=1:numcurvesplot
   plot(r(i,:),y1(i,:),'LineWidth',1.5,'Color',co(i,:)) 
end
xlabel('Radius (mm)','FontSize',fslabel)
ylabel('PDGFA (ng/mL)','FontSize',fslabel)
hold off

% LIF
ax4 = nexttile;
hold on
for i=1:numcurvesplot
   plot(r(i,:),y2(i,:),'LineWidth',1.5,'Color',co(i,:)) 
end
xlabel('Radius (mm)','FontSize',fslabel)
ylabel('LIF (ng/mL)','FontSize',fslabel)
hold off

axis([ax1,ax2,ax3,ax4],[0,5,-inf,inf]);

set(gcf,'Units','inches','Position',[1,1,8,6],'PaperPositionMode','auto');