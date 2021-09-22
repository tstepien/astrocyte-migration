titles_on = 'yes';
subpanels_on = 'no';

%%% plot times: day 0, 1, 2, 3, 4, 5, 6, 7
numcurvesplot = 8;
tplot = zeros(1,numcurvesplot);
for i=0:numcurvesplot-1
    tplot(i+1) = find(abs((t/24-i))==min(abs(t/24-i)));
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

% set up radius
radius_endo = zeros(length(t),1);
radius_ret = zeros(length(t),1);
for i=1:length(t)
    [~,~,radius_endo(i),radius_ret(i)] = thick_rad(t(i),rplot);
end

% set up PDGFA and LIF
rplot = zeros(numcurvesplot,nxpts+1);
y1 = zeros(numcurvesplot,nxpts+1);
y2 = zeros(numcurvesplot,nxpts+1);
for i=1:numcurvesplot
    for j = 1:nxpts
        rplot(i,j) = rmax * x(j);
        y1(i,j) = PDGFA(i*4-3,j);
        y2(i,j) = LIF(i*4-3,j);        
    end
    rplot(i,nxpts+1) = rplot(i,nxpts);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(subpanels_on,'yes')==1
    figure
    tiledlayout(2,2,'TileSpacing','compact','Padding','none')
    nexttile
else
    figure
end
hold on
for i=1:numcurvesplot
   plot(rplot(i,:),y1(i,:),'LineWidth',2.5,'Color',co(i,:)) 
end
hold off
xlabel('Radius (mm)','FontSize',fslabel,'Interpreter','latex')
ylabel('PDGFA (ng/mL)','FontSize',fslabel,'Interpreter','latex')
box off

h = legend('E15','E16','E17','E18','E19','E20','E21','E22/P0');
set(h,'FontSize',fslegend,'Position',[0.813 0.447 0.181 0.434]);

if strcmp(subpanels_on,'yes')==0
    set(gca,'XLim',[0,rmax],'FontSize',fsticks,'Position',[0.175 0.14 0.79 0.76])
    set(gcf,'Units','inches','Position',[2 2 7.75 5.75],'PaperPositionMode','auto')
else
    set(gca,'XLim',[0,rmax],'FontSize',fsticks)
end

if strcmp(titles_on,'yes')==1
    title('A                                             ',...
    'FontSize',fstitle)
end

pause(0.1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(subpanels_on,'yes')==1
    nexttile
else
    figure
end
hold on
for i=1:numcurvesplot
   plot(rplot(i,:),y2(i,:),'LineWidth',2.5,'Color',co(i,:)) 
end
hold off
xlabel('Radius (mm)','FontSize',fslabel,'Interpreter','latex')
ylabel('LIF (ng/mL)','FontSize',fslabel,'Interpreter','latex')
box off

h = legend('E15','E16','E17','E18','E19','E20','E21','E22/P0');
set(h,'FontSize',fslegend);

if strcmp(subpanels_on,'yes')==0
    set(gca,'XLim',[0,rmax],'FontSize',fsticks,'Position',[0.175 0.14 0.79 0.76])
    set(gcf,'Units','inches','Position',[2 2 7.75 5.75],'PaperPositionMode','auto')
else
    set(gca,'XLim',[0,rmax],'FontSize',fsticks)
end

if strcmp(titles_on,'yes')==1
    title('B                                             ',...
    'FontSize',fstitle)
end

pause(0.1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(subpanels_on,'yes')==1
    nexttile
else
    figure
end
hold on
plot(radius_ret(tplot),t(tplot)/24,'k','LineWidth',2.5)
scatter(radius_ret(tplot),t(tplot)/24,70,co,'filled')
hold off
ylabel('Time (days)','FontSize',fslabel,'Interpreter','latex')
xlabel('Retinal radius (mm)','FontSize',fslabel,'Interpreter','latex')
yticks(0:7)
box off

if strcmp(subpanels_on,'yes')==0
    set(gca,'XLim',[0,rmax],'YLim',[0,7],'FontSize',fsticks,'Position',[0.175 0.14 0.79 0.76])
    set(gcf,'Units','inches','Position',[2 2 7.75 5.75],'PaperPositionMode','auto')
else
    set(gca,'XLim',[0,rmax],'YLim',[0,7],'FontSize',fsticks)
end

if strcmp(titles_on,'yes')==1
    title('C                                             ',...
    'FontSize',fstitle)
end

pause(0.1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(subpanels_on,'yes')==1
    nexttile
else
    figure
end
hold on
plot(radius_endo(tplot),t(tplot)/24,'k','LineWidth',2.5)
scatter(radius_endo(tplot),t(tplot)/24,70,co,'filled')
hold off
ylabel('Time (days)','FontSize',fslabel,'Interpreter','latex')
xlabel('Endothelial cell radius (mm)','FontSize',fslabel,'Interpreter','latex')
yticks(0:7)
box off

if strcmp(subpanels_on,'yes')==0
    set(gca,'XLim',[0,rmax],'YLim',[0,7],'FontSize',fsticks,'Position',[0.175 0.14 0.79 0.76])
    set(gcf,'Units','inches','Position',[2 2 7.75 5.75],'PaperPositionMode','auto')
else
    set(gca,'XLim',[0,rmax],'YLim',[0,7],'FontSize',fsticks)
end

if strcmp(titles_on,'yes')==1
    title('D                                             ',...
    'FontSize',fstitle)
end

pause(0.1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(subpanels_on,'yes')==1
    set(gcf,'Units','inches','Position',[2,2,13,10],'PaperPositionMode','auto')
end