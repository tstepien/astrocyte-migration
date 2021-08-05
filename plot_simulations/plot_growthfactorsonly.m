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

fslabel = 22;
fsticks = 18;
flegend = 14;
fstitles = 24;

% set up radius
[~,~,radius_endo,radius_ret] = thick_rad(t,r);

% set up PDGFA and LIF
r = zeros(numcurvesplot,nxpts+1);
y1 = zeros(numcurvesplot,nxpts+1);
y2 = zeros(numcurvesplot,nxpts+1);
for i=1:numcurvesplot
    for j = 1:nxpts
        r(i,j) = rmax * x(j);
        y1(i,j) = PDGFA(i*4-3,j);
        y2(i,j) = LIF(i*4-3,j);        
    end
    r(i,nxpts+1) = r(i,nxpts);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(subpanels_on,'yes')==1
    figure
    subaxis(2,2,1,'MarginLeft',0.07,'MarginRight',0.05,'MarginTop',0.04,'MarginBottom',0.12)
else
    figure
end
hold on
for i=1:numcurvesplot
   plot(r(i,:),y1(i,:),'LineWidth',1.5,'Color',co(i,:)) 
end
hold off
xlabel('$r$ (mm)','FontSize',fslabel,'Interpreter','latex')
ylabel('PDGFA (ng/mL)','FontSize',fslabel,'Interpreter','latex')
box off

h = legend('E15','E16','E17','E18','E19','E20','E21','E22/P0');
set(h,'FontSize',flegend);

if strcmp(subpanels_on,'yes')==0
    set(gca,'XLim',[0,rmax],'FontSize',fsticks,'Position',[0.14 0.14 0.82 0.76])
    set(gcf,'Units','inches','Position',[2 2 7.75 5.75],'PaperPositionMode','auto')
else
    set(gca,'XLim',[0,rmax],'FontSize',fsticks)
end

if strcmp(titles_on,'yes')==1
    title('A                                                    ',...
    'FontSize',fstitles)
end

pause(0.1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(subpanels_on,'yes')==1
    subaxis(2,2,2,'MarginRight',0.01,'MarginLeft',0.07)
else
    figure
end
hold on
for i=1:numcurvesplot
   plot(r(i,:),y2(i,:),'LineWidth',1.5,'Color',co(i,:)) 
end
hold off
xlabel('$r$ (mm)','FontSize',fslabel,'Interpreter','latex')
ylabel('LIF (ng/mL)','FontSize',fslabel,'Interpreter','latex')
box off

h = legend('E15','E16','E17','E18','E19','E20','E21','E22/P0');
set(h,'FontSize',flegend);

if strcmp(subpanels_on,'yes')==0
    set(gca,'XLim',[0,rmax],'FontSize',fsticks,'Position',[0.14 0.14 0.82 0.76])
    set(gcf,'Units','inches','Position',[2 2 7.75 5.75],'PaperPositionMode','auto')
else
    set(gca,'XLim',[0,rmax],'FontSize',fsticks)
end

if strcmp(titles_on,'yes')==1
    title('B                                                    ',...
    'FontSize',fstitles)
end

pause(0.1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(subpanels_on,'yes')==1
    subaxis(2,2,3,'MarginLeft',0.07,'MarginRight',0.05,'MarginBottom',0.08,'MarginTop',0.08)
else
    figure
end
hold on
plot(radius_ret(tplot),t(tplot)/24,'k','LineWidth',1.5)
scatter(radius_ret(tplot),t(tplot)/24,70,co,'filled')
hold off
ylabel('$t$ (days)','FontSize',fslabel,'Interpreter','latex')
xlabel('Retinal radius (mm)','FontSize',fslabel,'Interpreter','latex')
yticks(0:7)
box off

if strcmp(subpanels_on,'yes')==0
    set(gca,'XLim',[0,rmax],'FontSize',fsticks,'Position',[0.14 0.14 0.82 0.76])
    set(gcf,'Units','inches','Position',[2 2 7.75 5.75],'PaperPositionMode','auto')
else
    set(gca,'XLim',[0,rmax],'YLim',[0,7],'FontSize',fsticks)
end

if strcmp(titles_on,'yes')==1
    title('C                                                    ',...
    'FontSize',fstitles)
end

pause(0.1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(subpanels_on,'yes')==1
    subaxis(2,2,4,'MarginRight',0.01,'MarginLeft',0.07)
else
    figure
end
hold on
plot(radius_endo(tplot),t(tplot)/24,'k','LineWidth',1.5)
scatter(radius_endo(tplot),t(tplot)/24,70,co,'filled')
hold off
ylabel('$t$ (days)','FontSize',fslabel,'Interpreter','latex')
xlabel('Endothelial cell radius (mm)','FontSize',fslabel,'Interpreter','latex')
yticks(0:7)
box off

if strcmp(subpanels_on,'yes')==0
    set(gca,'XLim',[0,rmax],'FontSize',fsticks,'Position',[0.14 0.14 0.82 0.76])
    set(gcf,'Units','inches','Position',[2 2 7.75 5.75],'PaperPositionMode','auto')
else
    set(gca,'XLim',[0,rmax],'YLim',[0,7],'FontSize',fsticks)
end

if strcmp(titles_on,'yes')==1
    title('D                                                    ',...
    'FontSize',fstitles)
end

pause(0.1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(subpanels_on,'yes')==1
    set(gcf,'Units','inches','Position',[2,2,13,10],'PaperPositionMode','auto')
end