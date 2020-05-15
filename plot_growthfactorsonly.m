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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(subpanels_on,'yes')==1
    figure
    subaxis(2,2,1,'MarginLeft',0.07,'MarginRight',0.05,'MarginTop',0.04,'MarginBottom',0.12)
else
    figure
end
hold on
j=0;
for i=tplot
    j=j+1;
    plot(r,q1_imp(i,:),'LineWidth',1.5,'Color',co(j,:))
end
hold off
xlabel('Radius (mm)','FontSize',fslabel,'Interpreter','latex')
ylabel('PDGFA (ng/mL)','FontSize',fslabel,'Interpreter','latex')
box on

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
stop

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(subpanels_on,'yes')==1
    subaxis(2,2,2,'MarginRight',0.01,'MarginLeft',0.07)
else
    figure
end
hold on
j=0;
for i=tplot
    j=j+1;
    plot(r,q2_imp(i,:),'LineWidth',1.5,'Color',co(j,:))
end
hold off
xlabel('Radius (mm)','FontSize',fslabel,'Interpreter','latex')
ylabel('LIF (ng/mL)','FontSize',fslabel,'Interpreter','latex')
box on

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
ylabel('Time (days)','FontSize',fslabel,'Interpreter','latex')
xlabel('Retina Radius (mm)','FontSize',fslabel,'Interpreter','latex')
yticks(0:7)
box on

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
ylabel('Time (days)','FontSize',fslabel,'Interpreter','latex')
xlabel('Endothelial Cell Radius (mm)','FontSize',fslabel,'Interpreter','latex')
yticks(0:7)
box on

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