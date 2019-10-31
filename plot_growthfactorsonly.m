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

fslabel = 16;
fsticks = 14;
flegend = 10;

figure
subaxis(2,2,1,'MarginLeft',0.07,'MarginRight',0.05,'MarginTop',0.04,'MarginBottom',0.12)
hold on
j=0;
for i=tplot
    j=j+1;
    plot(r,q1_imp(i,:),'LineWidth',1.5,'Color',co(j,:))
end
hold off
xlabel('r (mm)','FontSize',fslabel)
ylabel('PDGFA (ng/mL)','FontSize',fslabel)
set(gca,'XLim',[0,rmax],'FontSize',fsticks)

subaxis(2,2,2,'MarginRight',0.01,'MarginLeft',0.07)
hold on
j=0;
for i=tplot
    j=j+1;
    plot(r,q2_imp(i,:),'LineWidth',1.5,'Color',co(j,:))
end
hold off
xlabel('r (mm)','FontSize',fslabel)
ylabel('LIF (ng/mL)','FontSize',fslabel)
set(gca,'XLim',[0,rmax],'FontSize',fsticks)

h = legend('0 days (E15)','1 day (E16)','2 days (E17)','3 days (E18)',...
    '4 days (E19)','5 days (E20)','6 days (E21)','7 days (E22/P0)');
set(h,'FontSize',flegend,'Position',[0.8188 0.651 0.1635 0.3]);

subaxis(2,2,3,'MarginLeft',0.07,'MarginRight',0.05,'MarginBottom',0.08,'MarginTop',0.08)
hold on
plot(radius_ret(tplot),t(tplot)/24,'k','LineWidth',1.5)
scatter(radius_ret(tplot),t(tplot)/24,70,co,'filled')
hold off
ylabel('t (days)','FontSize',fslabel)
xlabel('retina radius (mm)','FontSize',fslabel)
yticks(0:7)
set(gca,'YLim',[0,7],'XLim',[0,rmax],'FontSize',fsticks)

subaxis(2,2,4,'MarginRight',0.01,'MarginLeft',0.07)
hold on
plot(radius_endo(tplot),t(tplot)/24,'k','LineWidth',1.5)
scatter(radius_endo(tplot),t(tplot)/24,70,co,'filled')
hold off
ylabel('t (days)','FontSize',fslabel)
xlabel('endothelial cell radius (mm)','FontSize',fslabel)
yticks(0:7)
set(gca,'YLim',[0,7],'XLim',[0,rmax],'FontSize',fsticks)

set(gcf,'Units','inches','Position',[2,2,10,8],'PaperPositionMode','auto')