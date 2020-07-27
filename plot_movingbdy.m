%%% max distance astrocytes spread
max_astrocytes = 2.67;

rad_APC = [0.17; 0.33; 0.5; 0.67; 1.67; 2.17; 2.67]; %not including E17
rad_days = [0; 1; 3; 4; 5; 6; 7];

%%% set up the plots

%%% variables
T = length(t);

%%% plot times: day 0, 1, 2, 3, 4, 5, 6, 7
numcurvesplot = 8;
plotind = zeros(1,numcurvesplot);
for i=0:numcurvesplot-1
    plotind(i+1) = find(abs((t/24-i))==min(abs(t/24-i)));
end

fslabel = 22;
fsticks = 18;

co = [0    0.4470    0.7410
    0.8500    0.3250    0.0980
    0.9290    0.6940    0.1250
    0.4940    0.1840    0.5560
    0.4660    0.6740    0.1880
    0.3010    0.7450    0.9330
    0.6350    0.0780    0.1840
    0         0.3725    0];

figure
hold on
    plot(t/24,mvgbdy,'k','LineWidth',1.5)
    scatter(t/24,mvgbdy,20,'k')
    scatter(rad_days,rad_APC,150,[0.5 0.5 0.5],'x','LineWidth',1.5)
hold off
xlabel('Time (days)','FontSize',fslabel,'Interpreter','latex')
ylabel('Moving Boundary (mm)','FontSize',fslabel,'Interpreter','latex')
ylim_mvgbdy = get(gca,'YLim');
if ylim_mvgbdy(2)<max_astrocytes
    set(gca,'FontSize',fsticks,'YLim',[0,max_astrocytes],'Position',[0.14 0.14 0.82 0.76])
else
    set(gca,'FontSize',fsticks,'XLim',[0,7],'Position',[0.14 0.14 0.82 0.76])
    xticks([0 1 2 3 4 5 6 7])
end

box on

title('B                                                ',...
    'FontSize',26)

set(gcf,'Units','inches','Position',[2 2 7.75 5.75],'PaperPositionMode','auto')