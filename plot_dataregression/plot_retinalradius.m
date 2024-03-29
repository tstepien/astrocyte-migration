clear variables;
clc;
close all;

% days since E15
timeindays = (0:7)';

timeindays_more = linspace(0,7,1000);

% radius data in microns
radiusdata = [1170
1500
1830
2170
2500
2830
3830
4000];

% regression lines
line_radius = 414.7*timeindays_more + 1029.2;

markersize = 100;

fstitle = 26;
fslabel = 26;
fsticks = 22;

figure
hold on
plot(timeindays_more,line_radius/1000,'LineWidth',2.5,'Color','k')

scatter(timeindays,radiusdata/1000,markersize,'filled',...
    'MarkerFaceColor','k');
hold off

box off
set(gca,'FontSize',fsticks,'Position',[0.175 0.08 0.79 0.82])

%xlabel('Days since E15','Interpreter','latex','FontSize',fslabel)
ylabel('Retinal radius (mm)','Interpreter','latex','FontSize',fslabel)

xlim([timeindays(1),timeindays(end)])
xticks(timeindays)
xticklabels({'E15','E16','E17','E18','E19','E20','E21','E22'})

title('D                                                ',...
    'FontSize',fstitle)

set(gcf,'Units','inches','Position',[2 2 7.75 5.75],'PaperPositionMode','auto')