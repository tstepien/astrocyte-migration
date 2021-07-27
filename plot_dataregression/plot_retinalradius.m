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

figure
hold on
plot(timeindays_more,line_radius,'LineWidth',2,'Color','k')

scatter(timeindays,radiusdata,markersize,'filled',...
    'MarkerFaceColor','k');
hold off

box on
set(gca,'FontSize',18,'Position',[0.14 0.08 0.82 0.82])

%xlabel('Days since E15','Interpreter','latex','FontSize',22)
ylabel('Retinal Radius ($\mu$m)','Interpreter','latex','FontSize',22)

xlim([timeindays(1),timeindays(end)])
xticks(timeindays)
xticklabels({'E15','E16','E17','E18','E19','E20','E21','E22'})

title('D                                                ',...
    'FontSize',26)

set(gcf,'Units','inches','Position',[2 2 7.75 5.75],'PaperPositionMode','auto')