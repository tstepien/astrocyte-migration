clear variables;
clc;
close all;

% days since E15
timeindays = (-2:6)';

timeindays_more = linspace(-2,6,1000);

% thickness data in microns
thickness_posterior = [75
87
101
101
110
161
161
163
188];

thickness_peripheral = [51
65
73
81
89
94
144
154
152];

% regression lines
line_posterior = 14.33*timeindays_more + 98.78;
line_peripheral = 13.77*timeindays_more + 72.8;

markersize = 100;
graycolor = [0.65 0.65 0.65];

figure
hold on
plot(timeindays_more,line_posterior,'LineWidth',2,'Color','k')
plot(timeindays_more,line_peripheral,'LineWidth',2,'Color',graycolor)

h1 = scatter(timeindays,thickness_posterior,markersize,'filled',...
    'MarkerFaceColor','k');
h2 = scatter(timeindays,thickness_peripheral,markersize,'^','filled',...
    'MarkerFaceColor',graycolor);
hold off

box on
set(gca,'FontSize',18,'Position',[0.14 0.08 0.82 0.82])

%xlabel('Days since E15','Interpreter','latex','FontSize',22)
ylabel('Retinal Thickness ($\mu$m)','Interpreter','latex','FontSize',22)

legend([h1 h2],{'Posterior','Peripheral'},'Location','northwest')

xticks(timeindays)
xticklabels({'E13','E14','E15','E16','E17','E18','E19','E20','E21'})

title('C                                                ',...
    'FontSize',26)

set(gcf,'Units','inches','Position',[2 2 7.75 5.75],'PaperPositionMode','auto')