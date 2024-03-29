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

fstitle = 26;
fslabel = 26;
fsticks = 22;
fslegend = 18;

figure
hold on
plot(timeindays_more,line_posterior/1000,'LineWidth',2.5,'Color','k')
plot(timeindays_more,line_peripheral/1000,'LineWidth',2.5,'Color',graycolor)

h1 = scatter(timeindays,thickness_posterior/1000,markersize,'filled',...
    'MarkerFaceColor','k');
h2 = scatter(timeindays,thickness_peripheral/1000,markersize,'^','filled',...
    'MarkerFaceColor',graycolor);
hold off

box off
set(gca,'FontSize',fsticks,'Position',[0.175 0.08 0.79 0.82])

%xlabel('Days since E15','Interpreter','latex','FontSize',fslabel)
ylabel('Retinal thickness (mm)','Interpreter','latex','FontSize',fslabel)

legend([h1 h2],{'Posterior','Peripheral'},'Location','northwest')

xticks(timeindays)
xticklabels({'E13','E14','E15','E16','E17','E18','E19','E20','E21'})

title('C                                                ',...
    'FontSize',fstitle)

set(gcf,'Units','inches','Position',[2 2 7.75 5.75],'PaperPositionMode','auto')