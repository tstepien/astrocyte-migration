clear variables;
clc;
close all;

%%% retinal ganglion cells

% days since E15
timeindays = (1:6)';

timeindays_more = linspace(1,6,1000);

% thickness data in microns
thickness_posterior = [0
27
46
31
32
31];

thickness_peripheral = [0
8
26
30
37
26];

% regression lines
line_posterior = -3.79*timeindays_more.^2 + 31.02*timeindays_more - 23.16;
line_peripheral = -2.49*timeindays_more.^2 + 23.81*timeindays_more - 24.12;

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
ylabel('Ganglion cell layer thickness (mm)','Interpreter','latex','FontSize',fslabel)

legend([h1 h2],{'Posterior','Peripheral'},'Location','northeast')

xticks(timeindays)
xticklabels({'E16','E17','E18','E19','E20','E21'})
ylim([0,0.05])

title('B                                                ',...
    'FontSize',fstitle)

set(gcf,'Units','inches','Position',[2 2 7.75 5.75],'PaperPositionMode','auto')