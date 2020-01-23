parameters_fixed

%%% APC, IPA, and retina radius (mm) for E15-E16, E18-E22/P0
rad_APC = [0.33; 0.33; 000; 0.5; 0.67; 1.67; 2.17; 2.67];

%%% set up the plots
%%% variables
T = length(t);

%%% only plot a subset of the times
% numcurvesplot = 7;
% if T<=numcurvesplot
%     plotind = 1:T;
%     numcurvesplot = length(plotind);
% else
%     plotind = floor(linspace(1,T,numcurvesplot));
% end

%%% plot times: day 0, 1, 2, 3, 4, 5, 6, 7
numcurvesplot = 8;
plotind = zeros(1,numcurvesplot);
for i=0:numcurvesplot-1
    plotind(i+1) = find(abs((t/24-i))==min(abs(t/24-i)));
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

%%% make moving boundary sharp on the plot
c1plot = zeros(T,R+1);
c2plot = zeros(T,R+1);
vel_cirplot = zeros(T,R+1);
vel_radplot = zeros(T,R+1);
rplot = zeros(T,R+1);
for i = 1:T
    c1plot(i,:) = [c1(i,1:j_init+(i-1)) , zeros(1,R+1-(j_init+(i-1)))];
    c2plot(i,:) = [c2(i,1:j_init+(i-1)) , zeros(1,R+1-(j_init+(i-1)))];
    vel_cirplot(i,:) = [vel_cir(i,1:j_init+(i-1)) , zeros(1,R+1-(j_init+(i-1)))];
    vel_radplot(i,:) = [vel_rad(i,1:j_init+(i-1)) , zeros(1,R+1-(j_init+(i-1)))];
    
    rplot(i,:) = [r(1:j_init+(i-1)) , r(j_init+(i-1)) , r(j_init+i:end)];
end

fslabel = 16;
fsticks = 14;

%% plot cell concentrations

figure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:7
subplot(2,4,i)
hold on
    plot(rplot(plotind(i),:),c1plot(plotind(i),:),...
        'LineWidth',1.5,'Color',co(1,:))
    plot(rplot(plotind(i),:),c2plot(plotind(i),:),...
        'LineWidth',1.5,'Color',co(2,:))
    ylim_sum = [0,5000];
    if i~=3
        line([rad_APC(i),rad_APC(i)],ylim_sum,'LineStyle','--',...
            'Color',[0.5,0.5,0.5],'LineWidth',1.25)
    end
box on
hold off
title([num2str(t(plotind(i))/24,3),' days (E',num2str(round(15+t(plotind(i))/24,1)),')']);
if i>=5
    xlabel('radius (mm)','FontSize',fslabel)
end
if mvgbdy(end)<1.5
    set(gca,'XLim',[0,mvgbdy(end)+5*dr],'YLim',ylim_sum,'FontSize',fsticks)
else
    set(gca,'XLim',[0,rmax],'YLim',ylim_sum,'FontSize',fsticks)
end
end

h = legend('APCs','IPAs');
set(h,'FontSize',fsticks,'Position',[0.713 0.08 0.13 0.22]);

set(gcf,'Units','inches','Position',[2,2,12,6],'PaperPositionMode','auto')