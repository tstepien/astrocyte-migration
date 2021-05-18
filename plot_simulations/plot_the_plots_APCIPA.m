parameters_fixed

%%% APC, IPA, and retina radius (mm) for E15-E16, E18-E22/P0
rad_APC = [0.17; 0.33; 000; 0.5; 0.67; 1.67; 2.17; 2.67];

%%% set up the plots
%%% variables
T = length(t);
R = length(r);
dr = r(2)-r(1);
rmax = m.rmax;
j_init = s0/dr+1;

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

%%% color order - MATLAB
% co = [0.8500    0.3250    0.0980 %red
%     0    0.4470    0.7410 %blue
%     ];

%%% color order - to match Chan-Ling color scheme for APC/IPA
co = [171/255 5/255 32/255;
    128/255 133/255 156/255];


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
for i=1:8
    if i==1
        subaxis(2,4,1,'MarginTop',0.07,'MarginBottom',0.175,'MarginLeft',0.06,'MarginRight',0.01)
    elseif i==5
        subaxis(2,4,5,'MarginTop',0.14,'MarginBottom',0.1)
    else
        subaxis(2,4,i)
    end
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
    if i==8
        title([num2str(t(plotind(i))/24,3),' days (E',num2str(round(15+t(plotind(i))/24,1)),'/P0)']);
    else
        title([num2str(t(plotind(i))/24,3),' days (E',num2str(round(15+t(plotind(i))/24,1)),')']);
    end
    if i>=5
        xlabel('Radius (mm)','FontSize',fslabel,'Interpreter','latex')
    end
    if mvgbdy(end)<1.5
        set(gca,'XLim',[0,mvgbdy(end)+5*dr],'YLim',ylim_sum,'FontSize',fsticks)
    else
        set(gca,'XLim',[0,rmax-1],'YLim',ylim_sum,'FontSize',fsticks)
        xticks([0 1 2 3 4])
    end
end

h = legend('APCs','IPAs');
set(h,'FontSize',fsticks,'Position',[0.155 0.8102 0.0946 0.1098]);

text(-16,12200,'A','FontSize',26,'FontWeight','bold')

set(gcf,'Units','inches','Position',[2,2,12,6.5],'PaperPositionMode','auto')