%%% max distance astrocytes spread
max_astrocytes = 2.67;

%%% set up the plots
[PO2,thickness] = oxygen(t,r);

T = length(t);
R = length(r);

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
subaxis(2,4,1,'MarginLeft',0.055,'MarginRight',0.01,'MarginTop',0.03,'MarginBottom',0.15)
ylim_sum = [0,5000];
hold on
for i=1:numcurvesplot
    plot(rplot(plotind(i),:),c1plot(plotind(i),:),'LineWidth',1.5,...
        'Color',co(i,:))
end
line([max_astrocytes,max_astrocytes],ylim_sum,'LineStyle','--',...
    'Color',[0.5,0.5,0.5],'LineWidth',1.25)
hold off
xlabel('radius (mm)','FontSize',fslabel)
ylabel('APCs (cells/mm^2)','FontSize',fslabel)
if mvgbdy(end)<1.5
    set(gca,'XLim',[0,mvgbdy(end)+5*dr],'FontSize',fsticks,'YLim',ylim_sum)
else
    set(gca,'XLim',[0,rmax],'FontSize',fsticks,'YLim',ylim_sum)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subaxis(2,4,2,'MarginLeft',0.06)
hold on
for i=1:numcurvesplot
    plot(rplot(plotind(i),:),c2plot(plotind(i),:),'LineWidth',1.5,...
        'Color',co(i,:))
end
line([max_astrocytes,max_astrocytes],ylim_sum,'LineStyle','--',...
    'Color',[0.5,0.5,0.5],'LineWidth',1.25)
hold off
xlabel('radius (mm)','FontSize',fslabel)
ylabel('IPAs (cells/mm^2)','FontSize',fslabel)
if mvgbdy(end)<1.5
    set(gca,'XLim',[0,mvgbdy(end)+5*dr],'FontSize',fsticks,'YLim',ylim_sum)
else
    set(gca,'XLim',[0,rmax],'FontSize',fsticks,'YLim',ylim_sum)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subaxis(2,4,3)
plot(t/24,mvgbdy,'-o')
xlabel('t (days)','FontSize',fslabel)
ylabel('Moving Boundary (mm)','FontSize',fslabel)
ylim_mvgbdy = get(gca,'YLim');
if ylim_mvgbdy(2)<max_astrocytes
    set(gca,'FontSize',fsticks,'YLim',[0,max_astrocytes])
else
    set(gca,'FontSize',fsticks,'XLim',[0,7])
    hold on
    line([0,7],[max_astrocytes,max_astrocytes],'LineStyle','--',...
        'Color',[0.5,0.5,0.5],'LineWidth',1.25)
    hold off
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subaxis(2,4,5,'MarginLeft',0.055,'MarginTop',0.08,'MarginBottom',0.1)
if size(thickness,1)==1 || size(thickness,2)==1
    plot(t/24,thickness,'-o','LineWidth',1.5)
    xlabel('t (days)','FontSize',fslabel)
else
    hold on
    for i=1:numcurvesplot
        plot(r,thickness(plotind(i),:),'LineWidth',1.5,'Color',co(i,:))
    end
    hold off
    xlabel('radius (mm)','FontSize',fslabel)
end
ylabel('Retinal Thickness (mm)','FontSize',fslabel)
set(gca,'XLim',[0,rmax],'FontSize',fsticks)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subaxis(2,4,6)
if size(thickness,1)==1 || size(thickness,2)==1
    plot(thickness,PO2,'-o','LineWidth',1.5)
    xlabel('total retinal thickness (mm)','FontSize',fslabel)
else
    hold on
    for i=1:numcurvesplot
        plot(r,PO2(plotind(i),:),'LineWidth',1.5,'Color',co(i,:))
    end
    hold off
    xlabel('radius (mm)','FontSize',fslabel)
end
ylabel('Partial Pressure of O_2 (mmHg)','FontSize',fslabel)
set(gca,'XLim',[0,rmax],'FontSize',fsticks)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subaxis(2,4,7)
hold on
for i=1:numcurvesplot
    plot(r,q1(plotind(i),:),'LineWidth',1.5,'Color',co(i,:))
end
hold off
xlabel('radius (mm)','FontSize',fslabel)
ylabel('Growth Factor: PDGFA (ng/mL)','FontSize',fslabel)
set(gca,'XLim',[0,rmax],'FontSize',fsticks)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subaxis(2,4,8)
hold on
for i=1:numcurvesplot
    plot(r,q2(plotind(i),:),'LineWidth',1.5,'Color',co(i,:))
end
hold off
xlabel('radius (mm)','FontSize',fslabel)
ylabel('Growth Factor: LIF (ng/mL)','FontSize',fslabel)
set(gca,'XLim',[0,rmax],'FontSize',fsticks)

h = legend('0 days','1 day','2 days','3 days','4 days','5 days','6 days','7 days');
set(h,'FontSize',fsticks,'Position',[0.85 0.57827 0.0833 0.372]);

set(gcf,'Units','inches','Position',[2,2,16,6],'PaperPositionMode','auto')
