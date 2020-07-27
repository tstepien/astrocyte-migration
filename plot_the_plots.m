parameters_fixed

%%% max distance astrocytes spread
max_astrocytes = 2.67;

rad_APC = [0.17; 0.33; 0.5; 0.67; 1.67; 2.17; 2.67]; %not including E17
rad_days = [0; 1; 3; 4; 5; 6; 7];

%%% set up the plots
%%% variables
R = length(r);
dr = r(2)-r(1);
rmax = m.rmax;
j_init = s0/dr+1;

%%% cell layer thickness and radius
[thickness_ret,~,~,~] = thick_rad(t,r);
    
%%% oxygen
PO2 = oxygen(r,thickness_ret,P0,Dalpha,M0);
% PO2frac = PO2./(Pm+PO2);

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
subaxis(3,3,3,'MarginLeft',0.07,'MarginRight',0.005,'MarginTop',0.03,'MarginBottom',0.15)
hold on
for i=1:numcurvesplot
    plot(rplot(plotind(i),:),c1plot(plotind(i),:)+c2plot(plotind(i),:),...
        'LineWidth',1.5,'Color',co(i,:))
end
% ylim_sum = get(gca,'YLim');
    ylim_sum = [0,5000];
line([max_astrocytes,max_astrocytes],ylim_sum,'LineStyle','--',...
    'Color',[0.5,0.5,0.5],'LineWidth',1.25)
hold off
xlabel('radius (mm)','FontSize',fslabel)
ylabel('APCs + IPAs (cells/mm^2)','FontSize',fslabel)
if mvgbdy(end)<1.5
    set(gca,'XLim',[0,mvgbdy(end)+5*dr],'FontSize',fsticks)
else
    set(gca,'XLim',[0,rmax],'FontSize',fsticks)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subaxis(3,3,1,'MarginLeft',0.055,'MarginRight',0.01)
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
subaxis(3,3,2,'MarginLeft',0.06,'MarginRight',0.01)
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
subaxis(3,3,4,'MarginTop',0.03,'MarginBottom',0.08)
hold on
plot(t/24,mvgbdy,'k','LineWidth',1.5)
scatter(t/24,mvgbdy,20,'k')
scatter(rad_days,rad_APC,150,[0.5 0.5 0.5],'x','LineWidth',1.5)
hold off
xlabel('t (days)','FontSize',fslabel)
ylabel('moving boundary (mm)','FontSize',fslabel)
ylim_mvgbdy = get(gca,'YLim');
if ylim_mvgbdy(2)<max_astrocytes
    set(gca,'FontSize',fsticks,'YLim',[0,max_astrocytes])
else
    set(gca,'FontSize',fsticks,'XLim',[0,7])
    xticks([0 1 2 3 4 5 6 7])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subaxis(3,3,5,'MarginTop',0.03)
if size(thickness_ret,1)==1 || size(thickness_ret,2)==1
    plot(t/24,thickness_ret,'-o','LineWidth',1.5)
    xlabel('t (days)','FontSize',fslabel)
else
    hold on
    for i=1:numcurvesplot
        plot(r,thickness_ret(plotind(i),:),'LineWidth',1.5,'Color',co(i,:))
    end
    hold off
    xlabel('radius (mm)','FontSize',fslabel)
end
ylabel('retinal thickness (mm)','FontSize',fslabel)
set(gca,'XLim',[0,rmax],'FontSize',fsticks)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subaxis(3,3,6,'MarginRight',0.005)
hy = hyaloid(r,P_hy,r_hy);
if size(thickness_ret,1)==1 || size(thickness_ret,2)==1
    plot(thickness_ret,hy+PO2,'-o','LineWidth',1.5)
    xlabel('total retinal thickness (mm)','FontSize',fslabel)
else
    hold on
    for i=1:numcurvesplot
        plot(r,hy+PO2(plotind(i),:),'LineWidth',1.5,'Color',co(i,:))
    end
    hold off
    xlabel('radius (mm)','FontSize',fslabel)
end
ylabel('Choroid+Hyaloid P_{O_2} (mmHg)','FontSize',fslabel)
set(gca,'XLim',[0,rmax],'FontSize',fsticks)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subaxis(3,3,7,'MarginBottom',0.06)
hold on
for i=1:numcurvesplot
    plot(r,q1(plotind(i),:),'LineWidth',1.5,'Color',co(i,:))
end
hold off
xlabel('radius (mm)','FontSize',fslabel)
ylabel('PDGFA (ng/mL)','FontSize',fslabel)
set(gca,'XLim',[0,rmax],'FontSize',fsticks)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subaxis(3,3,8,'MarginBottom',0.06)
hold on
for i=1:numcurvesplot
    plot(r,q2(plotind(i),:),'LineWidth',1.5,'Color',co(i,:))
end
hold off
xlabel('radius (mm)','FontSize',fslabel)
ylabel('LIF (ng/mL)','FontSize',fslabel)
set(gca,'XLim',[0,rmax],'FontSize',fsticks)

h = legend([num2str(t(plotind(1))/24),' days (E15)'],...
    [num2str(t(plotind(2))/24,3),' days (E',num2str(round(15+t(plotind(2))/24,1)),')'],...
    [num2str(t(plotind(3))/24,3),' days (E',num2str(round(15+t(plotind(3))/24,1)),')'],...
    [num2str(t(plotind(4))/24,3),' days (E',num2str(round(15+t(plotind(4))/24,1)),')'],...
    [num2str(t(plotind(5))/24,3),' days (E',num2str(round(15+t(plotind(5))/24,1)),')'],...
    [num2str(t(plotind(6))/24,3),' days (E',num2str(round(15+t(plotind(6))/24,1)),')'],...
    [num2str(t(plotind(7))/24,3),' days (E',num2str(round(15+t(plotind(7))/24,1)),')'],...
    [num2str(t(plotind(8))/24,3),' days (E',num2str(round(15+t(plotind(8))/24,1)),'/P0)']);
set(h,'FontSize',fsticks,'Position',[0.713 0.08 0.13 0.22]);

set(gcf,'Units','inches','Position',[2,2,16,10],'PaperPositionMode','auto')



%% plot velocities

figure
subaxis(2,2,1,'MarginLeft',0.05,'MarginRight',0.01,'MarginTop',0.03,'MarginBottom',0.05)
hold on
for i=1:numcurvesplot
    plot(rplot(plotind(i),:),vel_cirplot(plotind(i),:))
end
hold off
xlabel('r')
ylabel('circumferential velocity')
set(gca,'XLim',[0,mvgbdy(end)+5*dr])

subaxis(2,2,2,'MarginLeft',0.05,'MarginRight',0.01)
hold on
for i=1:numcurvesplot
    plot(rplot(plotind(i),:),vel_radplot(plotind(i),:))
end
hold off
xlabel('r')
ylabel('radial velocity')
set(gca,'XLim',[0,mvgbdy(end)+5*dr])

subaxis(2,2,3,'MarginTop',0.03,'MarginBottom',0.05)
hold on
for i=1:numcurvesplot
    plot(rplot(plotind(i),:),vel_cirplot(plotind(i),:)+vel_radplot(plotind(i),:))
end
hold off
xlabel('r')
ylabel('circumferential + radial velocity')
set(gca,'XLim',[0,mvgbdy(end)+5*dr])

subaxis(2,2,4,'MarginTop',0.03,'MarginBottom',0.05)
hold on
for i=1:numcurvesplot
    plot(rplot(plotind(i),:),vel_cirplot(plotind(i),:).*rplot(plotind(i),:))
end
hold off
xlabel('r')
ylabel('velocity')
set(gca,'XLim',[0,mvgbdy(end)+5*dr])


set(gcf,'Units','inches','Position',[2,2,12,8],'PaperPositionMode','auto')