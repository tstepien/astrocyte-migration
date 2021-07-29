%%% plot times: day 0, 1, 2, 3, 4, 5, 6, 7
numcurvesplot = 8;
rad_APC = [0.17; 0.33; 0.5; 0.67; 1.67; 2.17; 2.67]; %not including E17

%%% color order - to match Chan-Ling color scheme for APC/IPA
co = [171/255 5/255 32/255;
    128/255 133/255 156/255];

fslabel = 16;
fsticks = 14;
fslegend = 14;

figure(8);
tiledlayout(2,4,'TileSpacing','Compact','Padding','compact')
% extra points to make plots go to axis at outer edge
r = zeros(numcurvesplot,nxpts+1);
y1 = zeros(numcurvesplot,nxpts+1);
y2 = zeros(numcurvesplot,nxpts+1);

%%%%%%%%%%%%%%%%%%%%% top row: set up cell densities %%%%%%%%%%%%%%%%%%%%%%
rm = sol(:,nxpts,4);
for i=1:numcurvesplot
    for j = 1:nxpts
        r(i,j) = rm(i*4-3) * x(j);
        y1(i,j) = sol(i*4-3,j,1);
        y2(i,j) = sol(i*4-3,j,2);        
    end
    r(i,nxpts+1) = r(i,nxpts);
end

% find max value of APC+IPA and round up to nearest thousand
ylim_sum = ceil(max(max(y1+y2))/1000)*1000;

for i=1:numcurvesplot
    nexttile
    hold on
    plot(r(i,:),y1(i,:),'LineWidth',1.5,'Color',co(1,:))
    plot(r(i,:),y2(i,:),'LineWidth',1.5,'Color',co(2,:))
    
    if i<3
        line([rad_APC(i),rad_APC(i)],[0 ylim_sum],'LineStyle','--',...
            'Color',[0.5,0.5,0.5],'LineWidth',1.25)
    elseif i>3
        line([rad_APC(i-1),rad_APC(i-1)],[0 ylim_sum],'LineStyle','--',...
            'Color',[0.5,0.5,0.5],'LineWidth',1.25)
    end
    hold off
    
    if i>=5
        xlabel('Radius (mm)','FontSize',fslabel,'Interpreter','latex')
    end
    if i==1 || i==5
        ylabel('APCs + IPAs (cells/mm$^2$)','FontSize',fslabel,'Interpreter','latex')
    end
    set(gca,'XLim',[0,rmax],'YLim',[0,ylim_sum],'FontSize',fsticks)
    xticks(0:rmax)
end

h = legend('APCs','IPAs');
set(h,'FontSize',fsticks,'Position',[0.163 0.867 0.085 0.088]);

set(gcf,'Units','inches','Position',[1,1,12,6],'PaperPositionMode','auto');
