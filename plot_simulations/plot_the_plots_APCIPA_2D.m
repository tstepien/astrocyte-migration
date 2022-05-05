%%% plot times: day 0, 1, 2, 3, 4, 5, 6, 7
numcurvesplot = 8;
rad_APC = [0.17; 0.33; 0.5; 0.67; 1.67; 2.17; 2.67]; %not including E17

%%% color order - to match Chan-Ling color scheme for APC/IPA
co = [171/255 5/255 32/255;
    128/255 133/255 156/255];

load('chanlingBlueColormap.mat');
load('chanlingRedColormap.mat');

fslabel = 16;
fsticks = 14;
fslegend = 14;

% figure(9);
% tiledlayout(2,4,'TileSpacing','Compact','Padding','compact')
% extra points to make plots go to axis at outer edge
rplot = zeros(numcurvesplot,nxpts+1);
y1 = zeros(numcurvesplot,nxpts+1);
y2 = zeros(numcurvesplot,nxpts+1);

%%%%%%%%%%%%%%%%%%%%%%%%%% set up cell densities %%%%%%%%%%%%%%%%%%%%%%%%%%
rm = sol(:,nxpts,4);
for i=1:numcurvesplot
    for j = 1:nxpts
        rplot(i,j) = rm(i*4-3) * x(j);
        y1(i,j) = sol(i*4-3,j,1);
        y2(i,j) = sol(i*4-3,j,2);        
    end
    rplot(i,nxpts+1) = rplot(i,nxpts);
end

% find max value of APC+IPA and round up to nearest thousand
ylim_sum = ceil(max(max(y1+y2))/1000)*1000;

%%%%%%%%%%%%%%%%%%%%%%%% set up 2D cell densities %%%%%%%%%%%%%%%%%%%%%%%%%
gridpts = 100;
xgrid = linspace(-5,5,gridpts);
ygrid = linspace(-5,5,gridpts);

density1 = cell(numcurvesplot,1);
density2 = cell(numcurvesplot,1);
for k=1:numcurvesplot
    density1{k} = zeros(gridpts);
    density2{k} = zeros(gridpts);
    for i=1:gridpts
        for j=1:gridpts
            rr = sqrt(xgrid(i)^2+ygrid(j)^2);
            ind = find(rr<=rplot(k,:),1,'first');
            if isempty(ind)==0
                density1{k}(i,j) = y1(k,ind);
                density2{k}(i,j) = y2(k,ind);
            end
        end
    end
end

k = 7;

% first try
% figure
% tiledlayout(1,2)
% ax1 = nexttile;
% imagesc(density1{k})
% colormap(ax1,chanlingRedColormap)
% 
% ax2 = nexttile;
% imagesc(density2{k})
% colormap(ax2,chanlingBlueColormap)


%second try
figure
ax1 = axes;
im1 = imagesc(flipud(density1{k})); %imagesc(C,clims)
%alpha(im1,0.5);
axis xy
xlabel('x (mm)')
ylabel('y (mm)')

ax2 = axes;
im2 = imagesc(flipud(density2{k}));
alpha(im2,0.5);
axis xy

linkaxes([ax1,ax2])
ax2.Visible = 'off';
ax2.XTick = [];
ax2.YTick = [];
ax1.XLim = [0,gridpts];
ax1.XTick = linspace(0,gridpts,11);
ax1.XTickLabel = -5:5;
ax1.YLim = ax1.XLim;
ax1.YTick = ax1.XTick;
ax1.YTickLabel = ax1.XTickLabel;
colormap(ax1,chanlingRedColormap)
colormap(ax2,chanlingBlueColormap)

stop

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% rest of original APCIPA file
for i=1:numcurvesplot
    nexttile
    hold on
    plot(rplot(i,:),y1(i,:),'LineWidth',1.5,'Color',co(1,:))
    plot(rplot(i,:),y2(i,:),'LineWidth',1.5,'Color',co(2,:))
    
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
