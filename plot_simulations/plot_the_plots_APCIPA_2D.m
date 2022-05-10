%%% NOTE: must save as a .pdf file to get the correct coloring and size

%%% plot times: day 0, 1, 2, 3, 4, 5, 6, 7
numcurvesplot = 8;
rad_APC = [0.17; 0.33; 0.5; 0.67; 1.67; 2.17; 2.67]; %not including E17
legend_days = {'E15','E16','E17','E18','E19','E20','E21','E22/P0'};

%%% color order - to match Chan-Ling color scheme for APC/IPA
co = [171/255 5/255 32/255;
    128/255 133/255 156/255];

load('chanlingBlueColormap.mat');
load('chanlingRedColormap.mat');

fslabel = 16;
fsticks = 14;
fslegend = 14;
fscolorbar = 12;

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

% find max value of APC/IPA and round up to nearest 1500
ylim = max(max(y1(:)),max(y2(:)));

%%%%%%%%%%%%%%%%%%%%%%%% set up 2D cell densities %%%%%%%%%%%%%%%%%%%%%%%%%
gridpts = 2000;
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% plot panels %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
for i=1:numcurvesplot
    ax1 = axes;
    subplot(2,4,i,ax1)
    im1 = imagesc(ax1,flipud(density1{i}),[0,ylim]); % APC->red
    axis xy
    
    if i==4
        h1 = colorbar;
        set(h1,'FontSize',fscolorbar,...
            'Position',[0.921337093701508 0.583333333333333 0.011 0.342013888888889]);
        text(1.3*gridpts,1.15*gridpts,'APCs',...
            'HorizontalAlignment','center','FontSize',fscolorbar,...
            'Interpreter','latex')
        text(1.3*gridpts,1.05*gridpts,'(cells/mm$^2$)',...
            'HorizontalAlignment','center','FontSize',fscolorbar,...
            'Interpreter','latex')
    end

    ax2 = axes;
    subplot(2,4,i,ax2)
    im2 = imagesc(ax2,flipud(density2{i}),[0,ylim]); % IPA->blue
    alpha(im2,0.5);
    axis xy
    
    if i==8
        h2 = colorbar;
        set(h2,'FontSize',fscolorbar,...
            'Position',[0.92148899866361 0.109375 0.011 0.342013888888889]);
        text(1.3*gridpts,1.15*gridpts,'IPAs',...
            'HorizontalAlignment','center','FontSize',fscolorbar,...
            'Interpreter','latex')
        text(1.3*gridpts,1.05*gridpts,'(cells/mm$^2$)',...
            'HorizontalAlignment','center','FontSize',fscolorbar,...
            'Interpreter','latex')
    end

    % scale bar
    if i==1
        hold on
        plot([0.12*gridpts+0; 0.12*gridpts+0.1*gridpts], ...
            [0.15*gridpts; 0.15*gridpts],'-k','LineWidth', 5)
        hold off
        text(0.03*gridpts,0.07*gridpts, '1 mm', 'HorizontalAlignment','left','FontSize',fsticks)
    end

    text(0.5*gridpts,1.05*gridpts,legend_days{i},'HorizontalAlignment','center',...
        'FontSize',fslabel);

    linkaxes([ax1,ax2])
    colormap(ax1,chanlingRedColormap)
    colormap(ax2,chanlingBlueColormap)
    ax2.Visible = 'off';
    ax2.XTick = [];
    ax2.YTick = [];
    ax1.XLim = [0,gridpts];
    ax1.XTick = linspace(0,gridpts,11);
    ax1.XTickLabel = [];
    ax1.YLim = ax1.XLim;
    ax1.YTick = ax1.XTick;
    ax1.YTickLabel = ax1.XTickLabel;
end

set(gcf,'Units','inches','Position',[1,1,12,6],'PaperPositionMode','auto',...
    'PaperOrientation','landscape','PaperPosition',[-1.45 -0.5 12 6],...
    'PaperSize',[10.5,5.5]);