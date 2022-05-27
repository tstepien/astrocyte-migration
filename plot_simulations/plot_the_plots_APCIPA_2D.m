%%% NOTE: must save as a .pdf file to get the correct coloring and size

%%% need to plot using subplot (not tiledlayout) to get both colormaps on
%%% the same axes

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
rmax_plot = 3;
xgrid = linspace(-rmax_plot,rmax_plot,gridpts);
ygrid = linspace(-rmax_plot,rmax_plot,gridpts);

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

%%% the subplotposition matrix was populated by running the following code
%%% for gridpts = 100
%     outerpos = ax1.OuterPosition;
%     ti = ax1.TightInset;
%     left = outerpos(1) - lefttranslate(i)*ti(1);
%     bottom = outerpos(2) - bottomtranslate(i)*ti(2);
%     ax_width = outerpos(3) + 0.2*ti(3);
%     ax_height = outerpos(4) + 2.5*ti(4);
%     ax1.Position = [left bottom ax_width ax_height];
%     subplotposition(i,:) = [left bottom ax_width ax_height];
subplotposition = ...
    [0.00607143192418973	0.515376893594668	0.219515601048740	0.442070805247994
    0.236974166815948	0.515376893594668	0.219515601048740	0.442070805247994
    0.467876901707706	0.515376893594668	0.219526315334454	0.442070805247994
    0.698779636599463	0.515376893594668	0.219515601048740	0.442070805247994
    0.00607143192418973	0.0150181939821582	0.219515601048740	0.439319639609470
    0.236974166815948	0.0150181939821582	0.219515601048740	0.439319639609470
    0.467876901707706	0.0150181939821582	0.219526315334454	0.439319639609470
    0.698779636599463	0.0150181939821582	0.219515601048740	0.439319639609470];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% plot panels %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
lefttranslate = repmat([2,1.4,0.8,0.2],1,2);
bottomtranslate = [0.9*ones(1,4) , 1.4*ones(1,4)];

for i=1:numcurvesplot
    ax1 = axes;
    subplot(2,4,i,ax1)
    im1 = imagesc(ax1,flipud(density1{i}),[0,ylim]); % APC->red
    axis xy

    ax1.Position = subplotposition(i,:);
    
    if i==4
        h1 = colorbar;
        set(h1,'FontSize',fscolorbar,...
            'Position',[0.939558442449601 0.515625 0.0170138888888887 0.440972222222222]);
        text(1.36*gridpts,1.05*gridpts,'APCs (cells/mm$^2$)',...
            'HorizontalAlignment','right','FontSize',fscolorbar,...
            'Interpreter','latex')
    end

    ax2 = axes;
    subplot(2,4,i,ax2)
    im2 = imagesc(ax2,flipud(density2{i}),[0,ylim]); % IPA->blue
    alpha(im2,0.5);
    axis xy
    ax2.Position = subplotposition(i,:);
    
    if i==8
        h2 = colorbar;
        set(h2,'FontSize',fscolorbar,...
            'Position',[0.939558442449601 0.015625 0.0170138888888887 0.439236111111111]);
        text(1.36*gridpts,1.05*gridpts,'IPAs (cells/mm$^2$)',...
            'HorizontalAlignment','right','FontSize',fscolorbar,...
            'Interpreter','latex')
    end

    % scale bar
    if i==1
        hold on
        plot([0.12*gridpts+0; 0.12*gridpts+gridpts/(2*rmax_plot)], ...
            [0.15*gridpts; 0.15*gridpts],'-k','LineWidth', 5)
        hold off
        text(0.12*gridpts+gridpts/(4*rmax_plot),0.07*gridpts, '1 mm', 'HorizontalAlignment','center','FontSize',fsticks)
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
    ax1.XTick = linspace(0,gridpts,2*rmax_plot+1);
    ax1.XTickLabel = [];
    ax1.YLim = ax1.XLim;
    ax1.YTick = ax1.XTick;
    ax1.YTickLabel = ax1.XTickLabel;
end

set(gcf,'Units','inches','Position',[1,1,12,6],'PaperPositionMode','auto',...
    'PaperOrientation','landscape',...
    'PaperSize',[12,6]);