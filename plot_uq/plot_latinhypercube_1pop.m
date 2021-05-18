clear variables global;
clc;

load('parameter_analysis/1pop/latin_hypercube/latinhypercube_10000pts.mat')

param_names = {'$\mu$','$\alpha_{11}$','$\alpha_{12}$',...
    '$\gamma_1$','$T_e$','$P_\mathrm{hy}$','$r_\mathrm{hy}$'};


error_threshold = 0.5;

ind = (1:N)';
ind_good = ind(err_time < error_threshold);
ind_bad = ind(err_time >= error_threshold);

param_good = [mu(ind_good) alpha11(ind_good) alpha12(ind_good) ...
    gamma1(ind_good) Te(ind_good) P_hy(ind_good) r_hy(ind_good)];
param_bad = [mu(ind_bad) alpha11(ind_bad) alpha12(ind_bad) ...
    gamma1(ind_bad) Te(ind_bad) P_hy(ind_bad) r_hy(ind_bad)];

figure
for i=1:length(param_names)
    if i==1
        subaxis(2,4,i,'SpacingVert',0.15,'MarginLeft',0.05,'MarginRight',0,'MarginTop',0.15,'MarginBottom',0.08)
    elseif i==5 || i==9
        subaxis(2,4,i,'SpacingVert',0.15,'MarginLeft',0.05,'MarginRight',0,'MarginBottom',0.09)
    elseif i==4 || i==8
        subaxis(2,4,i,'SpacingVert',0.15,'MarginRight',0.02,'MarginLeft',0)
    else
        subaxis(2,4,i,'SpacingVert',0.15,'MarginLeft',0.04,'MarginRight',0.01)
    end
    
    h = histogram(param_good(:,i),'Normalization','probability',...
        'FaceColor','none','LineWidth',1.5);
%     h = histfit(param_good(:,i),[],'normal');
%     h(1).FaceColor = 'none';
%     h(2).Color = 'k';
    
    xlabel(param_names{i},'Interpreter','latex')
    if i==1 || i==5
        ylabel('Percentage','Interpreter','latex')
    end
    if i==1
        title('Adhesion constant','FontWeight','normal')
    elseif i==2
        title('APC proliferation rate wrt O_2','FontWeight','normal')
    elseif i==3
        title('APC proliferation rate wrt PDGFA','FontWeight','normal')
    elseif i==4
        title('APC apoptosis rate','FontWeight','normal')
    elseif i==5
        title('Edge tension','FontWeight','normal')
    elseif i==6
        title('Hyaloid artery maximum','FontWeight','normal')
    elseif i==7
        title('Hyaloid artery half-max value','FontWeight','normal')
    end
    xlim([0,bound(i,2)])
    
    set(gca,'FontSize',14)
end

sgtitle(strcat('\bf Time Error Threshold:',' ',num2str(error_threshold)))

set(gcf,'Units','inches','Position',[2,2,16,7],'PaperPositionMode','auto')

%% determine type of distribution

disttype = {'normal';'lognormal';'gamma';'exponential';'weibull';'logistic'};
%%% didn't use these distributions:
%%% 'beta';'birnbaumsaunders';'burr';'negative binomial';'extreme value';'kernel';
%%% 'generalized extreme value';'generalized pareto';'inversegaussian';
%%% 'nakagami';'loglogistic';'poisson';'rayleigh';'rician';'tlocationscale';
numDist = length(disttype);

param_dist = cell(numDist,numpar);
GoF_dist = zeros(numDist,numpar);

for i=1:numDist
    for j=1:numpar
        param_dist{i,j} = fitdist(param_good(:,j),disttype{i});
        GoF_dist(i,j) = chi2gof(param_good(:,j),'CDF',param_dist{i});
    end
end