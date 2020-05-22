clear variables global;
clc;

load('sensitivity analysis results/results_latinhypercube_5000pts.mat')

param_names = {'$\mu$','$\alpha_{11}$','$\alpha_{12}$',...
    '$\gamma_1$','$T_e$','$P_\mathrm{hy}$','$r_\mathrm{hy}$'};


error_threshold = 0.7;

ind = (1:N)';
ind_good = ind(err_time < error_threshold);
ind_bad = ind(err_time >= error_threshold);

param_good = [mu(ind_good) alpha11(ind_good) alpha12(ind_good) ...
    gamma1(ind_good) Te(ind_good) P_hy(ind_good) r_hy(ind_good)];
param_bad = [mu(ind_bad) alpha11(ind_bad) alpha12(ind_bad) ...
    gamma1(ind_bad) Te(ind_bad) P_hy(ind_bad) r_hy(ind_bad)];

figure
for i=1:7
    if i==1
        subaxis(2,4,i,'SpacingVert',0.1,'MarginLeft',0.045,'MarginRight',0,'MarginTop',0.07,'MarginBottom',0.08)
    elseif i==5 || i==9
        subaxis(2,4,i,'SpacingVert',0.12,'MarginLeft',0.045,'MarginRight',0,'MarginBottom',0.09)
    elseif i==4 || i==8
        subaxis(2,4,i,'SpacingVert',0.1,'MarginRight',0.02,'MarginLeft',0)
    else
        subaxis(2,4,i,'SpacingVert',0.12,'MarginLeft',0.04,'MarginRight',0.01)
    end
    
    h = histogram(param_good(:,i),'Normalization','probability',...
        'FaceColor','none','LineWidth',1.5);
    
    xlabel(param_names{i},'Interpreter','latex')
    if i==1 || i==5
        ylabel('Percentage','Interpreter','latex')
    end
    if i==1
        title('Adhesion constant')
    elseif i==2
        title('APC proliferation rate wrt O_2')
    elseif i==3
        title('APC proliferation rate wrt PDGFA')
    elseif i==4
        title('APC apoptosis rate')
    elseif i==5
        title('Edge tension')
    elseif i==6
        title('Hyaloid artery maximum')
    elseif i==7
        title('Hyaloid artery half-max value')
    end
    xlim([0,bound(i,2)])
    
    set(gca,'FontSize',14)
end

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