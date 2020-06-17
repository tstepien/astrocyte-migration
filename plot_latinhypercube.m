clear variables global;
clc;

load('parameter_analysis/latinhypercube_100000pts.mat')

param_names = {'$\mu$','$\alpha_{11}$','$\alpha_{12}$','$\alpha_{21}$',...
    '$\alpha_{22}$','$\beta_1$','$\beta_2$','$\beta_3$','$\gamma_1$',...
    '$\gamma_2$','$T_e$','$P_\mathrm{hy}$','$r_\mathrm{hy}$'};


% error_threshold_tot = 20;
error_threshold_time = 0.7;

ind = (1:N)';
% ind_good = ind(err_tot < error_threshold_tot & err_time < error_threshold_time);
% ind_bad = ind(err_tot >= error_threshold_tot & err_time >= error_threshold_time);
ind_good = ind(err_time < error_threshold_time);
ind_bad = ind(err_time >= error_threshold_time);

param_good = [mu(ind_good) alpha11(ind_good) alpha12(ind_good) ...
    alpha21(ind_good) alpha22(ind_good) beta1(ind_good) beta2(ind_good) ...
    beta3(ind_good) gamma1(ind_good) gamma2(ind_good) Te(ind_good) ...
    P_hy(ind_good) r_hy(ind_good)];
param_bad = [mu(ind_bad) alpha11(ind_bad) alpha12(ind_bad) ...
    alpha21(ind_bad) alpha22(ind_bad) beta1(ind_bad) beta2(ind_bad) ...
    beta3(ind_bad) gamma1(ind_bad) gamma2(ind_bad) Te(ind_bad) ...
    P_hy(ind_bad) r_hy(ind_bad)];

figure
tiledlayout(3,5,'TileSpacing','compact','Padding','compact')
for i=1:length(param_names)
    nexttile
    
    h = histogram(param_good(:,i),'Normalization','probability',...
        'FaceColor','none','LineWidth',1.5);
%     h = histfit(param_good(:,i),[],'normal');
%     h(1).FaceColor = 'none';
%     h(2).Color = 'k';
    
    xlabel(param_names{i},'Interpreter','latex')
    if i==1 || i==6 || i==11
        ylabel('Percentage','Interpreter','latex')
    end
    if i==1
        title('Adhesion constant','FontWeight','normal')
    elseif i==2
        title('APC prolif rate wrt O_2','FontWeight','normal')
    elseif i==3
        title('APC prolif rate wrt PDGFA','FontWeight','normal')
    elseif i==4
        title('IPA prolif rate wrt O_2','FontWeight','normal')
    elseif i==5
        title('IPA prolif rate wrt PDGFA','FontWeight','normal')
    elseif i==6
        title('Mass action rate','FontWeight','normal')
    elseif i==7
        title('Differentiation rate wrt O_2','FontWeight','normal')
    elseif i==8
        title('Differentiation rate wrt LIF','FontWeight','normal')
    elseif i==9
        title('APC apoptosis rate','FontWeight','normal')
    elseif i==10
        title('IPA apoptosis rate','FontWeight','normal')
    elseif i==11
        title('Edge tension','FontWeight','normal')
    elseif i==12
        title('Hyaloid artery maximum','FontWeight','normal')
    elseif i==13
        title('Hyaloid artery half-max value','FontWeight','normal')
    end
    xlim([0,bound(i,2)])
    
    set(gca,'FontSize',14)
end

% sgtitle(strcat('\bf Total Error Threshold:',' ',num2str(error_threshold_tot),...
%     ',   \bf Time Error Threshold:',' ',num2str(error_threshold_time)))
sgtitle(strcat('\bf Time Error Threshold:',' ',num2str(error_threshold_time)))

set(gcf,'Units','inches','Position',[2,2,16,7],'PaperPositionMode','auto')

%% determine type of distribution

% disttype = {'normal';'lognormal';'gamma';'exponential';'weibull';'logistic'};
% %%% didn't use these distributions:
% %%% 'beta';'birnbaumsaunders';'burr';'negative binomial';'extreme value';'kernel';
% %%% 'generalized extreme value';'generalized pareto';'inversegaussian';
% %%% 'nakagami';'loglogistic';'poisson';'rayleigh';'rician';'tlocationscale';
% numDist = length(disttype);
% 
% param_dist = cell(numDist,numpar);
% GoF_dist = zeros(numDist,numpar);
% 
% for i=1:numDist
%     for j=1:numpar
%         param_dist{i,j} = fitdist(param_good(:,j),disttype{i});
%         GoF_dist(i,j) = chi2gof(param_good(:,j),'CDF',param_dist{i});
%     end
% end