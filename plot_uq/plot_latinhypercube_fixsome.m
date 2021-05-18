clear variables global;
clc;

addpath emcee_mymod

percentholdon = 0.1;
what_set = 'maxthreshold'; %'maxthreshold' or 'maxmode'
fit_dist_plot = 'no'; % using percentholdon = 0.01 for distribution fits
titles_on = 'no';

load('parameter_analysis/2pop-someparamfixed/latin_hypercube/latinhypercube_100000pts.mat')

err_original = [err_dens err_rad err_time err_tot];
err_names = {'Density Error','Radius Error','Time Error','Total Error'};

% errDensity.original = err_dens;
% errRadius.original = err_rad;
% errTime.original = err_time;
% errTotal.original = err_tot;

param_original = [mu , alpha11 , alpha12 , alpha22 , beta1 , beta3 , ...
    gamma2 , P_hy , r_hy];

clear err_dens err_rad err_time err_tot mu alpha11 alpha12 alpha22 ...
    beta1 beta3 gamma2 P_hy r_hy;

param_names = {'$\mu$','$\alpha_{11}$','$\alpha_{12}$','$\alpha_{22}$',...
    '$\beta_1$','$\beta_3$','$\gamma_2$','$P_\mathrm{hy}$','$r_\mathrm{hy}$'};
num_param = length(param_names);
param_names_words = {'Adhesion constant','APC prolif rate wrt O_2',...
    'APC prolif rate wrt PDGFA','IPA prolif rate wrt O_2',...
    'Mass action rate','Differentiation rate wrt LIF',...
    'APC apoptosis rate','IPA apoptosis rate',...
    'Hyaloid artery maximum','Hyaloid artery half-max value'};

%% remove errors that were set to 10^4
maxthreshold = 10^4;

ind = (1:N)';
ind_maxthreshold = ind(err_original(:,4) < maxthreshold);
num_maxthreshold = length(ind_maxthreshold);

err_maxthreshold = err_original(ind_maxthreshold,:);

[~,ind_sort_maxthreshold] = sort(err_maxthreshold(:,4));
err_maxthreshold_sort = err_maxthreshold(ind_sort_maxthreshold,:);

fig1 = figure;
tiledlayout(2,2)
for i=1:4
    nexttile
    scatter(1:num_maxthreshold,err_maxthreshold_sort(:,i))
    xlim([0,num_maxthreshold])
    xlabel(err_names{i})
end
sgtitle(strcat(['Errors <10^4 (',num2str(num_maxthreshold),' parameter sets)']))

modes_error = zeros(1,4);
for i=1:4
    modes_error(i) = mode(err_maxthreshold(:,i));
end

%% look at errors that are smaller than the mode errors for density, radius, and time

% ind_maxmode = ind( err_original(:,1) < modes_error(1) ...
%     & err_original(:,2) < modes_error(2) ...
%     & err_original(:,3) < modes_error(3) );
% num_maxmode = length(ind_maxmode);
% 
% err_maxmode = err_original(ind_maxmode,:);
% 
% [~,ind_sort_maxmode] = sort(err_maxmode(:,4));
% err_maxmode_sort = err_maxmode(ind_sort_maxmode,:);
% 
% fig2 = figure;
% tiledlayout(2,2)
% for i=1:4
%     nexttile
%     scatter(1:num_maxmode,err_maxmode_sort(:,i))
%     xlim([0,num_maxmode])
%     xlabel(err_names{i})
% end
% sgtitle(strcat(['Errors < modes for density/radius/time (',num2str(num_maxmode),' parameter sets)']))

%% histograms of parameters

if strcmp(what_set,'maxmode')==1
    num_parametersets = num_maxmode;
    ind_parametersets = ind_maxmode;
    ind_sort = ind_sort_maxmode;
elseif strcmp(what_set,'maxthreshold')==1
    num_parametersets = num_maxthreshold;
    ind_parametersets = ind_maxthreshold;
    ind_sort = ind_sort_maxthreshold;
end

num_hold = ceil(percentholdon * num_parametersets);

param_sort = zeros(num_parametersets,num_param);
param_sort_hold = zeros(num_hold,num_param);
for i = 1:num_param
    param_sort(:,i) = param_original(ind_parametersets,i);
    param_sort_hold(:,i) = param_sort(ind_sort(1:num_hold),i);
end

fig3 = figure;
tiledlayout(2,5,'TileSpacing','compact','Padding','compact')
for i=1:num_param
    if i<5
        nexttile
    else
        nexttile(i+1)
    end
    
    if strcmp(fit_dist_plot,'no')==1
        histogram(param_sort_hold(:,i),'Normalization','probability',...
            'BinMethod','sturges','FaceColor','none','LineWidth',1.5);
    elseif strcmp(fit_dist_plot,'yes')==1
        distributionfit = {'exponential','normal','uniform','exponential','uniform',...
            'uniform','exponential','uniform','exponential','normal',...
            'normal','normal','normal'};
        if strcmp(distributionfit{i},'uniform')==0
            h = histfit(param_sort_hold(:,i),[],distributionfit{i});
            h(1).FaceColor = 'none';
            h(2).Color = 'k';
            box on
            yt = get(gca,'YTick');
            set(gca,'YTick',yt,'YTickLabel',round(yt/num_hold,2));
        else
            hold on
            histogram(param_sort_hold(:,i),'Normalization','probability',...
                'BinMethod','sqrt','FaceColor','none');
            plot(linspace(bound(i,1),bound(i,2),100),0.02*ones(1,100),'k',...
                'LineWidth',2.5)
            box on
            hold off
        end
    end
    
    xlabel(param_names{i},'Interpreter','latex')
    if strcmp(titles_on,'yes')==1
        title(param_names_words{i},'FontWeight','normal')
    end
    if i==1 || i==5
        ylabel('Percentage','Interpreter','latex')
    end
    xlim([0,bound(i,2)])
    
    set(gca,'FontSize',14)
end

if strcmp(titles_on,'yes')==1
    sgtitle(strcat(['Smallest ',num2str(percentholdon*100),'% Error (',num2str(num_hold),' parameter sets)']))
end

set(fig3,'Units','inches','Position',[2,2,16,4.74],'PaperPositionMode','auto')

%% determine type of distribution

disttype = {'normal';'lognormal';'gamma';'exponential';'weibull';'logistic'};
%%% didn't use these distributions:
%%% 'beta';'birnbaumsaunders';'burr';'negative binomial';'extreme value';'kernel';
%%% 'generalized extreme value';'generalized pareto';'inversegaussian';
%%% 'nakagami';'loglogistic';'poisson';'rayleigh';'rician';'tlocationscale';
num_dist = length(disttype);

param_dist = cell(num_dist,num_param);
GoF_dist = zeros(num_dist,num_param);

for i=1:num_dist
    for j=1:num_param
        param_dist{i,j} = fitdist(param_sort_hold(:,j),disttype{i});
        GoF_dist(i,j) = chi2gof(param_sort_hold(:,j),'CDF',param_dist{i});
    end
end

%% corner plot

num_keepscatter = 3; %number to keep for scatter plot
param_min = param_sort_hold(1:num_keepscatter,:);

param_mean = zeros(1,num_param);
param_mode = zeros(1,num_param);
for i=1:num_param
    param_mean(i) = mean(param_sort_hold(:,i));
    param_mode(i) = mode(param_sort_hold(:,i));
end

fig4 = figure;
ecornerplot(param_sort_hold,param_min,param_mean,bound,'names',param_names,'ks',true);
set(fig4,'Units','inches','Position',[2,2,10,8],'PaperPositionMode','auto')

%% correlation

corrmatrix = zeros(num_param,num_param);
for i=1:num_param
    for j=1:num_param
        if j<i
            corrmatrix(i,j) = corr(param_sort_hold(:,i),param_sort_hold(:,j));
        end
    end
end

max(abs(corrmatrix(:)))