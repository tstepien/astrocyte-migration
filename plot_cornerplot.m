clear variables;
clc;
close all;

addpath emcee_mymod

numparametersets = 1000000;
load('parameter_analysis/latin_hypercube/latinhypercube_1000000pts.mat')

param_names = {'$\mu$','$\alpha_{11}$','$\alpha_{12}$','$\alpha_{21}$',...
    '$\alpha_{22}$','$\beta_1$','$\beta_2$','$\beta_3$','$\gamma_1$',...
    '$\gamma_2$','$T_e$','$P_\mathrm{hy}$','$r_\mathrm{hy}$'};

p = [mu , alpha11 , alpha12 , alpha21 , alpha22 , beta1 , ...
    beta2 , beta3 , gamma1 , gamma2 , Te , P_hy , r_hy];

err_original = [err_dens err_rad err_time err_tot];
err_names = {'Density Error','Radius Error','Time Error','Total Error'};

clear err_dens err_rad err_time err_tot mu alpha11 alpha12 alpha21 alpha22 ...
    beta1 beta2 beta3 gamma1 gamma2 Te P_hy r_hy;

%% remove errors that were set to 10^4
maxthreshold = 10^4;

ind = (1:N)';
ind_maxthreshold = ind(err_original(:,4) < maxthreshold);
num_maxthreshold = length(ind_maxthreshold);

err_maxthreshold = err_original(ind_maxthreshold,:);

[~,ind_sort_maxthreshold] = sort(err_maxthreshold(:,4));
err_maxthreshold_sort = err_maxthreshold(ind_sort_maxthreshold,:);

%% corner plot
minnumkeep = 3; %number to keep for scatter plot
[sortedValues,sortIndex] = sort(err_tot);
minIndex = sortIndex(1:minnumkeep);

minparam = p(minIndex,:);
avgval = [mean(p(:,1)) , mean(p(:,2)) , mean(p(:,3)) , mean(p(:,4)), ...
    mean(p(:,5)) , mean(p(:,6)) , mean(p(:,7)) , mean(p(:,8)), ...
    mean(p(:,9)) , mean(p(:,10)) , mean(p(:,11)) , mean(p(:,12)), ...
    mean(p(:,13))];

figure
H = ecornerplot(p,minparam,avgval,'names',param_names,'ks',true);



set(gcf,'Units','inches','Position',[2,2,10,8],'PaperPositionMode','auto')