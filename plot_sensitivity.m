clear variables global;
clc;

load('sensitivity_analysis_1007.mat')

param_names = {'$\kappa$','$\mu$','$\alpha_1$','$\alpha_2$','$\beta$',...
    '$\gamma_1$','$\gamma_2$','$\bar{\xi}_1$','$\bar{\xi}_2$','$T_e$'};

figure
for i=1:10
    subaxis(2,5,i)
    h = plot(intrange(i,:),err(i,:),'-o');
    xlabel(param_names{i},'Interpreter','latex')
    if i==1 || i==6
        ylabel('error')
    end
end