clear variables global;
clc;

load('sensitivity analysis results/sensitivity_analysis_1011.mat')

param_names = {'$\mu$','$\alpha_1$','$\alpha_2$','$\beta$',...
    '$\gamma_1$','$\gamma_2$','$\bar{\xi}_1$','$\bar{\xi}_2$','$T_e$'};

figure
for i=1:9
    if i==1 || i==6
        subaxis(2,5,i,'SpacingVert',0.1,'MarginLeft',0.06,'MarginRight',0)
    elseif i==5
        subaxis(2,5,i,'SpacingVert',0.1,'MarginRight',0.01,'MarginLeft',0)
    else
        subaxis(2,5,i,'SpacingVert',0.1,'MarginLeft',0.04,'MarginRight',0.01)
    end
    h = plot(intrange(i,:),err(i,:),'-o');
    xlabel(param_names{i},'Interpreter','latex')
    if i==1 || i==6
        ylabel('error')
    end
    set(gca,'FontSize',14)
end

set(gcf,'Units','inches','Position',[2,2,16,8],'PaperPositionMode','auto')