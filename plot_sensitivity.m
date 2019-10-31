clear variables global;
clc;

load('sensitivity analysis results/sensitivity_analysis_1031_drp01_morepoints.mat')

param_names = {'$\mu$','$\alpha_1$','$\alpha_2$','$\beta$',...
    '$\gamma_1$','$\gamma_2$','$T_e$'};

figure
for i=1:7
    if i==1 || i==5
        subaxis(2,4,i,'SpacingVert',0.1,'MarginLeft',0.06,'MarginRight',0)
    elseif i==4
        subaxis(2,4,i,'SpacingVert',0.1,'MarginRight',0.01,'MarginLeft',0)
    else
        subaxis(2,4,i,'SpacingVert',0.1,'MarginLeft',0.04,'MarginRight',0.01)
    end
    h = plot(intrange(i,:),err(i,:),'-o');
    xlabel(param_names{i},'Interpreter','latex')
    if i==1 || i==5
        ylabel('error')
    end
    set(gca,'FontSize',14)
end

set(gcf,'Units','inches','Position',[2,2,16,8],'PaperPositionMode','auto')