clear variables global;
clc;

load('sensitivity analysis results/sensitivity_analysis_1115.mat')

param_names = {'$\mu$','$\alpha_1$','$\alpha_2$','$\beta$',...
    '$\gamma_1$','$\gamma_2$','$T_e$'};

[rownan,colnan] = find(isnan(err));

figure
for i=1:7
    if i==1
        subaxis(2,4,i,'SpacingVert',0.1,'MarginLeft',0.06,'MarginRight',0,'MarginTop',0.05)
    elseif i==5
        subaxis(2,4,i,'SpacingVert',0.12,'MarginLeft',0.06,'MarginRight',0,'MarginBottom',0.08)
    elseif i==4
        subaxis(2,4,i,'SpacingVert',0.1,'MarginRight',0.01,'MarginLeft',0)
    else
        subaxis(2,4,i,'SpacingVert',0.12,'MarginLeft',0.04,'MarginRight',0.01)
    end
    hold on
    h = plot(intrange(i,:),err(i,:),'-o','LineWidth',1.5);
    if sum(i==rownan)>0
        markplotx = intrange(i,colnan(i==rownan));
        markploty = min(err(i,:))*ones(size(markplotx));
        scatter(markplotx,markploty,100,'x','LineWidth',1.5,...
            'MarkerEdgeColor',[0.8500 0.3250 0.0980])
    end
    hold off
    xlabel(param_names{i},'Interpreter','latex')
    if i==1 || i==5
        ylabel('error')
    end
    if i==1
        title('Adhesion constant')
    elseif i==2
        title('APC proliferation rate')
    elseif i==3
        title('IPA proliferation rate')
    elseif i==4
        title('Differentiation rate')
    elseif i==5
        title('APC apoptosis rate')
    elseif i==6
        title('IPA apoptosis rate')
    elseif i==7
        title('Edge tension')
    end
    xlim(bound(i,:))
    set(gca,'FontSize',14)
end

set(gcf,'Units','inches','Position',[2,2,16,8],'PaperPositionMode','auto')