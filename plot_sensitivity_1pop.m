clear variables global;
clc;

load('paramter_analysis/sensitivity_analysis_0518_4.mat')

param_names = {'$\mu$','$\alpha_{11}$','$\alpha_{12}$',...
    '$\gamma_1$','$T_e$','$P_\mathrm{hy}$','$r_\mathrm{hy}$'};

[rownan,colnan] = find(isnan(err));

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
    hold on
    h = plot(intrange(i,:),err(i,:),'-o','LineWidth',1.5);
    if sum(i==rownan)>0
        markplotx = intrange(i,colnan(i==rownan));
        markploty = min(err(i,:))*ones(size(markplotx));
        scatter(markplotx,markploty,100,'x','LineWidth',1.5,...
            'MarkerEdgeColor',[0.8500 0.3250 0.0980])
    end
    
    if i==1
        paramofinterest = p.mu;
    elseif i==2
        paramofinterest = p.alpha11;
    elseif i==3
        paramofinterest = p.alpha12;
    elseif i==4
        paramofinterest = p.gamma1;
    elseif i==5
        paramofinterest = p.Te;
    elseif i==6
        paramofinterest = p.P_hy;
    elseif i==7
        paramofinterest = p.r_hy;
    end
    b1 = find(paramofinterest>=intrange(i,:),1,'last');
    b2 = find(paramofinterest<intrange(i,:),1,'first');
    baselineyval = interp1([intrange(i,b1),intrange(i,b2)],...
        [err(i,b1),err(i,b2)],paramofinterest);
    scatter(paramofinterest,baselineyval,100,'*','LineWidth',1.5,...
            'MarkerEdgeColor',[0.4940 0.1840 0.5560])
    
    hold off
    xlabel(param_names{i},'Interpreter','latex')
    if i==1 || i==5
        ylabel('error')
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
    if isequal(min(err(i,:)),max(err(i,:)))==0
        ylim([min(err(i,:)),max(err(i,:))])
    end
    set(gca,'FontSize',14)
end

set(gcf,'Units','inches','Position',[2,2,16,7],'PaperPositionMode','auto')