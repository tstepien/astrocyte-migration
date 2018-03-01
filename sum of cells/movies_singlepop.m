nFrames = length(t);
mov(1:nFrames) = struct('cdata',[],'colormap',[]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% movies %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
set(gcf,'Position',[440,378,800,400]);

h1 = [];
h2 = [];

for i = 1:length(t)
    set(gca,'nextplot','replacechildren');
    delete(h1);
    delete(h2);
    
%%%%%%% astrocytes
    subaxis(1,2,1,'MarginBottom',0.15,'MarginLeft',0.07,'MarginRight',0.05)
    plot(r,c(i,:),'LineWidth',1.5)
    set(gca,'FontSize',18,'ylim',[0,4000]);
    xlabel('$r$','FontSize',26,'Interpreter','latex')
    
    title('Astrocytes','FontSize',24,'FontWeight','bold')%,...
%         'Units','normalized','HorizontalAlignment', 'right',...
%         'Position',[0.325,1.02]);
    
%%%%%%% growth factors
    subaxis(1,2,2,'MarginRight',0.01)
    plot(r,p1(i,:),'LineWidth',1.5)
    set(gca,'FontSize',18,'ylim',[0,100]);
    xlabel('$r$','FontSize',26,'Interpreter','latex')
    
    title('PDGF-A','FontSize',24,'FontWeight','bold')%,...
%         'Units','normalized','HorizontalAlignment', 'right',...
%         'Position',[0.325,1.02]);
    
    set(gcf,'color','w');
    mov(i) = getframe(gcf);
end

v = VideoWriter('astrocytesmov.avi','Motion JPEG AVI');%,'Uncompressed AVI');
v.FrameRate = 15;
open(v);
writeVideo(v,mov)
close(v);