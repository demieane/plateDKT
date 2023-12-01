%% VISUALIZING LANGRANGE SHAPE FUNCTIONS
close all;clear all;
Nmesh = 100;
[xx,yy]=meshgrid(linspace(0,1,Nmesh),linspace(0,1,Nmesh));

% bx1, by1, ... (on the three-nodes)

nodal_b = zeros(1,6)';%column

figure; 
for kk = 1:6
    nodal_b(kk) = 1;
    for i=1:Nmesh
        for j = 1:Nmesh
            [ SF,~,~,~,~,~] = LNShapeFunDST(xx(i,j),yy(i,j));%row
            if (xx(i,j)+yy(i,j)>1)
                LLnew(i,j)=NaN;
            else
                LLnew(i,j) = SF*nodal_b;
            end
        end
    end
    nodal_b(kk) = 0;
    subplot(2,3,kk)
    surf(xx,yy,LLnew);
    shading interp;
    grid minor;
    colormap(viridis);
    xlabel('x','Interpreter','latex','FontSize',14);
    ylabel('y','Interpreter','latex','FontSize',14);
    %title('Lagrange Shape Function');
    set(gca,'FontSize',14);
end
