Nx = 10;
Ny = 10;
[xdummy,ydummy]=meshgrid(linspace(-1,1,Nx),linspace(-1,1,Ny));
s=0.2;
m=0.6;
fx=1./s/sqrt(2*pi).*exp(-0.5*((xdummy-m)./s).^2); 
a = chord;
b = span;
xx = a*0.5*(-xdummy+1);
yy = b*0.5*(ydummy+1);
        
figure
plot(xx,yy,'o');
hold on;
plot(xx(1,:),yy(1,:),'s');
xlabel('x');
ylabel('y');

        
% %         chord = 10;
% %         NACA = 0012;
% %         [TAU,EPSMAX,PTMAX] = indata(NACA);
% %         Npanels = 50;
% %         [Nodtot,X,Y] = setup(Npanels,TAU,EPSMAX,NACA,PTMAX);
% %         figure
% %         plot(X*chord, Y*chord,'b');
        
        
        xFoil = (X(2:Npanels/2)-0.5)./0.5; %xdummy 
        thickFoil = -Y(2:Npanels/2); %ydummy
        
        modfunc = (xFoil).^2;
        modfunc = modfunc./max(abs(modfunc));
        tx_foil_mod = thickFoil+ 1/100.*modfunc;
        
        figure
        plot(X,Y,'ro');
        hold on;grid on;
        plot(xFoil, thickFoil,'ks');
        plot(xFoil, tx_foil_mod,'b^');
        axis equal;
       
        thickdummy = xdummy.*0;
        for ii = 1:size(xdummy,2)
            thickdummy(ii,:) = interp1(-xFoil, tx_foil_mod, xdummy(ii,:),'linear','extrap');
        end
        
        figure
        surf(xdummy,ydummy,thickdummy);
       
        % thickness constant
        h0 = 0.001*chord;
        tx_modified = xx.*0 + h0;
        
        for ii = 1:size(xx,1)
            for jj = 1:size(xx,2)
               tx_modified(ii,jj) = thickdummy(ii,jj);
            end
        end
        
        figure;
        colormap(viridis);
        subplot(2,2,[ 1 2]);
        plot(xx(1,:),tx_modified(1,:),'ro');
        hold on;grid on;
%         h1=plot(xm,txx,'b^','MarkerSize',2);
%         plot(xc_fem_data(1,:),tx_linear(1,:),'k^-');
        legend("original", "interpolated")
        xlabel('x [m]');ylabel('y [m]');
        title('plate thickness spanwise (from mat-file)');
           subplot(2,2,3);
        hold on;grid on;
        surf(xx,yy,fx);
%         plot3(xm,ym,fxx,'b^', 'MarkerSize',2);
        title('bem data for fx');
        xlabel('x-axis');
        ylabel('y-axis');
        zlabel('P [Pa]');view([-25 25])

        subplot(2,2,4);hold on;grid on;
        surf(xx,yy,tx_modified);
%         plot3(xm,ym,txx,'b^', 'MarkerSize',2);
        title('bem data for tx');
        xlabel('x-axis');
        ylabel('y-axis');view([-25 25])

        