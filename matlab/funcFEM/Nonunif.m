function [fxx,txx,xx,yy]=Nonunif(x,y,IEN, p,e,t, chord, span, debugOn, importFromFile, fluid_dens, Uvel, h, d, triangleData)

% N1=50;
% N2=25;
% % [xx,yy]=meshgrid(linspace(-1,1,N),linspace(-1,1,N));
% % % [xx,yy]=meshgrid(linspace(0,chord,N),linspace(0,span,N));
% [xx,yy]=meshgrid(linspace(0,chord,N1),linspace(0,span,N2));
% xxx=xx(1,:);
% yyy=yy(:,1);
% s=0.2;
% m=0.5;
% s=0.05;
% m=0.25;
% fx=1./s/sqrt(2*pi).*exp(-0.5*((xx-m)./s).^2); 
% fx=1200*fx/max(max(fx));

if triangleData.ImportFromMatFile
    
    if importFromFile.toggle
        load(importFromFile.filename, 'xc_fem_data', 'yc_fem_data',...
            'th_fem_data', 'DCoefpres', 'inData');
    else
        error('Provide file to read from!')
    end
    xx=xc_fem_data;
    yy=yc_fem_data;

    fx=DCoefpres(:,:,d)*(0.5*fluid_dens*Uvel^2); %dimensionalize data [N]
    
%     Area = chord*span; % IS THIS OK?
%     fx=fx/Area;
%     fx=DCoefpres(:,:,d)*(0.5*fluid_dens*Uvel^2); %dimensionalize data
    
    a3=inData.a3;%0.0175;
    phase3=inData.phase3;%0;
    w3=inData.omega3;%10.92;
    dt=inData.dt;%0.0027;
    heave=a3*sin(w3*d*dt+phase3);
    dheave=w3*a3*cos(w3*d*dt+phase3);
    ddheave=-w3^2*a3*sin(w3*d*dt+phase3);
    
    fx = fx  - 0.*7850*h*ddheave; %FIX
    
%     tx_modified= 0.*th_fem_data + h; 
    tx_modified= th_fem_data + 2.5*chord.*(xc_fem_data-chord/2).^2;

    prepareLoad_csv(xx,yy,fx,'load_fx.csv',1); %prepare .csv file for ANSYS
    prepareLoad_csv(xx,yy,tx_modified,'sections_tx.csv',2); %prepare .csv file for ANSYS
    
    if debugOn
        figure
        subplot(2,2,[ 1 2]);
        plot(xc_fem_data(1,:),th_fem_data(1,:),'ro');
        hold on;grid on;
        plot(xc_fem_data(1,:),tx_modified(1,:),'b^');
        legend("original", "modified")
        xlabel('x [m]');ylabel('y [m]');
        title('plate thickness spanwise (from mat-file)');
        %
        subplot(2,2,3);
        hold on;grid on;
        surf(xx,yy,fx);
        title('bem data for fx');
        xlabel('x-axis');
        ylabel('y-axis');
        zlabel('P [Pa]');view([-25 25])
        %
        subplot(2,2,4);hold on;grid on;
        surf(xx,yy,tx_modified);
        title('bem data for tx');
        xlabel('x-axis');
        ylabel('y-axis');view([-25 25])
    end

    
elseif triangleData.VerificationPhD
    
   %% BENCHMARKS
   KATSIKADELIS = 0;
   SHUFFRIN = 0;
   DISTRIBUTED_LOAD = 1; % analytic formula
   
   if KATSIKADELIS 
       % KATSIKADELIS C-C-C-C
       % Katsikadelis & Sapoutzakis : A BEM solution to dynamic analysis of
       % plates 1991 Computational Mechanics, 7, 369-379
       % linear thickness along the x-axis 
       numSize = 150; 
       [xx,yy]=meshgrid(linspace(0,chord,numSize),linspace(0,span,numSize));
       alphaPachenari = 0.8;%0.6;%0.6;%0.8;
       alpha = (1-alphaPachenari);
       for i=1:numSize
           for j=1:numSize
               fx(i,j)=1.0;
               h1=0.001*span;
               tx_modified(i,j) = h1*(1+alpha*xx(i,j)/chord);
           end
       end
       taper_ratio_katsikadelis = max(max(tx_modified))/min(min(tx_modified))-1
       
   elseif SHUFFRIN == 1
       % SHUFFRIN C-F-F-F
       % linear thickness along the x-axis 
       numSize = 100; 
       [xx,yy]=meshgrid(linspace(0,chord,numSize),linspace(0,span,numSize));
       [xxNUM, yyNUM]=meshgrid(linspace(0,1,numSize),linspace(0,1,numSize));
       alphaShufrin = 0.5;% 0.25;%0.5;
       for i=1:numSize
           for j=1:numSize
               fx(i,j)=1.0;
               h0=0.1*span;%0.001*span;
               tx_modified(i,j) = h0*(1-alphaShufrin*xxNUM(i,j));%.*(1-alphaShufrin*yyNUM(i,j));
%                tx_modified(i,j) = h0*(1-alphaShufrin*xx(i,j)./chord);%.*(1-alphaShufrin*yyNUM(i,j));
           end
       end
       
   elseif DISTRIBUTED_LOAD == 1
       %% BENCHMARKS - DISTRIBUTED LOAD
        % load distribution
        Nx = 100;
        Ny = 100;
        [xdummy,ydummy]=meshgrid(linspace(-1,1,Nx),linspace(-1,1,Ny));
        s=0.2;
        m=0.6;
        fx=1./s/sqrt(2*pi).*exp(-0.5*((xdummy-m)./s).^2); 
        a = chord;
        b = span;
        xx = a*0.5*(-xdummy+1);
        yy = b*0.5*(ydummy+1);
        
        NACA = 0012;
        [TAU,EPSMAX,PTMAX] = indata(NACA);
        Npanels = 50;
        [Nodtot,X,Y] = setup(Npanels,TAU,EPSMAX,NACA,PTMAX);
        xFoil = (X(2:Npanels/2)-0.5)./0.5; %xdummy 
%         thickFoil = -Y(2:Npanels/2)*a/10; %ydummy (h/a<0.012)
        thickFoil = -Y(2:Npanels/2)*a; %ydummy (h/a<0.12)
        
        modfunc = (xFoil).^2;
        modfunc = modfunc./max(abs(modfunc));
        tx_foil_mod = thickFoil+ 1/100.*modfunc;
        
% %         figure
% %         plot(X,Y,'ro');
% %         hold on;grid on;
% %         plot(xFoil, thickFoil,'ks');
% %         plot(xFoil, tx_foil_mod,'b^');
% %         axis equal;
       
        thickdummy = xdummy.*0;
        for ii = 1:size(xdummy,2)
            thickdummy(ii,:) = interp1(-xFoil, tx_foil_mod, xdummy(ii,:),'linear','extrap');
        end
        
% %         figure
% %         surf(xdummy,ydummy,thickdummy);
       
        % thickness constant
        h0 = 0.001*chord;
        tx_modified = xx.*0 + h0;
        
        for ii = 1:size(xx,1)
            for jj = 1:size(xx,2)
               tx_modified(ii,jj) = 2*thickdummy(ii,jj);
            end
        end
        
        prepareLoad_csv(xx,yy,fx,'load_fx.csv',1); %prepare .csv file for ANSYS
        prepareLoad_csv(xx,yy,tx_modified,'sections_tx.csv',2); %prepare .csv file for ANSYS
   end
    
end

% Triangle center
xm=(1/3)*(x(IEN(1,:))+x(IEN(2,:))+x(IEN(3,:)));   % x barycentric coordinate
ym=(1/3)*(y(IEN(1,:))+y(IEN(2,:))+y(IEN(3,:)));   % y barycentric coordinate

% COMMENT: Use interpolation for this setup 
%   1. interp2 for rectangular domain
%   2. interpolation with scattered data for non-rectangular planforms
%
% COMMENT: If this is implemented in C matlap interpolation would not be
% available (discuss)

% % fxx=interp2(xx,yy,fx,xm,ym, 'spline');    
% % txx=interp2(xx,yy,tx_modified,xm,ym,'spline');%is this OK?

param = 2.55;
if importFromFile.toggle == 1 || DISTRIBUTED_LOAD == 1
%     param = 20.55;
    param = 5.55;
end

fxx=shepard_interp_2d(length(xx(:)),xx(:),yy(:),fx(:), param, length(xm(:)), xm, ym);
txx=shepard_interp_2d(length(xx(:)),xx(:),yy(:),tx_modified(:), param, length(xm(:)), xm, ym);
        
if debugOn
    %% figure 1
    figure;
    colormap(viridis);
    subplot(2,2,[ 1 2]);
    %
    if importFromFile.toggle     
        plot(xc_fem_data(1,:),th_fem_data(1,:),'ro');
    else
        plot(xx(1,:),tx_modified(1,:),'ro');
    end
    hold on;grid on;
    h1=plot(xm,txx,'b^','MarkerSize',2);
    legend("original", "interpolated")
    xlabel('x [m]');ylabel('y [m]');
    title('plate thickness spanwise (from mat-file)');
    %
    subplot(2,2,3);
    hold on;grid on;
    surf(xx,yy,fx);
    plot3(xm,ym,fxx,'b^', 'MarkerSize',2);
    title('bem data for fx');
    xlabel('x-axis');
    ylabel('y-axis');
    zlabel('P [Pa]');view([-25 25])
    %
    subplot(2,2,4);hold on;grid on;
    surf(xx,yy,tx_modified);
    plot3(xm,ym,txx,'b^', 'MarkerSize',2);
    title('bem data for tx');
    xlabel('x-axis');
    ylabel('y-axis');view([-25 25])

    %% figure 2
    FntSz = 16;
    figure;
    subplot(3,3,[1 2 3]);
    grid minor;hold on;
    xlabel('x (m)','Interpreter','Latex');
    ylabel('h (m)','Interpreter','Latex');
    h1=plot(xx(1,:),tx_modified(1,:),'k.-');
    set(gca,'FontSize',FntSz);
    set(gca,'TickLabelInterpreter','latex');
    %
    subplot(3,3,[ 4 5 6 7 8 9]);
    grid minor;hold on;colormap(viridis);
    xlabel('x (m)','Interpreter','Latex');
    ylabel('y (m)','Interpreter','Latex');
    zlabel('h (m)','Interpreter','Latex');
    %
    surf(xx,yy,tx_modified);
    set(gca,'FontSize',FntSz);
    set(gca,'TickLabelInterpreter','latex');
    view([-15 45])

end


end