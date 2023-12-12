function [] = animatorFunc(solution,Nt,ddt,pp,ee,tt,BBnodes, Name,inData)
% if animationSwitch = true, it produces a visualization of the 
% flapping foil movemennt. Giving information about the forcing
% CL and the pressure-type Kutta evaluation at each time step.

t = 0; % Initializing // corresponding to d=1
chord=inData.cRoot;%0.33;
span=inData.span;%1.0;

% animation settings      

hfig=figure;
set(hfig,'Position',[1 41 1366 652],'PaperPositionMode','auto');
hax = subplot(1,1,1,'parent',hfig,'XGrid','on','YGrid','on'); %It returns the objects identity
hold on;grid on;
xlabel(hax,'x');
ylabel(hax,'y');
zlabel(hax,'w');
htit = title(hax,['t=',num2str(t, '%5.3f sec')]);


writerObj = VideoWriter([Name,'.avi'],'Motion JPEG AVI');   
open(writerObj);

d=1;


w=solution.w(:,d);
F = scatteredInterpolant(pp(1,:)',pp(2,:)',w);
N1=50;
N2=25;
% [xx,yy]=meshgrid(linspace(0,chord,N1),linspace(-span/2,span/2,N2));
% [xx,yy]=meshgrid(linspace(0,chord,N1),linspace(-span/2,span/2,N2));
[xx,yy]=meshgrid(linspace(0,chord,N1),linspace(-span/2,0,N2));
% figure
% surf(xx,yy,ww);

set(hax,'ZLim',[-max(max(abs(solution.w))) max(max(abs(solution.w)))]);
% bx=solution.bx(:,d);
% by=solution.by(:,d);


hBCs = plot3(pp(1,BBnodes),pp(2,BBnodes),w(BBnodes),'ks','MarkerSize',3);
hPlate2=surf(xx,yy,F(xx,yy));
colormap(viridis);
% colorbar;shading interp;
% hPlate=plot3(pp(1,:),pp(2,:),w,'o','MarkerSize',3);
view([-25 25]);
% hPlate = pdeplot(pp,ee,tt,'XYData',w,"ZData",w,'colormap','jet');
% hPlate(1,1)

for d = 2:Nt%N+1
    %% Since I want to reproduce the final results here, I need to make
    % sure that in each time step I reveal the right amount of information
    % based on what was known until then.
    % Check out(!): ExistingUntilThatTime
    pause(0.05);
    
    t
    ddt
    
    t=(d-1)*ddt;
    
    w=solution.w(:,d);
    
    F2 = scatteredInterpolant(pp(1,:)',pp(2,:)',w);

% bx=solution.bx(:,d);
% by=solution.by(:,d);

    set(hBCs,'ZData',w(BBnodes));
   
    set(hPlate2,'Zdata',F2(xx,yy));
%     set(hPlate,'ZData',w);
%     hPlate = pdeplot(pp,ee,tt,'XYData',w,"ZData",w,'colormap','jet');

%     set(hPlate(1,1),'ZData',w);
    view([-25 25]);
    htit = title(hax,['t/T=',num2str(t/inData.T3, '%5.3f')]);


    Fig=getframe(hfig);
    writeVideo(writerObj,Fig);
end

close(writerObj);
end

