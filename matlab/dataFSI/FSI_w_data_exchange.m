%% prepare data for bem (half wing)

N1=size(DCoefpres,2)+1;%50;%BEM-mesh
N2=size(DCoefpres,1)+1;%25;%BEM-mesh
[xx,yy]=meshgrid(linspace(0,chord,N1),linspace(-span/2,span/2,N2));
%==========================================================================
for ii = 1:length(t)
    
    w = solution.w(:,ii);
    bx = solution.bx(:,ii);
    by = solution.by(:,ii);
    
    Fw = scatteredInterpolant(pp(1,:)',pp(2,:)',w);
    Fbx = scatteredInterpolant(pp(1,:)',pp(2,:)',bx);
    Fby = scatteredInterpolant(pp(1,:)',pp(2,:)',by);

    wInterp(:,:,ii) = Fw(xx,yy);
    bxInterp(:,:,ii) = Fbx(xx,yy);
    byInterp(:,:,ii) = Fby(xx,yy);
    
    for jj = 1:N2/2
%     for symmetry
        wInterp(N2-jj+1,:,ii) = Fw(xx(jj,:),yy(jj,:));
        bxInterp(N2-jj+1,:,ii) = Fbx(xx(jj,:),yy(jj,:));
        byInterp(N2-jj+1,:,ii) = Fby(xx(jj,:),yy(jj,:));
    end

    %==========================================================================
end


dtest = 200;
figure;
hold on;grid on;
hPlate2=surf(xx,yy,wInterp(:,:,dtest));

plot3(xx(1,:),yy(1,:),wInterp(1,:,dtest),'ks','MarkerFaceColor','r');
plot3(xx(N2,:),yy(N2,:),wInterp(1,:,dtest),'bs','MarkerFaceColor','r');
xlabel('x-axis (chord-wise)');
ylabel('y-axis (span-wise)');
title('w')
view([-25 25]);
% axis equal;
%     [fxx,txx,xx,yy]=Nonunif(x,y,IEN,pp,ee,tt, chord, span, 1, importFromFile, fluid_dens, U, d);
figure
hold on;grid on;
hPlate2=surf(xx,yy,bxInterp(:,:,dtest));
plot3(xx(1,:),yy(1,:),bxInterp(1,:,dtest),'ks','MarkerFaceColor','r');
xlabel('x-axis (chord-wise)');
ylabel('y-axis (span-wise)');
title('bx')
%
figure
hold on;grid on;
hPlate2=surf(xx,yy,byInterp(:,:,dtest));
plot3(xx(1,:),yy(1,:),byInterp(1,:,dtest),'ks','MarkerFaceColor','r');
xlabel('x-axis (chord-wise)');
ylabel('y-axis (span-wise)');
title('by')


% save femDATA_h05_r wInterp bxInterp byInterp xx yy
% save femDATA_h10_r wInterp bxInterp byInterp xx yy
% save femDATA_h15_r wInterp bxInterp byInterp xx yy
save femDATA_h182_r wInterp bxInterp byInterp xx yy

(inData.a3 + max(max(max(abs(wInterp)))))/inData.a3
