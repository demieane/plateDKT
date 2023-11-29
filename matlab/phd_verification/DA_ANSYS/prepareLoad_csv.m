%==========================================================================
% Preparation of .csv file for External Data ANSYS
%==========================================================================
close all;clear all;clc;
ANSYS = (-1);%convention, positive load toward the negative Z-axis
Nx = 70; 
Ny = 25;
[x,y]=meshgrid(linspace(-1,1,Nx),linspace(-1,1,Ny));
s=0.2;
m=0.6;
fx=ANSYS*1./s/sqrt(2*pi).*exp(-0.5*((x-m)./s).^2); 
a = 10;
b = 10;
x = a*0.5*(-x+1);
y = b*0.5*(y+1);
max(max(abs(fx)))
% fx=1200*fx/max(max(fx));
% % [x,y]=meshgrid(linspace(0,1,N),linspace(0,1,N));
% % s=0.02;
% % m=0.1;
% % fx=1./s/sqrt(2*pi).*exp(-0.5*((x-m)./s).^2); 
% % fx=1200*fx/max(max(fx));

%%
FntSz=16;
figure;hold on;
grid on;
xlabel('$$x (m)$$', 'Interpreter','latex','FontSize',FntSz);
ylabel('$$y (m)$$', 'Interpreter','latex','FontSize',FntSz);
zlabel('$$p_{o} (Pa)$$', 'Interpreter','latex','FontSize',FntSz);
surf(x,y,fx);
view([25 15]);
set(gca,'FontSize',FntSz);
c = colorbar;
c.Location='northoutside';
c.TickLabelInterpreter='latex';
set(gca,'TickLabelInterpreter','latex');

%%
z= 0;
A = zeros(Nx*Ny,3);
for ii=1:Ny
    for jj=1:Nx
        z=z+1;
        A(z,1)=x(ii,jj);
        A(z,2)=y(ii,jj);
        A(z,3)=fx(ii,jj);
    end
end
writecell({'x(m)','y(m)','pressure(Pa)'},'load.csv');
writematrix(A,'load.csv', 'WriteMode', 'append');

