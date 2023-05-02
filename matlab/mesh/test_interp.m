% test shepards interpolation method
clear all;close all;clc;
addpath('mesh/heathcote');

file1995 = 'bemDATA_h182_r';
load(file1995);%xc_fem_data, yc_fem_data, DCoefpres
xx=xc_fem_data;
yy=yc_fem_data;
d = 5;
fx = DCoefpres(:,:,d);
tx= th_fem_data;

load('mesh_h1_half');
pp=p;
tt=t;
ee=e;
% Nodal Coordinates
x=pp(1,:);
y=pp(2,:);
% IEN array
IEN=tt(1:3,:);

xm=(1/3)*(x(IEN(1,:))+x(IEN(2,:))+x(IEN(3,:)));   % x barycentric coordinate
ym=(1/3)*(y(IEN(1,:))+y(IEN(2,:))+y(IEN(3,:)));    % y barycentric coordinate

fxx=interp2(xx,yy,fx,xm,ym, 'spline');    
txx=interp2(xx,yy,tx,xm,ym, 'spline');%is this OK?

% fxx2=shepard_interp_2d(length(xm(:)),xm(:),ym(:),fxx(:), 2.22, length(xm(:)), xm, ym);
fxx2=shepard_interp_2d(length(xx(:)),xx(:),yy(:),fx(:), 5.22, length(xm(:)), xm, ym);
txx2=shepard_interp_2d(length(xx(:)),xx(:),yy(:),tx(:), 5.22, length(xm(:)), xm, ym);
%     
% error('we')
% FORCING
figure(1);
subplot(1,2,1);
hold on;grid on;
plot3(xx,yy,fx,'bo');
title('scattered data');
xlabel('x-axis');
ylabel('y-axis');view([-25 25])

subplot(1,2,2);
hold on;grid on;
h1=plot3(xm,ym,fxx,'rs');
h2=plot3(xx,yy,fx,'bo');
h3=plot3(xm,ym,fxx2,'k*');
legend('interp2','scattered data','shepard');
xlabel('x-axis');
xlabel('y-axis');
xlabel('z-axis');
view([-25 25]);
title('comparison')

        
% THICKNESS
figure(2);
subplot(1,2,1);
hold on;grid on;
plot3(xx,yy,tx,'bo');
title('scattered data');
xlabel('x-axis');
ylabel('y-axis');view([-25 25])

subplot(1,2,2);
hold on;grid on;
h1=plot3(xm,ym,txx,'rs');
h2=plot3(xx,yy,tx,'bo');
h3=plot3(xm,ym,txx2,'k*');
legend('interp2','scattered data','shepard');
xlabel('x-axis');
xlabel('y-axis');
xlabel('z-axis');
view([-25 25]);
title('comparison')

