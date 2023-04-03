% clear all;
% close all;
% clc;
addpath('mesh/heathcote');
load('mesh_h1_half');
file1995 = 'bemDATA_h182_r';
load(file1995);
%
chord=inData.cRoot;%0.33; %wing chord length
span=inData.span;%1; % wing span length
Uvel =inData.U;%2.52;%m/s
fluid_dens=1025;%kg/m3
%
m=7850;
E=210e9;
v=0.28;
h=1/1000;%0.12*chord;
% Boundary conditions 
%--------------------------------------------------------------------------
% The present features a polygonal domain with 4 support senarios
% controled by the choice of CC
% CC=1 , simply supported along selected edges
% CC=2 , fully clamped along selected edges
%--------------------------------------------------------------------------
CC=2; % boundary condition toggle
% Special Cantileaver Case CC=4
if CC==1
    disp(' SS case along selected edge (w=0)')
elseif CC==2
    disp(' Clamped along selected edge (w=bx=by=0)')
end
% Forcing
% 1- concetrated load, 2- uniform load, 3- distributed load (mapping func)
lll=3;%2; %loading case
importFromFile=struct('toggle',1,'filename',file1995);
%
P_load = 1; %[Pa] %pointing towards the Z-axis
% in ANSYS load pointing in the negative of Z-axis is positive
if lll==1
   Pxy=[5,5];%load position
end
% Generate mesh - ID, IEN, LM 
pp=p;
tt=t;
ee=e;

%% CREATE MODE.bin binary file for passing data 
precision = 'single';

% write to binary for communication with GPU executable
file = fopen('INDATA_FEM.bin', 'wb');
fwrite(file, chord, precision);
fwrite(file, span, precision);
fwrite(file, Uvel, precision);
fwrite(file, fluid_dens, precision);
fwrite(file, m, precision);
fwrite(file, E, precision);
fwrite(file, v, precision);
fwrite(file, h, precision);
fwrite(file, CC, 'int');
fwrite(file, lll,'int');
% fwrite(file, ee, precision);
fwrite(file, size(pp,1),'int');
fwrite(file, size(pp,2),'int');
fwrite(file, [pp(1,:), pp(2,:)], precision);
%
fwrite(file, size(tt,1),'int');
fwrite(file, size(tt,2),'int');
fwrite(file, [tt(1,:), tt(2,:), tt(3,:), tt(4,:)], 'int');
%
fwrite(file, size(ee,1),'int');
fwrite(file, size(ee,2),'int');
fwrite(file, [ee(1,:), ee(2,:), ee(3,:), ee(4,:), ee(5,:),...
    ee(6,:), ee(7,:)], precision);

fclose(file);

system('cp INDATA_FEM.bin ../c/INDATA_FEM.bin')


%% Understanding the local-to-global assembly
clr = viridis(10);
hfig=figure(1);hold on;grid on;
axis equal;
plot(pp(1,:),pp(2,:),'o', 'MarkerSize',2, 'Color',clr(1,:));
xlabel('x', 'interpreter','latex');
ylabel('y', 'interpreter','latex');
% view([25 45]);

for ii = 1:size(IEN,2)
    triangle=IEN(:,ii);
    plot(pp(1,triangle),pp(2,triangle),'s-','MarkerSize',3, 'Color',clr(1,:));
    plot([pp(1,triangle(end)) pp(1,triangle(1))],...
        [pp(2,triangle(end)) pp(2,triangle(1))],'s-','MarkerSize',3, 'Color',clr(1,:));
end

triangle=IEN(:,100);
plot(pp(1,triangle),pp(2,triangle),'o-','Color',clr(4,:), 'LineWidth',2);
plot([pp(1,triangle(end)) pp(1,triangle(1))],...
    [pp(2,triangle(end)) pp(2,triangle(1))],'-','Color',clr(4,:),'LineWidth',2);
set(gca,'FontSize',14);
saveas(hfig,'ch3_mesh.png');
% plot(pp(1,triangle),pp(2,triangle),'ro');