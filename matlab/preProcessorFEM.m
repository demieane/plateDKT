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
lll=2;%2; %loading case
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
% IEN array
IEN=tt(1:3,:);
% Number of elements
Nelem=length(IEN(1,:));
% Nodal Coordinates
% x=pp(1,:);
% y=pp(2,:);

% ID array
NN=max(max(IEN));
ID=zeros(3,NN);
%
ID(1,:)=1:3:3*NN-2;
ID(2,:)=2:3:3*NN-1;
ID(3,:)=3:3:3*NN;
% LM array degrees of freedom per element edof*elements
LM=zeros(9,Nelem);
for k=1:Nelem
   for i=1:3
       for j=1:3
           P=(3)*(j-1)+i; %the 9 dofs per triangle
           LM(P,k)=ID(i,IEN(j,k));
       end
   end
end

%% ASSUMING THAT THE PLATE SHAPE HAS 4-EDGES
Bound1=find(e(5,:)==1);
Bound2=find(e(5,:)==2);
Bound3=find(e(5,:)==3);
Bound4=find(e(5,:)==4);
%************************THIS IS THE ACTIVE BOUNDARY CONDITION*************
% COMMENT: The numbering is offered by the pdeModeler
Bnodes = [Bound4(1), Bound3];
%**************************************************************************
%
BBnodes = Bnodes.*0;
for i=1:length(Bnodes)
    BBnodes(i)=e(1,Bnodes(i));
end
BBnodes=sort(BBnodes);  
Bdofs1=ID(1,BBnodes);%w
Bdofs2=ID(2,BBnodes);%bx
Bdofs3=ID(3,BBnodes);%by
%
if CC == 1
    Bdofs_SS=sort(Bdofs1); %Dofs for w and theta
    Bdofs=Bdofs_SS;
elseif CC == 2
    Bdofs_CL=sort([Bdofs1 Bdofs2 Bdofs3]); %Dofs for w and theta
    Bdofs=Bdofs_CL;
end
%


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

fwrite(file, length(BBnodes), 'int');
fwrite(file, BBnodes, 'int');
fwrite(file, length(Bdofs), 'int');
fwrite(file, Bdofs, 'int');

fwrite(file, P_load, precision);
fclose(file);

system('cp INDATA_FEM.bin ../c/INDATA_FEM.bin')


%% read solution from binary file
fileID = fopen('../c/OUTDATA_FEM.bin','rb')
GEN_fromC = fread(fileID,1,'int')
rowsUsol = fread(fileID,1,'int')
colsUsol = fread(fileID,1,'int')

for i = 1:rowsUsol
    for j = 1:colsUsol
        Usol(i,j)=fread(fileID,1,'single');
    end
end

%%  VISUAL COMPARISON 
% load solMatlab
% u=U(1:GEN); % the vector of nodal unknowns (w1;bx1;by1;....wN;bxN;byN)
% %
% w=u(1:3:end);   % vertical displacement
% bx=u(2:3:end);  % rotation x
% by=u(3:3:end);  % rotation y
% 
% figure(1);
% hold on;grid on;
% subplot(1,3,[1 2]);
% plot3(pp(1,BBnodes),pp(2,BBnodes),w(BBnodes),'ks','MarkerSize',3);
% hh=pdeplot(pp,ee,tt,'XYData',w,"ZData",w,'colormap','jet');
% colorbar;shading interp;view([25 25]);%axis equal;
% zlim([-2.5*max(max(abs(w))) 2.5*max(max(abs(w)))])
% xlabel('x-axis');ylabel('y-axis');zlabel('w [m]');
% title('matlab')
% %     title('w displacement','FontWeight','normal');
% subplot(1,3,3);hold on;grid on;
% pdeplot(pp,ee,tt,'XYData',w,'colormap','jet','contour','on');
% colorbar;shading interp;
% xlabel('x-axis');ylabel('y-axis');
% title('(contour)','FontWeight','normal');

u_fromC=Usol(1:GEN_fromC); % the vector of nodal unknowns (w1;bx1;by1;....wN;bxN;byN)
%
w_fromC=u_fromC(1:3:end);   % vertical displacement
bx_fromC=u_fromC(2:3:end);  % rotation x
by_fromC=u_fromC(3:3:end);  % rotation y


figure(2);
hold on;grid on;
subplot(1,3,[1 2]);
plot3(pp(1,BBnodes),pp(2,BBnodes),w(BBnodes),'ks','MarkerSize',3);
hh=pdeplot(pp,ee,tt,'XYData',w_fromC,"ZData",w_fromC,'colormap','jet');
colorbar;shading interp;view([25 25]);%axis equal;
zlim([-2.5*max(max(abs(w_fromC))) 2.5*max(max(abs(w_fromC)))])
xlabel('x-axis');ylabel('y-axis');zlabel('w [m]');
title('sol from c')
%     title('w displacement','FontWeight','normal');
subplot(1,3,3);hold on;grid on;
pdeplot(pp,ee,tt,'XYData',w_fromC,'colormap','jet','contour','on');
colorbar;shading interp;
xlabel('x-axis');ylabel('y-axis');
title('(contour)','FontWeight','normal');

error('er')
%% connectivity in c



close all;

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

hbcs=plot(p(1,BBnodes),p(2,BBnodes),'rs','MarkerSize',5);
% title('boundary condition affected element nodes','FontWeight','normal');
legend([hbcs],'BCs')

saveas(hfig,'ch3_meshBCs.png');
