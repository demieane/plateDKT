%==========================================================================
% 
%           FSI INTERFACE
%
%==========================================================================
clear all;
close all;
clc;

debugOn = 1; 

MODAL_ANALYSIS = 1;
DYNAMIC_ANALYSIS = 0;
addpath('mesh/hydrofoil');%naca functions
% addpath('mesh/heathcote');
addpath('phd_verification');

addpath('phd_verification/rect_constant_thick');

%load('eig_rect_1');%334
load('eig_rect_2');%1336
%load('eig_rect_3');%5344
%
chord=10; %wing chord length
span=10;  % wing span length
Uvel = 1;%inData.U;%2.52;%m/s
fluid_dens=1025;%kg/m3

%% STEEL
m=7850;
E=210*10^9;%200e9;
v=0.3;
h=0.01;%0.12*chord;

%% Boundary conditions 
%--------------------------------------------------------------------------
% The present features a polygonal domain with 4 support senarios
% controled by the choice of CC
% CC=1 , simply supported along selected edges (Navier Solution static)
% CC=2 , fully clamped along selected edges
%--------------------------------------------------------------------------
CC=1; % boundary condition toggle
% Special Cantileaver Case CC=4
if CC==1
    disp(' SS case along selected edge (w=0)')
elseif CC==2
    disp(' Clamped along selected edge (w=bx=by=0)')
end
%
%% Forcing
% 1- concetrated load, 2- uniform load, 3- distributed load (mapping func)
lll=2; %loading case
% % importFromFile=struct('toggle',1,'filename',file1995);
%
if lll==1
    P_load = 100; %[N] %pointing towards the Z-axis
elseif lll==2
    P_load = 100/chord/span; %[Pa]
end
if lll==3
    P_load = [];
end
% in ANSYS load pointing in the negative of Z-axis is positive
if lll==1
   Pxy=[5,5];%load position
%    Pxy=[0.05,-0.15];%load position
end
%% Generate mesh - ID, IEN, LM 
pp=p;
tt=t;
ee=e;
% IEN array
IEN=tt(1:3,:);
% Number of elements
Nelem=length(IEN(1,:));
% Nodal Coordinates
x=pp(1,:);
y=pp(2,:);

%***************************ADDITION (FOR CONCENTRATED LOAD)***************
if lll==1
    tol=1;
    I=find(abs(x-Pxy(1))<tol & abs(y-Pxy(2))<tol);
    rr=sqrt((x(I)-Pxy(1)).^2+(y(I)-Pxy(2)).^2);
    [II,JJ]=min(rr);
    PNODE=I(JJ);
end
%**************************************************************************
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
    
GEN=max(max(LM)); %Total number of nodal unknowns taking into account the 
% connectivity between the triangles
%
%% Essential BCs (enforces on the displacements "w" and slopes "bx", "by")
% THE PLATE HAS 4 EDGES - RECTANGULAR
Bound1=find(e(5,:)==1);
Bound2=find(e(5,:)==2);
Bound3=find(e(5,:)==3);
Bound4=find(e(5,:)==4);
%************************THIS IS THE ACTIVE BOUNDARY CONDITION
% COMMENT: The numbering is offered by the pdeModeler
%Bnodes= [Bound4, Bound1(1)]; %FULL EDGE (x=0)
%Bnodes = [Bound4(1), Bound3]; %(x=a)
% Bnodes = [Bound4];
%Bnodes=Bound3; %for distributed load from function ANSYS 
Bnodes=[Bound1 Bound2 Bound3 Bound4];
%*************************************************************
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
%     Bdofs_SS=sort(Bdofs1); %Dofs for w and theta
    Bdofs_SS=sort([Bdofs1]); %Dofs for w and theta
    Bdofs=Bdofs_SS;
elseif CC == 2
    Bdofs_CL=sort([Bdofs1 Bdofs2 Bdofs3]); %Dofs for w and theta
    Bdofs=Bdofs_CL;
end
%
if debugOn
    figure;hold on;
    plot(p(1,:),p(2,:),'ro','MarkerSize',3);
    plot(p(1,BBnodes),p(2,BBnodes),'ks','MarkerSize',3);
    title('boundary condition affected element nodes','FontWeight','normal');
    xlabel('x-axis');
    ylabel('y-axis');
    if lll == 1
        plot(Pxy(1),Pxy(2),'bs');
    end
end

%% CREATE MODE.bin binary file for passing data 
modeFem = 1; %double:1, single(or mixed):2
% write to binary for communication with GPU executable
% file1 = fopen('MODE_FEM.bin', 'wb');
% fwrite(file1, modeFem, 'int');
% fclose(file1);
% system('cp MODE_FEM.bin ../c/MODE_FEM.bin');
%
if modeFem == 1
    precision = 'double';
    fileName = 'INDATA_FEM_double.bin';
elseif modeFem == 2
    precision = 'single';
    fileName = 'INDATA_FEM_single.bin';
end

% write to binary for communication with GPU executable
file = fopen(fileName, 'wb');
fwrite(file, modeFem, 'int');
fwrite(file, chord, precision);
fwrite(file, span, precision);
fwrite(file, Uvel, precision);
fwrite(file, fluid_dens, precision);
fwrite(file, m, precision);
fwrite(file, E, precision);
fwrite(file, v, precision);
fwrite(file, h, precision);%average thickness
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

if lll == 2 %uniform load
    fwrite(file, P_load, precision);
end
if lll==1%point load
    fwrite(file, P_load, precision);
    fwrite(file, Pxy(1), precision);
    fwrite(file, Pxy(2), precision);
    fwrite(file, PNODE, 'int');
end
if lll==3%distributed load & distributed thickness
    file1995 = 'bemDATA_h182_r';
    load(file1995); %xc_fem_data, yc_fem_data, DCoefpres
    xcp=xc_fem_data;
    ycp=yc_fem_data;
    tinstance = 100;
    fcp = DCoefpres(:,:,tinstance)*(0.5*fluid_dens*Uvel^2);
    tcp= th_fem_data.*0 + h;
%     tcp= th_fem_data;
    fwrite(file, length(xcp(:)),'int');
    % x-coords
    for ii = 1:length(xcp(:))
        fwrite(file, xcp(ii),precision);
    end
    % y-coords
    for ii = 1:length(ycp(:))
        fwrite(file, ycp(ii),precision);
    end
    % load at collocation point (to be updated)
    for ii = 1:length(fcp(:))
        fwrite(file, fcp(ii),precision);
    end
    % thickness at collocation point (to be updated) - TODO: increase the
    % values for the solution
    for ii = 1:length(tcp(:))
        fwrite(file, tcp(ii),precision);
    end
end
if DYNAMIC_ANALYSIS == 1
    fwrite(file, inData.omega3, precision);
    fwrite(file, inData.dt, precision);
    fwrite(file, inData.Nper, precision);
end
fclose(file);

if modeFem == 1
    system('cp INDATA_FEM_double.bin ../c/INDATA_FEM_double.bin');
elseif modeFem == 2
    system('cp INDATA_FEM_single.bin ../c/INDATA_FEM_single.bin');
end

% error('er')

%% RUN THE CODE

tstart = tic;
if modeFem == 1
    system('../c/./mainDKT_CPP2');
elseif modeFem == 2
    system('../c/./mainDKT_CPP');
%     system('../c/./mainDKT');
end
telapsed = toc(tstart)

%% Read solution from binary file

if modeFem == 1
    fileID = fopen('../c/OUTDATA_FEM_double.bin','rb')
elseif modeFem == 2
    fileID = fopen('../c/OUTDATA_FEM_single.bin','rb')
end
% fileID = fopen('../c/OUTDATA_FEM.bin','rb')
GEN_fromC = fread(fileID,1,'int')
rowsUsol = fread(fileID,1,'int')
colsUsol = fread(fileID,1,'int')

sizeM = rowsUsol/2;

for i = 1:rowsUsol
    for j = 1:colsUsol
        Usol(i,j)=fread(fileID,1,precision);
    end
end

u_fromC=Usol(1:GEN_fromC,:); % the vector of nodal unknowns (w1;bx1;by1;....wN;bxN;byN)

w_fromC=u_fromC(1:3:end,:);   % vertical displacement
bx_fromC=u_fromC(2:3:end,:);  % rotation x
by_fromC=u_fromC(3:3:end,:);  % rotation y

% % load('FEM_sol_h182_r_h2');
% % 
% % max(abs(Usol(1:sizeM,1:end-1) - u(1:sizeM,2:end) ))
% % 
% % 
% % 
% % uu_dot = u(1:GEN,:);

% uu_dot99 = Usol(1:GEN_fromC,:); 
% 
% uu99 = Usol(sizeM+1:sizeM+GEN_fromC,:);
% 
% w99=uu99(1:3:end,:);   % vertical displacement


%% DEBUG
SOLVE_FROM_C_MATRICES = 0;

if SOLVE_FROM_C_MATRICES == 1
    fileID = fopen('../c/OUTDATA_FEM_Kglob_Mglob_BCs.bin','rb')
    rowsUsol = fread(fileID,1,'int')
    colsUsol = fread(fileID,1,'int')

    for i = 1:rowsUsol
        for j = 1:colsUsol
            Kglob_aug2(i,j)=fread(fileID,1,precision);
        end
    end

    for i = 1:rowsUsol
        for j = 1:colsUsol
            Mglob_aug2(i,j)=fread(fileID,1,precision);
        end
    end

    for i = 1:rowsUsol
        for j = 1:1
            Fglob_aug2(i,j)=fread(fileID,1,precision);
        end
    end
      
    if MODAL_ANALYSIS == 1
        [XX,lamM,flag]=eigs(sparse(Kglob_aug2),sparse(Mglob_aug2),5,'smallestabs');
        cc=sort(diag(lamM));

        freq=sqrt(sort(diag(lamM),'ascend'))./(2*pi);

        freq'
    end
    
    % load singleMatrix
    U=Kglob_aug2\Fglob_aug2;%SOLVE SPARSE SYSTEM OF EQUATIONS
    u_fromC=U(1:GEN_fromC); % the vector of nodal unknowns (w1;bx1;by1;....wN;bxN;byN)

% % % %     u_fromC=Usol(1:GEN_fromC,:); % the vector of nodal unknowns (w1;bx1;by1;....wN;bxN;byN)

    w_fromC=u_fromC(1:3:end,:);   % vertical displacement
    bx_fromC=u_fromC(2:3:end,:);  % rotation x
    by_fromC=u_fromC(3:3:end,:);  % rotation y
    
    % save solFromC w_fromC

    % save singleMatrix Kglob_aug Fglob_aug Mglob_aug GEN
    % Kglob_aug_double = Kglob_aug;
    % Mglob_aug_double = Mglob_aug;
    % Fglob_aug_double = Fglob_aug;
    % save doubleMatrix Kglob_aug_double Fglob_aug_double Mglob_aug_double
end

%%  VISUAL COMPARISON 

GEN = size(pp,2)*3;
load solMatlab
u=Ustatic(1:GEN); % the vector of nodal unknowns (w1;bx1;by1;....wN;bxN;byN)
%
w=u(1:3:end);   % vertical displacement
bx=u(2:3:end);  % rotation x
by=u(3:3:end);  % rotation y

figure(1);
subplot(2,2,1);hold on;grid on;
% subplot(1,3,[1 2]);
plot3(pp(1,BBnodes),pp(2,BBnodes),w(BBnodes),'ks','MarkerSize',3);
hh=pdeplot(pp,ee,tt,'XYData',w,"ZData",w,'colormap',viridis);
colorbar;shading interp;view([25 25]);%axis equal;
zlim([-2.5*max(max(abs(w))) 2.5*max(max(abs(w)))])
xlabel('x-axis');ylabel('y-axis');zlabel('w [m]');
title('Matlab solver (double precision) [LEFT]')
%     title('w displacement','FontWeight','normal');
% subplot(1,3,3);
subplot(2,2,3);
hold on;grid on;
pdeplot(pp,ee,tt,'XYData',w,'colormap',viridis,'contour','on');
colorbar;shading interp;
xlabel('x-axis');ylabel('y-axis');
title('(contour)','FontWeight','normal');
% 
% u_fromC=Usol(1:GEN_fromC); % the vector of nodal unknowns (w1;bx1;by1;....wN;bxN;byN)
% %
% w_fromC=u_fromC(1:3:end);   % vertical displacement
% bx_fromC=u_fromC(2:3:end);  % rotation x
% by_fromC=u_fromC(3:3:end);  % rotation y

% figure(2);
subplot(2,2,2);hold on;grid on;
plot3(pp(1,BBnodes),pp(2,BBnodes),w(BBnodes),'ks','MarkerSize',3);
hh=pdeplot(pp,ee,tt,'XYData',w_fromC,"ZData",w_fromC,'colormap',viridis);
colorbar;shading interp;view([25 25]);%axis equal;
zlim([-2.5*max(max(abs(w_fromC))) 2.5*max(max(abs(w_fromC)))])
xlabel('x-axis');ylabel('y-axis');zlabel('w [m]');
title(['C solver (', num2str(precision),' precision) [RIGHT]'])
subplot(2,2,4);
hold on;grid on;
pdeplot(pp,ee,tt,'XYData',w_fromC,'colormap',viridis,'contour','on');
colorbar;shading interp;
xlabel('x-axis');ylabel('y-axis');
title('(contour)','FontWeight','normal');

check = [max(abs(w_fromC)), max(abs(w))]


error('er')
%%

load singleMatrix
load doubleMatrix


max(abs(Fglob_aug - Fglob_aug_double))
max(max(Kglob_aug - Kglob_aug_double))
max(max(Mglob_aug - Mglob_aug_double))


max(Kglob_aug(1:GEN,1:GEN) - Kglob_aug_double(1:GEN,1:GEN))

GEN = 759;
max(max(Kglob_aug(GEN:end,GEN:end) - Kglob_aug_double(GEN:end,GEN:end)))

max(max(Kglob_aug(GEN+1:end,1:GEN) - Kglob_aug_double(GEN+1:end,1:GEN)))

max(max(Kglob_aug(1:GEN,GEN+1:end) - Kglob_aug_double(1:GEN,GEN+1:end)))

max(max(Kglob_aug(GEN+1:end,GEN+1:end) - Kglob_aug_double(GEN+1:end,GEN+1:end)))






