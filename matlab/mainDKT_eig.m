% =========================================================================
% Discrete Kirchoff Triangle (PLATE ELEMENT)
% References : Batoz et al (1989), SYDENSTRICKER et al (1995)
%
% Small strains, Isotropic material, Homogeneous properties
%
% Written by: A. Karperaki
% Modified by: D. Anevlavi
%
% Search for the keyword "COMMENT" within the code.
%
% Discuss the limitations of the present methodology. w compared to
% the plate dimensions? when do i know it is an acceptable result?
%
% Which quantities I can predict with the weak coupling?
% What is the correct dimensionalization for the distributed load?
% Many questions have been answered in the 2-D version of the code!!

% =========================================================================

% COMMENT: Frequency non-dimensionalization is different in each work
% make sure we compared against the correct quantity!!!!!! 
% Typo in Pachenari Table 8, D: based on smallest thickness

clear all;
close all;
clc;

addpath('funcFEM');
addpath('dataFSI');

MODAL_ANALYSIS = 1;
DYNAMIC_ANALYSIS = 0;

tstart = tic;   

FntSz = 18;

%% MATERIAL PROPERTIES

% in docx -> at the limit of moderately thick plates 
% % % m=7850;%kg/m2 [mass distribution]
% % % E=210*10^9;%Pa [Young modulus]
% % % v=0.3;% [Poisson ratio]
% % % h=0.1;%m [thickness] 

% Comparison with Leissa 
% m=7850;%kg/m2 [mass distribution]
% E=210*10^6;%Pa [Young modulus]
% v=0.3;% [Poisson ratio]
% h=0.01;%m [thickness] 

% Comparison with Pachenary 2014 & Thesis 
m=7850;%kg/m2 [mass distribution]
E=210*10^9;%Pa [Young modulus]
v=0.3;% [Poisson ratio]
h=1;%0.01;%m [thickness] 
% h=1; % Shufrin

% G= E/(2*(1+v));
% Dplate=E*h^3/(12*(1-v^2));
% shear correction factor
% sc=5/6; 
% phi=E/(sc*G)*(h/10)^2;

%% MESH AND BOUNDARY CONDITIONS
%==========================================================================
% INFORMATION ABOUT MESH GENERATION USING pdeModeler
%
% Options -> Axis Limits (set to auto to view your sketch)
% Boundary -> Show Edge Labels
% Mesh -> Refine Mesh, -> Export Mesh
% Close pdeModeler or leave it open for easy mesh refinement
%--------------------------------------------------------------------------
% p (points, the mesh nodes) = [x1,x2,x3,x4...;y1,y2,y3]; nodal coordinates
%
% e (edges): 
%   e(1,k) is the index of the first point in mesh edge k
%   e(2,k) is the index of the second point in mesh edge k 
%   e(3,k) is the parameter value of the first point in edge k
%   e(4,k) is the parameter value at the second poiny of edge k
%   e(5,k) is the ID of the geometric edge containing the mesh edge
%   e(6,k),e(7,k) subdomain numbers
%
% t (triangles)
%==========================================================================

addpath('mesh');
% addpath('mesh/heathcote');
addpath('phd_verification/rect_constant_thick');
% load('mesh_h1_half');

chord = 10;
span = 10;
%load('eig_rect_1');%334
load('eig_rect_2');%1336
%load('eig_rect_3');%5344
%load('eig_rect_4');%21376
% load('eig_rect_5');%85504

%  8.9632   18.2788   18.2788   26.9483   32.7641
%  8.9699   18.3481   18.3538   27.1870   33.0950

%    36.0086
%    73.6559
%    73.6788
%   109.1383
%   132.8553
%   134.7860


%Create Rectangular plate
% pderect([0 10 0 10],'C1')

% Generate mesh - ID, IEN, LM 
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

%% BOUNDARY CONDITIONS
%--------------------------------------------------------------------------
% The present features a polygonal domain with support senarios
% controled by the choice of CC
% CC=1 , simply supported along selected edges
% CC=2 , fully clamped along selected edges
%--------------------------------------------------------------------------
%CC=1; %boundary condition toggle
CC=2;
% CC=4 Special cantilever case from Aggelinas other files
if CC==1
    disp(' SS case along selected edge (w=0)')
elseif CC==2
    disp(' Clamped along selected edge (w=bx=by=0)')
end

% Dependent on the mesh - PDE modeler for edge numbering 
% Essential BCs (enforces on the displacements "w" and slopes "bx", "by")   
Bound1=find(e(5,:)==1);
Bound2=find(e(5,:)==2);
Bound3=find(e(5,:)==3);
Bound4=find(e(5,:)==4);

%************************THIS IS THE ACTIVE BOUNDARY CONDITION*************
% COMMENT: The numbering is offered by the pdeModeler
%Bnodes= [Bound4, Bound1(1)]; %FULL EDGE
%Bnodes = Bound2; %for distributed load from function ANSYS
Bnodes = [Bound1 Bound2 Bound3 Bound4];
%Bnodes = [Bound1, Bound2(1)];
%**************************************************************************

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
    Bdofs_SS=sort([Bdofs1]); %Dofs for w and theta
    Bdofs=Bdofs_SS;
elseif CC == 2
    Bdofs_CL=sort([Bdofs1 Bdofs2 Bdofs3]); %Dofs for w and theta
    Bdofs=Bdofs_CL;
end

figure(1)
hold on;grid on;
h1=plot(p(1,:),p(2,:),'bo','MarkerSize',3);
h2=plot(p(1,BBnodes),p(2,BBnodes),'ks','MarkerSize',8);
% title('BC-affected element nodes','FontWeight','normal','Interpreter','latex',...
%     'FontSize',FntSz);
%title('Case-1: Mesh-1','FontWeight','normal','Interpreter','latex',...
%    'FontSize',FntSz);
xlabel('x (m)','Interpreter','latex', 'FontSize',FntSz);
ylabel('y (m)','Interpreter','latex', 'FontSize',FntSz);
grid minor;
set(gca,'FontSize',FntSz);
set(gca,'TickLabelInterpreter','latex');
pdeplot(pp,ee,tt);%,'colormap','copper')
legend([h2], 'BC-affected nodes','Interpreter','latex',...
    'FontSize',FntSz,'Location','northeast');


% error('er')

disp(['*=====================================================*']);
disp(['*                  PLATE 2D-FEM                       *']);
disp(['*=====================================================*']);
disp([' MESH SIZE (DKT Triangles 9-dofs)= ',num2str(size(t,2))]);
disp([' thick/width or thick/length < 0.1 (10%): thin plate'])
disp([' thick/width or thick/length > 0.1 (10%): moderately thick plate'])
disp([' > h/chord= ',num2str(h/chord)]);
disp([' > h/span = ',num2str(h/span)]);


%==========================================================================
%                             PROCESSOR
%==========================================================================
%% GAUSS POINTS AND WEIGHTS FOR TRIANGLES
% COMMENT: The code offers full integration. For particular cases of thin 
% plates with clamped boundary conditions, shear locking effects make
% the convergence of the code slow. Reduced integration could be a remedy
% for this issue. However, it is not included in the present version.
% When shear locking appears, the boundary conditions are not satisfied
% for sparse triangular meshes. 

% Number of GP
Ng=3; %TO-DO
[xw]=TriGaussPoints(Ng);

%==========================================================================
% BENDING STIFFNESS MATRIX (3x3) FOR EACH TRIANGLE
CONSTANT_THICK = 0;%1;
LINEAR_THICK = 1;
d=100;
if CONSTANT_THICK == 1
    thick=h*ones(1,Nelem);
    txxBEM = thick;
    loadFEM = [];
elseif LINEAR_THICK == 1
    importFromFile.toggle = 0; % NO-FILE to read
    [loadFEM,txxBEM]=Nonunif(x,y,IEN,p,e,t, chord, span, 1, importFromFile,...
        1025, 1, h, d);
end
% error('er')
[BeSt2]=BendingStiffness2(E,v,txxBEM,h); %[3,3] matrix

xg=xw(:,1);
yg=xw(:,2);
[l23,l31,l12,y12,y31,y23,x12,x31,x23,Area,a4,a5,a6,b4,b5,b6,...
    c4,c5,c6,d4,d5,d6,e4,e5,e6,C4,C5,C6,S4,S5,S6]=TrigElCoefsDKT(x,y,IEN);

% Calculation of shape functions and their derivatives at gauss points on
% the parent element
[ SF,DxsiSF,DetaSF,D2xsiSF,D2xsietaSF,D2etaSF ] = LNShapeFunDST(xg,yg);
[ SFm, DxsiSFm,DetaSFm] = LNShapeFunMassDST(xg,yg); 

%---> ax, ay, bx,by are constant for constant h (independednt of î,ç)
% %%
[GGDST,GGDKT] = matrixG();

GGin=inv(GGDST);
GGin2=inv(GGDKT);

% % %************************************** DISTRIBUTED LOAD VIA MAPPING*******
% % if lll==3
% %     [Fx,~]=Nonunif(x,y,IEN,pp,ee,tt, chord, span, 0, importFromFile,...
% %         fluid_dens, Uvel,h,  d);
% % %     Fx = -Fx; %for ansys
% % %     if debugOn
% % %         hold on;grid on;
% % %         plot3(pp(1,BBnodes),pp(2,BBnodes),pp(2,BBnodes).*0,'ks');
% % %         xlabel('x-axis');
% % %     end
% % end
% % if lll == 2 || lll == 1
% %     Fx = [];
% % end
%**************************************************************************

%==========================================================================
%       PREPARE Kg, Mg
%==========================================================================

for kk=1:Nelem %for each element (iS THIS TRIANGLE 1 IN t?)
  
    % Mass matrix
    [ Hm,HW] = massHmDKT(kk,l23,l31,l12,C4,C5,C6,S4,S5,S6,GGin);
    %% Integration
    kb=zeros(9,9); 
    mloc=zeros(9,9);
    mloc1=zeros(9,9);
    mloc2=zeros(9,9);

    [ Hxx,Hyy ] = rotationMass2(kk,a4,a5,a6,b4,b5,b6,c4,c5,c6,d4,d5,d6,e4,e5,e6 );
    [ Hxx1,Hyy1 ] = rotationMass(kk,l23,l31,l12,C4,C5,C6,S4,S5,S6);
    
    for ii=1:Ng % for each gauss point
    % Stiffness
      [Hx, Hy, Hx_xsi,Hx_eta,Hy_xsi, Hy_eta, Bb ] = ShapeFunDKT2(ii,kk,C4,C5,C6,S4,S5,S6,Area,...
          SF, DxsiSF,DetaSF,y31,x31,y12,x12,l23,l31,l12);
       
       [ Nm, HW,LW,L,HX3,~] = pseudoMassDKT(ii,kk,l23,l31,l12,y12,y31,x12,x31,Area,...
           Hm,GGin,GGin2,xg, yg,IEN,x,y,Hx,Hy,SFm,C4,C5,C6,S4,S5,S6,Hxx,Hyy,HW );

       % we give one more dimension 
       % to the BeSt [3x3] original matrix
       kb=kb+Area(kk)*xw(ii,3)*(Bb'*BeSt2(:,:,kk)*Bb);

       mloc=mloc+m*txxBEM(kk)*Area(kk)*xw(ii,3)*((HW'*(LW'*LW)*HW)+txxBEM(kk)^2/12*(Hx'*Hx)+txxBEM(kk)^2/12*(Hy'*Hy));
    end
    kloc=kb;

    Mg(:,kk)=[mloc(:,1);mloc(:,2);mloc(:,3);mloc(:,4);mloc(:,5);mloc(:,6);mloc(:,7);mloc(:,8);mloc(:,9)];

    Kg(:,kk)=[kloc(:,1);kloc(:,2);kloc(:,3);kloc(:,4);kloc(:,5);kloc(:,6);kloc(:,7);kloc(:,8);kloc(:,9)];

end

%% Global assembly (uses the LM)
% The sparse function accumulate the values that have identical
% subscripts. 
% S(i(k),j(k)) = v(k)

% disp('===================ASSEMBLY==============================')
iii=1:9; %Elnodes Number of element nodes
iii=repmat(iii',1,9);

rr=iii';
Ig=LM(iii(:),:);  %LM ARRAY (Dofs per Element)X(Elements)   -   (Nnodes*Ndof)X(Nel)
Jg=LM(rr(:),:);

Kglob=sparse(Ig(:),Jg(:),Kg(:),GEN,GEN);
% % spy(Kglob)
Mglob=sparse(Ig(:),Jg(:),Mg(:),GEN,GEN);
    
% COMMENT: The boundary conditions are enforced as extra equations
% in the sense of constraints in the present version

kkk=zeros(length(Bdofs),GEN);
mmm=zeros(length(Bdofs),GEN);
Dofs=length(pp)*3;

% Totbound=reshape(BBound',1,4*size(BBound,2));
%***************************ADDITION***************************************
for j=1:length(Bdofs)
    kkk(j,:)=[zeros(1,Bdofs(j)-1) 1 zeros(1,Dofs-Bdofs(j))];
end
%**************************************************************************
Kglob=[Kglob kkk'; kkk zeros(length(Bdofs))];
Mglob=[Mglob mmm'; mmm zeros(length(Bdofs))];

telapsed = toc(tstart);

% =========================================================================
%                           GENERAL COMMENTS
%
% COMMENT: Making the mesh more dense closer to the boundaries with clamped
% edges could solve the problem of fine grids?
% COMMENT: Discuss with Aggelina the bx, by fields in satisfying the BCs
%==========================================================================

[XX,lamM,flag]=eigs(Kglob,Mglob,15,'sm');
cc=sort(diag(lamM));

freq=sqrt(sort(diag(lamM),'ascend'))./(2*pi); %Hz
freq(1:5)

if CONSTANT_THICK == 1 || LINEAR_THICK == 1
%     h=min(txxBEM); %Katsikadelis
    h=max(txxBEM); % Shurfin & Eisenberger (from Pachenari)
    DD=E*h^3/(12*(1-v^2));
    cl=sqrt(m*h/DD);
else
    DD=E*mean(txxBEM)^3/(12*(1-v^2));
    cl=sqrt(m*mean(txxBEM)/DD);
end

NDfreq=cl*(freq*2*pi)*chord^2
XXX=XX(1:GEN,:);

rho =m*min(txxBEM);
DDk=E*min(txxBEM)^3/(12*(1-v^2));
Katsikadelis_freq=(sqrt(rho/DDk)*(freq*2*pi)*chord^2).^2

omega = (freq*2*pi);
d0 = E*max(txxBEM)^3/(12*(1-v^2));
shufrin_freq = omega*span^2*sqrt(m*max(txxBEM)/d0)/pi^2

%**************************************************************************
% Select to view the eigen-fuction of (1) w (2) w,x (3) w,y
YY=XXX(1:3:end,:); %w
% YY=XXX(3:3:end,:);  %theta_x = w,y
% YY=-XXX(2:3:end,:); %theta_y = -w,x
%**************************************************************************

% figure(2)
% pdeplot(pp,ee,tt,'xydata',YY(:,1),'contour','on');
% hold on;

figure
subplot(1,2,1)
pdeplot(pp,ee,tt,'zdata',YY(:,1));
xlabel('x');
subplot(1,2,2)
pdeplot(pp,ee,tt,'zdata',YY(:,2));
xlabel('x');
% subplot(4,2,3)
% pdeplot(pp,ee,tt,'zdata',YY(:,3))
% % subplot(4,2,4)
% % pdeplot(pp,ee,tt,'zdata',YY(:,4))
% % subplot(4,2,5)
% % pdeplot(pp,ee,tt,'zdata',YY(:,5))
% % subplot(4,2,6)
% % pdeplot(pp,ee,tt,'zdata',YY(:,6))
% % subplot(4,2,7)
% % pdeplot(pp,ee,tt,'zdata',YY(:,7))
% % subplot(4,2,8)
% % pdeplot(pp,ee,tt,'zdata',YY(:,8))

%% THESIS PLOTS
% eigModes = [1-8];

zlimSet = [-0.05 0.05];
caxisSet = [-0.03 0.03];
viewSet = [40 15];
% viewSet = [90 90];
FntSz = 14;
figure(3)
subplot(3,2,1); hold on; zlim(zlimSet); grid minor;
pdeplot(pp,ee,tt,'xydata',YY(:,1),'zdata',YY(:,1),'colormap','viridis');
colormap(viridis);%caxis([-0.1 0.1]);
view(viewSet);
hold on;
plot3(p(1,BBnodes),p(2,BBnodes),0.*p(2,BBnodes),'ks','MarkerSize',2);
xlabel('x', 'Interpreter','Latex','FontSize',FntSz-3);
ylabel('y', 'Interpreter','Latex','FontSize',FntSz-3);
zlabel('z', 'Interpreter','Latex','FontSize',FntSz-3);
title(['$$\lambda_{1}=$$',num2str(round(NDfreq(1),2))],...
    'Interpreter','Latex','FontSize',FntSz-3);
set(gca,'FontSize',FntSz);set(gca,'TickLabelInterpreter','latex');

subplot(3,2,2); hold on; zlim(zlimSet);grid minor;
pdeplot(pp,ee,tt,'xydata',YY(:,2),'zdata',YY(:,2),'colormap','viridis');
colormap(viridis);%caxis(caxisSet);
view(viewSet);
hold on;
plot3(p(1,BBnodes),p(2,BBnodes),0.*p(2,BBnodes),'ks','MarkerSize',2);
xlabel('x', 'Interpreter','Latex','FontSize',FntSz-3);
ylabel('y', 'Interpreter','Latex','FontSize',FntSz-3);
zlabel('z', 'Interpreter','Latex','FontSize',FntSz-3);
title(['$$\lambda_{2}=$$',num2str(round(NDfreq(2),2))],...
    'Interpreter','Latex','FontSize',FntSz-3);
set(gca,'FontSize',FntSz);set(gca,'TickLabelInterpreter','latex');

subplot(3,2,3); hold on; zlim(zlimSet);grid minor;
pdeplot(pp,ee,tt,'xydata',YY(:,3),'zdata',YY(:,3),'colormap','viridis');
colormap(viridis);%caxis(caxisSet);
view(viewSet);
hold on;
plot3(p(1,BBnodes),p(2,BBnodes),0.*p(2,BBnodes),'ks','MarkerSize',2);
xlabel('x', 'Interpreter','Latex','FontSize',FntSz-3);
ylabel('y', 'Interpreter','Latex','FontSize',FntSz-3);
zlabel('z', 'Interpreter','Latex','FontSize',FntSz-3);
title(['$$\lambda_{3}=$$',num2str(round(NDfreq(3),2))],...
    'Interpreter','Latex','FontSize',FntSz-3);
set(gca,'FontSize',FntSz);set(gca,'TickLabelInterpreter','latex');

subplot(3,2,4); hold on; zlim(zlimSet);grid minor;
pdeplot(pp,ee,tt,'xydata',YY(:,4),'zdata',YY(:,4),'colormap','viridis');
colormap(viridis);%caxis(caxisSet);
view(viewSet);
hold on;
plot3(p(1,BBnodes),p(2,BBnodes),0.*p(2,BBnodes),'ks','MarkerSize',2);
xlabel('x', 'Interpreter','Latex','FontSize',FntSz-3);
ylabel('y', 'Interpreter','Latex','FontSize',FntSz-3);
zlabel('z', 'Interpreter','Latex','FontSize',FntSz-3);
title(['$$\lambda_{4}=$$',num2str(round(NDfreq(4),2))],...
    'Interpreter','Latex','FontSize',FntSz-3);
set(gca,'FontSize',FntSz);set(gca,'TickLabelInterpreter','latex');

subplot(3,2,5); hold on; zlim(zlimSet);grid minor;
pdeplot(pp,ee,tt,'xydata',YY(:,5),'zdata',YY(:,5),'colormap','viridis');
colormap(viridis);%caxis(caxisSet);
view(viewSet);
hold on;
plot3(p(1,BBnodes),p(2,BBnodes),0.*p(2,BBnodes),'ks','MarkerSize',2);
xlabel('x', 'Interpreter','Latex','FontSize',FntSz-3);
ylabel('y', 'Interpreter','Latex','FontSize',FntSz-3);
zlabel('z', 'Interpreter','Latex','FontSize',FntSz-3);
title(['$$\lambda_{5}=$$',num2str(round(NDfreq(5),2))],...
    'Interpreter','Latex','FontSize',FntSz-3);
set(gca,'FontSize',FntSz);set(gca,'TickLabelInterpreter','latex');

subplot(3,2,6); hold on; zlim(zlimSet);grid minor;
pdeplot(pp,ee,tt,'xydata',YY(:,6),'zdata',YY(:,6),'colormap','viridis');
colormap(viridis);%caxis(caxisSet);
view(viewSet);
hold on;
plot3(p(1,BBnodes),p(2,BBnodes),0.*p(2,BBnodes),'ks','MarkerSize',2);
xlabel('x', 'Interpreter','Latex','FontSize',FntSz-3);
ylabel('y', 'Interpreter','Latex','FontSize',FntSz-3);
zlabel('z', 'Interpreter','Latex','FontSize',FntSz-3);
title(['$$\lambda_{6}=$$',num2str(round(NDfreq(6),2))],...
    'Interpreter','Latex','FontSize',FntSz-3);
set(gca,'FontSize',FntSz);set(gca,'TickLabelInterpreter','latex');

%% contour maps
% 
zlimSet = [-0.05 0.05];
caxisSet = [-0.03 0.03];
% viewSet = [40 15];
% viewSet = [90 90];
FntSz = 14;
figure%(3)
subplot(3,2,1); hold on; zlim(zlimSet); grid minor;
pdeplot(pp,ee,tt,'xydata',YY(:,1),'zdata',YY(:,1),'colormap','viridis',...
    'contour','on');
colormap(viridis);%caxis([-0.1 0.1]);
% view(viewSet);
hold on;
plot3(p(1,BBnodes),p(2,BBnodes),0.*p(2,BBnodes),'ks','MarkerSize',2);
xlabel('x', 'Interpreter','Latex','FontSize',FntSz-3);
ylabel('y', 'Interpreter','Latex','FontSize',FntSz-3);
zlabel('z', 'Interpreter','Latex','FontSize',FntSz-3);
title(['$$\lambda_{1}=$$',num2str(round(NDfreq(1),2))],...
    'Interpreter','Latex','FontSize',FntSz-3);
set(gca,'FontSize',FntSz);set(gca,'TickLabelInterpreter','latex');

subplot(3,2,2); hold on; zlim(zlimSet);grid minor;
pdeplot(pp,ee,tt,'xydata',YY(:,2),'zdata',YY(:,2),'colormap','viridis',...
    'contour','on');
colormap(viridis);%caxis(caxisSet);
% view(viewSet);
hold on;
plot3(p(1,BBnodes),p(2,BBnodes),0.*p(2,BBnodes),'ks','MarkerSize',2);
xlabel('x', 'Interpreter','Latex','FontSize',FntSz-3);
ylabel('y', 'Interpreter','Latex','FontSize',FntSz-3);
zlabel('z', 'Interpreter','Latex','FontSize',FntSz-3);
title(['$$\lambda_{2}=$$',num2str(round(NDfreq(2),2))],...
    'Interpreter','Latex','FontSize',FntSz-3);
set(gca,'FontSize',FntSz);set(gca,'TickLabelInterpreter','latex');

subplot(3,2,3); hold on; zlim(zlimSet);grid minor;
pdeplot(pp,ee,tt,'xydata',YY(:,3),'zdata',YY(:,3),'colormap','viridis',...
    'contour','on');
colormap(viridis);%caxis(caxisSet);
% view(viewSet);
hold on;
plot3(p(1,BBnodes),p(2,BBnodes),0.*p(2,BBnodes),'ks','MarkerSize',2);
xlabel('x', 'Interpreter','Latex','FontSize',FntSz-3);
ylabel('y', 'Interpreter','Latex','FontSize',FntSz-3);
zlabel('z', 'Interpreter','Latex','FontSize',FntSz-3);
title(['$$\lambda_{3}=$$',num2str(round(NDfreq(3),2))],...
    'Interpreter','Latex','FontSize',FntSz-3);
set(gca,'FontSize',FntSz);set(gca,'TickLabelInterpreter','latex');

subplot(3,2,4); hold on; zlim(zlimSet);grid minor;
pdeplot(pp,ee,tt,'xydata',YY(:,4),'zdata',YY(:,4),'colormap','viridis',...
    'contour','on');
colormap(viridis);%caxis(caxisSet);
% view(viewSet);
hold on;
plot3(p(1,BBnodes),p(2,BBnodes),0.*p(2,BBnodes),'ks','MarkerSize',2);
xlabel('x', 'Interpreter','Latex','FontSize',FntSz-3);
ylabel('y', 'Interpreter','Latex','FontSize',FntSz-3);
zlabel('z', 'Interpreter','Latex','FontSize',FntSz-3);
title(['$$\lambda_{4}=$$',num2str(round(NDfreq(4),2))],...
    'Interpreter','Latex','FontSize',FntSz-3);
set(gca,'FontSize',FntSz);set(gca,'TickLabelInterpreter','latex');

subplot(3,2,5); hold on; zlim(zlimSet);grid minor;
pdeplot(pp,ee,tt,'xydata',YY(:,5),'zdata',YY(:,5),'colormap','viridis',...
    'contour','on');
colormap(viridis);%caxis(caxisSet);
% view(viewSet);
hold on;
plot3(p(1,BBnodes),p(2,BBnodes),0.*p(2,BBnodes),'ks','MarkerSize',2);
xlabel('x', 'Interpreter','Latex','FontSize',FntSz-3);
ylabel('y', 'Interpreter','Latex','FontSize',FntSz-3);
zlabel('z', 'Interpreter','Latex','FontSize',FntSz-3);
title(['$$\lambda_{5}=$$',num2str(round(NDfreq(5),2))],...
    'Interpreter','Latex','FontSize',FntSz-3);
set(gca,'FontSize',FntSz);set(gca,'TickLabelInterpreter','latex');

subplot(3,2,6); hold on; zlim(zlimSet);grid minor;
pdeplot(pp,ee,tt,'xydata',YY(:,6),'zdata',YY(:,6),'colormap','viridis',...
    'contour','on');
colormap(viridis);%caxis(caxisSet);
% view(viewSet);
hold on;
plot3(p(1,BBnodes),p(2,BBnodes),0.*p(2,BBnodes),'ks','MarkerSize',2);
xlabel('x', 'Interpreter','Latex','FontSize',FntSz-3);
ylabel('y', 'Interpreter','Latex','FontSize',FntSz-3);
zlabel('z', 'Interpreter','Latex','FontSize',FntSz-3);
title(['$$\lambda_{6}=$$',num2str(round(NDfreq(6),2))],...
    'Interpreter','Latex','FontSize',FntSz-3);
set(gca,'FontSize',FntSz);set(gca,'TickLabelInterpreter','latex');


error('er')
%%
% eigModes = [1-8];

zlimSet = [-0.05 0.05];
caxisSet = [-0.03 0.03];
viewSet = [40 15];
% viewSet = [90 90];

figure(3)
subplot(4,2,1); hold on; zlim(zlimSet); grid minor;
pdeplot(pp,ee,tt,'xydata',YY(:,1),'zdata',YY(:,1));
colormap(viridis);%caxis([-0.1 0.1]);
hold on;
plot3(p(1,BBnodes),p(2,BBnodes),0.*p(2,BBnodes),'ks','MarkerSize',4);
view(viewSet);
xlabel('x', 'Interpreter','Latex','FontSize',FntSz-3);
ylabel('y', 'Interpreter','Latex','FontSize',FntSz-3);
zlabel('z', 'Interpreter','Latex','FontSize',FntSz-3);
title(['$$\Omega_{1}=$$',num2str(round(NDfreq(1),2))],...
    'Interpreter','Latex','FontSize',FntSz-3);

subplot(4,2,2); hold on; zlim(zlimSet);grid minor;
pdeplot(pp,ee,tt,'xydata',YY(:,2),'zdata',YY(:,2),'colormap','viridis');
colormap(viridis);%caxis(caxisSet);
view(viewSet);
hold on;
plot3(p(1,BBnodes),p(2,BBnodes),0.*p(2,BBnodes),'ks','MarkerSize',4);
xlabel('x', 'Interpreter','Latex','FontSize',FntSz-3);
ylabel('y', 'Interpreter','Latex','FontSize',FntSz-3);
zlabel('z', 'Interpreter','Latex','FontSize',FntSz-3);
title(['$$\Omega_{2}=$$',num2str(round(NDfreq(2),2))],...
    'Interpreter','Latex','FontSize',FntSz-3);

subplot(4,2,3); hold on; zlim(zlimSet);grid minor;
pdeplot(pp,ee,tt,'xydata',YY(:,3),'zdata',YY(:,3),'colormap','viridis');
colormap(viridis);%caxis(caxisSet);
view(viewSet);
hold on;
plot3(p(1,BBnodes),p(2,BBnodes),0.*p(2,BBnodes),'ks','MarkerSize',4);
xlabel('x', 'Interpreter','Latex','FontSize',FntSz-3);
ylabel('y', 'Interpreter','Latex','FontSize',FntSz-3);
zlabel('z', 'Interpreter','Latex','FontSize',FntSz-3);
title(['$$\Omega_{3}=$$',num2str(round(NDfreq(3),2))],...
    'Interpreter','Latex','FontSize',FntSz-3);

subplot(4,2,4); hold on; zlim(zlimSet);grid minor;
pdeplot(pp,ee,tt,'xydata',YY(:,4),'zdata',YY(:,4),'colormap','viridis');
colormap(viridis);%caxis(caxisSet);
view(viewSet);
hold on;
plot3(p(1,BBnodes),p(2,BBnodes),0.*p(2,BBnodes),'ks','MarkerSize',4);
xlabel('x', 'Interpreter','Latex','FontSize',FntSz-3);
ylabel('y', 'Interpreter','Latex','FontSize',FntSz-3);
zlabel('z', 'Interpreter','Latex','FontSize',FntSz-3);
title(['$$\Omega_{4}=$$',num2str(round(NDfreq(4),2))],...
    'Interpreter','Latex','FontSize',FntSz-3);

subplot(4,2,5); hold on; zlim(zlimSet);grid minor;
pdeplot(pp,ee,tt,'xydata',YY(:,5),'zdata',YY(:,5),'colormap','viridis');
colormap(viridis);%caxis(caxisSet);
view(viewSet);
hold on;
plot3(p(1,BBnodes),p(2,BBnodes),0.*p(2,BBnodes),'ks','MarkerSize',4);
xlabel('x', 'Interpreter','Latex','FontSize',FntSz-3);
ylabel('y', 'Interpreter','Latex','FontSize',FntSz-3);
zlabel('z', 'Interpreter','Latex','FontSize',FntSz-3);
title(['$$\Omega_{5}=$$',num2str(round(NDfreq(5),2))],...
    'Interpreter','Latex','FontSize',FntSz-3);

subplot(4,2,6); hold on; zlim(zlimSet);grid minor;
pdeplot(pp,ee,tt,'xydata',YY(:,6),'zdata',YY(:,6),'colormap','viridis');
colormap(viridis);%caxis(caxisSet);
view(viewSet);
hold on;
plot3(p(1,BBnodes),p(2,BBnodes),0.*p(2,BBnodes),'ks','MarkerSize',4);
xlabel('x', 'Interpreter','Latex','FontSize',FntSz-3);
ylabel('y', 'Interpreter','Latex','FontSize',FntSz-3);
zlabel('z', 'Interpreter','Latex','FontSize',FntSz-3);
title(['$$\Omega_{6}=$$',num2str(round(NDfreq(6),2))],...
    'Interpreter','Latex','FontSize',FntSz-3);

subplot(4,2,7); hold on; zlim(zlimSet);grid minor;
pdeplot(pp,ee,tt,'xydata',YY(:,7),'zdata',YY(:,7),'colormap','viridis');
colormap(viridis);%caxis(caxisSet);
view(viewSet);
hold on;
plot3(p(1,BBnodes),p(2,BBnodes),0.*p(2,BBnodes),'ks','MarkerSize',4);
xlabel('x', 'Interpreter','Latex','FontSize',FntSz-3);
ylabel('y', 'Interpreter','Latex','FontSize',FntSz-3);
zlabel('z', 'Interpreter','Latex','FontSize',FntSz-3);
title(['$$\Omega_{7}=$$',num2str(round(NDfreq(7),2))],...
    'Interpreter','Latex','FontSize',FntSz-3);

subplot(4,2,8); hold on; zlim(zlimSet);grid minor;
colormap(viridis);%caxis(caxisSet);
hold on;
plot3(p(1,BBnodes),p(2,BBnodes),0.*p(2,BBnodes),'ks','MarkerSize',4);
pdeplot(pp,ee,tt,'xydata',YY(:,8),'zdata',YY(:,8),'colormap','viridis');
xlabel('x', 'Interpreter','Latex','FontSize',FntSz-3);
ylabel('y', 'Interpreter','Latex','FontSize',FntSz-3);
zlabel('z', 'Interpreter','Latex','FontSize',FntSz-3);
title(['$$\Omega_{8}=$$',num2str(round(NDfreq(8),2))],...
    'Interpreter','Latex','FontSize',FntSz-3);
view(viewSet);


%% Including contours

% eigModes = [1, 4, 8, 15];


zlimSet = [-0.05 0.05];
caxisSet = [-0.03 0.03];
viewSet = [40 15];
decimals = 4;
% viewSet = [90 90];

figure(4)
subplot(4,2,1); hold on; zlim(zlimSet);grid minor;
pdeplot(pp,ee,tt,'xydata',YY(:,1),'zdata',YY(:,1),'colormap','viridis')
colormap(viridis);%caxis([-0.1 0.1]);
view(viewSet);
xlabel('x', 'Interpreter','Latex','FontSize',FntSz-3);
ylabel('y', 'Interpreter','Latex','FontSize',FntSz-3);
zlabel('z', 'Interpreter','Latex','FontSize',FntSz-3);
title(['$$\Omega_{1}=$$',num2str(round(NDfreq(1),decimals))],...
    'Interpreter','Latex','FontSize',FntSz-3);

    subplot(4,2,2); hold on; zlim(zlimSet);grid minor;
    pdeplot(pp,ee,tt,'xydata',YY(:,1),'colormap','viridis','contour','on')
    colormap(viridis);%caxis([-0.1 0.1]);
    % view(viewSet);
    xlabel('x', 'Interpreter','Latex','FontSize',FntSz-3);
    ylabel('y', 'Interpreter','Latex','FontSize',FntSz-3);
    zlabel('z', 'Interpreter','Latex','FontSize',FntSz-3);
    title(['$$\Omega_{1}=$$',num2str(round(NDfreq(1),decimals))],...
        'Interpreter','Latex','FontSize',FntSz-3);

subplot(4,2,3); hold on; zlim(zlimSet);
pdeplot(pp,ee,tt,'xydata',YY(:,4),'zdata',YY(:,4),'colormap','viridis')
colormap(viridis);%caxis(caxisSet);
view(viewSet);grid minor;
xlabel('x', 'Interpreter','Latex','FontSize',FntSz-3);
ylabel('y', 'Interpreter','Latex','FontSize',FntSz-3);
zlabel('z', 'Interpreter','Latex','FontSize',FntSz-3);
title(['$$\Omega_{4}=$$',num2str(round(NDfreq(4),decimals))],...
    'Interpreter','Latex','FontSize',FntSz-3);

    subplot(4,2,4); hold on; zlim(zlimSet);
    pdeplot(pp,ee,tt,'xydata',YY(:,4),'colormap','viridis','contour','on')
    colormap(viridis);%caxis(caxisSet);
    grid minor;
    xlabel('x', 'Interpreter','Latex','FontSize',FntSz-3);
    ylabel('y', 'Interpreter','Latex','FontSize',FntSz-3);
    zlabel('z', 'Interpreter','Latex','FontSize',FntSz-3);
    title(['$$\Omega_{4}=$$',num2str(round(NDfreq(4),decimals))],...
        'Interpreter','Latex','FontSize',FntSz-3);

subplot(4,2,5); hold on; zlim(zlimSet);grid minor;
pdeplot(pp,ee,tt,'xydata',YY(:,8),'zdata',YY(:,8),'colormap','viridis')
colormap(viridis);%caxis(caxisSet);
view(viewSet);
xlabel('x', 'Interpreter','Latex','FontSize',FntSz-3);
ylabel('y', 'Interpreter','Latex','FontSize',FntSz-3);
zlabel('z', 'Interpreter','Latex','FontSize',FntSz-3);
title(['$$\Omega_{8}=$$',num2str(round(NDfreq(8),decimals))],...
    'Interpreter','Latex','FontSize',FntSz-3);

    subplot(4,2,6); hold on; zlim(zlimSet);grid minor;
     pdeplot(pp,ee,tt,'xydata',YY(:,8),'colormap','viridis','contour','on')
    colormap(viridis);%caxis(caxisSet);
    xlabel('x', 'Interpreter','Latex','FontSize',FntSz-3);
    ylabel('y', 'Interpreter','Latex','FontSize',FntSz-3);
    zlabel('z', 'Interpreter','Latex','FontSize',FntSz-3);
    title(['$$\Omega_{8}=$$',num2str(round(NDfreq(8),decimals))],...
        'Interpreter','Latex','FontSize',FntSz-3);

subplot(4,2,7); hold on; zlim(zlimSet);grid minor;
pdeplot(pp,ee,tt,'xydata',YY(:,12),'zdata',YY(:,12),'colormap','viridis')
colormap(viridis);%caxis(caxisSet);
view(viewSet);
xlabel('x', 'Interpreter','Latex','FontSize',FntSz-3);
ylabel('y', 'Interpreter','Latex','FontSize',FntSz-3);
zlabel('z', 'Interpreter','Latex','FontSize',FntSz-3);
title(['$$\Omega_{12}=$$',num2str(round(NDfreq(12),decimals))],...
    'Interpreter','Latex','FontSize',FntSz-3);

subplot(4,2,8); hold on; zlim(zlimSet);grid minor;
colormap(viridis);%caxis(caxisSet);
pdeplot(pp,ee,tt,'xydata',YY(:,12),'colormap','viridis','contour','on')
xlabel('x', 'Interpreter','Latex','FontSize',FntSz-3);
ylabel('y', 'Interpreter','Latex','FontSize',FntSz-3);
zlabel('z', 'Interpreter','Latex','FontSize',FntSz-3);
title(['$$\Omega_{12}=$$',num2str(round(NDfreq(12),decimals))],...
    'Interpreter','Latex','FontSize',FntSz-3);

% % 
% % figure(5)
% % subplot(4,2,1)
% % pdeplot(pp,ee,tt,'xydata',YY(:,1),'zdata',YY(:,1),'colormap','viridis')
% % subplot(4,2,2)
% % pdeplot(pp,ee,tt,'xydata',YY(:,1),'colormap','viridis','contour','on')
% % subplot(4,2,3)
% % pdeplot(pp,ee,tt,'xydata',YY(:,4),'zdata',YY(:,4),'colormap','viridis')
% % subplot(4,2,4)
% % pdeplot(pp,ee,tt,'xydata',YY(:,4),'colormap','viridis','contour','on')
% % subplot(4,2,5)
% % pdeplot(pp,ee,tt,'xydata',YY(:,8),'zdata',YY(:,8),'colormap','viridis')
% % subplot(4,2,6)
% % pdeplot(pp,ee,tt,'xydata',YY(:,8),'colormap','viridis','contour','on')
% % subplot(4,2,7)
% % pdeplot(pp,ee,tt,'xydata',YY(:,15),'zdata',YY(:,15),'colormap','viridis')
% % subplot(4,2,8)
% % pdeplot(pp,ee,tt,'xydata',YY(:,15),'colormap','viridis','contour','on')

%% LATEX TABLE 1a, 1b, 1c RESULTS

%% CCCC constant thickness
% a/b = 1.0
anal_Leissa1973=[35.992, 73.413, 73.413, 108.27, 131.64, 132.24]';

% load('CCCC_mesh334_constant_thick.mat'); 
% NDfreq(1:6)
fem334 = [ 36.0121,   73.6732,   73.6961,  109.1788,  132.9158,  134.8479]'; %thin plate a/h=0.001
% fem334 = [ 36.0086,   73.6559,   73.6788,  109.1383,  132.8553,  134.7860]'; % thick plate a/h=0.1
Diff334 = (anal_Leissa1973 - fem334)./anal_Leissa1973*100

% load('CCCC_mesh1336_constant_thick.mat'); 
% NDfreq(1:6)

fem1336 = [ 35.9927,   73.4735,   73.4812,  108.4826,  131.9773,  132.9123]';
Diff1336 = (anal_Leissa1973 - fem1336)./anal_Leissa1973*100

% load('CCCC_mesh5344_constant_thick.mat'); 
% NDfreq(1:6)

fem5344 = [ 35.9870,   73.4138,   73.4159,  108.2833,  131.6832,  132.3826]';
Diff5344 = (anal_Leissa1973 - fem5344)./anal_Leissa1973*100

%% SSSS constant thickness
% a/b = 1.0
anal_Leissa1973=[19.7392, 49.3480, 49.3480, 78.9568, 98.6960, 98.6960]';

% load('SSSS_mesh334_constant_thick.mat'); 
% NDfreq(1:6)
fem334 = [  19.8320,   49.8485,   49.9212,   80.5741,  100.6417,  101.5332]'; %thin plate a/h=0.001
Diff334 = (anal_Leissa1973 - fem334)./anal_Leissa1973*100

% load('SSSS_mesh1336_constant_thick.mat'); 
% NDfreq(1:6)

fem1336 = [19.7626,    49.4771,   49.4947,   79.3722,   99.2081,   99.4196]';
Diff1336 = (anal_Leissa1973 - fem1336)./anal_Leissa1973*100

% load('SSSS_mesh5344_constant_thick.mat'); 
% NDfreq(1:6)

fem5344 = [19.7451,   49.3804,   49.3847,   79.0610,   98.8250,   98.8773]';
Diff5344 = (anal_Leissa1973 - fem5344)./anal_Leissa1973*100

%% CFFF constant thickness
% a/b = 1.0
anal_Leissa1973=[3.4917, 8.5246, 21.429, 27.331, 31.111, 54.443]';

% load('CFFF_mesh334_constant_thick.mat'); 0.5986
% NDfreq(1:6)0.
fem334 = [    3.4674,    8.5045,   21.2803,   27.1772,   31.0194,   54.6274]'; %thin plate a/h=0.001
Diff334 = (anal_Leissa1973 - fem334)./anal_Leissa1973*100

% load('CFFF_mesh1336_constant_thick.mat'); 
% NDfreq(1:6)

fem1336 = [  3.4702,    8.5059,   21.2835,   27.1947,   30.9730,   54.3050]';
Diff1336 = (anal_Leissa1973 - fem1336)./anal_Leissa1973*100

% load('CFFF_mesh5344_constant_thick.mat'); 
% NDfreq(1:6)

fem5344 = [   3.4708,    8.5062,   21.2838,   27.1977,   30.9590,   54.2140]';
Diff5344 = (anal_Leissa1973 - fem5344)./anal_Leissa1973*100

%% LINEAR TAPER (
% Shufrin vs present
% CFFF (taper=0.5)
shufrin=[0.3859 ,0.7563 ,1.8485 , 1.9438 ,2.4184 , 4.0317 ];
present = [0.3813,    0.7510,    1.7513,    1.9800,    2.3520,    3.9302];
present25=[0.3816,    0.7541,    1.7598,    1.9924,    2.3659,    3.9625];
diff = (shufrin-present25)./shufrin*100
% CCCC (taper=0.25)
shufrin = [3.1767, 6.4650, 6.4782, 9.5610, 11.5702, 11.6375];
present = [3.1521,    6.3541,    6.3662,    9.3067,   11.1978,   11.2804];
present25=[3.1630,    6.3801,    6.3920,    9.3505,   11.2599,   11.3580];
diff = (shufrin-present25)./shufrin*100

