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
% =========================================================================
clear all;
close all;
clc;
debugOn=0;
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
% addpath('mesh/uniform-load-ansys');
% addpath('mesh/distributed-load-ansys');
addpath('mesh/FSI-AUV');
% load("ANSYS_rect_v0.mat");
% load("ANSYS_rect_v1.mat");
% load("ANSYS_rect_v2.mat");
% load("ANSYS_rect_v3.mat");

% load("wingAUV_mesh.mat");
% load("wingAUV_mesh2.mat");
% load("wingAUV_mesh3.mat");

% %Create Rectangular plate
% pderect([-5 5 -5 5],'C1')

% load("eurodyn5");
load("eurodyn4");
% load("eurodyn3");
% load("eurodyn2");
% load("eurodyn1");

disp(['*=====================================================*']);
disp(['*                  PLATE 2D-FEM                       *']);
disp(['*=====================================================*']);
disp([' MESH SIZE (DKT Triangles 9-dofs)= ',num2str(size(t,2))]);
%==========================================================================
%                             PRE-PROCESSOR
%==========================================================================
chord=0.33; %wing chord length
span=1; % wing span length
U =2.52;%m/s
fluid_dens=1025;%kg/m3
%% Material properties (SI)
% mass distribution [kg/m3]
% plate Young modulus [Pa]
% plate Poisson ratio
% plate thickness 
%% STEEL
m=7850;
E=200e9;
v=0.3;
% h=0.01;
% h=0.12*chord/2;
h=0.12*chord;
% h=0.09*chord;

%% ALUMINIUM
% m=2689.8;
% v=0.37;
% E=68.3e9; 
% h=0.12*chord/2;

%% POLYURETHANE RUBBER (gives 10% of chord-length displacement)
% Info: https://matmatch.com/learn/material/polyurethane
% m=1700;
% E=0.03e9;%1.88e9;
% v=0.4;
% h=0.12*chord;%/2; %NACA0012 average thickness profile
%
disp([' thick/width or thick/length < 0.1 (10%): thin plate'])
disp([' thick/width or thick/length > 0.1 (10%): moderately thick plate'])
disp([' > h/chord= ',num2str(h/chord)]);
disp([' > h/span = ',num2str(h/span)]);
%
%% Boundary conditions 
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
%
%% Forcing
% 1- concetrated load, 2- uniform load, 3- distributed load (mapping func)
lll=3;%2; %loading case
importFromFile=struct('toggle',1,'filename','bemDATA.mat');
%
P_load = 1; %[Pa] %pointing towards the Z-axis
% in ANSYS load pointing in the negative of Z-axis is positive
if lll==1
   Pxy=[5,5];%load position
end
%
%% Rigidity
% G= E/(2*(1+v));
% Dplate=E*h^3/(12*(1-v^2));
% shear correction factor
% sc=5/6; 
% phi=E/(sc*G)*(h/10)^2;
%
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
               P=(3)*(j-1)+i;
               LM(P,k)=ID(i,IEN(j,k));
           end
       end
    end
    
GEN=max(max(LM)); %Total number of nodal unknowns taking into account the 
% connectivity between the triangles
%
%% Essential BCs (enforces on the displacements "w" and slopes "bx", "by")   
Bound1=find(e(5,:)==1);
Bound2=find(e(5,:)==2);
Bound3=find(e(5,:)==3);
Bound4=find(e(5,:)==4);
%************************THIS IS THE ACTIVE BOUNDARY CONDITION
% COMMENT: The numbering is offered by the pdeModeler
% Bnodes= [Bound4, Bound1(1)]; %FULL EDGE
Bnodes=Bound1; %for distributed load from function ANSYS
% Bnodes=[Bound1 Bound2 Bound3 Bound4];
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
    Bdofs_SS=sort(Bdofs1); %Dofs for w and theta
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
end
% error('r')
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
Ng=3;
[xw]=TriGaussPoints(Ng);
%==========================================================================
% BENDING STIFFNESS MATRIX (3x3) FOR EACH TRIANGLE

%  thick=h*ones(1,Nelem);
[~,txxBEM]=Nonunif(x,y,IEN,p,e,t, chord, span, debugOn, importFromFile,...
    fluid_dens, U);
[BeSt2]=BendingStiffness2(E,v,txxBEM,h); %[3,3] matrix

% error('er')
% [BeSt]=BendingStiffness(E,v,h); %[3,3] matrix
% Ds=sc*G*h*eye(2,2);

xg=xw(:,1);yg=xw(:,2);
[l23,l31,l12,y12,y31,y23,x12,x31,x23,Area,a4,a5,a6,b4,b5,b6,...
    c4,c5,c6,d4,d5,d6,e4,e5,e6,C4,C5,C6,S4,S5,S6]=TrigElCoefsDKT(x,y,IEN);

% Calculation of shape functions and their derivatives at gauss points on
% the parent element
[ SF,DxsiSF,DetaSF,D2xsiSF,D2xsietaSF,D2etaSF ] = LNShapeFunDST(xg,yg);
[ SFm, DxsiSFm,DetaSFm] = LNShapeFunMassDST(xg,yg); 

%
%---> ax, ay, bx,by are constant for constant h (independednt of �,�)
% %%
[ GGDST,GGDKT] = matrixG();
GGin=inv(GGDST);
GGin2=inv(GGDKT);


% for i=1:Ng;
%     
%     xg_m(i,:)=x(IEN(1,:))*(1-xg(i)-yg(i))+x(IEN(2,:))*xg(i)+x(IEN(3,:))*yg(i);
%     
%     yg_m(i,:)=y(IEN(1,:))*(1-xg(i)-yg(i))+y(IEN(2,:))*xg(i)+y(IEN(3,:))*yg(i);
%     
% end

%************************************** DISTRIBUTED LOAD VIA MAPPING*******
if lll==3
    [Fx,~]=Nonunif(x,y,IEN,p,e,t, chord, span, 0, importFromFile,...
        fluid_dens, U);
%     Fx = -Fx; %for ansys
    if debugOn
        hold on;grid on;
        plot3(p(1,BBnodes),p(2,BBnodes),p(2,BBnodes).*0,'ks');
        xlabel('x-axis');
    end
end
%**************************************************************************

Fglob=zeros(GEN,1);

for kk=1:Nelem %for each element (iS THIS TRIANGLE 1 IN t?)
  
% %% Mass matrix
    [ Hm,HW] = massHmDKT(kk,l23,l31,l12,C4,C5,C6,S4,S5,S6,GGin);
    %% Integration
    kb=zeros(9,9); 
    mloc=zeros(9,9);
    mloc1=zeros(9,9);
    mloc2=zeros(9,9);
    %************************** ADDITION
        floc1=zeros(10,1);
    %***********************************
    [ Hxx,Hyy ] = rotationMass2(kk,a4,a5,a6,b4,b5,b6,c4,c5,c6,d4,d5,d6,e4,e5,e6 );
    [ Hxx1,Hyy1 ] = rotationMass(kk,l23,l31,l12,C4,C5,C6,S4,S5,S6);
    
    for ii=1:Ng % for each gauss point
    %% Stiffness
%        [Hx, Hy, Hx_xsi, Hx_eta,Hy_xsi, Hy_eta, Bb ] = ShapeFunDKT(ii,kk,a4,a5,a6,b4,b5,b6,c4,c5,c6,d4,d5,d6,e4,e5,e6,Area,SF,DxsiSF,DetaSF,y31,x31,y12,x12) ;
      [Hx, Hy, Hx_xsi,Hx_eta,Hy_xsi, Hy_eta, Bb ] = ShapeFunDKT2(ii,kk,C4,C5,C6,S4,S5,S6,Area, SF, DxsiSF,DetaSF,y31,x31,y12,x12,l23,l31,l12);
%       [Hx1, Hy1, Hx_xsi1, Hx_eta1,Hy_xsi1, Hy_eta1, Bb1 ] = ShapeFunDKT(ii,kk,a4,a5,a6,b4,b5,b6,c4,c5,c6,d4,d5,d6,e4,e5,e6,Area,SF,DxsiSF,DetaSF,y31,x31,y12,x12) ;
       
       [ Nm, HW,LW,L,HX3,HY3] = pseudoMassDKT(ii,kk,l23,l31,l12,y12,y31,x12,x31,Area,Hm,GGin,GGin2,xg, yg,IEN,x,y,Hx,Hy,SFm,C4,C5,C6,S4,S5,S6,Hxx,Hyy,HW )   ;
% %        [ Nm,LW,L] = pseudoMassDKT(ii,kk,l23,l31,l12,y12,y31,x12,x31,Area,Hm,GGin,GGin2,xg, yg,IEN,x,y,Hx,Hy,SFm,C4,C5,C6,S4,S5,S6,HW)   ;

       % we give one more dimension 
       % to the BeSt [3x3] original matrix
       kb=kb+Area(kk)*xw(ii,3)*(Bb'*BeSt2(:,:,kk)*Bb);
%        kb=kb+Area(kk)*xw(ii,3)*(Bb'*BeSt*Bb);
%        mloc=mloc+m*h*Area(kk)*xw(ii,3)*(Nm'*Nm);
       mloc=mloc+m*txxBEM(kk)*Area(kk)*xw(ii,3)*(Nm'*Nm);
       %****************************ADDITION
       floc1=floc1+Area(kk)*xw(ii,3)*(LW');
       %************************************
%         mloc=mloc+m*txxBEM(kk)*Area(kk)*xw(ii,3)*((HW'*(LW'*LW)*HW)+txxBEM(kk)^2/12*(Hx'*Hx)+txxBEM(kk)^2/12*(Hy'*Hy));
%        mloc=mloc+m*txxBEM(kk)*Area(kk)*xw(ii,3)*((HW'*(LW'*LW)*HW));
    end
      
    kloc=kb;
    %************************** ADDITION
%         floc=P*HW'*floc1;
    if lll==2 % uniform load
        floc=Area(kk)*P_load/3*[1 0 0 1 0 0 1 0 0]';% lumped mass approach for the uniform load
    elseif lll==3 %distributed load via mapping func
        floc=Area(kk)*Fx(kk)/3*[1 0 0 1 0 0 1 0 0]';
    end
    %***********************************
    Mg(:,kk)=[mloc(:,1);mloc(:,2);mloc(:,3);mloc(:,4);mloc(:,5);mloc(:,6);mloc(:,7);mloc(:,8);mloc(:,9)];

    Kg(:,kk)=[kloc(:,1);kloc(:,2);kloc(:,3);kloc(:,4);kloc(:,5);kloc(:,6);kloc(:,7);kloc(:,8);kloc(:,9)];
    %************************** ADDITION
    for q=1:9
         Fglob(LM(q,kk))=Fglob(LM(q,kk))+floc(q);
    end 
    %***********************************
end

%% Global assembly (uses the LM)
% COMMENT: The boundary conditions are enforced as extra equations
% in the sense of constraints in the present version

iii=1:9; %Elnodes Number of element nodes
iii=repmat(iii',1,9);

rr=iii';
Ig=LM(iii(:),:);  %LM ARRAY (Dofs per Element)X(Elements)   -   (Nnodes*Ndof)X(Nel)
Jg=LM(rr(:),:);
Kglob=sparse(Ig(:),Jg(:),Kg(:),GEN,GEN);
spy(Kglob)
Mglob=sparse(Ig(:),Jg(:),Mg(:),GEN,GEN);

BBnodes_old=BBnodes; %DIMITRA
BBnodes=Bdofs;%DIMITRA

kkk=zeros(length(BBnodes),GEN);
mmm=zeros(length(BBnodes),GEN);
Dofs=length(p)*3;

% Totbound=reshape(BBound',1,4*size(BBound,2));
%***************************ADDITION***************************************
for j=1:length(BBnodes)
    kkk(j,:)=[zeros(1,Bdofs(j)-1) 1 zeros(1,Dofs-Bdofs(j))];
end
%**************************************************************************
Kglob=[Kglob kkk'; kkk zeros(length(BBnodes))];
Mglob=[Mglob mmm'; mmm zeros(length(BBnodes))];

% [XX,lamM,flag]=eigs(Kglob,Mglob,15,'sm');
% cc=sort(diag(lamM));
% 
% freq=sqrt(sort(diag(lamM),'ascend'))./(2*pi);
% 
% error('yy')

%FORCING AND SOLUTION
if lll==1
    Fglob=zeros(length(Kglob),1);
    Fglob(ID(1,PNODE))=P_load;
    U=Kglob\Fglob; %SOLVE SPARSE SYSTEM OF EQUATIONS
%hughe [ch.9] newmark - 2nd order
else
    Fglob1=[Fglob; zeros(length(BBnodes),1)];
    U=Kglob\Fglob1;%SOLVE SPARSE SYSTEM OF EQUATIONS
end
%==========================================================================
%                             POST-PROCESSOR
%==========================================================================
u=U(1:GEN); % the vector of nodal unknowns (w1;bx1;by1;....wN;bxN;byN)
%
w=u(1:3:end);   % vertical displacement
bx=u(2:3:end);  % rotation x
by=u(3:3:end);  % rotation y

BBnodes=BBnodes_old;%DIMITRA

%% PLOTS FOR THE NUMERICAL DISPLACEMENT
% CHECK IF BOUNDARY CONDITIONS ARE SATISFIED..
figure;
subplot(1,2,1);hold on;grid on;
plot3(pp(1,BBnodes),pp(2,BBnodes),bx(BBnodes),'ks','MarkerSize',3);
pdeplot(pp,ee,tt,'zdata',bx,'colormap','copper');
shading interp
xlabel('x-axis');ylabel('y-axis');zlabel('bx');
title('bx rotation','FontWeight','normal');
subplot(1,2,2);
% figure;
hold on;grid on;
plot3(pp(1,BBnodes),pp(2,BBnodes),by(BBnodes),'ks','MarkerSize',3);
pdeplot(pp,ee,tt,'zdata',by,'colormap','copper');
shading interp
xlabel('x-axis');ylabel('y-axis');zlabel('by');
title('by rotation','FontWeight','normal');

    figure;
    subplot(1,3,[1 2]);hold on;grid on;
    plot3(pp(1,BBnodes),pp(2,BBnodes),w(BBnodes),'ks','MarkerSize',3);
    pdeplot(pp,ee,tt,'XYData',w,"ZData",w,'colormap','jet');
    colorbar;shading interp;view([25 25]);%axis equal;
    zlim([-2.5*max(max(abs(w))) 2.5*max(max(abs(w)))])
    xlabel('x-axis');ylabel('y-axis');zlabel('w [m]');
%     title('w displacement','FontWeight','normal');
    subplot(1,3,3);hold on;grid on;
    pdeplot(pp,ee,tt,'XYData',w,'colormap','jet','contour','on');
    colorbar;shading interp;
    xlabel('x-axis');ylabel('y-axis');
    title('(contour)','FontWeight','normal');

figure;hold on;grid on;
plot3(pp(1,:),pp(2,:),w,'o','MarkerSize',3);
plot3(pp(1,BBnodes),pp(2,BBnodes),w(BBnodes),'ks','MarkerSize',3);
xlabel('x-axis');ylabel('y-axis');zlabel('w [m]');
title('w displacement (mesh & BCs)','FontWeight','normal');
% axis equal;
% zlim([-chord chord])
view([25 25])

maximum_w=max(max(abs(w)))
wdivh=max(max(abs(w)))/h
wdivchord=max(max(abs(w)))/chord*100

num_mesh123=[2.9922e-06;4.1277e-06];
ansys=[2.9440e-06;3.8859e-06];%6.2884e-06;8.1023e-06;1.2115e-05];%shell 181 cfff

(num_mesh123-ansys)./ansys*100


% =========================================================================
%                           GENERAL COMMENTS
%
% COMMENT: Making the mesh more dense closer to the boundaries with clamped
% edges could solve the problem of fine grids?
% COMMENT: Discuss with Aggelina the bx, by fields in satisfying the BCs
%==========================================================================

error('Up')
% 
P=P_load;
if lll==1 % concentrated load at the midddle of the rectangular plate
    %wanal=P/(16*pi*Dplate).*((a^2-rr1.^2)+2*rr1.^2.*log(rr1./a));
    b=10;
    a=10;
    summ=zeros(size(x));

    for m=1:2:150
        for n=1:2:150
    summ=summ+sin(m*pi/a.*x).*sin(n*pi/b.*y).*sin(m*pi/(a).*(a/2)).*sin(n*pi/(b).*(b/2))/((m^2/a^2+n^2/b^2))^2;
        end
    end
    wanal=4*P/(Dplate*a*b*pi^4)*summ;

else  % uniformly distributed load

    b=10;
    a=10;
    summ=zeros(size(x));
    for m=1:2:150
        for n=1:2:16
    summ=summ+sin((m*pi*x)/a).*sin((n*pi*y)/b)/(m*n*(m^2/a^2+n^2/b^2)^2);
        end
    end
    wanal=16*P/(Dplate*pi^6)*summ;
end

figure(2)
pdeplot(p,e,t,'zdata',wanal,'colormap','copper')
title('Concetrated Load Case')

figure(3)
subtitle('Concetrated Load Case')
subplot(2,1,1)
pdeplot(p,e,t,'xydata',w,'colormap','jet','contour','on')
%title('FEM')
axis equal
subplot(2,1,2)
pdeplot(p,e,t,'xydata',wanal,'colormap','jet','contour','on')
title('Navier Solution')
axis equal


MAx=max(max(abs(w)))
MAxAnal=max(max(wanal))
DEVpercent=(MAx-MAxAnal)/MAxAnal*100

error('er')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
% [XX,lam,flag]=eigs(Kglob,15,'sm');
% [XX,lam,flag]=eigs(Kglob,15,'sm');
% [XX,lam,flag]=eigs(Kglob,15,'sm');
[XX,lamM,flag]=eigs(Kglob,Mglob,15,'sm');
cc=sort(diag(lamM));

freq=sqrt(sort(diag(lamM),'ascend'))./(2*pi);


DD=E*mean(txxBEM)^3/(12*(1-v^2));
cl=sqrt(m*mean(txxBEM)/DD);
% DD=E*h^3/(12*(1-v^2));
% cl=sqrt(m*h/DD);
NDfreq=cl*freq*2*pi*100
XXX=XX(1:GEN,:);
YY=XXX(1:3:end,:);

figure(1)
pdeplot(pp,ee,tt,'xydata',YY(:,1),'contour','on')

figure(2)
subplot(4,2,1)
pdeplot(pp,ee,tt,'zdata',YY(:,1),'colormap','copper')
subplot(4,2,2)
pdeplot(pp,ee,tt,'zdata',YY(:,2),'colormap','copper')
subplot(4,2,3)
pdeplot(pp,ee,tt,'zdata',YY(:,3),'colormap','copper')
subplot(4,2,4)
pdeplot(pp,ee,tt,'zdata',YY(:,4),'colormap','copper')
subplot(4,2,5)
pdeplot(pp,ee,tt,'zdata',YY(:,5),'colormap','copper')
subplot(4,2,6)
pdeplot(pp,ee,tt,'zdata',YY(:,6),'colormap','copper')
subplot(4,2,7)
pdeplot(pp,ee,tt,'zdata',YY(:,7),'colormap','copper')
subplot(4,2,8)
pdeplot(pp,ee,tt,'zdata',YY(:,8),'colormap','copper')

figure(3)
subplot(4,2,1)
pdeplot(pp,ee,tt,'xydata',YY(:,1),'zdata',YY(:,1),'colormap','copper')
subplot(4,2,2)
pdeplot(pp,ee,tt,'xydata',YY(:,2),'zdata',YY(:,2),'colormap','copper')
subplot(4,2,3)
pdeplot(pp,ee,tt,'xydata',YY(:,3),'zdata',YY(:,3),'colormap','copper')
subplot(4,2,4)
pdeplot(pp,ee,tt,'xydata',YY(:,4),'zdata',YY(:,4),'colormap','copper')
subplot(4,2,5)
pdeplot(pp,ee,tt,'xydata',YY(:,5),'zdata',YY(:,5),'colormap','copper')
subplot(4,2,6)
pdeplot(pp,ee,tt,'xydata',YY(:,6),'zdata',YY(:,6),'colormap','copper')
subplot(4,2,7)
pdeplot(pp,ee,tt,'xydata',YY(:,7),'zdata',YY(:,7),'colormap','copper')
subplot(4,2,8)
pdeplot(pp,ee,tt,'xydata',YY(:,8),'zdata',YY(:,8),'colormap','copper')


figure(4)
subplot(4,2,1)
pdeplot(pp,ee,tt,'xydata',YY(:,1),'zdata',YY(:,1),'colormap','copper')
subplot(4,2,2)
pdeplot(pp,ee,tt,'xydata',YY(:,1),'colormap','copper','contour','on')
subplot(4,2,3)
pdeplot(pp,ee,tt,'xydata',YY(:,4),'zdata',YY(:,4),'colormap','copper')
subplot(4,2,4)
pdeplot(pp,ee,tt,'xydata',YY(:,4),'colormap','copper','contour','on')
subplot(4,2,5)
pdeplot(pp,ee,tt,'xydata',YY(:,8),'zdata',YY(:,5),'colormap','copper')
subplot(4,2,6)
pdeplot(pp,ee,tt,'xydata',YY(:,8),'colormap','copper','contour','on')
subplot(4,2,7)
pdeplot(pp,ee,tt,'xydata',YY(:,15),'zdata',YY(:,15),'colormap','copper')
subplot(4,2,8)
pdeplot(pp,ee,tt,'xydata',YY(:,15),'colormap','copper','contour','on')


