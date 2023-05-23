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
clear all;
close all;
clc;

MODAL_ANALYSIS = 1;
DYNAMIC_ANALYSIS = 1;

tstart = tic;   

debugOn=1;
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
addpath('funcFEM');
addpath('dataFSI');
addpath('mesh');
addpath('mesh/heathcote');

load('mesh_h1_half');

%Create Rectangular plate
% pderect([0 0.1 -0.3 0.3],'C1')
% pderect([0 0.1 -0.3 0],'C1')

% file1995 = 'bemDATA_h10_r';
% file1995 = 'bemDATA_h15_r';
file1995 = 'bemDATA_h182_r';
load(file1995);

disp(['*=====================================================*']);
disp(['*                  PLATE 2D-FEM                       *']);
disp(['*=====================================================*']);
disp([' MESH SIZE (DKT Triangles 9-dofs)= ',num2str(size(t,2))]);
%==========================================================================
%                             PRE-PROCESSOR
%==========================================================================
chord=inData.cRoot;%0.33; %wing chord length
span=inData.span;%1; % wing span length
Uvel =inData.U;%2.52;%m/s
fluid_dens=1025;%kg/m3

%% Material properties (SI)
% mass distribution [kg/m3]
% plate Young modulus [Pa]
% plate Poisson ratio
% plate thickness 
%% STEEL
m=7850;
E=210e9;
v=0.28;
h=1/1000;%0.12*chord;

%% ALUMINIUM
% m=2689.8;
% v=0.37;
% E=68.3e9; 
% h=0.001;%0.12*chord/2;

%% POLYURETHANE RUBBER (gives 10% of chord-length displacement)
% Info: https://matmatch.com/learn/material/polyurethane
% m=1700;
% E=0.5e7;%0.03e9;%1.88e9;
% v=0.4;
% h=0.12*chord/2;%/2; %NACA0012 average thickness profile
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
importFromFile=struct('toggle',1,'filename',file1995);
%
P_load = 1; %[Pa] %pointing towards the Z-axis
% in ANSYS load pointing in the negative of Z-axis is positive
if lll==1
%    Pxy=[5,5];%load position
   Pxy=[0.05,-0.15];%load position
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
           P=(3)*(j-1)+i; %the 9 dofs per triangle
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
% Bnodes = [Bound4(1), Bound3];
% Bnodes = [Bound2, Bound2];
Bnodes=Bound3; %for distributed load from function ANSYS
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
Ng=3; %TO-DO
[xw]=TriGaussPoints(Ng);

%==========================================================================
% BENDING STIFFNESS MATRIX (3x3) FOR EACH TRIANGLE
d=100;
%  thick=h*ones(1,Nelem);
[loadFem,txxBEM]=Nonunif(x,y,IEN,p,e,t, chord, span, 1, importFromFile,...
    fluid_dens, Uvel, h, d);
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

%
%---> ax, ay, bx,by are constant for constant h (independednt of î,ç)
% %%
[GGDST,GGDKT] = matrixG();

GGin=inv(GGDST);
GGin2=inv(GGDKT);

%************************************** DISTRIBUTED LOAD VIA MAPPING*******
if lll==3
    [Fx,~]=Nonunif(x,y,IEN,pp,ee,tt, chord, span, 0, importFromFile,...
        fluid_dens, Uvel,h,  d);
%     Fx = -Fx; %for ansys
%     if debugOn
%         hold on;grid on;
%         plot3(pp(1,BBnodes),pp(2,BBnodes),pp(2,BBnodes).*0,'ks');
%         xlabel('x-axis');
%     end
end
%**************************************************************************

% error('pp')
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
      [Hx, Hy, Hx_xsi,Hx_eta,Hy_xsi, Hy_eta, Bb ] = ShapeFunDKT2(ii,kk,C4,C5,C6,S4,S5,S6,Area,...
          SF, DxsiSF,DetaSF,y31,x31,y12,x12,l23,l31,l12);
%       [Hx1, Hy1, Hx_xsi1, Hx_eta1,Hy_xsi1, Hy_eta1, Bb1 ] = ShapeFunDKT(ii,kk,a4,a5,a6,b4,b5,b6,c4,c5,c6,d4,d5,d6,e4,e5,e6,Area,SF,DxsiSF,DetaSF,y31,x31,y12,x12) ;
       
       [ Nm, HW,LW,L,HX3,~] = pseudoMassDKT(ii,kk,l23,l31,l12,y12,y31,x12,x31,Area,...
           Hm,GGin,GGin2,xg, yg,IEN,x,y,Hx,Hy,SFm,C4,C5,C6,S4,S5,S6,Hxx,Hyy,HW );
% %        [ Nm,LW,L] = pseudoMassDKT(ii,kk,l23,l31,l12,y12,y31,x12,x31,Area,Hm,GGin,GGin2,xg, yg,IEN,x,y,Hx,Hy,SFm,C4,C5,C6,S4,S5,S6,HW)   ;

       % we give one more dimension 
       % to the BeSt [3x3] original matrix
       kb=kb+Area(kk)*xw(ii,3)*(Bb'*BeSt2(:,:,kk)*Bb);

%        kb=kb+Area(kk)*xw(ii,3)*(Bb'*BeSt*Bb);
%        mloc=mloc+m*h*Area(kk)*xw(ii,3)*(Nm'*Nm);
%        mloc=mloc+m*txxBEM(kk)*Area(kk)*xw(ii,3)*(Nm'*Nm);
       %****************************ADDITION
       floc1=floc1+Area(kk)*xw(ii,3)*(LW');
       %************************************
       mloc=mloc+m*txxBEM(kk)*Area(kk)*xw(ii,3)*((HW'*(LW'*LW)*HW)+txxBEM(kk)^2/12*(Hx'*Hx)+txxBEM(kk)^2/12*(Hy'*Hy));
%        mloc=mloc+m*txxBEM(kk)*Area(kk)*xw(ii,3)*((HW'*(LW'*LW)*HW));
    end
    kloc=kb;
    %************************** ADDITION
%     if lll==1
%         floc=P*HW'*floc1;
%         floc=P_load*HW'*floc1;
    if lll==2 % uniform load
        floc=Area(kk)*P_load/3*[1 0 0 1 0 0 1 0 0]';% lumped mass approach for the uniform load
    elseif lll==3 %distributed load via mapping func
        floc=Area(kk)*Fx(kk)/3*[1 0 0 1 0 0 1 0 0]';
    end
    %***********************************
    Mg(:,kk)=[mloc(:,1);mloc(:,2);mloc(:,3);mloc(:,4);mloc(:,5);mloc(:,6);mloc(:,7);mloc(:,8);mloc(:,9)];

    Kg(:,kk)=[kloc(:,1);kloc(:,2);kloc(:,3);kloc(:,4);kloc(:,5);kloc(:,6);kloc(:,7);kloc(:,8);kloc(:,9)];
    
    if lll == 2 || lll == 3
        %************************** ADDITION
        for q=1:9
             Fglob(LM(q,kk))=Fglob(LM(q,kk))+floc(q);
        end 
        %***********************************
    end
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
    
% error('er')
% A large condition number means that the matrix is close to being singular. 
% Kglobfull=full(Kglob);
% Mglobfull=full(Mglob);
% 
% cond(Kglobfull)
% 
% % error('er')
% clear Kglob_dense Mglob_dense
% Kglob_dense = zeros(GEN,GEN);
% Mglob_dense = zeros(GEN,GEN);
% for ii = 1:size(Ig,1)
%     for jj = 1:size(Ig,2)
% %         cnt = cnt + 1;
%         Kglob_dense(Ig(ii,jj),Jg(ii,jj)) = Kglob_dense(Ig(ii,jj),Jg(ii,jj)) + Kg(ii,jj);
%         Mglob_dense(Ig(ii,jj),Jg(ii,jj)) = Mglob_dense(Ig(ii,jj),Jg(ii,jj)) + Mg(ii,jj);
%     end
% end

% disp('===================Kglobfull=============================')
% Kglobfull(1:10,1:10)
% disp('===================Kglob_dense===========================')
% Kglob_dense(1:10,1:10)
% %
% disp('===================Mglobfull=============================')
% Mglobfull(1:10,1:10)
% disp('===================Mglob_dense===========================')
% Mglob_dense(1:10,1:10)

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

% Kglob_dense2=[Kglob_dense kkk'; kkk zeros(length(BBnodes))];
% Kglob_dense2(GEN+1:end,1:27)
% Kglob_dense2(1:27,GEN+1:end)

%FORCING AND SOLUTION
if lll==1
    Fglob(ID(1,PNODE))=P_load;
end

Fglob1=[Fglob; zeros(length(Bdofs),1)];

Fglob1(1:15)'

U = mldivide(Kglob,Fglob1); %or backslash

telapsed = toc(tstart);

% BBnodes=BBnodes_old;%DIMITRA

%==========================================================================
%                            POST-PROCESSOR
%==========================================================================
u=U(1:GEN); % the vector of nodal unknowns (w1;bx1;by1;....wN;bxN;byN)
%
w=u(1:3:end);   % vertical displacement
bx=u(2:3:end);  % rotation x
by=u(3:3:end);  % rotation y

figure;
subplot(1,3,[1 2]);hold on;grid on;
plot3(pp(1,BBnodes),pp(2,BBnodes),w(BBnodes),'ks','MarkerSize',3);
hh=pdeplot(pp,ee,tt,'XYData',w,"ZData",w);
colorbar;
shading interp;
colormap(viridis);
view([25 25]);%axis equal;
zlim([-2.5*max(max(abs(w))) 2.5*max(max(abs(w)))])
xlabel('x-axis');ylabel('y-axis');zlabel('w [m]');
%     title('w displacement','FontWeight','normal');
title(['Static solution: lll=',num2str(lll)]); 
subplot(1,3,3);hold on;grid on;
pdeplot(pp,ee,tt,'XYData',w,'colormap',viridis,'contour','on');
colorbar;shading interp;
xlabel('x-axis');ylabel('y-axis');
title('(contour)','FontWeight','normal');

max(abs(w))/inData.a3;

save solMatlab U

if MODAL_ANALYSIS == 1
    [XX,lamM,flag]=eigs(Kglob,Mglob,5,'sm');
    cc=sort(diag(lamM));

    freq=sqrt(sort(diag(lamM),'ascend'))./(2*pi);

    freq'
end

% error('er')

%% TIME-MARCHING
d_fromStatic = d;
% 
[FxDYN,~]=Nonunif(x,y,IEN,pp,ee,tt, chord, span, 0,....
    importFromFile,fluid_dens, Uvel, h, d_fromStatic);

if DYNAMIC_ANALYSIS == 1

    newmark = 0;
    implicitEuler = 0;
    crankNicolson = 1;
   
    T=2*pi/inData.omega3;%sec
    % wf=2*pi/T; %rad/s
    ddt=inData.dt;%T/100; %time-step
    t=[0:ddt:(inData.Nper)*T];%[0:h:2*T]; %time [sec]
    
    Ntimesteps = ceil((inData.Nper)*T/ddt)+1
    length(t)

    d=1; %starting point

    sizeM=size(Mglob,1);

    q=zeros(sizeM,length(t)); %displacement unknown vector (previous U)
    qdot=zeros(sizeM,length(t)); %velocity
    C = 0.*Mglob;
    % 
    [ C , res_Freq, a, b] = RayleighDamping( [], [], [], [], [], Kglob, Mglob, 1);
    % a
    % b
    Cfull=full(C);
    Cfull(1:10,1:10);

    
    % C=0.005*Mglob + 0.005*Kglob;   %a litte damping helps crank nicolson/newmark
% %     [Fx,~]=Nonunif(x,y,IEN,pp,ee,tt, chord, span, 0,....
% %             importFromFile,fluid_dens, Uvel, h, d);
% %     [Fglob_t] = createFglob(lll,GEN, Nelem,P_load, Fx,Area,LM,Bdofs);
    [Fglob_t] = createFglob(lll,GEN, Nelem,P_load, FxDYN*cos(inData.omega3*t(d)),Area,LM,Bdofs);
    Fm = [Fglob_t; zeros(sizeM,1)];
    size(Fm)

    if d==1
        G = zeros(length(Fm),length(t));
    end
    G(:,d) = Fm;
    
%     error('er')

    u=[qdot;q];

    solution = struct( 'w',[], 'bx',[], 'by', [],...
            'w_dot',[], 'bx_dot',[], 'by_dot', [],...
            'uu', [], 'uu_dot',[]);    

%     error('er')
    if newmark
        qdot2=zeros(sizeM,length(t)); %acceleration 
        beta = 0.25;
        gamma = 0.5;
        [Fglob_t] = createFglob(lll,GEN, Nelem,P_load,Fx,Area,LM,Bdofs);

        %initialization    
        AA=Mglob + gamma*ddt*C + ddt^2*beta*Kglob;
        BB=Fglob_t - C*qdot(:,1)- Kglob*q(:,1);
        qdot2(:,1)=AA\BB;
    end

    if newmark
        for d = 1:length(t)-1

            d
            % Update load vector
            [Fx,~]=Nonunif(x,y,IEN,pp,ee,tt, chord, span, 0, importFromFile,fluid_dens, Uvel, h, d+1);
    %         [Fx,~]=Nonunif(x,y,IEN,pp,ee,tt, chord, span, 0, importFromFile,fluid_dens, U, d);
            [Fglob_t] = createFglob(lll,GEN, Nelem,P_load,Fx,Area,LM,Bdofs);

            pr_vel = qdot(:,d)+(1-gamma)*ddt*qdot2(:,d);% + gamma*hhh*qdot2(:,d);
            pr_disp = q(:,d)+ddt*qdot(:,d)+ddt^2*(1/2-beta)*qdot2(:,d);%+hhh^2*beta*qdot2(:,d);

            AA = Mglob + gamma*ddt*C + ddt^2*beta*Kglob;
            BB=Fglob_t - C*pr_vel- Kglob*pr_disp;
            qdot2(:,d+1) = AA\BB; 

            q(:,d+1) = pr_disp+ddt^2*beta*qdot2(:,d+1);
            qdot(:,d+1) = pr_vel + gamma*ddt*qdot2(:,d+1);
            %
            u(:,d+1)=[qdot(:,d+1);q(:,d+1)];
            [solution] = solutionRetriever(GEN, sizeM, d+1, length(t), u, solution);%[w,bx,by]

        end
    %    
    %    
    else
    %
    %
    
%     error('er')
    % Am = [Mglob, sparse(sizeM,sizeM); sparse(sizeM,sizeM), speye(sizeM,sizeM)];
    % Bm = [C, Kglob; -speye(sizeM,sizeM), sparse(sizeM,sizeM)];
        for d = 1:10%length(t)-1
            d
            % Update load vector
% %             [Fx,~]=Nonunif(x,y,IEN,pp,ee,tt, chord, span, 0, importFromFile,fluid_dens, Uvel, h, d+1);
% %             
% %             [Fglob_t] = createFglob(lll,GEN, Nelem,P_load,Fx,Area,LM,Bdofs);
            [Fglob_t] = createFglob(lll,GEN, Nelem,P_load, FxDYN*cos(inData.omega3*t(d+1)),Area,LM,Bdofs);
            Fm = [Fglob_t; zeros(sizeM,1)];
            G(:,d+1) = Fm;

            theta = implicitEuler*(1) + crankNicolson*(1/2);
            u(:,d+1) = timeIntegration(u, d+1, GEN, Mglob, Kglob, C, G, ddt, theta); %[w,bx,by,lambda]
            usol=u(1:10,d+1)'

%             error('Testing');
            [solution] = solutionRetriever(GEN, sizeM, d+1, length(t), u, solution);%[w,bx,by]
        end
    %
    %
    %
    end
    
    
    G(1:10,1:10)
    u(1:10,1:10)

error('er')
    hmax = max(max(abs(solution.w)))

    (inData.a3 + hmax)/inData.a3

    chord

    w = solution.w(:,d);

    figure;
    subplot(1,3,[1 2]);hold on;grid on;
    plot3(pp(1,BBnodes),pp(2,BBnodes),w(BBnodes),'ks','MarkerSize',3);
    hh=pdeplot(pp,ee,tt,'XYData',w,"ZData",w,'colormap','jet');
    colorbar;shading interp;view([25 25]);%axis equal;
    % zlim([-2.5*max(max(abs(w))) 2.5*max(max(abs(w)))])
    xlabel('x-axis');ylabel('y-axis');zlabel('w [m]');
    %     title('w displacement','FontWeight','normal');
    subplot(1,3,3);hold on;grid on;
    pdeplot(pp,ee,tt,'XYData',w,'colormap',viridis,'contour','on');
    colorbar;shading interp;
    xlabel('x-axis');ylabel('y-axis');
    title('(contour)','FontWeight','normal');

    % save FEM_newmark

    % save FEM_sol_h15_r_h2
    save FEM_sol_h182_r_h2

%     error('er')
    % save FEM_sol_h05_r_h1

    debugOn=0;
    if debugOn
        Name = 'dynamicFEM_h182_r';
        animatorFunc(solution,length(t),ddt,pp,ee,tt,BBnodes,Name,inData)
    end

end

error('See below for comparisons with static case (analytic solution) & modal analysis!!!')


% =========================================================================
%                           GENERAL COMMENTS
%
% COMMENT: Making the mesh more dense closer to the boundaries with clamped
% edges could solve the problem of fine grids?
% COMMENT: Discuss with Aggelina the bx, by fields in satisfying the BCs
%==========================================================================
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




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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


