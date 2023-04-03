% =========================================================================
% Discrete Kirchoff Triangle (PLATE ELEMENT)
% References : Batoz et al (1989), SYDENSTRICKER et al (1995)
%
% Small strains, Isotropic material, Homogeneous properties
%
% Written by: A. Karperaki
% Modified by: D. Anevlavi
% =========================================================================
clear all;
close all;
clc;

addpath('mesh');

%% Boundary conditions 
CC=1; % boundary condition toggle

%% Forcing
lll=2; % 1- concetrated load, 2- distributed load

%% Material properties (SI)
% mass distribution [kg/m3]
m=7800;
% plate Young modulus [Pa]
E=200e9;
% plate Poisson ratio 
v=0.3;
% plate thickness 
h=0.01;
h/10
% Isotropic case
G= E/(2*(1+v));
Dplate=E*h^3/(12*(1-v^2));
% shear correction factor
sc=5/6; 
phi=E/(sc*G)*(h/10)^2;

%% Generate mesh - ID, IEN, LM 
%Create Rectangular plate

% pderect([-5 5 -5 5],'C1')
% save demyMesh p e t
% load Geometry1.mat
load mesh4%Sparse
%
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


%***************************ADDITION
tol=1;
I=find(abs(x-5)<tol & abs(y-5)<tol);
rr=sqrt((x(I)-5).^2+(y(I)-5).^2);
[II,JJ]=min(rr);
PNODE=I(JJ);
%***********************************

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


P = 1;

%% Essential BCs (enforces on the displacements "w" and slopes "bx", "by")   
Bound1=find(e(5,:)==1);
Bound2=find(e(5,:)==2);
Bound3=find(e(5,:)==3);
Bound4=find(e(5,:)==4);
Bnodes=[Bound1 Bound2 Bound3 Bound4];
 for i=1:length(Bnodes)
        BBnodes(i)=e(1,Bnodes(i));
 end
    BBnodes=sort(BBnodes);  
    Bdofs1=ID(1,BBnodes);
    Bdofs2=ID(2,BBnodes);
    Bdofs3=ID(3,BBnodes);
    Bdofs=sort([Bdofs1 Bdofs2 Bdofs3]); % % Dofs for w and theta

%% GAUSS POINTS AND WEIGHTS FOR TRIANGLES
% Number of GP
Ng=3;
[xw]=TriGaussPoints(Ng);
%==========================================================================
% BENDING STIFFNESS MATRIX (3x3) FOR EACH TRIANGLE

%  thick=h*ones(1,Nelem);
[BeSt]=BendingStiffness(E,v,h); %[3,3] matrix
Ds=sc*G*h*eye(2,2);

xg=xw(:,1);yg=xw(:,2);
[l23,l31,l12,y12,y31,y23,x12,x31,x23,Area,a4,a5,a6,b4,b5,b6,c4,c5,c6,d4,d5,d6,e4,e5,e6,C4,C5,C6,S4,S5,S6]=TrigElCoefsDKT(x,y,IEN);

% Calculation of shape functions and their derivatives at gauss points on
% the parent element
[ SF,DxsiSF,DetaSF,D2xsiSF,D2xsietaSF,D2etaSF ] = LNShapeFunDST(xg,yg);
[ SFm, DxsiSFm,DetaSFm] = LNShapeFunMassDST(xg,yg); 

%
%---> ax, ay, bx,by are constant for constant h (independednt of î,ç)
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

Fglob=zeros(GEN,1);

for kk=1:Nelem %for each element
  
% %% Mass matrix
    [ Hm,HW] = massHmDKT(kk,l23,l31,l12,C4,C5,C6,S4,S5,S6,GGin);
    %% Integration
    kb=zeros(9,9); 
    mloc=zeros(9,9);
    mloc1=zeros(9,9);
    mloc2=zeros(9,9);
    %************************** ADDITION
%         floc1=zeros(10,1);
    %***********************************
    [ Hxx,Hyy ] = rotationMass2(kk,a4,a5,a6,b4,b5,b6,c4,c5,c6,d4,d5,d6,e4,e5,e6 );
    [ Hxx1,Hyy1 ] = rotationMass(kk,l23,l31,l12,C4,C5,C6,S4,S5,S6);
    
    for ii=1:Ng % for each gauss point
    %% Stiffness
%        [Hx, Hy, Hx_xsi, Hx_eta,Hy_xsi, Hy_eta, Bb ] = ShapeFunDKT(ii,kk,a4,a5,a6,b4,b5,b6,c4,c5,c6,d4,d5,d6,e4,e5,e6,Area,SF,DxsiSF,DetaSF,y31,x31,y12,x12) ;
      [Hx, Hy, Hx_xsi,Hx_eta,Hy_xsi, Hy_eta, Bb ] = ShapeFunDKT2(ii,kk,C4,C5,C6,S4,S5,S6,Area, SF, DxsiSF,DetaSF,y31,x31,y12,x12,l23,l31,l12);
%       [Hx1, Hy1, Hx_xsi1, Hx_eta1,Hy_xsi1, Hy_eta1, Bb1 ] = ShapeFunDKT(ii,kk,a4,a5,a6,b4,b5,b6,c4,c5,c6,d4,d5,d6,e4,e5,e6,Area,SF,DxsiSF,DetaSF,y31,x31,y12,x12) ;
       
       [ Nm, HW,LW,L,HX3,HY3] = pseudoMassDKT(ii,kk,l23,l31,l12,y12,y31,x12,x31,Area,Hm,GGin,GGin2,xg, yg,IEN,x,y,Hx,Hy,SFm,C4,C5,C6,S4,S5,S6,Hxx,Hyy,HW )   ;
%        [ Nm,LW,L] = pseudoMassDKT(ii,kk,l23,l31,l12,y12,y31,x12,x31,Area,Hm,GGin,GGin2,xg, yg,IEN,x,y,Hx,Hy,SFm,C4,C5,C6,S4,S5,S6,HW)   ;

       kb=kb+Area(kk)*xw(ii,3)*(Bb'*BeSt*Bb);
       mloc=mloc+m*h*Area(kk)*xw(ii,3)*(Nm'*Nm);
%         mloc=mloc+m*h*Area(kk)*xw(ii,3)*((HW'*(LW'*LW)*HW)+h^2/12*(Hx'*Hx)+h^2/12*(Hy'*Hy));
       %mloc=mloc+m*h*Area(kk)*xw(ii,3)*((HW'*(LW'*LW)*HW));
    end
      
    kloc=kb;
    %************************** ADDITION
        floc=Area(kk)*P/3*[1 0 0 1 0 0 1 0 0]';% lumped mass approach for the uniform load
    %***********************************
    Mg(:,kk)=[mloc(:,1);mloc(:,2);mloc(:,3);mloc(:,4);mloc(:,5);mloc(:,6);mloc(:,7);mloc(:,8);mloc(:,9)];

    Kg(:,kk)=[kloc(:,1);kloc(:,2);kloc(:,3);kloc(:,4);kloc(:,5);kloc(:,6);kloc(:,7);kloc(:,8);kloc(:,9)];
    %************************** ADDITION
        for q=1:9
         Fglob(LM(q,kk))=Fglob(LM(q,kk))+floc(q);
        end 
    %***********************************
end

% error('er')
%% Global assembly (uses the LM)
iii=1:9; %Elnodes Number of element nodes
iii=repmat(iii',1,9);

rr=iii';
Ig=LM(iii(:),:);  %LM ARRAY (Dofs per Element)X(Elements)   -   (Nnodes*Ndof)X(Nel)
Jg=LM(rr(:),:);
Kglob=sparse(Ig(:),Jg(:),Kg(:),GEN,GEN);
Mglob=sparse(Ig(:),Jg(:),Mg(:),GEN,GEN);

% iii=1:9; %Elnodes Number of element nodes
% iii=repmat(iii',1,9);
% 
% rr=iii';
% Ig=LM(iii(:),:);  %LM ARRAY (Dofs per Element)X(Elements)   -   (Nnodes*Ndof)X(Nel)
% Jg=LM(rr(:),:);
% Kglob=sparse(Ig(:),Jg(:),Kg(:),GEN,GEN);
% Mglob=sparse(Ig(:),Jg(:),Mg(:),GEN,GEN);

if CC==1 %%
kkk=zeros(length(BBnodes),GEN);
mmm=zeros(length(BBnodes),GEN);
Dofs=length(p)*3;

% Totbound=reshape(BBound',1,4*size(BBound,2));
for j=1:length(BBnodes)
%     kkk(j,:)=[zeros(1,Bdofs1(j)-1) 1 zeros(1,Dofs-Bdofs1(j))];
    kkk(j,:)=[zeros(1,Bdofs(j)-1) 1 zeros(1,Dofs-Bdofs(j))];
end
Kglob=[Kglob kkk'; kkk zeros(length(BBnodes))];
Mglob=[Mglob mmm'; mmm zeros(length(BBnodes))];

end


% if CC==1 %%---> needs Work...debug the matrix augmentation procedure and consider mass matrix singularities
% kkk=zeros(length(BBnodes),GEN);
% mmm=zeros(length(BBnodes),GEN);
% Dofs=length(pp)*3;
% 
% % Totbound=reshape(BBound',1,4*size(BBound,2));
% for j=1:length(BBnodes)
%     kkk(j,:)=[zeros(1,Bdofs(j)-1) 1 zeros(1,Dofs-Bdofs(j))];
% end
% Kglob=[Kglob kkk'; kkk zeros(length(BBnodes))];
% Mglob=[Mglob mmm'; mmm zeros(length(BBnodes))];
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55

if lll==1
Fglob=zeros(length(Kglob),1);
% load TEST
Fglob(ID(1,PNODE))=P;
    U=Kglob\Fglob;
else
    Fglob1=[Fglob; zeros(length(BBnodes),1)];
    U=Kglob\Fglob1;
end

%--------------------------------------------
u=U(1:GEN); % 
w=u(1:3:end);

figure(1)
pdeplot(pp,ee,tt,'zdata',w,'colormap','copper')
shading interp

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


MAx=max(max(w))
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

DD=E*h^3/(12*(1-v^2));
cl=sqrt(m*h/DD);
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


