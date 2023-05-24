function [Fglob1] = createFglob(lll,GEN,Nelem,P_load,Fx,Area,LM,BBnodes)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

Fglob=zeros(GEN,1);
% GEN
% Nelem

% Fx_new = Fx.*sin(wf*t);
Fx_new = Fx;

for kk=1:Nelem %for each element (iS THIS TRIANGLE 1 IN t?)
  
% %% Mass matrix
%     [ Hm,HW] = massHmDKT(kk,l23,l31,l12,C4,C5,C6,S4,S5,S6,GGin);
    %% Integration
%     kb=zeros(9,9); 
%     mloc=zeros(9,9);
%     mloc1=zeros(9,9);
%     mloc2=zeros(9,9);
    %************************** ADDITION
%         floc1=zeros(10,1);
    %***********************************
%     [ Hxx,Hyy ] = rotationMass2(kk,a4,a5,a6,b4,b5,b6,c4,c5,c6,d4,d5,d6,e4,e5,e6 );
%     [ Hxx1,Hyy1 ] = rotationMass(kk,l23,l31,l12,C4,C5,C6,S4,S5,S6);
%     
%     for ii=1:Ng % for each gauss point
    %% Stiffness
%        [Hx, Hy, Hx_xsi, Hx_eta,Hy_xsi, Hy_eta, Bb ] = ShapeFunDKT(ii,kk,a4,a5,a6,b4,b5,b6,c4,c5,c6,d4,d5,d6,e4,e5,e6,Area,SF,DxsiSF,DetaSF,y31,x31,y12,x12) ;
% %       [Hx, Hy, Hx_xsi,Hx_eta,Hy_xsi, Hy_eta, Bb ] = ShapeFunDKT2(ii,kk,C4,C5,C6,S4,S5,S6,Area, SF, DxsiSF,DetaSF,y31,x31,y12,x12,l23,l31,l12);
%       [Hx1, Hy1, Hx_xsi1, Hx_eta1,Hy_xsi1, Hy_eta1, Bb1 ] = ShapeFunDKT(ii,kk,a4,a5,a6,b4,b5,b6,c4,c5,c6,d4,d5,d6,e4,e5,e6,Area,SF,DxsiSF,DetaSF,y31,x31,y12,x12) ;
       
% %        [ Nm, HW,LW,L,HX3,HY3] = pseudoMassDKT(ii,kk,l23,l31,l12,y12,y31,x12,x31,Area,Hm,GGin,GGin2,xg, yg,IEN,x,y,Hx,Hy,SFm,C4,C5,C6,S4,S5,S6,Hxx,Hyy,HW )   ;
% %        [ Nm,LW,L] = pseudoMassDKT(ii,kk,l23,l31,l12,y12,y31,x12,x31,Area,Hm,GGin,GGin2,xg, yg,IEN,x,y,Hx,Hy,SFm,C4,C5,C6,S4,S5,S6,HW)   ;

       % we give one more dimension 
       % to the BeSt [3x3] original matrix
% %        kb=kb+Area(kk)*xw(ii,3)*(Bb'*BeSt2(:,:,kk)*Bb);
%        kb=kb+Area(kk)*xw(ii,3)*(Bb'*BeSt*Bb);
%        mloc=mloc+m*h*Area(kk)*xw(ii,3)*(Nm'*Nm);
% %        mloc=mloc+m*txxBEM(kk)*Area(kk)*xw(ii,3)*(Nm'*Nm);
       %****************************ADDITION
%        floc1=floc1+Area(kk)*xw(ii,3)*(LW');
       %************************************
%         mloc=mloc+m*txxBEM(kk)*Area(kk)*xw(ii,3)*((HW'*(LW'*LW)*HW)+txxBEM(kk)^2/12*(Hx'*Hx)+txxBEM(kk)^2/12*(Hy'*Hy));
%        mloc=mloc+m*txxBEM(kk)*Area(kk)*xw(ii,3)*((HW'*(LW'*LW)*HW));
%     end
      
%     kloc=kb;
    %************************** ADDITION
%         floc=P*HW'*floc1;
    if lll==2 % uniform load
        floc=Area(kk)*P_load/3*[1 0 0 1 0 0 1 0 0]';% lumped mass approach for the uniform load
    elseif lll==3 %distributed load via mapping func
        floc=Area(kk)*Fx_new(kk)/3*[1 0 0 1 0 0 1 0 0]';
    end
%     floc
    %***********************************
%     Mg(:,kk)=[mloc(:,1);mloc(:,2);mloc(:,3);mloc(:,4);mloc(:,5);mloc(:,6);mloc(:,7);mloc(:,8);mloc(:,9)];
% 
%     Kg(:,kk)=[kloc(:,1);kloc(:,2);kloc(:,3);kloc(:,4);kloc(:,5);kloc(:,6);kloc(:,7);kloc(:,8);kloc(:,9)];
    %************************** ADDITION
    for q=1:9
%         LM(q,kk)
         Fglob(LM(q,kk))=Fglob(LM(q,kk))+floc(q);
    end 
    %***********************************
end


Fglob1=[Fglob; zeros(length(BBnodes),1)];
% size(Fglob1)

end

