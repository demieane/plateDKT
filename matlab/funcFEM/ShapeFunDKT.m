function [Hx, Hy, Hx_xsi, Hx_eta,Hy_xsi, Hy_eta, Bb ] = ShapeFunDKT(ii,kk,a4,a5,a6,b4,b5,b6,c4,c5,c6,d4,d5,d6,e4,e5,e6,Area, SF, DxsiSF,DetaSF,y31,x31,y12,x12) ;


%--- иy ич form
Hx(1)=1.5*(a6(kk)*SF(ii,6)-a5(kk)*SF(ii,5));%
Hx(2)=b5(kk)*SF(ii,5)+b6(kk)*SF(ii,6);%
Hx(3)=SF(ii,1)-c5(kk)*SF(ii,5)-c6(kk)*SF(ii,6);%

Hx(4)=1.5*(a4(kk)*SF(ii,4)-a6(kk)*SF(ii,6));%
Hx(5)=b6(kk)*SF(ii,6)+b4(kk)*SF(ii,4);%
Hx(6)=SF(ii,2)-c6(kk)*SF(ii,6)-c4(kk)*SF(ii,4);%

Hx(7)=1.5*(a5(kk)*SF(ii,5)-a4(kk)*SF(ii,4));%
Hx(8)=b4(kk)*SF(ii,4)+b5(kk)*SF(ii,5);%
Hx(9)=SF(ii,3)-c4(kk)*SF(ii,4)-c5(kk)*SF(ii,5);%


%%%
Hx_xsi(1)=1.5*(a6(kk)*DxsiSF(ii,6)-a5(kk)*DxsiSF(ii,5));%
Hx_xsi(2)=b5(kk)*DxsiSF(ii,5)+b6(kk)*DxsiSF(ii,6);%
Hx_xsi(3)=DxsiSF(ii,1)-c5(kk)*DxsiSF(ii,5)-c6(kk)*DxsiSF(ii,6);%

Hx_xsi(4)=1.5*(a4(kk)*DxsiSF(ii,4)-a6(kk)*DxsiSF(ii,6));%
Hx_xsi(5)=b6(kk)*DxsiSF(ii,6)+b4(kk)*DxsiSF(ii,4);%
Hx_xsi(6)=DxsiSF(ii,2)-c6(kk)*DxsiSF(ii,6)-c4(kk)*DxsiSF(ii,4);%

Hx_xsi(7)=1.5*(a5(kk)*DxsiSF(ii,5)-a4(kk)*DxsiSF(ii,4));%
Hx_xsi(8)=b4(kk)*DxsiSF(ii,4)+b5(kk)*DxsiSF(ii,5);%
Hx_xsi(9)=DxsiSF(ii,3)-c4(kk)*DxsiSF(ii,4)-c5(kk)*DxsiSF(ii,5);%


%%%
Hx_eta(1)=1.5*(a6(kk)*DetaSF(ii,6)-a5(kk)*DetaSF(ii,5));%
Hx_eta(2)=b5(kk)*DetaSF(ii,5)+b6(kk)*DetaSF(ii,6);%
Hx_eta(3)=DetaSF(ii,1)-c5(kk)*DetaSF(ii,5)-c6(kk)*DetaSF(ii,6);%

Hx_eta(4)=1.5*(a4(kk)*DetaSF(ii,4)-a6(kk)*DetaSF(ii,6));%
Hx_eta(5)=b6(kk)*DetaSF(ii,6)+b4(kk)*DetaSF(ii,4);%
Hx_eta(6)=DetaSF(ii,2)-c6(kk)*DetaSF(ii,6)-c4(kk)*DetaSF(ii,4);%

Hx_eta(7)=1.5*(a5(kk)*DetaSF(ii,5)-a4(kk)*DetaSF(ii,4));%
Hx_eta(8)=b4(kk)*DetaSF(ii,4)+b5(kk)*DetaSF(ii,5);%
Hx_eta(9)=DetaSF(ii,3)-c4(kk)*DetaSF(ii,4)-c5(kk)*DetaSF(ii,5);%




%%
Hy(1)=1.5*(d6(kk)*SF(ii,6)-d5(kk)*SF(ii,5));%
Hy(2)=-SF(ii,1)+e5(kk)*SF(ii,5)+e6(kk)*SF(ii,6);%
Hy(3)=-(b5(kk)*SF(ii,5)+b6(kk)*SF(ii,6));%

Hy(4)=1.5*(d4(kk)*SF(ii,4)-d6(kk)*SF(ii,6));%
Hy(5)=-SF(ii,2)+e4(kk)*SF(ii,4)+e6(kk)*SF(ii,6);%
Hy(6)=-(b4(kk)*SF(ii,4)+b6(kk)*SF(ii,6));%

Hy(7)=1.5*(d5(kk)*SF(ii,5)-d4(kk)*SF(ii,4));%
Hy(8)=-SF(ii,3)+e4(kk)*SF(ii,4)+e5(kk)*SF(ii,5);%
Hy(9)=-(b4(kk)*SF(ii,4)+b5(kk)*SF(ii,5));%
%%
Hy_xsi(1)=1.5*(d6(kk)*DxsiSF(ii,6)-d5(kk)*DxsiSF(ii,5));%
Hy_xsi(2)=-DxsiSF(ii,1)+e5(kk)*DxsiSF(ii,5)+e6(kk)*DxsiSF(ii,6);%
Hy_xsi(3)=-(b5(kk)*DxsiSF(ii,5)+b6(kk)*DxsiSF(ii,6));%

Hy_xsi(4)=1.5*(d4(kk)*DxsiSF(ii,4)-d6(kk)*DxsiSF(ii,6));%
Hy_xsi(5)=-DxsiSF(ii,2)+e4(kk)*DxsiSF(ii,4)+e6(kk)*DxsiSF(ii,6);%
Hy_xsi(6)=-(b4(kk)*DxsiSF(ii,4)+b6(kk)*DxsiSF(ii,6));%

Hy_xsi(7)=1.5*(d5(kk)*DxsiSF(ii,5)-d4(kk)*DxsiSF(ii,4));%
Hy_xsi(8)=-DxsiSF(ii,3)+e4(kk)*DxsiSF(ii,4)+e5(kk)*DxsiSF(ii,5);%
Hy_xsi(9)=-(b4(kk)*DxsiSF(ii,4)+b5(kk)*DxsiSF(ii,5));%



%%
Hy_eta(1)=1.5*(d6(kk)*DetaSF(ii,6)-d5(kk)*DetaSF(ii,5));%
Hy_eta(2)=-DetaSF(ii,1)+e5(kk)*DetaSF(ii,5)+e6(kk)*DetaSF(ii,6);%
Hy_eta(3)=-(b5(kk)*DetaSF(ii,5)+b6(kk)*DetaSF(ii,6));%

Hy_eta(4)=1.5*(d4(kk)*DetaSF(ii,4)-d6(kk)*DetaSF(ii,6));%
Hy_eta(5)=-DetaSF(ii,2)+e4(kk)*DetaSF(ii,4)+e6(kk)*DetaSF(ii,6);%
Hy_eta(6)=-(b4(kk)*DetaSF(ii,4)+b6(kk)*DetaSF(ii,6));%

Hy_eta(7)=1.5*(d5(kk)*DetaSF(ii,5)-d4(kk)*DetaSF(ii,4));%
Hy_eta(8)=-DetaSF(ii,3)+e4(kk)*DetaSF(ii,4)+e5(kk)*DetaSF(ii,5);
Hy_eta(9)=-(b4(kk)*DetaSF(ii,4)+b5(kk)*DetaSF(ii,5));%


Bb(1,:)=1/(2*Area(kk))*(y31(kk)*Hx_xsi+y12(kk)*Hx_eta);
Bb(2,:)=1/(2*Area(kk))*(-x31(kk)*Hy_xsi-x12(kk)*Hy_eta);
Bb(3,:)=1/(2*Area(kk))*(-x31(kk)*Hx_xsi-x12(kk)*Hx_eta+y31(kk)*Hy_xsi+y12(kk)*Hy_eta);
end