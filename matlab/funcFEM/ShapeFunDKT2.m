function [Hx, Hy, Hx_xsi, Hx_eta,Hy_xsi, Hy_eta, Bb ] = ShapeFunDKT2(ii,k,C4,C5,C6,S4,S5,S6,Area, SF, DxsiSF,DetaSF,y31,x31,y12,x12,l23,l31,l12)


%--- иy ич form


S4=S4(1,k);
S5=S5(1,k);
S6=S6(1,k);
C4=C4(1,k);
C5=C5(1,k);
C6=C6(1,k);

Hx(1)=S5*3/(2*l31(k))*SF(ii,5)-S6*3/(2*l12(k))*SF(ii,6);
Hx(2)=-SF(ii,5)*(3/4*S5*C5)-SF(ii,6)*(3/4*S6*C6);%
Hx(3)=SF(ii,1)+(0.5*C5^2-0.25*S5^2)*SF(ii,5)+(0.5*C6^2-0.25*S6^2)*SF(ii,6);%

Hx(4)=S6*3/(2*l12(k))*SF(ii,6)-S4*3/(2*l23(k))*SF(ii,4);
Hx(5)=-SF(ii,4)*(3/4*S4*C4)-SF(ii,6)*(3/4*S6*C6);%
Hx(6)=SF(ii,2)+(0.5*C4^2-0.25*S4^2)*SF(ii,4)+(0.5*C6^2-0.25*S6^2)*SF(ii,6);%

Hx(7)=S4*3/(2*l23(k))*SF(ii,4)-S5*3/(2*l31(k))*SF(ii,5);
Hx(8)=-SF(ii,4)*(3/4*S4*C4)-SF(ii,5)*(3/4*S5*C5);%
Hx(9)=SF(ii,3)+(0.5*C4^2-0.25*S4^2)*SF(ii,4)+(0.5*C5^2-0.25*S5^2)*SF(ii,5);%


%%%
Hx_xsi(1)=S5*3/(2*l31(k))*DxsiSF(ii,5)-S6*3/(2*l12(k))*DxsiSF(ii,6);
Hx_xsi(2)=-DxsiSF(ii,5)*(3/4*S5*C5)-DxsiSF(ii,6)*(3/4*S6*C6);%
Hx_xsi(3)=DxsiSF(ii,1)+(0.5*C5^2-0.25*S5^2)*DxsiSF(ii,5)+(0.5*C6^2-0.25*S6^2)*DxsiSF(ii,6);%

Hx_xsi(4)=S6*3/(2*l12(k))*DxsiSF(ii,6)-S4*3/(2*l23(k))*DxsiSF(ii,4);
Hx_xsi(5)=-DxsiSF(ii,4)*(3/4*S4*C4)-DxsiSF(ii,6)*(3/4*S6*C6);%
Hx_xsi(6)=DxsiSF(ii,2)+(0.5*C4^2-0.25*S4^2)*DxsiSF(ii,4)+(0.5*C6^2-0.25*S6^2)*DxsiSF(ii,6);%

Hx_xsi(7)=S4*3/(2*l23(k))*DxsiSF(ii,4)-S5*3/(2*l31(k))*DxsiSF(ii,5);
Hx_xsi(8)=-DxsiSF(ii,4)*(3/4*S4*C4)-DxsiSF(ii,5)*(3/4*S5*C5);%
Hx_xsi(9)=DxsiSF(ii,3)+(0.5*C4^2-0.25*S4^2)*DxsiSF(ii,4)+(0.5*C5^2-0.25*S5^2)*DxsiSF(ii,5);%



%%%
Hx_eta(1)=S5*3/(2*l31(k))*DetaSF(ii,5)-S6*3/(2*l12(k))*DetaSF(ii,6);
Hx_eta(2)=-DetaSF(ii,5)*(3/4*S5*C5)-DetaSF(ii,6)*(3/4*S6*C6);%
Hx_eta(3)=DetaSF(ii,1)+(0.5*C5^2-0.25*S5^2)*DetaSF(ii,5)+(0.5*C6^2-0.25*S6^2)*DetaSF(ii,6);%

Hx_eta(4)=S6*3/(2*l12(k))*DetaSF(ii,6)-S4*3/(2*l23(k))*DetaSF(ii,4);
Hx_eta(5)=-DetaSF(ii,4)*(3/4*S4*C4)-DetaSF(ii,6)*(3/4*S6*C6);%
Hx_eta(6)=DetaSF(ii,2)+(0.5*C4^2-0.25*S4^2)*DetaSF(ii,4)+(0.5*C6^2-0.25*S6^2)*DetaSF(ii,6);%

Hx_eta(7)=S4*3/(2*l23(k))*DetaSF(ii,4)-S5*3/(2*l31(k))*DetaSF(ii,5);
Hx_eta(8)=-DetaSF(ii,4)*(3/4*S4*C4)-DetaSF(ii,5)*(3/4*S5*C5);%
Hx_eta(9)=DetaSF(ii,3)+(0.5*C4^2-0.25*S4^2)*DetaSF(ii,4)+(0.5*C5^2-0.25*S5^2)*DetaSF(ii,5);%



%%
Hy(1)=-C5*3/(2*l31(k))*SF(ii,5)+C6*3/(2*l12(k))*SF(ii,6);
Hy(2)=-SF(ii,1)-(0.5*S5^2-0.25*C5^2)*SF(ii,5)-(0.5*S6^2-0.25*C6^2)*SF(ii,6);%
Hy(3)=SF(ii,5)*(3/4*S5*C5)+SF(ii,6)*(3/4*S6*C6);%

Hy(4)=-C6*3/(2*l12(k))*SF(ii,6)+C4*3/(2*l23(k))*SF(ii,4);
Hy(5)=-SF(ii,2)-(0.5*S4^2-0.25*C4^2)*SF(ii,4)-(0.5*S6^2-0.25*C6^2)*SF(ii,6);%
Hy(6)=SF(ii,4)*(3/4*S4*C4)+SF(ii,6)*(3/4*S6*C6);%

Hy(7)=-C4*3/(2*l23(k))*SF(ii,4)+C5*3/(2*l31(k))*SF(ii,5);
Hy(8)=-SF(ii,3)-(0.5*S4^2-0.25*C4^2)*SF(ii,4)-(0.5*S5^2-0.25*C5^2)*SF(ii,5);%
Hy(9)=SF(ii,4)*(3/4*S4*C4)+SF(ii,5)*(3/4*S5*C5);%


%%
Hy_xsi(1)=-C5*3/(2*l31(k))*DxsiSF(ii,5)+C6*3/(2*l12(k))*DxsiSF(ii,6);
Hy_xsi(2)=-DxsiSF(ii,1)-(0.5*S5^2-0.25*C5^2)*DxsiSF(ii,5)-(0.5*S6^2-0.25*C6^2)*DxsiSF(ii,6);%
Hy_xsi(3)=DxsiSF(ii,5)*(3/4*S5*C5)+DxsiSF(ii,6)*(3/4*S6*C6);%

Hy_xsi(4)=-C6*3/(2*l12(k))*DxsiSF(ii,6)+C4*3/(2*l23(k))*DxsiSF(ii,4);
Hy_xsi(5)=-DxsiSF(ii,2)-(0.5*S4^2-0.25*C4^2)*DxsiSF(ii,4)-(0.5*S6^2-0.25*C6^2)*DxsiSF(ii,6);%
Hy_xsi(6)=DxsiSF(ii,4)*(3/4*S4*C4)+DxsiSF(ii,6)*(3/4*S6*C6);%

Hy_xsi(7)=-C4*3/(2*l23(k))*DxsiSF(ii,4)+C5*3/(2*l31(k))*DxsiSF(ii,5);
Hy_xsi(8)=-DxsiSF(ii,3)-(0.5*S4^2-0.25*C4^2)*DxsiSF(ii,4)-(0.5*S5^2-0.25*C5^2)*DxsiSF(ii,5);%
Hy_xsi(9)=DxsiSF(ii,4)*(3/4*S4*C4)+DxsiSF(ii,5)*(3/4*S5*C5);%





%%
Hy_eta(1)=-C5*3/(2*l31(k))*DetaSF(ii,5)+C6*3/(2*l12(k))*DetaSF(ii,6);
Hy_eta(2)=-DetaSF(ii,1)-(0.5*S5^2-0.25*C5^2)*DetaSF(ii,5)-(0.5*S6^2-0.25*C6^2)*DetaSF(ii,6);%
Hy_eta(3)=DetaSF(ii,5)*(3/4*S5*C5)+DetaSF(ii,6)*(3/4*S6*C6);%

Hy_eta(4)=-C6*3/(2*l12(k))*DetaSF(ii,6)+C4*3/(2*l23(k))*DetaSF(ii,4);
Hy_eta(5)=-DetaSF(ii,2)-(0.5*S4^2-0.25*C4^2)*DetaSF(ii,4)-(0.5*S6^2-0.25*C6^2)*DetaSF(ii,6);%
Hy_eta(6)=DetaSF(ii,4)*(3/4*S4*C4)+DetaSF(ii,6)*(3/4*S6*C6);%

Hy_eta(7)=-C4*3/(2*l23(k))*DetaSF(ii,4)+C5*3/(2*l31(k))*DetaSF(ii,5);
Hy_eta(8)=-DetaSF(ii,3)-(0.5*S4^2-0.25*C4^2)*DetaSF(ii,4)-(0.5*S5^2-0.25*C5^2)*DetaSF(ii,5);%
Hy_eta(9)=DetaSF(ii,4)*(3/4*S4*C4)+DetaSF(ii,5)*(3/4*S5*C5);%


Bb(1,:)=1/(2*Area(k))*(y31(k)*Hx_xsi+y12(k)*Hx_eta);
Bb(2,:)=1/(2*Area(k))*(-x31(k)*Hy_xsi-x12(k)*Hy_eta);
Bb(3,:)=1/(2*Area(k))*(-x31(k)*Hx_xsi-x12(k)*Hx_eta+y31(k)*Hy_xsi+y12(k)*Hy_eta);
end