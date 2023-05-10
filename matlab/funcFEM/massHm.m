function [ Hm ] = massHm( k,HFX, HFY,l23,l31,l12,C4,C5,C6,S4,S5,S6 )
%UNTITLED4 Summary of this function goes here
%   {W} = Hm {u} --> 10 x1 = [10 x 9] x [9 x 1]



A12=4/27*l12(k)-2/27*l12(k);
B12=2/27*l12(k)-4/27*l12(k);

A23=4/27*l23(k)-2/27*l23(k);
B23=2/27*l23(k)-4/27*l23(k);

A31=4/27*l31(k)-2/27*l31(k);
B31=2/27*l31(k)-4/27*l31(k);

Hm(1,:)=[1 zeros(1,8)];
Hm(2,:)=[zeros(1,3) 1 zeros(1,5) ];
Hm(3,:)=[zeros(1,6) 1 zeros(1,2) ];


Hm(4,:)=[1/3 zeros(1,2) 1/3 zeros(1,2) 1/3 zeros(1,2) ];
%% 2-3

Hm(5,:)=[A23*(-S4(1,k)*HFX(1)+C4(1,k)*HFY(1)) A23*(-S4(1,k)*HFX(2)+C4(1,k)*HFY(2)) A23*(-S4(1,k)*HFX(3)+C4(1,k)*HFY(3)) ...
         A23*(-S4(1,k)*HFX(4)+C4(1,k)*HFY(4))+20/27 A23*(-S4(1,k)*HFX(5)+C4(1,k)*HFY(5))+4*l23(k)/27*S4(1,k) A23*(-S4(1,k)*HFX(6)+C4(1,k)*HFY(6))-4*l23(k)/27*C4(1,k) ...
         A23*(-S4(1,k)*HFX(7)+C4(1,k)*HFY(7))+7/27 A23*(-S4(1,k)*HFX(8)+C4(1,k)*HFY(8))-2*l23(k)/27*S4(1,k) A23*(-S4(1,k)*HFX(9)+C4(1,k)*HFY(9))+2*l23(k)/27*C4(1,k)];

Hm(6,:)=[B23*(-S4(1,k)*HFX(1)+C4(1,k)*HFY(1)) B23*(-S4(1,k)*HFX(2)+C4(1,k)*HFY(2)) B23*(-S4(1,k)*HFX(3)+C4(1,k)*HFY(3)) ...
         B23*(-S4(1,k)*HFX(4)+C4(1,k)*HFY(4))+7/27 B23*(-S4(1,k)*HFX(5)+C4(1,k)*HFY(5))+2*l23(k)/27*S4(1,k) B23*(-S4(1,k)*HFX(6)+C4(1,k)*HFY(6))-2*l23(k)/27*C4(1,k) ...
         B23*(-S4(1,k)*HFX(7)+C4(1,k)*HFY(7))+20/27 B23*(-S4(1,k)*HFX(8)+C4(1,k)*HFY(8))-4*l23(k)/27*S4(1,k) B23*(-S4(1,k)*HFX(9)+C4(1,k)*HFY(9))+4*l23(k)/27*C4(1,k)];

%% 3-1

Hm(7,:)=[A31*(-S5(1,k)*HFX(1)+C5(1,k)*HFY(1))+7/27 A31*(-S5(1,k)*HFX(2)+C5(1,k)*HFY(2))-2*l31(k)/27*S5(1,k) A31*(-S5(1,k)*HFX(3)+C5(1,k)*HFY(3))+2*l31(k)/27*C5(1,k) ...
         A31*(-S5(1,k)*HFX(4)+C5(1,k)*HFY(4)) A31*(-S5(1,k)*HFX(5)+C5(1,k)*HFY(5)) A31*(-S5(1,k)*HFX(6)+C5(1,k)*HFY(6)) ...
         A31*(-S5(1,k)*HFX(7)+C5(1,k)*HFY(7))+20/27 A31*(-S5(1,k)*HFX(8)+C5(1,k)*HFY(8))+4*l31(k)/27*S5(1,k) A31*(-S5(1,k)*HFX(9)+C5(1,k)*HFY(9))-4*l31(k)/27*C5(1,k)];

Hm(8,:)=[B31*(-S5(1,k)*HFX(1)+C5(1,k)*HFY(1))+20/27 B31*(-S5(1,k)*HFX(2)+C5(1,k)*HFY(2))-4*l31(k)/27*S5(1,k) B31*(-S5(1,k)*HFX(3)+C5(1,k)*HFY(3))+4*l31(k)/27*C5(1,k) ...
         B31*(-S5(1,k)*HFX(4)+C5(1,k)*HFY(4)) B31*(-S5(1,k)*HFX(5)+C5(1,k)*HFY(5)) B31*(-S5(1,k)*HFX(6)+C5(1,k)*HFY(6)) ...
         B31*(-S5(1,k)*HFX(7)+C5(1,k)*HFY(7))+7/27 B31*(-S5(1,k)*HFX(8)+C5(1,k)*HFY(8))+2*l31(k)/27*S5(1,k) B31*(-S5(1,k)*HFX(9)+C5(1,k)*HFY(9))-2*l31(k)/27*C5(1,k)];
%% 1-2

Hm(9,:)=[A12*(-S6(1,k)*HFX(1)+C6(1,k)*HFY(1))+20/27 A12*(-S6(1,k)*HFX(2)+C6(1,k)*HFY(2))+4*l12(k)/27*S6(1,k) A12*(-S6(1,k)*HFX(3)+C6(1,k)*HFY(3))-4*l12(k)/27*C6(1,k) ...
         A12*(-S6(1,k)*HFX(4)+C6(1,k)*HFY(4))+7/27 A12*(-S6(1,k)*HFX(5)+C6(1,k)*HFY(5))-2*l12(k)/27*S6(1,k) A12*(-S6(1,k)*HFX(6)+C6(1,k)*HFY(6))+2*l12(k)/27*C6(1,k) ...
         A12*(-S6(1,k)*HFX(7)+C6(1,k)*HFY(7)) A12*(-S6(1,k)*HFX(8)+C6(1,k)*HFY(8)) A12*(-S6(1,k)*HFX(9)+C6(1,k)*HFY(9))];

     
Hm(10,:)=[B12*(-S6(1,k)*HFX(1)+C6(1,k)*HFY(1))+7/27 B12*(-S6(1,k)*HFX(2)+C6(1,k)*HFY(2))+2*l12(k)/27*S6(1,k) B12*(-S5(1,k)*HFX(3)+C6(1,k)*HFY(3))-2*l12(k)/27*C6(1,k) ...
         B12*(-S6(1,k)*HFX(4)+C6(1,k)*HFY(4))+20/27 B12*(-S6(1,k)*HFX(5)+C6(1,k)*HFY(5))-4*l12(k)/27*S6(1,k) B12*(-S6(1,k)*HFX(6)+C6(1,k)*HFY(6))+4*l12(k)/27*C6(1,k) ...
         B12*(-S6(1,k)*HFX(7)+C6(1,k)*HFY(7)) B12*(-S6(1,k)*HFX(8)+C6(1,k)*HFY(8)) B12*(-S6(1,k)*HFX(9)+C6(1,k)*HFY(9))];




end
