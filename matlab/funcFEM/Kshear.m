function [ Bs, HFX, HFY ] = Kshear( Hs,asp,bsp,ayp,axp,byp,bxp )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

HFX(1)=asp(1)*Hs(1,1)+asp(2)*Hs(2,1)+asp(3)*Hs(3,1);
HFX(2)=asp(1)*Hs(1,2)+asp(2)*Hs(2,2)+asp(3)*Hs(3,2)+axp(1);
HFX(3)=asp(1)*Hs(1,3)+asp(2)*Hs(2,3)+asp(3)*Hs(3,3)+ayp(1);
HFX(4)=asp(1)*Hs(1,4)+asp(2)*Hs(2,4)+asp(3)*Hs(3,4);
HFX(5)=asp(1)*Hs(1,5)+asp(2)*Hs(2,5)+asp(3)*Hs(3,5)+axp(2);
HFX(6)=asp(1)*Hs(1,6)+asp(2)*Hs(2,6)+asp(3)*Hs(3,6)+ayp(2);
HFX(7)=asp(1)*Hs(1,7)+asp(2)*Hs(2,7)+asp(3)*Hs(3,7);
HFX(8)=asp(1)*Hs(1,8)+asp(2)*Hs(2,8)+asp(3)*Hs(3,8)+axp(3);
HFX(9)=asp(1)*Hs(1,9)+asp(2)*Hs(2,9)+asp(3)*Hs(3,9)+ayp(3);


HFY(1)=bsp(1)*Hs(1,1)+bsp(2)*Hs(2,1)+bsp(3)*Hs(3,1);
HFY(2)=bsp(1)*Hs(1,2)+bsp(2)*Hs(2,2)+bsp(3)*Hs(3,2)+bxp(1);
HFY(3)=bsp(1)*Hs(1,3)+bsp(2)*Hs(2,3)+bsp(3)*Hs(3,3)+byp(1);
HFY(4)=bsp(1)*Hs(1,4)+bsp(2)*Hs(2,4)+bsp(3)*Hs(3,4);
HFY(5)=bsp(1)*Hs(1,5)+bsp(2)*Hs(2,5)+bsp(3)*Hs(3,5)+bxp(2);
HFY(6)=bsp(1)*Hs(1,6)+bsp(2)*Hs(2,6)+bsp(3)*Hs(3,6)+byp(2);
HFY(7)=bsp(1)*Hs(1,7)+bsp(2)*Hs(2,7)+bsp(3)*Hs(3,7);
HFY(8)=bsp(1)*Hs(1,8)+bsp(2)*Hs(2,8)+bsp(3)*Hs(3,8)+bxp(3);
HFY(9)=bsp(1)*Hs(1,9)+bsp(2)*Hs(2,9)+bsp(3)*Hs(3,9)+byp(3);

Bs=[HFX;HFY];
end

