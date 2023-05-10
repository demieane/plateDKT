function [ ax,ay,bx,by] = CoefsDST( k,E,h,v,G,sc,D2xsiSF,D2xsietaSF,D2etaSF,b1,b2,c1,c2,phi)
% CoefsDST Function calculates the coefficients of Discrete Shear Triangle
% Coefficients correspond to the isotropic case 
% Second derivatives suggests that dependence on spatial coordinates is
% eliminated
   
ax=zeros(1,6);
ay=zeros(1,6);

bx=zeros(1,6);
by=zeros(1,6);

 AA=E.*h.^2./(12*(1-v^2)*sc*G);
%AA=E.*h.^3./(12*(1-v^2))*phi/10^2;

for i=1:6;  % ISOTROPIC CASE 
    % Cartesian derivatives of shape fucntions
    D2x=D2xsiSF(i)*b1(1,k).^2+2*D2xsietaSF(i)*b1(1,k)*b2(1,k)+D2etaSF(i)*b2(1,k).^2;
    D2y=D2xsiSF(i)*c1(1,k).^2+2*D2xsietaSF(i)*c1(1,k)*c2(1,k)+D2etaSF(i)*c2(1,k).^2;
    D2xy=D2xsiSF(i)*b1(1,k)*c1(1,k)+D2xsietaSF(i)*c2(1,k)*b1(1,k)+D2xsietaSF(i)*c1(1,k)*b2(1,k)+D2etaSF(i)*b2(1,k)*c2(1,k);
%     D2x=D2xsiSF(i);
%     D2y=D2etaSF(i);
%     D2xy=D2xsietaSF(i);
    % Assignment of coefficients (Batoz et al)
    ax(i)=AA*(D2x+(1-v)/2*D2y);     % FX
    ay(i)=AA*(v*D2xy+(1-v)/2*D2xy); %FX
    
    bx(i)=AA*(v*D2xy+(1-v)/2*D2xy);     %FY
    by(i)=AA*(D2y+(1-v)/2*D2x);         %FY

end


end

