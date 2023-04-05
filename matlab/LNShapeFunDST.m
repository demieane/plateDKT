function [ SF, DxsiSF,DetaSF,D2xsiSF,D2xsietaSF,D2etaSF ] = LNShapeFunDST(xg,yg) 
%UNTITLED Summary of this function goes here
%   %% Calculation of quadratic Lagrange Shape Functions 
% and their derivatives at the reference triangle


SF=zeros(length(xg),6);
DxsiSF=zeros(length(xg),6);
DetaSF=zeros(length(xg),6);

for i=1:length(xg)
    
    SF(i,1)=(1-xg(i)-yg(i))*(1-2*xg(i)-2*yg(i)); %checked
    SF(i,2)=xg(i)*(2*xg(i)-1);%checked
    SF(i,3)=yg(i)*(2*yg(i)-1);
        
    SF(i,4)=4*xg(i)*yg(i);
    SF(i,5)=4*yg(i)*(1-xg(i)-yg(i));
    SF(i,6)=4*xg(i)*(1-xg(i)-yg(i));
    
    
    % 1st derivative --> î
    DxsiSF(i,1)=-3+4*(xg(i)+yg(i));%checked
    DxsiSF(i,2)=4*xg(i)-1;%checked
    DxsiSF(i,3)=0;%checked
    DxsiSF(i,4)=4*yg(i);%checked                  
    DxsiSF(i,5)=-4*yg(i);%checked
    DxsiSF(i,6)=4-8*xg(i)-4*yg(i);
    
    
    %1st derivative ç
    DetaSF(i,1)=-3+4*(xg(i)+yg(i));
    DetaSF(i,2)=0;
    DetaSF(i,3)=4*yg(i)-1;
    DetaSF(i,4)=4*xg(i);        
    DetaSF(i,5)=4-4*xg(i)-8*yg(i);
    DetaSF(i,6)=-4*xg(i); 

end

% 2nd derivative--> î    %checked
    D2xsiSF(1)=4;
    D2xsiSF(2)=4;
    D2xsiSF(3)=0;
    D2xsiSF(4)=0;
    D2xsiSF(5)=0;
    D2xsiSF(6)=-8;
    % mixed derivative--> î,ç symmetric!%checked
    D2xsietaSF(1)=4;
    D2xsietaSF(2)=0;
    D2xsietaSF(3)=0;
    D2xsietaSF(4)=4;
    D2xsietaSF(5)=-4;
    D2xsietaSF(6)=-4;
    
        %2nd derivative ç    %checked
    D2etaSF(1)=4;
    D2etaSF(2)=0;
    D2etaSF(3)=4;
    D2etaSF(4)=0;
    D2etaSF(5)=-8;
    D2etaSF(6)=0;
    
end