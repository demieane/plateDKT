function [SF,DxsiSF,DetaSF]=SHAPE_FUN_LANTRIG(EL_TYPE,xg,yg)

% Function computes the shape functions and their partial derivatives
% at the parent element using barycetric coords xsi, eta
% Supports -Linear,-Quadratic and -Cubic Langrangian Triangles

if EL_TYPE==1;
    Nnodes=3;
elseif  EL_TYPE==2;
    Nnodes=6;
elseif  EL_TYPE==3;
    Nnodes=10;
end

SF=zeros(length(xg),Nnodes);
DxsiSF=zeros(length(xg),Nnodes);
DetaSF=zeros(length(xg),Nnodes);

%______________Linear Triangle
if Nnodes==3
    for i=1:length(xg);
    
    SF(i,1)=1-xg(i)-yg(i);
    SF(i,2)=xg(i);
    SF(i,3)=yg(i);
    
    DxsiSF(i,1)=-1;
    DxsiSF(i,2)=1;
    DxsiSF(i,3)=0;
    
    DetaSF(i,1)=-1;
    DetaSF(i,2)=0;
    DetaSF(i,3)=1;
    
    end
end
%______________Quadratic Trinagle
if Nnodes==6   
    for i=1:length(xg)
        
    SF(i,1)=(1-xg(i)-yg(i))*(1-2*xg(i)-2*yg(i));
    SF(i,2)=xg(i)*(2*xg(i)-1);
    SF(i,3)=yg(i)*(2*yg(i)-1);
    SF(i,4)=4*xg(i)*(1-xg(i)-yg(i));
    SF(i,5)=4*xg(i)*yg(i);
    SF(i,6)=4*yg(i)*(1-xg(i)-yg(i));
    
    DxsiSF(i,1)=-3+4*(xg(i)+yg(i));
    DxsiSF(i,2)=4*xg(i)-1;
    DxsiSF(i,3)=0;
    DxsiSF(i,4)=4-8*xg(i)-4*yg(i);
    DxsiSF(i,5)=4*yg(i);
    DxsiSF(i,6)=-4*yg(i);
    
    DetaSF(i,1)=-3+4*(xg(i)+yg(i));
    DetaSF(i,2)=0;
    DetaSF(i,3)=4*yg(i)-1;
    DetaSF(i,4)=-4*xg(i);
    DetaSF(i,5)=4*xg(i);
    DetaSF(i,6)=4-4*xg(i)-8*yg(i);
    
    end
end
%______________Qubic Triangle
if Nnodes==10
    for i=1:length(xg)
        
    SF(i,1)=9/2*(1-xg(i)-yg(i))*(1/3-xg(i)-yg(i))*(2/3-xg(i)-yg(i));
    SF(i,2)=9/2*xg(i)*(xg(i)-1/3)*(xg(i)-2/3);
    SF(i,3)=9/2*yg(i)*(yg(i)-1/3)*(yg(i)-2/3);
    SF(i,4)=27/2*(1-xg(i)-yg(i))*xg(i)*(2/3-xg(i)-yg(i));
    SF(i,5)=27/2*(1-xg(i)-yg(i))*xg(i)*(xg(i)-1/3);
    SF(i,6)=27/2*xg(i)*yg(i)*(xg(i)-1/3);    
    SF(i,7)=27/2*xg(i)*yg(i)*(yg(i)-1/3);  
    SF(i,8)=27/2*(1-xg(i)-yg(i))*yg(i)*(yg(i)-1/3);  
    SF(i,9)=27/2*(1-xg(i)-yg(i))*yg(i)*(2/3-xg(i)-yg(i)); 
    SF(i,10)=27*xg(i)*yg(i)*(1-xg(i)-yg(i)); 
    
    DxsiSF(i,1)=1/2*(-11-27*xg(i)^2+xg(i)*(36-54*yg(i))+36*yg(i)-27*yg(i)^2);
    DxsiSF(i,2)=1-9*xg(i)+27*xg(i)^2/2;
    DxsiSF(i,3)=0;
    DxsiSF(i,4)=9/2*(2+9*xg(i)^2-5*yg(i)+3*yg(i)^2+2*xg(i)*(-5+6*yg(i)));
    DxsiSF(i,5)=-9/2*(1+9*xg(i)^2-yg(i)+xg(i)*(-8+6*yg(i)));
    DxsiSF(i,6)=9/2*(-1+6*xg(i))*yg(i);  
    DxsiSF(i,7)=9/2*(-1+3*yg(i))*yg(i);  
    DxsiSF(i,8)=-9/2*yg(i)*(-1+3*yg(i));
    DxsiSF(i,9)=9/2*yg(i)*(-5+6*xg(i)+6*yg(i));
    DxsiSF(i,10)=-27*yg(i)*(-1+2*xg(i)+yg(i));
    
    DetaSF(i,1)=1/2*(-11-27*xg(i)^2+xg(i)*(36-54*yg(i))+36*yg(i)-27*yg(i)^2);
    DetaSF(i,2)=0;
    DetaSF(i,3)=1-9*yg(i)+27/2*yg(i)^2;
    DetaSF(i,4)=9/2*xg(i)*(-5+6*xg(i)+6*yg(i));
    DetaSF(i,5)=-9/2*xg(i)*(-1+3*xg(i));
    DetaSF(i,6)=9/2*xg(i)*(-1+3*xg(i));  
    DetaSF(i,7)=9/2*xg(i)*(-1+6*yg(i));   
    DetaSF(i,8)=-9/2*(1-8*yg(i)+9*yg(i)^2+xg(i)*(-1+6*yg(i)));
    DetaSF(i,9)=9/2*(2+3*xg(i)^2-10*yg(i)+9*yg(i)^2+xg(i)*(-5+12*yg(i)));
    DetaSF(i,10)=-27*xg(i)*(-1+xg(i)+2*yg(i)); 
    
    
    
    
    end


end