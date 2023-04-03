function [ GGDST,GGDKT ] = matrixG(  )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here


GGDST=zeros(10,10);
xsi=0; eta=0;
GGDST(1,:)=[1 xsi eta xsi*eta xsi^2 eta^2 xsi^2*eta eta^2*xsi xsi^3 eta^3];

xsi=1; eta=0;
GGDST(2,:)=[1 xsi eta xsi*eta xsi^2 eta^2 xsi^2*eta eta^2*xsi xsi^3 eta^3];

xsi=0; eta=1;
GGDST(3,:)=[1 xsi eta xsi*eta xsi^2 eta^2 xsi^2*eta eta^2*xsi xsi^3 eta^3];

xsi=1/3; eta=1/3;
GGDST(4,:)=[1 xsi eta xsi*eta xsi^2 eta^2 xsi^2*eta eta^2*xsi xsi^3 eta^3];

xsi=2/3; eta=1/3;
GGDST(5,:)=[1 xsi eta xsi*eta xsi^2 eta^2 xsi^2*eta eta^2*xsi xsi^3 eta^3];

xsi=1/3; eta=2/3;
GGDST(6,:)=[1 xsi eta xsi*eta xsi^2 eta^2 xsi^2*eta eta^2*xsi xsi^3 eta^3];

xsi=0; eta=2/3;
GGDST(7,:)=[1 xsi eta xsi*eta xsi^2 eta^2 xsi^2*eta eta^2*xsi xsi^3 eta^3];

xsi=0; eta=1/3;
GGDST(8,:)=[1 xsi eta xsi*eta xsi^2 eta^2 xsi^2*eta eta^2*xsi xsi^3 eta^3];

xsi=1/3; eta=0;
GGDST(9,:)=[1 xsi eta xsi*eta xsi^2 eta^2 xsi^2*eta eta^2*xsi xsi^3 eta^3];

xsi=2/3; eta=0;
GGDST(10,:)=[1 xsi eta xsi*eta xsi^2 eta^2 xsi^2*eta eta^2*xsi xsi^3 eta^3];



%%

xsi=0; eta=0;
GGDKT(1,:)=[1 xsi eta xsi*eta xsi^2 eta^2];

xsi=1; eta=0;
GGDKT(2,:)=[1 xsi eta xsi*eta xsi^2 eta^2 ];

xsi=0; eta=1;
GGDKT(3,:)=[1 xsi eta xsi*eta xsi^2 eta^2 ];

xsi=1/2; eta=1/2;
GGDKT(4,:)=[1 xsi eta xsi*eta xsi^2 eta^2 ];

xsi=0; eta=1/2;
GGDKT(5,:)=[1 xsi eta xsi*eta xsi^2 eta^2];

xsi=1/2; eta=0;
GGDKT(6,:)=[1 xsi eta xsi*eta xsi^2 eta^2 ];







end

