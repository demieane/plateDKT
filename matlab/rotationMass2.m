function [ Hxx,Hyy ] = rotationMass2( k,a4,a5,a6,b4,b5,b6,c4,c5,c6,d4,d5,d6,e4,e5,e6 )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
A6=a6(k);
A4=a4(k);
A5=a5(k);
C6=c6(k);
C4=c4(k);
C5=c5(k);
B6=b6(k);
B4=b4(k);
B5=b5(k);
D6=d6(k);
D4=d4(k);
D5=d5(k);
E6=e6(k);
E4=e4(k);
E5=e5(k);



Hxx(:,1)=[0 6*A6 -6*A5 6*A5-6*A6 -6*A6 6*A5]';
Hxx(:,2)=[0 4*B6 4*B5 -4*B5-4*B6 -4*B6 -4*B5]';  
Hxx(:,3)=[1 -4*C6-3 -3-4*C5 4+4*C5+4*C6 4*C6+2 4*C5+2]';
Hxx(:,4)=[0 -6*A6 0 6*A4+6*A6 6*A6 0]';
Hxx(:,5)=[0 4*B6 0 4*B4-4*B6 -4*B6 0]' ;
Hxx(:,6)=[0 -1-4*C6 0 4*C6-4*C4 2+4*C6 0]';
Hxx(:,7)=[0 0 6*A5 -6*A5-6*A4 0 -6*A5]';
Hxx(:,8)=[0 0 4*B5 4*B4-4*B5 0 -4*B5]';
Hxx(:,9)=[0 0 -4*C5-1 4*C5-4*C4 0 4*C5+2]';

Hyy(:,1)=[0 6*D6 -6*D5 6*D5-6*D6 -6*D6 6*D5]';
Hyy(:,2)=[-1 4*E6+3 4*E5+3 -4*E5-4*E6-4 -4*E6-2 -4*E5-2]';
Hyy(:,3)=[0 -4*B6 -4*B5 4*B5+4*B6 4*B6 4*B5]';
Hyy(:,4)=[0 -6*D6 0 6*D4+6*D6 6*D6 0]';
Hyy(:,5)=[0 1+4*E6 0 4*E4-4*E6 -4*E6-2 0]';%%
Hyy(:,6)=[0 -4*B6 0 -4*B4+4*B6 4*B6 0]';
Hyy(:,7)=[0 0 6*D5 -6*D5-6*D4 0 -6*D5]';
Hyy(:,8)=[0 0 1+4*E5 4*E4-4*E5 0 -4*E5-2]';
Hyy(:,9)=[0 0 -4*B5 -4*B4+4*B5 0 4*B5]';

end

