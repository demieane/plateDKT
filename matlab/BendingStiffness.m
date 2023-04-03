function [BeSt]=BendingStiffness(E,v,thick)

% calculates bending stiffness for constant properties

% BeSt=zeros(3,3,length(thick));

BeSt=zeros(3,3);

% ---------------------- ISOTROPIC MATERIAL -------------------------------
la=(E*thick^3/(12*(1-v^2)));

    
BeSt(1,1)=la;
BeSt(2,2)=la;
BeSt(3,3)=(1-v)*la/2;
BeSt(1,2)=v*la;
BeSt(2,1)=v*la;

clear la
end