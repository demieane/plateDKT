function [BeSt]=BendingStiffness2(E,v,thick,h)

% calculates bending stiffness for constant properties
% thick = thick.*0 + h;
% thick = thick.*0 + 1/1000;%min(thick);

% thick = thick + 0.0005;%h;
%******************************************************
% Make sure to add this function to supress very low
% values of thickness
% Monday from OMAE paper or IMAM
%******************************************************
max(h);

BeSt=zeros(3,3,length(thick));

% BeSt=zeros(3,3);

% ---------------------- ISOTROPIC MATERIAL -------------------------------
% la=(E*thick^3/(12*(1-v^2)));
% -------------------------------------------------------------------------
for kk = 1:length(thick)   
    la=(E*thick(kk)^3/(12*(1-v^2))); %[1,Nelem]
    BeSt(1,1,kk)=la;
    BeSt(2,2,kk)=la;
    BeSt(3,3,kk)=(1-v)*la/2;
    BeSt(1,2,kk)=v*la;
    BeSt(2,1,kk)=v*la;
%     clear la;
end

% clear la
end