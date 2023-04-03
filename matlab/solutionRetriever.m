function [solutionOUT] = solutionRetriever(GEN, Ndofs, d, Nt, u, solutionIN) %, time, T )
% This function retrives the data after the solution of the 
% system of ODE's has been completed

% Initialization
if d==2
    uu = zeros(GEN,Nt); %displacement
    uu_dot = zeros(GEN,Nt); %velocity
else
    %retrieve previous results
    uu=solutionIN.uu;
    uu_dot=solutionIN.uu_dot;
    
    solutionOUT.w=solutionIN.w;   % vertical displacement
    solutionOUT.bx=solutionIN.bx;  % rotation x
    solutionOUT.by=solutionIN.by;  % rotation y

    solutionOUT.w_dot=solutionIN.w_dot;   % vertical displacement
    solutionOUT.bx_dot=solutionIN.bx_dot;  % rotation x
    solutionOUT.by_dot=solutionIN.by_dot;  % rotation y
    
end

% RETRIEVING DATA
% u = [qdot; q] in my diploma thesis code it is the opposite

% apo to epilumeno
uu_dot(:,d) = u(1:GEN,d); 
uu(:,d) = u(Ndofs+1:Ndofs+GEN,d);

% w(:,d)=uu(1:3:end,d);   % vertical displacement
% bx(:,d)=uu(2:3:end,d);  % rotation x
% by(:,d)=uu(3:3:end,d);  % rotation y
% 
% w_dot(:,d)=uu_dot(1:3:end,d);   % vertical displacement
% bx_dot(:,d)=uu_dot(2:3:end,d);  % rotation x
% by_dot(:,d)=uu_dot(3:3:end,d);  % rotation y

solutionOUT.w(:,d)=uu(1:3:end,d);   % vertical displacement
solutionOUT.bx(:,d)=uu(2:3:end,d);  % rotation x
solutionOUT.by(:,d)=uu(3:3:end,d);  % rotation y

solutionOUT.w_dot(:,d)=uu_dot(1:3:end,d);   % vertical displacement
solutionOUT.bx_dot(:,d)=uu_dot(2:3:end,d);  % rotation x
solutionOUT.by_dot(:,d)=uu_dot(3:3:end,d);  % rotation y

solutionOUT.uu=uu;
solutionOUT.uu_dot=uu_dot;
end