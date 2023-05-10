function [ SFm, DxsiSFm,DetaSFm] = LNShapeFunMassDST(xg,yg) 
%UNTITLED Summary of this function goes here
%   %% Calculation of quadratic Lagrange Shape Functions 
% and their derivatives at the reference triangle


SFm=zeros(length(xg),3);
DxsiSFm=zeros(length(xg),3);
DetaSFm=zeros(length(xg),3);

for i=1:length(xg)
    
    SFm(i,1)=(1-xg(i)-yg(i)); %checked
    SFm(i,2)=xg(i);%checked
    SFm(i,3)=yg(i);
        

    % 1st derivative --> î
    DxsiSFm(i,1)=-1;%checked
    DxsiSFm(i,2)=1;%checked
    DxsiSFm(i,3)=0;%checked
    
     
    %1st derivative ç
    DetaSFm(i,1)=-1;
    DetaSFm(i,2)=0;
    DetaSFm(i,3)=1;

end


    