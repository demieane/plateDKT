function [ Bb ] = Kbend( Area,k,Hx_xsi, Hx_eta,Hy_xsi, Hy_eta,y12,y31,x12,x31,b1,b2,c1,c2 )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% Bb(1,:)=1/(2*Area(k))*(y31(k)*Hx_xsi+y12(k)*Hx_eta);
% Bb(2,:)=1/(2*Area(k))*(-x31(k)*Hy_xsi-x12(k)*Hy_eta);
% Bb(3,:)=1/(2*Area(k))*(-x31(k)*Hx_xsi-x12(k)*Hx_eta+y31(k)*Hy_xsi+y12(k)*Hy_eta);

Bb(1,:)=b1(k)*Hx_xsi+b2(k)*Hx_eta;
Bb(2,:)=c1(k)*Hy_xsi+c2(k)*Hy_eta;
Bb(3,:)=c1(k)*Hx_xsi+c2(k)*Hx_eta+b1(k)*Hy_xsi+b2(k)*Hy_eta;

end

