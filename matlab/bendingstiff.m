function [ Bb ] = bendingK( Area,k,Hx_xsi, Hx_eta,Hy_xsi, Hy_eta,b2,b3,c2,c3 )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

Bb(1)=1/(2*Area(k))*(b2(k)*Hx_xsi+b3(k)*Hx_eta);
Bb(2)=1/(2*Area(k))*(-c2(k)*Hy_xsi-c3(k)*Hy_eta);
Bb(3)=1/(2*Area(k))*(-c2(k)*Hx_xsi-c3(k)*Hx_eta+b2(k)*Hy_xsi+b3(k)*Hy_eta);


end

