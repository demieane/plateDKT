function [TAU,EPSMAX,PTMAX] = indata(NACA)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
IEPS=floor(NACA/1000);
IPTMAX=floor(NACA/100)-10*IEPS;
ITAU=NACA-1000*IEPS-100*IPTMAX;
EPSMAX=IEPS*0.01;
PTMAX=IPTMAX*0.1;
TAU=ITAU*0.01;
end

