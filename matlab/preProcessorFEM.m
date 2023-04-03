clear all;
close all;
clc;
addpath('mesh/heathcote');
load('mesh_h2_half');
file1995 = 'bemDATA_h182_r';
load(file1995);
%
chord=inData.cRoot;%0.33; %wing chord length
span=inData.span;%1; % wing span length
Uvel =inData.U;%2.52;%m/s
fluid_dens=1025;%kg/m3
%
m=7850;
E=210e9;
v=0.28;
h=1/1000;%0.12*chord;
% Boundary conditions 
%--------------------------------------------------------------------------
% The present features a polygonal domain with 4 support senarios
% controled by the choice of CC
% CC=1 , simply supported along selected edges
% CC=2 , fully clamped along selected edges
%--------------------------------------------------------------------------
CC=2; % boundary condition toggle
% Special Cantileaver Case CC=4
if CC==1
    disp(' SS case along selected edge (w=0)')
elseif CC==2
    disp(' Clamped along selected edge (w=bx=by=0)')
end
% Forcing
% 1- concetrated load, 2- uniform load, 3- distributed load (mapping func)
lll=3;%2; %loading case
importFromFile=struct('toggle',1,'filename',file1995);
%
P_load = 1; %[Pa] %pointing towards the Z-axis
% in ANSYS load pointing in the negative of Z-axis is positive
if lll==1
   Pxy=[5,5];%load position
end
% Generate mesh - ID, IEN, LM 
pp=p;
tt=t;
ee=e;

%% CREATE MODE.bin binary file for passing data 
precision = 'single';

% write to binary for communication with GPU executable
file = fopen('INDATA_FEM.bin', 'wb');
fwrite(file, chord, precision);
fwrite(file, span, precision);
fwrite(file, Uvel, precision);
fwrite(file, fluid_dens, precision);
fwrite(file, m, precision);
fwrite(file, E, precision);
fwrite(file, v, precision);
fwrite(file, h, precision);
fwrite(file, CC, 'int');
fwrite(file, lll,'int');
% fwrite(file, ee, precision);
fwrite(file, size(pp,1),'int');
fwrite(file, size(pp,2),'int');
fwrite(file, [pp(1,:), pp(2,:)], precision);
% fwrite(file, tt, precision);
fclose(file);

system('cp INDATA_FEM.bin ../c/INDATA_FEM.bin')