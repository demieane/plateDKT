% Eigenvalue problem
addpath('mesh/heathcote');
addpath('dataFSI');
load('mesh_h1_half');

m=7800;%kg/m^2
E=200e9;%Pa
v=0.3;
% average thickness
h=1/1000;%m

D = E*h^3/12/(1-v^2);

a = 10;
omega = 8.777;
lambda = omega*a^2*sqrt(m/D)
