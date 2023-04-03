

function [a,b] = rayleighDampCoef(om,m,n,zeta_1,zeta_m)

% 2*zeta_i = a + b*om_1^2
% om_i: i-th mode natural frequency, zeta_i: corresponding damping coef.

% Reference: Computation of Rayleigh Coefficients for large systems,
% I. Chowdhury, S. Dasgupta, 2003

% ---------------------- PROCEDURE -----------------------------
% 1. Select zeta_1
% 2. Select zeta_m, for m-th significant mode
% 3. Linear interpolation for intermediate modes
% 4. Linear extrapolation for higher modes, until 2.5m-th mode
% 5. b = (2*zeta_1*om_1 - 2*zeta_m*om_m)/(om_1^2 - om_m^2)
% 6. a = 2*zeta_1 - b*om_1^2
% 7. Calculate a,b based on om_1/zeta_1, om_2.5m,zeta_2.5m
% 8. Average values of a, b

% Input: natural frequencies vector, m, 2.5*m, corresponding zeta values

% 1. 1,...,m
b_1 = (2*zeta_1*om(1) - 2*zeta_m*om(m))/(om(1)^2 - om(m)^2);
a_1 = 2*zeta_1*om(1) - b_1*om(1)^2;

% Get zeta_max (nth - mode)
zeta_max = (zeta_m - zeta_1)*(om(n) - om(1))/(om(m) - om(1)) + zeta_1;

% 2. 1,...,n
b_2 = (2*zeta_1*om(1) - 2*zeta_max*om(n))/(om(1)^2 - om(n)^2);
a_2 = 2*zeta_1*om(1) - b_2*om(1)^2;


a = 0.75*(a_1) + 0.25*(a_2); % Could be mean value instead
b = 0.75*(b_1) + 0.25*(b_2); % Same as above

return;



