function uNEW = timeIntegration(u, d, GEN, Mglob, Kglob, C, G, dt, theta)
% 
%  A = [sparse(GEN,GEN) Mglob; speye(GEN,GEN) sparse(GEN,GEN)];
%  B = -[Kglob C; sparse(GEN,GEN) -speye(GEN,GEN)];


sizeM=size(Mglob,1);

A = [Mglob, sparse(sizeM,sizeM); sparse(sizeM,sizeM), speye(sizeM,sizeM)];
B = -[C, Kglob; -speye(sizeM,sizeM), sparse(sizeM,sizeM)];

%lamda (1: implicit Euler, 1/2: Crank - Nicolson)

% AA =  A + theta*dt*B;
% BB =  A - (1 - theta)*dt*B;

AA =  A - theta*dt*B;
BB =  A + (1 - theta)*dt*B;

full(AA(1:10,1:10))

full(BB(1:10,1:10))

error('debug');


Q = (1 - theta)*dt*G(:,d-1) + (theta)*dt*G(:,d);

u(:,d) = AA\(BB*u(:,d-1) + Q);

uNEW = u(:,d);

end