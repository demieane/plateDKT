% sparse matrix

A = zeros(3,3);
A(1,1) = 5;
A(1,2) = 4;
A(3,2) = 7;
A
As = sparse(A)

As(1)