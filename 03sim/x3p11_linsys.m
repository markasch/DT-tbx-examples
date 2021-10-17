% Compare methods for solving Ax=b,
% where A is a Hilbert matrix of order n
%       b is the RHS for x=[1 1 ... 1]^T
clear
n = 14;
A = hilb(n); x = ones(n,1); b = A * x;
% solve with LU
x_LU = A \ b;
% compute error
err_LU = norm(x-x_LU)/norm(x);
% solve with PCG
[x_PCG,flag,relres,iter,resvec] = pcg(A,b);
% compute error
err_CG = norm(x-x_PCG)/norm(x);
% output
fprintf('\n Condition number of A  = %0.4e \n', cond(A))
fprintf('\n Relative error for LU  = %0.4e \n', err_LU)
fprintf('\n Relative error for PCG = %0.4e \n', err_CG)

