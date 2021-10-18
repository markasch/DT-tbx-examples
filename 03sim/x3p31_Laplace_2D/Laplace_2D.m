%--- Laplace_2D.m ---%
% Solve  u_xx + u_yy = 1 in Omega
%        u = 0 in dOmega
% where Omega is an L-shaped region
% Parameters
clear
n = 32;  % number of mesh points
R = 'L'; % geometry of domain
% Geometry matrix
G = mynmgrid(R,n);
subplot(1,2,1); spy(G); axis square; % plot grid
% Discrete Laplacian matrix
A = mydelsq(G);
subplot(1,2,2);  spy(A); axis square; % sparse Laplacian matrix 
% Set up RHS
N = sum(G(:)>0);
b = ones(N,1);
% Solve the sparse system A U = b
u = A\b;
% Map solution back onto grid and plot
U = G;
U(G>0) = full(u(G(G>0)));
figure; clabel(contour(U));
prism; axis square ij