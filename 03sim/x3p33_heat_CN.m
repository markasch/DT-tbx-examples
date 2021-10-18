%--- heat_CN.m ---%
% Crank-Nicolson finite difference scheme for the 1D heat equation
%  u_t = u_xx in [0,1] x (0,T],
%  u(0,t) = u(1,t) = 0,
%  u(x,0) = f(x).
%
clear
% Physical parameters 
T = 10;     % final time 
N = 40;     % number of mesh points x 
% Numerical parameters
r = input(' Enter r = '); % r=k/h^2
h = 1/N; k = r*(h^2); ntmax = min(ceil(T/k),100);
x = [0:N].'/N;                    
N1 = N + 1; 

U = zeros(N1,ntmax+1); % storage of results

% Initial condition 
uo=(x < 0.5001).*(2*x) + (x > 0.5).*(2*(1-x)) ; 
uo(1) = 0; uo(N+1) = 0;  % impose zero BC for compatibility with IC
U(:,1) = uo; 
% Difference matrices
C1 = -r; C0 = 2 + 2*r; 
D1 = r;  D0 = 2 - 2*r;
E = ones(N-1,1); 
C = spdiags([C1*E C0*E C1*E], -1:1, N-1, N-1); 
D = spdiags([D1*E D0*E D1*E], -1:1, N-1, N-1); 
% Time loop
v = uo(2:N);
for n = 2:ntmax+1
    q = D*v;
    v = C\q;
    U(2:N,n) = v; 
end
% Plot the solution
hold on 
for n = 2:10:ntmax+1
    plot(x,U(:,n)), grid, xlabel('x'), ylabel('u(x,t)'); 
end 
hold off 