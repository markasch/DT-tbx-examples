% Neumann problem
%   u''=f 
%   u'(a)=alpha, u(b)=beta
% using 1st and 2nd order approximations of the
% Neumann boundary condition  
clear
% Physical parameters 
a = 0;                           % domain endpoints
b = 1;
alpha = 1;                       % Neumann condition
beta  = 2;                       % Dirichlet condition 
% Numerical Parameters
n = 8;                         % no. of discretization points
h = (b-a)/(n-1);               % mesh size
k = 1;
% Functions
f      = @(x) exp(x);  % right hand side
u_true = @(x) exp(x) + (alpha-exp(a))*(x - b) + beta - exp(b);  % true solution

% Loop for convergence estimation
for k = 1:5
  % Mesh
  x = linspace(a,b,n)';
  % Assemble matrix A (using sparse matrix storage):
  E = ones(n,1);
  A = spdiags([E -2*E E], [-1 0 1], n, n);
  % 1st order approximation
  A(1,1) = -1 ;
  A(1,2) = +1;
  % 2nd order approximation (uncomment)
  % A(1,1) = -3 ;
  % A(1,2) = +4;
  % A(1,3) = -1 ;
  A(n,n) = 1 ;
  A(n,n-1) = 0;
  % RHS
  F    = h^2 * f(x);
  F(1) = h*alpha;   % 1st order
  %F(1)= 2*h*alpha; % 2nd order
  F(n) = beta;
  % Solve
  u = A\F;
  % Error
  err(k,1) = h;
  err(k,2) = max(u - u_true(x)); 
  n = n*2 ;
  h = (b-a)/(n-1);    
end
% Output Results
figure
plot(x,u,'--ro',x,u_true(x))
figure
plot(log10(err(:,1)),log10(err(:,2)),'o-',log10(err(:,1)),2*log10(err(:,1)),'r-',log10(err(:,1)),log10(err(:,1)),'m-')
xlabel('log(h)');
ylabel('log(error)');
legend(['Discretization error';'Slope h^2';'Slope h'],4);


