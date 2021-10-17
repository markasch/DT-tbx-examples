function u=transpImplicit(u0,x,dx,dt)
% Solve the transport equation 
%   u_t + c u_x = 0. 
% Implicit scheme.
global c
lamb = c*dt/dx;
N=length(x);
e=ones(N,1);
Mat1 = diag(sparse(ones(N,1)),0) ...
    -lamb*diag(sparse(ones(N-1,1)),-1)/2 ...
    +lamb*diag(sparse(ones(N-1,1)),1)/2;
Mat1(1,end-1)=-lamb/2;
Mat1(end,2)=lamb/2;
F=u0;
u=Mat1\F;
end
