function u=transpExplicit(un,x,dx,dt)
% Solve the transport equation 
%   u_t + c u_x = 0. 
% Explicit scheme.
global c
if c < 0 then        
  disp('error: c must be positive')
end
lamb = c*dt/dx;
N=length(x);
u=zeros(N,1);
u(1)=un(end-1);
u(2:end)=un(2:end)-lamb*(un(2:end)-un(1:end-1));
end