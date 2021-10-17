function u=transpLW(un,x,dx,dt)
% Solve the transport equation 
%   u_t + c u_x = 0. 
% Lax-Wendroff scheme.
global c
lamb = c*dt/dx;
N = length(x);
u = zeros(N,1);
u = (1-lamb^2)*un ...
    + lamb*(lamb-1)*[un(2:end);un(2)]/2 ...
    + lamb*(1+lamb)*[un(end-1);un(1:end-1)]/2;            
end
