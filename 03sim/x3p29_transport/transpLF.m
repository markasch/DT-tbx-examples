function u=transpLF(un,x,dx,dt)
% Solve the transport equation u_t + c u_x = 0 
% Lax-Friedrichs scheme.
global c
lamb = c*dt/dx;
N=length(x);
t = 0;
u=zeros(N,1);
u(1)=un(end-1);
u(2:end)=(1-lamb)*[un(3:end);un(1)]/2+(1+lamb)*un(1:end-1)/2;
end
