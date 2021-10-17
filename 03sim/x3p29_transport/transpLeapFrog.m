function u=transpLeapFrog(un,unm1,x,dx,dt)
% Solve the transport equation 
%  u_t + c u_x = 0. 
% Leapfrog scheme.
global c
lambda=c*dt/dx;
u=unm1-lambda*([un(2:end);un(2)]-[un(end-1);un(1:end-1)]);
end
