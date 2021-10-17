function u=transpDec(un,x,dx,dt)
% Solve the transport equation 
%    u_t + c u_x = 0. 
% Downwind scheme.
global c
lamb = c*dt/dx;
if c<0 then        
    disp('error: c must be positive')
end
u(1)=un(end-1);
u(2:end)=un(2:end)-lamb*(un(2:end)-un(1:end-1));
end
