function en=enL2(unp1,un,N,dx,dt)
% calculate energy of solution at time t :
% E(t) = int |u_t|^2 + |u_x|^2

% Map the solution onto the grid 
Unp1 = reshape(unp1,N,N);

[ux,uy]=gradient(Unp1,dx);
ut=(unp1-un)/dt;

N=length(ux);
en2=sum(ux(:).^2) +sum(uy(:).^2) + sum(ut(:).^2);
% - .5*(ux(1).^2+ux(N).^2+ut(1).^2+ut(N).^2);
en=.5*(en2)*dx^2;