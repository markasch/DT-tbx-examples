function [u0, u1]=u0cyl(X,Y,dt)
% initial condition for wave eqn: E&M's expanding cylindrical
% wave [CPAM 23]

% parameters:
% r0=0.2; d=0.2;
r0=input(' r0 = '); d=input(' d = ');

r=sqrt(X.^2+Y.^2);
arg0=r-r0;
arg1=r-r0-dt;

p=find(r<d);
u0=zeros(size(X)); u1=zeros(size(X));
u0(p)=0.5*(cos(pi*arg0(p)/d)+1);
u1(p)=0.5*(cos(pi*arg1(p)/d)+1);

% plot :
subplot(2,1,1), mesh(X,Y,u0),title('u(x,y,0)'),view(62,36);
xlabel('x'),ylabel('y'),zlabel('u');
subplot(2,1,2), mesh(X,Y,u1),title('u(x,y,dt)'); view(32,34);
xlabel('x'),ylabel('y'),zlabel('u');
disp('Key to continue'); pause;