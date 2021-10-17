%--- transport.m ---%
% Solve u_t + c u_x = 0  on [a,b] with periodic boundary conditions,
% using the following methods with N interior points:
%
% -- Explicit
% -- Lax-Friedrichs
% -- Lax-Wendroff
% -- Leapfrog
% -- Implicit
%

global c
c  = 10;
x0 = 0.5;
Tf = 0.1;

u0 = @(x) exp(-500*((x-x0).^2));
uex = @(x,t) (u0(x-c*t)+u0(x+2*x0-c*t)+u0(x+4*x0-c*t)+...
     u0(x+6*x0-c*t)+u0(x+8*x0-c*t)+u0(x+10*x0-c*t)+ ...
     u0(x+12*x0-c*t)+u0(x+14*x0-c*t));

errorExp=[];
errorLF =[];
errorLW =[];
errorLP =[];
errorImp =[];
step = [];

N=50
for j=1:5
    x  = linspace(0,1,N);
    dx = x(2)-x(1);

    cfl = 0.99;
    dt = cfl*dx/c;
     
    T = 0:dt:Tf;
    M=length(T);

    uExp=u0(x');
    uLF=u0(x');
    uLW=u0(x');
    uLPm1=u0(x');    
    uLP=transpLW(uLPm1,x,dx,dt);
    uImp=u0(x');
    uExact=u0(x');

    errExp = 0;
    errLF  = 0;
    errLW  = 0;
    errLP  = 0;
    errImp = 0; 
      
    for i=1:M-1
        uExact = uex(x',T(i+1));
        uExp   = transpExplicit(uExp,x,dx,dt);
        uLF    = transpLF(uLF,x,dx,dt);
        uLW    = transpLW(uLW,x,dx,dt);
        uLPp1  = transpLeapFrog(uLP,uLPm1,x,dx,dt);
        uImp   = transpImplicit(uImp,x,dx,dt);
        
        errExp = errExp + sum((uExact-uExp).^2);
        errLF  = errLF  + sum((uExact-uLF).^2);
        errLW  = errLW  + sum((uExact-uLW).^2);
        errLP  = errLP  + sum((uExact-uLP).^2);
        errImp = errImp + sum((uExact-uImp).^2);
        uLPm1 = uLP; uLP = uLPp1;
    end
    
    errorExp=[errorExp;sqrt(dx*dt*errExp)];
    errorLF =[errorLF;sqrt(dx*dt*errLF)];
    errorLW =[errorLW;sqrt(dx*dt*errLW)];
    errorLP =[errorLP;sqrt(dx*dt*errLP)];
    errorImp=[errorImp;sqrt(dx*dt*errImp)];
    step=[step;dx];
    
    plot(x,uImp,x,uLF,x,uLW,x,uLP,x,uImp,x,uExact)
    h1=legend('Expl','LF','LW','Leap','Impl','Exact');
    pause
 
    N=N*2
end

figure();  hold on
loglog(step,(errorExp),'-o',step,(step/10),'-.')
loglog(step,(errorLF),'r-o')
loglog(step,(errorLW),'m-o',step,(step.^2*2),'-.')
loglog(step,(errorLP),'c-x')
loglog(step,(errorImp),'b-o')
h1=legend('Expl','h','LF','LW','h^2','Leap','Impl');
hold off
% Compute best fit to convergence 
LSerror('Explicit',5,step,errorExp);
LSerror('Implicit',5,step,errorImp);
LSerror('Lax-Friedrichs',5,step,errorLF);
LSerror('Leapfrog',5,step,errorLP);
LSerror('Lax-Wendroff',5,step,errorLW);

