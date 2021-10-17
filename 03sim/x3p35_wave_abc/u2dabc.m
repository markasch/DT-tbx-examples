%--- u2abc.m ---%
% Solves the 2-D wave-equation 
%      BOX(u) = 0, in Omega x (0,T)
%      u = f(x,y)  on dOmega, t > 0 (== Sigma)
%      u = u0(x,y,0) , du/dt = u1(x,y,0) in Omega
%
% by an explicit finite-difference scheme, (central-differences)
% for a given number of time-steps
%
% NOTES:
% 1. Neumann boundary conditions at x=0, y=0
% 2. Absorbing boundary conditions at x=1, y=1
% 3. Special corner treatment in 2-d
% 
clear; clf;
h  = input(' spatial increment, h = ');
rmax=1/sqrt(2);
disp(sprintf('Maximal value of r for stability, r_max = %g',rmax))
rf  = input('Ratio of r(=dt/dx) over CFL : ');
Tf = input('Final time = ');
Delt  = input('Delta(t) for plotting = ');
r=rf*rmax;
dt=h*r; k=dt; % r=(k/h);
Nt=ceil(Tf/dt); Nplot=round(Delt/dt);

mflag=input('Movie (1) or not(0) ? ');
if mflag==1; ifr=1; end

xmin = 0; xmax = 1; ymin = 0; ymax = 1;
x = xmin:h:xmax; y = ymin:h:ymax;
I = length(x); In = [2:I-1]; Bn = [1 I];

[X,Y]=meshgrid(xmin:h:xmax,ymin:h:ymax);

disp(sprintf('Space discretisation, h = %g',h))
disp(sprintf('Time discretisation, k  = %g',dt))
disp(sprintf('Number of time steps,N  = %g',Nt))
% Generate the grid 
G = numgrid('S',I+2);
% Generate the discrete Laplacian.
D = delsqr(G,r);
% coefficients for ABC
r1=(r-1)/(r+1);
r21=(-r^2+r-1)/(r+1);
r22=r^2/(2*(r+1));
r23=2/(r+1);
rc1=(3*r-2*sqrt(2))/(6*r+2*sqrt(2));
rc2=(3*r+2*sqrt(2))/(6*r+2*sqrt(2));
rc3=(sqrt(2)-3*r)/(sqrt(2)+3*r);
% Number of interior points
N = sum(G(:)>0);
% Numbering of boundary points:
No  =I:I:I^2; % North
So = 1:I:(I-1)*I+1; % South
We = 1:I; % West
Ea =(I-1)*I+1:I^2; % East
% For 2nd order condition:
dum = I^2+1; % fictitious point index
M = 1:I^2; % dim for U
NEm1 = No-1 + I; NEm1(I) = dum; % j+1,K-1
NWm1 = No-1 - I; NWm1(1) = dum; % j-1,K-1
NE = No + I; NE(I )= dum;       % j+1,K
NW = No - I; NW(1) = dum;       % j-1,K
ENm1 = Ea-I + 1; ENm1(I) = dum; % J-1,k+1
ESm1 = Ea-I - 1; ESm1(1) = dum; % J-1,k-1
EN = Ea + 1; EN(I) = dum;       % J,k+1
ES = Ea - 1; ES(1) = dum;       % J,k-1
JK = I^2; JKm1 = I^2-1; Jm1Km1 = I^2-I-1; JKm2 = I^2-2; 
Jm2K = I^2-2*I; Jm1K = I^2-I; 
J1 = I^2 - I + 1; K1 = I; 
% rhs = 0
rhs = zeros(N,1);
% Set-up initial conditions and RHS in interior
ic=menu('Initial condition','Exponential','Sinusoidal','Expanding wave');
if ic==1
  x0 = input(' x0 = '); alpha = input('alpha = ');
  xy0 = (X-x0).^2 + (Y-x0).^2;
  u0 = exp(-alpha*xy0); u1 = zeros(size(u0));
  vmax=.6; vmin=-.6; dv=.05; V=[vmin:dv:vmax];
elseif ic==2
  n = input('Wave number of initial condition (sin[n.pi.x]) = ');
  u0 = sin(n*pi*X).*sin(n*pi*Y); u1 = zeros(size(u0));
  vmin=-1; vmax=1; dv=.1; V=[vmin:dv:vmax];
elseif ic==3
  r0 = input(' r0 = '); d = input(' d = ');
  r = sqrt(X.^2+Y.^2);
  arg0 = r-r0; arg1 = r-r0-dt;
  p = find(r<d);
  u0 = zeros(size(X)); u1 = zeros(size(X));
  u0(p) = 0.5*(1+cos(pi*arg0(p)/d));
  u1(p) = 0.5*(1+cos(pi*arg1(p)/d));
  vmin=-.5; vmax=.5; dv=.05; V=[vmin:dv:-dv dv:dv:vmax];
end
% Choose order of absorbing BC 
abc = menu('Absorbing Boundary Condition','1st order','2nd order','None');
disp(sprintf('Using absorbing BC of order %g',abc));
if abc==3  % choose BC type
  bc = menu('Boundary Condition','Dirichlet','Neumann');
  if bc==1
    bctype = 'Dir'
  elseif bc==2
    bctype = 'Neu'
  end
end

% Indices for interior of the region
In=find(G);
% Convert matrix interior -> vector (for SQUARE)
u0=u0(:); u1=u1(:);
% Fictitious values
u0(dum) = 0; u1(dum) = 0; unm1(dum) = 0; un(dum) = 0; unp1(dum) = 0;
% First time step:
unm1 = u0; un = u0; unp1 = u0; t = 0;
if ic==3; 
   un = u1;   % for cylindrical wave IC, u^0, u^1 given
   unm1 = 2*u0 - u1;
else 
   un(M) = 0.5*D*u0(M) + k*u1(M);
   unm1(M) = un(M) - 2*k*u1(M);
end
unp1=un;
if abc==1
   unp1(No) = r1*un(No-1) - r1*u0(No) + u0(No-1); % y=1
   unp1(Ea) = r1*un(Ea-I) - r1*u0(Ea) + u0(Ea-I); % x=1
elseif abc==2
   unp1(No) = r21*(unp1(No-1)+unm1(No))+ ...
        r22*(unp1(NEm1)+unp1(NWm1)+unm1(NE)+unm1(NW))+...
        r23*(u0(No-1)+u0(No)) - unm1(No-1);
   unp1(Ea) = r21*(unp1(Ea-I)+unm1(Ea))+ ...
        r22*(unp1(ENm1)+unp1(ESm1)+unm1(ES)+unm1(EN))+...
        r23*(u0(Ea-I)+u0(Ea)) - unm1(Ea-I);
end
unm1=u0; % for next time-step

% calculate initial energy : E(0)
it=1;
if ic==3 
  energy(it)=enL2(un(M),unm1(M),I,h,dt); 
else
  energy(it)=enL2(u0(M),u0(M)-dt*u1(M),I,h,dt); % gives du/dt(t=0)=u1  
end;
% E(dt) :
it=2; energy(it)=enL2(un(M),unm1(M),I,h,dt);

% Map the solution onto the grid - including bdy points
Un=reshape(un(M),I,I); % valid for square ONLY !!!

%  display results
surf(X,Y,Un,'edgecolor','none','FaceLighting','gouraud')
shading interp; view(-15,50);
title(['u(x,y,t) for t=',num2str((+1)*k)])
axis tight manual
ax = gca;
ax.NextPlot = 'replaceChildren';

if mflag==1
  Mv(Nplot+1) = struct('cdata',[],'colormap',[]); 
  drawnow; Mv(ifr)=getframe; ifr=ifr+1; 
end
pause;
%% Loop on time:
disp('Looping on time ...')
for nt=1:Nt   
  t=(nt+1)*dt;
  % solve for t=n+1
  unp1(M) = D*un(M) - unm1(M); % + k^2*urhs(X,Y,G,nt*k);
  % boundary conditions for Unp1 
  unp1(So)=unp1(So+1); unp1(We)=unp1(We+I); % Neumann
  if abc==1
     unp1(No) = r1*unp1(No-1) - r1*un(No) + un(No-1); % absorb, y=1
     unp1(Ea) = r1*unp1(Ea-I) - r1*un(Ea) + un(Ea-I); % absorb, x=1
  else
    unp1(No)=r21*(unp1(No-1)+unm1(No))+ ...
	r22*(unp1(NEm1)+unp1(NWm1)+unm1(NE)+unm1(NW))+...
        r23*(un(No-1)+un(No)) - unm1(No-1);
    unp1(Ea)=r21*(unp1(Ea-I)+unm1(Ea))+ ...
        r22*(unp1(ENm1)+unp1(ESm1)+unm1(ES)+unm1(EN))+...
        r23*(un(Ea-I)+un(Ea)) - unm1(Ea-I);
    % fix corners (J,K-1), (J-1,K), (J,K).
    unp1(JKm1)=rc1*(unp1(Jm1Km1)+unp1(JKm2)) + ...
	rc2*(un(Jm1Km1)+un(JKm2)) + rc3*un(JKm1);
    unp1(Jm1K)=rc1*(unp1(Jm2K)+unp1(Jm1Km1)) + ...
	rc2*(un(Jm2K)+un(Jm1Km1)) + rc3*un(Jm1K);
    unp1(JK)=rc1*(unp1(Jm1K)+unp1(JKm1)) + ...
	rc2*(un(Jm1K)+un(JKm1)) + rc3*un(JK);
    % fix corners (J,1), (1,K) by putting the 1st order RBC
    unp1(K1)= r1*unp1(K1-1) - r1*un(K1) + un(K1-1); % absorb, y=1
    unp1(J1)= r1*unp1(J1-I) - r1*un(J1) + un(J1-I); % absorb, x=1
  end
  % Calculate the energy of the solution
  energy(nt+2)=enL2(unp1(M),un(M),I,h,dt);
  % update (except for last step)
  if nt < Nt
    unm1 = un  ;
    un   = unp1;
  end
  remt=rem(nt,Nplot);
  if remt==0
    % Map the solution onto the grid 
    Unp1=reshape(unp1(M),I,I);
    surf(X,Y,(Unp1),'edgecolor','none','FaceLighting','gouraud')
    shading interp
    axis([xmin xmax ymin ymax vmin vmax]);
    title(['u(x,y,t) for t=',num2str((nt+1)*k)]); view(-15,50);%view(-45,75);
    disp(sprintf('t = %g, Key to continue',t)); pause
    if mflag==1; drawnow; Mv(ifr)=getframe; ifr=ifr+1; end
  end
end

disp('Displaying results ...')
 %play the movie
 if mflag==1; movie(Mv,2,5); end

% Energy (inf-norm)
errmax=norm(unp1(M),inf);
disp(sprintf(' E(u(0)) = %0.5g',energy(1)));
disp(sprintf(' E(u(T)) = %0.5g',energy(Nt+2)));
disp(sprintf(' E(u(T))/E(u(0)) = %0.5g',energy(Nt+2)/energy(1)));
disp(sprintf(' |u(T)|_(inf) = %0.5g',errmax));

% plot E(t)
tvec=0:dt:t; Nt=length(energy);
figure; plot(tvec,energy/energy(1)); %axis([0 t 0 2]);
title('Energy'); xlabel('t'); ylabel('E(t)');