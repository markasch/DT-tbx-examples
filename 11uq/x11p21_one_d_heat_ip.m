%--- OneDHeatIP.m ---%
% Solve parameter estimation inverse problem 
% for the 1D heat equation
%
% Problem parameters
L = 10;
D = 0.5;
s = 0.03;
Tmax = 2;
xdim = 25; tdim = 75;
% Discretization
x = linspace(0,L,xdim);
t = linspace(0,Tmax,tdim);
dx = x(2)-x(1); dt = t(2)-t(1);
q = dt/dx^2;
r1 = 0.75*L; r2 = 0.8*L;
% Initial condition
u0 = zeros(1,xdim);
u0(x>=r1 & x<=r2) = 1;
xDat = 2:xdim-1;
tDat = tdim;
nxDat = length(xDat);
ntDat = length(tDat);
% Solve the heat equation
Z = heat(D,u0,q,tdim);
u = Z(tDat,xDat);
uDat = u + s*randn(ntDat,nxDat);
% Metropolis-Hastings MCMC
N = 10000; m = 100; XD = 1; X = zeros(1,N); X(1) = XD;
Z = heat(XD,u0,q,tdim);
u = Z(tDat,xDat);
oLLkd = sum(sum(-(u-uDat).^2))/(2*s^2);
LL = zeros(1,N); LL(1) = oLLkd;
w = 0.1;
% Main loop
for n = 2:N
   XDp = XD + w*(2*rand()-1);
   if XDp > 0
      Z = heat(XDp,u0,q,tdim);
      u = Z(tDat,xDat);
      nLLkd = sum(sum(-(u-uDat).^2))/(2*s^2);
      alpha = exp(nLLkd - oLLkd);
      if rand() < alpha
         XD = XDp;
         oLLkd = nLLkd;
         CZ = Z;
      end
   end
   X(n) = XD; LL(n) = oLLkd;
end
figure; plot(X);
%-----------------------------%
function Z = heat(D,u0,q,tdim)
xdim = length(u0);
Z = zeros(tdim,xdim);
Z(1,:) = u0;
for tin = 2:tdim
   tip = tin - 1;
   Z(tin,2:end-1) = Z(tip,2:end-1) + ...
       D*q*(Z(tip,1:end-2)-2*Z(tip,2:end-1)+Z(tip,3:end));
end