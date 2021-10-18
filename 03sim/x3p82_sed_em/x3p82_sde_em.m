%SDE_EM  Euler-Maruyama method for linear SDE
%
% dX   = b*X dt + sigma*X dW,   
% X(0) = X0.
%
% Brownian path discretized over [0,T] with step dt.
% Euler-Maruyama timestep Dt = k*dt.
%
randn('state',100)
b = 2; sigma = 1; X0 = 1;       % parameters
T = 1; N = 2^8; dt = T/N;         
dW = sqrt(dt)*randn(1,N);         % Brownian increments
W = cumsum(dW);                   % Discretized Brownian path
% Exact solution
Xtrue = X0*exp((b-0.5*sigma^2)*([dt:dt:T])+sigma*W); 
plot([0:dt:T],[X0,Xtrue],'r-'), hold on
k = 1; Dt = k*dt; J = N/k;        % J steps of EM with Dt = k*dt
Xem = zeros(1,J);                 % initialize
Xtemp = X0;
for j = 1:J
   Winc   = sum(dW(k*(j-1)+1:k*j)); 
   Xtemp  = Xtemp + Dt*b*Xtemp + sigma*Xtemp*Winc;
   Xem(j) = Xtemp;
end
% Plot numerical solution
plot([0:Dt:T],[Xzero,Xem],'b--*'), hold off
xlabel('t'); ylabel('X(t)')
% Calculate error at t=T
emerr = abs(Xem(end)-Xtrue(end))