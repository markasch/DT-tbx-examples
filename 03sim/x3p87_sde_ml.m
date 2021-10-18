%SDE_ML  Milstein method for nonlinear SDE
%
% dX = r*X*(K-X) dt + beta*X dW,   
% X(0) = X0,
% with r = 2, K= 1, beta = 1 and X0 = 0.5.
%
% Brownian path discretized over [0,1] with dt = 2^(-8).
% Milstein uses timestep = R*dt.
%
randn('state',100)
r = 2; K = 1; beta = 0.25; X0 = 0.5;   % parameters
T = 1; N = 2^(8); dt = T/N;         
dW = sqrt(dt)*randn(1,N);         % Brownian increments 
R = 2; Dt = R*dt; J = N/R;        % J steps of ML with Dt = R*dt
Xml = zeros(1,J);                 % preallocation for efficiency
Xtemp = X0;
for j = 1:J
   Winc   = sum(dW(R*(j-1)+1:R*j)); 
   Xtemp  = Xtemp + Dt*r*Xtemp.*(K-Xtemp) + beta*Xtemp.*Winc ...
            + 0.5*beta^2*Xtemp.*(Winc.^2 - Dt);
   Xml(j) = Xtemp;
end
% Plot solution
plot([0:Dt:T],[X0,Xml],'b--*')
xlabel('t'); ylabel('X(t)')