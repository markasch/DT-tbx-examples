%--- MBrown1.m ---%
% Simulate Brownian motion
randn('state', 100) % reference state
% Time interval, time step
T=1; N=500; dt=T/N;
% Initialize
dW = zeros(1,N);
W  = zeros(1,N);
% First approximation
dW(1) = sqrt(dt)*randn;
W(1) = dW(1);
% Loop over time
for j = 2:N
  dW(j) = sqrt(dt)*randn;
  W(j)  = W(j-1) + dW(j);
end
% Plot the path
plot([0:dt:T],[0,W],'r-')
xlabel('t'), ylabel('W(t)')
  
  
