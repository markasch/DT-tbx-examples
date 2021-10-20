% Kalman Filter for scalar Gaussian random walk
% Set parameters
sig_w = 1; sig_v = 0.5;
M = 1;
Q = sig_w^2;
H = 1;
R = sig_v^2;
% Initialize
m0 = 0;
P0 = 1;
% Simulate data
randn('state',1234);
steps = 100; T = [1:steps];
X = zeros(1,steps);
Y = zeros(1,steps);
x = m0;
for k=1:steps
  w = Q'*randn(1);
  x = M*x + w;
  y = H*x + sig_v*randn(1);
  X(k) = x;
  Y(k) = y;
end
plot(T,X,'-',T,Y,'.');
legend('Signal','Measurements');
xlabel('{k}'); ylabel('{x}_k');
% Kalman filter
m = m0;
P = P0;
for k=1:steps
  m = M*m;
  P = M*P*M' + Q;

  d = Y(:,k) - H*m;
  S = H*P*H' + R;
  K = P*H'/S;
  m = m + K*d;
  P = P - K*S*K';

  kf_m(k) = m;
  kf_P(k) = P;
end
% Plot   
clf; hold on
fill([T fliplr(T)],[kf_m+1.96*sqrt(kf_P) ...
  fliplr(kf_m-1.96*sqrt(kf_P))],1, ...
  'FaceColor',[.9 .9 .9],'EdgeColor',[.9 .9 .9])
plot(T,X,'-b',T,Y,'or',T, kf_m(1,:),'-g')
plot(T,kf_m+1.96*sqrt(kf_P),':r',T,kf_m-1.96*sqrt(kf_P),':r');
hold off
xlabel('{ k}'); ylabel('{ x}_k');