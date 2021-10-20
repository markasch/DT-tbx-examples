%--- UKfNav2D.m ---%
% Unscented Kalman Filter for 2D navigation
%------------------%
clear all
randn('state',123)
%% Parameters
dT = 0.1; % Time step
N  = 600; % Number of time steps for filter
xyA = [0 20]; % Station A coordinates
xyB = [20 0]; % Station B coordinates
% Step 1: Define UT Scaling parameters and weight vectors
n     = 4; % Dimension of state vector
alpha = 1; % Primary scaling parameter
beta  = 2; % Secondary scaling parameter (Gaussian assumption)
kappa = 0; % Tertiary scaling parameter
lambda = alpha^2*(n+kappa) - n;
wm = ones(2*n + 1,1)*1/(2*(n+lambda));
wc = wm;
wm(1) = lambda/(lambda+n);
wc(1) = lambda/(lambda+1) + 1 - alpha^2 + beta;
% Step 2: Define noise
Q = diag([0 0 4 4]); % Process
R = diag([1 1]);     % Measurement
% Step 3: Initialize state and covariance
x = zeros(4, N); % Initialize size of state estimate for all k
x(:,1) = [0; 0; 50; 50]; % Set initial state estimate
P0 = eye(4,4); % Set initial error covariance
% Simulate true state trajectory, measurement vector
w = sqrt(Q)*randn(4, N); % Generate random process noise from Q
v = sqrt(R)*randn(2, N); % Generate random measurement noise from R
xt = zeros(4, N); % Initialize size of true state for all k
xt(:,1)=[0; 0; 50; 50]+sqrt(P0)*randn(4,1); %Set true initial state
yt = zeros(2, N); % Initialize size of output vector for all k
M = [1 0 dT 0; 0 1 0 dT; 0 0 1 0; 0 0 0 1]; % Dynamics
for k = 2:N
  xt(:,k) = M*xt(:,k-1) + w(:,k-1);
  yt(:,k) = [sqrt((xt(1,k)-xyA(1))^2 + (xt(2,k)-xyA(2))^2); ...
  sqrt((xt(1,k)-xyB(1))^2 + (xt(2,k)-xyB(2))^2)] + v(:,k);
end
%% Initialize and run EKF for comparison
xe = zeros(4,N);
xe(:,1) = x(:,1);
P = P0;
M = [1 0 dT 0; 0 1 0 dT; 0 0 1 0; 0 0 0 1]; % Linear prediction
for k = 2:N
  % Prediction
  x_m = M*xe(:,k-1);
  P_m = M*P*M' + Q;
  % Observation
  dA = sqrt((x_m(1)-xyA(1)).^2 + (x_m(2)-xyA(2)).^2);
  dB = sqrt((x_m(1)-xyB(1)).^2 + (x_m(2)-xyB(2)).^2);
  y_m = [dA ; dB];
  H = [(x_m(1)-xyA(1))/dA, ...
  (x_m(2)-xyA(2))/dA, 0, 0; ...
  (x_m(1)-xyB(1))/dB, ...
  (x_m(2)-xyA(2))/dB, 0, 0];
  % Measurement Update
  K = P_m*H'/(H*P_m*H' + R); % Calculate Kalman gain
  xe(:,k) = x_m + K*(yt(:,k) - y_m); % Update state estimate
  P = (eye(4)-K*H)*P_m; % Update covariance estimate
end
%% Execute Unscented Kalman Filter
P = P0; % Set first value of P to the initial P0
for k = 2:N
  % Step 1: Generate the sigma-points
  sP = chol(P,'lower'); % Calculate square root of error covariance
  % chi_p = "chi previous" = chi(k-1)
  chi_p = [x(:,k-1), x(:,k-1)*ones(1,n)+sqrt(n+lambda)*sP, ...
  x(:,k-1)*ones(1,n)-sqrt(n+lambda)*sP];
  % Step 2: Prediction Transformation
  % Propagate each sigma-point through prediction
  % chi_m = "chi minus" = chi(k|k-1)
  chi_m = [1 0 dT 0; 0 1 0 dT; 0 0 1 0; 0 0 0 1]*chi_p;
  x_m = chi_m*wm; % Calculate mean of predicted state
  % Calculate covariance of predicted state
  P_m = Q;
  for i = 1:2*n+1
    P_m = P_m + wc(i)*(chi_m(:,i) - x_m)*(chi_m(:,i) - x_m)';
  end
  % Step 3: Observation Transformation
  % Propagate each sigma-point through observation
  psi_m = [sqrt((chi_m(1,:)-xyA(1)).^2 + (chi_m(2,:)-xyA(2)).^2); ...
  sqrt((chi_m(1,:)-xyB(1)).^2 + (chi_m(2,:)-xyB(2)).^2)];
  y_m = psi_m*wm; % Calculate mean of predicted output
  % Calculate covariance of predicted output
  % and cross-covariance between state and output
  Pyy = R;
  Pxy = zeros(n,2);
  for i = 1:2*n+1
    Pyy = Pyy + wc(i)*(psi_m(:,i) - y_m)*(psi_m(:,i) - y_m)';
    Pxy = Pxy + wc(i)*(chi_m(:,i) - x_m)*(psi_m(:,i) - y_m)';
  end
  % Step 4: Measurement Update
  K = Pxy/Pyy; % Calculate Kalman gain
  x(:,k) = x_m + K*(yt(:,k) - y_m); % Update state estimate
  P = P_m - K*Pyy*K'; % Update covariance estimate
end
%% Display results
t = dT*(1:N);
for i = 1:4
  subplot(2,2,i); plot(t,x(i,:),'b-',t,xe(i,:),'g-.',t,xt(i,:),'r-');
  xlabel('Time (s)'); ylabel(['x_',num2str(i)]); grid on; 
  legend('UKF','EKF','True');
end