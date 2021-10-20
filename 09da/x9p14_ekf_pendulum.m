%--- EKfPendulum.m ---%
% Extended Kalman Filter for noisy pendulum
%---------------------%
% Parameters
DT = 0.01;
g  = 9.81; L=1;
Q  = 0.01*[DT^3/3 DT^2/2; DT^2/2 DT];
R  = 0.1;
% Initialize
m0 = [1.8;0]; 
P0 = 0.1*eye(2);
Qsqrt = chol(Q,'lower');
steps = 500;
T = []; X = []; Y = [];
t = 0;
x = [1.5;0];
for k=1:steps
  x = [x(1)+x(2)*DT;
  x(2)-(g/L)*sin(x(1))*DT];
  w = Qsqrt * randn(2,1);
  x = x + w;
  y = sin(x(1)) + sqrt(R)*randn;
  t = t + DT;
  X = [X x];  Y = [Y y]; T = [T t];
end
% Plot the data
plot(T,Y,'g.',T,X(1,:),'r-');
pause
% EKF Filter
m = m0;
P = P0;
MM = zeros(size(m,1),length(Y));
PP = zeros(size(P,1),size(P,2),length(Y));
for k=1:length(Y)
  M  = [m(1)+m(2)*DT; m(2)-g*sin(m(1))*DT];
  Mx = [1 DT; -g*cos(m(1))*DT 1];
  m  = Mx;
  P  = Mx*P*Mx' + Q;

  h  = sin(m(1));
  Hx = [cos(m(1)) 0];
  S  = Hx*P*Hx' + R;
  K = P*Hx'/S;
  m = m + K*(Y(k) - h);
  P = P - K*S*K';

  MM(:,k) = m;
  PP(:,:,k) = P;
end
% Plot the filtering result
clf;
plot(T,X(1,:),'b',T,Y,'ro',T,MM(1,:),'g');
xlabel('Time t');
ylabel('Angle x_{1,k}') 