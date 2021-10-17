% ODE with instabilities:
%
%   y' = -1/t^2 + alpha*(y - 1/t), alpha>0, 1<t<10,
%   y(1) = 1
%
% Exact solution:  y = 1/t
%
% Compare the following methods:
% 1. RK4 - Explicit Runge-Kutta fourth order.
% 2. RKF45 - Explicit Adaptive Runge-Kutta-Fehlberg (use ode45)
% 3. BDF1 - Implicit Backward Euler with fixed time step
% 4. BDFn - Implicit Adaptive BDF (use ode15s)
clear all
% Define the initial conditions.
t0 = 1; tf = 10; y0 = 1;
% Solve the differential equation by explicit methods.
% 1. RKF45
[t1,y1] = ode45(@ode_f,[t0 tf],y0);
% 2. RK4
Nt = 100;
[t2,y2] = x3p25_rk4(@ode_f,[t0 tf],y0,Nt);
% Exact solution
ye1 = 1./t1; ye2 = 1./t2;
% Plot solutions
figure
subplot(1,2,2)
plot(t1,y1,t1,ye1), legend('numerical','exact')
xlabel('t'), ylabel('y(t)')
subplot(1,2,1)
plot(t2,y2,t2,ye2), legend('numerical','exact')
xlabel('t'), ylabel('y(t)')
% Solve the differential equation by implicit methods.
% 3. Backward Euler with fixed time step
Nh = 18; % h=0.5 gives better accuracy
[t3,y3] = beuler(@ode_f,[t0 tf],y0,Nh); % use fsolve
% 4. Adaptive step, adaptive order ode15s
[t4,y4] = ode15s(@ode_f,[t0 tf],y0);
% Plot implicit solutions
figure
subplot(1,2,2)
plot(t4,y4,t4,1./t4), legend('numerical','exact')
xlabel('t'), ylabel('y(t)')
subplot(1,2,1)
plot(t3,y3,t3,1./t3), legend('numerical','exact')
xlabel('t'), ylabel('y(t)')