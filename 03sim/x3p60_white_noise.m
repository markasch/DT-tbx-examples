%--- white_noise.m ---%
randn('state',0);   % ensure reproducibility
T = 0:0.001:1;      % define time interval
X = randn(size(T)); % Gaussian N(0,1) process
plot(T,X);  xlabel('Time, t'); ylabel('w(t)')
    
