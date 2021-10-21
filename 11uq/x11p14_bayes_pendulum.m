%--- BayesPendulum.m ---%
% Bayesian parameter estimation for noisy pendulum
clear
t     = 1:100;
theta = 0.2; % True value of parameter.
std   = 1; % standard deviation of noise
xt    = sin(theta*t) + std*randn(1,100); % Noisy measurements
figure
subplot(1,3,1); plot(t,xt); axis square
% Prior belief
ptheta = [0.275 0.15 0.275 0.025 0.05 0.225];
thetav = [0 .1 .2 .3 .4 .5]; % Values of theta
fac= 1/sqrt(2*pi*std^2);
% Compute posterior = prior x likelihood
for i = 1:6
   pr(i) = ptheta(i)*fac*prod(exp(-((xt-sin(thetav(i)*t)).^2)/(2*std^2)));
end
% Normalize
prs = pr./sum(pr);
% Plot results
subplot(1,3,2);bar(0:0.1:0.5,ptheta); xlim([-0.1 0.6]); ; axis square
subplot(1,3,3);bar(0:0.1:0.5,prs);    xlim([-0.1 0.6]); ; axis square
