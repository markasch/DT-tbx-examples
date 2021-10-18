%SDE_OU  Euler-Maruyama method for Ornstein-Uhlenbeck SDE
% Simulate trajectories from the OU process
% parameters 
lambda = 0.5; sigma2 = 1;
dt = 0.01;
T = (0:dt:1);
x0 = 4;
% mean and variance
M = exp(-lambda*T)*x0;
P = sigma2/(2*lambda)*(1 - exp(-2*lambda*T));
% initialize for 50 realizations
XX = zeros(50,length(T));
for n=1:size(XX,1)
   x = x0;
   for k=1:length(T)
      XX(n,k) = x;
      x = x - lambda * x * dt + sqrt(dt)*randn;
   end
end
figure(1); clf; hold on
% Shade the 95% quantiles
upper = M + 1.96*sqrt(P);
lower = M - 1.96*sqrt(P);
gray  = [.9 .9 .9];
fill([T fliplr(T)],[upper fliplr(lower)],gray,'EdgeColor',gray)
% Plot realizations
h1 = plot(T,XX,'-','Color',[.5 .5 .5],'LineWidth',0.5);
% Plot mean and quantiles
h2  = plot(T,M,'b-','LineWidth',1);
h34 = plot(T,upper,'--r',T,lower,'--r');
legend([h2(1) h34(1) h1(1)],'Mean','95% Quantiles','Realizations')
xlabel('Time, t'); ylabel('x(t)'); ylim([0 5]);
