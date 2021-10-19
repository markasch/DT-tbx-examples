close all
clear
clc

global gamma epsilon GAM OMEG

gamma=0.05;
epsilon=1.0;
OMEG=1.0;

GAM=7.5;

[t0 x0]=ode45(@duffing,[0 50],[3 4]);

figure(1)
plot(t0,x0(:,1),'b','LineWidth',2)
%axis tight
title('Duffing''s equation with 0.03% perturbation')
xlabel('t'), ylabel('x(t)')

% perturbation 0.03%

[t1 x1]=ode45(@duffing,[0 50],[3.01 4.01]);

hold on
plot(t1,x1(:,1),'--r','LineWidth',2)
legend('unperturbed', '0.03% perturbation')
hold off


figure(2)
plot(t0,x0(:,1),'b','LineWidth',2)
%axis tight
title('Duffing''s equation with 0.06% perturbation')
xlabel('t'), ylabel('x(t)')

% perturbation 0.06%

[t1 x1]=ode45(@duffing,[0 50],[3.02 4.02]);

hold on
plot(t1,x1(:,1),'--r','LineWidth',2)
legend('unperturbed', '0.06% perturbation')
hold off
