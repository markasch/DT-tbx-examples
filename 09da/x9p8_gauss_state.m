% Simulate a Guassian random walk.
% initialize
randn('state',123)
R=1; Q=1; K=100;
% simulate
X_init = sqrt(Q)*randn(K,1);
X = cumsum(X_init);
W = sqrt(R)*randn(K,1);
Y = X + W;
% plot
plot(1:K,X,1:K,Y(1:K,1),'ro')
xlabel('k'), ylabel('x_k')