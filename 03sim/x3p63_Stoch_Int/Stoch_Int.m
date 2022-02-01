%--- StochInt.m ---%
% Approximations of Stochastic Integrals
%
% Ito and Stratonovich integrals for W.dW
%
randn('state',123) % reference state
% time interval, time step                
T = 1; N = 500; dt = T/N;
% paths
dW = sqrt(dt)*randn(1,N);               % increments
W  = cumsum(dW);                        % cumulative sum
% approximations of stochastic integrals
ito   = sum([0,W(1:end-1)].*dW)  
strat = sum((0.5*([0,W(1:end-1)]+W) + 0.5*sqrt(dt)*randn(1,N)).*dW) 
% approximation errors
itoerr   = abs(ito   - 0.5*(W(end)^2-T))
straterr = abs(strat - 0.5*W(end)^2)
