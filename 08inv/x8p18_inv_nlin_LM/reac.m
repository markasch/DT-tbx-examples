function yhat = reac(beta,x) 
%REAC simple model for reaction kinetics. 
%   YHAT = REAC(BETA,X) gives the predicted values of 
%   reaction rate, YHAT, as a function of the vector of  
%   parameters, BETA, and the vector of data, X. 
% 
%   The model form is: 
%   y = (b1*x)./(b2 + x)
b1 = beta(1); b2 = beta(2);
yhat = (b1*x)./(b2 + x);