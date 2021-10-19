% Nonlinear regression for rate of enzyme reaction
% 
% Load the measurements 
obs=[0.038 0.194 0.425 0.626 1.253 2.500 3.740; 
  0.050 0.127 0.094 0.2122 0.2729 0.2665 0.3317]; 
reactants = obs(1,:); 
rate = obs(2,:); 
beta = [0.9 0.2]; 
% Levenberg-Marquardt method
options=statset('Display','iter'); 
[betahat,resid,J] = nlinfit(reactants,rate,@reac,beta,options);

