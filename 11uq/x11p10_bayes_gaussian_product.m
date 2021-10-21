%--- BayesGaussianProduct.m ---%
% Bayes posterior as product of two Gaussians
prior_mean = 0; prior_sd   = sqrt(1.21);
obs_mean   = 2; obs_err_sd = sqrt(0.64);
n = 100;
% Get the prior and observational error variance
prior_var = prior_sd^2;
obs_err_var = obs_err_sd^2;
% Compute the posterior variance
if (n==2 )
  post_var = 1. / (1. / prior_var + 2. / obs_err_var) % for n=2 obs
else
  post_var = 1. / (1./prior_var + 1./obs_err_var)  
end
post_sd = sqrt(post_var);
% Compute the posterior mean
if (n==2 )
  K=1./(obs_err_var/(2*prior_var)+1); % for n=2 obs
  post_mean = prior_mean + K*(obs_mean-prior_mean)% for n=2 obs
else
  post_mean = post_var * (prior_mean/prior_var + obs_mean/obs_err_var)
end
% plot
figure; hold on;
p_prior = plot_gaussian(prior_mean, prior_sd,2,'--','r');
p_likel = plot_gaussian(obs_mean, obs_err_sd,2,'-.','g');
p_poste = plot_gaussian(post_mean, post_sd,2,'-','b');
hold off
legend('Prior','Likelihood','Posterior')