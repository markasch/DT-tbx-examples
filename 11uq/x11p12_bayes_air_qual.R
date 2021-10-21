## New York Air Quality Measurements
help("airquality")
## Estimating average wind speeds
wind = airquality$Wind
## Plot histogram and compute statistics
hist(wind, col = "gray", border = "white", xlab = "Wind Speed (MPH)")
n = length(wind)
ybar = mean(wind)
ybar
# >> [1] 9.957516 ## "frequentist" estimate
tau = 1/sd(wind)
## Based on research, we believe tha average wind speeds are 
## closer to 12 mph, but probably no greater than 15.
## Then a potential prior would be N(12, 2)
a = 12
b = 2
## The posterior mean and standard deviation are
postmean = 1/(1 + n*tau) * a + n*tau/(1 + n*tau) * ybar
postsd = 1/(1 + n*tau)
## Sample from the posterior
set.seed(123)
posterior_sample = rnorm(n = 10000, mean = postmean, sd = postsd)
hist(posterior_sample, col = "gray", border = "white", xlab = "Wind Speed (MPH)")
abline(v = median(posterior_sample))
abline(v = ybar, lty = 3)
## Compute posterior, Bayesian statistics 
median(posterior_sample)
# >> [1] 10.00324
## confidence intervals
quantile(x = posterior_sample, probs = c(0.025, 0.975)) 
# >> 2.5%     97.5% 
#  >> 9.958984 10.047404 