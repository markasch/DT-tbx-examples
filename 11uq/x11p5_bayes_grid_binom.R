# Bayesian inference by Grid Evaluation
# define grid
d <- 80
p_grid <- seq( from=0 , to=1 , length.out=d )
# define 3 priors:
# 1. uniform/flat/uninformative
prior1 <- rep( 1 , d )
# 2. step function
prior2 <- ifelse( p_grid < 0.5 , 0 , 1 )
# 3. peaked
prior3 <- exp( -5*abs( p_grid - 0.5 ) )
# compute likelihood at each point in the grid
likelihood <- dbinom( 6 , size=9 , prob=p_grid )
# compute product of likelihood and prior
unstd.posterior1 <- likelihood * prior1
unstd.posterior2 <- likelihood * prior2
unstd.posterior3 <- likelihood * prior3
# normalize the posterior, so it sums to 1
posterior1 <- unstd.posterior1 / sum(unstd.posterior1)
# plot the posteriors
plot( p_grid , posterior1 , type="b" ,
      xlab="probability of data" , ylab="posterior probability" )
mtext( "80 points-flat prior" )
#
posterior2 <- unstd.posterior2 / sum(unstd.posterior2)
plot( p_grid , posterior2 , type="b" ,
      xlab="probability of data" , ylab="posterior probability" )
mtext( "80 points-step prior" )
#
posterior3 <- unstd.posterior3 / sum(unstd.posterior3)
plot( p_grid , posterior3 , type="b" ,
      xlab="probability of data" , ylab="posterior probability" )
mtext( "80 points-peaked prior" )