##### ---- Parameters and data
sig2 <- 1 ;          # likelihood variance
tau2 <- 10 ; mu <- 5 # Gaussian prior
theta<-0 ; delta<-2  # Gaussian proposal
N<-10000 ;           # length of chain
# observations
y<-c(9.37, 10.18, 9.16, 11.60, 10.33)
# initialize
THETA<-NULL ; set.seed(1)
#### ---- Metropolis algorithm
for(n in 1:N)
{
  theta.star <- rnorm(1,theta,sqrt(delta)) # candidate
  log.r <-( sum(dnorm(y,theta.star,sqrt(sig2),log=TRUE)) +
           dnorm(theta.star,mu,sqrt(tau2),log=TRUE) ) -
          ( sum(dnorm(y,theta,sqrt(sig2),log=TRUE)) +
           dnorm(theta,mu,sqrt(tau2),log=TRUE) )
  # acceptance test
  if(log(runif(1))<log.r) { theta<-theta.star }
  THETA<-c(THETA,theta) # update chain
}
##### ---- Output
# compute exxact sample posterior
n<-5
y<-round(rnorm(n,10,1),2)
mu.n<-( mean(y)*n/sig2 + mu/tau2 )/( n/sig2+1/tau2) 
tau2.n<-1/(n/sig2+1/tau2)
# set up
par(mfrow=c(1,2))
Nkeep<-seq(10,N,by=10)
# plot Markov Chain
plot(Nkeep,THETA[Nkeep],type="l",xlab="iteration",ylab=expression(theta))
# plot histogram of posterior distribution
hist(THETA[-(1:50)],prob=TRUE,main="",xlab=expression(theta),ylab="density")
th<-seq(min(THETA),max(THETA),length=100)
lines(th,dnorm(th,mu.n,sqrt(tau2.n)),col='blue' )
