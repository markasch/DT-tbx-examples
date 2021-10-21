#--- BayesBetaBinomEpidemic.R ---#
# Compute some prior probabilities from a beta(a,b) law
a<-2 ; b<-20
a/(a+b)
(a-1)/(a-1+b-1)
pbeta(.20,a,b) - pbeta(.05,a,b)
pbeta(.10,a,b)
# Bayes estimation:
par(mar=c(3,3,1,1),mgp=c(1.75,.75,0))
par(mfrow=c(1,2))
# Likelihood for Y=0
dbinom(0,20,.05)
n<-20
x<-0:n
delta<-.3
# plot three binomial likelihoods
plot( range(x-delta), c(0,.4),xlab="number infected",
     ylab="probability",type="n")
points( x-delta,dbinom(x,n,.05),type="h",col=1,lwd=3)
points( x,dbinom(x,n,.10),type="h",col=2,lwd=3)
points( x+delta,dbinom(x,n,.20),type="h",col=3,lwd=3)
legend(10,.35,legend=c(
    expression(paste(theta,"=0.05",sep="")), 
    expression(paste(theta,"=0.10",sep="")),
    expression(paste(theta,"=0.20",sep="")) ),
     lwd=c(3,3,3), 
    col=c(1:3) ,bty="n") 
# Parameters for the conjugate prior and posterior
a<-2 ; b<-20
y<-0 ; n<-20
# Plot the posterior
theta<-seq(0,1,length=500)
plot(theta, dbeta(theta,a+y,b+n-y), col=2,
     type="l",
     xlab="proportion infected in population",
     ylab="", lwd=2, ylim=c(0,16)
     )
lines(theta, dbeta(theta,a,b),col="blue",lwd=2) # prior
legend(.5,14,legend=c( expression(paste(italic("p"),"(",theta,")",sep="")), 
     expression(paste(italic("p"),"(",theta,"|",italic("y"),")",sep=""))  ), 
  bty="n", lwd=c(2,2),col=c("blue",2))
# Compute posterior probabilities
(a+y)/(b+n-y)
(a+y-1)/(a+y-1+b+n-y-1)
pbeta(.20,a+y,b+n-y) - pbeta(.05,a+y,b+n-y)
pbeta(.10,a+y,b+n-y)
#
#Sensitivity analysis
#
par(mar=c(3,3,1,1),mgp=c(1.75,.75,0))
par(mfrow=c(1,2))
g<-50
th0<-seq(0.01,.5,length=g)
nu0<-seq(1,25,length=g) 
# Loop and compute contours
PP10<-PM<-PLQ<-PUQ<-matrix(0,g,g)
for(i in 1:g) 
  {
  for(j in 1:g) 
    {
      a<-nu0[i]*th0[j]
      b<-nu0[i]*(1-th0[j]) 
   
     PM[i,j]<- (a+y)/(a+y+b+n-y) 
     PP10[i,j]<- pbeta(.10,a+y,b+n-y)
     PLQ[i,j]<- qbeta(.05,a+y,b+n-y) 
     PUQ[i,j]<- qbeta(.95,a+y,b+n-y) 
    } 
}
cols = rev(colorRampPalette(c('darkred','red','blue','lightblue'))(14))
contour(nu0,th0,PM,xlab=expression(italic(w)), ylab=expression(theta[0]),col=cols)
cols = rev(colorRampPalette(c('darkred','red','blue','lightblue'))(8))
contour(nu0,th0,PP10,xlab=expression(italic(w)),col=cols, 
   levels=c(0.1,0.3,.5,.70,.90,.975) )
#### Adjusted Wald 95% interval of confidence
a<-2 ; b<-2 
th<-  (y+a)/(n+a+b)
th+c(-1,1)*1.96*sqrt(th*(1-th)/n)

qbeta(c(.025,.975),a+y,b+n-y)
