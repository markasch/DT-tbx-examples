#--- BayesRegDiabetes.R ---#
#### Diabetes example
load("diabetes.RData")
# Response variable
yf<-diabetes$y
# Center and scale
yf<-(yf-mean(yf))/sd(yf)
# Covariates
Xf<-diabetes$X
# Center and scale
Xf<-t( (t(Xf)-apply(Xf,2,mean))/apply(Xf,2,sd))
## set up training and test data
n<-length(yf)
set.seed(1)
# random indices
i.test  <- sample(1:n,100)
i.train <- (1:n)[-i.test]
y<-yf[i.train] ; y.test <- yf[i.test]
X<-Xf[i.train,]; X.test <- Xf[i.test,]
# fit an LS regression
olsfit <- lm(y~-1+X)
y.test.ols <- X.test%*%olsfit$coef
# plot coefficients
plot(olsfit$coef,type="h",lwd=2,xlab="regression coefficient index",ylab=expression(hat(beta)[LS]))
# plot test values versus their LS estimation
plot(y.test,y.test.ols,pch=1,col="blue",xlab=expression(italic(y)[test]),
ylab=expression(hat(italic(y))[test])) ; abline(0,1)
# compute the mean prediction error
sprintf(' Mean prediction error  = %4.3f',mean( (y.test-y.test.ols )^2 ))