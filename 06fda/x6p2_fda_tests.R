# fda examples
library(fda)
# spline and Fourier basis functions
spline.basis <- create.bspline.basis(rangeval=c(0,10), nbasis=5)
plot(spline.basis, lty=1, lwd=2)
fourier.basis=create.fourier.basis(rangeval=c(0,10), nbasis=5)
plot(fourier.basis, lty=1, lwd=2)
# smooth interpolation of Brownian motion
set.seed(1234)
BM <- cumsum(rnorm(10000)) # BM trajectory on [0,N], N=10ˆ4
plot.ts(BM, xlab="", ylab="")
B25.basis <- create.bspline.basis(rangeval=c(0,10000), nbasis=25)
BM.fd     <- smooth.basis(y=BM, fdParobj=B25.basis)
lines(BM.fd, lwd=3,col=2)
# generate an ensemble of BM trajectories
# Generate an ensemble of BM trajectories
N <- 100
BM.mat <- matrix(0, ncol=N, nrow=10000)
for(n in 1:N){
	BM.mat[, n] <- cumsum(rnorm(10000))/100
    }
B25.basis <- create.bspline.basis(rangeval=c(0,10000), nbasis=25)
BM.fd <- smooth.basis(y=BM.mat, fdParobj=B25.basis)
# Compute functional statistics
BW.mean <- mean.fd(BM.fd$fd)
BM.sd   <- std.fd(BM.fd$fd)
BM.cov  <- var.fd(BM.fd$fd)
lines(BM.sd, lwd=3, col=3); lines(BM.mean, lty=2, lwd=3, col=4)
# compute PCA
# compute PCA
BM.pca <- pca.fd(BM.fd$fd, nharm=4)
plot(BM.pca$harmonics)
BM.pca$varprop