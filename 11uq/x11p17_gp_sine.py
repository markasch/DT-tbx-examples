#--- GP-sine.py ---#
import numpy as np
import matplotlib.pyplot as pl
# The true unknown function we are trying to approximate
f = lambda x: np.sin(0.9*x).flatten()
# Define the SE kernel function
def kernel(x, xp, ell):
    # GP squared exponential kernel
    sqdist = np.sum(x**2,1).reshape(-1,1) + 
             np.sum(xp**2,1) - 2*np.dot(x, xp.T)
    return np.exp(-.5 * (1/ell) * sqdist)
# Set parameters for regression
N = 5          # number of training points.
n = 50         # number of test points.
s = 0.00005    # noise variance.
# Sample input points and get noisy, training observations
Xtrain = np.random.uniform(-5, 5, size=(N,1))
ytrain = f(Xtrain) + s*np.random.randn(N)
# Compute the covariance (kernel) and its Cholesky factorization
param = 0.1
K = kernel(Xtrain, Xtrain, param)
L = np.linalg.cholesky(K + s*np.eye(N))
# Choose test points for posterior predictions
Xtest = np.linspace(-5, 5, n).reshape(-1,1)
# Compute posterior mean at test points.
Ks  = kernel(Xtrain, Xtest, param)
LKs = np.linalg.solve(L, Ks)
mu  = np.dot(LKs.T, np.linalg.solve(L, ytrain)).reshape((n,))
# Compute posterior variance at test points.
Kss = kernel(Xtest, Xtest, param)
s2  = np.diag(Kss) - np.sum(LKs**2, axis=0)
s   = np.sqrt(s2)
# PLOTS:
pl.figure(1)
pl.clf()
pl.plot(Xtrain, ytrain, 'bs', ms=10)
pl.plot(Xtest, f(Xtest), 'b-')
pl.gca().fill_between(Xtest.flat, mu-2*s, mu+2*s, color="#dddddd")
pl.plot(Xtest, mu, 'r--', lw=2)
pl.title('Mean predictions plus 2 std.deviations')
pl.axis([-5, 5, -3, 3])
# Draw samples from the posterior at test points.
L = np.linalg.cholesky(Kss + 1e-6*np.eye(n) - np.dot(LKs.T, LKs))
f_post = mu.reshape(-1,1) + np.dot(L, np.random.normal(size=(n,4)))
pl.figure(2)
pl.clf()
pl.plot(Xtest, f_post)
pl.plot(Xtrain, ytrain, 'bs', ms=10)
pl.plot(Xtest, f(Xtest), 'b-')
pl.gca().fill_between(Xtest.flat, mu-2*s, mu+2*s, color="#dddddd")
pl.plot(Xtest, mu, 'r--', lw=2)
#pl.title('4 samples from the GP posterior')
pl.axis([-5, 5, -3, 3])