# Simple Linear  Regression  
x <- c(7, 3, 4, 6, 10, 9)
y <- c(276, 43, 82, 136, 417, 269)
SLRmodel <- lm(y~x)
SLRmodel
coefs <- coef(SLRmodel)
# Plot
plot(x, y, pch=20,col="red", xlab="Number ", ylab="Time (sec.)")
abline(coefs[1],coefs[2])
# Diagnostics
summary(SLRmodel)