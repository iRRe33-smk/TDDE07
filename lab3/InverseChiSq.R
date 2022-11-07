library(LaplacesDemon)
x <- dinvchisq(1,1,1)
x <- rinvchisq(10,1)

sigma_sq_0 = 1.8
v0 = 100

#Plot Probability Functions
x <- seq(from=0.01, to=5, by=0.01)
plot(x, dinvchisq(x,v0,sigma_sq_0), type="l", main="Probability Function",
     ylab="density", col="red")

x[which.max(dinvchisq(x,v0,sigma_sq_0))]
abline(v=x[which.max(dinvchisq(x,v0,sigma_sq_0))])
