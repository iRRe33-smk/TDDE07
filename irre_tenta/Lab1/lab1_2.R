y = c(33,24,48,32,55,74,23,76, 17)
mu = 3.5
n = 9

# 2.a)

tau_sq = sum((log(y)-mu)^2)/n   
X_draws = rchisq(1E5,n-1)

sigma_sq = (n-1)*tau_sq / X_draws   # Simulating draws from Inv−χ2(n,τ) distribution

# Theoretical density function for Inv−χ2(n,τ) distribution 
inversechisq <- function(x, v, s){
  density = (((v/2)^(v/2))/gamma(v/2)) * (s^v) * (x^(-(v/2+1))) * (exp((-v*s^2)/(2*x)))
  return (density)
}
# Plotting posterior draws from sampling 
plot(density(sigma_sq), lwd=3, main="Simulated and theoretical posteriors for sigma^2") 
# Plotting posterior distribution using therotical density function 
lines(seq(0.1,15,by=0.1), inversechisq(seq(0.1,15,by=0.1), n, sqrt(tau_sq)), type = "l", col =  "red", lty=2, lwd=3)  
legend("topright", inset=.02, c("Simulated Posterior Draws", "Theoretical Posterior Density"), fill=c("black", "red"))

# 2.b)
# Calculating and plotting Gini coefficient from posterior Draws
distribution_gini = 2 * pnorm(sqrt(sigma_sq/2)) - 1
plot(1:length(distribution_gini), sort(distribution_gini), type = "l")
plot(density(distribution_gini), main = "Posterior Distribution Gini Coefficent") 

# 2.c)
quants = quantile(distribution_gini, c(0.025,0.975))
quants
abline(v=quants,col="red")


# 2.d)

# Finding density of gini distribution
dense = density(sorted)
# storing coordinate and desnity values in a dataframe
df = data.frame(x=dense$x,y=dense$y)
# sorting dataframe based on density values
df = df[order(df$y),]
# Finding cumulative sum of sorted density values
df$cumsum = cumsum(df$y)
# finding all points at which the cumulative sum is > 5% if the total density ( 95% Highest density posterior interval)
df = df[df$cumsum>.05 * sum(df$y),c("x","y","cumsum")]
head(df)

# verifying if we have found the correct interval
sum(df$y) / sum(dense$y)

# plotting the HDPI
points(df,col="red")
legend("topright", inset=.01, c("95% Equal Tail Confidence Intervel", "Highest Density Posterior Intervel"), fill=c("green", "red"))
#black lines indicate equal tail credible interval 95%
#red dots indicate 95% highest posterior density

