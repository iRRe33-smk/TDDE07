
library(manipulate)

sucess <- 13
n_trails <- 50
failures <- 37

alpha <- 5
beta <- 5

# Task 1. a)

# Function to calculate and plot mean of posterior Distribution Random Draws
PosteriorDraws_MeanPlot <- function(alpha, beta, n, sucess_proportions){
  mean_posterior <-  list()
  for(i in 1:n){
    posterior_draws = rbeta(i, alpha+sucess, beta+failures)
    mean_posterior[i] <- mean(posterior_draws)
  }
  print(mean_posterior[n])
  plot(c(1:n), mean_posterior, type = 'l', lwd = 1, col = "blue", xlab = "number of Random Draws", 
       ylab = 'mean of Draws', main = 'Mean plot Posteriors')
}

# Function to calculate and plot SD of posterior Distribution Random Draws
PosteriorDraws_SDPlot <- function(alpha, beta, n, sucess_proportions){
  sd_posterior <- list()
  for(i in 1:n){
    posterior_draws = rbeta(i, alpha+sucess, beta+failures)
    sd_posterior[i] <- sd(posterior_draws)
  }
  print(sd_posterior[n])
  plot(c(1:n), sd_posterior, type = 'l', lwd = 1, col = "red", xlab = "number of Random Draws",  ylab = 'SD of Draws', main = 'SD plot Posteriors')
  
}

manipulate(
  PosteriorDraws_MeanPlot(alpha, beta, n, sucess/n_trails),
  n = slider(1, 100000, step=100, initial = 10000, label = "Number of trails(Random Draws)")
)

# True mean of beta distribution -> alpha/(alpha+beta), so true mean of posterior
true_mean = (alpha+sucess)/((alpha+sucess)+(beta+failures))
legend("bottomright", inset=.02, c("Mean posterior Draws", "Theoretical Mean"), fill=c("blue", "green"))
abline(h=true_mean, col = "green", type="l", lty=2)
print(true_mean)

manipulate(
  PosteriorDraws_SDPlot(alpha, beta, n, sucess/n_trails),
  n = slider(1, 100000, step=100, initial = 10000, label = "Number of trails(Random Draws)")
)

# True variance of beta distribution -> (alpha*beta)/((alpha+beta+1) * (alpha + beta)^2)
true_sd = sqrt(((alpha+sucess)*(beta+failures))/((alpha+sucess+beta+failures+1) * (alpha + sucess + beta + failures)^2))
print(true_sd)
legend("bottomright", inset=.02, c("Standard Deviation posterior Draws", "Theoretical Standard Deviation"), fill=c("red", "green"))
abline(h=true_sd, col = "green", type="l", lty=2)

# Task 1. b)
posterior_draws <-  rbeta(10000, alpha+sucess, beta+failures)              # taking 10000 random draws from the posterior beta distribution
count_density_less <- length(posterior_draws[posterior_draws<0.3])         # counting number of draws with value less than 0.3
posterior_probablity_less <- count_density_less/10000                    # calculating probability of draws with value < 0.3  
posterior_probablity_less_exact <- pbeta(0.3, alpha+sucess, beta+failures) # using pbeta to find the same probability using cumulative frequency of beta distribution
posterior_probablity_less_exact
posterior_probablity_less

# Task 1. c)
posterior_draws <-  rbeta(10000, alpha+sucess, beta+failures)              # taking 10000 random draws from the posterior beta distribution
posterior_logodds <- log(posterior_draws/(1-posterior_draws))              # converting to log odds using simulated draws from beta posterior
hist(posterior_logodds, freq = FALSE)                                                    # plotting posterior distribution of log-odds
lines(density(posterior_logodds))                                           # plotting posterior distribution of log-odds





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
num = length(distribution_gini)
sorted = sort(distribution_gini)

# Finding cumulative sum of densities 
c_sum = cumsum(sorted)
total = c_sum[length(c_sum)]
# Normalizing
c_sum_norm <- c_sum / total


lower <- 0
upper <- 0
N = length(c_sum)
for (i in 1:as.integer(N/2)){
  # Finding lower bound of 95% equal tail credible interval as the point which contains 2.5% of the distribution 
  if (c_sum_norm[i] < 0.025){
    lower <- i
  }
  # Finding upper bound of 95% equal tail credible interval as the point which contains 97.5% of the distribution 
  if (c_sum_norm[N-i] > 0.975){
    upper <- N-i
  }
  
  
}

# Plotting Confidence Intervals 
plot(density(distribution_gini), main = "95% Equal Trail Confidence Intervel - Posterior Distribution Gini Coefficent") 
abline(v = sorted[lower], col = "green")
abline(v = sorted[upper], col = "green")
legend("topright", inset=.02, c("Confidence Intervel"), fill=c("green"))


# Checking if found confidence interval is correct
sum(sorted[lower:upper]/sum(sorted))
sum(sorted[1:lower]/sum(sorted))
sum(sorted[upper:N]/sum(sorted))

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






#3.a

# Function for the posterior distribution with n data points
p_cond <- function(y,mu,k){
  n = length(y)
  res = (1/(2*pi*besselI(k,0))^n) * exp(sum(k*cos(y-mu))-k)
  return(res) 
}



y= c(1.83, 2.02, 2.33, -2.79, 2.07, 2.02, -2.44, 2.14, 2.54, 2.23)

# Simulating 15000 posterior distribution over a k grid of 0-10 
N=1500
res = c()
used_k = c()
for (k in seq(0,10*N)){
  res = c(res,p_cond(y,2.51,k/N))
  used_k = c(used_k,k/N)
}

# Normalizing density values by dividing by total sum of densities
res_norm = res/sum(res)
# Plotting k values against densities
plot(used_k,res_norm, main = "Posterior Distribution of k Values")

# The mode of a distribution with a discrete random variable is the value of the term that occurs the most often
# i.e the value of k with the highest probability density
used_k[res==max(res)]
abline(v=used_k[res==max(res)], col = "red", type="l", lty=2)
legend("topright", inset=.01, c("Posterior mode of k"), fill=c("red"), )

#dense = density(res)







