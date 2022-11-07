
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
