library(mvtnorm)
data <- read.table("eBayNumberOfBidderData.dat", sep = "" , header = T , na.strings ="", stringsAsFactors= F)
data = as.data.frame(data)

col_names = colnames(data)

# lab 3.2.1

model <- glm(nBids ~ . , data[,-2], family = poisson(link = "log"))
model
summary(model)

# lab 3.2.2
Y = data$nBids
X = as.matrix(data[,-1])
n = dim(data)[1]
mu = rep(0,9)


Sigma = 100*solve(t(X)%*%X)
beta_prior = rmvnorm(1,mean=mu,sigma=Sigma)
dim(beta_prior)


LogPostPoisson <- function(betas,Y,X,mu,Sigma){
  linPred <- betas%*%t(X);
  # finding the log likelihood 
  logLik <- sum(Y*linPred - exp(linPred) - log(factorial(Y)));
  #print(logLik)
  # finding the log logPrior using the parameters given
  logPrior <- dmvnorm(betas, mu, Sigma, log=TRUE);
  #print(logPrior)
  # returning the log posterior as the sum of loglikihood and logPrior
  return(logLik + logPrior)
} 
OptimRes <- optim(beta_prior,LogPostPoisson,gr=NULL,Y,X,mu,Sigma,method=c("BFGS"),control=list(fnscale=-1),hessian=TRUE)

# Getting the hessian matrix from the results returned from the optim function
beta_mode_hessian = OptimRes$hessian
# Getting the mode β values as the parameters returned from the optim function
beta_mode = OptimRes$par
colnames(beta_mode) = col_names[-1]
# Getting the J(beta_mode) as the negative of the hessian returned 
J = -beta_mode_hessian
# Calculating inverse of J
J_inverse = solve(-beta_mode_hessian)

posterior_beta_draw = rmvnorm(1,beta_mode, J_inverse)

# lab 3.2.3
num_random_walk = 1E4
# Function for metropolis random walk
MetropolisRandomWalk <- function(beta_initial, c, J_inverse, log_density_function, mu, Sigma, Y, X){
  metropolis_betas = matrix(data = NA, nrow = num_random_walk, ncol = dim(X)[2])
  # Initializing previous betas as the initial value
  prev_beta = beta_initial
  for(i in 1:num_random_walk){
    # Sampling betas from the multivariate normal distribution 
    sample = rmvnorm(1, mean = prev_beta, sigma = c * J_inverse)
    # Calculating the log of the density values for the sampled betas
    log_density_sample = log_density_function(sample, Y, X, mu, Sigma)
    # Calculating the log of the density values for the previous samples of accepted betas
    log_density_prev_sample = log_density_function(prev_beta, Y, X, mu, Sigma) 
    # Calculating alpha values
    alpha = min(1, exp((log_density_sample - log_density_prev_sample)))
    # Generating a random number between 0 and 1
    rand = runif(n=1, min = 0, max = 1)
    # Accepting sample based on the probability alpha
    if(rand<alpha){
      metropolis_betas[i,] = sample
      prev_beta = sample
    }
    else{
      metropolis_betas[i,] = prev_beta
    }
    
  }
  return(metropolis_betas)
}
# Sampling betas based using Metropolis Random Walk
metropolis_posterior_betas = MetropolisRandomWalk(posterior_beta_draw, .6, J_inverse, LogPostPoisson, mu, Sigma, Y, X)
colnames(metropolis_posterior_betas) = col_names[-1]

# plotting Sampled value to access the convergence of betas
par(mfrow=c(3,3))
for(i in 1:dim(X)[2]){
  plot(metropolis_posterior_betas[,i], type = "l", ylim = c(0,2*model$coefficients[i]),ylab = col_names[1+i])
  abline(h=model$coefficients[i], col="red")
}

# lab 3.2.4
par(mfrow=c(1,1))
# Creation data for test case 
sample = matrix(c(1, 1, 0, 1, 0, 1,0,1.2,0.8))
# Getting predictions for the test data based on sampled betas as the mean of the Poisson Distribution is parameters passed it
# The plot of the parameters passed to it (exp(betas*X)) would be good plot of the predictive distribution
prediction = exp(metropolis_posterior_betas %*% sample)

#print(exp(model$coefficients%*%sample))
hist(prediction, main = "Predictive distribution", sub =paste("glm predicts: ",exp(model$coefficients%*%sample)))

# Calculating probability of no bidders
prob_no_bids = mean(rpois(num_random_walk,prediction)==0)
prob_no_bids



