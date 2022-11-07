library(mvtnorm)

# Lab 2.2.a)

# importing data 
data <- read.table("WomenAtWork.dat", sep = "" , header = T , na.strings ="", stringsAsFactors= F)
col_names = colnames(data)
# Finding the number of covariates in the data the -1 is done as the data imported contains the response y
n_coefs <- dim(data)[2]-1

# Initializing parameters of prior
tau = 5
mu = rep(0,n_coefs)
Sigma = tau*tau*diag(n_coefs)

# sampling one element from the prior of β to use as initial value of β when using optim
prior = rmvnorm(1,mean=mu,sigma=Sigma)

# initializing y as a matrix of the response variable from the data set
y = as.matrix(data$Work)
# initializing x as a matrix the covariates from the data set
X = as.matrix(data[,-1])

# Function to get the log posterior of the logistic regression 
LogPostLogistic <- function(betas,y,X,mu,Sigma){
  
  linPred <- X%*%betas;
  # finding the log likelihood 
  logLik <- sum( linPred*y - log(1 + exp(linPred)));
  # finding the log logPrior using the parameters given
  logPrior <- dmvnorm(betas, mu, Sigma, log=TRUE);
  
  # returning the log posterior as the sum of loglikihood and logPrior
  return(logLik + logPrior)
}

# Using the optim function to find the values of β for which the logposterior is maximum. 
# This would give us the mode of β's in the posterior as the probability densities would be maximum only at the mode of the parameters in the distribution
OptimRes <- optim(prior,LogPostLogistic,gr=NULL,y,X,mu,Sigma,method=c("BFGS"),control=list(fnscale=-1),hessian=TRUE)

# Getting the hessian matrix from the results returned from the optim function
beta_mode_hessian = OptimRes$hessian
# Getting the mode β values as the parameters returned from the optim function
beta_mode = OptimRes$par
# Getting the J(beta_mode) as the negative of the hessian returned 
J = -beta_mode_hessian
# Calculating inverse of J
J_inverse = solve(-beta_mode_hessian)

# Presenting values of beta_mode and inverse J(beta_mode)
print(beta_mode)
print(J_inverse)

# Verifying if the estimated results are reasonable using glmModel
glmModel <- glm(Work ~ 0 + ., data = data, family = binomial)
# The values are found to be quite similar
glmModel

# Sampling 1000 draws from the posterior of the logistic regression 
N = 1000
posterior_beta_draws = rmvnorm(N,beta_mode, J_inverse)

# Plotting the density of the regression coefficients of the variable NSmallChild (6th covariate)
plot(main=col_names[-1][6],density(posterior_beta_draws[,6]))
print(summary(posterior_beta_draws[,6]))
# Getting sampled values for the beta paramters for NSmallChild(6th covariate)
small_kids_draw = posterior_beta_draws[,6]
# Finding the 95% equal tail credible interval from the sampled values for NSmallChild and adding it to the plot
q = quantile(small_kids_draw, probs=c(0.025,0.975))
print(q)
q = unname(q)
q_lo = q[1]
q_hi = q[2]
abline(v=q_lo, col="red")
abline(v=q_hi, col="red")
legend("topleft", inset=.02, c("95% equal tailed credible intervel"), fill=c("red"))



# Lab 2.2.b)

# Creating a data point for a women of age 43, with 2 small children, 12 years of education, 8 years experience and with husband income 20
X_test = c(1,20,12, 8, 43, 2, 0)

# Getting predictions that women is working based for the test case created above based on the sampled beta parameters from the posterior of logistic regression
res = X_test %*% t(posterior_beta_draws) 
yhat_test = exp(res)/(1+exp(res))

# Plotting the distribution of the predicted values of the test case
plot(main = "Density plot of probablity of the women working", density(yhat_test))
hist(yhat_test, main = "Histogram of Prediction of probablity of the women working")

# Lab 2.2.c)
# Finding the distribution for the number of working women out of the 11 of the same type specified above 
# Since the binomial distribution is the sum of bernolli trial
# Hence random draws from the binomial distribution with number of trails as 11 (considering 11 women of same specification)
# and probability of success (working) found for a single women can be get the samples of number of women working out of 11
num_of_working_women <- c()
# looping through all the beta values sampled for posterior
for (ind in 1:dim(posterior_beta_draws)[1]){
  # Considering p as the probability of success for 1 women found using the beta parameter
  p  = yhat_test[ind]
  # Getting number of women working as random draws from the binomial distribution with n = 11 and p = probability of success for one women found using posterior betas
  num_of_working_women = c(num_of_working_women, rbinom(1, 11, p))
}
# Plotting posterior predictive distribution of number of women working out of 11
hist(num_of_working_women, freq = FALSE, probability = TRUE, main = "Histogram of probablity of number of out of 11 women working")
plot(density(num_of_working_women), main = "Density plot of probablity of number of out of 11 women working")

