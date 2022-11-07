library(mvtnorm)
library(scales)

# Lab 2.1.a)

# Initializing hyper parameters of the linear regression model
mu0 = c(-10, 100, -100)
omega0 = 0.12 * diag(3)
v0 = 3
sigma0_sq = 2

# Number of elements to sample for σ^2 and β
N <- 1E2

# Reading data from the given txt file
data <- read.table("TempLambohov.txt", sep = "" , header = T , na.strings ="", stringsAsFactors= F)

# Creating a matrix of 1, time, time^2 from the data read
X = data.frame(const = 1, t = data$time, t2 = data$time^2)

# creating a sampling of 201 values between 0 and 1 for the time variable
time_sample = data.frame(const=1,t = seq(0,1,0.005),t2 =seq(0,1,0.005)^2)

# function to calculate temperate based on the linear regression formula β0 + β1 * time + β2 *time^2
pred_temp <- function(betas){
  pred_temp = as.matrix(time_sample)%*%t(betas)
}

# Sampling for prior σ^2 from Inv−χ2(v0,σ0^2) distribution
X_draws = rchisq(N,v0)                # Taking N samples from the χ2(v0)
sigma_sq = (v0)*sigma0_sq / X_draws   # Simulating draws from Inv−χ2(v0,σ0^2) distribution

# Plotting data points
plot(data,ylim=c(-25, 60))


betas_prior = matrix(nrow=N,ncol=3 )
x_plot = seq(0,1,0.005) 
for (i in 1:N){
  # Simulating draws for prior β from N(μ0, σ^*2Ω0^−1)
  beta = rmvnorm(1, mean = mu0, sigma = sigma_sq[i]*solve(omega0))
  betas_prior[i,]=beta
  # Plotting temperatures predicted from the prior draws
  points(x_plot,t(pred_temp(beta)),col=alpha("red",0.07))
}

# Lab 2.1.b)

# Function to simulate for σ^2 and β from the joint posterior distribution
posterior <- function(X,y, omega0,  mu_0, v0, n, sigma0_sq ){
  posterior_parameters = matrix(nrow=N,ncol=4)
  
  # calculating beta_hat = ((X′X)^−1)X′y 
  beta_hat = solve(t(X)%*%X)%*%t(X)%*%y
  # calculating μn = ((X′X+Ω0)^−1)(X′X*beta_hat + Ω0*μ0)
  mu_n = solve(t(X)%*%X + omega0) %*% (t(X)%*%X%*%beta_hat + omega0%*%mu_0)
  # calculating Ωn = X′X+Ω0
  omega_n = t(X)%*%X + omega0
  # calculating νn = ν0 + n
  v_n = v0 + n
  # calculating v_s = vn*σn^2 = v0 * σ0^2 + (y′y + μ0′* Ω0 * μ0 − μn′* Ωn * μn)
  v_s = v0*sigma0_sq + (t(y)%*%y + t(mu_0)%*%omega0%*%mu_0 - t(mu_n)%*%omega_n%*%mu_n)
  
  # Sampling for posterior σ^2 from Inv−χ2(vn,σn^2) distribution
  X_draws = rchisq(N,v_n)    # Taking N samples from the χ2(vn)
  sigma_sq = v_s / X_draws   # Simulating draws from Inv−χ2(vn,σn^2) distribution
  
  posterior_parameters[,4] = sigma_sq
  temps = matrix(nrow=N,ncol=201)
  for (i in 1:N){
    # Simulating draws for posterior β from N(μn, σ^*2Ωn^−1)
    beta = rmvnorm(1, mean = mu_n, sigma = sigma_sq[i]*solve(omega_n))
    posterior_parameters[i,1:3]=beta
    
    # Predicting temperature values from the posterior β
    temps[i,] = t(pred_temp(beta))
    
    # plotting predcited temperature values
    points(x_plot,t(pred_temp(beta)),col=alpha("blue",0.01))
    
  }
  ret = list(posterior_parameters=posterior_parameters, temps=temps)
  return(ret)
}
legend("topleft", inset=.02, c("Tempature Data", "Calulated Prior Temperature", "Caulated posterior Temperature"), fill=c("black", "red", "blue"))
ret = posterior(as.matrix(X),as.matrix(data$temp),omega0, mu0,v0,dim(data)[1],sigma0_sq)
betas = ret$posterior_parameters[,1:3]
sigma_sq = ret$posterior_parameters[,4]
temps = ret$temps

# # 2.1.b.i)
# # Plotting histogram for each of the marginal posterior parameters
hist(betas[,1], freq = FALSE, main = "Histogram of Beta 0")
hist(betas[,2], freq = FALSE, main = "Histogram of Beta 1")
hist(betas[,3], freq = FALSE, main = "Histogram of Beta 2")
hist(sigma_sq, freq = FALSE, main = "Histogram of Sigma^2")



# 2.1.b.ii)

# Scatter plot of temperature value
plot(data,ylim=c(-25, 60))

post_medians = c()
quants = matrix(nrow=dim(temps)[2],ncol = 2)
for (i in 1:dim(temps)[2]){
  # Finding the median of the predicted temperatures
  med = median(temps[,i])
  # finding the 95% equal tail posterior probablity interval of predicted terperature
  q = quantile(temps[,i],probs=c(.025,0.975))
  q = unname(q)
  
  quants[i,] = c(q[1],q[2])
  
  post_medians = c(post_medians,med)
}
lines(x_plot,post_medians, col = "green")
lines(x_plot,quants[,1], col = "blue")
lines(x_plot,quants[,2], col = "red")

# 2.1.c)

# The time at which the temperature is maximum can be found by solveing the function f'(time) = 0 and solving for time
# This points can be considered as the maximas as f''(time) = 2 * β2 < 0 for all values in β2
# hence time at which temperature is max = -β1/(2*β2)


max_temp_time = -betas[,2]/(2*betas[,3])
abline(v=max_temp_time, col=alpha("violet",0.1))
legend("topleft", inset=.02, c("Tempature Data", "Upper limit 95% equal tail credible intervel", 
                               "lower limit 95% equal tail credible intervel", "Median Posterior Temperature", "Max Temperature time period"), 
       fill=c("black", "red", "blue", "green", "violet"))




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
