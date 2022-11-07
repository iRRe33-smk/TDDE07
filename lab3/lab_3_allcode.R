library(scales)
data = readRDS("Precipitation.rds")

# 3.1.a)
y = log(data)
n = length(y)
y_mean = mean(y)

# initializing hyper parameters
mu_0 = mean(y)
t_sq_0 = var(y)
sigma_sq_0 = 1.8
v0 = 100


# Plotting prior to check if it is similar to the data given
prior_sigmasq = c()
N_Draws = 2000
prior_mu = rnorm(N_Draws,mu_0,t_sq_0)
for(i in 1:N_Draws){
  X = rchisq(1,v0)
  prior_sigmasq = c(prior_sigmasq, v0*sigma_sq_0 / X)
}

# Drawing form the simulated prior mean and sigma^2
prior_datapoints = c()
for(i in 1:N_Draws){
  prior_datapoints = c(prior_datapoints, rnorm(1, prior_mu[i],prior_sigmasq[i]))
}

# Plotting histogram of data points
hist(prior_datapoints)
mean(prior_datapoints)
var(prior_datapoints)

# Initial Setting 
N_Draws = n
mu = rnorm(1,mu_0,t_sq_0)
X = rchisq(1,v0)
sigmasq = v0*sigma_sq_0 / X

Gibbs_posterior_mu = c()
Gibbs_posterior_sigma = c()

# Posterior form Gibbs Sampling
for( i in 1:N_Draws){
  
  w = (n/sigmasq) / ( (n/sigmasq) + (1 / t_sq_0))
  mu_n = w*y_mean + (1-w)*mu_0
  t_sq_n = 1/((n/sigmasq) + (1/t_sq_0))
  
  mu = rnorm(1,mu_n,t_sq_n)
  Gibbs_posterior_mu = c(Gibbs_posterior_mu, mu)
  
  v_n = v0 + n
  X = rchisq(1,v_n)
  thing = (v0*sigma_sq_0 + sum((mu - y)^2))/(n + v0)
  sigmasq = v_n*thing/X
  Gibbs_posterior_sigma = c(Gibbs_posterior_sigma, sigmasq)
}

acf_mu = acf(Gibbs_posterior_mu)
acf_sigma = acf(Gibbs_posterior_sigma)

# Calculating Inefficiency Factors for mu and sigma^2
if_mu = 1+2*sum(acf_mu$acf[-1])
if_sigma = 1+2*sum(acf_sigma$acf[-1])

hist(Gibbs_posterior_mu)
hist(Gibbs_posterior_sigma)
plot(Gibbs_posterior_mu,type="l", xlab = "iteration", ylab = "mu", main="Sampled Markov chain for Mu")
plot(Gibbs_posterior_sigma,type="l", xlab = "iteration", ylab = "sigma", main="Sampled Markov chain for sigma")

# Direct Draws Done for verification
mu_0 = mean(y)
k0 = 1
sigma_sq_0 = 1.8
v0 = 100

# Draws from Piror
nDraws = 2000

X = rchisq(nDraws,v0)
prior_direct_sigmasq = v0*sigma_sq_0 / X

prior_direct_mu = c()
for(i in 1:nDraws){
  prior_direct_mu = c(prior_direct_mu, rnorm(1,mu_0,prior_direct_sigmasq[i]/k0))
}
hist(mu)
hist(prior_direct_sigmasq)

# Posterior from Direct Sampling
nDraws = 500
mu_n = ((k0/(k0+n))*mu_0) + ((n/(k0+n))*y_mean)
kn = k0 + n
vn = v0 + n 
s_sq = sum((y-y_mean)^2)/(n-1)
sigmasq_n = (v0*sigma_sq_0 + (n-1)*s_sq + (((k0*n)/(k0+n))*(y_mean - mu_0)^2))/vn

X = rchisq(nDraws,vn)
direct_sigmasq = vn*sigmasq_n / X

direct_mu = c()
for(i in 1:nDraws){
  direct_mu = c(direct_mu, rnorm(1,mu_n,direct_sigmasq[i]/kn))
}

hist(direct_mu)
hist(direct_sigmasq)
plot(direct_mu,type="l")
plot(direct_sigmasq,type="l")

# lab 3.1.b

hist(data)
draws_posterior = c()
for(i in 1:n){
  draws_posterior = c(draws_posterior, rnorm(1,Gibbs_posterior_mu[i], Gibbs_posterior_sigma[i]))
}
hist(exp(draws_posterior), breaks = 100, xlim=range(0,80))

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
  plot(metropolis_posterior_betas[,i], type = "l", ylim = c(0,2*model$coefficients[i]))
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

mu = 13
sigma_sq = 3
t_max = 300
fi = 1
# lab 3.3.1

# Function to simulate Auto Regressive Function
AR_Process <- function(mu, FI, sigma_sq)
{
  path = matrix(data = NA, nrow = length(FI), ncol = t_max)
  i = 1
  for (fi in FI){
    path[i,1] = mu
    for (t in 2:t_max){
      path[i,t] = mu + fi*(path[i,t-1]-mu) + rnorm(1,0,sigma_sq)
    } 
    i = i+1
  }
  return(path)
}

# Simulating Draws from the AR1 process with phi value -0.8, 0, 0.8
FI = seq(-.8,.8,.8)
path = AR_Process(mu, FI, sigma_sq)
par(mfrow=c(3,2))
for (fi in 1:length(FI)){
  plot(path[fi,],main=paste("phi = ",FI[fi]),type="line",col="blue")  
  acf(path[fi,],main=paste("ACF phi = ",FI[fi]))
}

# lab 3.3.2 

# Creatind Data to pass to stan model
fi_1 = .2
fi_2 = .95

FI = c(fi_1, fi_2)
path = AR_Process(mu, FI, sigma_sq)
par(mfrow=c(2,2))
for (fi in 1:length(FI)){
  plot(path[fi,],main=paste("phi = ",FI[fi]),type="line",col="blue")  
  acf(path[fi,],main=paste("ACF phi = ",FI[fi]))
}

library(rstan)

# Creating Stan Model to simulate the two AR(1)-processes

StanModel = '
data {
  matrix[300,2] M; // synthetic data
  
}
parameters {
  real mu;
  real<lower=0> sigma2;
  real phi_1;
  real phi_2;
}
model {
  mu ~ uniform(-100,200); // Normal with mean 0, st.dev. 100
  sigma2 ~ uniform(0,100); // Scaled-inv-chi2 with nu 1,sigma 2
  phi_1 ~ uniform(-10,10);
  phi_2 ~ uniform(-10,10);
  
  for(i in 2:300){
    M[i,1] ~ normal(mu + phi_1*(M[i-1,1]-mu),sigma2);
    M[i,2] ~ normal(mu + phi_2*(M[i-1,2]-mu),sigma2);   
  }
}'

data <- list(M=t(path))
warmup <- 10
niter <- 2000
fit <- stan(model_code=StanModel,data=data, warmup=warmup,iter=niter,chains=4)

# Print the fitted model
print(fit,digits_summary=3)

# Extract posterior samples
postDraws <- extract(fit)

# Do traceplots of the first chain
par(mfrow = c(1,1))

# Do automatic traceplots of all chains
traceplot(fit)
# Plotting joint posterior for mu and phi1
plot(postDraws$mu[1:(niter-warmup)], postDraws$phi_1[1:(niter-warmup)], xlab = "mu", ylab = "phi_1", main = "Joint Posterior of mu and phi 1")
# Plotting joint posterior for mu and phi2
plot(postDraws$mu[1:(niter-warmup)], postDraws$phi_2[1:(niter-warmup)], xlab = "mu", ylab = "phi_1", main = "Joint Posterior of mu and phi 2")
# Bivariate posterior plots
pairs(fit)






