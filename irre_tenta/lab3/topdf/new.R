#-----------------------------------------------------------
#-------------------------------LAB3------------------------
#-----------------------------------------------------------

library(scales)
data = readRDS("Precipitation.rds")

#---------------------------------3.1--------------------

y = log(data)
n = length(y)
y_mean = mean(y)

# initializing hyper parameters 
#these are not given and no method for choosing provided
#have chosen reasonable values and checked with draws from posterior
mu_0 = mean(y)*1.2 
t_sq_0 = var(y)*.8
sigma_sq_0 = .5
v0 = 100*.1


# Plotting prior to check if it is at all similar to the data given
#simulated draws from inverss chi^2
prior_sigmasq = c()
N_Draws = 2000
prior_mu = rnorm(N_Draws,mu_0,t_sq_0)
for(i in 1:N_Draws){
  X = rchisq(1,v0)
  prior_sigmasq = c(prior_sigmasq, v0*sigma_sq_0 / X)
}

# Drawing form the simulated prior mean and sigma^2
draws_from_prior = c()
for(i in 1:N_Draws){
  draws_from_prior = c(draws_from_prior, rnorm(1, prior_mu[i],prior_sigmasq[i]))
}



# Plotting histogram of data points
hist(draws_from_prior)
mean(draws_from_prior)
var(draws_from_prior)

#plotting given data
hist(y)
mean(y)
var(y)



#------------------gibbs sampling------------------

N_Draws = n
# initializing hyper parameters 
mu_0 = mean(y)*1.2 
t_sq_0 = var(y)*.8
sigma_sq_0 = .5
v0 = 100*.1
v_n = v0 + n
# Initial Setting 

mu = rnorm(1,mu_0,t_sq_0) # draw using our hyperparams

#draw from inv_chi^2
X = rchisq(1,v0) 
sigmasq = v0*sigma_sq_0 / X

Gibbs_posterior_mu = c()
Gibbs_posterior_sigma = c()


#mu is mean of the normal distribution for y
#mu_n is the mean of the normal distribution fro mu
#t_sq_n is the variance of the normal distribution for mu

# v0, v_n and sigma_sq_0 are used in drawing for sigma_sq

for( i in 1:N_Draws){
  
  w = (n/sigmasq) / ( (n/sigmasq) + (1 / t_sq_0))
  mu_n = w*y_mean + (1-w)*mu_0
  t_sq_n = 1/((n/sigmasq) + (1/t_sq_0))
  
  mu = rnorm(1,mu_n,t_sq_n)
  Gibbs_posterior_mu = c(Gibbs_posterior_mu, mu)
  
  X = rchisq(1,v_n)
  thing = (v0*sigma_sq_0 + sum((mu - y)^2))/(n + v0)
  sigmasq = v_n*thing/X
  Gibbs_posterior_sigma = c(Gibbs_posterior_sigma, sigmasq)
}
#autocorrelations of sampled params
acf_mu = acf(Gibbs_posterior_mu)
acf_sigma = acf(Gibbs_posterior_sigma)

#inefficiencty factors
if_mu = 1+2*sum(acf_mu$acf[-1])
if_sigma = 1+2*sum(acf_sigma$acf[-1])

hist(Gibbs_posterior_mu)
hist(Gibbs_posterior_sigma)
plot(Gibbs_posterior_mu,type="l")
plot(Gibbs_posterior_sigma,type="l")

#--------------verification, through direct draws, not part of question -----------------------
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

#-----------------------------------------------3.1.b--------------------------

hist(data)
draws_posterior = c()
for(i in 1:n){
  draws_posterior = c(draws_posterior, rnorm(1,Gibbs_posterior_mu[i], Gibbs_posterior_sigma[i]))
}
hist(exp(draws_posterior), breaks = 200, xlim=range(0,80))

#-----------------------------------------------------------------------------
#----------------------------------------------3.2----------------------------
#-----------------------------------------------------------------------------

#----------------------------------------------3.2.a--------------------------


library(mvtnorm)
data <- read.table("eBayNumberOfBidderData.dat", sep = "" , header = T , na.strings ="", stringsAsFactors= F)
data = as.data.frame(data)

col_names = colnames(data)

#poisson distributed, with exponentail in expression
model <- glm("nBids ~ . - const" , data, family = poisson(link = "log"))
model
summary(model)

#----------------------------------3.2.b--------------------------------
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

#-----------------------------------3.2c---------------------------------
num_random_walk = 1E5
MetropolisRandomWalk <- function(beta_initial, c, J_inverse, log_density_function, mu, Sigma, Y, X){
  
  metropolis_betas = matrix(data = NA, nrow = num_random_walk, ncol = dim(X)[2])
  prev_beta = beta_initial
  
  for(i in 1:num_random_walk){
    
    sample = rmvnorm(1, mean = prev_beta, sigma = c * J_inverse)
    
    log_density_sample = log_density_function(sample, Y, X, mu, Sigma)
    log_density_prev_sample = log_density_function(prev_beta, Y, X, mu, Sigma) 
    
    alpha = min(1, exp((log_density_sample - log_density_prev_sample)))
    #comparing the sum if log-prior and log likelihood of the two samples
    #update betas if new one is sufficiently 'better', based on random draw
    
    rand = runif(n=1, min = 0, max = 1)
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

metropolis_posterior_betas = MetropolisRandomWalk(posterior_beta_draw, .6, J_inverse, LogPostPoisson, mu, Sigma, Y, X)
colnames(metropolis_posterior_betas) = col_names[-1]



par(mfrow=c(3,3))
for(i in 1:dim(X)[2]){
  plot(metropolis_posterior_betas[,i], type = "l",ylab = col_names[1+i])# ylim = c(0,2*model$coefficients[i]))
  abline(h=model$coefficients[i], col="red")
}

#-------------------------------------------3.2.d------------------------------------

par(mfrow=c(1,1))
sample = matrix(c(1, 1, 0, 1, 0, 1,0,1.2,0.8))
intensities = exp(metropolis_posterior_betas %*% sample)

poisson_draws = matrix(nrow=num_random_walk,ncol=1)
for (i in 1:num_random_walk){
  poisson_draws[i,1] = rpois(1,intensities[i])
}

#print(exp(model$coefficients%*%sample))
hist(poisson_draws, main = "Predictive distribution", sub =paste("glm predicts: ",exp(model$coefficients%*%sample)))

prob_no_bids = mean(poisson_draws==0)
prob_no_bids


#----------------------------------------------------------------------------------------
#------------------------------------------3.3-------------------------------------------
#----------------------------------------------------------------------------------------

mu = 13
sigma_sq = 3
t_max = 300
fi = 1
# ------------------ 3.3.a-------------------------------------------------------------
par(mfrow=c(3,2))
for (fi in seq(-.8,.8,.8)){
  path = matrix(nrow = t_max)
  path[1] = mu
  for (t in 2:t_max){
    path[t] = mu + fi*(path[t-1]-mu) + rnorm(1,0,sigma_sq)
  } 
  plot(path,main=fi,type="line",col="blue")  
  acf(path,main=paste("ACF ",fi))
}


#----------------3.3.2 create data---------------------------------------------------------

fi_1 = .2
mu_1 = 15
sigma_sq_1 = 4


fi_2 = .95
mu_2 = 2
sigma_sq_2 = 8


par(mfrow=c(1,1))

path = matrix(nrow = t_max,ncol=2)
path[1,] = c(mu_1,mu_2)
for (t in 2:t_max){
  path[t,1] = mu_1 + fi_1*(path[t-1,1]-mu_1) + rnorm(1,0,sqrt(sigma_sq_1))
  path[t,2] = mu_2 + fi_2*(path[t-1,2]-mu_2) + rnorm(1,0,sqrt(sigma_sq_2))
} 
plot(path[,2],type="line",col="blue")  
lines(path[,1],col="red")#acf(path,main=paste("ACF ",fi))



write.table(path[,1],file="path1.txt",header=T)
write.table(path[,2],file="path2.txt",header=T)
paths = path

library(rstan)
#y=c(4,5,6,4,0,2,5,3,8,6,10,8)
#N=length(y)

StanModel = '
data {
  vector[300] M; // synthetic data
  //matrix[300,2] M
  
  
}
parameters {
  real mu;
  real<lower=0> sigma2;
  real phi;
}
model {
  mu ~ uniform(-100,100); // Normal with mean 0, st.dev. 100
  sigma2 ~ uniform(0,100); // Scaled-inv-chi2 with nu 1,sigma 2
  phi ~ uniform(-10,10);
  
  
  for(i in 2:300){
                  
    M[i] ~ normal(mu + phi*(M[i-1]-mu),sqrt(sigma2));
  }
}'

data1 <- list(M=paths[,1])
data2 <- list(M=paths[,2])

################################################################################
#fitting for first series
niter <- 5000
warmup <- as.integer(niter*.15)

fit <- stan(model_code=StanModel,data=data1, warmup=warmup,iter=niter,chains=4)
# Print the fitted model
print(fit,digits_summary=3)
# Extract posterior samples
postDraws <- extract(fit)
# Do traceplots of the first chain
par(mfrow = c(1,1))
plot(postDraws$mu[1:(niter-warmup)],type="l",ylab="mu",main="Traceplot")
# Do automatic traceplots of all chains
traceplot(fit)
# Bivariate posterior plots
pairs(fit)

################################################################################
#fitting for second series
niter <- 5000
warmup <- as.integer(niter*.15)
fit <- stan(model_code=StanModel,data=data2, warmup=warmup,iter=niter,chains=4)
# Print the fitted model
print(fit,digits_summary=3)
# Extract posterior samples
postDraws <- extract(fit)
# Do traceplots of the first chain
par(mfrow = c(1,1))
plot(postDraws$mu[1:(niter-warmup)],type="l",ylab="mu",main="Traceplot")
# Do automatic traceplots of all chains
traceplot(fit)
# Bivariate posterior plots
pairs(fit)
































