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

#--------------verification, not part of question -----------------------
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

