library(scales)
data = readRDS("Precipitation.rds")

# 3.1.a)
y = log(data)
n = length(y)
y_mean = mean(y)

# initializing hyper parameters
mu_0 = mean(y)*1.2
t_sq_0 = var(y)*.8
sigma_sq_0 = 1.8*2
v0 = 100*.1


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

if_mu = 1+2*sum(acf_mu$acf[-1])
if_sigma = 1+2*sum(acf_sigma$acf[-1])

hist(Gibbs_posterior_mu)
hist(Gibbs_posterior_sigma)
plot(Gibbs_posterior_mu,type="l")
plot(Gibbs_posterior_sigma,type="l")

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

# 1.b

hist(data)
draws_posterior = c()
for(i in 1:n){
  draws_posterior = c(draws_posterior, rnorm(1,Gibbs_posterior_mu[i], Gibbs_posterior_sigma[i]))
}
hist(exp(draws_posterior), breaks = 200, xlim=range(0,80))

