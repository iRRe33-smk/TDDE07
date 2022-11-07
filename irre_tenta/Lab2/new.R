#-----------------------------------------------------------
#-------------------------------LAB2------------------------
#-----------------------------------------------------------


library(mvtnorm)
#--------------------------------2.1a----------------------
dat = read.table("TempLambohov.txt",header = T )
summary(dat)

plot(dat)

y = dat[,2]
X = data.frame(dat[,1])
X$time_sq = dat[,1]^2

colnames(X) = c("time", "time_sq")
X$const = 1

X = X[,c(3,1,2)]
head(X)

#hyperparams (for beta)
mu_0 = c(-10, 100,-100)
omega_0 = 0.1 * diag(3)
nu_0 = 3
sigma_sq_0 = 2


#draw sigma from inv_chgi_sq
numdraws = 50
Xdraws = rchisq(numdraws,nu_0)

#as per slide in L3 but with other params
sigma_sq = nu_0*sigma_sq_0/Xdraws
#sigma_sq

#draw beta from N
betas = matrix(nrow=numdraws,ncol=3)
for (i in 1:numdraws){
  tmp = rmvnorm(1,mu_0,sigma_sq[i]*solve(omega_0)) # L5
  #print(tmp)
  betas[i,] = tmp 
}

X_sample = data.frame(const = 1, time = seq(0,1,1/364),time_sq = seq(0,1,1/364)^2)
plot(dat,ylim=c(-100,100), main="Measured data and simulated draws.")
draws_predict =betas%*%t(X_sample) 
for (i in 1:numdraws){
  points(X_sample$time,draws_predict[i,],col="#F7000010")
}

#-----------------------------------2.1b.i---------------------------------
"
create all the '_n' varialbes
draw some number of sigmas.
use the new sigmas to draw betas
show marginalö posterior for betas
X,y, omega0,  mu_0, v0, n, sigma0_sq
"
posterior <- function(numdraws, X, y, omega_0, mu_0, nu_0,sigma_sq_0 ){
  n = length(y)
  XtX = t(X)%*%X
  mu_n = 0
  #mu_n
  beta_hat = solve(t(X)%*%X)%*%t(X)%*%y # (X.t * X)^-1 * X.t * y
  
  mu_n = solve(XtX + omega_0) %*% (XtX%*%beta_hat + omega_0%*%mu_0)
  
  #omega_n
  omega_n = XtX + omega_0
  
  #nu_n
  nu_n = nu_0 + n
  #vs_n
  vs_n = nu_0*sigma_sq_0 + (t(y)%*%y + t(mu_0)%*%omega_0%*%mu_0 - t(mu_n)%*%omega_n%*%mu_n)
  
  
  #draw sigmas
  X_draws = rchisq(numdraws, nu_n)
  sigma_sq = vs_n / X_draws
  
  #print(sigma_sq)
  
  parameters = matrix(nrow = numdraws,ncol=length(mu_0)+1)
  om_inv = solve(omega_n)
  
  for (i in 1:numdraws){
    parameters[i,1:length(mu_0)] = rmvnorm(1,mean= mu_n,sigma = sigma_sq[i]*om_inv)
  }
  parameters[,length(mu_0)+1] = sigma_sq
  
  return(parameters)
}

parameters = posterior(1000, as.matrix(X), as.matrix(y), omega_0, mu_0, nu_0, sigma_sq_0 )

print(head(parameters))



hist(parameters[,1],freq = F,main="beta0- constant")#,xlim=c(-100,300))
hist(parameters[,2],freq = F, main ="beta1- time")#,xlim=c(-100,300))
hist(parameters[,3],freq = F, main = "beta2- time^2")#,xlim=c(-100,300))
hist(parameters[,4],freq = F, main = "sigma^2- variance")

#------------------------------------2.1b.ii--------------------------


plot(dat)

betas = parameters[,1:3]

n = dim(X_sample)[1]
preds = betas %*%t(X_sample)

median_preds = matrix(ncol=n)
quants = matrix(nrow=2,ncol=n)
for (i in  1:n){
  median_preds[i] = median(preds[,i])
  quants[,i] = quantile(preds[,i],c(.025,0.975))

}

points(X_sample$time,median_preds, col = "#ff000023")
points(X_sample$time,quants[1,], col ="#0000ff23")
points(X_sample$time,quants[2,], col ="#0000ff23")
#points(X_sample$time,median_preds, col = alpha("red",.8))

#--------------------------------2.1c---------------------


beta_derivs = matrix(nrow = dim(betas)[1], ncol = 2)
zeros = matrix(rep(0,n),nrow=n)
beta_derivs[,1] = betas[,2]
beta_derivs[,2] = 2*betas[,3]
highest_temp = -beta_derivs[,1]/beta_derivs[,2]
head(highest_temp)
hist(highest_temp,freq=F, main="Posterior distribtuiton for x~, time of year with highest expected temp")

#-----------------------------2.2a------------------------------
library(mvtnorm)

#data
dat = read.table("WomenAtWork.dat",header=T)
head(dat)

X = as.matrix(dat[,-1])
y = as.matrix(dat[,1])



#priors

#mean (for betas)
num_covs = dim(X)[2]
betas_prior = as.matrix(rep(0,num_covs)) # ~N(0,tau^2 * I)

tau = 5
Sigma_prior = tau^2 * diag(num_covs)


LogPostLogistic <- function(betas,y,X,mu,Sigma){
  linPred <- X%*%betas;
  logLik <- sum( linPred*y - log(1 + exp(linPred)) );
  #if (abs(logLik) == Inf) logLik = -20000; # Likelihood is not finite, stear the optimizer away from here!
  logPrior <- dmvnorm(t(betas), t(mu), Sigma, log=TRUE);
  
  return(logLik + logPrior)
}

LogPostLogistic(betas_prior,y,X,betas_prior,Sigma_prior)

optres = optim(betas_prior,LogPostLogistic, gr = NULL,y,X,betas_prior, Sigma_prior,method=c("BFGS"),control=list(fnscale=-1),hessian=TRUE)
optres

J = optres$hessian
beta_mode = optres$par



beta_posterior_draws = rmvnorm(1E4,mean = beta_mode, sigma = solve(-J))

small_child_draws = beta_posterior_draws[,6]

hist(small_child_draws, freq=F,sub="red lines indicate boundries of 95% equal tail")
abline(v=quantile(small_child_draws,probs=c(0.025,0.975)), col="red")
#the interval is fully negative, indicating we number of small children is significant at least to that level of certainty

glm(formula = "Work ~ 0 + .",data = dat, family = "binomial")
#---------------------------------2.2b-----------------------

posterior_draw_predictive <- function(X, beta_mode, Sigma, numdraws = 1){
  
  betas = rmvnorm(numdraws,beta_mode,sigma=Sigma)
  
  linpred  = matrix(nrow=numdraws,ncol=1)
  
  for (i in 1:numdraws){
    linpred[i] = X %*% betas[i,]
  }  
  
  tmp =exp(linpred)
  prob = tmp / ( 1 + tmp)
  
  return(prob)
  
  
}

N = 1E4
X_test = matrix(data = c(1, 20, 12, 8, 43, 0,2),ncol = 7, nrow=1 )

probs = posterior_draw_predictive(X_test,beta_mode,solve(-J), numdraws = N)

hist(probs, freq=F, main= "Posterior draws for probabilty that the woman is working")
lines(density(probs))



#---------------------------2.2c-----------------------------

posterior_draw_predictive_11 <- function(X, beta_mode, Sigma, numdraws = 1){
  
  betas = rmvnorm(numdraws,beta_mode,sigma=Sigma)
  
  num_women  = matrix(nrow=numdraws,ncol=1)
  
  for (i in 1:numdraws){
    linpred = X %*% betas[i,]
    tmp =exp(linpred)
    prob = tmp / ( 1 + tmp)
    num = rbinom(1,11,prob)
    num_women[i] = num
    }  

  return(num_women)
  
  
}

N = 1E4
X_test = matrix(data = c(1, 20, 12, 8, 43, 0,2),ncol = 7, nrow=1 )

num_women = posterior_draw_predictive_11(X_test,beta_mode,solve(-J), numdraws = N)

hist(num_women, freq = F)


