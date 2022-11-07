#-------------------1.a-------------------
#prior
a = 16
b = 24

n=100
choice_a = 38
choice_b = 27
choice_c = 35


#posterior 
post_a = a + choice_a
post_b = b + choice_b + choice_c

val = .4
prob_less_than = pbeta(val,shape1 = post_a, shape2 = post_b,lower.tail = FALSE)
print(prob_less_than)


N = 10000
random_draws = rbeta(N,shape1 = post_a, shape2 = post_b)

hist(1-random_draws,probability=TRUE, main = "probability distribution of 1-Theta_a", sub = paste("prob 1-theta_a < ",val, " = ",prob_less_than,sep=""))

abline(v=1-val, col="red")
#-----------------------------1.b------------------

rat = (1-random_draws)/random_draws

hist(rat)

q=quantile(x=rat,probs=c(0.025,0.975))
abline(v=q, col="red")

1/(1+q)


#---------------------1.c----------------

#marginal likelihood
#posterior over prior. When conjugate prior

marg= beta(a+choice_a,a+b+n) / beta(a,b)


#-----------------------1.d----------------

alpha = 60*c(1,1,1)/3 
y= c(38,27,35)

n=1E5
rDirichlet = function(n,alpha,y){
  z = matrix(nrow=n,ncol=length(alpha))
  for (j in 1:n){
    x = matrix(nrow=1,ncol=length(alpha))
    for (i in 1:length(alpha)){
      x[1,i] = rgamma(n=1,shape =a+y[i])
    } 
    
    z[j,] = x / sum(x)
  }  
  return(z)
  
  
  
}
z = rDirichlet(n,alpha,y)

posterior_prob = mean(z[,1] > z[,3])
print(posterior_prob)


#--------------------------------3.a--------------------
library(geoR)
library(mvtnorm)
head(X)

#priors
#sigma_sq_0 = rinvchisq(n = 1,df = 1,scale = 2^2)
ncovs = dim(X)[2]
mu_0 = rep(0,ncovs)
omega_0 = 1/25*diag(ncovs)
nu_0 = 1

sigma_sq_0 = 2^2



#BayesLinReg(y,X,mu_0,omega_0,nu_0,sigma_sq_0,10000)

draws = BayesLinReg(y = y,X = X,mu_0 = mu_0,Omega_0 = omega_0,v_0 = nu_0,sigma2_0 = sigma_sq_0,nIter = 10000)

intervals = matrix(nrow=ncovs,ncol=3)
colnames(intervals) = c("2.5%","mean","97.5%")
rownames(intervals) = colnames(X)
for (i in 1:ncovs){
  draw_i = draws$betaSample[,i]
  intervals[i,c(1,3)] = quantile(draw_i,probs=c(0.025,0.975))
  intervals[i,2] = mean(draw_i)
}
intervals

#interval for b1 /x1 is positive indicating that the variable x1 is significant to our data.

#----------------------3.b---------------
posterior_median = median(draws$sigma2Sample)
posterior_median

#----------------------3.c----------------
intervals[c("x1","x1*x3","x1*x4"),]

beta1 = draws$betaSample[2,]

beta5 = draws$betaSample[6,]

beta6 = draws$betaSample[7,]


b1_B = beta1 + beta5
b1_C = beta1 + beta6

diff = b1_B - b1_C

plot(density(diff))

