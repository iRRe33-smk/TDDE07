#-----------------------------------------------------------
#-------------------------------LAB1------------------------
#-----------------------------------------------------------


#---------------------------------1.1a----------------------
#priors
prior_alpha = 5
prior_beta = 5

#observed data
trials = 50
sucesses = 13
fails = trials - sucesses

#posteriors  for beta distr
posterior_alpha = (sucesses+prior_alpha)
posterior_beta = (fails + prior_beta)


N = 5000

drawn_mean = c()
drawn_sd = c()
for ( n in seq(1,N)){
  draw = rbeta(n,prior_alpha+sucesses,prior_beta+fails)
  
  drawn_mean = c(drawn_mean , mean(draw))
  drawn_sd = c(drawn_sd , sd(draw))   
}

#true vals
#mean / expected value of Beta-distr
true_mean = (sucesses+prior_alpha) / (trials + prior_alpha + prior_beta) 

#sd of Beta distr
true_sd = sqrt((posterior_alpha*posterior_beta) / ((posterior_alpha+posterior_beta)^2 * ( posterior_alpha + posterior_beta +1))) 


plot(drawn_mean, col="blue", xlab="Number of drawn samples")
abline(h=true_mean, lty = 2, lwd = 2, col ="red")

plot(drawn_sd,col="blue", xlab = "number of drawn samples")
abline(h=true_sd, lty = 2, lwd = 2, col ="red")



#----------------------------------1.1b-------------------


simulated_posterior_prob = mean(rbeta(10000, posterior_alpha,posterior_beta) <.3 ) #draw from posterior distribution, and check if smaller than 0.3

exact_posterior_prob = pbeta(.3, posterior_alpha, posterior_beta) #get the exact value using pbeta, 

print(paste("simulated:",simulated_posterior_prob,"exact:",exact_posterior_prob,sep=" "))



#---------------------------------1.1c------------------

log_odds <- function(p){
  return(log ( p/(1-p)))
} 


ndraws = 10000
posterior_theta_draws = rbeta(ndraws,posterior_alpha, posterior_beta) # again, draw from posterior distribution

posterior_logodds = log_odds(posterior_theta_draws)

hist(posterior_logodds, freq=F, main = "posterior distribution of phi = logodds(theta)")
lines(density(posterior_logodds))



#-----------------------------2a--------------------------

#priors
mu = 3.5
#sigma unknown with prior p(s^2) ??? 1/s^2

#data
y = c(33, 24, 48, 32, 55, 74, 23, 76 , 17)

n = length(y)
numdraws = 1E5

tau_sq = sum((log(y)-mu)^2)/n   
X_draws = rchisq(numdraws,n-1)

sigma_sq = (n-1)*tau_sq / X_draws   # Simulating draws from Inv?????2(n,??) distribution

# Theoretical density function for Inv?????2(n,??) distribution 
inversechisq <- function(x, v, s){
  density = (((v/2)^(v/2))/gamma(v/2)) * (s^v) * (x^(-(v/2+1))) * (exp((-v*s^2)/(2*x)))
  return (density)
}

# Plotting posterior draws from sampling 
plot(density(sigma_sq), lwd=3, main="Simulated and theoretical posteriors for sigma^2") 
# Plotting posterior distribution using therotical density function 
lines(seq(0.1,15,by=0.1), inversechisq(seq(0.1,15,by=0.1), n, sqrt(tau_sq)), type = "l", col =  "red", lty=2, lwd=3)  
legend("topright", inset=.02, c("Simulated Posterior Draws", "Theoretical Posterior Density"), fill=c("black", "red"))


#------------------------------2b-------------------------

gini = 2*pnorm(sqrt(sigma_sq)/sqrt(2))-1 #formula given in instructions
plot(density(gini),main="posterior distribution of Gini-coefficient")


#-----------------------------2c--------------------------

quants = quantile(gini,probs=c(0.025,0.975))

plot(density(gini),main="posterior distribution of Gini-coefficient")
abline(v=quants, col="green","lwd"=2,lty=2)


#----------------------------2d--------------------------

sorted = sort(gini)
# 2.d)
# Finding density of gini distribution
dense = density(sorted)
# storing coordinate and desnity values in a dataframe
df = data.frame(x=dense$x,y=dense$y)
# sorting dataframe based on density values
df = df[order(df$y),]
# Finding cumulative sum of sorted density values
df$cumsum = cumsum(df$y)
# finding all points at which the cumulative sum is > 5% if the total density ( 95% Highest posterior interval)
TOTAL = sum(df$y)
df_2 = df[df$cumsum>TOTAL*0.05,c("x","y")]

points(df_2,col="red")
legend("topright", inset=.01, c("95% Equal Tail Confidence Interval", "Highest Density Posterior
Interval"), fill=c("green", "red"))


#-----------------------------3a&b---------------------------



y = c(1.83, 2.02, 2.33, -2.79, 2.07, 2.02, -2.44, 2.14, 2.54, 2.23)


p_cond <- function(y,k,mu=2.51){
  n = length(y)
  res = (1/(2*pi*besselI(k,0))^n) * exp(sum(k*cos(y-mu))-k)
  return(res) 
}


probs = c()
part_length = 0.01
tested_k_values = seq(0,10,part_length) 
max_prob = 0
k_max_prob = 0
for (k in tested_k_values){
  p =p_cond(y,k)
  probs = c(probs,p)
  
  if (p > max_prob){
    max_prob = p
    k_max_prob = k
  }
}

plot(x =tested_k_values, y = probs / sum(probs * part_length),ylab = "density", xlab = "k", main = "posterior distribution of k, with mode plotted.", sub = paste("k_mode  =",k_max_prob))
abline(v=k_max_prob, col ="red", lty=2)

