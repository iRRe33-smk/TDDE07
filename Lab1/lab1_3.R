

#3.a

# Function for the posterior distribution with n data points
p_cond <- function(y,mu,k){
  n = length(y)
  res = (1/(2*pi*besselI(k,0))^n) * exp(sum(k*cos(y-mu))-k)
  return(res) 
}



y= c(1.83, 2.02, 2.33, -2.79, 2.07, 2.02, -2.44, 2.14, 2.54, 2.23)

# Simulating 15000 posterior distribution over a k grid of 0-10 
N=1500
res = c()
used_k = c()
for (k in seq(0,10*N)){
  res = c(res,p_cond(y,2.51,k/N))
  used_k = c(used_k,k/N)
}

# Normalizing density values by dividing by total sum of densities
res_norm = res/sum(res)
# Plotting k values against densities
plot(used_k,res_norm, main = "Posterior Distribution of k Values")

# The mode of a distribution with a discrete random variable is the value of the term that occurs the most often
# i.e the value of k with the highest probability density
used_k[res==max(res)]
abline(v=used_k[res==max(res)], col = "red", type="l", lty=2)
legend("topright", inset=.01, c("Posterior mode of k"), fill=c("red"), )

#dense = density(res)





