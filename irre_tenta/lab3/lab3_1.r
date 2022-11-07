setwd("baysian-learning/lab3")
data = read.table("Precipitation.rds",stringsAsFactors =T, header=T, sep="",na.strings="")

data = readRDS("Precipitation.rds")

y = log(data)
n = length(y)

mu_0 = mean(y)
t_sq_0 = var(y)
sigma_sq_0 = var(y)
v0 = 1.5
#v0 = -2*sigma_sq_0/(t_sq_0-sigma_sq_0)
#v*t/v-2 = var(y)

mu_draws = c()
sigma_draws = c()

mu = rnorm(1,mu_0,t_sq_0)

X = rchisq(1,v0)
s_2 = v0*sigma_sq_0 / X+

y_mean = mean(y)
#sample_s = sqrt(var(y))
for( i in 2:5000){
  mu_draws = c(mu_draws,mu)
  sigma_draws = c(sigma_draws, s_2)
  
  w = (n/s_2) / ( (n/s_2) + (1 / t_sq_0))
  
  mu_n = w*y_mean + (1-w)*mu_0
  
  v_n = v0 + n

  t_sq_n = 1/(n / s_2) + (1 / t_sq_0)
  

  #mu
  mu = rnorm(1,mu_n,t_sq_n)
  
  X = rchisq(1,v_n)
  thing = (v0*sigma_sq_0 + sum((mu - y)^2))/(n + v0)
  s_2 = v_n*thing/X
  

}



plot(mu_draws,type="l")
plot(sigma_draws,type="l")