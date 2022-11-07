y = c(33,24,48,32,55,74,23,76, 17)
var_measure = var(y)
mu = 3.5
plot(density(rnorm(1000,mu,sqrt(var_measure))))

n = 9
tau_sq = sum((log(y)-mu)^2)/n

X_draws = c()
theta_draws = c()
print("start loop")
for (tmp in 1:1E5)
{

  X = rchisq(1,n-1)
  X_draws = c(X_draws,X)
  
  sigma_sq = (n-1)*tau_sq / X
  theta_draw = rnorm(1,mean(X_draws),sigma_sq/n)
  theta_draws = c(theta_draws,theta_draw)
  }
print("end loop")



plot(density(theta_draws))

#salary_draws = rnorm(1E5,mu,sqrt(8))
#hist(exp(salary_draws))
#mean(salary_draws)
