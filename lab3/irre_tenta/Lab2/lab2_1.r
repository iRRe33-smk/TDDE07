library(mvtnorm)
library(scales)

# Lab 2.1.a)

# Initializing hyper parameters of the linear regression model
mu0 = c(-10, 100, -100)
omega0 = 0.12 * diag(3)
v0 = 3
sigma0_sq = 2

# Number of elements to sample for σ^2 and β
N <- 1E2

# Reading data from the given txt file
data <- read.table("TempLambohov.txt", sep = "" , header = T , na.strings ="", stringsAsFactors= F)

# Creating a matrix of 1, time, time^2 from the data read
X = data.frame(const = 1, t = data$time, t2 = data$time^2)

# creating a sampling of 201 values between 0 and 1 for the time variable
time_sample = data.frame(const=1,t = seq(0,1,0.005),t2 =seq(0,1,0.005)^2)

# function to calculate temperate based on the linear regression formula β0 + β1 * time + β2 *time^2
pred_temp <- function(betas){
  pred_temp = as.matrix(time_sample)%*%t(betas)
}

# Sampling for prior σ^2 from Inv−χ2(v0,σ0^2) distribution
X_draws = rchisq(N,v0)                # Taking N samples from the χ2(v0)
sigma_sq = (v0)*sigma0_sq / X_draws   # Simulating draws from Inv−χ2(v0,σ0^2) distribution

# Plotting data points
plot(data,ylim=c(-25, 60))


betas_prior = matrix(nrow=N,ncol=3 )
x_plot = seq(0,1,0.005) 
for (i in 1:N){
  # Simulating draws for prior β from N(μ0, σ^*2Ω0^−1)
  beta = rmvnorm(1, mean = mu0, sigma = sigma_sq[i]*solve(omega0))
  betas_prior[i,]=beta
  # Plotting temperatures predicted from the prior draws
  points(x_plot,t(pred_temp(beta)),col=alpha("red",0.07))
}

# Lab 2.1.b)

# Function to simulate for σ^2 and β from the joint posterior distribution
posterior <- function(X,y, omega0,  mu_0, v0, n, sigma0_sq ){
  posterior_parameters = matrix(nrow=N,ncol=4)
  
  # calculating beta_hat = ((X′X)^−1)X′y 
  beta_hat = solve(t(X)%*%X)%*%t(X)%*%y
  # calculating μn = ((X′X+Ω0)^−1)(X′X*beta_hat + Ω0*μ0)
  mu_n = solve(t(X)%*%X + omega0) %*% (t(X)%*%X%*%beta_hat + omega0%*%mu_0)
  # calculating Ωn = X′X+Ω0
  omega_n = t(X)%*%X + omega0
  # calculating νn = ν0 + n
  v_n = v0 + n
  # calculating v_s = vn*σn^2 = v0 * σ0^2 + (y′y + μ0′* Ω0 * μ0 − μn′* Ωn * μn)
  v_s = v0*sigma0_sq + (t(y)%*%y + t(mu_0)%*%omega0%*%mu_0 - t(mu_n)%*%omega_n%*%mu_n)
  
  # Sampling for posterior σ^2 from Inv−χ2(vn,σn^2) distribution
  X_draws = rchisq(N,v_n)    # Taking N samples from the χ2(vn)
  sigma_sq = v_s / X_draws   # Simulating draws from Inv−χ2(vn,σn^2) distribution
  
  posterior_parameters[,4] = sigma_sq
  temps = matrix(nrow=N,ncol=201)
  for (i in 1:N){
    # Simulating draws for posterior β from N(μn, σ^*2Ωn^−1)
    beta = rmvnorm(1, mean = mu_n, sigma = sigma_sq[i]*solve(omega_n))
    posterior_parameters[i,1:3]=beta
    
    # Predicting temperature values from the posterior β
    temps[i,] = t(pred_temp(beta))
    
    # plotting predcited temperature values
    points(x_plot,t(pred_temp(beta)),col=alpha("blue",0.01))
    
  }
  ret = list(posterior_parameters=posterior_parameters, temps=temps)
  return(ret)
}
legend("topleft", inset=.02, c("Tempature Data", "Calulated Prior Temperature", "Caulated posterior Temperature"), fill=c("black", "red", "blue"))
ret = posterior(as.matrix(X),as.matrix(data$temp),omega0, mu0,v0,dim(data)[1],sigma0_sq)
betas = ret$posterior_parameters[,1:3]
sigma_sq = ret$posterior_parameters[,4]
temps = ret$temps

# # 2.1.b.i)
# # Plotting histogram for each of the marginal posterior parameters
hist(betas[,1], freq = FALSE, main = "Histogram of Beta 0")
hist(betas[,2], freq = FALSE, main = "Histogram of Beta 1")
hist(betas[,3], freq = FALSE, main = "Histogram of Beta 2")
hist(sigma_sq, freq = FALSE, main = "Histogram of Sigma^2")



# 2.1.b.ii)

# Scatter plot of temperature value
plot(data,ylim=c(-25, 60))

post_medians = c()
quants = matrix(nrow=dim(temps)[2],ncol = 2)
for (i in 1:dim(temps)[2]){
  # Finding the median of the predicted temperatures
  med = median(temps[,i])
  # finding the 95% equal tail posterior probablity interval of predicted terperature
  q = quantile(temps[,i],probs=c(.025,0.975))
  q = unname(q)

  quants[i,] = c(q[1],q[2])

  post_medians = c(post_medians,med)
}
lines(x_plot,post_medians, col = "green")
lines(x_plot,quants[,1], col = "blue")
lines(x_plot,quants[,2], col = "red")

# 2.1.c)

# The time at which the temperature is maximum can be found by solveing the function f'(time) = 0 and solving for time
# This points can be considered as the maximas as f''(time) = 2 * β2 < 0 for all values in β2
# hence time at which temperature is max = -β1/(2*β2)


max_temp_time = -betas[,2]/(2*betas[,3])
abline(v=max_temp_time, col=alpha("violet",0.1))
legend("topleft", inset=.02, c("Tempature Data", "Upper limit 95% equal tail credible intervel", 
                                   "lower limit 95% equal tail credible intervel", "Median Posterior Temperature", "Max Temperature time period"), 
       fill=c("black", "red", "blue", "green", "violet"))


