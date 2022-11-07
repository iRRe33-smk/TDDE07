mu = 13
sigma_sq = 3
t_max = 300
fi = 1
# lab 3.3.1

# Function to simulate Auto Regressive Function
AR_Process <- function(mu, FI, sigma_sq)
{
  path = matrix(data = NA, nrow = length(FI), ncol = t_max)
  i = 1
  for (fi in FI){
    path[i,1] = mu
    for (t in 2:t_max){
      path[i,t] = mu + fi*(path[i,t-1]-mu) + rnorm(1,0,sigma_sq)
    } 
    i = i+1
  }
  return(path)
}

# Simulating Draws from the AR1 process with phi value -0.8, 0, 0.8
FI = seq(-.8,.8,.8)
path = AR_Process(mu, FI, sigma_sq)
par(mfrow=c(3,2))
for (fi in 1:length(FI)){
  plot(path[fi,],main=paste("phi = ",FI[fi]),type="line",col="blue")  
  acf(path[fi,],main=paste("ACF phi = ",FI[fi]))
}

# lab 3.3.2 

# Creatind Data to pass to stan model
fi_1 = .2
fi_2 = .95

FI = c(fi_1, fi_2)
path = AR_Process(mu, FI, sigma_sq)
par(mfrow=c(2,2))
for (fi in 1:length(FI)){
  plot(path[fi,],main=paste("phi = ",FI[fi]),type="line",col="blue")  
  acf(path[fi,],main=paste("ACF phi = ",FI[fi]))
}

library(rstan)

# Creating Stan Model to simulate the two AR(1)-processes

StanModel = '
data {
  matrix[300,2] M; // synthetic data
  
}
parameters {
  real mu;
  real<lower=0> sigma2;
  real phi_1;
  real phi_2;
}
model {
  mu ~ uniform(-100,200); // Normal with mean 0, st.dev. 100
  sigma2 ~ uniform(0,100); // Scaled-inv-chi2 with nu 1,sigma 2
  phi_1 ~ uniform(-10,10);
  phi_2 ~ uniform(-10,10);
  
  for(i in 2:300){
    M[i,1] ~ normal(mu + phi_1*(M[i-1,1]-mu),sigma2);
    M[i,2] ~ normal(mu + phi_2*(M[i-1,2]-mu),sigma2);   
  }
}'

data <- list(M=t(path))
warmup <- 10
niter <- 2000
fit <- stan(model_code=StanModel,data=data, warmup=warmup,iter=niter,chains=4)

# Print the fitted model
print(fit,digits_summary=3)

# Extract posterior samples
postDraws <- extract(fit)

# Do traceplots of the first chain
par(mfrow = c(1,1))

# Do automatic traceplots of all chains
traceplot(fit)
# Plotting joint posterior for mu and phi1
plot(postDraws$mu[1:(niter-warmup)], postDraws$phi_1[1:(niter-warmup)], xlab = "mu", ylab = "phi_1", main = "Joint Posterior of mu and phi 1")
# Plotting joint posterior for mu and phi2
plot(postDraws$mu[1:(niter-warmup)], postDraws$phi_2[1:(niter-warmup)], xlab = "mu", ylab = "phi_1", main = "Joint Posterior of mu and phi 2")
# Bivariate posterior plots
pairs(fit)
























