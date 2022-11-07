



mu = 13
sigma_sq = 3
t_max = 300
fi = 1
# ------------------ 3.3.1
#num
par(mfrow=c(3,2))
for (fi in seq(-.8,.8,.8)){
  path = matrix(nrow = t_max)
  path[1] = mu
  for (t in 2:t_max){
    path[t] = mu + fi*(path[t-1]-mu) + rnorm(1,0,sigma_sq)
  } 
  plot(path,main=fi,type="line",col="blue")  
  acf(path,main=paste("ACF ",fi))
}


#----------------3.3.2 create data

fi_1 = .2
fi_2 = .95
par(mfrow=c(3,1))
for (fi in seq(-.8,.8,.8)){
  path = matrix(nrow = t_max,ncol=2)
  path[1,] = c(mu,mu)
  for (t in 2:t_max){
    path[t,1] = mu + fi_1*(path[t-1,1]-mu) + rnorm(1,0,sigma_sq)
    path[t,2] = mu + fi_2*(path[t-1,2]-mu) + rnorm(1,0,sigma_sq)
    } 
  plot(path[,1],type="line",col="blue")  
  lines(path[,2],col="red")#acf(path,main=paste("ACF ",fi))
}


write.table(path,file="paths.txt",header=T)


library(rstan)
#y=c(4,5,6,4,0,2,5,3,8,6,10,8)
#N=length(y)

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
  mu ~ uniform(-100,100); // Normal with mean 0, st.dev. 100
  sigma2 ~ uniform(0,100); // Scaled-inv-chi2 with nu 1,sigma 2
  phi_1 ~ uniform(-10,10);
  phi_2 ~ uniform(-10,10);
  
  for(i in 2:300){
    M[i,1] ~ normal(mu * phi_1*(M[i-1,1]-mu),sqrt(sigma2));
    M[i,2] ~ normal(mu * phi_2*(M[i-1,2]-mu),sqrt(sigma2));   
  }
}'

data <- list(M=path)
warmup <- 1000
niter <- 2000
fit <- stan(model_code=StanModel,data=data, warmup=warmup,iter=niter,chains=4)
# Print the fitted model
print(fit,digits_summary=3)
# Extract posterior samples
postDraws <- extract(fit)
# Do traceplots of the first chain
par(mfrow = c(1,1))
plot(postDraws$mu[1:(niter-warmup)],type="l",ylab="mu",main="Traceplot")
# Do automatic traceplots of all chains
traceplot(fit)
# Bivariate posterior plots
pairs(fit)
























