require(mvtnorm)

nsim <- 1000
betas <- c(1,2)
nbetas <- length(betas)
nsamples <- 100
estimates <- matrix(nrow=nsim, ncol=nbetas)

par(mfcol=c(1,3))

####################################################
# simulation 1. Noise is fixed. Design matrix varies

noise <-scale(rnorm(nsamples))

for (sim in 1:nsim) 
{
   # X = cbind(scale(rnorm(nsamples)), scale(rnorm(nsamples)))
   X= rmvnorm(nsamples, mean = c(0, 0), sigma = matrix(c(1,.8,.8,1),nrow=2),method="chol")
   y <- X %*% betas + rnorm(nsamples)
   estimates[sim,] <- coef(lm(y~X))[-1]
}
boxplot(estimates,main="Simulation 1",ylim=c(0,3))
grid()


#######################################################
# simulation 2. Noise is varied. Design matrix constant

X = cbind(scale(rnorm(nsamples)), scale(rnorm(nsamples)))

for (sim in 1:nsim) 
{
  noise <-rnorm(nsamples)  
  y <- X %*% betas + scale(rnorm(nsamples))
  estimates[sim,] <- coef(lm(y~X))[-1]
}
boxplot(estimates,main="Simulation 2",ylim=c(0,3))
grid()

#######################################################
# Simulation 3. Both noise & design matrix vary


for (sim in 1:nsim) 
{
  noise <-scale(rnorm(nsamples))
  X = cbind(scale(rnorm(nsamples)), scale(rnorm(nsamples)))
  
  y <- X %*% betas + rnorm(nsamples)
  estimates[sim,] <- coef(lm(y~X))[-1]
}
boxplot(estimates,main="Simulation 3",ylim=c(0,3))
grid()


