library(mvtnorm)

n = 250

#This can be varied to simulate different signals
alpha = c(rep(0.5,10), rep(0, 90))
beta = c(rep(0.5,5), rep(0, 195))

X1 = list()
X2 = list()

i=1
for(i in 1:n)
{
  print(i)
  #Simulate  times of measurement
  times = 1:10

  #random effect
  ui = rnorm(1,0,0.25)

  #Shared canonical vector (with some time effect)
  Zi = times * 0.5 + ui

  #Simulate data and add some noise
  X1i = sapply(alpha, function(a) rnorm(length(times), Zi * a, 0.5))
  X2i = sapply(beta, function(a) rnorm(length(times), Zi * a, 0.5))
  colnames(X1i) = paste0("X", 1:ncol(X1i))
  colnames(X2i) = paste0("Y", 1:ncol(X2i))

  #Check the simulated cross correlation
  #image(cor(X1i, X2i))

  #Remove some observations
  p_observed = 1
  X1i = cbind(i=i, times=times, X1i)[rbinom(length(times),1,p_observed)==1,]
  X2i = cbind(i=i, times=times, X2i)[rbinom(length(times),1,p_observed)==1,]

  X1[[i]] = X1i
  X2[[i]] = X2i
}

X1 = do.call("rbind", X1)
X2 = do.call("rbind", X2)
image(cor(X1,X2))

