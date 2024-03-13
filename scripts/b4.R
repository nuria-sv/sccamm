
# library -----------------------------------------------------------------
library(runjags)
library(ggplot2)
library(reshape2)
library(toscca)
library(RColorBrewer)
library(hrbrthemes)

# model -------------------------------------------------------------------
dpcca.mod = "model
{
  # ########## #
  #   latent   #
  ##############

  ## priors ##

  v.inv ~ dgamma(0.001,0.001); # for precision
  v <- 1/v.inv;                 # variance ar lantent err
  phi ~ dunif(-1,1);            # one time coeff
  #a ~ dnorm(0, 0.01);           # level coeff
  b ~ dnorm(0, 1);           # det time coeff

  for(i in 1:n) {
    # mu[i] ~ dnorm(0, 0.01);
    z0[i] ~ dnorm(0, 0.001); # initial state
  }



  ## likelihood ##
  for(i in 1:n) {
      Z[1,i] ~ dnorm(phi*z0[i], v.inv);

  }

  for(s in 2:t) {
    for(i in 1:n) {
      Z[s,i] ~ dnorm(phi*Z[s-1, i] + b*(s/t), v.inv);

    }
  }

  # ################ #
  #   observations   #
  ####################

  ## priors ##
  for(j in 1:p) {
    w1[j] ~ dnorm(0, prec1[j]);
    prec1[j] ~ dgamma(1/2, gamprec1[j]);
    gamprec1[j] ~ dgamma(10^4,10^4);
    S1.inv[j] ~ dgamma(0.001, 0.001);
    S1[j] <- 1/S1.inv[j];
  }

  for(j in 1:q) {
    w2[j] ~ dnorm(0, prec2[j]);
    prec2[j] ~ dgamma(1/2, gamprec2[j]);
    gamprec2[j] ~ dgamma(10^4,10^4);
    S2.inv[j] ~ dgamma(0.001, 0.001);
    S2[j] <- 1/S2.inv[j];
  }

  ## likelihood ##
  for(i in 1:n) {
    for(s in 1:t) {
      for(j in 1:p) {
        X1[s, i, j] ~ dnorm(Z[s,i]*w1[j], S1[j]);
      }
      for(j in 1:q) {
        X2[s, i,j] ~ dnorm(Z[s,i]*w2[j], S2[j]);
      }
    }
  }


}"

  # simulatios --------------------------------------------------------------


  n = 100
  p = 100
  q = 100
  K = 1
  t = 10

  p.sel = sample(1:p, 2, replace = F)
  q.sel = sample(1:q, 1)
  w1 = c(-0.9, 1.3)
  w2 = 2.5

  # noise
  X1 = array(sapply(1:p, function(j) rnorm(n*t)), c(t, n, p))
  X2 = array(sapply(1:q, function(j) rnorm(n*t)), c(t, n, q))

  # latent variable
  z0  = 0                                 # initial value
  a   = 0                                 # level
  b   = 0.2                               # det. trend
  phi = 0.0                               # ar(1)
  eps = matrix(rnorm(t*n, 0, 0.01), t, n)          # error term
  mu  = matrix(0, t, n)                   # sapply(1:n, function(i) runif(1, -1, 1))

  Z  = matrix(z0, t, n)
  for (s in 2:t) {

    # Z[s, ] = tr*(s/t) + b*Z[s-1]

    for (i in 1:n) {
      mu[s,i] =  a + b*(s/t)
      Z[s, i] = phi*Z[s -1 ,i] + mu[s,i] + eps[s,i]

    }


    X1[s, , p.sel[1]] = w1[1]*Z[s,] + rnorm(1, 0, 0.01)
    X1[s, , p.sel[2]] = w1[2]*Z[s,] + rnorm(1, 0, 0.01)
    X2[s, , q.sel[1]] = w2[1]*Z[s,] + rnorm(1, 0, 0.01)
  }

  data = list(X1 = X1, X2 = X2, t = t, n = n, p = p, q = q)


  rownames(Z) = paste0("t", 1:t)
  colnames(Z) = paste0("i", 1:n)
  Z = data.frame(Z)
  Z$t =  1:t
  Z.melt = melt(Z, id.var = "t"); colnames(Z.melt) <- c("t", "i", "z")

  ggplot(Z.melt, aes(x=t,y=z,group=i,colour=i)) +
    geom_line() +
    scale_x_discrete(limits = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10")) +
    ggtitle("Latent process") +
    xlab("Time") +
    ylab("z") +
    theme_ipsum()


  # results ----------------------------------------------------------------------
  mod = run.jags(data = data, model=dpcca.mod, monitor=c("phi", "b", "w1", "w2", "S1", "S2"), n.chains=3,
                 method="rjparallel", burnin = 100, sample = 1000, adapt = 100)



  # plots -------------------------------------------------------------------

  aa = do.call("rbind", mod$mcmc)

  colp = rep(1, p); colp[p.sel] <- 2
  colq = rep(1, q); colq[q.sel] <- 2


  par(mfrow=c(1,2))

  plot(apply(aa, 2, median)[3:(p+3)], xlab = "p", ylab = "W1", col = colp, pch = colp)

  plot(apply(aa, 2, median)[(p+3):(p+q + 2)], xlab = "q", ylab = "W2", col = colq, pch = colq)

  apply(aa, 2, median)[-(1:(2*p+2*q+2))]



  colnames(aa)

  ind_param = c(1,2, p.sel+2, p +q.sel+2, p+9+2) # -((p+q):(p+q+2))
  true_param = abs(c(phi, b,  w1, w2, 0))

  bb1 = abs(mod$mcmc[[1]][,ind_param])

  bb2 = abs(mod$mcmc[[2]][,ind_param])
  bb3 = abs(mod$mcmc[[3]][,ind_param])

  par(mfrow=c(3,2), ask=TRUE)

  for(i in 1:length(ind_param))

  {

    plot(as.matrix(bb1[,i]), type="p", main=colnames(bb1)[i], ylim = c(min(bb1[,i], bb2[,i], bb3[,i], true_param[i], 0), max(bb1[,i], bb2[,i], bb3[,i], true_param[i], 0.95)))

    points(as.matrix(bb2[,i]), col=2)
    points(as.matrix(bb3[,i]), col=3)
    abline(h=true_param[i], col = "blue")

  }

  Z_hat  = matrix(0, t, n); mu_hat = matrix(NA, t, n)


  Z_hat  = matrix(0, t, n); mu_hat = matrix(NA, t, n)
  for (s in 2:t) {
    for (i in 1:n) {
      mu_hat[s,i] =  0 + par_hat[1]*(s/t)
      Z_hat[s, i] =  + mu_hat[s,i]
    }
  }
