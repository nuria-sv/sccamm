# simulations lower dim



# clean up and library ----------------------------------------------------

rm(list = ls())

library(MASS)
library(ggplot2)
library(corrplot)
library(ellipse)
library(rgl)

# set.seed(12345)
# simulations -------------------------------------------------------------

# environment variables
p1 = 3                                                                          # variables in X1
p2 = 2                                                                          # variables in X2
n1 = 10                                                                          # individuals in X1
n2 = n1                                                                         # individuals in X2
t1 = 4                                                                          # measurements in X1
t2 = t1                                                                         # mearuements in X2
k  = 1                                                                          # components
n  = max(n1, n2)
t = max(t1, t2)

p1.sel = sample(1:p1, 2, replace = F)
p2.sel = sample(1:p2, 2, replace = F)

# distribution variables
mu1 = rep(0, p1)
mu2 = mu1
c1 = matrix(runif(p1^2)*2-1, ncol=p1)
S1 = t(c1)%*%c1
c2 = matrix(runif(p2^2)*2-1, ncol=p2)
S2 = t(c2)%*%c2

S = magic::adiag(S1, S2)                                                        # block diagonal matrix of S1 and S2
# weights
w1 = rnorm(p1); w1[-p1.sel] = 0
w2 = rnorm(p2); w2[-p2.sel] = 0

W = matrix(c(w1, w2), nrow = (p1+p2), ncol = k)

# X1 = MASS::mvrnorm(n = n1, mu = mu1, Sigma = S1)



# latent ------------------------------------------------------------------
S_time = matrix(c(1, 0.7, 0.4, 0,
                  0.7, 1, 0.65, 0.4,
                  0.4, 0.65, 1, 0.65,
                  0, 0.4, 0.65, 1), ncol = t1)
z = MASS::mvrnorm(n = n, mu = rep(0, ncol(S_time)), Sigma = S_time)

matplot(t(z), type = "l")
lines(colMeans(z), col = "tomato", lty = 3, lwd = 3)

# observations ------------------------------------------------------------
X =  array(sapply(1:(p1 + p2), function(j) rnorm(n*t)), c(t, n, (p1 + p2)))
for (s in 1:t) {
  X[s,,] = MASS::mvrnorm(n = n, mu = W%*%z[,s], Sigma = W%*%(3)%*%t(W) + S)
}


par(mfrow = c(3, 4))
for (i in 1:(p1 + p2)) {
  matplot(t(X[,,i]), type = "l", ylab = ifelse(i<=p1, paste0("X1 ", i), paste0("X2 ", i - p1)), main = paste0("Var. ", ifelse(i<=p1, paste0("X1 ", i), paste0("X2 ", i - p1))))
}

# plots -------------------------------------------------------------------
# X1.mvt = mvtnorm::rmvnorm(n1, mean = mu1, sigma = S1)
ggplot(as.data.frame(X1), aes(x=V1, y=V2))+
  geom_point(alpha = .2) +
  geom_density_2d()+
  theme_bw()


corrplot(cor(X1),
         method="ellipse",
         tl.pos="n",
         title="Matrix Correlations")

plot(as.data.frame(X1),
     col = "black",
     main = "Bivariate Normal with Confidence Intervals")


ggplot(as.data.frame(X1), aes(x=as.data.frame(X1)[,1], y=as.data.frame(X1)[,2]))+
  geom_point(alpha = .5) +
  geom_bin2d() +
  scale_fill_viridis_c()+
  theme_bw()

ggplot(as.data.frame(X1), aes(x=as.data.frame(X1)[,1], y=as.data.frame(X1)[,2]))+
  geom_point(alpha = .2) +
  geom_density_2d()+
  theme_bw()

ggplot(as.data.frame(X1), aes(x=as.data.frame(X1)[,1], y=as.data.frame(X1)[,3]))+
  geom_point()+
  stat_ellipse(type = "norm", linetype = 2) +
  stat_ellipse(type = "t")+
  geom_smooth()+
  theme_bw()

bivn.kde <- kde2d(as.data.frame(X1)[,1], as.data.frame(X1)[,3], n = 50)
image(bivn.kde)
contour(bivn.kde, add = TRUE)

# Classic Bivariate Normal Diagram
dev.off()
rho <- cor(as.data.frame(X1))
y_on_x <- lm(as.data.frame(X1)[,2] ~ as.data.frame(X1)[,1])    # Regressiion Y ~ X
x_on_y <- lm(as.data.frame(X1)[,1] ~ as.data.frame(X1)[,2])    # Regression X ~ Y
plot(as.data.frame(X1)[,c(1,2)],
     col = "grey",
     main = "Bivariate Normal with Confidence Intervals")
lines(ellipse(rho), col="red")       # ellipse() from ellipse package
lines(ellipse(rho, level = .99), col="green")
lines(ellipse(rho, level = .90), col="blue")
abline(y_on_x)
abline(x_on_y, col="brown")


# Basic perspective plot
persp(bivn.kde, phi = 10, theta = 30, shade = .1, border = NA) # from base graphics package
# RGL interactive plot
col2 <- heat.colors(length(bivn.kde$z))[rank(bivn.kde$z)]
persp3d(x=bivn.kde, col = col2)

