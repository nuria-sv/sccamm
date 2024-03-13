# toscca ME linear
toscca.core = function(alphaInit, A, B, nonzero_a, nonzero_b, iter = 20, tol = 10^(-6), silent = FALSE, lmeformula = " ~ -1 + time + (1|id)")
{

  # checks
  if(ncol(B) <= max(nonzero_b)) {
    message("At least one of the nonzero options for B is not sparse. Changing to meet criteria")
    nonzeroB = nonzero_b[nonzero_b < (ncol(B) - 2)]
  }

  if(ncol(A) <= max(nonzero_a)) {
    message("At least one of the nonzero options for A is not sparse. Changing to meet criteria")
    nonzero_a = nonzero_a[nonzero_a < (ncol(A) - 2)]
  }


  #Create the matrix A
  alpha = sapply(nonzero_a, function(x) c(alphaInit))

  varTol1 = matrix(0, nrow = nrow(A), ncol = length(nonzero_a))
  varTol2 = matrix(0, nrow = nrow(B), ncol = length(nonzero_b))
  i = 0
  e = 10


  # format data
  id_a = A[,"id"]
  id_b = B[,"id"]
  time_a = A[,"time"]
  time_b = B[,"time"]

  A = as.matrix(A[, 3:ncol(A)])
  B = as.matrix(B[, 3:ncol(B)])

  while (e > tol & i <= iter) {
    i = i +1

    # refresh
    if(i > 1) varTol1 = gamma
    if(i > 1) varTol2 = zeta


    gamma =  A %*% alpha
    dist  = sqrt(colSums(gamma^2))
    gamma = sweep(gamma, 2, dist, "/")


    me = lmer(as.formula(paste("gamma", lmeformula)), data = data.frame(gamma = gamma, time = time_a, id = id_a), REML = TRUE)
    pred_me = predict(me, newdata = data.frame(time = time_b, id = id_b), allow.new.levels = TRUE, re.form = NULL)

    beta = t(B) %*% pred_me

    rm(me, pred_me)

    beta = apply(rbind(beta,nonzero_b), 2, function(x)
    {
      nonzero1 = x[length(x)]
      y = x[-length(x)]
      thres = abs(y)[order(abs(y), decreasing=TRUE)[nonzero1+1]]
      tmp = (abs(y) - thres)
      tmp[tmp<=0] = 0
      sign(y) * tmp
    })

    zeta = B %*% beta
    dist = sqrt(colSums(zeta^2))
    zeta = sweep(zeta, 2, dist, "/")


    me = lmer(as.formula(paste("zeta", lmeformula)), data = data.frame(zeta = zeta, time = time_b, id = id_b), REML = TRUE)
    pred_me = predict(me, newdata = data.frame(time = time_a, id = id_a), allow.new.levels = TRUE, re.form = NULL)

    alpha = t(A) %*% pred_me

    rm(me, pred_me)

    alpha = apply(rbind(alpha,nonzero_a), 2, function(x)
    {
      nonzero1 = x[length(x)]
      y = x[-length(x)]
      thres = abs(y)[order(abs(y), decreasing=TRUE)[nonzero1+1]]
      tmp = (abs(y) - thres)
      tmp[tmp<=0] = 0
      sign(y) * tmp
    })


    if(length(nonzero_a) == 1) e = mean(abs(gamma - varTol1)) + mean(abs(zeta - varTol2))
    if(length(nonzero_a) > 1) e  = mean(colMeans(abs(gamma - varTol1))) + mean(colMeans(abs(zeta - varTol2)))

    textSCCA = paste0(" Common convergence error: ", round(e, 5), " & Iterations: ", i)
    if(isFALSE(silent) & (e<= tol || i > iter)) cat(textSCCA, "\r")

  }
#
#   dist = sqrt(colSums(alpha^2))
#   alpha = sweep(alpha, 2, dist, "/")
#
#   dist = sqrt(colSums(beta^2))
#   beta = sweep(beta, 2, dist, "/")

  return(list(a = alpha, b = beta, conv = e, iter = i))
}
