
####################################################################################
####################################################################################
###                                                                              ###
### penalized CCA of X and Y-matrices with sets of x- and y-variables measured   ###
###  at varying numbers of time-points with the criss-cross algorithm            ###
###                                                                              ###
####################################################################################
####################################################################################
scalar1 <- function(x) {x / sqrt(sum(x^2))}
# read the simulated data, here all X and Y are still available in all n for all time-points
# load("/users/ahzwinderman/desktop/nuria/herhaalde metingen syntaxen/nuria1.Rdata")
#rm(input2)

# re-organize d2$X and d2$Y in long-format and throw away some rows to simulate a situation of n patients measured varying numbers of times
dim(d2$X)
dim(d2$Y)
T = dim(d2$X)[1]
n = dim(d2$X)[2]
XX = data.frame(id=1:n, tijd=1, d2$X[1,,])
YY = data.frame(id=1:n, tijd=1, d2$Y[1,,])
for (t in 2:T) {
   XX=rbind(XX, data.frame(id=1:n, tijd=t, d2$X[t,,]))
   YY=rbind(YY, data.frame(id=1:n, tijd=t, d2$Y[t,,]))
}
dim(XX)
dim(YY)
class(XX)
class(YY)
names(XX)
names(YY)
table(XX$tijd)
table(YY$tijd)
table(XX$id)
table(YY$id)

verwijder_X = sample(1:nrow(XX), 0.2*nrow(XX), replace=FALSE)
verwijder_Y = sample(1:nrow(YY), 0.3*nrow(YY), replace=FALSE)
table(verwijder_X)
table(verwijder_Y)
XX=XX[-verwijder_X,]
YY=YY[-verwijder_Y,]

dim(XX)
dim(YY)
table(XX$tijd)            # nog steeds 1 ... 10
table(YY$tijd)            # nog steeds 1 ... 10
table(XX$id)
length(table(XX$id))   # nog steeds 100
table(YY$id)
length(table(YY$id))   # nog steeds 100

rm(n, t, T, verwijder_X, verwijder_Y)
rm(d2)

ls()
c(dim(XX), dim(YY))
table(XX$tijd)
table(YY$tijd)
table(XX$id)
table(YY$id)

### data are ready; 2 data.frames; XX (800 x 2502) and YY (700 x 1302)
XX2=XX
YY2=YY
# XX2=scale(XX2)  # scale sets? How?
# YY2=scale(YY2)

colnames(XX2)[1:2] <- c("id", "time")
colnames(YY2)[1:2] <- c("id", "time")

#########################################################################################################################################################################################################


### HMP data: variables from db4 and variables from metabolomics2
# datasets XX and YY are special variants of X and Y, such that only time-points are selected with X- and Y-measurements at the time-point in the same subject (so XX and YY are paired, if you like)
# load("/users/ahzwinderman/desktop/nuria/herhaalde metingen syntaxen/HMP_db4vars_metabolomics2vars.Rdata")
c(dim(X), dim(XX), dim(Y), dim(YY))

#########################################################################################################################################################################################################

###### First run the functions defined below from line 167 onwards

#X=scale(X)   # scale the 2 matrices ?? If so, How?
#Y=scale(Y)
#XX=scale(XX) # same pats in XX and YY
#YY=scale(YY)

##################
## simulated data
##################

# estimate random-effects of each x- and y-variable by separate LME's, eigenlijk alleen zinvol als dezelfde n pats in X en Y zitten
date()
RE_XX = do_separate_lmes(XX2, lmeformule=" ~ poly(time,1) + (1+time|id)")
RE_YY = do_separate_lmes(YY2, lmeformule=" ~ poly(time,1) + (1+time|id)")
date()    # circa two minutes

# calculate cross-correlations of the estimated random-effects
ff1=cor(RE_XX, RE_YY)

# plot the histogram of the cross-correlations and also determine the histogram of ncol(RE_XX)*ncol(RE_YY) correlations drawn from a null-distribution
ff1a=as.numeric(ff1)                                          # 2*2500 * 2*1300 = 13 million correlations, but 1573926 NAs. Only 20*10 = 100 correlations are truly unequal to zero
ff2=hist(ff1a,plot=FALSE)
ff3=rnorm(sum(!is.na(ff1a)),0,1/sqrt(nrow(RE_XX)-3))          # draw about 13M Fisher-transforms of correlations from the normal distribution with zero mean and variance 1/(n-3)
ff3=(exp(2*ff3)-1) / (exp(2*ff3)+1)                           # calculate the inverse Fisher-transforms
ff4=hist(ff3, plot=FALSE,breaks=ff2$breaks)                   # histogram of the randomly selected correlations
kk=ncol(RE_XX)*ncol(RE_YY)
alpha=0.05
plot(0,0, type="n", xlim=c(min(ff2$breaks),max(ff2$breaks)), ylim=c(-max(ff2$density,ff4$density),max(ff2$density,ff4$density)), xlab="", ylab="density")#, yaxt="n")
i=2
for (i in 2:length(ff2$breaks)) {
   polygon(x = c(ff2$breaks[(i-1)], ff2$breaks[(i-1)], ff2$breaks[i], ff2$breaks[i],ff2$breaks[(i-1)]), y=c(0,ff2$density[(i-1)],ff2$density[(i-1)],0,0), col="grey")
   polygon(x = c(ff4$breaks[(i-1)], ff4$breaks[(i-1)], ff4$breaks[i], ff4$breaks[i],ff4$breaks[(i-1)]), y=c(0,-ff4$density[(i-1)],-ff4$density[(i-1)],0,0), col="pink")
}
abline(h=0, col=1)
text(x=min(ff2$breaks),y=max(ff2$density,ff4$density),"histogram of observed cross-correlations", col=1, adj=0)
text(x=min(ff2$breaks),y=-max(ff2$density,ff4$density),"histogram of cross-correlations under the null-hypothesis of no association", col="red", adj=0, cex=0.7)

# some stats of the distribution of correlations (compared to stats of the null-distribution)
round(c(mean(ff1a,na.rm=TRUE), sd(ff1a,na.rm=TRUE), quantile(ff1a, probs=c(0.025,0.25,0.5,0.75,0.0975),na.rm=TRUE)),6)
round(c(mean(ff3), sd(ff3,na.rm=TRUE), quantile(ff3, probs=c(0.025,0.25,0.5,0.75,0.0975),na.rm=TRUE)),6)
which(abs(as.numeric(ff1))>0.9)

# do a standard penalized CCA with the criss-cross/NIPALS algorithm on RE_XX and RE_YY
date()
res0 = estimate0_ab(RE_XX, RE_YY, alpha=1, lambda=0.00001, eps = 0.001, maxiter=1000)
date()   # about 2 minutes and 30 seconds
c(res0$iter)
makeplots0(RE_XX, RE_YY, res0, weigths=TRUE, loadings=TRUE)

# do the criss-cross algorithm on the original X- and Y-matrices
descriptives(XX2, YY2, plotit=FALSE, k=16)
res4a = estimate4_ab(XX2, YY2, lmeformule=" ~ -1 + poly(time,3) + (1+time|id)", alpha=0.99, lambda=0.001, eps = 0.001, maxiter=1000)
c(res4a$conv, res4a$iter)
makeplots(XX2, YY2, res4a, lmeformule=" ~ poly(time,3) + (1+time|id)", weigths=TRUE, loadings=TRUE, change=TRUE);dev.new()
plotsofpatientsof1variable(XX2, set="X", number=1, res4a, k=16, lmeformule=" ~ poly(time,3) + (1|id)"); dev.new()
plotsofselectedvariables(XX2, set="X", variables=1:20, res4a, lmeformule=" ~ poly(time,3) + (1|id)"); dev.new()
plotsofselectedvariables(YY2, set="Y", variables=1:10, res4a, lmeformule=" ~ poly(time,3) + (1|id)"); dev.new()

# TOSCCA
source("C:/Users/PC/OneDrive/github/ccalb/scripts/toscca_me.R")
res_toscca = toscca.core(alphaInit = runif(ncol(XX2)-2), XX2, YY2, 50, 50, lmeformula = " ~ -1 + poly(time,3) + (1|id)")
c(res_toscca$conv, res_toscca$iter)
makeplots(XX2, YY2, res_toscca, lmeformule=" ~ poly(time,3) + (1+time|id)", weigths=TRUE, loadings=TRUE, change=TRUE); dev.new()
plotsofpatientsof1variable(XX2, set="X", number=1, res_toscca, k=16, lmeformule=" ~ poly(time,3) + (1|id)"); dev.new()
plotsofselectedvariables(XX2, set="X", variables=1:20, res_toscca, lmeformule=" ~ poly(time,3) + (1|id)"); dev.new()
plotsofselectedvariables(YY2, set="Y", variables=1:10, res_toscca, lmeformule=" ~ poly(time,3) + (1|id)")

##################
## HMP data
##################

# NURIA'S HMP DATA
load("D:/long_cca/HMP/data_restored_id.RData")
X.original = X2
Y.original = X1
id_unique = data.frame(id_name = unique(X.original$d2_SubjectID), id= 1:length(unique(X.original$d2_SubjectID)))

X.temp = merge(X.original, id_unique, by.x = "d2_SubjectID", by.y = "id_name"); X.temp$id_name <- NULL
Y.temp = merge(Y.original, id_unique, by.x = "d2_SubjectID", by.y = "id_name"); Y.temp$id_name <- NULL
X.temp$time <- ave(X.temp$id, X.temp$id, FUN=seq_along)
Y.temp$time <- ave(Y.temp$id, Y.temp$id, FUN=seq_along)

X.temp$id = as.numeric(X.temp$id); X.temp$time = as.numeric(X.temp$time); X.temp$d2_SubjectID <- NULL; X.temp$VisitID <- NULL
Y.temp$id = as.numeric(Y.temp$id); Y.temp$time = as.numeric(Y.temp$time); Y.temp$d2_SubjectID <- NULL; Y.temp$VisitID <- NULL

X.temp = na.omit(X.temp); Y.temp = na.omit(Y.temp)

X = X.temp[X.temp$id %in% Y.temp$id,c("id", "time", colnames(X.temp)[1:(ncol(X.temp) - 2)])]; rm(X.temp)
Y = Y.temp[Y.temp$id %in% X$id,c("id", "time", colnames(Y.temp)[1:(ncol(Y.temp) - 2)])]; rm(Y.temp)


# calculate a few descriptive statistics and plots of k randomly selected X- and Y-variables with lme-estimated means
descriptives(X, Y, plotit=FALSE, k=16)

# estimate the weights for the X- and Y-variables using the criss-cross/NIPALS algorithm using some lme-model to link the latent variables of the X- or Y-variables
date()
res4 = estimate4_ab(X, Y, lmeformule=" ~ -1 + poly(time,3) + (1|id)", alpha=0.99, lambda=0.0001, eps = 0.001, maxiter=1000)
date()   # 22 seconds for the HMP data
c(res4$conv, res4$iter)   # 362 iterations for the HMP-data

# plot the weights for the X- and Y-variables and the correlations of the X- and Y-variables with the latent variables
#       and the estimated means of the latent variables
makeplots(X, Y, res4, lmeformule=" ~ poly(time,3) + (1|id)", weigths=TRUE, loadings=TRUE, change=TRUE, scale = T); dev.new()

# plot the observed and predicted values of a (X-/Y-)variable of a bunch of k (max. 16) randomly selected patients
which(abs(res4$a) > 0.09)  # 22 for HMP_data
plotsofpatientsof1variable(X, set="X", number=22, res4, k=16, lmeformule=" ~ poly(time,3) + (1|id)"); dev.new()     # scaling van pred_thetaas?

# plots of selected variables
plotsofselectedvariables(X, set="X", variables=which(abs(res4$a) > 0.09), res4, lmeformule=" ~ poly(time,3) + (1|id)"); dev.new()


source("C:/Users/PC/OneDrive/github/ccalb/scripts/toscca_me.R")
res_toscca = toscca.core(alphaInit = runif(ncol(X)-2), X, Y, 100, 50, lmeformula = " ~ poly(time,3) + (1|id)")
c(res_toscca$conv, res_toscca$iter)
makeplots(X, Y, res_toscca, lmeformule=" ~ poly(time,3) + (1|id)", weigths=TRUE, loadings=TRUE, change=TRUE, scale = T); dev.new()
plotsofpatientsof1variable(X, set="X", number=1, res_toscca, k=16, lmeformule=" ~ poly(time,3) + (1|id)"); dev.new()
plotsofselectedvariables(X, set="X", variables=1:20, res_toscca, lmeformule=" ~ poly(time,3) + (1|id)"); dev.new()
plotsofselectedvariables(Y, set="Y", variables=1:10, res_toscca, lmeformule=" ~ poly(time,3) + (1|id)")


#########################################################################################################################################################################################################

plotsofselectedvariables=function(Z, set="X", variables=1:16, resx, lmeformule=NA) {
   library(lme4)
   rangorde=order(Z[,"id"], Z[,"time"])
   Z=Z[rangorde,]
   id=Z[,"id"]
   time=Z[,"time"]
   a = resx$a
   setx=set
   if (set != "X") {
      a=resx$b
      setx="Y"
   }
   if (is.na(lmeformule)) {lmeformule=" ~ time + (1|id)"}
   k = length(variables)
   if (k > 16) {
      k=16
      variables=sort(sample(variables,k,replace=FALSE))
   }
   ff1=c(4,4) * (k > 9) + c(3,3) * (k <= 9) * (k > 4) + c(2,2) * (k <= 4) * (k > 1) + c(1,1) * (k==1)
   par(mfrow=ff1)
   library(lme4)
   i = variables[1]
   for (i in variables) {
      zvar=Z[,(i+2)]
      ff2 = lmer(as.formula(paste("zvar",lmeformule)))
      ff3 = predict(ff2, newdata=data.frame(time=sort(unique(time)), id = -1), allow.new.levels=TRUE, re.form=NA)   # re.form=NA because we use the population-means here
      plot(time, zvar, xlab="time", ylab=paste(setx,"-variable:",i,sep=""), ylim=c(min(zvar,ff3),max(zvar,ff3)))
      lines(sort(unique(time)), ff3, lwd=2, col=2)
      theta = as.matrix(Z[,-c(1,2)]) %*% a
      ff4 = lmer(as.formula(paste("theta",lmeformule)))
      predtheta = predict(ff4)
      ff5 = predict(ff4, newdata=data.frame(time=sort(unique(time)), id = -1), allow.new.levels=TRUE, re.form=NA)   # re.form=NA because we use the population-means here
      lines(sort(unique(time)), predict(lm(zvar ~ predtheta),newdata=data.frame(predtheta=ff5)), col=3, lwd=2)
   }
 }

#########################################################################################################################################################################################################

plotsofpatientsof1variable = function(Z, set="X", number, resx, k=16, lmeformule=NA) {
   library(lme4)
   rangorde=order(Z[,"id"], Z[,"time"])
   Z=Z[rangorde,]
   id=Z[,"id"]
   time=Z[,"time"]
   a = resx$a
   setx=set
   if (set != "X") {
      a=resx$b
      setx="Y"
   }
   if (is.na(lmeformule)) {lmeformule=" ~ time + (1|id)"}
   theta=as.matrix(Z[,3:(ncol(Z))]) %*% a
   gg1=lmer(as.formula(paste("theta", lmeformule)))
   predLV=predict(gg1)#, newdata=data.frame(time=time, id=id), re.form=NULL)
   zvar=Z[,(number+2)]
   ff1=lmer(as.formula(paste("zvar", lmeformule)))
   predictzvar=predict(ff1)
   predzvar_LV=predict(lm(zvar ~ predLV))
   patnrs=sort(unique(id))
   i=1
   selectedpatnrs=sort(sample(patnrs, k, replace=FALSE))
   ff1=c(4,4) * (k > 9) + c(3,3) * (k <= 9) * (k > 4) + c(2,2) * (k <= 4) * (k > 1) + c(1,1) * (k==1)
   par(mfrow=ff1)
   i=selectedpatnrs[1]
   for (i in selectedpatnrs) {
      og=min(zvar[id==i],predictzvar[id==i], predzvar_LV[id==i])
      bg=max(zvar[id==i],predictzvar[id==i], predzvar_LV[id==i])
      plot(time[id==i], zvar[id==i], type="b", xlab="time", ylab="variable value", ylim=c(og,bg),main=paste("patient:",i))
      lines(time[id==i], predictzvar[id==i], type="l", col=2)
      lines(time[id==i], predzvar_LV[id==i], type="l", col=3)
   }
   title(paste(setx, "-variable: ", number,sep=""), outer=TRUE, line=-1)
}

#########################################################################################################################################################################################################

makeplots = function(X, Y, resx, lmeformule=NA, weigths=TRUE, loadings=FALSE, change=FALSE, scale = FALSE) {

   if (is.na(lmeformule)) {lmeformule=" ~ time + (1|id)"}

   library(lme4)
   if(weigths) {
      dev.new()
      par(mfrow=c(1,2))
      plot(1:length(resx$a), resx$a, xlab="variable number", ylab="weights of X-variables")
      abline(h=0)
      plot(1:length(resx$b), resx$b, xlab="variable number", ylab="weights of Y-variables")
      abline(h=0)
   }

   theta = (as.matrix(X[,-c(1,2)]) %*% resx$a); if(scale == TRUE) {theta = scalar1(theta)}
   xi = (as.matrix(Y[,-c(1,2)]) %*% resx$b); if(scale == TRUE) {xi = scalar1(xi)}

  if(loadings) {
      dev.new()
      par(mfrow=c(1,2))
      plot(1:length(resx$a) , cor(X[,-c(1,2)], as.numeric(theta)), xlab="variable number", ylab="loadings of X-variables on theta")
      abline(h=0)
      plot(1:length(resx$b) , cor(Y[,-c(1,2)], as.numeric(xi))   , xlab="variable number", ylab="loadings of Y-variables on xi")
      abline(h=0)
   }

   if(change) {
      dev.new()
      rangorde = order(X[,"time"], X[,"id"])
      X=X[rangorde,]
      id = X[,"id"]
      time = X[,"time"]
      ff1 = aggregate(theta, by=list(time), FUN=mean)
      ff3 = lmer(as.formula(paste("theta", lmeformule)))
      ff4 = predict(ff3, newdata=data.frame(time=sort(unique(time)), id=-1), re.form=NA, allow.new.levels=TRUE)
      par(mfrow=c(1,2))
      ff1 = aggregate(theta, by=list(time), FUN=mean)
      ff3 = lmer(as.formula(paste("theta", lmeformule)), data = data.frame(X))
      ff4 = predict(ff3, newdata=data.frame(time=sort(unique(time)), id=-1), re.form=NA, allow.new.levels=TRUE)
      plot(ff1[,1], ff1[,2], type="b", xlab="time-point", ylab="mean of the latent variable of the X-variables", ylim=c(min(ff1[,2],ff4),max(ff1[,2],ff4)))
      lines(sort(unique(time)), ff4, col=2, lwd=2)
      rangorde = order(Y[,"time"], Y[,"id"])
      Y=Y[rangorde,]
      id = Y[,"id"]
      time = Y[,"time"]
      ff2 = aggregate(xi, by=list(time), FUN=mean)
      ff5 = lmer(as.formula(paste("xi", lmeformule)), data = data.frame(Y))
      ff6 = predict(ff5, newdata=data.frame(time=sort(unique(time)), id=-1), re.form=NA, allow.new.levels=TRUE)
      plot(ff2[,1], ff2[,2], type="b", xlab="time-point", ylab="mean of the latent variable of the Y-variables", ylim=c(min(ff2[,2],ff6),max(ff2[,2],ff6)))
      lines(sort(unique(time)), ff6, col=3, lwd=2)
      title("estimated change of the mean of the latent variables with time", outer=T, line=-1)
   }
}

#########################################################################################################################################################################################################

estimate4_ab = function(X, Y, lmeformule=NA, alpha=1, lambda=0.0005, eps = 0.001, maxiter=100) {
   # criss-cross algorithm without assuming X and Y are matrices with the same numer of rows

   if (is.na(lmeformule)) {lmeformule=" ~ -1 + time + (1|id)"}

   library(lme4)
   library(glmnet)

   id_x = X[,"id"]
   id_y = Y[,"id"]
   time_x = X[,"time"]
   time_y = Y[,"time"]

   # step 0: initialize weight-vector a
   a = rep(1, (ncol(X)-2))/(ncol(X)-2)
   b = rep(0, (ncol(Y)-2))/(ncol(Y)-2)

   iter = 0
   conv = 1
   while (conv > eps & iter <= maxiter) {

      iter = iter + 1

      # step 1: calculate theta's using current value of the a-vector
      theta = as.matrix(X[,-c(1,2)]) %*% a
      ff1 = data.frame(theta=theta, time=time_x, id=id_x)
      ff2 = lmer(as.formula(paste("theta", lmeformule)), data=ff1, REML=TRUE)
      ff3 = predict(ff2, newdata=data.frame(time=time_y, id=id_y), allow.new.levels=TRUE, re.form=NULL)

      # step 2: regress predicted-theta on Y and obtain the b-vector
      bnew=glmnet(y = ff3, x = Y[,-c(1,2)], family="gaussian", alpha=alpha, lambda=lambda, intercept=FALSE)$beta[,1]
      rm(ff1, ff2, ff3)

      # step 3: calculate xi's using current value of the b-vector and stack and scale xi's
      xi = as.matrix(Y[, -c(1,2)]) %*% bnew
      ff1 = data.frame(xi=xi, time=time_y, id=id_y)
      ff2 = lmer(as.formula(paste("xi", lmeformule)), data=ff1, REML=TRUE)
      ff3 = predict(ff2, newdata=data.frame(time=time_x, id=id_x), allow.new.levels=TRUE, re.form=NULL)

      # step 4: re-calculate the a-vector a by regression of stacked_xi on stacked_X
      anew = glmnet(y = ff3, x = X[,-c(1,2)], family="gaussian", alpha=alpha, lambda=lambda, intercept=FALSE)$beta[,1]
      rm(ff1, ff2, ff3)

      # step 5: check convergence
      conv = max(c(abs(a-anew), abs(b-bnew)))
      a = anew
      b = bnew
   }
   return(list(a=anew, b=bnew, conv=conv, iter=iter))
}

#########################################################################################################################################################################################################

descriptives = function(X, Y, plotit=TRUE, k=16, lmeformule=NA) {
   ggmin1="                    "
   gg0=data.frame(kol1=" ", kol1a=ggmin1, kol2=" ")
   gg1=data.frame(kol1=paste("number of records: ", nrow(X),sep=""), kol1a=ggmin1, kol2=paste("number of records: ", nrow(Y), sep=""))
   id_x = unique(X[,1])
   id_y = unique(Y[,1])
   nxy=which(id_x %in% id_y)
   nyx=which(id_x %in% id_y)
   gg2=data.frame(kol1=paste("subjects: ", length(id_x),sep=""), kol1a=ggmin1, kol2=paste("subjects: ", length(id_y),sep=""))
   gg2a=aggregate(rep(1,nrow(X)), by=list(X[,1]), FUN=sum)[,2]
   gg2b=aggregate(rep(1,nrow(Y)), by=list(Y[,1]), FUN=sum)[,2]
   gg2c=data.frame(kol1=paste("mean (SD) number of repeated measures/subject: ",round(mean(gg2a),3), " (",round(sd(gg2a),3),")", sep=""), kol1a=ggmin1,
         kol2=paste("mean (SD) number of repeated measures/subject: ",round(mean(gg2b),3), " (",round(sd(gg2b),3),")", sep=""))
   gg3=data.frame(kol1=paste("number of X-subjects in Y: ", length(nxy), sep=""), kol1a=ggmin1, kol2=paste("number of Y-subjects in X: ", length(nyx),sep=""))

   p = ncol(X) - 2
   q = ncol(Y) - 2
   gg4=data.frame(kol1=paste("numer of x-variables: ", p,sep=""), kol1a="     ", kol2=paste("number of y-variables: ", q, sep=""))

   time_x = X[,2]
   time_y = X[,2]
   kx=unique(time_x)
   ky=unique(time_y)
   gg5=data.frame(kol1=paste("number of time-points: ", length(kx),sep=""), kol1a=ggmin1, kol2=paste("number of time-points: ", length(ky),sep=""))
   gg6=data.frame(kol1=paste("mean of time-points: ", round(mean(X[,2]),3), sep=""), kol1a=ggmin1, kol2=paste("mean of time-points: ", round(mean(Y[,2]),3),sep=""))
   gg7=data.frame(kol1=paste("SD of time-points: ", round(sd(X[,2]),3),sep=""), kol1a=ggmin1, kol2=paste("SD of time-points: ", round(sd(Y[,2]),3),sep=""))

   gg10=rbind(gg0, gg1, gg2,gg2c, gg3, gg0,gg4,gg0,gg5,gg6,gg7,gg0)
   names(gg10)[c(1,2,3)]=c("X-variables", "" , "Y-variables")
   row.names(gg10)=1:nrow(gg10)
   print(noquote(gg10))
   if (length(kx) <= 10 & length(ky) <= 10) {
      ff1=table(X[,2])
      ff2=table(Y[,2])
      ff1=data.frame(tp=as.numeric(names(ff1)),aantalx=as.numeric(ff1))
      ff2=data.frame(tp=as.numeric(names(ff2)),aantaly=as.numeric(ff2))
      names(ff1)=c("time-point","occurring in X")
      names(ff2)=c("time-point","occurring in Y")
      ff3=merge(ff1,ff2)
      print(ff3)
   }
   if (length(kx) > 10 | length(ky) > 10) {
      dev.new()
      par(mfrow=c(1,2))
      hist(X[,2], xlab="time-points", main="x-variables")
      hist(Y[,2], xlab="time-points", main="y-variables")
   }

   if (plotit) {
      if(is.na(lmeformule)) {lmeformule="~ time + (1|id)"}
      if (k > 16) {k=16}
      ff1=c(4,4) * (k > 9) + c(3,3) * (k <= 9) * (k > 4) + c(2,2) * (k <= 4) * (k > 1) + c(1,1) * (k==1)
      dev.new()
      par(mfrow=ff1)
      library(lme4)
      selectie_x = sort(sample(1:p,k,replace=FALSE))
      i = selectie_x[1]
      for (i in selectie_x) {
         id=X[,1]
         time =X[,2]
         zvar=X[,(i+2)]
         ff2 = lmer(as.formula(paste("zvar",lmeformule)))
         ff3 = predict(ff2, newdata=data.frame(time=sort(unique(c(time_x, time_y))), id = -1), allow.new.levels=TRUE, re.form=NA)
         plot(time, zvar, xlab="time", ylab=paste("x-variable:",i), ylim=c(min(zvar,ff3),max(zvar,ff3)), xlim=c(min(sort(unique(c(time_x, time_y)))), max(sort(unique(c(time_x, time_y))))))
         lines(sort(unique(c(time_x, time_y))), ff3, lwd=2, col=3)
      }
      dev.new()
      par(mfrow=ff1)
      selectie_y = sort(sample(1:q,k,replace=FALSE))
      i = selectie_y[1]
      for (i in selectie_y) {
         id=Y[,1]
         time =Y[,2]
         zvar=Y[,(i+2)]
         ff2 = lmer(as.formula(paste("zvar",lmeformule)))
         ff3 = predict(ff2, newdata=data.frame(time=sort(unique(c(time_x, time_y))), id = -1), allow.new.levels=TRUE, re.form=NA)
         plot(time, zvar, xlab="time", ylab=paste("y-variable:",i), ylim=c(min(zvar,ff3),max(zvar,ff3)), xlim=c(min(sort(unique(c(time_x, time_y)))), max(sort(unique(c(time_x, time_y))))))
         lines(sort(unique(c(time_x, time_y))), ff3, lwd=2, col=3)
      }
   }
}

#########################################################################################################################################################################################################

estimate0_ab = function(X, Y, alpha=1, lambda=0.0005, eps = 0.001, maxiter=100) {
   # criss-cross algorithm assuming that X and Y are matrices with the same numbers of rows
   library(glmnet)
   a = rep(1, ncol(X))/ncol(X)
   b = rep(0, ncol(Y))/ncol(Y)
   iter = 0
   conv = 1
   while (conv > eps & iter <= maxiter) {
      iter = iter + 1
      # step 1: calculate theta's using current value of the a-vector
      theta = as.matrix(X) %*% a
       # step 2: regress theta on Y to obtain the b-vector
      bnew = glmnet(y = (theta), x = Y, family="gaussian", alpha=alpha, lambda=lambda, intercept=FALSE)$beta[,1]
      # step 3: calculate xi's using current value of the b-vector and stack and scale xi's
      xi = as.matrix(Y) %*% bnew
      # step 4: re-calculate the a-vector a by regression of xi on X
      anew = glmnet(y = (xi), x = X, family="gaussian", alpha=alpha, lambda=lambda, intercept=FALSE)$beta[,1]
      # step 5: check convergence
      conv = max(c(abs(a-anew), abs(b-bnew)))
      a = anew
      b = bnew
   }
   return(list(a=anew, b=bnew, conv=conv, iter=iter))
}

#########################################################################################################################################################################################################

do_separate_lmes=function(Z, lmeformule=NA) {
   library(lme4)
   if (is.na(lmeformule)) {lmeformule=" ~ 1 + time + (1|id)"}
   id=Z[,1]
   time=Z[,2]
   j=1
   zvar=Z[,(j+2)]
   ff1=lmer(as.formula(paste("zvar",lmeformule)))
   ff2=ranef(ff1)$id
   randomeffects=as.data.frame(matrix(NA, nrow=nrow(ff2), ncol=ncol(ff2)*(ncol(Z)-2)))
   randomeffects[,1:ncol(ff2)]=ff2
   j=2
   for (j in 2:(ncol(Z)-2)) {
      zvar=Z[,(j+2)]
      ff1=lmer(as.formula(paste("zvar",lmeformule)))
      ff2=ranef(ff1)$id
      randomeffects[,((ncol(ff2)*(j-1)+1) : (ncol(ff2)*(j-1)+ncol(ff2)) )]=ff2
      if ((j||100)==0) {print(paste("variable:",j))}
   }
   return(randomeffects)
}

#########################################################################################################################################################################################################

makeplots0 = function(X, Y, resx, weigths=TRUE, loadings=FALSE) {
   if(weigths) {
      dev.new()
      par(mfrow=c(1,2))
      plot(1:length(resx$a), resx$a, xlab="variable number", ylab="weights of X-variables")
      abline(h=0)
      plot(1:length(resx$b), resx$b, xlab="variable number", ylab="weights of Y-variables")
      abline(h=0)
   }
  if(loadings) {
      theta = as.matrix(X) %*% resx$a
      xi = as.matrix(Y) %*% resx$b
      dev.new()
      par(mfrow=c(1,2))
      plot(1:length(resx$a) , cor(X, as.numeric(theta)), xlab="variable number", ylab="loadings of X-variables on theta")
      abline(h=0)
      plot(1:length(resx$b) , cor(Y, as.numeric(xi))   , xlab="variable number", ylab="loadings of Y-variables on xi")
      abline(h=0)
   }
}


