### R code from vignette source 'Predict.Rnw'

###################################################
### code chunk number 1: Predict.Rnw:38-42
###################################################
options(continue="  ", width=60)
options(SweaveHooks=list(fig=function() par(mar=c(3.1, 4.1, 3.1, 1.1))))
pdf.options(pointsize=8) #text in graph about the same as regular text
library(EpiILM, quietly=TRUE)


###################################################
### code chunk number 2: Predict.Rnw:45-52
###################################################
set.seed(101)
n <- 100
contact <- matrix(0, n, n)
for(i in 1:(n-1)) {
   contact[i, ((i+1):n)] <- rbinom((n-i), 1, 0.05)
   contact[((i+1):n), i] <- contact[i, ((i+1):n)]
}


###################################################
### code chunk number 3: Predict.Rnw:63-66
###################################################
set.seed(101)
netdat <- epidata( type = "SI", n = 100, tmax = 25,  sus.par = 0.1, 
                   contact = contact)


###################################################
### code chunk number 4: Predict.Rnw:71-75
###################################################
t_end <- max(netdat$inftime)
mcmcout_net <- epimcmc(object = netdat,  tmax = t_end,
        niter = 50000, sus.par.ini = 0.01,pro.sus.var = 0.01,
        prior.sus.dist = "uniform", prior.sus.par = c(0, 10000))


###################################################
### code chunk number 5: nettrace
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(mcmcout_net, partype = "parameter", start = 10001, density = FALSE)


###################################################
### code chunk number 6: Predict.Rnw:85-92
###################################################
set.seed(1001)
pred.net.15 <- pred.epi(netdat, xx = mcmcout_net, tmin = 15, 
                       burnin = 1000, criterion = "newly infectious", 
                       n.samples = 500)
pred.net.20 <- pred.epi(netdat, xx = mcmcout_net, tmin = 20, 
                       burnin = 1000, criterion = "newly infectious", 
                       n.samples = 500)


###################################################
### code chunk number 7: t15
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(pred.net.15, col = "red", lwd = 2, pch = 19)


###################################################
### code chunk number 8: t20
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(pred.net.20, col = "red", lwd = 2, pch = 19)


###################################################
### code chunk number 9: Predict.Rnw:115-123
###################################################
# simulate spatial locations
set.seed(101)
x <- runif(100, 0, 10)
y <- runif(100, 0, 10)
# simulate covariate, number of animals on each farm
A <- round(rexp(100,1/50))
SI.cov <- epidata(type = "SI", n = 100, tmax = 25, x = x, y = y,
                 Sformula = ~A, sus.par = c(0.5, 0.5), beta = 6)


###################################################
### code chunk number 10: tM1
###################################################
getOption("SweaveHooks")[["fig"]]()
t_end <- max(SI.cov$inftime)
prior_par <- matrix(rep(1, 4), ncol = 2, nrow = 2)
mcmcout_M1 <- epimcmc(SI.cov, Sformula = ~A, tmax = t_end, niter = 50000, 
                      sus.par.ini = c(0.001, 0.001), beta.ini = 0.01, 
                      pro.sus.var = c(0.08, 0.4), pro.beta.var = 0.5, 
                      prior.sus.dist = c("gamma","gamma"), 
                      prior.sus.par = prior_par, 
                      prior.beta.dist = "uniform", prior.beta.par = c(0, 10000) )
summary(mcmcout_M1, start = 10001)
#  MCMC traceplot for the estimation of the model (2) parameters
plot(mcmcout_M1, partype = "parameter", start = 10001, density = FALSE)


###################################################
### code chunk number 11: Predict.Rnw:151-157
###################################################
set.seed(101)
mcmcout_M2 <- epimcmc(SI.cov, tmax = t_end, niter = 50000, sus.par.ini = 0.01, 
                      beta.ini = 0.01,  pro.sus.var = 0.1, pro.beta.var = 0.5, 
                      prior.sus.dist = "uniform",
                      prior.sus.par = c(0, 10000), prior.beta.dist  = "uniform",
                      prior.beta.par = c(0, 10000))


###################################################
### code chunk number 12: tM2
###################################################
getOption("SweaveHooks")[["fig"]]()
summary(mcmcout_M2,  start = 10001)  
# MCMC traceplot for the estimation of the model (3) parameters
plot(mcmcout_M2, partype = "parameter", start = 10001, density = FALSE)


###################################################
### code chunk number 13: Predict.Rnw:167-175
###################################################
set.seed(101)
pred.model1 <- pred.epi(SI.cov, Sformula = ~A, xx = mcmcout_M1, 
                    criterion = "newly infectious",  n.samples = 500)
# convert ILM model (2) into epidata object
model2.data <- as.epidata(type = "SI", n = 100, x = x, y = y, 
                          inftime = SI.cov$inftime)
pred.model2 <- pred.epi(model2.data, xx = mcmcout_M2, 
                     criterion = "newly infectious", n.samples = 500)


###################################################
### code chunk number 14: predictM1
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(pred.model1, col = "red", type = "b", lwd = 2)


###################################################
### code chunk number 15: predictM2
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(pred.model2, col = "red", type = "b", lwd = 2)


###################################################
### code chunk number 16: Predict.Rnw:193-198
###################################################
loglike <- epilike(SI.cov, tmax = t_end, Sformula = ~A, sus.par = c(0.597, 1.071),
                                         beta = 6.606)
dic1 <- epidic(burnin = 10000, niter = 50000, LLchain = mcmcout_M1$Loglikelihood,
                                         LLpostmean = loglike)
dic1


###################################################
### code chunk number 17: Predict.Rnw:202-205
###################################################
loglike <- epilike(model2.data, tmax = t_end,  sus.par = 6.210, beta = 4.942)
dic2 <- epidic(burnin = 10000, niter = 50000, LLchain = mcmcout_M2$Loglikelihood,
                                          LLpostmean = loglike)


