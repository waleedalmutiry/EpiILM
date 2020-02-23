library("EpiILM")

set.seed(789)
# generating the XY coordinates of individuals:
x  <-  runif(256,  0, 100)
y  <-  runif(256,  0, 100)

# generating the sus. covariate:
A <- round(rexp(256, 1/100))
 
# simulating an epidemic:
out_cov <- epidata(type = "SI", n = 256, tmax = 10, x = x, y = y, Sformula = ~A, sus.par = c(0.01, 0.05), beta = 2)
out_cov

# plotting some characteristics of the epidmeic:
par(mfrow=c(2, 2))
plot(out_cov, plottype = "curve", curvetype = "complete")
plot(out_cov, plottype = "curve", curvetype = "susceptible")
plot(out_cov, plottype = "curve", curvetype = "totalinfect")
plot(out_cov, plottype = "curve", curvetype = "newinfect")
plot(out_cov, plottype = "spatial")

# performing the MCMC using the epimcmc function:
t_end <- max(out_cov$inftime)
unif_range <- matrix(c(0, 0, 1, 1), nrow =  2, ncol =  2)

mcmcout_M8 <- epimcmc(out_cov, Sformula = ~A,
    tmax = t_end, niter = 50000,
    sus.par.ini = c(0.03, 0.005), beta.ini = 2,
    pro.sus.var = c(0.005, 0.005), pro.beta.var = 0.01,
    prior.sus.dist = c("uniform", "uniform"), prior.sus.par = unif_range,
    prior.beta.dist = "uniform", prior.beta.par = c(0, 10),
    adapt = TRUE, acc.rate = 0.5)
            
plot(mcmcout_M8, partype = "parameter", density = FALSE )

summary(mcmcout_M8, start = 1001)


mcmcout_M9 <- epimcmc(out_cov,
    tmax = t_end, niter = 50000, sus.par.ini = 0.01,
    beta.ini = 2, pro.sus.var = 0.1, pro.beta.var = 0.5,
    prior.sus.dist = "uniform", prior.sus.par = c(0, 3),
    prior.beta.dist = "uniform", prior.beta.par = c(0, 10),
    adapt = TRUE, acc.rate = 0.5)

plot(mcmcout_M9, partype = "parameter", density = FALSE )

summary(mcmcout_M9, start = 1001)


#set.seed(23456)
predepi1<-pred.epi(object = out_cov, xx = mcmcout_M8, criterion = "newly infectious", n.samples = 500, burnin = 1001, tmin = 1, Sformula = ~A)

par(mfrow=c(1,2))
plot(predepi1, col = "red", type = "b", lwd = 2, pch = 20)



loglike1 <- epilike(object = out_cov, tmax = t_end, Sformula = ~A, sus.par = c(0.08806, 0.04421), beta = 1.96839)
loglike2 <- epilike(object = out_cov, tmax = t_end, sus.par = 0.735, beta = 1.554)

dic1 <- epidic(burnin = 10000, niter = 50000, LLchain = mcmcout_M8$Loglikelihood, LLpostmean = loglike1)
dic1

dic2 <- epidic(burnin = 10000, niter = 50000, LLchain = mcmcout_M9$Loglikelihood, LLpostmean = loglike2)
dic2

