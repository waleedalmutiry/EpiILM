# Load the EpiILM package
library("EpiILM")
set.seed(101)

# sample size
n <- 500

# Generating an undirected binary contact network:
contact <- matrix(0, n, n)
for(i in 1:(n-1)) {
	contact[i, ((i+1):n)] <- rbinom((n-i), 1, 0.05)
	contact[((i+1):n), i] <- contact[i, ((i+1):n)]
	}
	
# Generating the susceptibility and transmissibilty covariates: 
X1 <- round(rexp(n, 1/100))
X2 <- round(rgamma(n, 50, 0.5))

# Simulate epidemic form SIR network-based ILM

infp <- rep(3, n)
SIR.net <- epidata(type = "SIR", n = 500, tmax = 50, 
                                     sus.par = c(0.003, 0.002), trans.par = c(0.0003, 0.0002), 
                                     contact = contact, infperiod = infp,
                                     Sformula = ~ -1 + X1 + X2, Tformula = ~ -1 + X1 + X2)
SIR.net

# Epidemic curve for SIR.net
plot(SIR.net, plottype = "curve", curvetype = "complete")

# epimcmc function to estimate the model parameters:
t_end <- max(SIR.net$inftime)
prior_par <- matrix(rep(1, 4), ncol = 2, nrow = 2)

# This took 305.7 seconds on an Apple MacBook Pro with i5-core Intel 2.4 GHz 
# processors with 8 GB of RAM.
mcmcout_SIR.net <- epimcmc(SIR.net, tmax = t_end, niter = 20000,
                                      Sformula = ~-1 + X1 + X2, Tformula = ~-1 + X1 + X2,  
                                      sus.par.ini = c(0.003, 0.01), trans.par.ini = c(0.01, 0.01),
                                      pro.sus.var = c(0.0, 0.1), pro.trans.var = c(0.1, 0.1),
                                      prior.sus.dist = c("gamma", "gamma"), prior.trans.dist = c("gamma", "gamma"), 
                                      prior.sus.par = prior_par, prior.trans.par = prior_par,
                                      adapt = TRUE, acc.rate = 0.5)

# Summary of MCMC results
summary(mcmcout_SIR.net, start = 10001) 
plot(mcmcout_SIR.net, partype = "parameter", start = 10001, density = FALSE)

# Posterior predictive forecasting
# Posterior prediction starting from time point 5
pred.SIR.net.point.5 <- pred.epi(SIR.net, xx = mcmcout_SIR.net, tmin = 5, 
                                               Sformula = ~-1 + X1 + X2, Tformula = ~-1 + X1 + X2,  
                                               burnin = 1000, criterion = "newly infectious", 
                                               n.samples = 500)

# Posterior prediction starting from time point 14
pred.SIR.net.point.8 <- pred.epi(SIR.net, xx = mcmcout_SIR.net, tmin = 14, 
                                                  Sformula = ~-1 + X1 + X2, Tformula = ~-1 + X1 + X2,  
                                                  burnin = 1000, criterion = "newly infectious", 
                                                  n.samples = 500)

# plot predictions:
par(mfrow = c(2,1))
plot(pred.SIR.net.point.5, col = "red", lwd = 2, pch = 19)
plot(pred.SIR.net.point.8, col = "red", lwd = 2, pch = 19)
