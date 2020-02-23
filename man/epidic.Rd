\name{epidic}
\alias{epidic}

\title{
 Deviance Information Criterion (DIC)
}
\description{
Computes the Deviance Information Criterion for individual level models
}
\usage{
epidic (burnin, niter, LLchain, LLpostmean)
}

\arguments{
  \item{burnin}{
Burnin period for MCMC
}
  \item{niter}{
Number of MCMC iterations
}
  \item{LLchain}{
Loglikelihood values from the MCMC output 
}
  \item{LLpostmean}{
Loglikelihood value of the model with posterior mean of estimates
}
}

\examples{

## Example 1:  spatial SI model
# generate 100 individuals

x <- runif(100, 0, 10)

y <- runif(100, 0, 10)

covariate <- runif(100, 0, 2)

out1 <- epidata(type = "SI", n = 100, Sformula = ~covariate, tmax = 15,
              sus.par = c(0.1, 0.3), beta = 5.0, x = x, y = y)

unif_range <- matrix(c(0, 0, 10000, 10000), nrow = 2, ncol = 2)

# estimate parameters
mcmcout <- epimcmc(out1, tmax = 15, niter = 1500,
                    Sformula = ~covariate,   
                    sus.par.ini = c(0.003, 0.01), beta.ini =0.01,
                    pro.sus.var = c(0.1, 0.1),pro.beta.var = 0.5,
                    prior.sus.par = unif_range,
                    prior.sus.dist = c("uniform","uniform"), prior.beta.dist = "uniform",
                    prior.beta.par = c(0, 10000), adapt = TRUE, acc.rate = 0.5 )
# store the estimates
sus.parameters = c(mean(unlist(mcmcout$Estimates[1])), mean(unlist(mcmcout$Estimates[2])))
beta.par = mean(unlist(mcmcout$Estimates[3]))

# likelihood value
loglike <- epilike(out1, tmax = 15, Sformula = ~covariate, sus.par = sus.parameters, 
                   beta = beta.par)

# deviance information criterion calculation for the above epidemic
dic <- epidic(burnin = 500, niter = 1500, LLchain = mcmcout$Loglikelihood,
              LLpostmean = loglike)
dic

}

\references{
Spiegelhalter, D., Best, N., Carlin, B., Van der Linde, A. (2002). Bayesian Measures of Model Complexity and Fit. \emph{Journal of the Royal Statistical Society. Series B (Statistical Methodology)}, 64(4), 583-639.
}

