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


\references{
Spiegelhalter, D., Best, N., Carlin, B., Van der Linde, A. (2002). Bayesian Measures of Model Complexity and Fit. \emph{Journal of the Royal Statistical Society. Series B (Statistical Methodology)}, 64(4), 583-639.
}

