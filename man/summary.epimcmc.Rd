\name{summary.epimcmc}
\alias{summary.epimcmc}
\title{Summary method for \code{epimcmc} objects}
\description{Summarize a \code{\link{epimcmc}} object and return an object of class \code{summary.epimcmc}.
}
\usage{
\method{summary}{epimcmc}(object, \dots)
}
\arguments{
  \item{object}{an S3 object of class \code{epimcmc} (i.e. the output of the \code{epimcmc} function).}
  \item{\dots}{potential further arguments (require by generic).}
}

\seealso{
  \code{\link{epimcmc}}, \code{\link{plot.epimcmc}}.
}

\examples{


## Example:  spatial SI model
# generate 100 individuals

x <- runif(100, 0, 10)

y <- runif(100, 0, 10)

covariate <- runif(100, 0, 2)

out1 <- epidata(type = "SI", n = 100, Sformula = ~covariate, tmax = 15,
               sus.par = c(0.1, 0.3), beta = 5.0, x = x, y = y)

alphapar1 <- matrix(c(1, 1, 1, 1), ncol = 2, nrow = 2)

betapar1 <- c(1, 1)

epi <- epimcmc(object = out1, tmin = 1, tmax = 15,
               niter = 1000, sus.par.ini = c(1, 1), beta.ini = 1,
               Sformula = ~covariate, pro.sus.var = c(0.5, 0.3), pro.beta.var = 0.1,
               prior.sus.dist = c("gamma", "gamma"), prior.beta.dist = "gamma",
               prior.sus.par = alphapar1, prior.beta.par = betapar1,
               adapt = TRUE, acc.rate = 0.5)

# summary of mcmc output
summary(epi)

}