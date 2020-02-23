\name{plot.epimcmc}
\alias{plot.epimcmc}
\title{Plot the output of \code{epimcmc} object}
\description{
  \code{plot.epimcmc} is an S3 method that plots the output of an S3 object of class \code{epimcmc}. 
}
\usage{
\method{plot}{epimcmc}(x, partype, start = 1, end = NULL, thin = 1, ...)
}
\arguments{
  \item{x}{An S3 object of class \code{epimcmc} (i.e. the output of the \code{epimcmc} function).}
  \item{partype}{Determines which of two options to plot the output of the \code{epimcmc} function are used: ``parameter'' produces trace plots for each of the model parameters, and ``loglik'' produces trace plot of the log-likelihood of the MCMC samples.}
  \item{start, end, thin}{options for creating \code{\link{mcmc}} object.}
  \item{...}{additional arguments that are passed to the generic \code{plot} function.}
}
\value{
plot.
}
\seealso{
  \code{\link{epimcmc}}, \code{\link{summary.epimcmc}}, \code{\link{mcmc}}, \code{\link{plot.mcmc}}.
}
\keyword{plot}

\examples{


## Example :  spatial SI model
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

# plot estimates
plot(epi, partype = "parameter")
}