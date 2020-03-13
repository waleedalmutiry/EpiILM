\name{EpiILM-package}
\alias{EpiILM-package}
\alias{EpiILM}
\docType{package}

\title{\pkg{EpiILM}: \packageTitle{EpiILM}}

\description{
The \R package \pkg{EpiILM} is provided for simulating from, and carrying out Bayesian MCMC-based statistical inference for spatial and/or network-based individual-level modelling framework. The package allows for the incorporation of individual-level susceptibility and transmissibility covariates in models, and provides various methods of summarizing epidemic data sets.
}

\details{
The \R package \pkg{EpiILM} can be used to carry out simulation of epidemics, estimate the basic reproduction number, plot various epidemic summary graphics, calculate the log-likelihood, carry out Bayesian inference using Metropolis-Hastings MCMC, and implement posterior predictive checks and model selection for a given data set and model. The key functions for this package are detailed in the value section. One of the important functions \code{\link{epimcmc}} depends heavily on the \code{\link[adaptMCMC]{MCMC}} from the \pkg{adaptMCMC} package for performing the MCMC analysis. This function implements the robust adaptive Metropolis sampler of Vihola (2012) for tuning the covariance matrix of the (normal) jump distribution adaptively to achieve the desired acceptance rate. The package has other features for making predictions or forecasting for a specific model via the \code{pred.epi} function. The main functions, including for epidemic simulation (\code{\link{epidata}}) and likelihood calculation (\code{\link{epilike}}) are coded in Fortran in order to achieve the goal of agile implementation.
}
\value{
	Key functions for this package:
	\item{\code{\link{epidata}}}{Simulates epidemics for the specified model type and parameters.}
	\item{\code{\link{epilike}}}{Calculates the log-likelihood for the specified model and data set.}
	\item{\code{\link{epimcmc}}}{Runs an MCMC algorithm for the estimation of specified model parameters.}
	\item{\code{\link{pred.epi}}}{Computes posterior predictions for a specified epidemic model.}
}

\author{
Vineetha Warriyar. K. V., Waleed Almutiry, and Rob Deardon \cr Maintainer: \packageMaintainer{EpiILM}
}

\examples{
\dontrun{
demo(EpiILM.spatial)
demo(EpiILM.network)
}
}
\references{
Deardon, R., Brooks, S. P., Grenfell, B. T., Keeling, M. J., Tildesley, M. J., Savill, N. J., Shaw, D. J., and Woolhouse, M. E. (2010). Inference for individual level models of infectious diseases in large populations. \emph{Statistica Sinica}, 20, 239-261.

Vihola, M. (2012) Robust adaptive Metropolis algorithm with coerced acceptance rate. \emph{Statistics and Computing}, 22(5), 997-1008. doi:10.1007/s11222-011-9269-5.
}

\keyword{package}
