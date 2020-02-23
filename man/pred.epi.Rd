\name{pred.epi}
\alias{pred.epi}

\title{
Posterior predictive check.
}
\description{
Computing the posterior predictive check based on different summary statistics.
}
\usage{
pred.epi (object, xx, criterion , n.samples, burnin = NULL, tmin = NULL,

		  Sformula = NULL, Tformula = NULL, showProgressBar = interactive())
}

\arguments{
  \item{object}{
An object of class \code{epidata} that can be the output of \code{\link{epidata}} or \code{\link{as.epidata}}
}
  \item{xx}{
  		An object of class \code{epimcmc} that is the output of \code{\link{epimcmc}}.
}
  \item{criterion}{
The (multivariate) statistical criteria used in the posterior predictive check. It has three options: ``newly infectious'' which is a multivariate statistics represented by the number of newly infectious individuals over time, ``epidemic length'' represents the length of epidemic, and ``peak time'' represents the time of the peak of epidemic.
}
  \item{n.samples}{
The number of epidemics that needs to be simulated in the posterior predictive check procedure.
}
  \item{burnin}{
A scalar value which represents the number of samples needs to be discarded from the MCMC output.
}
  \item{tmin}{
	The first time point at which the infection occurs, default value is one.
}
\item{Sformula}{
An object of class formula. See \link{formula}.

Individual-level covariate information passable through this argument. An expression of the form \code{ ~ model} is interpreted as a specification that the susceptibility function, \eqn{\Omega_s(i) }  is modelled by a linear predictor specified symbolically by the model term. Such a model consists of a series of terms separated by \code{+} and \code{-} operators. If there is no covariate information, \code{Sformula} is null.
}
\item{Tformula}{
An object of class formula. See \link{formula}.

Individual-level covariate information passable through this argument. An expression of the form \code{ ~ model} is interpreted as a specification that the transmissibility function, \eqn{\Omega_t(i) }  is modelled by a linear predictor specified symbolically by the model term. Such a model consists of a series of terms separated by \code{+} and \code{-} operators. If there is no covariate information, \code{Tformula} is null.
}

\item{showProgressBar}{
logical. If TRUE a progress bar is shown.
}
}

\value{
An object of class \code{pred.epi} that contains the following:
   \describe{
    \item{type:}{The compartmental framework model used in the analysis.}
    \item{criterion:}{The (multivariate) statistical criteria used in the posterior predictive check.}
    \item{crit.sim:}{The output of the evaluated \code{criterion} on the simulated epidemics.}
    \item{crit.obs:}{The output of the evaluated \code{criterion} on the observed epidemics.}
    \item{tmax:}{The last time point at which data is observed.}
    \item{n.samples:}{The number of simulated epidemics used in the posterior predictive check procedure.}
  }

}

\seealso{
\code{\link{epimcmc}}, \code{\link{epidata}}, \code{\link{epilike}}, \code{\link{plot.pred.epi}}.
}



\examples{
## Example 1:  spatial SI model
# generate 100 individuals

x <- runif(100, 0, 10)

y <- runif(100, 0, 10)

covariate <- cbind(runif(100, 0, 2), rbinom(100, 1, 0.5))

out <- epidata(type = "SI",n = 100, Sformula = ~covariate, tmax = 15,
               sus.par = c(0.1, 0.3, 0.01), beta = 5.0, x = x, y = y)

alphapar2 <- matrix(c(1, 1, 1, 1, 1, 1), ncol = 2, nrow = 3)

betapar2 <- c(1, 1)

epi<-epimcmc(object = out, tmin = 1, tmax = 15,
	       niter = 500, sus.par.ini = c(1, 1, 1), beta.ini = 1,
	       Sformula = ~covariate,
				 pro.sus.var = c(0.5, 0.3, 0.2), pro.beta.var = 0.1,
				 prior.sus.dist = c("gamma", "gamma", "gamma"),
				 prior.beta.dist = "gamma",
				 prior.sus.par = alphapar2, prior.beta.par = betapar2,
	       adapt = TRUE, acc.rate = 0.5)

epipred1 <- pred.epi (object = out, xx = epi,
				criterion = "newly infectious",
			  n.samples = 100, burnin = 200, tmin = 1,
			  Sformula = ~covariate)

plot(epipred1, col = "red", type = "b", lwd = 2)

epipred2 <- pred.epi (object = out, xx = epi,
				criterion = "peak time",
			  n.samples = 100, burnin = 200, tmin = 1,
			  Sformula = ~covariate)

plot(epipred2, col = "dark gray")

}


\references{

Deardon R, Brooks, S. P., Grenfell, B. T., Keeling, M. J., Tildesley, M. J., Savill, N. J., Shaw, D. J.,  Woolhouse, M. E. (2010). Inference for individual level models of infectious diseases in large populations. \emph{Statistica Sinica}, 20, 239-261.
}
