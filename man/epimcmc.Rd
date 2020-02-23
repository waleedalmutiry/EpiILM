\name{epimcmc}
\alias{epimcmc}

\title{
Monte Carlo Simulation
}
\description{
Runs an MCMC algorithm for the estimation of specified model parameters
}
\usage{
epimcmc (object, tmin = NULL, tmax,

         niter, sus.par.ini, trans.par.ini = NULL, beta.ini = NULL, spark.ini = NULL,

         Sformula = NULL, Tformula = NULL,

         pro.sus.var, pro.trans.var = NULL, pro.beta.var = NULL, pro.spark.var = NULL,

         prior.sus.dist, prior.trans.dist = NULL, prior.beta.dist = NULL,

         prior.spark.dist = NULL, prior.sus.par, prior.trans.par, prior.beta.par = NULL,

         prior.spark.par = NULL, adapt = FALSE, acc.rate = NULL)
}

\arguments{
\item{object}{
An object of class \code{epidata} that can be the output of \code{\link{epidata}} or \code{\link{as.epidata}}.
}
\item{tmin}{
	The first time point at which the infection occurs, default value is one.
	}
\item{tmax}{
	The last time point at which data is observed.
	}
\item{niter}{
	Number of MCMC iterations.
	}
\item{sus.par.ini}{
Initial value(s) of the susceptibility parameter(s) (>0).
}
\item{trans.par.ini}{
Initial value(s) of the transmissibility parameter(s) (>0).
}
\item{beta.ini}{
	Initial value(s) of the spatial parameter(s) (>0) or the network parameter(s) (>0) if contact network is used.
	}
\item{spark.ini}{
	Initial value of the spark parameter (>=0).
	}
\item{Sformula}{
	An object of class formula. See \link{formula}

Individual-level covariate information associated with susceptibility can be passed through this argument. An expression of the form \code{ ~ model} is interpreted as a specification that the susceptibility function, \eqn{\Omega_S(i) }  is modelled by a linear predictor specified symbolically by the model term. Such a model consists of a series of terms separated by \code{+} and \code{-} operators. If there is no susceptibility covariate information, \code{Sformula} is null.
}
\item{Tformula}{
	An object of class formula. See \link{formula}

Individual-level covariate information associated with transmissibility can be passed through this argument. An expression of the form \code{ ~ -1+model} is interpreted as a specification that the transmissibility function, \eqn{\Omega_T(j) }  is modelled by a linear predictor specified symbolically by the model terms without the incorporation of the intercept term. Such a model consists of a series of terms separated by \code{+} and \code{-} operators. If there is no transmissibility covariate information, \code{Tformula} is null.
}
\item{pro.sus.var}{
	Proposal density variance(s) for susceptibility parameter(s). If a zero value is assigned to the proposal variance of any parameter, the parameter is considered fixed to its \code{sus.par.ini} value.
}
\item{pro.trans.var}{
	Proposal density variance(s) for transmissibility parameter(s). If a zero value is assigned to the proposal variance of any parameter, the parameter is considered fixed to its \code{sus.par.ini} value.
}
\item{pro.beta.var}{
	Proposal density variance(s) for beta parameter(s). If a zero value is assigned to the proposal variance of any parameter, the parameter is considered fixed to its \code{sus.par.ini} value.
}
\item{pro.spark.var}{
	Proposal density variance for the spark parameter.
}
\item{prior.sus.dist}{
Select the prior distribution(s) for the susceptibility parameter(s) with the choice of "halfnormal" for positive half normal distribution, "gamma"  for gamma distribution and "uniform" for uniform distribution
}
\item{prior.trans.dist}{
Select the prior distribution(s) for the transmissibility parameter(s) with the choice of "halfnormal" for positive half normal distribution, "gamma"  for gamma distribution and "uniform" for uniform distribution
}
\item{prior.beta.dist}{
	Select the prior distribution(s) for the beta parameter(s) with the choice of "halfnormal" for half normal distribution, "gamma"  for gamma distribution and "uniform" for uniform distribution
	}

\item{prior.spark.dist}{
	Select the prior distribution for the spark parameter with the choice of "halfnormal" for half normal distribution, "gamma"  for gamma distribution and "uniform" for uniform distribution
	}
\item{prior.sus.par}{
A vector (matrix) of the prior distribution parameters for updating the susceptibility parameter(s).
}
\item{prior.trans.par}{
A vector (matrix) of the prior distribution parameters for updating the transmissibility parameter(s).
}

\item{prior.beta.par}{
	A vector (matrix) of the prior distribution parameters for updating the kernel parameter(s).
	}
\item{prior.spark.par}{
	A vector of the prior distribution parameters for updating the spark parameter.
	}
\item{adapt}{
	To enable the adapted MCMC method in the \code{\link{MCMC}} function, default is FALSE.
	}
\item{acc.rate}{
	To set an acceptance rate. This option will be ignored if \code{adapt} = FALSE. See \code{\link{MCMC}} for more details.
 	}
}

\value{
Returns an object of class \code{epimcmc} that contains:
   \describe{
    \item{type:}{the compartmental framework model used in the analysis.}
    \item{kernel.type:}{the used \code{kernel.type} in the function (distance-based or network-based).}
    \item{Estimates:}{the MCMC output of the updated model parameters.}
    \item{Loglikelihood:}{the loglikelihood of the updated model parameters.}
    \item{Fullsamples:}{the MCMC output of all the model parameters (including fixed parameters).}
    \item{n.sus.par:}{the number of parameters in the susceptibility function.}
    \item{n.trans.par:}{the number of parameters in the transmissibility function.}
    \item{n.ker.par:}{the number of parameters in the kernel function.}
  }
}


\details{

Independent Gaussian random walks are used as the Metropolis-Hastings MCMC proposal for all parameters. The \code{epimcmc} function depends on the \code{\link{MCMC}} function from the \pkg{\link{adaptMCMC}} package.

}
\seealso{
  \code{\link{summary.epimcmc}}, \code{\link{plot.epimcmc}}, \code{\link{epidata}}, \code{\link{epilike}}, \code{\link{pred.epi}}.
}

\examples{


## Example 1:  spatial SI model
# generate 100 individuals

x <- runif(100, 0, 10)

y <- runif(100, 0, 10)

covariate <- runif(100, 0, 2)

out1 <- epidata(type = "SI", n = 100, Sformula = ~covariate, tmax = 15,
               sus.par = c(0.1, 0.3), beta = 5.0, x = x, y = y)

alphapar1 <- matrix(c(1, 1, 1, 1), ncol = 2, nrow = 2)

betapar1 <- c(10, 2)

epi <- epimcmc(object = out1, tmin = 1, tmax = 15,
               niter = 1000, sus.par.ini = c(1, 1), beta.ini = 1,
               Sformula = ~covariate, pro.sus.var = c(0.5, 0.3), pro.beta.var = 0.1,
               prior.sus.dist = c("gamma", "gamma"), prior.beta.dist = "gamma",
               prior.sus.par = alphapar1, prior.beta.par = betapar1,
               adapt = TRUE, acc.rate = 0.5)
epi

## Example 2:  spatial SIR model

lambda <- rep(3, 100)

out2 <- epidata(type = "SIR", n = 100, tmax = 15, sus.par = 0.3, beta = 5.0, infperiod = lambda,
        x = x, y = y)

alphapar2 <- c(1, 1)
betapar2 <- c(1, 1)

epi2 <- epimcmc(object = out2, tmin = 1, tmax = 15,
               niter = 1000, sus.par.ini = 1, beta.ini = 1,
               Sformula = NULL, pro.sus.var = 0.3, pro.beta.var = 0.1,
               prior.sus.dist = "gamma", prior.beta.dist = "gamma",
               prior.sus.par = alphapar2, prior.beta.par = betapar2,
               adapt = FALSE, acc.rate = NULL)
epi2

}


\references{
Rob Deardon, Xuan Fang, and Grace P. S. Kwong (2014). Statistical modelling of spatio-temporal infectious disease tranmission in Analyzing and Modeling Spatial and Temporal Dynamics of Infectious Diseases, \emph{(Ed: D. Chen, B. Moulin, J. Wu), John Wiley & Sons.}. Chapter 11.
}
