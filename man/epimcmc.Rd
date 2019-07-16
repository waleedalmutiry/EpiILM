\name{epimcmc}
\alias{epimcmc}

\title{
Monte Carlo Simulation
}
\description{
Runs an MCMC algorithm forthe estimation of specified model parameters
}
\usage{
epimcmc (type, x = NULL, y = NULL, inftime, tmin = NULL, tmax, infperiod = NULL, 

         niter, alphaini, betaini, sparkini = NULL, Sformula = NULL, 

         contact = NULL, pro.var.a, pro.var.b, pro.var.sp = NULL, prioralpha, 

         halfnorm.var.a = NULL, gamma.par.a = NULL, unif.range.a = NULL, priorbeta, 

         halfnorm.var.b = NULL, gamma.par.b = NULL, unif.range.b = NULL, 

         priorsp = NULL, halfnorm.var.sp = NULL, gamma.par.sp = NULL, 

         unif.range.sp = NULL, tempseed = NULL)
}

\arguments{
\item{type}{
 Type of compartment framework, with the choice of "SI" for
Susceptible-Infectious diseases and "SIR" for
Susceptible-Infectious-Removed
}
  \item{x}{X coordinates of individuals}
 \item{y }{Y coordinates of individuals}
 \item{inftime }{Times at which individuals are infected to initialize epidemic simulation}
\item{tmin}{ The first time point at which the infection occurs, default value is one}
 \item{tmax}{ The last time point at which data is observed}
 \item{infperiod}{Length of infectious period for each individual }
\item{niter }{Number of MCMC iterations}
  \item{alphaini}{Initial value of susceptibility parameter(s)(>0)}
 \item{betaini}{Initial value of spatial parameter(s) (>0) or network parameter (s) (>0) if contact is used}
 \item{sparkini}{Initial value of spark parameter (>=0)}

\item{Sformula}{
An object of class formula. See \link{formula}

Individual-level covariate information associated with susceptibility can be passed through this argument. An expression of the form \code{ ~ model} is interpreted as a specification that the susceptibility function, \eqn{\Omega_s(i) }  is modelled by a linear predictor specified symbolically by the model term. Such a model consists of a series of terms separated by \code{+} and \code{-} operators. If there is no covariate information, \code{Sformula} is null

}

\item{contact}{Contact network matrix (matrices) }


\item{pro.var.a}{Proposal density variance for alpha parameter(s)}

\item{pro.var.b}{Proposal density variance for beta parameter(s) }
\item{pro.var.sp}{Proposal density variance for spark parameter }

\item{prioralpha}{

Select the prior distribution for alpha parameter(s) with the choice of "halfnormal" for positive half normal distribution, "gamma"  for gamma distribution and "uniform" for uniform distribution }

\item{halfnorm.var.a}{Half normal prior variance for alpha}

\item{gamma.par.a}{Gamma  prior: shape and rate parameters for alpha}

\item{unif.range.a}{Uniform prior: Minimum and maximum range for alpha }

\item{priorbeta}{

Select the prior distribution for beta parameter(s) with the choice of "halfnormal" for half normal distribution, "gamma"  for gamma distribution and "uniform" for uniform distribution }

\item{halfnorm.var.b}{Half normal prior variance for beta}

\item{gamma.par.b}{Gamma  prior: shape and rate parameters for beta}

\item{unif.range.b}{Uniform prior: Minimum and maximum range for beta}


\item{priorsp}{

Select the prior distribution for spark parameter with the choice of "halfnormal" for half normal distribution, "gamma"  for gamma distribution and "uniform" for uniform distribution }

\item{halfnorm.var.sp}{Half normal prior variance for spark}

\item{gamma.par.sp}{Gamma  prior: shape and rate parameters for spark}

\item{unif.range.sp}{Uniform prior: Minimum and maximum range for spark}

  \item{tempseed}{
Integer seed value to initialize the (Fortran) random number generator, default value is  a random seed. 
}

}

\value{
A list is returned with the following components:
 \item{Estimates }{MCMC output as \pkg{coda} object} 
\item{Loglikelihood }{Log likelihood value of each posterior estimate}

}


\details{

Independent Gaussian random walks are used as the Metropolis-Hastings MCMC proposal for all parameters

}

\examples{


## Example 1:  spatial SI model 
# generate 100 individuals
\dontrun{
x <- runif(100, 0, 10)

y <- runif(100, 0, 10)

covariate <- runif(100, 0, 2)

out <- epidata(type = "SI", n = 100, Sformula = ~covariate, tmax = 15, 
              alpha = c(0.1, 0.3), beta = 5.0, x = x, y = y)


mcout <- epimcmc(type = "SI", x = x, y = y, inftime = out$inftime, Sformula = ~covariate, 
                tmax = 15, niter = 1000, alphaini = c(0.01, 0.01), betaini = 0.01,
                pro.var.a = c(0.01, 0.005), pro.var.b = 0.05, 
                prioralpha = "halfnormal", halfnorm.var.a = c(10**5, 10**5), 
                priorbeta = "halfnormal", halfnorm.var.b = 10**5)
}
}


\references{
Rob Deardon, Xuan Fang, and Grace P. S. Kwong (2014).
 Statistical modelling of spatio-temporal infectious disease tranmission in Analyzing and Modeling Spatial and Temporal Dynamics of Infectious Diseases,
 \emph{(Ed: D. Chen, B. Moulin, J. Wu), John Wiley & Sons.}. Chapter 11.
   
   
}






