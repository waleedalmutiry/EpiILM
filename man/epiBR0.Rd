\name{epiBR0}
\alias{epiBR0}

\title{
Basic reproduction number (R0)
}
\description{
Calculates the basic reproduction number for a specified SIR model and data set through simulation
}
\usage{
epiBR0 (x = NULL, y = NULL, contact = NULL, sus.par, trans.par = NULL, beta,

        spark = NULL, infperiod, Sformula = NULL, Tformula = NULL, tmax,

        niter)
}

\arguments{
 \item{x}{
X coordinates of individuals
}
  \item{y}{
 Y coordinates of individuals
}
  \item{contact}{
Contact network(s)
}

\item{sus.par}{
Susceptibility parameter(>0)
}
\item{trans.par}{
Transmissibility parameter(>0)
}
  \item{beta}{
Spatial parameter(s) (>0) or network parameter (s) (>0) if contact is used
}
  \item{spark}{
Sparks parameter (>=0), representing infections unexplained by other partsof the model or infections coming in from outside the observed population, default value is zero
}

  \item{infperiod}{
Length of infectious period for each individual
}

  \item{Sformula}{
An object of class formula. See \link{formula}

Individual-level covariate information associated with susceptibility can be passed through this argument. An expression of the form \code{ ~ model} is interpreted as a specification that the susceptibility function, \eqn{\Omega_S(i) }  is modelled by a linear predictor specified symbolically by the model term. Such a model consists of a series of terms separated by \code{+} and \code{-} operators. If there is no susceptibility covariate information, \code{Sformula} is null.
}
\item{Tformula}{
An object of class formula. See \link{formula}

Individual-level covariate information associated with transmissibility can be passed through this argument. An expression of the form \code{ ~ -1+model} is interpreted as a specification that the transmissibility function, \eqn{\Omega_T(j) }  is modelled by a linear predictor specified symbolically by the model terms without the incorporation of the intercept term. Such a model consists of a series of terms separated by \code{+} and \code{-} operators. If there is no transmissibility covariate information, \code{Tformula} is null.
}
  \item{tmax}{
The last time point of simulation
}

  \item{niter}{
Number of epidemic simulations to calculate basic reproduction number
}
 }

\value{
A list is returned with the following components:
\item{BasicR0}{
The basic reproduction number value
}
\item{simulated_BR0}{
Number of infections per simulation
}
}

\examples{
# generate 100 X-Y coordinates for a distance-based ILM

x <- runif(100, 0, 10)

y <- runif(100, 0, 10)

# Suppose we know the length of infectious period for each individual. Also, assume 
# susceptibility parameter = 1.5 and spatial parameter = 5 for this SIR model

infperiod <- rep(3, 100)

# For a 1000 iteration with a last observed time 
# point 15, we can calculate the basic reproduction number 

out <- epiBR0 (x = x, y = y,  sus.par = 1.5,  beta = 5, infperiod= infperiod,  
     
               tmax = 15, niter = 1000)

out$BasicR0


}