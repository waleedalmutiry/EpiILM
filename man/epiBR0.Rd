\name{epiBR0}
\alias{epiBR0}

\title{
Basic reproduction number (R0)
}
\description{
Calculates the basic reproduction number for the specified SIR model and data set
}
\usage{
epiBR0 (x = NULL, y = NULL, contact = NULL, alpha, beta, spark = NULL, infperiod, 

        Sformula = NULL, tmax, niter, tempseed = NULL)
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
 
  \item{alpha}{
Susceptibility parameter(>0)
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

Individual-level covariate information passable through this argument. An expression of the form \code{ ~ model} is interpreted as a specification that the susceptibility function, \eqn{\Omega_s(i) }  is modelled by a linear predictor specified symbolically by the model term. Such a model consists of a series of terms separated by \code{+} and \code{-} operators. If there is no covariate information, \code{Sformula} is null.
}

  \item{tmax}{
The last time point of simulation
}  

  \item{niter}{
Number of epidemic simulations to calculate basic reproduction number 
}
  \item{tempseed}{
Integer seed value to initialize the (Fortran) random number generator, default value is  a random seed. 
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
