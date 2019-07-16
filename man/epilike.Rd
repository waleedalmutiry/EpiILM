\name{epilike}
\alias{epilike}

\title{
Calculates the log likelihood 
}
\description{
Calculates the log likelihood for thespecified individual level model and data set
}
\usage{
epilike (type, x = NULL, y = NULL, inftime, infperiod = NULL, tmin = NULL, tmax, alpha, 

         beta, spark = NULL, Sformula = NULL, contact = NULL)
}

\arguments{
  \item{type}{
 Type of compartment framework, with the choice of "SI" for
Susceptible-Infectious diseases and "SIR" for
Susceptible-Infectious-Removed
}
  \item{x}{
 X coordinates of individuals
}
  \item{y}{
 Y coordinates of individuals
}
  \item{inftime}{
Infection times 
}
  \item{infperiod}{
Length of infectious period for each individual
}
  \item{tmin}{
The first time point at which data is observed, default value is one
}
  \item{tmax}{
The last time point at which data is observed
}
  \item{alpha}{
Susceptibility parameter(>0)
}
  \item{beta}{
Spatial parameter(s) (>0) or network parameter (s) (>0) if contact is used
}
  \item{spark}{
Sparks parameter(>=0), representing infections unexplained by other partsof the model or infections coming in from outside the observed population, default value is zero
}

   \item{Sformula}{
An object of class formula. See \link{formula}

Individual-level covariate information passable through this argument. An expression of the form \code{ ~ model} is interpreted as a specification that the susceptibility function, \eqn{\Omega_s(i) }  is modelled by a linear predictor specified symbolically by the model term. Such a model consists of a series of terms separated by \code{+} and \code{-} operators. If there is no covariate information, \code{Sformula} is null

}

\item{contact}{
Contact network matrix (matrices)
}
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

epilike(type = "SI", x = x, y = y, inftime = out$inftime, tmax = 15, 
         alpha = c(0.1, 0.3), beta = 5, Sformula = ~covariate)
}
## Example 2: spatial SIR model 
# generate infectious period (=3) for 100 individuals 
\dontrun{
lambda <- rep(3, 100)

out <- epidata(type = "SIR", n = 100, tmax = 15, alpha  =0.3, beta = 5.0,
              infperiod = lambda, x = x, y = y)

epilike(type = "SIR", x = x, y = y, infperiod = lambda, inftime = out$inftime, 
        tmax = 15, alpha = 0.3, beta = 5.0)

}
}


\references{

Deardon R, Brooks, S. P., Grenfell, B. T., Keeling, M. J., Tildesley, M. J., Savill, N. J., Shaw, D. J.,  Woolhouse, M. E. (2010). Inference for individual level models of infectious diseases in large populations. \emph{Statistica Sinica}, 20, 239-261.
}
