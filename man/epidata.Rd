\name{epidata}
\alias{epidata}

\title{
Simulates epidemic for the specifiedmodel type and parameters 
}
\description{
This function allows the user to simulate epidemics under different models and scenarios
}
\usage{
epidata (type, n, tmin = NULL, tmax, alpha, beta, spark = NULL, Sformula = NULL, 

         x = NULL, y = NULL, inftime = NULL, infperiod = NULL, contact = NULL, 
 
         tempseed = NULL)
}

\arguments{
  \item{type}{
 Type of compartment framework, with the choice of "SI" for
Susceptible-Infectious diseases and "SIR" for
Susceptible-Infectious-Removed
}
  \item{n}{
Population size
}
  \item{tmin}{
The time point at which simulation begins, default value is one
}
  \item{tmax}{
The last time point of simulation
}
  \item{alpha}{
Susceptibility parameter (>0)
}
  \item{beta}{
Spatial parameter(s) (>0) or network parameter (s) (>0) if contact is used}

\item{spark}{
Sparks parameter (>=0), representing infections unexplained by other partsof the model (eg. infections coming in from outside the observed population), default value is zero
}
  
  \item{Sformula}{
An object of class formula. See \link{formula}

Individual-level covariate information associated with susceptibility can be passed through this argument. An expression of the form \code{ ~ model} is interpreted as a specification that the susceptibility function, \eqn{\Omega_s(i) }  is modelled by a linear predictor specified symbolically by the model term. Such a model consists of a series of terms separated by \code{+} and \code{-} operators. If there is no covariate information, \code{Sformula} is null

}
  \item{x}{
 X coordinates of individuals
}
  \item{y}{
Y coordinates of individuals
}
  \item{inftime}{
Times at which individuals are infected to initialize epidemic simulation
}
  \item{infperiod}{
Length of infectious period for each individual 
}
  \item{contact}{
Contact network matrix (matrices)
}
 \item{tempseed}{
Integer seed value to initialize the (Fortran) random number generator, default value is  a random seed. 
}


}


\details{
We consider following two individual level models: 

{\bold{Spatial model:}}

\deqn{P(i,t) =1- \exp\{-\Omega_s(i) \sum_{j \in I(t)}{d_{ij}^{-\beta}- \varepsilon}\}}

{\bold{Network model:}} 

\deqn{P(i,t) =1- \exp\{-\Omega_s(i) \sum_{j \in I(t)}{(\beta_1 C^{(1)}_{ij}} + \dots + \beta_n  C^{(n)}_{ij} )- \varepsilon\}}

where \eqn{P(i,t)} is the probability that susceptible individual i is infected at time point t, becoming infectious at time t+1; and  \eqn{\Omega_s(i) } is  a susceptibility function which accommodates potential risk factors associated with susceptible individual i contracting the disease.


}
\value{

\item{inftime }{Times at which individuals become infected/infectious}

\item{removaltime}{Times at which individuals are removed, only when \code{type = "SIR"} is specified.} 
}

\examples{


## Example 1:  spatial SI model 
# generate 100 individuals

x <- runif(100, 0, 10)

y <- runif(100, 0, 10)

covariate <- runif(100, 0, 2)

out <- epidata(type = "SI",n = 100, Sformula = ~covariate, tmax = 15, 
               alpha = c(0.1, 0.3), beta = 5.0, x = x, y = y)

# Plots of epidemic progression (optional)

epispatial(type = "SI", x = x, y = y, inftime = out$inftime)

epicurve(type = "SI", inftime = out$inftime, plottype = "newinfect")

## Example 2: spatial SIR model 
# generate infectious period(=3) for 100 individuals 

lambda <- rep(3, 100)

epidata(type = "SIR", n = 100, tmax = 15, alpha = 0.3, beta = 5.0, infperiod = lambda, 
        x = x, y = y)

## Example 3:   SI network model 
\dontrun{
contact1 <- matrix(rbinom(10000, 1, 0.1), nrow = 100, ncol = 100)

contact2 <- matrix(rbinom(10000, 1, 0.1), nrow = 100, ncol = 100)

diag(contact1[,] ) <- 0

diag(contact2[,] ) <- 0

contact <- array(c(contact1, contact2), dim = c(100, 100, 2))

epidata(type = "SI", n = 100, tmax = 15, alpha = 0.3, beta = c(3.0, 5.0), 
        contact = contact)

}
}
\references{

Deardon R, Brooks, S. P., Grenfell, B. T., Keeling, M. J., Tildesley, M. J., Savill, N. J., Shaw, D. J.,  Woolhouse, M. E. (2010). Inference for individual level models of infectious diseases in large populations. \emph{Statistica Sinica}, 20, 239-261.

Rob Deardon, Xuan Fang, and Grace P. S. Kwong (2014).
 Statistical modelling of spatio-temporal infectious disease transmission in analyzing and modeling Spatial and temporal dynamics of infectious diseases,
 \emph{(Ed: D. Chen, B. Moulin, J. Wu), John Wiley & Sons.} Chapter 11.
   
   
}
