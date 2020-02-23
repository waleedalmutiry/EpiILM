\name{epidata}
\alias{epidata}

\title{
Simulates epidemic for the specifiedmodel type and parameters
}
\description{
This function allows the user to simulate epidemics under different models and scenarios
}
\usage{
epidata (type, n, tmin = NULL, tmax, sus.par, trans.par = NULL, beta = NULL, spark = NULL,

          Sformula = NULL, Tformula = NULL, x = NULL, y = NULL,

          inftime = NULL, infperiod = NULL, contact = NULL)

          
}

\arguments{
  \item{type}{
 Type of compartment framework, with the choice of "SI" for
Susceptible-Infectious diseases and "SIR" for
Susceptible-Infectious-Removed.
}
  \item{n}{
Population size
}
  \item{tmin}{
The time point at which simulation begins, default value is one.
}
  \item{tmax}{
The last time point of simulation.
}
\item{sus.par}{
Susceptibility parameter (>0).
}
\item{trans.par}{
Transmissibility parameter (>0).
}
  \item{beta}{
Spatial parameter(s) (>0) or network parameter (s) (>0) if contact network is used.
}

\item{spark}{
Sparks parameter (>=0), representing infections unexplained by other partsof the model (eg. infections coming in from outside the observed population), default value is zero.
}

\item{Sformula}{
An object of class formula. See \link{formula}.

Individual-level covariate information associated with susceptibility can be passed through this argument. An expression of the form \code{ ~ model} is interpreted as a specification that the susceptibility function, \eqn{\Omega_s(i) }  is modelled by a linear predictor specified symbolically by the model term. Such a model consists of a series of terms separated by \code{+} and \code{-} operators. If there is no susceptibility covariate information, \code{Sformula} is null.

}
\item{Tformula}{
An object of class formula. See \link{formula}.

Individual-level covariate information associated with transmissibility can be passed through this argument. An expression of the form \code{ ~ -1+model} is interpreted as a specification that the transmissibility function, \eqn{\Omega_t(j) }  is modelled by a linear predictor specified symbolically by the model terms without the incorporation of the intercept term. Such a model consists of a series of terms separated by \code{+} and \code{-} operators. If there is no transmissibility covariate information, \code{Tformula} is null.

}
  \item{x}{
 X coordinates of individuals.
}
  \item{y}{
Y coordinates of individuals.
}
  \item{inftime}{
Times at which individuals are infected to initialize epidemic simulation.
}
  \item{infperiod}{
Length of infectious period for each individual.
}
  \item{contact}{
  A contact network matrix or an array of contact network matrices.
}


}


\details{
We consider following two individual level models:

{\bold{Spatial model:}}

\deqn{P(i,t) =1- \exp\{-\Omega_s(i) \sum_{j \in I(t)}{\Omega_t(j)d_{ij}^{-\beta}- \varepsilon}\}}

{\bold{Network model:}}

\deqn{P(i,t) =1- \exp\{-\Omega_s(i) \sum_{j \in I(t)}{\Omega_t(j)(\beta_1 C^{(1)}_{ij}} + \dots + \beta_n  C^{(n)}_{ij} )- \varepsilon\}}

where \eqn{P(i,t)} is the probability that susceptible individual i is infected at time point t, becoming infectious at time t+1; 
\eqn{\Omega_s(i)} is  a susceptibility function which accommodates potential risk factors associated with susceptible individual i contracting the disease; \eqn{\Omega_t(j)} is  a transmissibility function which accommodates potential risk factors associated with infectious individual j; \eqn{\varepsilon} is a sparks term which represents infections originating from outside the population being observed or some other unobserved infection mechanism. 

The susceptibility function can incorporate any individual-level covariates  of interest and \eqn{\Omega_s(i)} is treated as a linear function of the covariates, i.e., \eqn{\Omega_s(i) = \alpha_0 + \alpha_1 X_1(i) + \alpha_2 X_2 (i) + \dots  + 
\alpha_{n_s} X_{n_s} (i)}, where \eqn{X_1(i), \dots, X_{n_s} (i)} denote  \eqn{n_s}covariates associated with susceptible individual $i$, along with susceptibility parameters \eqn{\alpha_0,\dots,\alpha_{n_s} >0}. If the model does not contain any susceptibility covariates then  \eqn{\Omega_s(i) = \alpha_0} is used. In a similar way, the transmissibility function can incorporate any individual-level covariates of interest associated with infectious individual. \eqn{\Omega_t(j)} is also treated as a linear function of the covariates, but without the intercept term, i.e., \eqn{\Omega_t(j) = \phi_1 X_1(j) + \phi_2 X_2 (j) + \dots  + \phi_{n_t} X_{n_t} (j)}, where \eqn{X_1(j), \dots, X_{n_t} (j)} denote  the \eqn{n_t} covariates associated with  infectious individual j, along with transmissibility parameters \eqn{\phi_1,\dots,\phi_{n_t} >0}. If the model does not contain any transmissibility covariates then  \eqn{\Omega_t(j)  = 1} is used.


}
\value{
An object of class \code{epidata} is returned containing the following:
\describe{
\item{type}{Type of compartment framework, with the choice of "SI" for
Susceptible-Infectious diseases and "SIR" for
Susceptible-Infectious-Removed
}
\item{XYcoordinates}{The XY-coordinates of individuals.}
\item{contact}{Contact network matrix.}
\item{inftime}{The infection times of individuals.}
\item{remtime}{The removal times of individuals when \code{type} = ``SIR''.}
}
}

\seealso{
  \code{\link{plot.epidata}}, \code{\link{epimcmc}}, \code{\link{epilike}}, \code{\link{pred.epi}}.
}

\examples{


## Example 1:  spatial SI model
# generate 100 individuals

x <- runif(100, 0, 10)

y <- runif(100, 0, 10)

covariate <- runif(100, 0, 2)

out1 <- epidata(type = "SI",n = 100, Sformula = ~covariate, tmax = 15,
               sus.par = c(0.1, 0.3), beta = 5.0, x = x, y = y)

# Plots of epidemic progression (optional)

plot(out1, plottype = "spatial")
plot(out1, plottype = "curve", curvetype = "newinfect")

## Example 2: spatial SIR model
# generate infectious period(=3) for 100 individuals

lambda <- rep(3, 100)

out2 <- epidata(type = "SIR", n = 100, tmax = 15, sus.par = 0.3, beta = 5.0, infperiod = lambda,
        x = x, y = y)

plot(out2, plottype = "spatial")
plot(out2, plottype = "curve", curvetype = "newinfect")

## Example 3:   SI network model

contact1 <- matrix(rbinom(10000, 1, 0.1), nrow = 100, ncol = 100)

contact2 <- matrix(rbinom(10000, 1, 0.1), nrow = 100, ncol = 100)

diag(contact1[,] ) <- 0

diag(contact2[,] ) <- 0

contact <- array(c(contact1, contact2), dim = c(100, 100, 2))

out3 <- epidata(type = "SI", n = 100, tmax = 15, sus.par = 0.3, beta = c(3.0, 5.0),
        contact = contact)
plot(out3, plottype = "curve", curvetype = "complete")
plot(out3, plottype = "curve", curvetype = "susceptible")
plot(out3, plottype = "curve", curvetype = "newinfect")
plot(out3, plottype = "curve", curvetype = "totalinfect")

}
\references{

Deardon R, Brooks, S. P., Grenfell, B. T., Keeling, M. J., Tildesley, M. J., Savill, N. J., Shaw, D. J.,  Woolhouse, M. E. (2010). Inference for individual level models of infectious diseases in large populations. \emph{Statistica Sinica}, 20, 239-261.

Rob Deardon, Xuan Fang, and Grace P. S. Kwong (2014).
 Statistical modelling of spatio-temporal infectious disease transmission in analyzing and modeling Spatial and temporal dynamics of infectious diseases, \emph{(Ed: D. Chen, B. Moulin, J. Wu), John Wiley & Sons.} Chapter 11.
}
