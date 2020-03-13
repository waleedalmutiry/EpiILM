\name{as.epidata}
\alias{as.epidata}

\title{
Discrete Time level information of an Epidemic
}
\description{
This function allows the user to generate objects of class \code{"epidata"}. The output of this function provides information required to be used in the other functions in the package.
}
\usage{
as.epidata (type, n, x = NULL, y = NULL, inftime, infperiod = NULL, contact = NULL)
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

\value{
An object of class \code{epidata} is returned containing the following:
\describe{
\item{type}{Type of compartment framework, with the choice of "SI" for
Susceptible-Infectious diseases and "SIR" for
Susceptible-Infectious-Removed.
}
\item{XYcoordinates}{The XY-coordinates of individuals if distance-based ILM is used, otherwise is NULL.}
\item{contact}{Contact network matrix if network-based ILM is used, otherwise is NULL.}
\item{inftime}{The infection times of individuals.}
\item{remtime}{The removal times of individuals when \code{type} = ``SIR''.}
}
}

\seealso{
  \code{\link{epidata}}, \code{\link{plot.epidata}}.
}


\examples{
# generate 100 X-Y coordinates for a distance-based ILM

x <- runif(100, 0, 10)

y <- runif(100, 0, 10)

# suppose we know the infection times for a spatial SI model based on above X-Y coordinates

inftime  <- rpois(100, 8) 

# Now we can convert above information to an epidata object with

out <- as.epidata(type = "SI", n = 100, x = x, y = y, inftime = inftime ) 
out

}