\name{epispatial}
\alias{epispatial}

\title{ Spatial Plot  }

\description{Produce spatial plots of epidemicprogression over time}

\usage{epispatial (type, x, y, inftime, removaltime = NULL, time = NULL, tmin = NULL)
}

\arguments{
\item{type}{
 Type of compartment framework, with the choice of "SI" for
Susceptible-Infectious diseases and "SIR" for
Susceptible-Infectious-Removed
}
 
\item{x}{ X coordinates of individuals}

\item{y }{Y coordinates of individuals}
\item{inftime}{Times at which individuals become infected/infectious }
\item{removaltime }{Times at which individuals are removed}
\item{time}{ Specify time points at which the spatial square is plotted}

\item{tmin}{ Initial time point at which infection occurs, default value is one}

}






 



