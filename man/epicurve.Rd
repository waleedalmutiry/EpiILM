\name{epicurve}
\alias{epicurve}

\title{
Epidemic curves }
\description{
Produces graphs of various epidemiccurves}
\usage{
epicurve (type, plottype, inftime, removaltime = NULL, tmin = NULL, timepoints = NULL)
}

\arguments{
  \item{type}{
 Type of compartment framework, with the choice of "SI" for
Susceptible-Infectious diseases and "SIR" for
Susceptible-Infectious-Removed 
}
\item{plottype}{
Plot type to be plotted, with the following choices:
"complete" for plotting the number of susceptible, infected and removed individuals at each time point, "susceptible" for  plotting  the number of susceptible individuals at each time point,
"newinfect" for plotting  the number of newly infected individuals at each time point, and
"totalinfect" for  plotting  the cumulative number of infected individuals at each time point

}
  \item{inftime}{
Times at which individuals become infected/infectious 
}
  \item{removaltime}{
Time points where the individual is removed
}
  \item{tmin}{
Initial time point at which infection occurs, default value is one}
 
\item{timepoints}{
 Specify time points at which the curve is plotted}

}




