\name{plot.epidata}
\alias{plot.epidata}
\title{S3 method to provide some graphics of epidemic.}
\description{
Produces various graphs summarizing epidemic of class \code{epidata}.
}

\usage{
\method{plot}{epidata}(x, plottype, curvetype = NULL, time_id = NULL, tmin = NULL, timepoints = NULL, ...)
}

\arguments{
	\item{x}{
An object of class \code{epidata} that can be the output of \code{\link{epidata}} or \code{\link{as.epidata}}
}
	\item{plottype}{
		Provide two types of plots. When \code{plottype} = "curve", graphs of various epidemic curves can be produced, and when \code{plottype} = "spatial", spatial plots of epidemic progression over time is produced.
		}
	\item{curvetype}{
		It has four options: "complete", "susceptible", "newinfect", and "totalinfect". See details for more information.
		}
	\item{time_id}{
		Specify time points at which the spatial square is plotted.
		}
	\item{tmin}{
		Initial time point at which infection occurs, default value is one.
		}
	\item{timepoints}{
		Specify time points at which the curve is plotted
		}
	\item{...}{
		........
		}
}
\details{
	The argument \code{plottype} has two options. When \code{plottype}=``spatial'' spatial plots are produced for the epidemic progression over time, and when \code{plottype}=``curve'', the argument \code{curvetype} has to be specified to one of the four available options: "complete" for plotting the number of susceptible, infected and removed individuals at each time point, "susceptible" for  plotting  the number of susceptible individuals at each time point, "newinfect" for plotting  the number of newly infected individuals at each time point, and "totalinfect" for  plotting  the cumulative number of infected individuals at each time point.
}
\value{
plot
}
\seealso{
\code{\link{epidata}}, \code{\link{plot.epimcmc}}, \code{\link{plot.pred.epi}}.
}
\examples{


## Example :  spatial SI model
# generate 100 individuals

x <- runif(100, 0, 10)

y <- runif(100, 0, 10)

covariate <- runif(100, 0, 2)

out1 <- epidata(type = "SI",n = 100, Sformula = ~covariate, tmax = 15,
               sus.par = c(0.1, 0.3), beta = 5.0, x = x, y = y)

# Plots of epidemic progression 

plot(out1, plottype = "spatial")
plot(out1, plottype = "curve", curvetype = "newinfect")

}

