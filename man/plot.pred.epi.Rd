\name{plot.pred.epi}
\alias{plot.pred.epi}
\title{S3 method to provide plots of posterior predictive check.}
\description{
Produces various graphs for the output of the posterior predictive check of class \code{pred.epi}.
}

\usage{
\method{plot}{pred.epi}(x, ...)
}

\arguments{
	\item{x}{
		An object of class \code{pred.epi} which is the output of the \code{pred.epi} function.
	}
	\item{...}{
		additional arguments to be passed to \code{plot} generic function.
	}
}
\value{
plot
}
\seealso{
 \code{\link{pred.epi}}, \code{\link{plot.epidata}}, \code{\link{plot.epimcmc}}.
}
