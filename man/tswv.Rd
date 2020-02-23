\name{tswv}
\docType{data}
\alias{tswv}
\title{Tomato Spotted Wilt Virus (TSWV) data}
\description{
  A real data example for the use of ILM model in cereal crops is described and the data is extracted from plots presented in Hughes et al. (1997). The data show the locations of healthy plants and plants with TSWV diseases in a field plot on 18 December 1929.
}
\usage{data(tswv)}
\format{
  A data frame with following variables
   \describe{
    \item{x}{ X coordinate}
    \item{y}{ Y coordinate}
    \item{inftime}{Infection times}
    \item{removaltime}{Times at which individuals are removed}
  }
}

\examples{
data("tswv")

x <- tswv$x 
y <- tswv$y 
inftime <- tswv$inftime 
removaltime <- tswv$removaltime 
infperiod <- rep(3, length(x))  

# change to epilate object
 epidat.tswv <- as.epidata(type = "SIR", n = 520, x = x, y = y, 
                            inftime = inftime, infperiod = infperiod)

# plot
plot(epidat.tswv, plottype = "spatial", tmin = 2)
}

\references{
Hughes, G., McRoberts,N., Madden, L.V., Nelson, S. C. (1997). Validating mathematical models of plant disease progress in space and time. \emph{IMA Journal of Mathematics Applied in Medicine and Biology}, 14, 85-112.
}