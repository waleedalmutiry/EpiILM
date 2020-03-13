\name{tswv}
\docType{data}
\alias{tswv}
\title{Tomato Spotted Wilt Virus (TSWV) data}
\description{
  Data extracted from Hughes et al. (1997). Data obtained from a field experiment as the spatial dynamics of tomato spotted wilt virus (tswv). 
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