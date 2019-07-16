################################################################################
# Part of the R/EpiILM package
#
# AUTHORS:
#         Vineetha Warriyar. K. V. <vineethawarriyar.kod@ucalgary.ca> and
#         Rob Deardon <robert.deardon@ucalgary.ca>
#
# HISTORY:
#          Version 1.0: 2017-04-14
#          Version 1.1: 2017-04-17
#
# Free software under the terms of the GNU General Public License, version 2,
# a copy of which is available at http://www.r-project.org/Licenses/.
################################################################################

epicurve <- function(type, plottype, inftime, removaltime = NULL, tmin = NULL, timepoints = NULL) {
  
  # Error checks for input arguments
  if (is.null(type) || !(type %in% c("SI", "SIR"))) {
    stop("epicurve: Specify type as \"SI\" or \"SIR\" ", call. = FALSE)
  }
  
  if (is.null(plottype) || !(plottype %in% c("complete","susceptible","totalinfect","newinfect"))) {
    stop("epicurve: Specify plottype as \"complete\" , \"susceptible\",\"totalinfect\" or  \"newinfect\"", call. = FALSE)
  }
  
  n <- length(inftime)
  
  if (is.null(removaltime) && type == "SIR") {
    stop(' epicurve: Specify removal time')
  }
  
  if (!is.null(removaltime)) {
      if (length(removaltime) != n) {
        stop('epicurve: Length of removaltime is not compatible')
      }
      if (type == "SI") {
        stop('epicurve: Type must be "SIR"')
      }
  }
  
  if (is.null(tmin)) {
  tmin <- 1
  }
  
  # initializations
  totalinf <- rep(0)
  sus      <- rep(0)
  newinf   <- rep(0)
  removed  <- rep(0)
  tmax     <- max(inftime)
  time     <- rep(tmin:tmax)
  
  # plot for Susceptible-Infectious (SI)
  if (type == "SI") {
    for (i in tmin:tmax) {
      newinf[i] <- length(inftime[inftime==i])
      xc   <- subset(inftime, inftime <= i & inftime != 0)
      totalinf[i] <- length(xc)
      sus[i] <- n - totalinf[i]
    }
    if (tmin > 1) {
      newinf <- newinf[tmin:tmax]
      totalinf <- totalinf[tmin:tmax]
      sus <- sus[tmin:tmax]
    }
   
   # plot for infected and susceptible individuals
    if (plottype == "complete") {
    plot(x    = time,
         y    = sus,
         xlim = c(tmin, tmax),
         ylim = c(1, n),
         main = "Epidemic Curves",
         ylab = "Number of individuals ",
         xlab = "time",
         type = "l",
         lwd  = 2,
         cex  = 0.5,
         xaxt = "n")
    lines(time, totalinf, col = "red",lwd = 2)
    axis(1, at=1:tmax)
    legend("topright", inset = .001, cex = 0.8, bty = "n", c("Infected", "Susceptible"),
           horiz = TRUE, lty = c(1, 1), lwd = c(2, 2),col = c("red", "black"))
    }
   
   # plot for infected and susceptible individuals with specified time points
    if ((plottype == "complete") & (!is.null(timepoints))) {
      plot(x    = time,
           y    = sus,
           xlim = c(timepoints[1], timepoints[2]),
           ylim = c(1,n),
           main = "Epidemic Curves",
           ylab = "Number of individuals",
           xlab = "time",
           type = "l",
           lwd  = 2,
           cex  = 0.5,
           xaxt = "n")
    lines(time, totalinf, col = "red", lwd = 2)
    axis(1, at=1:timepoints[2])
    legend("topright", inset = .001, cex = 0.8, bty = "n", c("Infected", "Susceptible"),
           horiz = TRUE, lty = c(1, 1), lwd = c(2, 2), col = c("red", "black"))
    }
  # end if - SI model
  }
  
  # Plot for Susceptible-Infectious-Removed (SIR)
  if (type == "SIR") {
    dat <- data.frame(inftime, removaltime)
    for (i in tmin:tmax) {
      xcc <- subset(dat, inftime <= i & inftime != 0 & i<removaltime)
      totalinf[i] <- length(xcc$inftime)
    }
    for (i in tmin:tmax) {
      newinf[i] <- length(inftime[inftime==i])
    }
    for (i in tmin:tmax) {
      xcc <- subset(dat, inftime<=i & inftime != 0)
      xc <- subset(xcc, i >= removaltime)
      removed[i] <- length(xc$inftime)
      sus[i] <- n - length(xcc$inftime)
    }
    if (tmin>1) {
      newinf   <- newinf[tmin:tmax]
      totalinf <- totalinf[tmin:tmax]
      removed  <- removed[tmin:tmax]
      sus      <- sus[tmin:tmax]
    }
    
    # plot for infected and susceptible individuals
    if (plottype == "complete") {
      plot(x    = time,
           y    = sus,
           xlim = c(tmin, tmax),
           ylim = c(1, n),
           main = "Epidemic Curve",
           ylab = "Number of individuals",
           xlab = "time",
           type = "l",
           lwd  = 2,
           cex  = 0.5,
           xaxt = "n")
      lines(time, totalinf, col = "red", lwd = 2)
      lines(time, removed, col = "blue", lwd = 2)
      axis(1, at=1:tmax)
      legend("topright", inset = .001, cex = 0.8, bty = "n", c("Infected", "Susceptible", "Removed"),
             horiz = TRUE, lty = c(1, 1, 1), lwd = c(2, 2, 2), col = c("red", "black", "blue"))
    }
    
    # plot for infected and susceptible individuals with specified time points
    if ((plottype == "complete") & (!is.null(timepoints))) {
      plot(x    = time,
           y    = sus,
           xlim = c(timepoints[1], timepoints[2]),
           ylim = c(1,n),
           main = "Epidemic Curve",
           ylab = "Number of individuals",
           xlab = "time",
           type ="l",
           lwd  = 2,
           cex  = 0.5,
           xaxt = "n")
      lines(time, totalinf, col = "red", lwd = 2)
      lines(time, removed, col = "blue", lwd = 2)
      axis(1, at=1:timepoints[2])
      legend("topright", inset = .001, bty = "n", c("Infected", "Susceptible", "Removed"),
             horiz = TRUE, lty = c(1, 1, 1), lwd = c(2, 2, 2), col = c("red", "black", "blue"))
    }
  # end if - SIR model
  }
  
  # plot for total infected individuals
  if (plottype == "totalinfect") {
    plot(x    = time,
         y    = totalinf,
         xlim = c(min(time), max(time)),
         ylim = c(0, (max(totalinf)+1)),
         main = "Epidemic Curve",
         ylab = "Number of infected individuals",
         xlab = "time",
         type = "b",
         pch  = 20,
         lwd  = 2,
         xaxt = "n")
    axis(1, at=1:max(time))
  }
  
  # plot for newly infected individuals
  if (plottype == "newinfect") {
    plot(x    = time,
         y    = newinf,
         xlim = c(min(time), max(time)),
         ylim = c(0, (max(newinf)+1)),
         main = "Epidemic Curve",
         ylab = "Number of new infections",
         xlab = "time",
         type = "b",
         pch  = 20,
         lwd  = 2,
         xaxt = "n")
    axis(1, at=1:max(time))
  }
  
  # plot for susceptible individuals
  if (plottype == "susceptible") {
    plot(x    = time,
         y    = sus,
         xlim = c(min(time), max(time)),
         ylim = c(0, (max(sus)+1)),
         main = "Epidemic Curve",
         ylab = "Number of susceptibles",
         xlab = "time",
         type = "b",
         pch  = 20,
         lwd  = 2,
         xaxt = "n")
    axis(1, at=1:max(time))
  }
  
  # with specified time points
  if (!is.null(timepoints)) {
  
  # plot for total infected individuals
  if (plottype == "totalinfect") {
    plot(x    = time,
         y    = totalinf,
         xlim = c(timepoints[1], timepoints[2]),
         ylim = c(1, n),
         main = "Epidemic Curve",
         ylab = "Number of infected individuals",
         xlab = "time",
         type = "b",
         pch  = 20,
         lwd  = 2,
         xaxt = "n")
    axis(1, at=1:timepoints[2])
  }
  
  # plot for newly infected individuals
  if (plottype == "newinfect") {
    plot(x    = time,
         y    = newinf,
         xlim = c(timepoints[1], timepoints[2]),
         ylim = c(0, (max(newinf)+1)),
         main = "Epidemic Curve",
         ylab = "Number of new infections",
         xlab = "time",
         type = "b",
         pch  = 20,
         lwd  = 2,
         xaxt = "n")
    axis(1, at=1:timepoints[2])
  }
  
  # plot for susceptible individuals
  if (plottype == "susceptible") {
    plot(x    = time,
         y    = sus,
         xlim = c(timepoints[1], timepoints[2]),
         ylim = c(1, n),
         main = "Epidemic curve",
         ylab = "Number of susceptibles",
         xlab = "time",
         type = "b",
         pch  = 20,
         lwd  = 2,
         xaxt = "n")
    axis(1, at=1:timepoints[2])
  }
# end if for timepoints specification
  }
# End of function
}




