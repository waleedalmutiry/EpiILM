################################################################################
# Part of the R/EpiILM package
#
# AUTHORS:
#         Vineetha Warriyar. K. V. <vineethawarriyar.kod@ucalgary.ca> and
#         Rob Deardon <robert.deardon@ucalgary.ca>
#
#
# Free software under the terms of the GNU General Public License, version 2,
# a copy of which is available at http://www.r-project.org/Licenses/.
################################################################################

epispatial <- function(type, x, y, inftime, removaltime = NULL, time = NULL, tmin = NULL) {
 
# input argument- error checks
  if (is.null(type) || !(type %in% c("SI", "SIR"))) {
    stop("epiplot: Specify type as \"SI\" or \"SIR\" ", call. = FALSE)
  }
  
  n <- length(inftime)
  
  if ((length(y) != n) || (length(x) != n)) {
    stop('epispatial: Length of x or y is not compatible ')
  }
  if (is.null(removaltime) && type == "SIR") {
    stop(' epispatial: Specify removaltime ')
  }
  if (!is.null(removaltime)) {
      if (length(removaltime) != n) {
        stop('epispatial: Length of removaltime is not compatible')
      }
      if (type == "SI") {
        stop('epispatial: Type must be "SIR"')
      }
  }
  if (is.null(tmin)) {
    tmin <- 1
  }
  
  old.par <- par(mfrow = c(3,3))
  par(cex = 0.5)
  par(xpd = TRUE)
  
# plots for susceptible- infected (SI)
  if (type == "SI") {
    dat <- data.frame(x, y, inftime)
    # no specific time point (s)
    if (is.null(time)) {
      for (i in tmin:max(inftime)) {
        xcc <- subset(dat, inftime <= i & inftime != 0)
        plot(x    = x,
             y    = y,
             xlim = c(min(x), max(x)),
             ylim = c(min(y), max(y)),
             sub  = paste("time ", i))
        points(xcc$x, xcc$y, pch = 16, col = "red")
        op <- par(usr = c(0, 1, 0, 1), xpd = NA)
        legend(-.03, 1.15, c("Infected", "Susceptible"), col = c("red", "black"),
               pch = c(16, 21), bty = "n", horiz = TRUE)
      }
    }
    # specific time points
    if (!is.null(time)) {
        for (i in 1:length(time)) {
            xcc <- subset(dat, inftime <= time[i] & inftime != 0)
            plot(x    = x,
                 y    = y,
                 xlim = c(min(x), max(x)),
                 ylim = c(min(y), max(y)),
                 sub  = paste("time ", time[i]))
            points(xcc$x, xcc$y, pch = 16, col = "red")
            op <- par(usr = c(0, 1, 0, 1), xpd = NA)
            legend(-.03, 1.15, c("Infected", "Susceptible"), col = c("red", "black"),
                   pch = c(16,21), bty = "n", horiz = TRUE)
        }
    }
  # end if - SI model
  }

# plots for susceptible- infected-removed  (SIR)
  if (type == "SIR") {
    dat <- data.frame(x, y, inftime, removaltime)
    if (is.null(time)) {
      for(i in tmin:max(inftime)) {
        xcc <- subset(dat, inftime <= i & inftime != 0)
        plot(x    = x,
             y    = y,
             xlim = c(min(x), max(x)),
             ylim = c(min(y), max(y)),
             pch  = 21,
             sub  = paste("time ", i))
        xred <- subset(xcc, i < xcc$removaltime)
        points(xred$x, xred$y, pch = 16, col = "red")
        xblue <- subset(xcc, i >= xcc$removaltime)
        points(xblue$x, xblue$y, pch = 4, col = "blue")
        op <- par(usr = c(0, 1, 0, 1), xpd = NA)
        legend(-.03, 1.15, c("Infected","Susceptible", "Removed"), col = c("red", "black", "blue"),
               pch = c(16, 21, 4), bty = "n", horiz = TRUE)
      }
    }

# with specific time points
  if (!is.null(time)) {
    for (i in 1:length(time)) {
      xcc <- subset(dat, inftime <= time[i] & inftime != 0)
      plot(x    = x,
           y    = y,
           xlim = c(min(x), max(x)),
           ylim = c(min(y), max(y)),
           pch  = 21,
           sub  = paste("time ", time[i]))
      xred <- subset(xcc, time[i] < xcc$removaltime)
      points(xred$x, xred$y, pch = 16, col = "red")
      xblue <- subset(xcc, time[i] >= xcc$removaltime)
      points(xblue$x, xblue$y, pch = 4, col = "blue")
      op <- par(usr = c(0, 1, 0, 1), xpd = NA)
      legend(-.03, 1.15, c("Infected", "Susceptible", "Removed"), col = c("red", "black", "blue"),
             pch = c(16, 21, 4), bty = "n", horiz = TRUE)
    }
  }
  # end if - SIR model
  }
  par(old.par)
# End function
}

