################################################################################
# Part of the R/EpiILM package
#
# AUTHORS:
#         Waleed Almutiry <wkmtierie@qu.edu.sa>,
#         Vineetha Warriyar. K. V. <vineethawarriyar.kod@ucalgary.ca> and
#         Rob Deardon <robert.deardon@ucalgary.ca>
#
# Free software under the terms of the GNU General Public License, version 2,
# a copy of which is available at http://www.r-project.org/Licenses/.
################################################################################

plot.epidata <- function(x, plottype, curvetype = NULL, time_id = NULL, tmin = NULL, timepoints = NULL, ...) {

    if (class(x) != "epidata") {
        stop("The x must be in a class of \"epidata\"", call. = FALSE)
    } else {

        n <- length(x$inftime)

        if (is.null(tmin)) {
        tmin <- 1
        }

        if (plottype == "curve") {
            # Error checks for input arguments
            if (is.null(curvetype) | !(curvetype %in% c("complete","susceptible","totalinfect","newinfect"))) {
                stop("epicurve: Specify plottype as \"complete\" , \"susceptible\",\"totalinfect\" or  \"newinfect\"", call. = FALSE)
            }

            # initializations
            totalinf <- rep(0)
            sus      <- rep(0)
            newinf   <- rep(0)
            removed  <- rep(0)
            tmax     <- max(x$inftime)
            timerange     <- rep(tmin:tmax)

            # plot for Susceptible-Infectious (SI)
            if (x$type == "SI") {
                for (i in tmin:tmax) {
                    newinf[i] <- length(x$inftime[x$inftime==i])
                    xc   <- subset(x$inftime, x$inftime <= i & x$inftime != 0)
                    totalinf[i] <- length(xc)
                    sus[i] <- n - totalinf[i]
                }
                if (tmin > 1) {
                    newinf <- newinf[tmin:tmax]
                    totalinf <- totalinf[tmin:tmax]
                    sus <- sus[tmin:tmax]
                }


                # plot for infected and susceptible individuals
                if (curvetype == "complete") {

                    plot(timerange, sus,
                        xlim = c(tmin, tmax), ylim = c(1, n),
                        main = "Epidemic Curves", ylab = "Number of individuals ", xlab = "time",
                        type = "l", lwd  = 2, cex  = 0.5, xaxt = "n")
                    lines(timerange, totalinf, col = "red",lwd = 2)
                    axis(1, at=1:tmax)
                    legend("topright", inset = .001, cex = 0.8, bty = "n", c("Infected", "Susceptible"),
                    horiz = TRUE, lty = c(1, 1), lwd = c(2, 2),col = c("red", "black"))

                } else if ((curvetype == "complete") & (!is.null(timepoints))) {

                    # plot for infected and susceptible individuals with specified time points
                    plot(timerange, sus,
                        xlim = c(timepoints[1], timepoints[2]), ylim = c(1,n),
                        main = "Epidemic Curves", ylab = "Number of individuals", xlab = "time",
                        type = "l", lwd  = 2, cex  = 0.5, xaxt = "n")
                    lines(timerange, totalinf, col = "red", lwd = 2)
                    axis(1, at=1:timepoints[2])
                    legend("topright", inset = .001, cex = 0.8, bty = "n", c("Infected", "Susceptible"),
                    horiz = TRUE, lty = c(1, 1), lwd = c(2, 2), col = c("red", "black"))

                } else if ((curvetype == "totalinfect") & (is.null(timepoints))) {

                    plot(timerange, totalinf,
                        xlim = c(min(timerange), max(timerange)), ylim = c(0, (max(totalinf)+1)),
                        main = "Epidemic Curve", ylab = "Number of infected individuals", xlab = "time",
                        type = "b", pch  = 20, lwd  = 2, xaxt = "n")
                    axis(1, at=1:max(timerange))

                } else if ((curvetype == "totalinfect") & (!is.null(timepoints))) {

                    plot(timerange, totalinf,
                        xlim = c(timepoints[1], timepoints[2]), ylim = c(1, n),
                        main = "Epidemic Curve", ylab = "Number of infected individuals", xlab = "time",
                        type = "b", pch  = 20, lwd  = 2, xaxt = "n")
                    axis(1, at=1:timepoints[2])

                } else if ((curvetype == "newinfect") & (is.null(timepoints))) {

                    plot(timerange, newinf,
                       xlim = c(min(timerange), max(timerange)), ylim = c(0, (max(newinf)+1)),
                       main = "Epidemic Curve", ylab = "Number of new infections", xlab = "time",
                       type = "b", pch  = 20, lwd  = 2, xaxt = "n")
                    axis(1, at=1:max(timerange))

                } else if ((curvetype == "newinfect") & (!is.null(timepoints))) {

                    plot(timerange, newinf,
                        xlim = c(timepoints[1], timepoints[2]), ylim = c(0, (max(newinf)+1)),
                        main = "Epidemic Curve", ylab = "Number of new infections", xlab = "time",
                        type = "b", pch  = 20, lwd  = 2, xaxt = "n")
                    axis(1, at=1:timepoints[2])

                } else if ((curvetype == "susceptible") & (is.null(timepoints))) {

                    plot(timerange, sus,
                        xlim = c(min(timerange), max(timerange)), ylim = c(0, (max(sus)+1)),
                        main = "Epidemic Curve", ylab = "Number of susceptibles", xlab = "time",
                        type = "b", pch  = 20, lwd  = 2, xaxt = "n")
                    axis(1, at=1:max(timerange))

                } else if ((curvetype == "susceptible") & (!is.null(timepoints))) {

                    plot(timerange, sus,
                        xlim = c(timepoints[1], timepoints[2]), ylim = c(1, n),
                        main = "Epidemic Curve", ylab = "Number of susceptibles", xlab = "time",
                        type = "b", pch  = 20, lwd  = 2, xaxt = "n")
                    axis(1, at=1:timepoints[2])

                }

            } else if (x$type == "SIR") {
                # Plot for Susceptible-Infectious-Removed (SIR)
                dat <- data.frame(x$inftime, x$remtime)

                for (i in tmin:tmax) {
                    xcc <- subset(dat, x$inftime <= i & x$inftime != 0 & i < x$remtime)
                    totalinf[i] <- length(xcc$x.inftime)
                }

                for (i in tmin:tmax) {
                    newinf[i] <- length(x$inftime[x$inftime==i])
                }

                for (i in tmin:tmax) {
                    xcc <- subset(dat, x$inftime<=i & x$inftime != 0)
                    xc <- subset(xcc, i >= xcc$x.remtime)
                    removed[i] <- length(xc$x.remtime)
                    sus[i] <- n - length(xcc$x.inftime)
                }
                if (tmin>1) {
                    newinf   <- newinf[tmin:tmax]
                    totalinf <- totalinf[tmin:tmax]
                    removed  <- removed[tmin:tmax]
                    sus      <- sus[tmin:tmax]
                }

                # plot for infected and susceptible individuals
                if ((curvetype == "complete") & (is.null(timepoints))) {

                    plot(timerange, sus,
                        xlim = c(tmin, tmax), ylim = c(1, n),
                        main = "Epidemic Curve", ylab = "Number of individuals", xlab = "time",
                        type = "l", lwd  = 2, cex  = 0.5, xaxt = "n")
                    lines(timerange, totalinf, col = "red", lwd = 2)
                    lines(timerange, removed, col = "blue", lwd = 2)
                    axis(1, at=1:tmax)
                    legend("topright", inset = .001, cex = 0.8, bty = "n", c("Infected", "Susceptible", "Removed"),
                    horiz = TRUE, lty = c(1, 1, 1), lwd = c(2, 2, 2), col = c("red", "black", "blue"))

                } else if ((curvetype == "complete") & (!is.null(timepoints))) {
                # plot for infected and susceptible individuals with specified time points

                    plot(timerange, sus,
                        xlim = c(timepoints[1], timepoints[2]), ylim = c(1,n),
                        main = "Epidemic Curve", ylab = "Number of individuals", xlab = "time",
                        type ="l", lwd  = 2, cex  = 0.5, xaxt = "n")
                    lines(timerange, totalinf, col = "red", lwd = 2)
                    lines(timerange, removed, col = "blue", lwd = 2)
                    axis(1, at=1:timepoints[2])
                    legend("topright", inset = .001, bty = "n", c("Infected", "Susceptible", "Removed"),
                    horiz = TRUE, lty = c(1, 1, 1), lwd = c(2, 2, 2), col = c("red", "black", "blue"))

                } else if ((curvetype == "totalinfect") & (is.null(timepoints))) {

                    plot(timerange, totalinf,
                        xlim = c(min(timerange), max(timerange)), ylim = c(0, (max(totalinf)+1)),
                        main = "Epidemic Curve", ylab = "Number of infected individuals", xlab = "time",
                        type = "b", pch  = 20, lwd  = 2, xaxt = "n")
                    axis(1, at=1:max(timerange))

                } else if ((curvetype == "totalinfect") & (!is.null(timepoints))) {

                    plot(timerange, totalinf,
                        xlim = c(timepoints[1], timepoints[2]), ylim = c(0, n),
                        main = "Epidemic Curve", ylab = "Number of infected individuals", xlab = "time",
                        type = "b", pch  = 20, lwd  = 2, xaxt = "n")
                    axis(1, at=1:timepoints[2])

                } else if ((curvetype == "newinfect") & (is.null(timepoints))) {

                    plot(timerange, newinf,
                        xlim = c(min(timerange), max(timerange)), ylim = c(0, (max(newinf)+1)),
                        main = "Epidemic Curve", ylab = "Number of new infections", xlab = "time",
                        type = "b", pch  = 20, lwd  = 2, xaxt = "n")
                    axis(1, at=1:max(timerange))

                } else if ((curvetype == "newinfect") & (!is.null(timepoints))) {

                    plot(timerange, newinf,
                        xlim = c(timepoints[1], timepoints[2]), ylim = c(0, (max(newinf)+1)),
                        main = "Epidemic Curve", ylab = "Number of new infections", xlab = "time",
                        type = "b", pch  = 20, lwd  = 2, xaxt = "n")
                    axis(1, at=1:timepoints[2])

                } else if ((curvetype == "susceptible") & (is.null(timepoints))) {

                    plot(timerange, sus,
                        xlim = c(min(timerange), max(timerange)), ylim = c(0, (max(sus)+1)),
                        main = "Epidemic Curve", ylab = "Number of susceptibles", xlab = "time",
                        type = "b", pch  = 20, lwd  = 2, xaxt = "n")
                    axis(1, at=1:max(timerange))

                } else if ((curvetype == "susceptible") & (!is.null(timepoints))) {

                    plot(timerange, sus,
                        xlim = c(timepoints[1], timepoints[2]), ylim = c(1, n),
                        main = "Epidemic Curve", ylab = "Number of susceptibles", xlab = "time",
                        type = "b", pch  = 20, lwd  = 2, xaxt = "n")
                    axis(1, at=1:timepoints[2])

                }

            }
        } else if (plottype == "spatial") {

			op1 <- par(no.readonly = TRUE)

            old.par <- par(mfrow = c(3,3))
            par(cex = 0.5)
            par(xpd = TRUE)

    # plots for susceptible- infected (SI)
            if (x$type == "SI") {
                dat <- data.frame(x$XYcoordinates, x$inftime)
                # no specific time point (s)
                if (is.null(time_id)) {
                    for (i in tmin:max(x$inftime)) {
                        xcc <- subset(dat, dat$x.inftime <= i & dat$x.inftime != 0)
                        xx    = x$XYcoordinates[,1]
                        yy    = x$XYcoordinates[,2]
                        plot(xx, yy,
                        xlim = c(min(xx), max(xx)),
                        ylim = c(min(yy), max(yy)),
                        sub  = paste("time ", i))
                        points(xcc[,1], xcc[,2], pch = 16, col = "red")
                        op <- par(usr = c(0, 1, 0, 1), xpd = NA)
                        legend(-.03, 1.15, c("Infected", "Susceptible"), col = c("red", "black"),
                        pch = c(16, 21), bty = "n", horiz = TRUE)
                    }
                } else if (!is.null(time_id)) {
                    for (i in 1:length(time_id)) {
                        xcc <- subset(dat, dat$x.inftime <= time_id[i] & dat$x.inftime != 0)
                        xx    = x$XYcoordinates[,1]
                        yy    = x$XYcoordinates[,2]
                        plot(xx, yy,
                        xlim = c(min(xx), max(xx)),
                        ylim = c(min(yy), max(yy)),
                        sub  = paste("time ", time_id[i]))
                        points(xcc[,1], xcc[,2], pch = 16, col = "red")
                        op <- par(usr = c(0, 1, 0, 1), xpd = NA)
                        legend(-.03, 1.15, c("Infected", "Susceptible"), col = c("red", "black"),
                        pch = c(16,21), bty = "n", horiz = TRUE)
                    }
                }
            # end if - SI model
            } else if (x$type == "SIR") {
                dat <- data.frame(x$XYcoordinates, x$inftime, x$remtime)
                if (is.null(time_id)) {
                    for(i in tmin:max(x$inftime)) {
                        xcc <- subset(dat, dat$x.inftime <= i & dat$x.inftime != 0)
                        xx    = x$XYcoordinates[,1]
                        yy    = x$XYcoordinates[,2]
                        plot(xx, yy,
                        xlim = c(min(xx), max(xx)),
                        ylim = c(min(yy), max(yy)),
                        pch  = 21,
                        sub  = paste("time ", i))
                        xred <- subset(xcc, i < xcc[,4])
                        points(xred[,1], xred[,2], pch = 16, col = "red")
                        xblue <- subset(xcc, i >= xcc[,4])
                        points(xblue[,1], xblue[,2], pch = 16, col = "blue")
                        op <- par(usr = c(0, 1, 0, 1), xpd = NA)
                        legend(-.03, 1.15, c("Infected","Susceptible", "Removed"), col = c("red", "black", "blue"),
                        pch = c(16, 21, 16), bty = "n", horiz = TRUE)
                    }
                } else if (!is.null(time_id)) {
                    for (i in 1:length(time_id)) {
                        xcc <- subset(dat, dat$x.inftime <= time_id[i] & dat$x.inftime != 0)
                        xx    = x$XYcoordinates[,1]
                        yy    = x$XYcoordinates[,2]
                        plot(xx, yy,
                        xlim = c(min(xx), max(xx)),
                        ylim = c(min(yy), max(yy)),
                        pch  = 21,
                        sub  = paste("time ", time_id[i]))
                        xred <- subset(xcc, time_id[i] < xcc[,4])
                        points(xred[,1], xred[,2], pch = 16, col = "red")
                        xblue <- subset(xcc, time_id[i] >= xcc[,4])
                        points(xblue[,1], xblue[,2], pch = 16, col = "blue")
                        op <- par(usr = c(0, 1, 0, 1), xpd = NA)
                        legend(-.03, 1.15, c("Infected", "Susceptible", "Removed"), col = c("red", "black", "blue"),
                        pch = c(16, 21, 16), bty = "n", horiz = TRUE)
                    }
                }
            # end if - SIR model
            }

            par(old.par)

			on.exit(par(op1))

    # End function
        } else {
            stop("The plottype option should be either \"curve\" or \"spatial\"", call. = FALSE)
        }
    }
    # End of function
}
