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

    if (!is(x, "epidata")) {
        stop("The x must be in a class of \"epidata\"", call. = FALSE)
    } else {

        n <- length(x$inftime)

        if (is.null(tmin)) {
        tmin <- 1
        }

        if (plottype == "curve") {
            # Error checks for input arguments
            if (any(is.null(curvetype) | !(curvetype %in% c("complete","susceptible","totalinfect","newinfect"))) == TRUE) {
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
                if (all((curvetype == "complete") & (is.null(timepoints))) == TRUE) {

				    op1 <- par(no.readonly = TRUE)
					op <- par(mar = c(5.9, 4.0, 1.5, 0.5), omi = c(0.2, 0.05, 0.15, 0.15), 
						      mfrow = c(1,1))
				    
                    plot(timerange, sus,
                        xlim = c(tmin, tmax), ylim = c(1, n),
                        main = "Epidemic Curves", ylab = "Number of individuals ", xlab = "Time",
                        type = "l", lwd  = 2, cex  = 1, xaxt = "n")
                    lines(timerange, totalinf, col = "red", lwd = 2)
                    axis(1, at = 1:tmax)
					opar <- par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0),
					mar = c(0, 0, 0, 0), new = TRUE)
					on.exit(par(opar), add = TRUE)
					plot(0, 0, type = 'n', bty = 'n', xaxt = 'n', yaxt = 'n')
					legend("bottom", c("Infected", "Susceptible"), col = c("red", "black"),
					lty = c(1, 1), lwd = c(2, 2), bty = "n", horiz = TRUE, cex = 1)
					op <- par(mar = c(5.9, 4.0, 1.5, 0.5), omi = c(0.2, 0.05, 0.15, 0.15), mfrow = c(1,1))
					on.exit(par(op), add = TRUE)

                } else if (all((curvetype == "complete") & (!is.null(timepoints))) == TRUE) {

				    op1 <- par(no.readonly = TRUE)
					op <- par(mar = c(5.9, 4.0, 1.5, 0.5), omi = c(0.2, 0.05, 0.15, 0.15), 
						      mfrow = c(1,1))
                    # plot for infected and susceptible individuals with specified time points
                    plot(timerange, sus,
                        xlim = c(timepoints[1], timepoints[2]), ylim = c(1,n),
                        main = "Epidemic Curves", ylab = "Number of individuals", xlab = "time",
                        type = "l", lwd  = 2, cex  = 1, xaxt = "n")
                    lines(timerange, totalinf, col = "red", lwd = 2)
                    axis(1, at=1:timepoints[2])
					opar <- par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0),
					mar = c(0, 0, 0, 0), new = TRUE)
					on.exit(par(opar), add = TRUE)
					plot(0, 0, type = 'n', bty = 'n', xaxt = 'n', yaxt = 'n')
					legend("bottom", c("Infected", "Susceptible"), col = c("red", "black"),
					lty = c(1, 1), lwd = c(2, 2), bty = "n", horiz = TRUE, cex = 1)
					op <- par(mar = c(5.9, 4.0, 1.5, 0.5), omi = c(0.2, 0.05, 0.15, 0.15), mfrow = c(1,1))
					on.exit(par(op), add = TRUE)

                } else if (all((curvetype == "totalinfect") & (is.null(timepoints))) == TRUE) {

				    op1 <- par(no.readonly = TRUE)
					op <- par(mar = c(5.9, 4.0, 1.5, 0.5), omi = c(0.2, 0.05, 0.15, 0.15), 
						      mfrow = c(1,1))
                    plot(timerange, totalinf,
                        xlim = c(min(timerange), max(timerange)), ylim = c(0, (max(totalinf)+1)),
                        main = "Epidemic Curve", ylab = "Number of infected individuals", xlab = "time",
                        type = "b", pch  = 20, lwd  = 2, xaxt = "n")
                    axis(1, at=1:max(timerange))

                } else if (all((curvetype == "totalinfect") & (!is.null(timepoints))) == TRUE) {

				    op1 <- par(no.readonly = TRUE)
					op <- par(mar = c(5.9, 4.0, 1.5, 0.5), omi = c(0.2, 0.05, 0.15, 0.15), 
						      mfrow = c(1,1))
                    plot(timerange, totalinf,
                        xlim = c(timepoints[1], timepoints[2]), ylim = c(1, n),
                        main = "Epidemic Curve", ylab = "Number of infected individuals", xlab = "time",
                        type = "b", pch  = 20, lwd  = 2, xaxt = "n")
                    axis(1, at=1:timepoints[2])

                } else if (all((curvetype == "newinfect") & (is.null(timepoints))) == TRUE) {

				    op1 <- par(no.readonly = TRUE)
					op <- par(mar = c(5.9, 4.0, 1.5, 0.5), omi = c(0.2, 0.05, 0.15, 0.15), 
						      mfrow = c(1,1))
                    plot(timerange, newinf,
                       xlim = c(min(timerange), max(timerange)), ylim = c(0, (max(newinf)+1)),
                       main = "Epidemic Curve", ylab = "Number of new infections", xlab = "time",
                       type = "b", pch  = 20, lwd  = 2, xaxt = "n")
                    axis(1, at=1:max(timerange))

                } else if (all((curvetype == "newinfect") & (!is.null(timepoints))) == TRUE) {

				    op1 <- par(no.readonly = TRUE)
					op <- par(mar = c(5.9, 4.0, 1.5, 0.5), omi = c(0.2, 0.05, 0.15, 0.15), 
						      mfrow = c(1,1))
                    plot(timerange, newinf,
                        xlim = c(timepoints[1], timepoints[2]), ylim = c(0, (max(newinf)+1)),
                        main = "Epidemic Curve", ylab = "Number of new infections", xlab = "time",
                        type = "b", pch  = 20, lwd  = 2, xaxt = "n")
                    axis(1, at=1:timepoints[2])

                } else if (all((curvetype == "susceptible") & (is.null(timepoints))) == TRUE) {

				    op1 <- par(no.readonly = TRUE)
					op <- par(mar = c(5.9, 4.0, 1.5, 0.5), omi = c(0.2, 0.05, 0.15, 0.15), 
						      mfrow = c(1,1))
                    plot(timerange, sus,
                        xlim = c(min(timerange), max(timerange)), ylim = c(0, (max(sus)+1)),
                        main = "Epidemic Curve", ylab = "Number of susceptibles", xlab = "time",
                        type = "b", pch  = 20, lwd  = 2, xaxt = "n")
                    axis(1, at=1:max(timerange))

                } else if (all((curvetype == "susceptible") & (!is.null(timepoints))) == TRUE) {

				    op1 <- par(no.readonly = TRUE)
					op <- par(mar = c(5.9, 4.0, 1.5, 0.5), omi = c(0.2, 0.05, 0.15, 0.15), 
						      mfrow = c(1,1))
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
                if (all((curvetype == "complete") & (is.null(timepoints))) == TRUE) {

				    op1 <- par(no.readonly = TRUE)
					op <- par(mar = c(5.9, 4.0, 1.5, 0.5), omi = c(0.2, 0.05, 0.15, 0.15), 
						      mfrow = c(1,1))
                    plot(timerange, sus,
                        xlim = c(tmin, tmax), ylim = c(1, n),
                        main = "Epidemic Curve", ylab = "Number of individuals", xlab = "time",
                        type = "l", lwd  = 2, cex  = 1, xaxt = "n")
                    lines(timerange, totalinf, col = "red", lwd = 2)
                    lines(timerange, removed, col = "blue", lwd = 2)
                    axis(1, at=1:tmax)
					opar <- par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0),
					mar = c(0, 0, 0, 0), new = TRUE)
					on.exit(par(opar), add = TRUE)
					plot(0, 0, type = 'n', bty = 'n', xaxt = 'n', yaxt = 'n')
					legend("bottom", c("Infected", "Susceptible", "Removed"), col = c("red", "black","blue"),
					lty = c(1, 1, 1), lwd = c(2, 2, 2), bty = "n", horiz = TRUE, cex = 1)
					op <- par(mar = c(5.9, 4.0, 1.5, 0.5), omi = c(0.2, 0.05, 0.15, 0.15), mfrow = c(1,1))
					on.exit(par(op), add = TRUE)

                } else if (all((curvetype == "complete") & (!is.null(timepoints))) == TRUE) {
                # plot for infected and susceptible individuals with specified time points

				    op1 <- par(no.readonly = TRUE)
					op <- par(mar = c(5.9, 4.0, 1.5, 0.5), omi = c(0.2, 0.05, 0.15, 0.15), 
						      mfrow = c(1,1))
                    plot(timerange, sus,
                        xlim = c(timepoints[1], timepoints[2]), ylim = c(1,n),
                        main = "Epidemic Curve", ylab = "Number of individuals", xlab = "time",
                        type ="l", lwd  = 2, cex  = 0.5, xaxt = "n")
                    lines(timerange, totalinf, col = "red", lwd = 2)
                    lines(timerange, removed, col = "blue", lwd = 2)
                    axis(1, at=1:timepoints[2])
					opar <- par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0),
					mar = c(0, 0, 0, 0), new = TRUE)
					on.exit(par(opar), add = TRUE)
					plot(0, 0, type = 'n', bty = 'n', xaxt = 'n', yaxt = 'n')
					legend("bottom", c("Infected", "Susceptible", "Removed"), col = c("red", "black","blue"),
					lty = c(1, 1, 1), lwd = c(2, 2, 2), bty = "n", horiz = TRUE, cex = 1)
					op <- par(mar = c(5.9, 4.0, 1.5, 0.5), omi = c(0.2, 0.05, 0.15, 0.15), mfrow = c(1,1))
					on.exit(par(op), add = TRUE)

                } else if (all((curvetype == "totalinfect") & (is.null(timepoints))) == TRUE) {

				    op1 <- par(no.readonly = TRUE)
					op <- par(mar = c(5.9, 4.0, 1.5, 0.5), omi = c(0.2, 0.05, 0.15, 0.15), 
						      mfrow = c(1,1))
                    plot(timerange, totalinf,
                        xlim = c(min(timerange), max(timerange)), ylim = c(0, (max(totalinf)+1)),
                        main = "Epidemic Curve", ylab = "Number of infected individuals", xlab = "time",
                        type = "b", pch  = 20, lwd  = 2, xaxt = "n")
                    axis(1, at=1:max(timerange))

                } else if (all((curvetype == "totalinfect") & (!is.null(timepoints))) == TRUE) {

				    op1 <- par(no.readonly = TRUE)
					op <- par(mar = c(5.9, 4.0, 1.5, 0.5), omi = c(0.2, 0.05, 0.15, 0.15), 
						      mfrow = c(1,1))
                    plot(timerange, totalinf,
                        xlim = c(timepoints[1], timepoints[2]), ylim = c(0, n),
                        main = "Epidemic Curve", ylab = "Number of infected individuals", xlab = "time",
                        type = "b", pch  = 20, lwd  = 2, xaxt = "n")
                    axis(1, at=1:timepoints[2])

                } else if (all((curvetype == "newinfect") & (is.null(timepoints))) == TRUE) {

				    op1 <- par(no.readonly = TRUE)
					op <- par(mar = c(5.9, 4.0, 1.5, 0.5), omi = c(0.2, 0.05, 0.15, 0.15), 
						      mfrow = c(1,1))
                    plot(timerange, newinf,
                        xlim = c(min(timerange), max(timerange)), ylim = c(0, (max(newinf)+1)),
                        main = "Epidemic Curve", ylab = "Number of new infections", xlab = "time",
                        type = "b", pch  = 20, lwd  = 2, xaxt = "n")
                    axis(1, at=1:max(timerange))

                } else if (all((curvetype == "newinfect") & (!is.null(timepoints))) == TRUE) {

				    op1 <- par(no.readonly = TRUE)
					op <- par(mar = c(5.9, 4.0, 1.5, 0.5), omi = c(0.2, 0.05, 0.15, 0.15), 
						      mfrow = c(1,1))
                    plot(timerange, newinf,
                        xlim = c(timepoints[1], timepoints[2]), ylim = c(0, (max(newinf)+1)),
                        main = "Epidemic Curve", ylab = "Number of new infections", xlab = "time",
                        type = "b", pch  = 20, lwd  = 2, xaxt = "n")
                    axis(1, at=1:timepoints[2])

                } else if (all((curvetype == "susceptible") & (is.null(timepoints))) == TRUE) {

				    op1 <- par(no.readonly = TRUE)
					op <- par(mar = c(5.9, 4.0, 1.5, 0.5), omi = c(0.2, 0.05, 0.15, 0.15), 
						      mfrow = c(1,1))
                    plot(timerange, sus,
                        xlim = c(min(timerange), max(timerange)), ylim = c(0, (max(sus)+1)),
                        main = "Epidemic Curve", ylab = "Number of susceptibles", xlab = "time",
                        type = "b", pch  = 20, lwd  = 2, xaxt = "n")
                    axis(1, at=1:max(timerange))

                } else if (all((curvetype == "susceptible") & (!is.null(timepoints))) == TRUE) {

				    op1 <- par(no.readonly = TRUE)
					op <- par(mar = c(5.9, 4.0, 1.5, 0.5), omi = c(0.2, 0.05, 0.15, 0.15), 
						      mfrow = c(1,1))
                    plot(timerange, sus,
                        xlim = c(timepoints[1], timepoints[2]), ylim = c(1, n),
                        main = "Epidemic Curve", ylab = "Number of susceptibles", xlab = "time",
                        type = "b", pch  = 20, lwd  = 2, xaxt = "n")
                    axis(1, at=1:timepoints[2])

                }

            }
        } else if (plottype == "spatial") {

			op1 <- par(no.readonly = TRUE)

    # plots for susceptible- infected (SI)
            if (x$type == "SI") {
                dat <- data.frame(x$XYcoordinates, x$inftime)
                # no specific time point (s)
                if (is.null(time_id)) {
					ntimes <- max(x$inftime) - tmin + 1
					mfrow1 <- switch(min(ntimes,13), c(1,1), c(1,2), c(2,2), c(2,2), c(3,2), c(3,2), 
					c(3,3), c(3,3), c(3,3), c(3,2), c(3,2), c(3,2), c(3,3))
					
					sepwindow <- seq(prod(mfrow1), prod(mfrow1)*ceiling(ntimes/prod(mfrow1)), 
					by = prod(mfrow1))
					
					op <- par(mar = c(5.9, 4.0, 1.5, 0.5), omi = c(0.2, 0.05, 0.15, 0.15), 
						      mfrow = mfrow1)
					u <- 1
					
                    for (i in tmin:max(x$inftime)) {
                        xcc <- subset(dat, dat$x.inftime <= i & dat$x.inftime != 0)
                        xx    = x$XYcoordinates[,1]
                        yy    = x$XYcoordinates[,2]
                        plot(xx, yy,
                        xlim = c(min(xx), max(xx)),
                        ylim = c(min(yy), max(yy)),
                        main  = paste("time ", i), cex = 1, ...)
                        points(xcc[,1], xcc[,2], pch = 16, col = "red")
						
                        if (any(u == sepwindow) | (u == ntimes)) {
                            opar <- par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0),
                            mar = c(0, 0, 0, 0), new = TRUE)
                            on.exit(par(opar), add = TRUE)
                            plot(0, 0, type = 'n', bty = 'n', xaxt = 'n', yaxt = 'n')
                            legend("bottom", c("Infected", "Susceptible"), col = c("red", "black"),
							pch = c(16, 21), bty = "n", horiz = TRUE, cex = 1.5)
							op <- par(mar = c(5.9, 4.0, 1.5, 0.5), omi = c(0.2, 0.05, 0.15, 0.15), mfrow = mfrow1)
						   	on.exit(par(op), add = TRUE)
						}
						u <- u + 1
                    }

                } else if (!is.null(time_id)) {
					ntimes <- length(time_id)
					mfrow1 <- switch(min(ntimes,13), c(1,1), c(1,2), c(2,2), c(2,2), c(3,2), c(3,2), 
					c(3,3), c(3,3), c(3,3), c(3,2), c(3,2), c(3,2), c(3,3))
					
					sepwindow <- seq(prod(mfrow1), prod(mfrow1)*ceiling(length(time_id)/prod(mfrow1)), 
					by = prod(mfrow1))
										
					op <- par(mar = c(5.9, 4.0, 1.5, 0.5), omi = c(0.2, 0.05, 0.15, 0.15), 
						      mfrow = mfrow1)

                    for (i in 1:length(time_id)) {
                        xcc <- subset(dat, dat$x.inftime <= time_id[i] & dat$x.inftime != 0)
                        xx    = x$XYcoordinates[,1]
                        yy    = x$XYcoordinates[,2]
                        plot(xx, yy,
                        xlim = c(min(xx), max(xx)),
                        ylim = c(min(yy), max(yy)),
                        main  = paste("time ", time_id[i]), cex = 1, ...)
                        points(xcc[,1], xcc[,2], pch = 16, col = "red")
                        if (any(i == sepwindow) | (i == length(time_id))) {
                            opar <- par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0),
                            mar = c(0, 0, 0, 0), new = TRUE)
                            on.exit(par(opar), add = TRUE)
                            plot(0, 0, type = 'n', bty = 'n', xaxt = 'n', yaxt = 'n')
                              legend("bottom", c("Infected", "Susceptible"), col = c("red", "black"),
							 pch = c(16, 21), bty = "n", horiz = TRUE, cex = 1.5)
							op <- par(mar = c(5.9, 4.0, 1.5, 0.5), omi = c(0.2, 0.05, 0.15, 0.15), mfrow = mfrow1)
						   	on.exit(par(op), add = TRUE)
						}
                    }
                }
            # end if - SI model
            } else if (x$type == "SIR") {
                dat <- data.frame(x$XYcoordinates, x$inftime, x$remtime)
                if (is.null(time_id)) {
					ntimes <- max(x$inftime) - tmin + 1
					mfrow1 <- switch(min(ntimes,13), c(1,1), c(1,2), c(2,2), c(2,2), c(3,2), c(3,2), 
					c(3,3), c(3,3), c(3,3), c(3,2), c(3,2), c(3,2), c(3,3))
					
					sepwindow <- seq(prod(mfrow1), prod(mfrow1)*ceiling(ntimes/prod(mfrow1)), 
					by = prod(mfrow1))
					
					op <- par(mar = c(5.9, 4.0, 1.5, 0.5), omi = c(0.2, 0.05, 0.15, 0.15), 
						      mfrow = mfrow1)
					u <- 1
					
                    for(i in tmin:max(x$inftime)) {
                        xcc <- subset(dat, dat$x.inftime <= i & dat$x.inftime != 0)
                        xx    = x$XYcoordinates[,1]
                        yy    = x$XYcoordinates[,2]
                        plot(xx, yy,
                        xlim = c(min(xx), max(xx)),
                        ylim = c(min(yy), max(yy)),
                        pch  = 21,
                        main  = paste("time ", i), cex = 1)#, ...)
                        xred <- subset(xcc, i < xcc[,4])
                        points(xred[,1], xred[,2], pch = 16, col = "red")
                        xblue <- subset(xcc, i >= xcc[,4])
                        points(xblue[,1], xblue[,2], pch = 16, col = "blue")

                        if (any(u == sepwindow) | (u == ntimes)) {
                            opar <- par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0),
                            mar = c(0, 0, 0, 0), new = TRUE)
                            on.exit(par(opar), add = TRUE)
                            plot(0, 0, type = 'n', bty = 'n', xaxt = 'n', yaxt = 'n')
                            legend("bottom", c("Infected", "Susceptible", "Removed"), 
                            col = c("red", "black", "blue"), pch = c(16, 21, 16), 
                            bty = "n", horiz = TRUE, cex = 1.5)
						    op <- par(mar = c(5.9, 4.0, 1.5, 0.5), omi = c(0.2, 0.05, 0.15, 0.15), 
						    mfrow = mfrow1)
						   	on.exit(par(op), add = TRUE)
						}
						u <- u + 1
                    }
                } else if (!is.null(time_id)) {
					ntimes <- length(time_id)
					mfrow1 <- switch(min(ntimes,13), c(1,1), c(1,2), c(2,2), c(2,2), c(3,2), c(3,2), 
					c(3,3), c(3,3), c(3,3), c(3,2), c(3,2), c(3,2), c(3,3))
					
					sepwindow <- seq(prod(mfrow1), prod(mfrow1)*ceiling(length(time_id)/prod(mfrow1)), 
					by = prod(mfrow1))
										
					op <- par(mar = c(5.9, 4.0, 1.5, 0.5), omi = c(0.2, 0.05, 0.15, 0.15), 
						      mfrow = mfrow1)

                    for (i in 1:length(time_id)) {
                        xcc <- subset(dat, dat$x.inftime <= time_id[i] & dat$x.inftime != 0)
                        xx    = x$XYcoordinates[,1]
                        yy    = x$XYcoordinates[,2]
                        plot(xx, yy,
                        xlim = c(min(xx), max(xx)),
                        ylim = c(min(yy), max(yy)),
                        pch  = 21,
                        main  = paste("time ", time_id[i]), cex = 1, ...)
                        xred <- subset(xcc, time_id[i] < xcc[,4])
                        points(xred[,1], xred[,2], pch = 16, col = "red")
                        xblue <- subset(xcc, time_id[i] >= xcc[,4])
                        points(xblue[,1], xblue[,2], pch = 16, col = "blue")
                        if (any(i == sepwindow) | (i == length(time_id))) {
                            opar <- par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0),
                            mar = c(0, 0, 0, 0), new = TRUE)
                            on.exit(par(opar), add = TRUE)
                            plot(0, 0, type = 'n', bty = 'n', xaxt = 'n', yaxt = 'n')
	                        legend(-.03, 1.15, c("Infected", "Susceptible", "Removed"), 
	                        col = c("red", "black", "blue"), pch = c(16, 21, 16), 
	                        bty = "n", horiz = TRUE, cex = 1)
							op <- par(mar = c(5.9, 4.0, 1.5, 0.5), omi = c(0.2, 0.05, 0.15, 0.15), mfrow = mfrow1)
						   	on.exit(par(op), add = TRUE)
						}
                    }
                }
            # end if - SIR model
            }

#            par(old.par)

			on.exit(par(op1))

    # End function
        } else {
            stop("The plottype option should be either \"curve\" or \"spatial\"", call. = FALSE)
        }
    }
    # End of function
}
