################################################################################
# Part of the R/EpiILM package
#
# AUTHORS:
#         Waleed Almutiry <wkmtierie@qu.edu.sa>,
#         Vineetha Warriyar. K. V. <vineethawarriyar.kod@ucalgary.ca>, and
#         Rob Deardon <robert.deardon@ucalgary.ca>
#
# Algorithm based on:
#                Deardon R, Brooks, S. P., Grenfell, B. T., Keeling, M. J., Tildesley,
#                M. J., Savill, N. J., Shaw, D. J.,  Woolhouse, M. E. (2010).
#                Inference for individual level models of infectious diseases in large
#                populations. Statistica Sinica, 20, 239-261.
#
# Free software under the terms of the GNU General Public License, version 2,
# a copy of which is available at http://www.r-project.org/Licenses/.
################################################################################

pred.epi <- function (object, xx, criterion , n.samples, burnin = NULL, tmin = NULL, Sformula = NULL, Tformula = NULL,  showProgressBar = interactive()) {

	if (class(object) != "epidata") {
		stop("The object of the epidemic data must be in a class of \"epidata\"", call. = FALSE)
	}

	if (class(xx) != "epimcmc") {
		stop("The object of the MCMC samples must be in a class of \"epimcmc\"", call. = FALSE)
	}

	n <- length(object$inftime)

	if (is.null(burnin)) {
		burnin <- 1
	}

	if (is.null(tmin)){
		tmin = 1
	}

	t_end <- max(object$inftime)

	ns <- xx$n.sus.par
	nt <- xx$n.trans.par
	ni <- xx$n.ker.par

    #	if (is.null(tempseed)) {

    #	temp <- 0

    #} else {

    #	temp <- tempseed
    #	set.seed(temp)

    #}

    if (object$type == "SIR") {
        infperiod = object$remtime - object$inftime
    } else {
        infperiod = NULL
    }


	mcmcout_posterior <- matrix(xx$Fullsample[sample(seq(burnin,length(xx$Fullsample[,1])), n.samples, replace = FALSE), ],
                            nrow = n.samples,
                            ncol = ncol(xx$Fullsample))

														xx$Fullsample[sample(seq(burnin,length(xx$Fullsample[,1])), n.samples, replace = FALSE), ]
	newinftime <- replace(object$inftime, object$inftime > tmin, 0)
	out1 <- vector(mode = "list", length = n.samples)

	if (is.null(object$contact)) {

		x <- object$XYcoordinates[,1]
		y <- object$XYcoordinates[,2]

		if ( (nt == 0) & (ns + ni) == ncol(xx$Fullsample)) {

			cat("generate", n.samples, "epidemics \n")
			if (showProgressBar) {
	    pb <- txtProgressBar(min = 0, max = n.samples, style = 3)
			}
			update.step <- max(5, floor(n.samples/100))

			for (i in 1:n.samples) {
				if (showProgressBar && i %% update.step == 0) {
		      setTxtProgressBar(pb, i)
		    }

				out1[[i]] <- epidata(type = object$type, n = n, tmax = t_end, Sformula = Sformula, tmin = tmin,
				sus.par = c(mcmcout_posterior[i, 1: ns]),
				beta = c(mcmcout_posterior[i, (ns + 1): (ns + ni)]),
				x = x, y = y, inftime = newinftime, infperiod = infperiod)

			}
			if (showProgressBar) {
				close(pb)
			}
		} else if ( (nt > 0) & (ns + nt + ni) == ncol(xx$Fullsample)) {

			cat("  generate", n.samples, "epidemics \n")
			if (showProgressBar) {
			pb <- txtProgressBar(min = 0, max = n.samples, style = 3)
			}
			update.step <- max(5, floor(n.samples/100))

			for (i in 1:n.samples) {

				if (showProgressBar && i %% update.step == 0) {
					setTxtProgressBar(pb, i)
				}

				out1[[i]] <- epidata(type = object$type, n = n, tmax = t_end, Sformula = Sformula, Tformula = Tformula, tmin = tmin,
				sus.par = c(mcmcout_posterior[i, 1:ns]), trans.par = c(mcmcout_posterior[i, (ns+1): (ns+nt)]),
				beta = c(mcmcout_posterior[i, (ns+nt+1): (ns+nt+ni)]),
				x = x, y = y, inftime = newinftime, infperiod = infperiod)

			}
			if (showProgressBar) {
				close(pb)
			}

		} else if ( (nt == 0) & (ns + ni) < ncol(xx$Fullsample)) {

			cat("  generate", n.samples, "epidemics \n")
			if (showProgressBar) {
			pb <- txtProgressBar(min = 0, max = n.samples, style = 3)
			}
			update.step <- max(5, floor(n.samples/100))

			for (i in 1:n.samples) {

				if (showProgressBar && i %% update.step == 0) {
					setTxtProgressBar(pb, i)
				}

				out1[[i]] <- epidata(type = object$type, n = n, tmax = t_end, Sformula = Sformula, tmin = tmin,
				sus.par = c(mcmcout_posterior[i, 1: ns]), trans.par = 1,
				beta = c(mcmcout_posterior[i, (ns + 1): (ns + ni)]),
				spark = mcmcout_posterior[i, (ns + ni + 1)],
				x = x, y = y, inftime = newinftime, infperiod = infperiod)
               
			}
			if (showProgressBar) {
				close(pb)
			}

		} else if ( (nt > 0) & (ns + nt + ni) < ncol(xx$Fullsample)) {

			cat("  generate", n.samples, "epidemics \n")
			if (showProgressBar) {
			pb <- txtProgressBar(min = 0, max = n.samples, style = 3)
			}
			update.step <- max(5, floor(n.samples/100))

			for (i in 1:n.samples) {

				if (showProgressBar && i %% update.step == 0) {
					setTxtProgressBar(pb, i)
				}

				out1[[i]] <- epidata(type = object$type, n = n, tmax = t_end, Sformula = Sformula, Tformula = Tformula, tmin = tmin,
				sus.par = c(mcmcout_posterior[i, 1:ns]), trans.par = c(mcmcout_posterior[i, (ns+1): (ns+nt)]),
				beta = c(mcmcout_posterior[i, (ns+nt+1): (ns+nt+ni)]),
				spark = c(mcmcout_posterior[i, (ns+nt+ni+1)]),
				x = x, y = y, inftime = newinftime, infperiod = infperiod)
								

			}
			if (showProgressBar) {
				close(pb)
			}
		}

	} else {

		contact <- object$contact

		if ( (nt == 0) & (ns + ni) == ncol(xx$Fullsample)) {

			cat("  generate", n.samples, "epidemics \n")
			if (showProgressBar) {
			pb <- txtProgressBar(min = 0, max = n.samples, style = 3)
			}
			update.step <- max(5, floor(n.samples/100))

			for (i in 1:n.samples) {

				if (showProgressBar && i %% update.step == 0) {
					setTxtProgressBar(pb, i)
				}

        if (ni != 0) {
  				out1[[i]] <- epidata(type = object$type, n = n, tmax = t_end, Sformula = Sformula, tmin = tmin,
    				sus.par = c(mcmcout_posterior[i, 1: ns]),
    				beta = c(mcmcout_posterior[i, (ns + 1): (ns + ni)]),
    				contact = contact, inftime = newinftime, infperiod = infperiod)
          
        } else if (ni == 0) {
  				out1[[i]] <- epidata(type = object$type, n = n, tmax = t_end, Sformula = Sformula, tmin = tmin,
    				sus.par = c(mcmcout_posterior[i, 1: ns]),
    				contact = contact, inftime = newinftime, infperiod = infperiod)
            
        }
			}
			if (showProgressBar) {
				close(pb)
			}

		} else if ( (nt > 0) & (ns + nt + ni) == ncol(xx$Fullsample)) {

			cat("  generate", n.samples, "epidemics \n")
			if (showProgressBar) {
			pb <- txtProgressBar(min = 0, max = n.samples, style = 3)
			}
			update.step <- max(5, floor(n.samples/100))

			for (i in 1:n.samples) {

				if (showProgressBar && i %% update.step == 0) {
					setTxtProgressBar(pb, i)
				}

        if (ni != 0) {
  				out1[[i]] <- epidata(type = object$type, n = n, tmax = t_end, Sformula = Sformula, Tformula = Tformula, tmin = tmin,
    				sus.par = c(mcmcout_posterior[i, 1:ns]), trans.par = c(mcmcout_posterior[i, (ns+1): (ns+nt)]),
    				beta = c(mcmcout_posterior[i, (ns+nt+1): (ns+nt+ni)]),
    				contact = contact, inftime = newinftime, infperiod = infperiod)
            
        } else if (ni == 0) {
          out1[[i]] <- epidata(type = object$type, n = n, tmax = t_end, Sformula = Sformula, Tformula = Tformula, tmin = tmin,
            sus.par = c(mcmcout_posterior[i, 1:ns]), trans.par = c(mcmcout_posterior[i, (ns+1): (ns+nt)]),
            contact = contact, inftime = newinftime, infperiod = infperiod)
            
        }
			}
			if (showProgressBar) {
				close(pb)
			}

		} else if ( (nt == 0) & (ns + ni) < ncol(xx$Fullsample)) {

			cat("  generate", n.samples, "epidemics \n")
			if (showProgressBar) {
			pb <- txtProgressBar(min = 0, max = n.samples, style = 3)
			}
			update.step <- max(5, floor(n.samples/100))

			for (i in 1:n.samples) {

				if (showProgressBar && i %% update.step == 0) {
					setTxtProgressBar(pb, i)
				}

        if (ni != 0) {
  				out1[[i]] <- epidata(type = object$type, n = n, tmax = t_end, Sformula = Sformula, tmin = tmin,
    				sus.par = c(mcmcout_posterior[i, 1: ns]),
    				beta = c(mcmcout_posterior[i, (ns + 1): (ns + ni)]),
    				spark = mcmcout_posterior[i, (ns + ni + 1)],
    				contact = contact, inftime = newinftime, infperiod = infperiod)
           
        } else if (ni == 0) {
          out1[[i]] <- epidata(type = object$type, n = n, tmax = t_end, Sformula = Sformula, tmin = tmin,
            sus.par = c(mcmcout_posterior[i, 1: ns]),
            spark = mcmcout_posterior[i, (ns + ni + 1)],
            contact = contact, inftime = newinftime, infperiod = infperiod)
         
        }
			}

			if (showProgressBar) {
				close(pb)
			}

		} else if ( (nt > 0) & (ns + nt + ni) < ncol(xx$Fullsample)) {

			cat("  generate", n.samples, "epidemics \n")
			if (showProgressBar) {
			pb <- txtProgressBar(min = 0, max = n.samples, style = 3)
			}
			update.step <- max(5, floor(n.samples/100))

			for (i in 1:n.samples) {

				if (showProgressBar && i %% update.step == 0) {
					setTxtProgressBar(pb, i)
				}

        if (ni != 0) {
          out1[[i]] <- epidata(type = object$type, n = n, tmax = t_end, Sformula = Sformula, Tformula = Tformula, tmin = tmin,
            sus.par = c(mcmcout_posterior[i, 1:ns]), trans.par = c(mcmcout_posterior[i, (ns+1): (ns+nt)]),
            beta = c(mcmcout_posterior[i, (ns+nt+1): (ns+nt+ni)]),
            spark = c(mcmcout_posterior[i, (ns+nt+ni+1)]),
            contact = contact, inftime = newinftime, infperiod = infperiod)
            
        } else if (ni == 0) {
          out1[[i]] <- epidata(type = object$type, n = n, tmax = t_end, Sformula = Sformula, Tformula = Tformula, tmin = tmin,
            sus.par = c(mcmcout_posterior[i, 1:ns]), trans.par = c(mcmcout_posterior[i, (ns+1): (ns+nt)]),
            spark = c(mcmcout_posterior[i, (ns+nt+ni+1)]),
            contact = contact, inftime = newinftime, infperiod = infperiod)
           
        }
			}
			if (showProgressBar) {
				close(pb)
			}
		}

	}

	if (criterion == "newly infectious") {

		crit <- matrix(0, ncol = t_end, nrow = n.samples)

		for (j in 1:n.samples) {
			newinf <- rep(0)
			for (i in 1: t_end) {
				newinf[i] <- sum(out1[[j]]$inftime == i)
			}
			crit[j,] <- newinf
		}



		crit_true <- rep(0)
		for (i in 1: t_end) {
			crit_true[i] <- length(object$inftime[object$inftime == i])
		}

		result <- list(type = object$type, criterion = criterion, crit.sim = crit, crit.obs = crit_true, tmin = tmin, tmax = t_end, n.samples = n.samples)

	} else if (criterion == "epidemic length") {


		if (object$type == "SI") {

			crit <- NULL

			for (j in 1:n.samples) {
				crit[j] <- max(out1[[j]]$inftime) - min(out1[[j]]$inftime)
			}

			crit_true <- max(object$inftime) - min(object$inftime)

		} else if (object$type == "SIR") {

			crit <- NULL

			for (j in 1:n.samples) {
				crit[j] <- max(out1[[j]]$remtime) - min(out1[[j]]$inftime)
			}

			crit_true <- max(object$remtime) - min(object$inftime)

		}
			result <- list(type = object$type, criterion = criterion, crit.sim = crit, crit.obs = crit_true, tmin = tmin, tmax = t_end, n.samples = n.samples)


			} else if (criterion == "peak time") {

		    if (object$type == "SI") {

		      crit <- NULL

		      for (j in 1:n.samples) {

		        oot <- out1[[j]]$inftime[which(out1[[j]]$inftime!=0)]

		        dff <- as.data.frame(table(oot))

		        res <- dff[which(dff[, 2] == max(dff[, 2])), 1]

		        if (length(res) > 1) {

		          crit[j] <- max(as.numeric(as.character(res)))

		        } else {

		          crit[j] <- as.numeric(as.character(res))

		        }

		      }

		      oot <- object$inftime[which(object$inftime!=0)]

		      dff <- as.data.frame(table(oot))

		      crit_true <- max(as.numeric(as.character(dff[which(dff[, 2] == max(dff[, 2])), 1])))

		      result <- list(type = object$type, criterion = criterion, crit.sim = crit, crit.obs = crit_true, tmin = tmin, tmax = t_end, n.samples = n.samples)

		    } else if (object$type == "SIR") {

		      crit_inf <- NULL

		      crit_rem <- NULL

		      for (j in 1:n.samples) {

		        oot <- out1[[j]]$inftime[which(out1[[j]]$inftime!=0)]

		        oot1 <- out1[[j]]$remtime[which(out1[[j]]$remtime!=0)]

		        dff <- as.data.frame(table(oot))

		        dff1 <- as.data.frame(table(oot1))

		        res <- dff[which(dff[, 2] == max(dff[, 2])), 1]

		        res1<- dff1[which(dff1[, 2] == max(dff1[, 2])), 1]

		        if (length(res) > 1) {

		          crit_inf[j] <- max(as.numeric(as.character(res)))

		        } else {

		          crit_inf[j] <- as.numeric(as.character(res))

		        }

		        if (length(res1) > 1) {

		          crit_rem[j] <- max(as.numeric(as.character(res1)))

		        } else {

		          crit_rem[j] <- as.numeric(as.character(res1))

		        }

		      }

		      oot <- object$inftime[which(object$inftime!=0)]

		      oot1 <- object$remtime[which(object$remtime!=0)]

		      dff <- as.data.frame(table(oot))

		      dff1 <- as.data.frame(table(oot1))

		      crit_true_inf <- max(as.numeric(as.character(dff[which(dff[, 2] == max(dff[, 2])), 1])))

		      crit_true_rem <- max(as.numeric(as.character(dff1[which(dff1[, 2] == max(dff1[, 2])), 1])))

		      result <- list(type = object$type, criterion = criterion, crit.sim = list(crit_inf, crit_rem),
		                     crit.obs = list(crit_true_inf, crit_true_rem), tmin = tmin, tmax = t_end, n.samples = n.samples)

		    }


		  }

		  class(result) <- "pred.epi"

		  result

}

plot.pred.epi <- function (x, ...) {

	if (class(x) != "pred.epi") {
		stop("The object must be in a class of \"pred.epi\"", call. = FALSE)
	}

	if (x$criterion == "newly infectious") {

		lowerq <- NULL
		upperq <- NULL

		for (i in 1: x$tmax) {
			lowerq[i] <- quantile(x$crit.sim[,i], 0.025)
			upperq[i] <- quantile(x$crit.sim[,i], 0.975)
		}

		time <- rep(1: x$tmax)
		mx <- max(upperq)
		mn <- min(lowerq)

		plot (time, x$crit.obs, xlim = c(min(time), max(time)), ylim = c(mn, mx),
			ylab = "Newly infections", xlab =" Time", xaxt = "n", col = "black", type = "o", pch = 20)

		axis(1, at = 1:max(time))

		for (i in 1: x$n.samples) {
            lines(time[x$tmin:x$tmax], x$crit.sim[i,x$tmin:x$tmax], col = "grey", lwd = 0.2)
		}

        lines(time[x$tmin:x$tmax], apply(x$crit.sim[,x$tmin:x$tmax], 2, mean), lty = 1, ...)
        lines(time, x$crit.obs, col = "black", type = "o", pch = 20, lwd = 2)

		axis(1, at = 1:max(time))
		lines(time[x$tmin:x$tmax], lowerq[x$tmin:x$tmax], lty = 2, ...)
		lines(time[x$tmin:x$tmax], upperq[x$tmin:x$tmax], lty = 2, ...)

	} else if (x$criterion == "epidemic length") {

		hist(x$crit.sim, main = "", xlab = "The length of epidemic", ylab = "Time", ...)
		abline(v = quantile(x$crit.sim, c(0.025, 0.975)), col = "red", lty = 2)
		abline(v = mean(x$crit.sim), col = "red")
		abline(v = x$crit.obs, col = "blue")

	} else if (x$criterion == "peak time") {

		if(x$type == "SI") {

			hist(x$crit.sim, main = "", xlab = "The time of the peak of infection times", ylab = "Time", ...)
			abline(v = quantile(x$crit.sim, c(0.025, 0.975)), col = "red", lty = 2)
			abline(v = mean(x$crit.sim), col = "red")
			abline(v = x$crit.obs, col = "blue")

		} else if (x$type == "SIR") {

      par(mfrow = c(1, 2))

			hist(x$crit.sim[[1]], main = "", xlab = "The time of the peak of infection times", ylab = "Time", ...)
			abline(v = quantile(x$crit.sim[[1]], c(0.025, 0.975)), col = "red", lty = 2)
			abline(v = mean(x$crit.sim[[1]]), col = "red")
			abline(v = x$crit.obs[[1]], col = "blue")

			hist(x$crit.sim[[2]], main = "", xlab = "The time of the peak of removal times", ylab = "Time", ...)
			abline(v = quantile(x$crit.sim[[2]], c(0.025, 0.975)), col = "red", lty = 2)
			abline(v = mean(x$crit.sim[[2]]), col = "red")
			abline(v = x$crit.obs[[2]], col = "blue")

		}
	}
}
