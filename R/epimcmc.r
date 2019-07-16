################################################################################
# Part of the R/EpiILM package
#
# AUTHORS:
#         Vineetha Warriyar. K. V. <vineethawarriyar.kod@ucalgary.ca> and
#         Rob Deardon <robert.deardon@ucalgary.ca>
#
# Based on:
#                Deardon R, Brooks, S. P., Grenfell, B. T., Keeling, M. J., Tildesley,
#                M. J., Savill, N. J., Shaw, D. J.,  Woolhouse, M. E. (2010).
#                Inference for individual level models of infectious diseases in large
#                populations. Statistica Sinica, 20, 239-261.
#
#
# Free software under the terms of the GNU General Public License, version 2,
# a copy of which is available at http://www.r-project.org/Licenses/.
################################################################################

epimcmc <- function (type, x = NULL, y = NULL, inftime, tmin = NULL, tmax,
               infperiod = NULL, niter, alphaini, betaini, sparkini = NULL,
               Sformula = NULL, contact = NULL, pro.var.a, pro.var.b, pro.var.sp = NULL,
               prioralpha, halfnorm.var.a = NULL, gamma.par.a = NULL, unif.range.a = NULL,
               priorbeta, halfnorm.var.b = NULL, gamma.par.b = NULL, unif.range.b = NULL,
               priorsp = NULL, halfnorm.var.sp = NULL, gamma.par.sp = NULL, unif.range.sp = NULL,
               tempseed = NULL) {
  
  
  if (is.null(tempseed)) {
      tempseed <- 0
  }
 
 # error checks for input arguments
  if (is.null(type) || !(type %in% c("SI", "SIR"))) {
    stop("epimcmc: Specify type as \"SI\" or \"SIR\" ", call. = FALSE)
  }
  
  ns <- length(alphaini)
  ni <- length(betaini)
  
  if (!is.vector(inftime)) {
    stop('epimcmc: inftime is not a vector')
  }
  
  n <- length(inftime)
  
  if (is.null(infperiod) && type == "SIR") {
    stop(' epimcmc: Specify removal distance, infperiod ')
  }
  
  if (!is.null(infperiod)) {
      if (length(infperiod) != n) {
        stop('epimcmc: Length of infperiod is not compatible')
      }
      if (type == "SI") {
        stop('epimcmc: Type must be "SIR"')
      }
  } else {
   infperiod <- rep(0, n)
  }
  
  if (is.null(contact) &  (is.null(x) || is.null(y))) {
      stop('epimcmc: Specify contact network or x, y coordinates')
  }
  
  if (!is.null(x)) {
      if ((length(y) != n) || (length(x) != n)) {
        stop('epimcmc: Length of x or y is not compatible ')
      }
  }
  if (is.null(tmin)) {
    tmin <- 1
  }
  if (!is.null(sparkini)) {
    flag <- 1
    if (is.null(pro.var.sp)) {
    stop('epimcmc: Specify proposal variance for spark')
    }
    if (is.null(priorsp)) {
        stop('epimcmc: Specify prior for spark')
    }
  } else {
    sparkini <- 0
    flag     <- 0
  }

# formula for susceptibility function
  if (!is.null(Sformula)) {
    covmat <- model.matrix(Sformula)
    
    if ((ncol(covmat) == length(all.vars(Sformula))) & (ns != length(all.vars(Sformula)))) {
      stop('epimcmc: Check Sformula (no intercept term) and the dimension of alpha')
    }
    
    if ((ncol(covmat) > length(all.vars(Sformula))) & (ns != ncol(covmat))) {
      stop('epimcmc: Check Sformula (intercept term) and the dimension of alpha')
    }
  } else {
    if (ns == 1) {
      covmat <- matrix(1.0, nrow = n, ncol = ns)
    }
  }

# likelihood selection for SI and SIR
  if (type == "SI") {
      tnum <- 1
  }
  if (type == "SIR") {
      tnum <- 2
  }

# input prior - error check
  if (is.null(prioralpha) || !(prioralpha %in% c("halfnormal", "gamma","uniform"))) {
    stop("epimcmc: Specify prior for alpha as \"halfnormal\" ,  \"gamma\" or \"uniform\"  ", call. = FALSE)
  }

  if (is.null(priorbeta) || !(priorbeta %in% c("halfnormal", "gamma","uniform"))) {
    stop("epimcmc: Specify prior for beta as \"halfnormal\" ,  \"gamma\" or \"uniform\"  ", call. = FALSE)
  }
  if (!is.null(priorsp)) {
    if (!(priorsp %in% c("halfnormal", "gamma","uniform"))) {
        stop("epimcmc: Specify prior for spark term as \"halfnormal\" ,  \"gamma\" or \"uniform\"  ", call. = FALSE)
    }
  }

  # initialization
  unifmin <- vector(mode = "numeric", length = ns)
  unifmax <- vector(mode = "numeric", length = ns)
  gshape  <- vector(mode = "numeric", length = ns)
  gscale  <- vector(mode = "numeric", length = ns)
  halfvar <- vector(mode = "numeric", length = ns)

# prior selection for alpha
  if (prioralpha == "gamma") {
    anum <- 1
    if (ns == 1) {
        if (is.null(gamma.par.a)) {
            stop('epimcmc: Specify gamma.par.a')
        } else {
            gshape <- gamma.par.a[1]
            gscale <- gamma.par.a[2]
        }
    }
    if (ns > 1) {
        if (is.null(gamma.par.a) ||  !is.matrix(gamma.par.a)) {
            stop('epimcmc: Specify gamma.par.a as a matrix with each row corresponds to each alpha')
        } else {
            for (i in 1:ns) {
                gshape[i] <- gamma.par.a[i,1]
               gscale[i] <- gamma.par.a[i,2]
            }
        }
    }
  # end if - gamma prior
  }
  
  if (prioralpha == "halfnormal") {
    anum <- 2
    if (ns == 1) {
        if (is.null(halfnorm.var.a)) {
            stop('epimcmc: Specify halfnorm.var.a ')
        } else {
            halfvar <- halfnorm.var.a
        }
    }
    if (ns > 1) {
        if (is.null(halfnorm.var.a) || length(halfnorm.var.a) != ns) {
            stop('epimcmc: Specify halfnorm.var.a as a vector' )
        } else {
          for (i in 1:ns) {
            halfvar[i] <- halfnorm.var.a[i]
          }
        }
    }
  # end if - half normal prior
  }

if (prioralpha == "uniform") {
    anum <- 3
    if (ns == 1) {
        if (is.null(unif.range.a)) {
            stop('epimcmc: Specify unif.range.a ')
        } else {
           unifmin <- unif.range.a[1]
           unifmax <- unif.range.a[2]
        }
    }
    if (ns > 1) {
        if (is.null(unif.range.a) || !is.matrix(unif.range.a)) {
            stop('epimcmc: Specify unif.range.a as a matrix with each row corresponds to each alpha')
        } else {
            for (i in 1:ns) {
                unifmin[i] <- unif.range.a[i,1]
                unifmax[i] <- unif.range.a[i,2]
            }
        }
    }
  # end if - uniform prior
  }

  # initialization
  unifminb <- vector(mode = "numeric", length = ni)
  unifmaxb <- vector(mode = "numeric", length = ni)
  gshapeb  <- vector(mode = "numeric", length = ni)
  gscaleb  <- vector(mode = "numeric", length = ni)
  halfvarb <- vector(mode = "numeric", length = ni)

# prior selection for beta
if (priorbeta == "gamma") {
    bnum <- 1
    if (ni == 1) {
        if (is.null(gamma.par.b)) {
            stop('epimcmc: Specify gamma.par.b ')
        } else {
            gshapeb <- gamma.par.b[1]
            gscaleb <- gamma.par.b[2]
        }
    }
    if (ni > 1) {
        if (is.null(gamma.par.b) || !is.matrix(gamma.par.b)) {
            stop('epimcmc: Specify gamma.par.b as a matrix with each row corresponds to each beta')
        } else {
            for (i in 1:ni) {
                gshapeb[i] <- gamma.par.b[i,1]
                gscaleb[i] <- gamma.par.b[i,2]
            }
        }
    }
  # end if - gamma prior
  }

if  (priorbeta == "halfnormal") {
    bnum <- 2
    if (ni == 1) {
        if (is.null(halfnorm.var.b)) {
            stop('epimcmc: Specify halfnorm.var.b ')
        } else {
            halfvarb <- halfnorm.var.b
        }
    }
    if (ni > 1) {
        if (is.null(halfnorm.var.b) || length(halfnorm.var.b) != ni) {
            stop('epimcmc: Specify halfnorm.var.b as a vector' )
        } else {
            for (i in 1:ni) {
                halfvarb[i] <- halfnorm.var.b[i]
            }
        }
    }
  # end id - half normal prior
  }

if (priorbeta == "uniform") {
    bnum <- 3
    if (ni == 1) {
        if (is.null(unif.range.b)) {
            stop('epimcmc: Specify unif.range.b ')
        } else {
            unifminb <- unif.range.b[1]
            unifmaxb <- unif.range.b[2]
        }
    }
    if (ni > 1) {
        if (is.null(unif.range.b) || !is.matrix(unif.range.b)) {
            stop('epimcmc: Specify unif.range.b as a matrix with each row corresponds to each beta')
        } else {
            for (i in 1:ni) {
                unifminb[i] <- unif.range.b[i,1]
                unifmaxb[i] <- unif.range.b[i,2]
            }
        }
    }
  # end if - uniform prior
  }

  # initializations
  unifminsp <- 0.0
  unifmaxsp <- 0.0
  gshapesp  <- 0.0
  gscalesp  <- 0.0
  halfvarsp <- 0.0

# Prior selection for spark if it exist in the model
  if (!is.null(priorsp)) {
    if (priorsp == "gamma") {
     snum <- 1
         if (is.null(gamma.par.sp)) {
             stop('epimcmc: Specify gamma.par.sp ')
         } else {
             gshapesp <- gamma.par.sp[1]
             gscalesp <- gamma.par.sp[2]
         }
    }
   if (priorsp == "halfnormal") {
     snum <- 2
     if (is.null(halfnorm.var.sp)) {
         stop('epimcmc: Specify halfnorm.var.sp ')
     } else {
        halfvarsp <- halfnorm.var.sp
     }
   }
   if (priorsp == "uniform") {
     snum <- 3
     if (is.null(unif.range.sp)) {
         stop('epimcmc: Specify unif.range.sp ')
     } else {
         unifminsp <- unif.range.sp[1]
         unifmaxsp <- unif.range.sp[2]
     }
   }
  } else {
     snum <- 0
  }

# proposal variance for alpha
  prostda <- vector(mode = "numeric", length = ns)
  if (ns == 1) {
    prostda <- sqrt(pro.var.a)
  }
  if (ns > 1) {
      if (length(pro.var.a) != ns) {
        stop('epimcmc: Specify proposal variance for each alpha parameter')
      }
    for (i in 1:ns) {
            prostda[i] <- sqrt(pro.var.a[i])
        }
  }

# proposal variance for  beta
  prostdb <- vector(mode = "numeric", length = ni)
  if (ni == 1) {
    prostdb <- sqrt(pro.var.b)
  }
  if (ni > 1) {
      if (length(pro.var.b) != ni) {
        stop('epimcmc: Specify proposal variance for each beta parameter')
      }
    for (i in 1:ni) {
            prostdb[i] <- sqrt(pro.var.b[i])
        }
  }

# proposal variance for  spark
  if (!is.null(pro.var.sp)) {
    prostdsp <- sqrt(pro.var.sp)
  } else {
     prostdsp <- 0.0
  }
 
 # Calling fortran subroutine for Purely Spatial models
  if (is.null(contact)) {
    tmp <- .Fortran("mcmc",
                     tnum = as.integer(tnum),
                     x = as.numeric(x),
                     y = as.numeric(y),
                     tau = as.integer(inftime),
                     n = as.integer(n),
                     lambda = as.integer(infperiod),
                     tmin = as.integer(tmin),
                     tmax = as.integer(tmax),
                     nsim = as.integer(niter),
                     aalpha = as.numeric(alphaini),
                     ns = as.integer(ns),
                     ni = as.integer(ni),
                     bbeta = as.double(betaini),
                     covmat = as.vector(covmat),
                     prostda = as.double(prostda),
                     prostdb = as.double(prostdb),
                     anum = as.integer(anum),
                     bnum = as.integer(bnum),
                     halfvar = as.double(halfvar),
                     unifmin = as.double(unifmin),
                     unifmax = as.double(unifmax),
                     gshape = as.double(gshape),
                     gscale = as.double(gscale),
                     halfvarb = as.double(halfvarb),
                     unifminb = as.double(unifminb),
                     unifmaxb = as.double(unifmaxb),
                     gshapeb = as.double(gshapeb),
                     gscaleb = as.double(gscaleb),
                     simalpha = matrix(0, nrow = niter, ncol = ns),
                     simbeta = matrix(0, nrow = niter, ncol = ni),
                     sspark = as.numeric(sparkini),
                     flag = as.integer(flag),
                     prostdsp = as.double(prostdsp),
                     snum = as.integer(snum),
                     halfvarsp = as.double(halfvarsp),
                     unifminsp = as.double(unifminsp),
                     unifmaxsp = as.double(unifmaxsp),
                     gshapesp = as.double(gshapesp),
                     gscalesp = as.double(gscalesp),
                     simspark = as.double(rep(0,niter)),
                     llikeval = as.double(rep(0,niter)),
                     tempseed = as.integer(tempseed))
                     
    if (flag == 0) {
        result        <- data.frame(ALPHA = tmp$simalpha, BETA = tmp$simbeta)
        Loglikelihood <- tmp$llikeval
    } else {
      result        <- data.frame(ALPHA = tmp$simalpha, BETA = tmp$simbeta, SPARK = tmp$simspark)
      Loglikelihood <- tmp$llikeval
    }
  }

# Calling fortran subroutines for contact network models
  if (!is.null(contact)) {
      if (length(contact)/(n * n) != ni) {
        stop('epimcmc:  Dimension of beta  and the number of contact networks are not matching')
      }
    network <- array(contact, c(n, n, ni))
    
    tmp <- .Fortran("conmcmc",
                    tnum = as.integer(tnum),
                    tau = as.integer(inftime),
                    n = as.integer(n),
                    lambda = as.integer(infperiod),
                    tmin = as.integer(tmin),
                    tmax = as.integer(tmax),
                    nsim = as.integer(niter),
                    aalpha = as.numeric(alphaini),
                    ns = as.integer(ns),
                    ni = as.integer(ni),
                    bbeta = as.double(betaini),
                    covmat = as.vector(covmat),
                    network = as.vector(network),
                    prostda = as.double(prostda),
                    prostdb = as.double(prostdb),
                    anum = as.integer(anum),
                    bnum = as.integer(bnum),
                    halfvar = as.double(halfvar),
                    unifmin = as.double(unifmin),
                    unifmax = as.double(unifmax),
                    gshape = as.double(gshape),
                    gscale = as.double(gscale),
                    halfvarb = as.double(halfvarb),
                    unifminb = as.double(unifminb),
                    unifmaxb = as.double(unifmaxb),
                    gshapeb = as.double(gshapeb),
                    gscaleb = as.double(gscaleb),
                    simalpha = matrix(0, nrow = niter, ncol = ns),
                    simbeta = matrix(0, nrow = niter, ncol = ni),
                    sspark = as.numeric(sparkini),
                    flag = as.integer(flag),
                    prostdsp = as.double(prostdsp),
                    snum = as.integer(snum),
                    halfvarsp = as.double(halfvarsp),
                    unifminsp = as.double(unifminsp),
                    unifmaxsp = as.double(unifmaxsp),
                    gshapesp = as.double(gshapesp),
                    gscalesp = as.double(gscalesp),
                    simspark = as.double(rep(0, niter)),
                    llikeval = as.double(rep(0, niter)),
                    tempseed = as.integer(tempseed))
                    
    if (flag == 0) {
        result       <- data.frame(ALPHA = tmp$simalpha, BETA = tmp$simbeta)
       Loglikelihood <- tmp$llikeval
    } else {
        result        <- data.frame(ALPHA = tmp$simalpha, BETA = tmp$simbeta, SPARK = tmp$simspark)
        Loglikelihood <- tmp$llikeval
    }
  }
  # mcmc result as a coda object
  result <- coda::mcmc(result)
  list(Estimates = result, Loglikelihood = Loglikelihood)


# End of function
}










