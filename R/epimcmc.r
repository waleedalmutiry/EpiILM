################################################################################
# Part of the R/EpiILM package
#
# AUTHORS:
#         Waleed Almutiry <wkmtierie@qu.edu.sa>,
#         Vineetha Warriyar. K. V. <vineethawarriyar.kod@ucalgary.ca> and
#         Rob Deardon <robert.deardon@ucalgary.ca>
#
# Based on:
#                Deardon R, Brooks, S. P., Grenfell, B. T., Keeling, M. J., Tildesley,
#                M. J., Savill, N. J., Shaw, D. J.,  Woolhouse, M. E. (2010).Inference
#                for individual level models of infectious diseases in large
#                populations. Statistica Sinica, 20, 239-261.
#
#
# Free software under the terms of the GNU General Public License, version 2,
# a copy of which is available at http://www.r-project.org/Licenses/.
################################################################################

epimcmc <- function (object, tmin = NULL, tmax, niter,
            sus.par.ini, trans.par.ini = NULL, beta.ini = NULL, spark.ini = NULL,
            Sformula = NULL, Tformula = NULL,
            pro.sus.var, pro.trans.var = NULL, pro.beta.var = NULL, pro.spark.var = NULL,
            prior.sus.dist, prior.trans.dist = NULL, prior.beta.dist = NULL, prior.spark.dist = NULL,
            prior.sus.par, prior.trans.par = NULL, prior.beta.par = NULL, prior.spark.par = NULL,
            adapt = FALSE, acc.rate = NULL) {

  if (!is(object, "epidata")) {
    	stop("The object must be in a class of \"epidata\"", call. = FALSE)
  } else {

    # error checks for input arguments

    if (any(is.null(object$type) | !(object$type %in% c("SI", "SIR"))) == TRUE) {
      stop("epimcmc: Specify type as \"SI\" or \"SIR\" ", call. = FALSE)
    }

    if (!is.vector(object$inftime)) {
      stop('epimcmc: inftime is not a vector')
    }

    n <- length(object$inftime)

    if (object$type == "SIR") {
      infperiod = object$remtime - object$inftime
    }

    if (all(is.null(object$contact) &  is.null(object$XYcoordinates)) == TRUE) {
      stop("epimcmc: Specify contact network or x, y coordinates")
    }

    if (is.null(tmin)) {
      tmin <- 1
    }


# susceptibility parameters:

    ns <- length(sus.par.ini)

    # formula for susceptibility function
    if (!is.null(Sformula)) {
      covmat.sus <- model.matrix(Sformula)
      if (all((ncol(covmat.sus) == length(all.vars(Sformula))) & (ns != length(all.vars(Sformula)))) == TRUE) {
        stop("epimcmc: Check Sformula (no intercept term) and the dimension of susceptibility parameters", call. = FALSE)
      } else if (all((ncol(covmat.sus) > length(all.vars(Sformula))) & (ns != ncol(covmat.sus))) == TRUE) {
        stop("epimcmc: Check Sformula (intercept term) and the dimension of susceptibility parameters", call. = FALSE)
      }
    } else {
      if (ns == 1) {
        covmat.sus <- matrix(1.0, nrow = n, ncol = ns)
      }
    }

    if (any(is.null(prior.sus.dist) | !(all(prior.sus.dist %in% c("halfnormal", "gamma","uniform")))) == TRUE) {
      stop("epimcmc: Specify prior for susceptibility parameters as \"halfnormal\" ,  \"gamma\" or \"uniform\"  ", call. = FALSE)
    }

    prostda <- vector(mode = "numeric", length = ns)
    if (ns == 1) {
      prostda <- pro.sus.var
    } else if (ns > 1) {
      if (length(pro.sus.var) != ns) {
        stop('epimcmc: Specify proposal variance for each susceptibility parameter')
      }
      for (i in 1:ns) {
        prostda[i] <- pro.sus.var[i]
      }
    }

    alphapar <- list(NULL)
    len.alphapar <- ns - sum(prostda == 0)
    if (len.alphapar != 0){
      j <- 1
      for (i in 1:ns) {
        if (prostda[i] != 0){
          alphapar[[j]] <- list(NULL)
          alphapar[[j]][[1]] <- prior.sus.dist[i]
          j <- j + 1
        } else {
          j <- j
        }
      }
    } else {
      alphapar[[1]] <- list(NULL)
      alphapar[[1]][[1]] <- "uniform" # just to fill the list and it is not being used at all.
    }

    if (all(ns == 1 & len.alphapar != 0) == TRUE) {
      if (any((prior.sus.dist == "gamma") | (prior.sus.dist == "uniform")) == TRUE) {
        if (is.null(prior.sus.par)) {
          stop('epimcmc: Specify prior.sus.par as a vector or a 1 by 2 matrix.')
        } else {
          alphapar[[1]][[2]] <- vector()
          alphapar[[1]][[2]][1] <- prior.sus.par[1]
          alphapar[[1]][[2]][2] <- prior.sus.par[2]
        }
      } else if (prior.sus.dist == "halfnormal") {
        if (is.null(prior.sus.par)) {
          stop('epimcmc: Specify the variance of the half normal prior distribution for prior.sus.par.')
        } else {
          alphapar[[1]][[2]] <- vector()
          alphapar[[1]][[2]][1] <- prior.sus.par[1]
          alphapar[[1]][[2]][2] <- 0
        }
      }
    } else if (all(ns == 1 & len.alphapar == 0) == TRUE) {
      alphapar[[1]][[2]] <- vector()
      alphapar[[1]][[2]][1] <- 0
      alphapar[[1]][[2]][2] <- 0
    } else if (all(ns > 1 & len.alphapar != 0) == TRUE) {
      if (any(is.null(prior.sus.par) |  !is.matrix(prior.sus.par)) == TRUE) {
        stop('epimcmc: Specify prior.sus.par as a matrix with each row corresponds to each susceptibility parameter')
      } else {
        j <- 1
        for (i in 1:ns) {
          if (prostda[i] != 0){
            alphapar[[j]][[2]] <- vector()
            alphapar[[j]][[2]][1] <- prior.sus.par[i,1]
            alphapar[[j]][[2]][2] <- prior.sus.par[i,2]
            j <- j + 1
          } else {
            j <- j
          }
        }
      }
    } else if (all(ns > 1 & len.alphapar == 0) == TRUE) {
      alphapar[[1]][[2]] <- vector()
      alphapar[[1]][[2]][1] <- 0
      alphapar[[1]][[2]][2] <- 0
    }# end if - ns

# transmissibility parameters:

    if (!is.null(trans.par.ini)) {
      nt <- length(trans.par.ini)
      flag.trans <- 1
    } else {
      nt <- 1
      flag.trans <- 0
      phipar <- 1
    }

    # formula for transmissibility function
    if (flag.trans == 1) {
      if (is.null(Tformula)) {
        stop("epimcmc: Tformula is missing. It has to be specified with no intercept term and number of columns equal to the length of trans.par", call. = FALSE)
      } else if (!is.null(Tformula)) {
        covmat.trans <- model.matrix(Tformula)
        if (all((ncol(covmat.trans) == length(all.vars(Tformula))) & (nt != length(all.vars(Tformula)))) == TRUE) {
          stop("epimcmc: Check Tformula. It has to be with no intercept term and number of columns equal to the length of trans.par", call. = FALSE)
        }
      }
    } else if (flag.trans == 0) {
      covmat.trans <- matrix(1.0, nrow = n, ncol = nt)
    }

    if (all(!is.null(prior.trans.dist) & !(all(prior.trans.dist %in% c("halfnormal", "gamma","uniform")))) == TRUE) {
      stop("epimcmc: Specify prior for transmissibility parameters as \"halfnormal\" ,  \"gamma\" or \"uniform\"  ", call. = FALSE)
    }

    if (all(is.null(prior.trans.dist) & flag.trans == 1) == TRUE) {
      stop("epimcmc: Specify prior for transmissibility parameters as \"halfnormal\" ,  \"gamma\" or \"uniform\"  ", call. = FALSE)
    }

    # proposal variance for transmissibility parameters
    prostdt <- vector(mode = "numeric", length = nt)
    if (all(flag.trans == 1 & nt == 1) == TRUE) {
      if (is.null(pro.trans.var)) {
        stop('epimcmc: Specify proposal variance for the transmissibility parameter')
      }
      prostdt <- pro.trans.var
    } else if (all(flag.trans == 1 & nt > 1) == TRUE) {
      if (is.null(pro.trans.var)) {
        stop('epimcmc: Specify proposal variance for each transmissibility parameter')
      }
      if (length(pro.trans.var) != nt) {
        stop('epimcmc: Specify proposal variance for each transmissibility parameter')
      }
      for (i in 1:nt) {
        prostdt[i] <- pro.trans.var[i]
      }
    } else {
      prostdt <- NULL
    }

    if (flag.trans == 1){
      phipar <- list(NULL)
      len.phipar <- nt - sum(prostdt == 0)
      if (len.phipar != 0){
        j <- 1
        for (i in 1:nt) {
          if (prostdt[i] != 0){
            phipar[[j]] <- list(NULL)
            phipar[[j]][[1]] <- prior.trans.dist[i]
            j <- j + 1
          } else {
            j <- j
          }
        }
      } else {
        phipar[[1]] <- list(NULL)
        phipar[[1]][[1]] <- "uniform" # just to fill the list and it is not being used at all.
      }
    } else {
      phipar <- 1
    }

    # prior selection for transmissibility parameters
    if (flag.trans == 1) {

      if (nt == 1 & len.phipar != 0) {
        if (any((prior.trans.dist == "gamma") | (prior.trans.dist == "uniform")) == TRUE) {
          if (is.null(prior.trans.par)) {
            stop('epimcmc: Specify prior.trans.par as a vector or a 1 by 2 matrix.')
          } else {
            phipar[[1]][[2]] <- vector()
            phipar[[1]][[2]][1] <- prior.trans.par[1]
            phipar[[1]][[2]][2] <- prior.trans.par[2]
          }
        } else if (prior.trans.dist == "halfnormal") {
          if (is.null(prior.trans.par)) {
            stop('epimcmc: Specify the variance of the half normal prior distribution for prior.trans.par.')
          } else {
            phipar[[1]][[2]] <- vector()
            phipar[[1]][[2]][1] <- prior.trans.par[1]
            phipar[[1]][[2]][2] <- 0
          }
        }
      } else if (all(nt == 1 & len.phipar == 0) == TRUE) {
        phipar[[1]][[2]] <- vector()
        phipar[[1]][[2]][1] <- 0
        phipar[[1]][[2]][2] <- 0
      } else if (all(nt > 1 & len.phipar != 0) == TRUE) {
        if (any(is.null(prior.trans.par) |  !is.matrix(prior.trans.par)) == TRUE) {
          stop('epimcmc: Specify prior.trans.par as a matrix with each row corresponds to each transmissibility parameter')
        } else {
          j <- 1
          for (i in 1:nt) {
            if (prostdt[i] != 0){
              phipar[[j]][[2]] <- vector()
              phipar[[j]][[2]][1] <- prior.trans.par[i,1]
              phipar[[j]][[2]][2] <- prior.trans.par[i,2]
              j <- j + 1
            } else {
              j <- j
            }
          }
        }
      } else if (all(nt > 1 & len.phipar == 0) == TRUE) {
        phipar[[1]][[2]] <- vector()
        phipar[[1]][[2]][1] <- 0
        phipar[[1]][[2]][2] <- 0
      }# end if - nt
    }# end if - flag.trans

# BETA:

    if (!is.null(object$contact)) {
       if (length(object$contact)/(n * n) == 1) {
         ni <- 1

         if (any(!is.null(beta.ini) | !is.null(prior.beta.dist) | !is.null(pro.beta.var) | !is.null(prior.beta.par)) == TRUE) {
          stop("As the model has only one contact network, The model does not have a network parameter beta, beta.ini, prior.dist.beta, prior.beta.var, and pro.beta.var must be assigned to their default values = NULL", .call = FALSE)
         }

         flag.beta <- 0

         betapar <- list(NULL)
         betapar[[1]] <- "uniform" # just to fill the list and it is not being used at all.
         betapar[[2]] <- vector()
         betapar[[2]][1] <- 0
         betapar[[2]][2] <- 0

       } else if (length(object$contact)/(n * n) > 1) {

        if (is.null(beta.ini)) {
          stop("epimcmc: The initial values beta.ini of the network parameters has to specified", call. = FALSE)
        }
         ni <- length(beta.ini)

         if (length(object$contact)/(n * n) != ni) {
             stop('epimcmc:  Dimension of beta  and the number of contact networks are not matching')
         }

         flag.beta <- 1

         if (any(is.null(prior.beta.dist) | !(all(prior.beta.dist %in% c("halfnormal", "gamma","uniform")))) == TRUE) {
           stop("epimcmc: Specify prior for beta as \"halfnormal\" ,  \"gamma\" or \"uniform\"  ", call. = FALSE)
         }


         if (length(prior.beta.dist) != ni) {
           stop('epimcmc: Specify prior distributions for each of the network parameters beta, prior.beta.dist')
         }

         # proposal variance for  beta
         prostdb <- vector(mode = "numeric", length = ni)
         if (length(pro.beta.var) != ni) {
           stop('epimcmc: Specify proposal variance for each beta parameter')
         }
         for (i in 1:ni) {
           prostdb[i] <- pro.beta.var[i]
         }


         betapar <- list(NULL)
         if (flag.beta == 0) {
           len.betapar <- 0
         } else {
           len.betapar <- ni - sum(prostdb == 0)
         }

         if (len.betapar != 0){
           j <- 1
           for (i in 1:ni) {
             if (prostdb[i] != 0){
               betapar[[j]] <- list(NULL)
               betapar[[j]][[1]] <- prior.beta.dist[i]
               j <- j + 1
             } else {
               j <- j
             }
           }
         } else {
           betapar[[1]] <- list(NULL)
           betapar[[1]][[1]] <- "uniform" # just to fill the list and it is not being used at all.
         }
         # prior selection for kernel parameters
         if (len.betapar != 0) {
           if (any(is.null(prior.beta.par) |  !is.matrix(prior.beta.par)) == TRUE) {
             stop('epimcmc: Specify prior.beta.par as a matrix with each row corresponds to each beta parameter')
           } else {
             j <- 1
             for (i in 1:ni) {
               if (prostdb[i] != 0){
                 betapar[[j]][[2]] <- vector()
                 betapar[[j]][[2]][1] <- prior.beta.par[i,1]
                 betapar[[j]][[2]][2] <- prior.beta.par[i,2]
                 j <- j + 1
               } else {
                 j <- j
               }
             }
           }
         } else if (len.betapar == 0) {
           betapar[[1]][[2]] <- vector()
           betapar[[1]][[2]][1] <- 0
           betapar[[1]][[2]][2] <- 0
         }# end if - ni
        }

    } else if (is.null(object$contact)) {

        if (is.null(object$XYcoordinates)) {
          stop("epimcmc: Specify contact network or x, y coordinates")
        }

        ni <- length(beta.ini)

        if (ni > 1) {
          stop("epimcmc: The input of beta.ini has more than one value while the considered distance-based ILM needs only one spatial parameter, beta.")
        }

        if (any(is.null(prior.beta.dist) | !(all(prior.beta.dist %in% c("halfnormal", "gamma","uniform")))) == TRUE) {
          stop("epimcmc: Specify prior for beta as \"halfnormal\" ,  \"gamma\" or \"uniform\"  ", call. = FALSE)
        }

        if (length(pro.beta.var) != ni) {
          stop('epimcmc: Specify only one proposal variance for the spatial parameter beta, pro.beta.var')
        }

        if (length(prior.beta.dist) != ni) {
          stop('epimcmc: Specify only one prior distribution for the spatial parameter beta, prior.beta.dist')
        }

        flag.beta <- 1

        prostdb <- pro.beta.var

        betapar <- list(NULL)
        len.betapar <- ni - sum(prostdb == 0)
        betapar[[1]] <- list(NULL)
        if (len.betapar != 0){
          betapar[[1]][[1]] <- prior.beta.dist
        } else {
          betapar[[1]][[1]] <- "uniform" # just to fill the list and it is not being used at all.
        }
        # prior selection for kernel parameters
        if (len.betapar != 0) {
          if (any((prior.beta.dist == "gamma") | (prior.beta.dist == "uniform")) == TRUE) {
            if (is.null(prior.beta.par)) {
              stop('epimcmc: Specify prior.beta.par as a vector or a 1 by 2 matrix.')
            } else {
              betapar[[1]][[2]] <- vector()
              betapar[[1]][[2]][1] <- prior.beta.par[1]
              betapar[[1]][[2]][2] <- prior.beta.par[2]
            }
          } else if (prior.beta.dist == "halfnormal") {
            if (is.null(prior.beta.par)) {
              stop('epimcmc: Specify the variance of the half normal prior distribution for prior.beta.par.')
            } else {
              betapar[[1]][[2]] <- vector()
              betapar[[1]][[2]][1] <- prior.beta.par[1]
              betapar[[1]][[2]][2] <- 0
            }
          }

        } else if (len.betapar == 0) {
          betapar[[1]][[2]] <- vector()
          betapar[[1]][[2]][1] <- 0
          betapar[[1]][[2]][2] <- 0
        }
    }

# SPARK parameter:

    if (!is.null(spark.ini)) {
      flag <- 1
      if (is.null(pro.spark.var)) {
        stop("epimcmc: Specify proposal variance for spark", call. = FALSE)
      }
      if (is.null(prior.spark.dist)) {
        stop("epimcmc22: Specify prior for spark", call. = FALSE)
      }
      if (any(is.null(prior.spark.dist) | !(prior.spark.dist %in% c("halfnormal", "gamma","uniform"))) == TRUE) {
        stop("epimcmc2: Specify prior for spark term as \"halfnormal\" ,  \"gamma\" or \"uniform\"  ", call. = FALSE)
      }
      sparkpar <- list(NULL)
      length(sparkpar) <- 2
      sparkpar[[1]] <- prior.spark.dist
      if (any((prior.spark.dist == "gamma") | (prior.spark.dist == "uniform")) == TRUE) {
        if (is.null(prior.spark.par)) {
          stop('epimcmc: Specify prior.spark.par as a vector or a 1 by 2 matrix.')
        } else {
          sparkpar[[1]][[2]] <- vector()
          sparkpar[[1]][[2]][1] <- prior.spark.par[1]
          sparkpar[[1]][[2]][2] <- prior.spark.par[2]
        }
      } else if (prior.spark.dist == "halfnormal") {
        if (is.null(prior.spark.par)) {
          stop('epimcmc: Specify the variance of the half normal prior distribution for prior.spark.par.')
        } else {
          sparkpar[[1]][[2]] <- vector()
          sparkpar[[1]][[2]][1] <- prior.spark.par[1]
          sparkpar[[1]][[2]][2] <- 0
        }
      }

    } else if (is.null(spark.ini)) {
      flag <-  0
      sparkpar <- list(NULL)
      sparkpar[[2]] <- vector()
      sparkpar[[2]][1] <- 0
      sparkpar[[2]][2] <- 0
    }

    # proposal variance for  spark

    if (flag == 1) {
      prostdsp <- pro.spark.var
    } else {
      prostdsp <- 0.0
    }



    if (all(flag == 0 & flag.trans == 0 & flag.beta == 0) == TRUE) {
      mcmcscale <- c(prostda)
    } else if (all(flag == 0 & flag.trans == 0 & flag.beta == 1) == TRUE) {
      mcmcscale <- c(prostda, prostdb)
    } else if (all(flag == 0 & flag.trans == 1 & flag.beta == 0) == TRUE) {
      mcmcscale <- c(prostda, prostdt)
    } else if (all(flag == 0 & flag.trans == 1 & flag.beta == 1) == TRUE) {
      mcmcscale <- c(prostda, prostdt, prostdb)
    } else if (all(flag == 1 & flag.trans == 0 & flag.beta == 0) == TRUE) {
      mcmcscale <- c(prostda, prostdsp)
    } else if (all(flag == 1 & flag.trans == 0 & flag.beta == 1) == TRUE) {
      mcmcscale <- c(prostda, prostdb, prostdsp)
    } else if (all(flag == 1 & flag.trans == 1 & flag.beta == 0) == TRUE) {
      mcmcscale <- c(prostda, prostdt, prostdsp)
    } else if (all(flag == 1 & flag.trans == 1 & flag.beta == 1) == TRUE) {
      mcmcscale <- c(prostda, prostdt, prostdb, prostdsp)
    }

    scalemc <- mcmcscale[mcmcscale !=0]


    if (all(flag == 0 & flag.trans == 0 & flag.beta == 0) == TRUE) {
      mcmcinit <- c(sus.par.ini)
    } else if (all(flag == 0 & flag.trans == 0 & flag.beta == 1) == TRUE) {
      mcmcinit <- c(sus.par.ini, beta.ini)
    } else if (all(flag == 0 & flag.trans == 1 & flag.beta == 0) == TRUE) {
      mcmcinit <- c(sus.par.ini, trans.par.ini)
    } else if (all(flag == 0 & flag.trans == 1 & flag.beta == 1) == TRUE) {
      mcmcinit <- c(sus.par.ini, trans.par.ini, beta.ini)
    } else if (all(flag == 1 & flag.trans == 0 & flag.beta == 0) == TRUE) {
      mcmcinit <- c(sus.par.ini, spark.ini)
    } else if (all(flag == 1 & flag.trans == 0 & flag.beta == 1) == TRUE) {
      mcmcinit <- c(sus.par.ini, beta.ini, spark.ini)
    } else if (all(flag == 1 & flag.trans == 1 & flag.beta == 0) == TRUE) {
      mcmcinit <- c(sus.par.ini, trans.par.ini, spark.ini)
    } else if (all(flag == 1 & flag.trans == 1 & flag.beta == 1) == TRUE) {
      mcmcinit <- c(sus.par.ini, trans.par.ini, beta.ini, spark.ini)
    }


    initmc <- mcmcinit[mcmcscale!=0]

    p <- function(xx) {

      ms <- which(mcmcscale[1: ns] == 0)
      nns <- ns - length(ms)
      if (length(ms) != 0) {
        xs <- vector("double", length = ns)
        xs[ms] <- mcmcinit[ms]
        xs[-ms] <- xx[1: nns]
      } else {
        xs <- xx[1: nns]
      }

      if (flag.trans == 1) {
        mt <- which(mcmcscale[(1+ns): (ns+nt)] == 0)
        nnt <- nt - length(mt)
        if (length(mt) != 0) {
          xt <- vector("double", length = nt)
          xt[mt] <- mcmcinit[mt+ns]
          xt[-mt] <- xx[(1+nns): (nns+nnt)]
        } else {
          xt <- xx[(1+nns): (nns+nnt)]
        }
      } else {
        nnt <- 0
        xt <- NULL
      }

      if (flag.beta == 1) {
        if (flag.trans == 1) {
          mi <- which(mcmcscale[(1+ns+nt): (ns+nt+ni)] == 0)
        } else {
          mi <- which(mcmcscale[(1+ns): (ns+ni)] == 0)
        }
        nni <- ni - length(mi)

        if (flag.trans == 1) {
          if (length(mi) != 0) {
            xi <- vector("double", length = ni)
            xi[mi] <- mcmcinit[mi+ns+nt]
            xi[-mi] <- xx[(1+nns+nnt): (nns+nnt+nni)]
          } else {
            xi <- xx[(1+nns+nnt): (nns+nnt+nni)]
          }
        } else {
          if (length(mi) != 0) {
            xi <- vector("double", length = ni)
            xi[mi] <- mcmcinit[mi+ns+nnt]
            xi[-mi] <- xx[(1+nns+nnt): (nns+nnt+nni)]
          } else {
            xi <- xx[(1+nns+nnt): (nns+nnt+nni)]
          }
        }
      } else {
        nni <- 0
        xi <- NULL
      }

      if (flag == 1) {
        if(length(xx) < nns+nnt+nni) {
          xsp <- xx[nns+nnt+nni+1]
        } else if (all(length(xx) == nns+nnt+nni & mcmcscale[ns+nt+ni+1] == 0) == TRUE) {
          xsp <- mcmcinit[ns+nt+ni+1]
        }
      } else {
        xsp <- NULL
      }


      loglikeILM <- epilike(object = object, tmin = tmin, tmax = tmax, sus.par = xs,
                            trans.par = xt, beta = xi, spark = xsp,
                            Sformula = Sformula, Tformula = Tformula)

      if (nns != 0) {
        palpha = 0
        for (i in 1:nns) {
          if (alphapar[[i]][[1]] == "uniform") {
            palpha <- palpha + dunif(xx[i], min = alphapar[[i]][[2]][1],
                                     max = alphapar[[i]][[2]][2], log = TRUE)
          } else if (alphapar[[i]][[1]] == "gamma") {
            palpha <- palpha + dgamma(xx[i], shape = alphapar[[i]][[2]][1],
                                      rate = alphapar[[i]][[2]][2], log = TRUE)
          } else if (alphapar[[i]][[1]] == "halfnormal") {
            palpha <- palpha + dhalfnorm(xx[i], scale = alphapar[[i]][[2]][1], log = TRUE)
          }
        }
      } else {
        palpha = 0
      }

      if (all(flag.trans == 1 & nnt != 0) == TRUE) {
        pphi = 0
        for (i in 1:nnt) {
          if (phipar[[i]][[1]] == "uniform") {
            pphi <- pphi + dunif(xx[i+nns], min = phipar[[i]][[2]][1],
                                 max = phipar[[i]][[2]][2], log = TRUE)
          } else if (phipar[[i]][[1]] == "gamma") {
            pphi <- pphi + dgamma(xx[i+nns], shape = phipar[[i]][[2]][1],
                                  rate = phipar[[i]][[2]][2], log = TRUE)
          } else if (phipar[[i]][[1]] == "halfnormal") {
            pphi <- pphi + dhalfnorm(xx[i+nns], scale = phipar[[i]][[2]][1], log = TRUE)
          }
        }
      } else {
        pphi = 0
      }

      if (nni != 0) {
        pbeta = 0
        for (i in 1:nni) {
          if (betapar[[i]][[1]] == "uniform") {
            pbeta <- pbeta + dunif(xx[i+nns+nnt], min = betapar[[i]][[2]][1],
                                   max = betapar[[i]][[2]][2], log = TRUE)
          } else if (betapar[[i]][[1]] == "gamma") {
            pbeta <- pbeta + dgamma(xx[i+nns+nnt], shape = betapar[[i]][[2]][1],
                                    rate = betapar[[i]][[2]][2], log = TRUE)
          } else if (betapar[[i]][[1]] == "halfnormal") {
            pbeta <- pbeta + dhalfnorm(xx[i+nns+nnt], scale = betapar[[i]][[2]][1], log = TRUE)
          }
        }
      } else {
        pbeta = 0
      }


      if (flag == 1) {
        if (sparkpar[[1]] == "uniform") {
          pspark <- dunif(xx[nns+nnt+nni+1], min = sparkpar[[2]][1],
                          max = sparkpar[[2]][2], log = TRUE)
        } else if (sparkpar[[1]] == "gamma") {
          pspark <- dgamma(xx[nns+nnt+nni+1], shape = sparkpar[[2]][1],
                           rate = sparkpar[[2]][2], log = TRUE)
        } else if (sparkpar[[1]] == "halfnormal") {
          pspark <- pspark + dhalfnorm(xx[nns+nnt+nni+1], scale = sparkpar[[2]][1], log = TRUE)
        }
      } else {
        pspark = 0.0
      }

      return(loglikeILM + palpha + pphi + pbeta + pspark)
    }

    if (is.null(object$contact)) {
        kernel.type <- "distance-based discrete-time ILM"
        result <- MCMC(p = p, n = niter, init = initmc, scale = scalemc,
        adapt = adapt, acc.rate = acc.rate)
    } else if (!is.null(object$contact)) {
        kernel.type <- "network-based discrete-time ILM"
        result <- MCMC(p = p, n = niter, init = initmc, scale = scalemc,
        adapt = adapt, acc.rate = acc.rate)
    }

    alphanames <- c()
    for (i in 1:ns) {
        alphanames <- c(alphanames, paste("alpha.", i, sep = ""))
    }

    if (flag.trans == 1) {
        phinames <- c()
        for (i in 1:nt) {
            phinames <- c(phinames, paste("phi.", i, sep = ""))
        }
    }

    if (flag.beta == 1) {
      betanames <- c()
      for (i in 1:ni) {
          betanames <- c(betanames, paste("beta.", i, sep = ""))
      }
    }

    if (all(flag == 0 & flag.trans == 0 & flag.beta == 0) == TRUE) {
        nnames <- c(alphanames)
    } else if (all(flag == 0 & flag.trans == 0 & flag.beta == 1) == TRUE) {
        nnames <- c(alphanames, betanames)
    } else if (all(flag == 1 & flag.trans == 0 & flag.beta == 0) == TRUE) {
        nnames <- c(alphanames, "spark")
    } else if (all(flag == 1 & flag.trans == 0 & flag.beta == 1) == TRUE) {
        nnames <- c(alphanames, betanames, "spark")
    } else if (all(flag == 0 & flag.trans == 1 & flag.beta == 0) == TRUE) {
        nnames <- c(alphanames, phinames)
    } else if (all(flag == 0 & flag.trans == 1 & flag.beta == 1) == TRUE) {
        nnames <- c(alphanames, phinames, betanames)
    } else if (all(flag == 1 & flag.trans == 1 & flag.beta == 0) == TRUE) {
        nnames <- c(alphanames, phinames, "spark")
    } else if (all(flag == 1 & flag.trans == 1 & flag.beta == 1) == TRUE) {
        nnames <- c(alphanames, phinames, betanames, "spark")
    }

    finalnames <- nnames[mcmcscale!=0]

  	Estimates <- data.frame(result$samples)
  	colnames(Estimates) <- finalnames

  	Loglikelihood <- c(result$log.p)
#  	colnames(Loglikelihood) <- "Loglikelihood"

    ki <- match(nnames,finalnames)

    nki <- which(nnames %in% finalnames)

    fullsamples <- matrix(0, ncol = length(nnames), nrow = niter)
    fullsamples[,nki] <- result$samples
    nnki <- which(is.na(ki))

    if (length(nnki) == 1) {
      fullsamples[,nnki] <- rep(mcmcinit[nnki], niter)
    } else if (length(nnki) > 1){
      for (i in 1:length(nnki)) {
        fullsamples[,nnki[i]] <- rep(mcmcinit[nnki[i]], niter)
      }
    }
  }

  # mcmc result as a coda object

  if (flag.trans == 1) {
    n.trans.par <- nt
  } else {
    n.trans.par <- 0
  }
  if (flag.beta == 1) {
    n.ker.par <- ni
  } else {
    n.ker.par <- 0
  }

  output <- list(type = object$type, kernel.type = kernel.type, Estimates = Estimates,
                 Loglikelihood = Loglikelihood, Fullsample = fullsamples,
                 n.sus.par = ns, n.trans.par = n.trans.par, n.ker.par = n.ker.par)


  class(output) <- "epimcmc"
  output

} # End of function

summary.epimcmc <- function(object, ...) {
    if (!is(object, "epimcmc")) {
        stop("The object has to be of \"epimcmc\" class", call. = FALSE)
    }
    if (object$type == "SI") {
        cat(paste("Model:", object$type, object$kernel.type,"\n"))
        cat(paste("Method: Markov chain Monte Carlo (MCMC) \n"))
        summary(window(mcmc(object$Estimates), ...))
    } else if (object$type == "SIR") {
        cat(paste("Model:", object$type, object$kernel.type,"\n"))
        cat(paste("Method: Markov chain Monte Carlo (MCMC) \n"))
        summary(window(mcmc(object$Estimates), ...))
    }
}


plot.epimcmc <- function(x, partype, start = 1, end = NULL, thin = 1, ...) {
    if (is.null(end)) {
        end = nrow(x$Estimates)
    }
    if (!is(x, "epimcmc")) {
        stop("The object x has to be of \"epimcmc\" class", call. = FALSE)
    }
    if (partype == "parameter") {
        plot(window(as.mcmc(x$Estimates), start = start, end = end, thin = thin), ...)
    } else if (partype == "loglik") {
        plot(window(as.mcmc(x$Loglikelihood), start = start, end = end, thin = thin), ...)
    }
}
