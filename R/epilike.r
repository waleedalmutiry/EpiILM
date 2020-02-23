################################################################################
# Part of the R/EpiILM package
#
# AUTHORS:
#         Vineetha Warriyar. K. V. <vineethawarriyar.kod@ucalgary.ca>,
#         Waleed Almutiry <wkmtierie@qu.edu.sa>, and
#         Rob Deardon <robert.deardon@ucalgary.ca>
#
#
# Algorithm based on:
#                Deardon R, Brooks, S. P., Grenfell, B. T., Keeling, M. J., Tildesley,
#                M. J., Savill, N. J., Shaw, D. J.,  Woolhouse, M. E. (2010).
#                Inference for individual level models of infectious diseases in large
#                populations. Statistica Sinica, 20, 239-261.
#
#
# Free software under the terms of the GNU General Public License, version 2,
# a copy of which is available at http://www.r-project.org/Licenses/.
################################################################################

epilike <- function(object, tmin = NULL, tmax, sus.par, trans.par = NULL,
                    beta = NULL, spark = NULL, Sformula = NULL, Tformula = NULL) {

if (class(object) != "epidata") {
   stop("The object must be in a class of \"epidata\"", call. = FALSE)
} else {

  # Error checks for input arguments
  if (is.null(object$type) | !(object$type %in% c("SI", "SIR"))) {
    stop("epilike: Specify type as \"SI\" or \"SIR\" ", call. = FALSE)
  }

  ns <- length(sus.par)
  if (!is.null(trans.par)){
    nt <- length(trans.par)
    flag.trans <- 1
    phi <- trans.par
  } else if (is.null(trans.par)) {
    nt <- 1
    flag.trans <- 0
    phi <- 1
    trans.par <- 1
  }

  if (!is.vector(object$inftime)) {
    stop('epilike: inftime is not a vector')
  }
  n <- length(object$inftime)

  if (is.null(object$contact) &  is.null(object$XYcoordinates)) {
      stop('epilike: Specify contact network or x, y coordinates')
  }

  if (is.null(object$remtime) & object$type == "SIR") {
    stop(' epilike: Specify removal times')
  }

  if(object$type == "SIR") {
      infperiod <- object$remtime - object$inftime
  }

  if (is.null(tmin)) {
    tmin <- 1
  }
  if (is.null(spark)) {
    spark <- 0
  }

  if (is.null(object$contact)) {
    ni <- length(beta)
    if (is.null(beta)) {
      stop("epilike: A scalar value of the spatial parameter beta must be specified", call. = FALSE)
    }
    if (ni != 1) {
      stop("epilike: The input of beta has more than one value while the considered distance-based ILM needs only one spatial parameter, beta.", call. = FALSE)
    }
  } else {
    if (length(object$contact)/(n * n) == 1) {
      ni <- 1
      if (is.null(beta)) {
        beta <- 1
      } else {
        stop("epilike: As the model has only one contact network, The model does not need a network parameter beta, beta must be assigned to its default values = NULL", call. = FALSE)
      }
    } else if (length(object$contact)/(n * n) > 1) {
      ni <- length(beta)
      if (is.null(beta)) {
        stop("epilike: A vector values of the network parameters beta must be specified", call. = FALSE)
      }
      if (length(object$contact)/(n * n) != ni) {
        stop('epilike: Dimension of beta and the number of contact networks are not matching')
      }
    }
    network <- array(object$contact, c(n, n, ni))
  }

  # formula for susceptibility function
  	if (!is.null(Sformula)) {
  		covmat.sus <- model.matrix(Sformula)

  		if ((ncol(covmat.sus) == length(all.vars(Sformula))) & (ns != length(all.vars(Sformula)))) {
  			stop("epilike: Check Sformula (no intercept term) and the dimension of susceptibility parameters", call. = FALSE)
  		} else if ((ncol(covmat.sus) > length(all.vars(Sformula))) & (ns != ncol(covmat.sus))) {
  			stop("epilike: Check Sformula (intercept term) and the dimension of susceptibility parameters", call. = FALSE)
  		}
  	} else {
  		if (ns == 1) {
  			covmat.sus <- matrix(1.0, nrow = n, ncol = ns)
  		}
  	}

    # formula for transmissibility function
    if (flag.trans == 1) {
      if (is.null(Tformula)) {
        stop("epilike: Tformula is missing. It has to be specified with no intercept term and number of columns equal to the length of trans.par", call. = FALSE)
      } else if (!is.null(Tformula)) {
    		covmat.trans <- model.matrix(Tformula)

    		if ((ncol(covmat.trans) == length(all.vars(Tformula))) & (nt != length(all.vars(Tformula)))) {
    			stop("epilike: Check Tformula. It has to be with no intercept term and number of columns equal to the length of trans.par", call. = FALSE)
    		}
      }
    } else if (flag.trans == 0) {
    	covmat.trans <- matrix(1.0, nrow = n, ncol = nt)
    }

# Calling fortran subroutines for Purely Spatial models : SI and SIR

  if ((object$type == "SI") & is.null(object$contact)) {

    tmp1 <- .Fortran("like",
                   x = as.vector(object$XYcoordinates[,1], mode = "double"),
                   y = as.vector(object$XYcoordinates[,2], mode = "double"),
                   tau = as.vector(object$inftime, mode = "integer"),
                   n = as.integer(n),
                   tmin = as.integer(tmin),
                   tmax = as.integer(tmax),
                   ns = as.integer(ns),
                   nt = as.integer(nt),
                   ni = as.integer(ni),
                   alpha = as.vector(sus.par, mode = "double"),
                   phi = as.vector(trans.par, mode = "double"),
                   beta = as.vector(beta, mode = "double"),
                   spark = as.double(spark),
                   covmatsus = matrix(as.double(covmat.sus), ncol = ncol(covmat.sus), nrow = n),
                   covmattrans = matrix(as.double(covmat.trans), ncol = ncol(covmat.trans), nrow = n),
                   val = as.double(0))

  } else if ((object$type == "SIR") & is.null(object$contact)) {

    tmp1 <- .Fortran("likesir",
                    x = as.vector(object$XYcoordinates[,1], mode = "double"),
                    y = as.vector(object$XYcoordinates[,2], mode = "double"),
                    tau = as.vector(object$inftime, mode = "integer"),
                    lambda = as.vector(infperiod, mode = "integer"),
                    n = as.integer(n),
                    tmin = as.integer(tmin),
                    tmax = as.integer(tmax),
                    ns = as.integer(ns),
                    nt = as.integer(nt),
                    ni = as.integer(ni),
                    alpha = as.vector(sus.par, mode = "double"),
                    phi = as.vector(trans.par, mode = "double"),
                    beta = as.vector(beta, mode = "double"),
                    spark = as.double(spark),
                    covmatsus = matrix(as.double(covmat.sus), ncol = ncol(covmat.sus), nrow = n),
                    covmattrans = matrix(as.double(covmat.trans), ncol = ncol(covmat.trans), nrow = n),
                    val = as.double(0))

  } else if ((object$type == "SI") & !is.null(object$contact)) {
    tmp1 <- .Fortran("likecon",
                    tau = as.vector(object$inftime, mode = "integer"),
                    n = as.integer(n),
                    ns = as.integer(ns),
                    nt = as.integer(nt),
                    ni = as.integer(ni),
                    tmin = as.integer(tmin),
                    tmax = as.integer(tmax),
                    alpha = as.vector(sus.par, mode = "double"),
                    phi = as.vector(trans.par, mode = "double"),
                    beta = as.vector(beta, mode = "double"),
                    spark = as.double(spark),
                    covmatsus = matrix(as.double(covmat.sus), ncol = ncol(covmat.sus), nrow = n),
                    covmattrans = matrix(as.double(covmat.trans), ncol = ncol(covmat.trans), nrow = n),
                    network = as.vector(network),
                    val = as.double(0))

  } else if ((object$type == "SIR") & !is.null(object$contact)) {

    # Calling fortran subroutines for Contact network models: SI and SIR

    if (length(object$contact)/(n * n) != ni) {
      stop('epilike:  Dimension of beta  and the number of contact networks are not matching')
    }
    network <- array(object$contact, c(n, n, ni))

    tmp1 <- .Fortran("likeconsir",
                      tau = as.vector(object$inftime, mode = "integer"),
                      lambda = as.vector(infperiod, mode = "integer"),
                      n = as.integer(n),
                      ns = as.integer(ns),
                      nt = as.integer(nt),
                      ni = as.integer(ni),
                      tmin = as.integer(tmin),
                      tmax = as.integer(tmax),
                      alpha = as.vector(sus.par, mode = "double"),
                      phi = as.vector(trans.par, mode = "double"),
                      beta = as.vector(beta, mode = "double"),
                      spark = as.double(spark),
                      covmatsus = matrix(as.double(covmat.sus), ncol = ncol(covmat.sus), nrow = n),
                      covmattrans = matrix(as.double(covmat.trans), ncol = ncol(covmat.trans), nrow = n),
                      network = as.vector(network),
                      val = as.double(0))
  }
  return(tmp1$val)
}# End of function
}
