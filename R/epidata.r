################################################################################
# Part of the R/EpiILM package
#
# AUTHORS:
#         Vineetha Warriyar. K. V. <vineethawarriyar.kod@ucalgary.ca>,
#         Waleed Almutiry <wkmtierie@qu.edu.sa>, and
#         Rob Deardon <robert.deardon@ucalgary.ca>
#
# Algorithm based on:
#         Deardon R, Brooks, S. P., Grenfell, B. T., Keeling, M. J., Tildesley,
#         M. J., Savill, N. J., Shaw, D. J.,  Woolhouse, M. E. (2010).
#         Inference for individual level models of infectious diseases in large
#         populations. Statistica Sinica, 20, 239-261.
#
# Free software under the terms of the GNU General Public License, version 2,
# a copy of which is available at http://www.r-project.org/Licenses/.
################################################################################

epidata <- function(type, n, tmin = NULL, tmax, sus.par, trans.par = NULL, beta = NULL, spark = NULL,
                    Sformula = NULL, Tformula = NULL, x = NULL, y = NULL,
                    inftime = NULL, infperiod = NULL, contact = NULL) {

  # Error checks for input arguments
  if (any(is.null(type) | !(type %in% c("SI", "SIR"))) == TRUE) {
       stop("Specify type as \"SI\" or \"SIR\" ", call. = FALSE)
  }

  ns <- length(sus.par)
  if  (!is.null(trans.par)) {

    nt <- length(trans.par)
    flag.trans <- 1
    phi <- trans.par

  } else {

    nt <- 1
    flag.trans <- 0
    phi <- 1

  }

  if (all(is.null(contact) &  (is.null(x) | is.null(y))) == TRUE) {
      stop('epidata: Specify contact network or x, y coordinates')
  }

  if (is.null(contact)) {
    ni <- length(beta)
    if (ni != 1) {
      stop("epidata: The input of beta has more than one value while the considered distance-based ILM needs only one spatial parameter, beta.", call. = FALSE)
    }
    if (!is.null(x)) {
      if (any((length(y) != n) | (length(x) != n)) == TRUE) {
        stop('epidata: Length of x or y is not compatible ')
      }
    }
  } else {
    if (length(contact)/(n * n) == 1) {
      ni <- 1
      if (is.null(beta)) {
        beta <- 1
      } else {
        stop("epidata: As the model has only one contact network, The model does not need a network parameter beta, beta must be assigned to its default values = NULL", call. = FALSE)
      }
    } else if (length(contact)/(n * n) > 1) {
      ni <- length(beta)
      if (length(contact)/(n * n) != ni) {
        stop('epidata: Dimension of beta and the number of contact networks are not matching')
      }
    }
    network <- array(contact, c(n, n, ni))
  }

  if (is.null(tmin)){
    tmin <- 1
  }

  if (is.null(spark)) {
    spark <- 0
  }

  if (!is.null(inftime)) {
      if ((length(inftime) != n)) {
        stop('epidata: Length of inftime is not compatible ')
      }
  } else {
        inftime <- rep(0, n)
  }

  if (all(is.null(infperiod) & type == "SIR") == TRUE) {
    stop(' epidata: Specify removal, infperiod')
  }

  if (!is.null(infperiod)) {
      if (length(infperiod) != n) {
        stop('epidata: Length of infperiod is not compatible')
      }
      if (type == "SI") {
        stop('epidata: Type must be "SIR"')
      }
       remt <- rep(0, n)
  }

  # formula for susceptibility (covariate) function
   if (!is.null(Sformula)) {
     covmat.sus <- model.matrix(Sformula)

     if (all((ncol(covmat.sus) == length(all.vars(Sformula))) & (ns != length(all.vars(Sformula)))) == TRUE) {
       stop('epidata: Check Sformula (no intercept term) and the dimension of sus.par')
     }

     if (all((ncol(covmat.sus) > length(all.vars(Sformula))) & (ns != ncol(covmat.sus))) == TRUE) {
       stop('epidata: Check Sformula (intercept term) and the dimension of sus.par')
     }
   } else {
       if (ns == 1) {
         covmat.sus <- matrix(1.0, nrow = n, ncol = ns)
       }
       if (ns > 1) {
         stop('epidata: Please specify the susceptibility covariate')
       }
   }

   # formula for transmissibility (covariate) function
   if (flag.trans == 1) {
     if (is.null(Tformula)) {
       stop("epidata: Tformula is missing. It has to be specified with no intercept term and number of columns equal to the length of trans.par", call. = FALSE)
     } else if (!is.null(Tformula)) {
   		covmat.trans <- model.matrix(Tformula)

   		if (all((ncol(covmat.trans) == length(all.vars(Tformula))) & (nt != length(all.vars(Tformula)))) == TRUE) {
   			stop("epidata: Check Tformula. It has to be with no intercept term and number of columns equal to the length of trans.par", call. = FALSE)
   		}
    }
   } else if (flag.trans == 0) {
   	covmat.trans <- matrix(1.0, nrow = n, ncol = nt)
   }

 # Calling fortran subroutine - Purely Spatial: Susceptible-Infectious(SI)
  if (all((type == "SI") & is.null(contact)) == TRUE) {
    tmp <- .Fortran("dataxy",
                    x = as.vector(x, mode = "double"),
                    y = as.vector(y, mode = "double"),
                     n = as.integer(n),
                     tmin = as.integer(tmin),
                     tmax = as.integer(tmax),
                     ns = as.integer(ns),
                     nt = as.integer(nt),
                     ni = as.integer(ni),
                     alpha = as.vector(sus.par, mode = "double"),
                     phi = as.vector(phi, mode = "double"),
                     beta =as.vector(beta, mode = "double"),
                     spark = as.numeric(spark),
                     covmatsus = matrix(as.double(covmat.sus), ncol = ncol(covmat.sus), nrow = n),
                     covmattrans = matrix(as.double(covmat.trans), ncol = ncol(covmat.trans), nrow = n),
                     tau = as.vector(inftime, mode = "integer")
                     )

    result1 <- list(type = type, XYcoordinates = cbind(x, y), contact = NULL, inftime = tmp$tau)
  }

 # Calling fortran subroutine - Purely Spatial: Susceptible-Infectious-Removed (SIR)
  if (all((type == "SIR") & is.null(contact)) == TRUE) {
    tmp <- .Fortran("dataxysir",
                     n = as.integer(n),
                     tmin = as.integer(tmin),
                     tmax = as.integer(tmax),
                     ns = as.integer(ns),
                     nt = as.integer(nt),
                     ni = as.integer(ni),
                     alpha = as.vector(sus.par, mode = "double"),
                     phi = as.vector(phi, mode = "double"),
                     beta =as.vector(beta, mode = "double"),
                     spark = as.numeric(spark),
                     covmatsus = matrix(as.double(covmat.sus), ncol = ncol(covmat.sus), nrow = n),
                     covmattrans = matrix(as.double(covmat.trans), ncol = ncol(covmat.trans), nrow = n),
                     lambda = as.vector(infperiod, mode = "integer"),
                     x = as.vector(x, mode = "double"),
                     y = as.vector(y, mode = "double"),
                     tau = as.vector(inftime, mode = "integer"),
                     remt = as.vector(remt, mode = "integer")
                     )

    result1 <- list(type = type, XYcoordinates = cbind(x, y), contact = NULL, inftime = tmp$tau, remtime = tmp$remt)
  }

  # Calling fortran subroutine - Contact networks: Susceptible-Infectious (SI)
  if (all((type == "SI") & !is.null(contact)) == TRUE) {
    tmp <- .Fortran("datacon",
                     n = as.integer(n),
                     tmin = as.integer(tmin),
                     tmax = as.integer(tmax),
                     ns = as.integer(ns),
                     nt = as.integer(nt),
                     ni = as.integer(ni),
                     alpha = as.numeric(sus.par),
                     phi = as.numeric(phi),
                     beta = as.numeric(beta),
                     spark = as.numeric(spark),
                     covmatsus = matrix(as.double(covmat.sus), ncol = ncol(covmat.sus), nrow = n),
                     covmattrans = matrix(as.double(covmat.trans), ncol = ncol(covmat.trans), nrow = n),
                     network = as.vector(network),
                     tau = as.integer(inftime)
                     )

    result1 <- list(type = type, XYcoordinates = cbind(x, y), contact = contact, inftime = tmp$tau)
  } else if (all((type == "SIR") & !is.null(contact)) == TRUE) {
 # Calling fortran subroutine - Contact networks: Susceptible-Infectious-Removed (SIR)
    tmp <- .Fortran("dataconsir",
                     n = as.integer(n),
                     tmin = as.integer(tmin),
                     tmax = as.integer(tmax),
                     ns = as.integer(ns),
                     nt = as.integer(nt),
                     ni = as.integer(ni),
                     lambda = as.integer(infperiod),
                     alpha = as.numeric(sus.par),
                     phi = as.numeric(phi),
                     beta = as.numeric(beta),
                     spark = as.numeric(spark),
                     covmatsus = matrix(as.double(covmat.sus), ncol = ncol(covmat.sus), nrow = n),
                     covmattrans = matrix(as.double(covmat.trans), ncol = ncol(covmat.trans), nrow = n),
                     network = as.vector(network),
                     tau = as.integer(inftime),
                     remt = as.integer(remt)
                     )

    result1 <- list(type = type, XYcoordinates = cbind(x, y), contact = contact, inftime = tmp$tau, remtime = tmp$remt)
  }
  class(result1) <- "epidata"

  result1
  # End of function
}
