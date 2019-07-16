################################################################################
# Part of the R/EpiILM package
#
# AUTHORS:
#         Vineetha Warriyar. K. V. <vineethawarriyar.kod@ucalgary.ca> and
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

epilike <- function(type, x = NULL, y = NULL, inftime, infperiod = NULL, tmin = NULL, tmax, alpha,
                   beta, spark = NULL, Sformula = NULL, contact = NULL) {


  # Error checks for input arguments
  if (is.null(type) || !(type %in% c("SI", "SIR"))) {
    stop("elike: Specify type as \"SI\" or \"SIR\" ", call. = FALSE)
  }
  
  ns <- length(alpha)
  ni <- length(beta)
  
  if (!is.vector(inftime)) {
    stop('epilike: inftime is not a vector')
  }
  n <- length(inftime)
  
  if (is.null(contact) &  (is.null(x) || is.null(y))) {
      stop('epiBR0: Specify contact network or x, y coordinates')
  }
  if (!is.null(x)) {
      if ((length(y) != n) || (length(x) != n)) {
        stop('epilike: Length of x/y is not compatible')
      }
  }
  if (is.null(infperiod) && type == "SIR") {
    stop(' epilike: Specify removal distance,infperiod ')
  }
  if (!is.null(infperiod)) {
      if (length(infperiod) != n) {
        stop('epilike: Length of infperiod is not compatible')
      }
      if(type == "SI") {
        stop('epilike: Type must be "SIR"')
      }
  }
  if (is.null(tmin)) {
    tmin <- 1
  }
  if (is.null(spark)) {
    spark <- 0
  }

# formula for susceptibility function
  if(!is.null(Sformula)) {
    covmat <- model.matrix(Sformula)
    
    if ((ncol(covmat) == length(all.vars(Sformula))) & (ns != length(all.vars(Sformula)))) {
      stop('epilike: Check Sformula (no intercept term) and the dimension of alpha')
    }
    if ((ncol(covmat) > length(all.vars(Sformula))) & (ns != ncol(covmat))) {
      stop('epilike: Check Sformula (intercept term) and the dimension of alpha')
    }
  } else {
      if (ns == 1) {
        covmat <- matrix(1.0, nrow = n, ncol = ns)
      }
      if (ns > 1) {
        stop('epilike: Please specify covariate')
      }
  }
 
 val <- 0.0
# Calling fortran subroutines for Purely Spatial models : SI and SIR
  if ((type == "SI") && is.null(contact)) {
    tmp1 <- .Fortran("like",
                     x = as.numeric(x),
                     y = as.numeric(y),
                     tau = as.integer(inftime),
                     n = as.integer(n),
                     tmin = as.integer(tmin),
                     tmax = as.integer(tmax),
                     ns = as.integer(ns),
                     ni = as.integer(ni),
                     alpha = as.double(alpha),
                     beta = as.double(beta),
                     spark = as.double(spark),
                     covmat = as.vector(covmat),
                     val = as.double(val))
  }
  if ((type == "SIR") && is.null(contact)) {
    tmp1 <- .Fortran("likesir",
                      x = as.numeric(x),
                      y = as.numeric(y),
                      tau = as.integer(inftime),
                      lambda = as.integer(infperiod),
                      n = as.integer(n),
                      tmin = as.integer(tmin),
                      tmax = as.integer(tmax),
                      ns = as.integer(ns),
                      ni = as.integer(ni),
                      alpha = as.double(alpha),
                      beta = as.double(beta),
                      spark = as.double(spark),
                      covmat = as.vector(covmat),
                      val = as.double(val))
  }

# Calling fortran subroutines for Contact network models: SI and SIR
  if (!is.null(contact)) {
    if (length(contact)/(n * n) != ni) {
      stop('epilike:  Dimension of beta  and the number of contact networks are not matching')
    }
    network <- array(contact, c(n, n, ni))
  }
  if ((type == "SI") && !is.null(contact)) {
    tmp1 <- .Fortran("likecon",
                     tau = as.integer(inftime),
                     n = as.integer(n),
                     ns = as.integer(ns),
                     ni = as.integer(ni),
                     tmin = as.integer(tmin),
                     tmax = as.integer(tmax),
                     alpha = as.numeric(alpha),
                     beta = as.numeric(beta),
                     spark = as.double(spark),
                     covmat = as.vector(covmat),
                     network = as.vector(network),
                     val = as.double(val))
  }
  if ((type == "SIR") && !is.null(contact)) {
    tmp1 <- .Fortran("likeconsir",
                      tau = as.integer(inftime),
                      lambda = as.integer(infperiod),
                      n = as.integer(n),
                      ns = as.integer(ns),
                      ni = as.integer(ni),
                      tmin = as.integer(tmin),
                      tmax = as.integer(tmax),
                      alpha = as.numeric(alpha),
                      beta = as.numeric(beta),
                      spark = as.double(spark),
                      covmat = as.vector(covmat),
                      network = as.vector(network),
                      val = as.double(val))
  }
  return(tmp1$val)
# End of function
}

