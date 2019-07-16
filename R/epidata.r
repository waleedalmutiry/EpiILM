################################################################################
# Part of the R/EpiILM package
#
# AUTHORS:
#         Vineetha Warriyar. K. V. <vineethawarriyar.kod@ucalgary.ca> and
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

epidata <- function(type, n, tmin = NULL, tmax, alpha, beta, spark = NULL, Sformula = NULL,
                    x = NULL, y = NULL, inftime = NULL, infperiod = NULL, contact = NULL,
                    tempseed = NULL) {

  # Error checks for input arguments
  if (is.null(type) || !(type %in% c("SI", "SIR"))) {
       stop("Specify type as \"SI\" or \"SIR\" ", call. = FALSE)
  }
  
  if (is.null(tempseed)) {
      tempseed <- 0
  }
  
  ns <- length(alpha)
  ni <- length(beta)
  
  if (is.null(tmin)){
    tmin <- 1
  }
  
  if (is.null(spark)) {
    spark <- 0
  }
  
  if (is.null(contact) &  (is.null(x) || is.null(y))) {
      stop('epidata: Specify contact network or x, y coordinates')
  }
  
  if (!is.null(x)) {
    if ((length(y) != n) || (length(x) != n)) {
      stop('epidata: Length of x or y is not compatible ')
    }
  }
  
  if (!is.null(inftime)) {
      if ((length(inftime) != n)) {
        stop('epidata: Length of inftime is not compatible ')
      }
  } else {
        inftime <- rep(0, n)
  }
  
  if (is.null(infperiod) && type == "SIR") {
    stop(' epidata: Specify removal distance, infperiod ')
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
    covmat <- model.matrix(Sformula)
    
    if ((ncol(covmat) == length(all.vars(Sformula))) & (ns != length(all.vars(Sformula)))) {
      stop('epidata: Check Sformula (no intercept term) and the dimension of alpha')
    }
    
    if ((ncol(covmat) > length(all.vars(Sformula))) & (ns != ncol(covmat))) {
      stop('epidata: Check Sformula (intercept term) and the dimension of alpha')
    }
  } else {
      if (ns == 1) {
        covmat <- matrix(1.0, nrow = n, ncol = ns)
      }
      if (ns > 1) {
        stop('epidata: Please specify covariate')
      }
  }
 
 # Calling fortran subroutine - Purely Spatial: Susceptible-Infectious(SI)
  if ((type == "SI") && is.null(contact)) {
    tmp <- .Fortran("dataxy",
                     x = as.double(x),
                     y = as.double(y),
                     n = as.integer(n),
                     tmin = as.integer(tmin),
                     tmax = as.integer(tmax),
                     ns = as.integer(ns),
                     ni = as.integer(ni),
                     alpha = as.numeric(alpha),
                     beta = as.numeric(beta),
                     spark = as.numeric(spark),
                     covmat = as.vector(covmat),
                     tau = as.integer(inftime),
                     tempseed = as.integer(tempseed)
                     )
                     
    result1 <- list(inftime = tmp$tau)
  }
 
 # Calling fortran subroutine - Purely Spatial: Susceptible-Infectious-Removed (SIR)
  if ((type == "SIR") && is.null(contact)) {
    tmp <- .Fortran("dataxysir",
                     n = as.integer(n),
                     tmin = as.integer(tmin),
                     tmax = as.integer(tmax),
                     ns = as.integer(ns),
                     ni = as.integer(ni),
                     alpha = as.numeric(alpha),
                     beta = as.numeric(beta),
                     spark = as.numeric(spark),
                     covmat = as.vector(covmat),
                     lambda = as.integer(infperiod),
                     x = as.double(x),
                     y = as.double(y),
                     tau = as.integer(inftime),
                     remt = as.integer(remt),
                     tempseed = as.integer(tempseed)
                     )
                     
    result1 <- list(inftime = tmp$tau, removaltime = tmp$remt)
  }
  
  # Calling fortran subroutine - Contact networks: Susceptible-Infectious (SI)
  if (!is.null(contact)) {
      if (length(contact)/(n * n) != ni) {
        stop('epidata:  Dimension of beta  and the number of contact networks are not matching')
      }
    network <- array(contact, c(n, n, ni))
  }
  if ((type == "SI") && !is.null(contact)) {
    tmp <- .Fortran("datacon",
                     n = as.integer(n),
                     tmin = as.integer(tmin),
                     tmax = as.integer(tmax),
                     ns = as.integer(ns),
                     ni = as.integer(ni),
                     alpha = as.numeric(alpha),
                     beta = as.numeric(beta),
                     spark = as.numeric(spark),
                     covmat = as.vector(covmat),
                     network = as.vector(network),
                     tau = as.integer(inftime),
                     tempseed = as.integer(tempseed)
                     )
                     
    result1 <- list(inftime = tmp$tau)
  }
 
 # Calling fortran subroutine - Contact networks: Susceptible-Infectious-Removed (SIR)
  if ((type == "SIR") && !is.null(contact)) {
    tmp <- .Fortran("dataconsir",
                     n = as.integer(n),
                     tmin = as.integer(tmin),
                     tmax = as.integer(tmax),
                     n = as.integer(ns),
                     ni = as.integer(ni),
                     lambda = as.integer(infperiod),
                     alpha = as.numeric(alpha),
                     beta = as.numeric(beta),
                     spark = as.numeric(spark),
                     covmat = as.vector(covmat),
                     network = as.vector(network),
                     tau = as.integer(inftime),
                     remt = as.integer(remt),
                     tempseed = as.integer(tempseed)
                     )
                    
    result1 <- list(inftime = tmp$tau, removaltime = tmp$remt)
  }
  return(result1)
  # End of function
}



