################################################################################
# Part of the R/EpiILM package
#
# AUTHORS:
#         Vineetha Warriyar. K. V. <vineethawarriyar.kod@ucalgary.ca>,
#         Waleed Almutiry <wkmtierie@qu.edu.sa>, and
#         Rob Deardon <robert.deardon@ucalgary.ca>
#
# Free software under the terms of the GNU General Public License, version 2,
# a copy of which is available at http://www.r-project.org/Licenses/.
################################################################################

epiBR0 <- function(x = NULL, y = NULL, contact = NULL, sus.par,  trans.par = NULL, beta,
                  spark = NULL, infperiod, Sformula = NULL, Tformula = NULL, tmax, niter) {

  # Error checks for input arguments
  if (is.null(contact) &  (is.null(x) | is.null(y))) {
    stop('epiBR0: Specify contact network or x, y coordinates')
  }

  if (is.null(spark)) {
       spark <- 0
  }

  ns     <- length(sus.par)
  nt     <- length(trans.par)
  ni     <- length(beta)

  if (!is.null(x)) {
    n <- length(x)

    if ((length(y) != n) | (length(x) != n)) {
      stop('epiBR0: Length of x or y is not compatible')
    }
  }

  if (!is.null(contact)) {
      if (is.matrix(contact)) {
        if (dim(contact)[1] != dim(contact)[2]) {
          stop('epiBR0:  The contact network matrix is not a square matrix.')
        }
        n <- dim(contact)[1]
        network <- array(contact, c(n, n, ni))
      } else if (is.array(contact)) {
        if (dim(contact)[1] != dim(contact)[2] | length(contact)/(sqrt(dim(contact)[1]) * sqrt(dim(contact)[2])) != dim(contact)[3]) {
          stop('epiBR0:  One or all of the contact network matrix are not a square matrix.')
        }
        n <- sqrt(dim(contact)[1])
        network <- contact #array(contact, c(n, n, ni))

      } else {
        stop('epiBR0:  The contact network must be specified as an n by n square matrix or an array of n by n square matrices.')
      }
   }

  if (length(infperiod) != n) {
    stop('epiBR0: Length of infperiod is not compatible')
  }

  # formula for susceptibility function
  if (!is.null(Sformula)) {
    covmat.sus <- model.matrix(Sformula)

    if ((ncol(covmat.sus) == length(all.vars(Sformula))) & (ns != length(all.vars(Sformula)))) {
      stop('epiBR0: Check Sformula (no intercept term) and the dimension of sus.par')
    }

    if ((ncol(covmat.sus) > length(all.vars(Sformula))) & (ns != ncol(covmat.sus))) {
      stop('epiBR0: Check Sformula (intercept term) and the dimension of sus.par')
    }
  } else {
      if (ns == 1) {
        covmat.sus <- matrix(1.0, nrow = n, ncol = ns)
      }
      if (ns > 1) {
        stop('epiBR0: Please specify the susceptibility covariate')
      }
  }

  # formula for transmissibility function

  if (!is.null(Tformula)) {
    covmat.trans <- model.matrix(Tformula)

    if ((ncol(covmat.trans) == length(all.vars(Tformula))) & (nt != length(all.vars(Tformula)))) {
      stop('epiBR0: Check Tformula (no intercept term) and the dimension of trans.par')
    }
  } else {
      if (nt == 0) {
        covmat.trans <- matrix(1.0, nrow = n, ncol = 1)
      }
      if (nt > 0) {
        stop('epiBR0: Please specify the transmissibility covariate')
      }
  }

  val <- 0.0
  # Calling Fortran subroutines : contact network model
  if (!is.null(contact)) {
    tmp <- .Fortran("rconsir",
                    n = as.integer(n),
                    tmax = as.integer(tmax),
                    ns = as.integer(ns),
                    nt = as.integer(nt),
                    ni = as.integer(ni),
                    lambda = as.integer(infperiod),
                    suspar = as.numeric(sus.par),
                    transpar = as.numeric(trans.par),
                    beta = as.numeric(beta),
                    spark = as.numeric(spark),
                    covmatsus = as.vector(covmat.sus),
                    covmattrans = as.vector(covmat.trans),
                    network = as.vector(network),
                    sim = as.integer(niter),
                    val = as.double(val),
                    countinf = as.integer(rep(0,niter))
                    )

  result1 <- list(BasicR0 = tmp$val,
                  simulated_BR0 = tmp$countinf)
  } else {
    # Calling Fortran subroutines : spatial model
    tmp <- .Fortran("rxysir",
                     n = as.integer(n),
                     tmax = as.integer(tmax),
                     ns = as.integer(ns),
                     nt = as.integer(nt),
                     ni = as.integer(ni),
                     suspar = as.numeric(sus.par),
                     transpar = as.numeric(trans.par),
                     beta = as.numeric(beta),
                     spark = as.numeric(spark),
                     covmatsus = as.vector(covmat.sus),
                     covmattrans = as.vector(covmat.trans),
                     lambda = as.integer(infperiod),
                     x = as.double(x),
                     y = as.double(y),
                     sim = as.integer(niter),
                     val = as.double(val),
                     countinf = as.integer(rep(0,niter))
                     )

  result1 <- list(BasicR0 = tmp$val,
                  simulated_BR0 = tmp$countinf)
  }
  return(result1)
  }
