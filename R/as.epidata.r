################################################################################
# Part of the R/EpiILM package
#
# AUTHORS:
#         Waleed Almutiry <wkmtierie@qu.edu.sa>,
#         Vineetha Warriyar. K. V. <vineethawarriyar.kod@ucalgary.ca>, and
#         Rob Deardon <robert.deardon@ucalgary.ca>
#
# Free software under the terms of the GNU General Public License, version 2,
# a copy of which is available at http://www.r-project.org/Licenses/.
################################################################################

as.epidata <- function(type, n, x = NULL, y = NULL, inftime, infperiod = NULL, contact = NULL) {

  # Error checks for input arguments
  if (is.null(type) | !(type %in% c("SI", "SIR"))) {
       stop("as.epidata: Specify type as \"SI\" or \"SIR\".", call. = FALSE)
  }

  if (is.null(n)) {
       stop("as.epidata: The number of individuals \"n\" has to be specified.", call. = FALSE)
  }

  if (is.null(contact) &  (is.null(x) | is.null(y))) {
      stop('as.epidata: Specify either contact network or x, y coordinates.')
  }

  if (!is.null(contact) &  (!is.null(x) | !is.null(y))) {
      stop('as.epidata: Specify either contact network or x, y coordinates.')
  }

  if (!is.null(x) & !is.null(y) ) {
    if ((length(y) != n) | (length(x) != n)) {
      stop('as.epidata: Length of x or y is not compatible.')
    }
    XYcoordinates <- cbind(x,y)
    contact <- NULL
  }

  if (!is.null(contact)) {
      if (is.matrix(contact)) {
        if (length(contact)/(n*n) != 1) {
          stop('as.epidata:  The contact network matrix is not an n by n square matrix.')
        }
      } else if (is.array(contact)) {
        if (length(contact)/(n*n) != dim(contact)[3]) {
          stop('as.epidata:  One or all of the contact network matrix are not an n by n square matrix.')
        }
      } else {
        stop('as.epidata:  The contact network must be specified as an n by n square matrix or an array of n by n square matrices.')
      }
    XYcoordinates <- NULL
  }

  if (!is.null(inftime)) {
      if ((length(inftime) != n)) {
        stop('as.epidata: Length of inftime is not compatible.')
      }
  } else {
      stop('as.epidata: The inftime has to be specified as a vector of length \"n\".')
  }

  if (type == "SIR") {
    if (is.null(infperiod)) {
      stop('as.epidata: The infectious period has to be specified as a vector of length \"n\" via the option \"infperiod\".')
    }
    if (!is.null(infperiod)) {
        if (!is.vector(infperiod)) {
          stop('as.epidata: The infectious period has to be specified as a vector of mode = \"integer\" and with length \"n\" via the option \"infperiod\".')
        }
        if (length(infperiod) != n) {
          stop('as.epidata: Length of the infectious period vector \"infperiod\" is not compatible')
        }
        remtime <- inftime + infperiod
    }
  } else if (type == "SI") {
    if (!is.null(infperiod)) {
      stop('as.epidata: There is conflict inputs between the option type =\"SI\" and \"infperiod\".')
    }
    remtime <- NULL
  }

  result1 <- list(type = type, XYcoordinates = XYcoordinates, contact = contact, inftime = inftime, remtime = remtime)

  class(result1) <- "epidata"

  result1
  # End of function
}
