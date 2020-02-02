#' Applying the Box-Cox Transformation.
#'
#' \code{bct} returns the Box-Cox transformed numeric vector (Box and Cox,
#' 1964).
#'
#' @param y a positive real number vector.
#' @param lambda a scalar transformation parameter.
#'
#' @return bct returns the Box-Cox transformed numeric vector,
#' \code{z = log(y)} for \code{lambda = 0},
#' \code{z = (y ^ lambda - 1) / lambda} for \code{lambda ne 1}.
#'
#' @references Box, G.E.P. and Cox, D.R. (1964). An analysis of transformations
#' (with discussion).
#' \emph{Journals of the Royal Statistical Society, Series B}, 26,
#' 211-246, \doi{10.1111/j.2517-6161.1964.tb00553.x}.
#'
#' @examples
#'   y <- exp(rnorm(10))
#'   z <- bct(y, 0) #log transformation
#'
#' @export
bct <- function(y, lambda) {
  if (sum(y < 0)) {
    stop("All valuesin y must be positive.")
  }
  if (lambda == 0) {
    z <- log(y)
  }
  else {
    z <- (y ^ lambda - 1) / lambda
  }
  return(z)
}
