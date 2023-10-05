#' Applying the Box-Cox Transformation to a numeric vector.
#'
#' \code{bct.v} returns the Box-Cox transformed numeric vector (Box and Cox,
#' 1964).
#'
#' @param y a positive real number vector.
#' @param lmdint a vector containing the end-points of the interval to be
#'   searched for a transformation parameter. Default is c(-3, 3).
#'
#' @return \code{bct.v} returns the Box-Cox transformed numeric vector and
#' estimated transformation parameter.
#' \describe{
#'   \item{\code{transformed}}{The Box-Cox transformed numeric vector.}
#'   \item{\code{lambda}}{a numeric value of the estimate of the transformation
#'         parameter.}
#' }
#'
#' @references Box, G.E.P. and Cox, D.R. (1964). An analysis of transformations
#' (with discussion).
#' \emph{Journals of the Royal Statistical Society, Series B}, 26,
#' 211-246, \doi{10.1111/j.2517-6161.1964.tb00553.x}.
#'
#' @examples
#'   y <- exp(rnorm(50))
#'   bct.v(y)
#'
#' @export
bct.v <- function(y, lmdint = c(-3, 3)) {
  if (!is.numeric(y)) {
    stop("The argument must be numeric.")
  }
  if (sum(y <= 0) != 0) {
    stop("The argument includes non-positive value(s)")
  }
  if (sum(is.na(y)) > 0) {
    y <- y[!is.na(y)]
    warning("Missing values are excluded.")
  }
  lambda <- bcreg(y ~ 1, data.frame(y = y), lmdint)$lambda
  z <- bct(y, lambda)
  return(list(transformed = z, lambda = lambda))
}
