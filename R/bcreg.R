#' Linear regression model with the Box-Cox Transformation.
#'
#' \code{bcreg} returns the maximum likelihood estimates for parameters of
#' the linear regression models with the Box-Cox transformation (Box and Cox,
#' 1964).
#'
#' @param formula a two-sided linear formula object describing the model, with
#'   the response on the left of a ~ operator and the terms, separated by +
#'   operators, on the right.
#' @param data a data frame in which to interpret the variables named in the
#'   formula.
#' @param lmdint a vector containing the end-points of the interval to be
#'   searched for a transformation parameter. Default is c(-3, 3).
#'
#' @return bcreg returns a list including following components:
#' \describe{
#'   \item{\code{lambda}}{a numeric value with the estimate of
#'         the transformation parameter.}
#'   \item{\code{beta}}{a vector with the estimates of the regression
#'         parameters.}
#'   \item{\code{sigma}}{a numeric value with the estimate of the scale
#'         parameter.}
#'   \item{\code{betainf}}{a data frame with inference results for \code{beta}
#'         under the assumption that \code{lambda} is known.}
#'   \item{\code{lik}}{a numeric value with the maximized likelihood.}
#'   \item{\code{lmObject}}{an object of "\code{lm}" containing the results of
#'         \code{lm} function on the transformed scale}
#' }
#'
#' @references Box, G.E.P. and Cox, D.R. (1964). An analysis of transformations
#' (with discussion).
#' \emph{Journals of the Royal Statistical Society, Series B}, 26,
#' 211-246, \doi{10.1111/j.2517-6161.1964.tb00553.x}.
#'
#' @seealso \code{\link{lm}}
#'
#' @examples
#'   data(aidscd4)
#'   #Transformation of baseline observation for aid.cd4 data
#'   bcreg(cd4.bl ~ 1, aidscd4[aidscd4$weekc == 8, ])
#'
#' @importFrom MASS ginv
#' @importFrom stats lm model.matrix optimize
#'
#' @export
bcreg <- function(formula, data, lmdint = c(-3, 3)) {
  formula <- formula(formula)
  data <- as.data.frame(data)
  y <- data[, as.character(formula)[2]]
  if (sum(y < 0, na.rm = TRUE) > 0L) {
    stop("outcome must be positive.")
  }
  options(na.action = "na.pass")
  X <- model.matrix(formula, data = data)
  if (sum(is.na(X)) > 0L) {
    stop("There are missing values in explanatory variables.")
  }
  flg.na <- !is.na(y)
  y <- y[flg.na]
  X <- X[flg.na, ]
  lbc <- function(l, y, X) {
    if (l == 0) z <- log(y)
    else z <- (y ^ l - 1) / l
    n <- length(y)
    beta <- ginv(t(X) %*% X) %*% t(X) %*% z
    zm <- z - X %*% beta
    sgm <- sqrt(sum(zm ^ 2) / n)
    lik <- -n / 2 * (log(2 * pi) + 1) - n * log(sgm) + (l - 1) * sum(log(y))
    return(lik)
  }
  res <- optimize(lbc, interval = lmdint, maximum = TRUE, y = y, X = X)
  lmd <- res$maximum
  z <- (y ^ lmd - 1) / lmd
  lik <- res$objective
  data <- as.data.frame(data[flg.na, ])
  data$z <- as.numeric(z)
  evars <- as.character(formula)[3]
  formula2 <- formula(paste("z~", evars))
  lmres <- lm(formula2, data)
  beta <- coef(lmres)
  sgm <- summary(lmres)$sigma
  names(sgm) <- "Residual"
  bcres <- as.data.frame(summary(lmres)$coefficients)
  list(lambda = lmd, beta = beta, sigma = sgm,
            betainf = bcres, logLik = lik, lmObject = lmres)
}
