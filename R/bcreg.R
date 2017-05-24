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
#'
#' @return LmeMargBc returns a list including following conponents:
#' \describe{
#'   \item{lambda}{estimate of a transformation parameter.}
#'   \item{beta}{estimate of a regression parameter vector.}
#'   \item{sigma}{estimate of a scale parameter.}
#'   \item{transformed}{inference results for beta under the
#'         assumption that lambda is known.}
#'   \item{lik}{maximized liklihood.}
#' }
#'
#' @references Box, G.E.P. and Cox, D.R. (1964). An analysis of transformations
#' (with discussion). Journals of the Royal Statistical Society, Series B, 26,
#' 211-246.
#'
#' @examples
#'   data(aidscd4)
#'   #Transformation of baseline observation for aid.cd4 data
#'   bcreg(cd4.bl ~ 1, aidscd4[aidscd4$weekc == 8, ])
#'
#' @importFrom MASS ginv
#'
#' @export
bcreg <- function(formula, data){
  formula <- formula(formula)
  data <- as.data.frame(data)
  y <- data[, as.character(formula)[2]]
  if (sum(y < 0, na.rm = T) > 0L){
    stop("outcome must be positive.")
  }
  options(na.action = "na.pass")
  X <- model.matrix(formula, data = data)
  if (sum(is.na(X)) > 0L){
    stop("There are missing values in explanatory variables.")
  }
  flg.na <- !is.na(y)
  y <- y[flg.na]
  X <- X[flg.na, ]
  n <- length(y)
  lbc <- function(l, y, X){
    if (l == 0) z <- log(y)
    else z <- (y ^ l - 1) / l
    n <- length(y)
    beta <- ginv(t(X) %*% X) %*% t(X) %*% z
    zm <- z - X %*% beta
    sgm <- sqrt(sum(zm ^ 2) / n)
    lik <- -n / 2 * (log(2 * pi) + 1) - n * log(sgm) + (l - 1) * sum(log(y))
    return(lik)
  }
  res <- optimize(lbc, c(-3, 3), y = y, X = X, maximum = T)
  lmd <- res$maximum
  z <- (y ^ lmd - 1) / lmd
  beta <- ginv(t(X) %*% X) %*% t(X) %*% z
  zm <- z - X %*% beta
  sgm <- sqrt(sum(zm ^ 2) / n)
  lik <- res$objective
  data <- as.data.frame(data[flg.na, ])
  data$z <- as.numeric(z)
  evars <- as.character(formula)[3]
  formula2 <- formula(paste("z~", evars))
  bcres <- summary(lm(formula2, data))$coefficients
  res <- list(lambda = lmd, beta = as.numeric(beta), sigma = sgm,
              transformed = bcres, lik = lik)
  return(res)
}
