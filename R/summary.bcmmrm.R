#' Summarize a bcmmrm Object.
#'
#' Additional information about the Box-Cox transformed MMRM analysis
#' represented by \code{object} is extracted and included as components
#' of \code{object}.
#'
#' @param object an object inheriting from class "\code{bcmmrm}", representing
#' the Box-Cox transformed MMRM analysis.
#' @param robust an optional logical value used to specify whether to apply
#' the robust inference. Default is \code{TRUE}.
#' @param ssadjust an optional logical value used to specify whether to apply
#' the empirical small sample adjustment. Default is \code{TRUE}.
#' @param ... some methods for this generic require additional arguments.
#' None are used in this method.
#'
#' @return an object inheriting from class \code{summary.bcmmrm} with all
#' components included in \code{object} (see \code{\link{glsObject}} for
#' a full description of the components) plus the following components:
#' \describe{
#'   \item{\code{median}}{a list including inference results of the model median
#'         for specified values of \code{robust} and \code{ssadjust}.}
#'   \item{\code{meddif}}{a list including inference results of the model median
#'         difference for specified values of \code{robust} and \code{ssadjust}.}
#'   \item{\code{robust}}{a specified value of \code{robust}.}
#'   \item{\code{ssadjust}}{a specified value of \code{ssadjust}.}
#' }
#'
#' @seealso \code{\link{bcmmrm}}, \code{\link{bcmmrmObject}}, \code{\link{summary}}
#'
#' @examples
#'  data(aidscd4)
#'  resar <- bcmarg(cd4 ~ as.factor(treatment), aidscd4, weekc, id, "AR(1)")
#'  summary(resar)
#'
#' @export
summary.bcmmrm <- function(object, robust = T, ssadjust = T, ...) {
  if (robust & ssadjust){
    med <- object$median.rob.adj
    dif <- object$meddif.rob.adj
  }
  if (!robust & ssadjust){
    med <- object$median.mod.adj
    dif <- object$meddif.mod.adj
  }
  if (robust & !ssadjust){
    med <- object$median.rob
    dif <- object$meddif.rob
  }
  if (!robust & !ssadjust){
    med <- object$median.mod
    dif <- object$meddif.mod
  }
  structure(class = "summary.bcmmrm",
            c(object, list(median = med,
                 meddif = dif,
                 robust = robust,
                 ssadjust = ssadjust
                 )))
}

#' @export
print.summary.bcmmrm <-
  function(x, digits = 3, ...)
  {
    mCall <- x$call
    covstr <- mCall$structure
    if (is.null(covstr) & !is.null(mCall$id)){
      covstr <- "UN"
    }
    cat("Model median inference based on MMRM with Box-Cox transformation\n")
    cat("\nData and variable information:\n")
    cat("  Outcome:", deparse(mCall$outcome), "\n")
    cat("  Group:", deparse(mCall$group), "\n")
    cat("  Time:", deparse(mCall$time), "\n")
    cat("  ID:", deparse(mCall$id), "\n")
    cat("  Covariate(s):", deparse(mCall$covv), "\n")
    cat("  Data:", deparse( mCall$data ), "\n")
    cat("\nAnalysis information:\n")
    cat("  Covariance structure:", deparse(covstr), "\n")
    cat("  Robust inference:", x$robust, "\n")
    cat("  Empirical small sample adjustment:", x$ssadjust, "\n")
    cat("\nAnalysis results:\n")
    cat("  Estimated transformation parameter: ",
        format(x$lambda, digits = digits), "\n")


    nt <- length(x$median)
    ng <- nrow(x$median[[1]])
    group0 <- x$meddif[[1]]$group0
    glabel0 <- as.character(group0)
    group1 <- x$meddif[[1]]$group1
    glabel1 <- as.character(group1)
    for (g in 1:ng){
      if (sum(group0 == g) != 0){
        glabel0[group0 == g] <- as.character(x$group.tbl$label[g])
      }
      if (sum(group1 == g) != 0){
        glabel1[group1 == g] <- as.character(x$group.tbl$label[g])
      }
    }

    for (t in 1:nt){
      medt <- x$median[[t]]
      names(medt)[1] <- deparse(mCall$group)
      medt[, 1] <- x$group.tbl$label
      cat("\n", "\nModel median inferences for", deparse(x$call$time), "=" ,
          x$time.tbl$label[t], "\n","\n")
      print(medt, digits = digits)
    }
    for (t in 1:nt){
      dift <- x$meddif[[t]]
      names(dift)[1] <- paste(deparse(mCall$group), "1")
      names(dift)[2] <- paste(deparse(mCall$group), "0")
      dift[, 1] <- glabel1
      dift[, 2] <- glabel0
      cat("\n", "\nInferences of model median difference between groups (",
          names(dift)[1], "-", names(dift)[2], ") for", deparse(x$call$time),
          "=", x$time.tbl$label[t], "\n","\n")
      dift[, 8] <- round(dift[, 8], digits = digits)
      print(dift, digits = digits)
    }
    invisible(x)
  }
