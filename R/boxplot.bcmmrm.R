#' Box-whisker plot for a bcmmrm Object.
#'
#' Box-whisker plot for the transformed residuals of each treatment groups
#' at a specified time point with error bar plot (mean +- SD)  using
#' \code{bcmmrmObject}.
#'
#' @param x an object inheriting from class "\code{bcmmrm}", representing
#' the Box-Cox transformed MMRM analysis.
#' @param timepoint an numeric value of a specified level of \code{time}
#' variable. The default is \code{NULL} and the last level is specified.
#' @param xlab a title for the x axis. The default is \code{NULL} and
#'  the name of \code{time} or \code{group} variable is used.
#' @param ylab a title for the y axis. The default is \code{NULL} and
#'  the default title is "Transformed residuals".
#' @param main a main title for the plot. The default is \code{TRUE} and
#'  default title is "Box-whisker plot for transformed residuals".
#' @param sub a sub title for the plot. The default is \code{NULL}.
#' @param  verbose a logical optional value specifying to print the detailed
#' plot information in the console. The default is \code{FALSE}.
#' @param ... some methods for this generic require additional arguments.
#'
#' @return a box-whisker plot for transformed residual.
#'
#' @seealso \code{\link{bcmmrm}}, \code{\link{bcmmrmObject}}
#'
#' @examples
#'  data(aidscd4)
#'  lmd.bl <- bcreg(cd4.bl ~ 1, data = aidscd4[aidscd4$weekc == 8, ])$lambda
#'  aidscd4$cd4.bl.tr <- (aidscd4$cd4.bl ^ lmd.bl - 1) / lmd.bl
#'  resar <- bcmmrm(outcome = cd4, group = treatment, data = aidscd4,
#'                 time = weekc, id = id, covv = c("cd4.bl.tr", "sex"),
#'                 cfactor = c(0, 1), structure = "AR(1)", glabel =
#'                 c("Zid/Did", "Zid+Zal", "Zid+Did", "Zid+Did+Nev"))
#'  boxplot(resar, xlab = "Treatment", col = 1:4)
#' @importFrom graphics boxplot axis abline points arrows lines
#' @importFrom stats sd
#' @export
boxplot.bcmmrm <- function(x, timepoint = NULL, xlab = NULL,
                           ylab = NULL, main = TRUE, sub = NULL,
                           verbose = FALSE, ...) {
  if (is.null(timepoint)) {
    nt <- length(x$median.mod)
    timepoint <- x$time.tbl$label[nt]
  }
  dgdat <- x$outdata[x$outdata[, deparse(x$call$time)] ==
                        timepoint, ]
  rgm <- ceiling(max(abs(dgdat$res.tr), na.rm = TRUE))
  group <- deparse(x$call$group)
  if (is.null(xlab)) {
    xlab <- group
  }
  if (is.null(ylab)) {
    ylab <- "Transformed residual"
  }
  if (!is.null(main)) {
    if (main) {
      main <- paste("Box-whisker plot for transformed residuals at",
                    deparse(x$call$time), "=", timepoint)
    }
  }
  bp <- boxplot(dgdat$res.tr ~ dgdat[, group],
          ylim = c(-rgm, rgm), main = main, sub = sub,
          ylab = ylab, xlab = xlab, xaxt = "n")
  glabel <- x$group.tbl$label
  axis(1, at = x$group.tbl$code, labels = glabel)
  abline(h = 0)
  xi <- 0.2 + seq(bp$n)
  mn.t <- tapply(dgdat$res.tr, dgdat[, group], mean, na.rm = TRUE)
  sd.t <- tapply(dgdat$res.tr, dgdat[, group], sd, na.rm = TRUE)
  points(xi, mn.t, pch = 4)
  arrows(xi, mn.t - sd.t, xi, mn.t + sd.t,
         code = 3, angle = 75, length = .1)
  if (verbose) {
    cat("Plot information:\n")
    cat("\nError bar plot: mean +- SD\n")
    cat("\nTimepoint:", deparse(x$call$time), "=", timepoint)
  }
}
