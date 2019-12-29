#' Plot a bcmmrm Object.
#'
#' Plot for the model medians of each treatment groups with the 95 percent
#' confidence intervals stored in \code{bcmmrmObject}.
#'
#' @param x an object inheriting from class "\code{bcmmrm}", representing
#' the Box-Cox transformed MMRM analysis.
#' @param robust an optional logical value used to specify whether to apply
#' the robust inference. Default is \code{TRUE}.
#' @param ssadjust an optional logical value used to specify whether to apply
#' the empirical small sample adjustment. Default is \code{TRUE}.
#' @param tnom a optional logical value indicating the scale of x axis
#' of the longitudinal median plot. When \code{tnom} is \code{TRUE}, nominal
#' scale is used and widths between any combinations of neighbor time points
#' are same. When \code{tnom} is \code{FALSE}, actual scale of \code{time}
#' variable is used. Default is \code{TRUE}.
#' @param dt an numeric value indicating shift length between groups
#'  in the longitudinal median plot. A multiplying factor for the default
#'  settings specified (e.g. if \code{2} is specified,
#'  shift length is twice longer than that for the default setting).
#'  Default is \code{1}.
#' @param timepoint an numeric value of a specified level of \code{time}
#' variable at which median plot is created. When \code{timepoint} is
#' \code{NULL} and number of time points is not 1, longitudinal median plot
#' (x axis is \code{time}) is created. Otherwise, median plot where
#' x axis is \code{group} is created. Default is \code{NULL}.
#' @param xlab a title for the x axis. Default is \code{NULL} and
#'  the name of \code{time} or \code{group} variable is used.
#' @param ylab a title for the y axis. Default is \code{NULL} and
#'  the name of \code{outcome} variable is used.
#' @param xlim a numeric vector with length of 2 indicating limits of x
#'  axis. Default is \code{NULL} and limits are calculated automatically.
#' @param ylim a numeric vector with length of 2 indicating limits of y
#'  axis. Default is \code{NULL} and limits are calculated automatically.
#' @param lwd an optional positive numeric value indicating line width.
#' Default is \code{2}.
#' @param col an integer or a character vector with length of the number of
#'  groups indicating colors of lines for each treatment group. Default is
#'  \code{NULL} and all of colors are black.
#' @param lty an optional integer or a character vector with length of the
#'  number of groups indicating line types of lines for each treatment group.
#'  Default is \code{NULL} and \code{1:ng} is used, where ng is number of
#'  groups.
#' @param main a main title for the plot. Default is \code{TRUE} and
#'  default title is "(Longitudunal) Plot for median of each group".
#' @param sub a sub title for the plot. Default is \code{NULL}.
#' @param legend a logical optional value specifying to add legends to
#'  plots. When \code{legend} is \code{TRUE} legends are added to the plot.
#'  Otherwise, legends are not added.  Default is \code{TRUE}.
#' @param loc a character value indicating the location of the legends.
#' The location can be specified by setting \code{loc} to a single keyword
#' from the list \code{"bottomright"}, \code{"bottom"}, \code{"bottomleft"},
#' \code{"left"}, \code{"topleft"}, \code{"top"}, \code{"topright"},
#' \code{"right"}, and \code{"center"}. Default is \code{"topright"}.
#' @param ... some methods for this generic require additional arguments.
#'
#' @return a median plot.
#'
#' @seealso \code{\link{bcmmrm}}, \code{\link{bcmmrmObject}}
#'
#' @examples
#'  data(aidscd4)
#'  resar <- bcmmrm(outcome = cd4, group = treatment, data = aidscd4,
#'                 time = weekc, id = id, structure = "AR(1)", glabel =
#'                 c("Zid/Did", "Zid+Zal", "Zid+Did", "Zid+Did+Nev"))
#'  plot(resar, xlab = "Week", ylab = "CD4+1", col = 1:4, main = NULL)
#'  plot(resar, timepoint = 32, xlab = "Treatment", ylab = "CD4+1")
#' @importFrom graphics plot grid lines arrows axis legend
#' @export
plot.bcmmrm <- function(x, robust = T, ssadjust = T, dt = 1,
                        timepoint = NULL, tnom = T, xlab = NULL,
                        ylab = NULL, xlim = NULL, ylim = NULL,
                        lwd = 2, col = NULL, lty = NULL,
                        main = T, sub = NULL, legend = T, loc = "topright",
                        ...) {
  x <- summary(x, robust = robust, ssadjust = ssadjust)
  med <- x$median
  if (is.null(ylab)) {
    ylab <- deparse(x$call$outcome)
  }
  ng <- nrow(med[[1]])
  if (is.null(col)) {
    col <- rep("black", ng)
  }
  glabel <- x$group.tbl$label
  if (length(med) == 1 | !is.null(timepoint)) {
    if (!is.null(timepoint)) {
      med <- med[[x$time.tbl$code[x$time.tbl$label == timepoint]]]
    } else {
      med <- med[[1]]
    }

    Group <- 1:ng
    if (is.null(ylim)) {
      upl <- max(med$upper.CL) * 1.2
      lol <- min(med$lower.CL) * 0.8
      ylim <- c(lol, upl)
    }
    if (is.null(xlim)) {
      xlim <- c(0.5, ng + 0.5)
    }

    if (is.null(xlab)) {
      xlab <- deparse(x$call$group)
    }
    if (!is.null(main)) {
      if (main == T) {
        main <- "Plot for model median of each group"
      }
    }

    plot(1, 1, xlim = xlim, ylim = ylim, type = "n", ylab = ylab, xlab = xlab,
         xaxt = "n",
         panel.first = grid(NA, NULL, lty = 2, col = "#E9DECA", ...),
         main = main, sub = sub)
    axis(1, at = Group, labels = glabel)

    points(Group, med$median, lwd = lwd, col = col)
    arrows(Group, med$upper.CL, Group, med$lower.CL,
           code = 3, angle = 90, length = .07, lwd = lwd, col = col)
    cat("Analysis information:\n")
    cat("  Covariance structure:", deparse(x$call$structure), "\n")
    cat("  Robust inference:", x$robust, "\n")
    cat("  Empirical small sample adjustment:", x$ssadjust, "\n")
    cat("\nError bar: 95% confidence interval\n")
    cat("\nTimepoint:", deparse(x$call$time), "=", timepoint)
  } else {
    Time2 <- x$time.tbl$label
    Time1 <- x$time.tbl$code
    nt <- max(Time1)
    if (tnom == 1) {
      Time0 <- Time1
      hh <- 1
    } else {
      Time0 <- Time2
      hh <- (Time2[nt] - Time2[1]) / nt
    }
    Time <- list()
    est <- list()
    upper <- list()
    lower <- list()
    for (i in 1:ng) {
      Time[[i]] <- Time0 + 0.05 * dt * hh * (2 * i - ng - 1)
      est0 <- c()
      upper0 <- c()
      lower0 <- c()
      for (j in 1:nt) {
        est0 <- c(est0, med[[j]]$median[i])
        upper0 <- c(upper0, med[[j]]$upper.CL[i])
        lower0 <- c(lower0, med[[j]]$lower.CL[i])
      }
      est[[i]] <- est0
      upper[[i]] <- upper0
      lower[[i]] <- lower0
    }
    up <- c()
    lo <- c()
    for (i in 1:nt) {
      up <- c(up, med[[i]]$upper.CL)
      lo <- c(lo, med[[i]]$lower.CL)
    }
    if (is.null(ylim)) {
      upl <- max(up) * 1.2
      lol <- min(lo) * 0.8
      ylim <- c(lol, upl)
    }
    if (is.null(xlim)) {
      if (tnom == 1) {
        xlim <- c(0.5, nt + 0.5)
      } else {
        lx <- min(Time2)
        ux <- max(Time2)
        dux <- lx  * 0.1
        xlim <- c(min(c(lx - dux, lx - 1)), max(c(ux + dux, ux + 1)))
      }
    }
    if (is.null(xlab)) {
      xlab <- deparse(x$call$time)
    }
    if (!is.null(main)) {
      if (main == T) {
        main <- "Longitudinal plot for model median of each group"
      }
    }
    plot(1, 1, xlim = xlim, ylim = ylim, type = "n", ylab = ylab, xlab = xlab,
         xaxt = "n",  panel.first = grid(NA, NULL, lty = 2, col = "#E9DECA"),
         main = main, sub = sub, ...)
    if (tnom == 1) {
      axis(1, at = Time0, labels = Time2)
    } else {
      axis(1, at = Time0)
    }
    if (is.null(lty)) {
      lty <- 1:ng
    }
    for (i in 1:ng) {
      lines(Time[[i]], est[[i]], lty = lty[i], lwd = lwd, col = col[i])
      arrows(Time[[i]], upper[[i]], Time[[i]], lower[[i]],
             code = 3, angle = 90, length = .07 * dt, lwd = lwd,
             col = col[i])
    }
    if (legend) {
      legend(loc, legend = glabel, lty = lty, col = col, lwd = lwd)
    }
    cat("Analysis information:\n")
    cat("  Covariance structure:", deparse(x$call$structure), "\n")
    cat("  Robust inference:", x$robust, "\n")
    cat("  Empirical small sample adjustment:", x$ssadjust, "\n")
    cat("\nError bar: 95% confidence interval\n")
  }
}
