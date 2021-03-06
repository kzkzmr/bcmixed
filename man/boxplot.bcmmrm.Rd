% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/boxplot.bcmmrm.R
\name{boxplot.bcmmrm}
\alias{boxplot.bcmmrm}
\title{Box-whisker plot for a bcmmrm Object.}
\usage{
\method{boxplot}{bcmmrm}(
  x,
  timepoint = NULL,
  xlab = NULL,
  ylab = NULL,
  main = TRUE,
  sub = NULL,
  verbose = FALSE,
  ...
)
}
\arguments{
\item{x}{an object inheriting from class "\code{bcmmrm}", representing
the Box-Cox transformed MMRM analysis.}

\item{timepoint}{an numeric value of a specified level of \code{time}
variable. The default is \code{NULL} and the last level is specified.}

\item{xlab}{a title for the x axis. The default is \code{NULL} and
the name of \code{time} or \code{group} variable is used.}

\item{ylab}{a title for the y axis. The default is \code{NULL} and
the default title is "Transformed residuals".}

\item{main}{a main title for the plot. The default is \code{TRUE} and
default title is "Box-whisker plot for transformed residuals".}

\item{sub}{a sub title for the plot. The default is \code{NULL}.}

\item{verbose}{a logical optional value specifying to print the detailed
plot information in the console. The default is \code{FALSE}.}

\item{...}{some methods for this generic require additional arguments.}
}
\value{
a box-whisker plot for transformed residual.
}
\description{
Box-whisker plot for the transformed residuals of each treatment groups
at a specified time point with error bar plot (mean +- SD)  using
\code{bcmmrmObject}.
}
\examples{
 data(aidscd4)
 lmd.bl <- bcreg(cd4.bl ~ 1, data = aidscd4[aidscd4$weekc == 8, ])$lambda
 aidscd4$cd4.bl.tr <- (aidscd4$cd4.bl ^ lmd.bl - 1) / lmd.bl
 resar <- bcmmrm(outcome = cd4, group = treatment, data = aidscd4,
                time = weekc, id = id, covv = c("cd4.bl.tr", "sex"),
                cfactor = c(0, 1), structure = "AR(1)", glabel =
                c("Zid/Did", "Zid+Zal", "Zid+Did", "Zid+Did+Nev"))
 boxplot(resar, xlab = "Treatment", col = 1:4)
}
\seealso{
\code{\link{bcmmrm}}, \code{\link{bcmmrmObject}}
}
