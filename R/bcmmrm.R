#' Model Median Inference for Longitudinal Data in Randomized Clinical Trials.
#'
#' \code{bcmmrm} provides inference results for the model median differences
#' between treatment groups proposed by Maruo et al. (2017), which focuses on
#' continuous and positive longitudinally observed outcomes and a situation
#' where the efficacy of some treatments is compared based on a randomized,
#' parallel group clinical trial. If \code{time} and \code{id} are not
#' specified, inference results reduce to the results for the context of
#' linear regression model provided by Maruo et al. (2015).
#'
#' @param outcome a name of positive outcome (dependent) variable included in
#'   \code{data}.
#' @param group  a name of treatment group variable included in \code{data}.
#' @param data a data frame that may include \code{outcome}, \code{group},
#'   \code{time}, \code{id}, and specified covariate variables.
#' @param time a name of time variable for repeated measurements included
#'   in \code{data}. The default is \code{NULL}.
#' @param id a name of subject id variable for repeated measurements included
#'   in \code{data}. The default is \code{NULL}.
#' @param covv a character vector with names of covariate variables included
#'   in \code{data}. The default is \code{NULL}.
#' @param cfactor an integer vector including nominal variable indicators for
#'   covariate variables. Nominal variable: \code{1}, continuous variable:
#'   \code{0}. The default is \code{NULL}.
#' @param structure specify the covariance structure from \code{c("UN", "CS",
#'   "AR(1)")}. The default is \code{"UN"}.
#' @param conf.level a numeric value of the confidence level for the
#'   confidence intervals. The default is 0.95.
#' @param lmdint a vector containing the end-points of the interval to be
#'   searched for a transformation parameter. The default is \code{c(-3, 3)}.
#' @param glabel a vector of length number of treatment groups containing
#'   the labels of \code{group} variable. The default is \code{NULL} and the
#'   levels of \code{group} variable in \code{data} are used.
#' @param tlabel a vector of length number of repeated measures containing
#'   the labels of \code{time} variable. The default is \code{NULL} and the
#'   levels of \code{time} variable in \code{data} are used.
#'
#' @return an object of class "\code{bcmmrm}" representing the results of model
#'   median inference based on the Box-Cox transformed MMRM model.
#'   Generic functions such as \code{print}, \code{plot}, and \code{summary}
#'   have methods to show the results of the fit. See \code{\link{bcmmrmObject}}
#'   for the components of the fit.
#'
#' @note If baseline observation for the outcome variable is available, Box-Cox
#'   transformed baseline should be included as a covariate for accuracy of
#'   estimation.\cr Although this function can be applied to non-randomized
#'   trial data, performances of the above approach have not evaluated in
#'   context of non-randomized trials.
#'
#' @references \itemize{
#'   \item Maruo, K., Isogawa, N., Gosho, M. (2015). Inference of median
#'         difference based on the Box-Cox model in randomized clinical trials.
#'         \emph{Statistics in Medicine}, 34, 1634-1644.\cr
#'   \item Maruo, K., Yamaguchi, Y., Noma, H., Gosho, M. (2017). Interpretable
#'         inference on the mixed effect model with the Box-Cox transformation.
#'         \emph{Statistics in Medicine}, 36, 2420-2434.
#' }
#'
#' @seealso \code{\link{bcmarg}}, \code{\link{bcmmrmObject}}
#'
#' @examples
#' data(aidscd4)
#'
#' # covariate: Box-Cox transformed baseline (continuous) and sex (nominal),
#' # covariance structure: AR(1) structure
#' # *Note: The UN structure is preferred although the AR(1)
#' #        structure is used in this example to reduce calculation time
#'
#' # Box-Cox transformation for the baseline
#' lmd.bl <- bcmarg(cd4.bl ~ 1, data = aidscd4[aidscd4$weekc == 8, ])$lambda
#' aidscd4$cd4.bl.tr <- bct(aidscd4$cd4.bl, lmd.bl)
#'
#' # Median inference for each group and week
#' bcmmrm(outcome = cd4, group = treatment, data = aidscd4, time = weekc,
#'        id = id, covv = c("cd4.bl.tr", "sex"), cfactor = c(0, 1),
#'        structure = "AR(1)", glabel = c("Zid/Did", "Zid+Zal", "Zid+Did",
#'        "Zid+Did+Nev"))
#'
#' @importFrom stats model.matrix pnorm pt qnorm qt
#'
#' @export
bcmmrm <- function(outcome, group, data, time = NULL, id = NULL,
                   covv = NULL, cfactor = NULL, structure = "UN",
                   conf.level = 0.95, lmdint = c(-3, 3),
                   glabel = NULL, tlabel = NULL) {
  if (conf.level <= 0 | conf.level >= 1) {
    stop("conf.level must be within the range of c(0,1)")
  }
  Call <- match.call()
  y <- deparse(substitute(outcome))
  data0 <- data
  group <- deparse(substitute(group))
  data0$group <- data0[, group]
  data0$y <- data0[, y]
  covform <- ""
  cct <- 0
  if (!is.null(covv)) {
    cct <- numeric(length(covv)) + 1
    for (i in seq_len(length(covv))) {
      if (cfactor[i] == 0) {
        covform <- paste(covform, "+", covv[i])
      } else {
        covform <- paste(covform, "+as.factor(", covv[i], ")")
        cct[i] <- length(unique(data0[, covv[i]])) - 1
      }
    }
  }
  ng <- length(unique(data0$group))
  if (is.null(glabel)) {
    glabel <- names(table(as.factor(data0$group)))
  }
  group.tbl <- data.frame(code = 1:ng, label = glabel)

  if (deparse(substitute(time)) != "NULL" &
      deparse(substitute(id)) != "NULL") {
    if (is.null(Call$structure)) {
      Call$structure <- "UN"
    }
    time <- deparse(substitute(time))
    id <- deparse(substitute(id))
    data0$time0 <- data0[, time]
    time.tbl <- sort(unique(data0$time0))
    for (i in seq_len(length(time.tbl))) {
      data0[data0$time0 == time.tbl[i], "time"] <- i
    }
    data0$id <- data0[, id]
    formula <- formula(paste("y ~ as.factor(group) + as.factor(time) +
                             as.factor(group):as.factor(time)", covform))
    nt <- length(unique(data0$time))
    if (is.null(tlabel)) {
      tlabel <- time.tbl
    }
    time.tbl <- data.frame(code = 1:nt, label = tlabel)
    msflg <- table(data0$id, is.na(data0$y))[, 1]
    N <- sum(msflg != 0)
    omis <- names(msflg)[msflg == 0]
    if (length(omis) > 0L) {
      for (i in seq_len(length(omis))) {
        data0 <- data0[data0$id != omis[i], ]
      }
    }
    rm(id)
    rm(time)
    try1 <- bcmarg(formula, data0, time, id, structure, lmdint)
    fitted <- try1$glsObject$fitted
  } else {
    formula <- formula(paste("y ~ as.factor(group)", covform))
    data0$time <- 1
    nt <- 1
    try1 <- bcmarg(formula, data0)
    time.tbl <- data.frame(code = 1, label = 1)
    fitted <- try1$glsObject$fitted.values
  }

  le <- try1$lambda
  lik <- try1$logLik
  beta <- try1$beta
  alp <- try1$alpha
  V <- try1$V
  RR <- V / (sqrt(diag(V)) %*% t(sqrt(diag(V))))
  nc <- sum(cct)
  nb <- length(beta)
  ns <- length(alp)
  iII <- try1$Vtheta.mod
  iIIr <- try1$Vtheta.rob
  options(na.action = "na.pass")
  X <- model.matrix(formula, data = data0)
  if (nt != 1) {
    X <- X[!duplicated(data0$id), ]
  }
  dimnames(X) <- NULL
  xcm <- 0
  if (!is.null(covv)) {
    if (nc == 1)  {
      xcm <- mean(X[, ng + nt])
    } else {
      xcm <- c(colMeans(X[, (ng + nt):(ng + nt + nc - 1)]))
    }
  }
  adj.prm <- try1$adj.prm
  N <- adj.prm[1]
  Ncc <- adj.prm[2]
  n.data <- adj.prm[3]
  n.na <- adj.prm[4]
  if (structure == "UN" & nt != 1) {
    nu <- Ncc - nt
    nu2 <- sqrt(Ncc / nu)
  } else if (nt == 1) {
    nu <- n.data - nb
    nu2 <- sqrt(n.data / nu)
  } else{
    nu <- (N - ng) * (nt - 1) - n.na
    nu2 <- sqrt(n.data / (n.data - nb))
  }
  cbn <- choose(ng, 2)
  comb <- matrix(0, cbn, 2)
  count <- 1
  for (ii in 1:ng) {
    for (jj in 1:ng) {
      if (ii < jj) {
        comb[count, 1] <- ii
        comb[count, 2] <- jj
        count <- count + 1
      }
    }
  }
  median.mod <- list()
  median.rob <- list()
  median.mod.adj <- list()
  median.rob.adj <- list()
  meddif.mod <- list()
  meddif.rob <- list()
  meddif.mod.adj <- list()
  meddif.rob.adj <- list()
  for (tp0 in 1:nt) {
    xi <- numeric(ng)
    bt <- 0
    dbt <- numeric(nt - 1)
    if (tp0 != 1) {
      bt <- beta[ng + tp0 - 1]
      dbt[tp0 - 1] <- 1
    }
    bc <- 0
    if (!is.null(covv)) {
      bc <- beta[(ng + nt):(ng + nt + nc - 1)]
    }
    dtv <- numeric(ns)
    SExi <- numeric(ng)
    SExir <- numeric(ng)
    Dt <- matrix(0, 1 + nb + ns, ng)
    for (i in 1:ng) {
      bg <- 0
      dbg <- numeric(ng - 1)
      if (i != 1) {
        bg <- beta[i]
        dbg[i - 1] <- 1
      }
      bgt <- 0
      dbgt <- numeric((ng - 1) * (nt - 1))
      if (i != 1 & tp0 != 1) {
        bgti <- (tp0 - 2) * (ng - 1) + i - 1
        bgt <- beta[ng + nt + nc - 1 + bgti]
        dbgt[bgti] <- 1
      }
      xi[i] <- (le * (beta[1] + bg + bt + bgt + sum(xcm * bc)) + 1) ^ (1 / le)
      dtl <- le ^ (-2) * xi[i] * (1 - le * log(xi[i]) - xi[i] ^ (-le))
      if (!is.null(covv)) {
        dtb <- xi[i] ^ (1 - le) * c(1, dbg, dbt, xcm, dbgt)
      } else {
        dtb <- xi[i] ^ (1 - le) * c(1, dbg, dbt, dbgt)
      }
      Dt[, i] <- c(dtl, dtb, dtv)
      SExi[i] <- sqrt(c(t(Dt[, i]) %*% iII %*% Dt[, i]))
      SExir[i] <- sqrt(c(t(Dt[, i]) %*% iIIr %*% Dt[, i]))
    }
    delta <- numeric(cbn)
    SEd <- numeric(cbn)
    SEdr <- numeric(cbn)
    for (i in 1:cbn) {
      delta[i] <- xi[comb[i, 2]] - xi[comb[i, 1]]
      SEd[i] <- sqrt(c(t(Dt[, comb[i, 2]] - Dt[, comb[i, 1]]) %*% iII %*%
                         (Dt[, comb[i, 2]] - Dt[, comb[i, 1]])))
      SEdr[i] <- sqrt(c(t(Dt[, comb[i, 2]] - Dt[, comb[i, 1]]) %*% iIIr %*%
                          (Dt[, comb[i, 2]] - Dt[, comb[i, 1]])))
    }
    SExi2 <- SExi * nu2
    SExir2 <- SExir * nu2
    SEd2 <- SEd * nu2
    SEdr2 <- SEdr * nu2
    sig.level <- 1 - conf.level
    tt <- qnorm(1 - sig.level / 2)
    tvalue <- delta / SEd
    pvalue <- (1 - pnorm(abs(tvalue))) * 2
    tvaluer <- delta / SEdr
    pvaluer <- (1 - pnorm(abs(tvaluer))) * 2
    tvalue2 <- delta / SEd2
    tvaluer2 <- delta / SEdr2
    tt2 <- qt(1 - sig.level / 2, nu)
    pvalue2 <- (1 - pt(abs(tvalue2), nu)) * 2
    pvaluer2 <- (1 - pt(abs(tvaluer2), nu)) * 2
    lowerm <- xi - SExi * tt
    upperm <- xi + SExi * tt
    lowermr <- xi - SExir * tt
    uppermr <- xi + SExir * tt
    lowerm2 <- xi - SExi2 * tt2
    upperm2 <- xi + SExi2 * tt2
    lowermr2 <- xi - SExir2 * tt2
    uppermr2 <- xi + SExir2 * tt2
    lowerd <- delta - SEd * tt
    upperd <- delta + SEd * tt
    lowerdr <- delta - SEdr * tt
    upperdr <- delta + SEdr * tt
    lowerd2 <- delta - SEd2 * tt2
    upperd2 <- delta + SEd2 * tt2
    lowerdr2 <- delta - SEdr2 * tt2
    upperdr2 <- delta + SEdr2 * tt2
    median.mod[[tp0]] <- data.frame(group = 1:ng, median = xi, SE = SExi,
                                    lower.CL = lowerm, upper.CL = upperm)
    median.rob[[tp0]] <- data.frame(group = 1:ng, median = xi, SE = SExir,
                                    lower.CL = lowermr, upper.CL = uppermr)
    median.mod.adj[[tp0]] <- data.frame(group = 1:ng, median = xi, SE = SExi2,
                                        lower.CL = lowerm2, upper.CL = upperm2)
    median.rob.adj[[tp0]] <- data.frame(group = 1:ng, median = xi, SE = SExir2,
                                        lower.CL = lowermr2,
                                        upper.CL = uppermr2)
    meddif.mod[[tp0]] <- data.frame(group1 = comb[, 2], group0 = comb[, 1],
                                    delta = delta, SE = SEd, lower.CL = lowerd,
                                    upper.CL = upperd, t.value = tvalue,
                                    p.value = pvalue)
    meddif.rob[[tp0]] <- data.frame(group1 = comb[, 2], group0 = comb[, 1],
                                    delta = delta, SE = SEdr,
                                    lower.CL = lowerdr, upper.CL = upperdr,
                                    t.value = tvaluer, p.value = pvaluer)
    meddif.mod.adj[[tp0]] <- data.frame(group1 = comb[, 2], group0 = comb[, 1],
                                        delta = delta, SE = SEd2,
                                        lower.CL = lowerd2, upper.CL = upperd2,
                                        t.value = tvalue2, p.value = pvalue2)
    meddif.rob.adj[[tp0]] <- data.frame(group1 = comb[, 2], group0 = comb[, 1],
                                        delta = delta, SE = SEdr2,
                                        lower.CL = lowerdr2,
                                        upper.CL = upperdr2,
                                        t.value = tvaluer2, p.value = pvaluer2)
  }
  data$ytr <- (data0$y ^ le - 1) / le
  data$ytr.fitted <- data$ytr
  data$ytr.fitted[!is.na(data$ytr.fitted)] <- fitted
  data$res.tr <- data$ytr - data$ytr.fitted

  structure(class = "bcmmrm",
            list(call = Call, median.mod = median.mod,
                 median.rob = median.rob, median.mod.adj = median.mod.adj,
                 median.rob.adj = median.rob.adj,
                 meddif.mod = meddif.mod, meddif.rob = meddif.rob,
                 meddif.mod.adj = meddif.mod.adj,
                 meddif.rob.adj = meddif.rob.adj, lambda = le, R = RR,
                 logLik = lik, betainf = try1$betainf, time.tbl = time.tbl,
                 group.tbl = group.tbl, inf.marg = try1, outdata = data,
                 conf.level = conf.level))
}

#' @export
logLik.bcmmrm <- function(object, REML = F, ...) {
  if (REML) {
    stop("REML method can not be used in bcmmrm.")
  }
  structure(object$logLik,
            df = 1 + length(c(object$inf.marg$beta, object$inf.marg$alpha)),
            class = "logLik")
}

#' @export
coef.bcmmrm <- function(object, ...) {
  coef(object$inf.marg$glsObject)
}

#' @export
print.bcmmrm <- function(x, ...) {
    mCall <- x$call
    covstr <- mCall$structure
    if (is.null(covstr) & !is.null(mCall$id)) {
      covstr <- "UN"
    }
    cat("Model median estimation based on MMRM with Box-Cox transformation\n")
    cat("  Outcome:", deparse(mCall$outcome), "\n")
    cat("  Group:", deparse(mCall$group), "\n")
    cat("  Time:", deparse(mCall$time), "\n")
    cat("  ID:", deparse(mCall$id), "\n")
    cat("  Covariate(s):", deparse(mCall$covv), "\n")
    cat("  Covariance structure:", deparse(covstr), "\n")
    cat("  Data:", deparse(mCall$data), "\n")
    cat("  Estimated transformation parameter: ",
        format(x$lambda, digits = 3), "\n")
    cat("  Log-likelihood: ", format(x$logLik), "\n", sep = "")

    nt <- length(x$median.mod)
    med <- c()
    for (i in 1:nt) {
      med <- cbind(med, x$median.mod[[i]]$median)
    }
    med <- cbind(x$group.tbl$label, as.data.frame(med))
    if (nt > 1) {
      names(med) <- c(paste(mCall$group, "|", mCall$time), x$time.tbl$label)
      cat("\nModel median estimates (row: group, col: time):\n")
      print(med, digits = 3, ...)
    } else {
      names(med) <- c(paste(mCall$group), "Estimate")
      cat("\nModel median estimates:\n")
      print(med, digits = 3, ...)
    }
    invisible(x)
}
