#' Marginal Model of the Mixed Effect Model with the Box-Cox Transformation.
#'
#' \code{bcmarg} returns the inference results the parameters of
#'   the marginal model of the linear mixed effect model with the Box-Cox
#'   transformation proposed by Maruo et al. (2017). If time and id are not
#'   specified, inference results reduce to the results for the context of
#'   linear regression model provided by Maruo et al. (2015).
#'
#' @param formula a two-sided linear formula object describing the model, with
#'   the response on the left of a ~ operator and the terms, separated by +
#'   operators, on the right.
#' @param data a data frame containing the variables used in the model.
#' @param time time variable name for repeated measurements. The default is
#'   NULL.
#' @param id subject id variable name for repeated measurements. The default is
#'   NULL.
#' @param structure specify the covariance structure from c("UN", "CS",
#'   "AR(1)"). The default is "UN".
#' @param lmdint a vector containing the end-points of the interval to be
#'   searched for a transformation parameter. The default is c(-3, 3).
#'
#' @return an object of class "\code{bcmarg}". Objects of this class have
#' methods for the generic functions  \code{coef}, \code{logLik},
#' \code{print}, and \code{summary}. The object includes following components
#' for the marginal model parameter inference:
#' \describe{
#'   \item{\code{lambda}}{a numeric value of the estimate of the transformation
#'         parameter.}
#'   \item{\code{beta}}{a vector with the estimates of the regression
#'         parameters.}
#'   \item{\code{alpha}}{a vector with the estimates of the covariance
#'         parameters.}
#'   \item{\code{V}}{variance-covariance matrix for any subject with no missing
#'         values.}
#'   \item{\code{betainf}}{a matrix containing the inference results for
#'        \code{beta} under the assumption that lambda is known.
#'        Note that standard errors might be underestimated
#'        although statistical tests would be asymptotically valid.}
#'   \item{\code{Vtheta.mod}}{model-based variance-covariance matrix for MLE of
#'         the vector of all parameters: \code{c(lambda, beta, alpha)}.}
#'   \item{\code{Vtheta.rob}}{robust variance-covariance matrix for MLE of
#'         the vector of all parameters.}
#'   \item{\code{logLik}}{a numeric value of the maximized likelihood.}
#'   \item{\code{adj.prm}}{a vector with parameters used for the empirical
#'         small sample adjustment in \code{\link{bcmmrm}}:
#'         c(number of subjects, number of completed subjects, number of
#'         outcome observations, number of missing observations).}
#'   \item{\code{glsObject}}{an object of "\code{gls}" (or "\code{lm}" when
#'         \code{time} and \code{id} are not specified) containing results of
#'         \code{gls} (or \code{lm}) function on the transformed scale. }
#' }
#'
#' @references \itemize{
#'   \item Maruo, K., Isogawa, N., Gosho, M. (2015). Inference of median
#'         difference based on the Box-Cox model in randomized clinical trials.
#'         \emph{Statistics in Medicine}, 34, 1634-1644,
#'         \doi{10.1002/sim.6408}.\cr
#'   \item Maruo, K., Yamaguchi, Y., Noma, H., Gosho, M. (2017). Interpretable
#'         inference on the mixed effect model with the Box-Cox transformation.
#'         \emph{Statistics in Medicine}, 36, 2420-2434,
#'         \doi{10.1002/sim.7279}.
#' }
#'
#' @seealso \code{\link{bcmmrm}} \code{\link[nlme]{gls}}
#'
#' @examples
#'  data(aidscd4)
#'  bcmarg(cd4 ~ as.factor(treatment) * as.factor(weekc) + age,
#'         data = aidscd4, time = weekc, id = id, structure = "AR(1)")
#'
#' @importFrom nlme gls glsControl corSymm varIdent corCompSymm corAR1
#'             corMatrix
#' @importFrom MASS ginv
#' @importFrom stats coef ftable model.matrix na.omit optimize xtabs predict
#'
#' @export
bcmarg <- function(formula, data, time = NULL, id = NULL, structure = "UN",
                   lmdint = c(-3, 3)) {
  backup_options <- options()
  Call <- match.call()
  tr <- function(X) sum(diag(X))
  covcalcUN <- function(glsob, nt) {
    flg <- 0
    j <- 0
    while (flg == 0) {
      j <- j + 1
      corm <- corMatrix(glsob$modelStruct$corStruct)[[j]]
      if (nrow(corm) == nt) flg <- 1
    }
    varstruct <- glsob$modelStruct$varStruct
    varests <- coef(varstruct, uncons = FALSE, allCoef = TRUE)
    covm <- corm * glsob$sigma ^ 2 * t(t(varests)) %*% t(varests)
    return(covm)
  }
  covcalcCSAR <- function(glsob, nt) {
    flg <- 0
    j <- 0
    while (flg == 0) {
      j <- j + 1
      corm <- corMatrix(glsob$modelStruct$corStruct)[[j]]
      if (nrow(corm) == nt) flg <- 1
    }
    covm <- corm * glsob$sigma ^ 2
    return(covm)
  }
  dec2bin <- function(num, digit=0) {
    if (num <= 0 && digit <= 0) {
      return(NULL)
    }else{
      return(append(Recall(num %/% 2, digit - 1), num %% 2))
    }
  }
  likn <- function(formula, data, structure, lmd) {
    y <- data$y
    if (lmd == 0) data$z <- log(y)
    else data$z <- (y ^ lmd - 1) / lmd
    if (structure == "UN") {
      RS <- gls(model = formula, data = data,
                correlation = corSymm(form = ~ time | id),
                weights = varIdent(form = ~ 1 | time), method = "ML",
                control = glsControl(msMaxIter = 100), na.action = na.omit)
    }
    if (structure == "CS") {
      RS <- gls(model = formula, data = data,
                correlation = corCompSymm(form = ~ time | id),
                method = "ML", control = glsControl(msMaxIter = 100),
                na.action = na.omit)
    }
    if (structure == "AR(1)") {
      RS <- gls(model = formula, data = data,
                correlation = corAR1(form = ~ time | id),
                method = "ML", control = glsControl(msMaxIter = 100),
                na.action = na.omit)
    }
    return(RS$logLik[1] + (lmd - 1) * sum(log(y), na.rm = TRUE))
  }

  formula <- formula(formula)
  yc <- as.character(formula)[2]

  options(na.action = "na.pass")
  X <- model.matrix(formula, data = data)

  if (sum(is.na(X)) > 0L) {
    stop("There are missing values in explanatory variables.")
  }
  if ((deparse(substitute(time)) == "NULL" &
       deparse(substitute(id)) != "NULL") |
      (deparse(substitute(time)) != "NULL" &
       deparse(substitute(id)) == "NULL")) {
    stop("If either time or id is specified,
         both of time and id must be specified.")
  }
  data$y <- data[, yc]

  if (sum(data$y < 0, na.rm = TRUE) > 0L) {
    stop("outcome must be positive.")
  }
  if (deparse(substitute(time)) == "NULL" &
      deparse(substitute(id)) == "NULL") {
    try1 <- bcreg(formula, data)
    nt <- 1
    data$id <- as.character(seq_len(nrow(data)))
    msflg <- table(data$id, is.na(data$y))[, 1]
    N <- sum(msflg != 0)
  } else {
    evars <- as.character(formula)[3]
    timec <- deparse(substitute(time))
    data$id <- as.character(data[, deparse(substitute(id))])
    ttbl <- as.numeric(names(table(data[, timec])))
    nt <- length(ttbl)
    if (!structure %in% c("UN", "CS", "AR(1)")) {
      stop("Select structure from c(\"UN\", \"CS\", \"AR(1)\")")
    }
    time2 <- data[, timec]
    for (i in 1:nt) {
      time2[data[, timec] == ttbl[i]] <- i
    }
    data$time <- time2
    msflg <- table(data$id, is.na(data$y))[, 1]
    N <- sum(msflg != 0)
    omis <- names(msflg)[msflg == 0]
    if (length(omis) > 0L) {
      for (i in seq_len(length(omis))) {
        data <- data[data$id != omis[i], ]
      }
    }
    data0 <- data[!duplicated(data$id), ]
    data0$y <- NULL
    data01 <- data0
    data1 <- c()
    for (i in 1:nt) {
      data01$time <- i
      data1 <- rbind(data1, data01)
    }
    data2 <- data[, c("id", "time", "y")]
    data <- merge(data2, data1, by = c("id", "time"), all = TRUE)
    for (i in 1:nt) {
      data[data$time == i, timec] <- ttbl[i]
    }

    if (nt == 1) {
      try1 <- bcreg(formula, data, lmdint)
    } else {
      formula2 <- formula(paste("z~", evars))
      try1 <- optimize(likn, interval = lmdint, maximum = TRUE,
                       formula = formula(paste("z~", evars)), data = data,
                       structure = structure)
    }
  }

  if (nt == 1) {
    le <- try1$lambda
    beta <- try1$beta
    alp <- try1$sigma ^ 2
    names(alp) <- "Residual"
    V <- as.matrix(alp)
    data$z <- (data$y ^ le - 1) / le
    bcres <- try1$betainf
    structure <- "UN"
    lik <- try1$logLik
    RS <- try1$lmObject
  } else {
    le <- try1$maximum
    lik <- try1$objective

    data$z <- (data$y ^ le - 1) / le

    if (structure == "UN") {
      RS <- gls(model = formula2, data = data,
                correlation = corSymm(form = ~ time | id),
                weights = varIdent(form = ~ 1 | time), method = "ML",
                control = glsControl(msMaxIter = 100), na.action = na.omit)
      RS$call$model <- eval(formula2)
      V <- covcalcUN(RS, nt)
      alp <- c()
      nmalp <- c()
      for (i in 1:nt) {
        for (j in i:nt) {
          alp <- c(alp, V[i, j])
          nmalp <- c(nmalp, paste("UN(", i, ",", j, ")", sep = ""))
        }
        names(alp) <- nmalp
      }
    }
    if (structure == "CS") {
      RS <- gls(model = formula2, data = data,
                correlation = corCompSymm(form = ~ time | id),
                method = "ML", control = glsControl(msMaxIter = 100),
                na.action = na.omit)
      RS$call$model <- eval(formula2)
      V <- covcalcCSAR(RS, nt)
      alp <- c(V[1, 2], V[1, 1] - V[1, 2])
      names(alp) <- c("ID", "Residual")
    }
    if (structure == "AR(1)") {
      RS <- gls(model = formula2, data = data,
                correlation = corAR1(form = ~ time | id),
                method = "ML", control = glsControl(msMaxIter = 100),
                na.action = na.omit)
      RS$call$model <- eval(formula2)
      V <- covcalcCSAR(RS, nt)
      alp <- c(V[1, 1], V[1, 2] / V[1, 1])
      names(alp) <- c("Variance", "Correlation")
    }
    beta <- RS$coefficients
    bcres <- summary(RS)$tTable
    time.tbl <- sort(unique(data[, deparse(substitute(time))]))
    rownames(V) <- time.tbl
    colnames(V) <- time.tbl
  }
  ns <- length(alp)
  msflg <- table(data$id, is.na(data$y))[, 1]
  N <- sum(msflg != 0)
  omis <- names(msflg)[msflg == 0]
  idt <- cbind((1:N) %x% (numeric(nt) + 1), (numeric(N) + 1) %x% (1:nt))

  if (length(omis) > 0L) {
    for (i in seq_len(length(omis))) {
      data <- data[data$id != omis[i], ]
    }
  }
  Ncc <- sum(msflg == nt)
  flg.na <- !is.na(data$y)
  n.data <- sum(flg.na)
  X <- model.matrix(formula, data = data)
  y <- data$y[flg.na]
  z <- data$z[flg.na]
  n.na <- length(flg.na) - n.data
  adj.prm <- c(N, Ncc, n.data, n.na)
  nb <- length(beta)
  Xb <- list()
  idt <- cbind((1:N) %x% (numeric(nt) + 1), (numeric(N) + 1) %x% (1:nt))
  for (j in 1:nb) {
    dum <- as.data.frame(cbind(idt, X[, j]))
    names(dum) <- c("a", "b", "c")
    Xb[[j]] <- as.matrix(ftable(xtabs(c ~ ., dum)))
  }
  dimnames(X) <- NULL
  nti <- table(data$id[!is.na(data$y)])
  ntic <- cumsum(nti)
  ntic <- c(0, ntic)
  Xflg <- (1:N) %x% (numeric(nt) + 1)
  X <- X[flg.na, ]
  Xflg <- Xflg[flg.na]
  ly <- log(y)
  dz <- le ^ (-2) * (y ^ le * (le * ly - 1) + 1)
  ddz <- le ^ (-1) * y ^ le * ly ^ 2 - 2 * le ^ (-2) * y ^ le * ly +
    2 * le ^ (-3) * (y ^ le - 1)
  if (!is.matrix(X)) {
    X <- as.matrix(X)
  }
  re <- z - X %*% beta
  dum <- as.data.frame(cbind(idt[flg.na, ], re))
  names(dum) <- c("a", "b", "c")
  rt <- as.matrix(xtabs(c ~ ., dum))
  dum <- as.data.frame(cbind(idt[flg.na, ], dz))
  names(dum) <- c("a", "b", "c")
  dzt <- as.matrix(xtabs(c ~ ., dum))
  dum <- as.data.frame(cbind(idt[flg.na, ], ddz))
  names(dum) <- c("a", "b", "c")
  ddzt <- as.matrix(xtabs(c ~ ., dum))
  dum <- as.data.frame(cbind(idt[flg.na, ], ly))
  names(dum) <- c("a", "b", "c")
  lyt <- as.matrix(xtabs(c ~ ., dum))
  lysum <- rowSums(lyt)
  dp <- numeric(N)
  t2 <- 2 ^ nt - 1
  nsf <- c()
  for (j1 in 1:nt) {
    for (j2 in j1:nt) {
      nsf <- rbind(nsf, t(c(j1, j2)))
    }
  }
  colf <- list()
  for (j in 1:(t2)) {
    mv <- dec2bin(j - 1, nt)
    aprt <- as.matrix(apply(rt == 0, 1, "==", mv))
    if (t2 == 1) {
      aprt <- t(aprt)
    }
    dp[colSums(aprt) == nt] <- j
    colf[[j]] <- (1:nt)[mv == 0]
  }
  ndp <- numeric(t2)
  dpf <- list()
  for (j in 1:t2) {
    ndp[j] <- sum(dp == j)
    dpf[[j]] <- (dp == j)
  }
  Hl <- 0
  Jl <- 0
  Hb <- matrix(0, nb, nb)
  Hs <- matrix(0, ns, ns)
  Hlb <- matrix(0, 1, nb)
  Hls <- matrix(0, 1, ns)
  Hbs <- matrix(0, nb, ns)
  Jb <- matrix(0, nb, nb)
  Js <- matrix(0, ns, ns)
  Jlb <- matrix(0, 1, nb)
  Jls <- matrix(0, 1, ns)
  Jbs <- matrix(0, nb, ns)
  for  (j in 1:t2) {
    if (ndp[j] != 0) {
      rj <- rt[dp == j, colf[[j]], drop = FALSE]
      nj <- nrow(rj)
      lysj <- lysum[dp == j]
      S <- V[colf[[j]], colf[[j]], drop = FALSE]
      nSj <- nrow(S)
      cholS <- chol(S)
      mat0 <- matrix(0, nSj, nSj)
      dS <- list()
      dS2 <- list()
      if (structure == "UN") {
        for (j1 in 1:ns) {
          dS[[j1]] <- 1 * (S == alp[j1])
          dS2[[j1]] <- list()
          for (j2 in 1 : ns) {
            dS2[[j1]][[j2]] <- mat0
          }
        }
      }
      if (structure == "CS") {
        dS[[1]] <- diag(nSj)
        dS[[2]] <- matrix(1, nSj, nSj)
        for (j1 in 1:2) {
          dS2[[j1]] <- list()
          for (j2 in 1:2) {
            dS2[[j1]][[j2]] <- mat0
          }
        }
      }
      if (structure == "AR(1)") {
        dS[[1]] <- S / S[1, 1]
        dS[[2]] <- mat0
        for (j1 in 1:2) {
          dS2[[j1]] <- list()
          for (j2 in 1:2) {
            dS2[[j1]][[j2]] <- mat0
          }
        }
        for (j1 in 1:nSj) {
          for (j2 in 1:nSj) {
            rp <- abs(j1 - j2)
            dS[[2]][j1, j2] <- rp * alp[1] * alp[2] ^ (rp - 1)
            dS2[[1]][[2]][j1, j2] <- rp * alp[2] ^ (rp - 1)
            dS2[[2]][[1]][j1, j2] <- rp * alp[2] ^ (rp - 1)
            dS2[[2]][[2]][j1, j2] <- rp * (rp - 1) * alp[1] *
              alp[2] ^ (rp - 2)
          }
        }
      }
      iS <- solve(S)
      sSr <- backsolve(cholS, t(rj), transpose = TRUE)
      sxb <- list()
      sdz <- list()
      sddz <- list()
      xjb <- list()
      cs_xSr <- list()
      dzj <- dzt[dp == j, colf[[j]], drop = FALSE]
      ddzj <- ddzt[dp == j, colf[[j]], drop = FALSE]
      sdz <- backsolve(cholS, t(dzj), transpose = TRUE)
      sddz <- backsolve(cholS, t(ddzj), transpose = TRUE)
      cs_zSr <- colSums(sdz * sSr)
      for (jb in 1:nb) {
        xj <- Xb[[jb]][dp == j, colf[[j]], drop = FALSE]
        xjb[[jb]] <- xj
        sxb0 <- backsolve(cholS, t(xj), transpose = TRUE)
        sxb[[jb]] <- sxb0
        cs_xSr[[jb]] <- colSums(sxb0 * sSr)
      }
      qAr <- list()
      cs_rAr <- list()
      rAA <- list()
      AA <- list()
      tris <- numeric(ns)
      for (js in 1:ns) {
        j1 <- nsf[js, 1]
        j2 <- nsf[js, 2]
        if (sum(colf[[j]] == j1) == 1 & sum(colf[[j]] == j2) == 1) {
          AA0 <- -iS %*% dS[[js]] %*% iS
          qr0 <- qr(AA0)
          AA[[js]] <- AA0
          rAA[[js]] <- qr.R(qr0)
          qx0 <- t(qr.Q(qr0)) %*% t(rj)
          rx0 <- qr.R(qr0) %*% t(rj)
          tris[js] <- tr(iS %*% dS[[js]])
          qAr[[js]] <- qx0
          cs_rAr[[js]] <- colSums(qx0 * rx0)
        } else {
          qAr[[js]] <- NA
        }
      }
      Hl <- Hl - sum(sddz * sSr) - sum(sdz ^ 2)
      Jl <- Jl + sum(lysj ^ 2) + sum(cs_zSr ^ 2) - 2 * sum(lysj * cs_zSr)
      for (j1 in 1:nb) {
        for (j2 in j1:nb) {
          if (!is.na(sxb[[j1]][1]) & !is.na(sxb[[j2]][1])) {
            Hb[j1, j2] <- Hb[j1, j2] - sum(sxb[[j1]] * sxb[[j2]])
            Jb[j1, j2] <- Jb[j1, j2] + sum(cs_xSr[[j1]] * cs_xSr[[j2]])
          }
        }
      }
      for (j1 in 1:ns) {
        for (j2 in j1:ns) {
          if (!is.na(qAr[[j1]][1]) & !is.na(qAr[[j2]][1])) {
            AA2 <- iS %*% (2 * dS[[j1]] %*% iS %*% dS[[j2]] -
                             dS2[[j1]][[j2]]) %*% iS
            qr0 <- qr(AA2)
            qx0 <- t(qr.Q(qr0)) %*% t(rj)
            rx0 <- qr.R(qr0) %*% t(rj)
            Hs[j1, j2] <- Hs[j1, j2] - 0.5 * nj *
              tr(AA[[j1]] %*% dS[[j2]] + iS %*% dS2[[j1]][[j2]]) -
              0.5 * sum(qx0 * rx0)
            Js[j1, j2] <- Js[j1, j2] + 0.25 * nj * tris[j1] * tris[j2] +
              0.25 * sum(cs_rAr[[j1]] * cs_rAr[[j2]]) +
              0.25 * sum(tris[j1] * cs_rAr[[j2]]) +
              0.25 * sum(tris[j2] * cs_rAr[[j1]])
          }
        }
      }
      for (j2 in 1:nb) {
        Hlb[1, j2] <- Hlb[1, j2] + sum(sdz * sxb[[j2]])
        Jlb[1, j2] <- Jlb[1, j2] + sum(lysj * cs_xSr[[j2]]) -
          sum(cs_zSr * cs_xSr[[j2]])
      }
      for (j2 in 1:ns) {
        if (!is.na(qAr[[j2]][1])) {
          dzrA <- rAA[[j2]] %*% t(dzj)
          Hls[1, j2] <- Hls[1, j2] - sum(dzrA * qAr[[j2]])
          Jls[1, j2] <- Jls[1, j2] - 0.5 * sum(lysj * tris[j2]) -
            0.5 * sum(lysj * cs_rAr[[j2]]) +
            0.5 * sum(cs_zSr * tris[j2]) + 0.5 * sum(cs_zSr * cs_rAr[[j2]])
        }
      }
      for (j1 in 1:nb) {
        for (j2 in 1:ns) {
          if (!is.na(qAr[[j2]][1])) {
            xrA <- rAA[[j2]] %*% t(xjb[[j1]])
            Hbs[j1, j2] <- Hbs[j1, j2] + sum(xrA * qAr[[j2]])
            Jbs[j1, j2] <- Jbs[j1, j2] -
              0.5 * sum(tris[[j2]] * cs_xSr[[j1]]) -
              0.5 * sum(cs_xSr[[j1]] * cs_rAr[[j2]])
          }
        }
      }
    }
    for (j1 in 1:nb) {
      for (j2 in 1:j1) {
        Hb[j1, j2] <- Hb[j2, j1]
        Jb[j1, j2] <- Jb[j2, j1]
      }
    }
    for (j1 in 1:ns) {
      for (j2 in 1:j1) {
        Hs[j1, j2] <- Hs[j2, j1]
        Js[j1, j2] <- Js[j2, j1]
      }
    }
  }
  H <- rbind(cbind(Hl, Hlb, Hls), cbind(t(Hlb), Hb, Hbs),
             cbind(t(Hls), t(Hbs), Hs))
  J <- rbind(cbind(Jl, Jlb, Jls), cbind(t(Jlb), Jb, Jbs),
             cbind(t(Jls), t(Jbs), Js))
  iII <- ginv(-H)
  iIIr <- iII %*% J %*% iII
  Vname <- c("lambda", paste("beta:", names(beta)),
             paste("alpha:", names(alp)))
  colnames(iII) <- Vname
  row.names(iII) <- Vname
  colnames(iIIr) <- Vname
  row.names(iIIr) <- Vname
  options(backup_options)
  structure(class = "bcmarg",
            list(call = Call, formula = formula, lambda = le, beta = beta,
                 alpha = alp, V = V, betainf = bcres, Vtheta.mod = iII,
                 Vtheta.rob = iIIr, logLik = lik, adj.prm = adj.prm,
                 glsObject = RS))
}

#' @export
logLik.bcmarg <- function(object, REML = FALSE, ...) {
  if (REML) {
    stop("REML method can not be used in bcmmrm.")
  }
  structure(object$logLik,
            df = 1 + length(c(object$beta, object$alpha)),
            class = "logLik")
}

#' @export
coef.bcmarg <- function(object, ...) {
  coef(object$glsObject)
}


#' @export
print.bcmarg <- function(x, digits = 3, ...) {
    mCall <- x$call
    covstr <- mCall$structure
    if (is.null(covstr) & !is.null(mCall$id)) {
      covstr <- "UN"
    }
    cat("Box-Cox transformed mixed model analysis\n")
    cat("  Formula:", deparse(x$formula), "\n")
    cat("  Time:", deparse(mCall$time), "\n")
    cat("  ID:", deparse(mCall$id), "\n")
    cat("  Covariance structure:", deparse(covstr), "\n")
    cat("  Data:", deparse(mCall$data), "\n")
    cat("  Estimated transformation parameter: ",
        format(x$lambda, digit = digits), "\n")
    cat("  Log-likelihood: ", format(x$logLik), "\n", sep = "")
    cat("\nCoefficients:\n")
    print(x$beta, digit = digits, ...)
    invisible(x)
  }

#' @export
summary.bcmarg <- function(object, ...) {
  structure(class = "summary.bcmarg", object)
}

#' @export
print.summary.bcmarg <- function(x, digits = 3, ...) {
    mCall <- x$call
    covstr <- mCall$structure
    if (is.null(covstr) & !is.null(mCall$id)) {
      covstr <- "UN"
    }
    V <- x$V
    Sd <- sqrt(diag(V))
    R <- V / (Sd %*% t(Sd))
    betainf <- x$betainf
    betainf[, 4] <- round(betainf[, 4], digits = digits)

    cat("Box-Cox transformed mixed model analysis\n")
    cat("  Formula:", deparse(x$formula), "\n")
    cat("  Time:", deparse(mCall$time), "\n")
    cat("  ID:", deparse(mCall$id), "\n")
    cat("  Covariance structure:", deparse(covstr), "\n")
    cat("  Data:", deparse(mCall$data), "\n")
    cat("  Log-likelihood: ", format(x$logLik), "\n", sep = "")
    cat("  Estimated transformation parameter: ",
        format(x$lambda, digits = digits), "\n")
    cat("\nCoefficients on the transformed scale:\n")
    print(betainf, digits = digits, ...)
    cat("\nCovariance parameters on the transformed scale:\n")
    print(x$alp, digits = digits)
    cat("\nCorrelations on the transformed scale:\n")
    print(R, digits = digits, ...)
    invisible(x)
}

#' @export
vcov.bcmarg <- function(object, robust = TRUE, ...) {
  if (robust) {
    vcov <- object$Vtheta.rob
  } else {
    vcov <- object$Vtheta.mod
  }
  return(vcov)
}


#' @export
fitted.bcmarg <- function(object, transformed = FALSE, ...) {
  z <- fitted(object$glsObject)
  if  (transformed) {
    fitted <- z
  } else {
    lmd <- object$lambda
    if (lmd == 0) {
      fitted <- exp(z)
    } else {
      fitted <- (lmd * z + 1) ^ (1 / lmd)
    }
  }
  return(fitted)
}


#' @export
predict.bcmarg <- function(object, newdata = NULL, transformed = FALSE, ...) {
  if (missing(newdata)) {
    z <- predict(object$glsObject)
  } else {
    z <- predict(object$glsObject, newdata)
  }
  if  (transformed) {
    pred<- z
  } else {
    lmd <- object$lambda
    if (lmd == 0) {
      pred<- exp(z)
    } else {
      pred <- (lmd * z + 1) ^ (1 / lmd)
    }
  }
  return(pred)
}
