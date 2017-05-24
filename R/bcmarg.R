#' Marginal Model of the Mixed Effect Model with the Box-Cox Transformation.
#'
#' \code{bcmarg} returns the inference resutls the parameters of
#'   the marginal model of the linear mixed effect model with the Box-Cox
#'   transformation proposed by Maruo et al. (2017). If time and id are not
#'   specified, inference results reduce to the results for the context of
#'   linear regression model proveded by Maruo et al. (2015).
#'
#' @param formula a two-sided linear formula object describing the model, with
#'   the response on the left of a ~ operator and the terms, separated by +
#'   operators, on the right.
#' @param data a data frame in which to interpret the variables named in the
#'   formula.
#' @param time time variable name for repeated measurements. Default is NULL.
#' @param id subject id variable name for repeated measurements. Default is
#'   NULL.
#' @param structure specify the covariance structure from c("UN", "CS",
#'   "AR(1)"). Default is "UN".
#'
#' @return bcmarg returns a list including following conponents for the
#' marginal model parameter inference:
#' \describe{
#'   \item{\code{lambda}}{estimate of a transformation parameter.}
#'   \item{\code{beta}}{estimate of a regression parameter vector.}
#'   \item{\code{alpha}}{estimate of a scale parameter vector.}
#'   \item{\code{V}}{variace-covariance matrix for any subject with no missing
#'         values.}
#'   \item{\code{transformed}}{inference results for beta under the assumption
#'         that lambda is known.}
#'   \item{\code{Vtheta.mod}}{model-based variance-covariance matrix for MLE of
#'         a vector of all parameters.}
#'   \item{\code{Vtheta.rob}}{robust variance-covariance matrix for MLE of
#'         a vector of all parameters.}
#'   \item{\code{lik}}{maximized liklihood.}
#'   \item{\code{adj.prm}}{small-sample adjustment parameter used in
#'         \code{\link{bcmmrm}}: c(number of subjects, number of
#'         completed subjects, number of outcome observations,
#'         number of missing observations)}
#' }
#'
#' @references \itemize{
#'   \item Maruo, K., Isogawa, N., Gosho, M. (2015). Inference of median
#'         difference based on the Box-Cox model in randomized clinical trials.
#'         \emph{Statistics in Medicine}, 34, 1634-1644.\cr
#'   \item Maruo, K., Yamaguchi, Y., Noma, H., Gosho, M. (2017). Interpretable
#'         inference on the mixed effect model with the Box-Cox transformation.
#'         \emph{Statistics in Medicine}, early view (DOI: 10.1002/sim.7279).
#' }
#'
#' @seealso \code{\link{bcmmrm}}
#'
#' @examples
#'  data(aidscd4)
#'  bcmarg(cd4 ~ as.factor(treatment) + as.factor(weekc) + as.factor(treatment):as.factor(weekc), aidscd4, weekc, id, "AR(1)")
#'
#' @importFrom nlme gls
#' @importFrom MASS ginv
#'
#' @export
bcmarg <- function(formula, data, time = NULL, id = NULL, structure = "UN"){
  tr <- function(X) sum(diag(X))
  covcalcUN <- function(glsob, nt){
    flg <- 0
    j <- 0
    while (flg == 0){
      j <- j + 1
      corm <- corMatrix(glsob$modelStruct$corStruct)[[j]]
      if (nrow(corm) == nt) flg <- 1
    }
    varstruct <- glsob$modelStruct$varStruct
    varests <- coef(varstruct, uncons = F, allCoef = T)
    covm <- corm * glsob$sigma ^ 2 * t(t(varests)) %*% t(varests)
    return(covm)
  }
  covcalcCSAR <- function(glsob, nt){
    flg <- 0
    j <- 0
    while (flg == 0){
      j <- j + 1
      corm <- corMatrix(glsob$modelStruct$corStruct)[[j]]
      if (nrow(corm) == nt) flg <- 1
    }
    covm <- corm * glsob$sigma ^ 2
    return(covm)
  }
  dec2bin <- function(num, digit=0){
    if (num <= 0 && digit <= 0) {
      return(NULL)
    }else{
      return(append(Recall(num %/% 2, digit - 1), num %% 2))
    }
  }
  likn <- function(formula, data, structure, lmd){
    y <- data$y
    if (lmd == 0) data$z <- log(y)
    else data$z <- (y ^ lmd - 1) / lmd
    if (structure == "UN"){
      RS <- gls(model = formula, data = data,
                correlation = corSymm(form = ~ time | id),
                weights = varIdent(form = ~ 1 | time), method = "ML",
                control = glsControl(msMaxIter = 100), na.action = na.omit)
    }
    if (structure == "CS"){
      RS <- gls(model = formula, data = data,
                correlation = corCompSymm(form = ~ time | id),
                method = "ML", control = glsControl(msMaxIter = 100),
                na.action = na.omit)
    }
    if (structure == "AR(1)"){
      RS <- gls(model = formula, data = data,
                correlation = corAR1(form = ~ time | id),
                method = "ML", control = glsControl(msMaxIter = 100),
                na.action = na.omit)
    }
    return(RS$logLik[1] + (lmd - 1) * sum(log(y), na.rm = T))
  }

  formula <- formula(formula)
  yc <- as.character(formula)[2]

  options(na.action = "na.pass")
  X <- model.matrix(formula, data = data)

  if (sum(is.na(X)) > 0L){
    stop("There are missing values in explanatory variables.")
  }
  if ((deparse(substitute(time)) == "NULL" &
       deparse(substitute(id)) != "NULL") |
      (deparse(substitute(time)) != "NULL" &
       deparse(substitute(id)) == "NULL")){
    stop("If either time or id is specified,
         both of time and id must be specified.")
  }
  data$y <- data[, yc]

  if (sum(data$y < 0, na.rm = T) > 0L) {
    stop("outcome must be positive.")
  }
  if (deparse(substitute(time)) == "NULL" & deparse(substitute(id)) == "NULL"){
    try1 <- try(bcreg(formula, data))
    nt <- 1
    data$id <- as.character(1:nrow(data))
  } else {
    evars <- as.character(formula)[3]
    timec <- deparse(substitute(time))
    data$id <- as.character(data[, deparse(substitute(id))])
    ttbl <- as.numeric(names(table(data[, timec])))
    nt <- length(ttbl)
    time2 <- data[, timec]
    for (i in 1:nt){
      time2[data[, timec] == ttbl[i]] <- i
    }
    data$time <- time2
    nt <- length(unique(data$time))
    if (nt == 1){
      try1 <- try(bcreg(formula, data))
    } else {
      formula2 <- formula(paste("z~", evars))
      try1 <- try(optimize(likn, lower = -2, upper = 2, maximum = T,
                           formula = formula2, data = data,
                           structure = structure))
      if (class(try1) != "try-error"){
        if ((abs(try1$maximum) - 2) < 0.01){
          try1 <- try(optimize(likn, lower = -3, upper = 3, maximum = T,
                               formula = formula2, data = data,
                               structure = structure))
        }
      }
    }
  }
  if (class(try1) != "try-error"){
    if (nt == 1){
      le <- try1$lambda
      beta <- try1$beta
      alp <- try1$sigma ^ 2
      V <- as.matrix(alp)
      data$z <- (data$y ^ le - 1) / le
      bcres <- try1$transformed
      structure <- "UN"
      lik <- try1$lik
    } else {
      le <- try1$maximum
      lik <- try1$objective

      data$z <- (data$y ^ le - 1) / le

      if (structure == "UN"){
        RS <- gls(model = formula2, data = data,
                  correlation = corSymm(form = ~ time | id),
                  weights = varIdent(form = ~ 1 | time), method = "ML",
                  control = glsControl(msMaxIter = 100), na.action = na.omit)
        V <- covcalcUN(RS, nt)
        alp <- c()
        for (i in 1:nt){
          for (j in i:nt){
            alp <- c(alp, V[i, j])
          }
        }
      }
      if (structure == "CS"){
        RS <- gls(model = formula2, data = data,
                  correlation = corCompSymm(form = ~ time | id),
                  method = "ML", control = glsControl(msMaxIter = 100),
                  na.action = na.omit)
        V <- covcalcCSAR(RS, nt)
        alp <- c(V[1, 1] - V[1, 2], V[1, 2])
      }
      if (structure == "AR(1)"){
        RS <- gls(model = formula2, data = data,
                  correlation = corAR1(form = ~ time | id),
                  method = "ML", control = glsControl(msMaxIter = 100),
                  na.action = na.omit)
        V <- covcalcCSAR(RS, nt)
        alp <- c(V[1, 1], V[1, 2] / V[1, 1])
      }
      beta <- as.numeric(RS$coefficients)
      bcres <- summary(RS)$tTable
    }

    ns <- length(alp)
    msflg <- table(data$id, is.na(data$y))[, 1]
    N <- sum(msflg != 0)
    omis <- names(msflg)[msflg == 0]
    if (length(omis) > 0L){
      for (i in 1:length(omis)){
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
    for (j in 1:nb){
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
    for (j1 in 1:nt){
      for (j2 in j1:nt){
        nsf <- rbind(nsf, t(c(j1, j2)))
      }
    }
    colf <- list()
    for (j in 1:(t2)){
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
    for  (j in 1:t2){
      if (ndp[j] != 0){
        rj <- rt[dp == j, colf[[j]], drop = F]
        nj <- nrow(rj)
        lysj <- lysum[dp == j]
        S <- V[colf[[j]], colf[[j]], drop = F]
        nSj <- nrow(S)
        cholS <- chol(S)
        mat0 <- matrix(0, nSj, nSj)
        dS <- list()
        dS2 <- list()
        if (structure == "UN"){
          for (j1 in 1:ns){
            dS[[j1]] <- 1 * (S == alp[j1])
            dS2[[j1]] <- list()
            for (j2 in 1 : ns){
              dS2[[j1]][[j2]] <- mat0
            }
          }
        }
        if (structure == "CS"){
          dS[[1]] <- diag(nSj)
          dS[[2]] <- matrix(1, nSj, nSj)
          for (j1 in 1:2){
            dS2[[j1]] <- list()
            for (j2 in 1:2){
              dS2[[j1]][[j2]] <- mat0
            }
          }
        }
        if (structure == "AR(1)"){
          dS[[1]] <- S / S[1, 1]
          dS[[2]] <- mat0
          for (j1 in 1:2){
            dS2[[j1]] <- list()
            for (j2 in 1:2){
              dS2[[j1]][[j2]] <- mat0
            }
          }
          for (j1 in 1:nSj){
            for (j2 in 1:nSj){
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
        sSr <- backsolve(cholS, t(rj), transpose = T)
        sxb <- list()
        sdz <- list()
        sddz <- list()
        xjb <- list()
        cs_xSr <- list()
        dzj <- dzt[dp == j, colf[[j]], drop = F]
        ddzj <- ddzt[dp == j, colf[[j]], drop = F]
        sdz <- backsolve(cholS, t(dzj), transpose = T)
        sddz <- backsolve(cholS, t(ddzj), transpose = T)
        cs_zSr <- colSums(sdz * sSr)
        for (jb in 1:nb){
          xj <- Xb[[jb]][dp == j, colf[[j]], drop = F]
          xjb[[jb]] <- xj
          sxb0 <- backsolve(cholS, t(xj), transpose = T)
          sxb[[jb]] <- sxb0
          cs_xSr[[jb]] <- colSums(sxb0 * sSr)
        }
        qAr <- list()
        cs_rAr <- list()
        rAA <- list()
        AA <- list()
        tris <- numeric(ns)
        for (js in 1:ns){
          j1 <- nsf[js, 1]
          j2 <- nsf[js, 2]
          if (sum(colf[[j]] == j1) == 1 & sum(colf[[j]] == j2) == 1){
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
        for (j1 in 1:nb){
          for (j2 in j1:nb){
            if (!is.na(sxb[[j1]][1]) & !is.na(sxb[[j2]][1])){
              Hb[j1, j2] <- Hb[j1, j2] - sum(sxb[[j1]] * sxb[[j2]])
              Jb[j1, j2] <- Jb[j1, j2] + sum(cs_xSr[[j1]] * cs_xSr[[j2]])
            }
          }
        }
        for (j1 in 1:ns){
          for (j2 in j1:ns){
            if (!is.na(qAr[[j1]][1]) & !is.na(qAr[[j2]][1])){
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
        for (j2 in 1:nb){
          Hlb[1, j2] <- Hlb[1, j2] + sum(sdz * sxb[[j2]])
          Jlb[1, j2] <- Jlb[1, j2] + sum(lysj * cs_xSr[[j2]]) -
            sum(cs_zSr * cs_xSr[[j2]])
        }
        for (j2 in 1:ns){
          if (!is.na(qAr[[j2]][1])){
            dzrA <- rAA[[j2]] %*% t(dzj)
            Hls[1, j2] <- Hls[1, j2] - sum(dzrA * qAr[[j2]])
            Jls[1, j2] <- Jls[1, j2] - 0.5 * sum(lysj * tris[j2]) -
              0.5 * sum(lysj * cs_rAr[[j2]]) +
              0.5 * sum(cs_zSr * tris[j2]) + 0.5 * sum(cs_zSr * cs_rAr[[j2]])
          }
        }
        for (j1 in 1:nb){
          for (j2 in 1:ns){
            if (!is.na(qAr[[j2]][1])){
              xrA <- rAA[[j2]] %*% t(xjb[[j1]])
              Hbs[j1, j2] <- Hbs[j1, j2] + sum(xrA * qAr[[j2]])
              Jbs[j1, j2] <- Jbs[j1, j2] -
                0.5 * sum(tris[[j2]] * cs_xSr[[j1]]) -
                0.5 * sum(cs_xSr[[j1]] * cs_rAr[[j2]])
            }
          }
        }
      }
      for (j1 in 1:nb){
        for (j2 in 1:j1){
          Hb[j1, j2] <- Hb[j2, j1]
          Jb[j1, j2] <- Jb[j2, j1]
        }
      }
      for (j1 in 1:ns){
        for (j2 in 1:j1){
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
    res <- list(lambda = le, beta = beta, alpha = alp, V = V,
                transformed = bcres, Vtheta.mod = iII, Vtheta.rob = iIIr,
                lik = lik, adj.prm = adj.prm)
    return(res)
  }
}
