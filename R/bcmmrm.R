#' Model Median Inference for Longitudinal Data in Randomized Clinical Trials.
#'
#' \code{bcmmrm} provides inference results for the model median differences
#' between treatment groups proposed by Maruo et al. (2017), which focuses on
#' continuous and positive longitudinally observed outcomes and a situation
#' where the efficacy of some treatments is compared based on a randomized,
#' parallel group clinical trial. If time and id are not specified, inference
#' results reduce to the results for the context of linear regression model
#' proveded by Maruo et al. (2015).
#'
#' @param outcome positive outcome (dependent) variable.
#' @param group treatment group variable.
#' @param data a data frame which includes outcome, group, time, id, and
#'   specified covariate variables.
#' @param time time variable for repeated measurements.
#' @param id subject id variable for repeated measurements.
#' @param timepoint analysis time point.
#' @param covv covariate variables. default is NULL.
#' @param cfactor nominal variable indicator for covaraite variables. nominal
#'   variable: 1, continuous variable: 0. default is NULL.
#' @param structure specify the covariance structure from c("UN", "CS",
#'   "AR(1)"). default is "UN"
#'
#' @return bcmmrm returns a list including following components for the model
#'   median inference. \describe{
#'     \item{median}{inrefence results for the model medians for all groups.}
#'     \item{meddif.mod}{model-based inference results for the model median
#'           differences between all pairs of groups.}
#'     \item{meddif.rob}{robust inference results for the model median
#'           differences between all pairs of groups.}
#'     \item{meddif.mod.adj}{model-based inference results for the model
#'           median differences between all pairs of groups with the empirical
#'           adjustment.}
#'     \item{meddif.rob.adj}{robust inference results for the model median
#'           differences between all pairs of groups with the empirical
#'           adjustment.}
#'     \item{lambda}{estimate of a transformation parameter.}
#'     \item{R}{correlation matrix for any subject with no missing values.}
#'     \item{transformed}{inference results for beta under the assumption
#'           that lambda is known.}
#'     \item{inf.marg}{result of \code{\link{bcmarg}} function.}
#'   }
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
#'         \emph{Statistics in Medicine}, early view (DOI: 10.1002/sim.7279).
#' }
#'
#' @seealso \code{\link{bcmarg}}
#'
#' @examples
#' data(aidscd4)
#' # covariates: no covariate, covariance structure: CS structure
#' bcmmrm(cd4, treatment, aidscd4, weekc, id, 32, structure="CS")
#'
#' # covariate: Box-Cox transformed baseline,age (continuous variables),
#' # and sex (nominal variable), covariance structure: UN structure
#' lmd.bl <- bcreg(cd4.bl ~ 1, aidscd4[aidscd4$weekc==8, ])$lambda
#' aidscd4$cd4.bl.tr <- (aidscd4$cd4.bl^lmd.bl-1)/lmd.bl
#' bcmmrm(cd4, treatment, aidscd4, weekc, id, 32, c("cd4.bl.tr"), c(0), structure = "AR(1)")
#'
#' @importFrom stats model.matrix pnorm pt qnorm qt
#'
#' @export
bcmmrm <- function(outcome, group, data, time = NULL, id = NULL,
                   timepoint = NULL, covv = NULL, cfactor = NULL,
                   structure = "UN"){
  y <- deparse(substitute(outcome))
  group <- deparse(substitute(group))
  data$group <- data[, group]
  data$y <- data[, y]
  covform <- ""
  cct <- 0
  if (!is.null(covv)){
    cct <- numeric(length(covv)) + 1
    for (i in 1:length(covv)){
      if (cfactor[i] == 0) {
        covform <- paste(covform, "+", covv[i])
      } else {
        covform <- paste(covform, "+as.factor(", covv[i], ")")
        cct[i] <- length(unique(data[, covv[i]])) - 1
      }
    }
  }

  if (deparse(substitute(time)) != "NULL" & deparse(substitute(id)) != "NULL"){
    time <- deparse(substitute(time))
    id <- deparse(substitute(id))
    data$time0 <- data[, time]
    time.tbl <- sort(unique(data$time0))
    for (i in 1:length(time.tbl)){
      data[data$time0 == time.tbl[i], "time"] <- i
      if (time.tbl[i] == timepoint) {
        tp0 <- i
      }
    }
    data$id <- data[, id]
    formula <- formula(paste("y ~ as.factor(group) + as.factor(time) +
                             as.factor(group):as.factor(time)", covform))
    nt <- length(unique(data$time))
    rm(id)
    rm(time)
    try1 <- try(bcmarg(formula, data, time, id, structure))
  } else {
    formula <- formula(paste("y ~ as.factor(group)", covform))
    data$time <- 1
    nt <- 1
    tp0 <- 1
    try1 <- try(bcmarg(formula, data))
  }

  if (class(try1) != "try-error"){
    le <- try1$lambda
    lik <- try1$lik
    beta <- try1$beta
    alp <- try1$alpha
    V <- try1$V
    RR <- V / (sqrt(diag(V)) %*% t(sqrt(diag(V))))
    ng <- length(unique(data$group))
    nc <- sum(cct)
    nb <- length(beta)
    ns <- length(alp)
    iII <- try1$Vtheta.mod
    iIIr <- try1$Vtheta.rob
    options(na.action = "na.pass")
    X <- model.matrix(formula, data = data)
    dimnames(X) <- NULL
    xcm <- 0
    if (!is.null(covv)){
      if (nc == 1)  {
        xcm <- mean(X[data$time == unique(data$time)[1], ng + nt])
      } else {
        xcm <- c(colMeans(X[data$time == unique(data$time)[1],
                            (ng + nt):(ng + nt + nc - 1)]))
      }
    }
    xi <- numeric(ng)
    bt <- 0
    dbt <- numeric(nt - 1)
    if (tp0 != 1){
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
    for (i in 1:ng){
      bg <- 0
      dbg <- numeric(ng - 1)
      if ( i != 1) {
        bg <- beta[i]
        dbg[i - 1] <- 1
      }
      bgt <- 0
      dbgt <- numeric((ng - 1) * (nt - 1))
      if (i != 1 & tp0 != 1){
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

    cbn <- choose(ng, 2)
    comb <- matrix(0, cbn, 2)
    count <- 1
    for (ii in 1:ng){
      for (jj in 1:ng){
        if (ii < jj) {
          comb[count, 1] <- ii
          comb[count, 2] <- jj
          count <- count + 1
        }
      }
    }
    delta <- numeric(cbn)
    SEd <- numeric(cbn)
    SEdr <- numeric(cbn)
    adj.prm <- try1$adj.prm
    N <- adj.prm[1]
    Ncc <- adj.prm[2]
    n.data <- adj.prm[3]
    n.na <- adj.prm[4]
    for (i in 1:cbn){
      delta[i] <- xi[comb[i, 2]] - xi[comb[i, 1]]
      SEd[i] <- sqrt(c(t(Dt[, comb[i, 2]] - Dt[, comb[i, 1]]) %*% iII %*%
                         (Dt[, comb[i, 2]] - Dt[, comb[i, 1]])))
      SEdr[i] <- sqrt(c(t(Dt[, comb[i, 2]] - Dt[, comb[i, 1]]) %*% iIIr %*%
                          (Dt[, comb[i, 2]] - Dt[, comb[i, 1]])))
    }
    if (structure == "UN" & nt != 1){
      nu <- Ncc - nt
      SEd2 <- SEd * sqrt(Ncc / nu)
      SEdr2 <- SEdr * sqrt(Ncc / nu)
    } else if (nt == 1) {
      nu <- n.data - nb
      SEd2 <- SEd * sqrt(n.data / nu)
      SEdr2 <- SEdr * sqrt(n.data / nu)
    } else{
      nu <- (N - ng) * (nt - 1) - n.na
      SEd2 <- SEd * sqrt(n.data / (n.data - nb))
      SEdr2 <- SEdr * sqrt(n.data / (n.data - nb))
    }
    sig.level <- 0.05
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
    lower <- delta - SEd * tt
    upper <- delta + SEd * tt
    lowerr <- delta - SEdr * tt
    upperr <- delta + SEdr * tt
    lower2 <- delta - SEd2 * tt2
    upper2 <- delta + SEd2 * tt2
    lowerr2 <- delta - SEdr2 * tt2
    upperr2 <- delta + SEdr2 * tt2
    medres <- data.frame(group = 1:ng, median = xi, SE.naive = SExi,
                         SE.robust = SExir)
    meddif.mod <- data.frame(group1 = comb[, 2], group0 = comb[, 1],
                             delta = delta, SE = SEd, lower.CL = lower,
                             upper.CL = upper, t.value = tvalue,
                             p.value = pvalue)
    meddif.rob <- data.frame(group1 = comb[, 2], group0 = comb[, 1],
                             delta = delta, SE = SEdr, lower.CL = lowerr,
                             upper.CL = upperr, t.value = tvaluer,
                             p.value = pvaluer)
    meddif.mod.adj <- data.frame(group1 = comb[, 2], group0 = comb[, 1],
                                 delta = delta, SE = SEd2, lower.CL = lower2,
                                 upper.CL = upper2, t.value = tvalue2,
                                 p.value = pvalue2)
    meddif.rob.adj <- data.frame(group1 = comb[, 2], group0 = comb[, 1],
                                 delta = delta, SE = SEdr2, lower.CL = lowerr2,
                                 upper.CL = upperr2, t.value = tvaluer2,
                                 p.value = pvaluer2)
    result <- list(median = medres, meddif.naive = meddif.mod,
                   meddif.rob = meddif.rob, meddif.mod.adj = meddif.mod.adj,
                   meddif.rob.adj = meddif.rob.adj, lambda = le, R = RR,
                   lik = lik, transformed = try1$transformed, inf.marg = try1)
  } else {
    result <- "Not converged"
  }
  return(result)
}
