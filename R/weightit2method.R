#Propensity score estimation with regression
weightit2ps <- function(formula, data, s.weights, estimand, subset, stabilize, ps, ...) {
  #exact.factor should be a factor vector of values where each uniue value gets its own weights
  A <- list(...)
  if (length(A$link) == 0) A$link <- "logit"
  if (length(A$family) == 0) A$family <- binomial(link = A$link)

  if (length(ps) > 0) {
    mf <- model.frame(formula, data[subset,])
    t <- model.response(mf)
  }
  else {
    fit <- glm(formula, data = data.frame(data, .s.weights = s.weights)[subset,],
               weights = .s.weights,
               family = A$family,
               control = list(),
               ...)
    ps <- fit$fitted.values
    t <- fit$y
  }

  #Computing weights
  if (toupper(estimand) == "ATE") {
    w <- t/ps + (1-t)/(1-ps)
  }
  else if (toupper(estimand) == "ATT") {
    w <- t + (1-t)*ps/(1-ps)
  }
  else if (toupper(estimand) == "ATC") {
    w <- (1-t) + t*(1-ps)/ps
  }
  else if (toupper(estimand) == "ATO") {
    w <- t*(1-ps) + (1-t)*ps
  }
  else w <- NULL
  if (stabilize) {
    num.fit <- glm(t ~ 1, family = A$family)
    num.ps <- num.fit$fitted.values
    w <- w*ifelse(t == 1, mean(t), 1-mean(t))
  }

  obj <- list(ps = ps,
              w = w)
  return(obj)

}
weightit2ps.multi <- function(formula, data, s.weights, subset, estimand, focal, stabilize, ps, ...) {
  A <- list(...)
  if (length(A$link) == 0) A$link <- "logit"
  else {
    acceptable.links <- c("logit", "probit", "bayes.probit")
    which.link <- acceptable.links[pmatch(A$link, acceptable.links, nomatch = 0)][1]
    if (length(which.link) == 0) {
      A$link <- "logit"
      warning("Only  \"logit\",\"probit\" and \"bayes.probit\" are allowed as links for multinomial treatments. Using link = \"logit\".",
            call. = FALSE)
    }
    else A$link <- which.link
  }

  tt <- terms(formula)
  mf <- model.frame(tt, data[subset,], drop.unused.levels = TRUE)
  t <- model.response(mf)

  if (length(ps) > 0) {
  }
  else {
    if (A$link %in% c("logit", "probit")) {
      check.package("mlogit")
      covs <- model.matrix(tt, data=mf)[,-1]
      mult <- mlogit::mlogit.data(data.frame(.t = t, covs, .s.weights = s.weights[subset]), varying = NULL, shape = "wide", sep = "", choice = ".t")
      fit <- mlogit::mlogit(as.formula(paste0(".t ~ 1 | ", paste(colnames(covs), collapse = " + "),
                                        " | 1")), data = mult, estimate = TRUE,
                            probit = ifelse(A$link == "probit", TRUE, FALSE),
                            weights = .s.weights, ...)
      ps <- fitted(fit, outcome = FALSE)
    }
    else if (A$link == "bayes.probit") {
      check.package("MNP")
      fit <- MNP::mnp(formula, data[subset,], verbose = TRUE)
      ps <- MNP::predict.mnp(fit, type = "prob")$p
    }
  }
  #ps should be matrix of probs for each treat
  #Computing weights
  if (toupper(estimand) == "ATE") {
    w <- rep(0, nrow(ps))
    for (i in seq_len(nunique(t))) {
      w[t == levels(t)[i]] <- 1/ps[t == levels(t)[i], i]
    }
  }
  else if (toupper(estimand) == "ATT") {
    w <- rep(0, nrow(ps))
    for (i in seq_len(nunique(t))) {
      w[t == levels(t)[i]] <- 1/ps[t == levels(t)[i], i]
    }
    w <- w*ps[, which(levels(t) == focal)]
  }
  else if (toupper(estimand) == "ATO") {
    w <- rep(0, nrow(ps))
    for (i in seq_len(nunique(t))) {
      w[t == levels(t)[i]] <- 1/ps[t == levels(t)[i], i]
    }
    w <- w*apply(ps, 1, prod)
  }
  else w <- NULL

  #ps <- sapply(seq_along(t), function(i) ps[i, which(levels(t) == t[i])])

  obj <- list(w = w)
  return(obj)
}
weightit2ps.cont <- function(formula, data, s.weights, subset, stabilize, ps, ...) {
  A <- list(...)
  if (length(A$link) == 0) A$link <- "identity"
  if (length(A$family) == 0) A$family <- gaussian(link = A$link)

  mf <- model.frame(formula, data[subset,])
  t <- model.response(mf)

  stabilize <- TRUE

  if (length(ps) > 0) {
  }
  else {
    fit <- glm(formula, data = data.frame(data, .s.weights = s.weights)[subset,],
               weights = .s.weights,
               family = A$family,
               control = list(),
               ...)
    p.denom <- fit$fitted.values
    den.denom <- dnorm(t, p.denom, sqrt(summary(fit)$dispersion))

    if (stabilize) {
      if (length(A$num.formula) == 0) A$num.formula <- ~ 1
      num.fit <- glm(update.formula(A$num.formula, .t ~ .),
                     data = data.frame(.t =t, data, .s.weights = s.weights)[subset,],
                     weights = .s.weights,
                     family = A$family,
                     control = list(), ...)
      p.num <- num.fit$fitted.values
      den.num <- dnorm(t, p.num, sqrt(summary(num.fit)$dispersion))
      w <- den.num/den.denom
    }
    else {
      w <- 1/den.denom
    }
  }

  obj <- list(ps = p.denom,
              w = w)
  return(obj)
}

#Generalized boosted modeling with twang
weightit2gbm <- function(formula, data, s.weights, estimand, subset, stabilize, verbose, ...) {
  A <- list(...)
  if (length(A$stop.method) == 0) {
    warning("No stop.method was provided. Using \"es.mean\".",
            call. = FALSE, immediate. = TRUE)
    A$stop.method <- "es.mean"
  }
  else if (length(A$stop.method) > 1) {
    warning("Only one stop.method is allowed at a time. Using just the first stop.method (\"", A$stop.method[1],"\").",
            call. = FALSE, immediate. = TRUE)
    A$stop.method <- A$stop.method[1]
  }

  tt <- terms(formula)
  mf <- model.frame(tt, data[subset,], drop.unused.levels = TRUE)
  treat <- model.response(mf)
  covs <- model.matrix(tt, data=mf)[,-1]
  new.data <- data.frame(treat = treat, covs)

  if (estimand == "ATC") {
    treat <- 1 - treat
    estimand <- "ATT"
  }

  check.package("twang")
  fit <- do.call(twang::ps, c(list(formula = formula(new.data),
                                   data = new.data,
                                   estimand = estimand, sampw = s.weights[subset],
                                   verbose = verbose, print.level = 2), A))

  s <- names(fit$ps)[1]

  w <- cobalt::get.w(fit, stop.method = s)[[1]]

  obj <- list(w = w,
              ps = fit$ps[,s])

  return(obj)
}
weightit2gbm.multi <- function(formula, data, s.weights, estimand, focal, subset, stabilize, verbose, ...) {
  A <- list(...)
  if (length(A$stop.method) == 0) {
    warning("No stop.method was provided. Using \"es.mean\".",
            call. = FALSE, immediate. = TRUE)
    A$stop.method <- "es.mean"
  }
  else if (length(A$stop.method) > 1) {
    warning("Only one stop.method is allowed at a time. Using just the first stop.method.",
            call. = FALSE, immediate. = TRUE)
    A$stop.method <- A$stop.method[1]
  }

  tt <- terms(formula)
  mf <- model.frame(tt, data[subset,], drop.unused.levels = TRUE)
  treat <- model.response(mf)
  covs <- model.matrix(tt, data=mf)[,-1]
  new.data <- data.frame(treat = treat, covs)

  check.package("twang")
  fit <- do.call(twang::mnps, c(list(formula = formula(new.data),
                                     data = new.data,
                                     estimand = estimand, sampw = s.weights[subset],
                                     verbose = verbose, print.level = 2,
                                     treatATT = focal), A))

  s <- fit$stopMethods[1]

  w <- cobalt::get.w(fit, stop.method = s)[[1]]

  out <- list(w = w)
}

#CBPS
weightit2cbps <- function(formula, data, subset, estimand, verbose, s.weights, ...) {
  A <- list(...)

  tt <- terms(formula)
  mf <- model.frame(tt, data[subset,], drop.unused.levels = TRUE)
  treat <- model.response(mf)
  covs <- model.matrix(tt, data=mf)[,-1]
  new.data <- data.frame(treat = treat, covs)

  check.package("CBPS")
  fit <- CBPS::CBPS(formula(new.data), data = new.data, ATT = switch(estimand, ATT = 1, ATC = 2, ATE = 0),
                    method = ifelse(length(A$over) == 0 || isTRUE(A$over), "over", "exact"),
                    standardize = FALSE, sample.weights = s.weights[subset], ...)


  w <- cobalt::get.w(fit, estimand = switch(estimand, ATE = "ate", "att"))

  obj <- list(w = w,
              ps = fit$fitted.values)

  return(obj)
}
weightit2cbps.multi <- function(formula, data, subset, s.weights, verbose, ...) {
  A <- list(...)

  tt <- terms(formula)
  mf <- model.frame(tt, data[subset,], drop.unused.levels = TRUE)
  treat <- model.response(mf)
  covs <- model.matrix(tt, data=mf)[,-1]
  new.data <- data.frame(treat = treat, covs)

  check.package("CBPS")
  fit <- CBPS::CBPS(formula(new.data), data = new.data,
                    method = ifelse(length(A$over) == 0 || isTRUE(A$over), "over", "exact"),
                    standardize = FALSE, sample.weights = s.weights[subset], ...)


  w <- cobalt::get.w(fit)

  obj <- list(w = w)
}
weightit2cbps.cont <- weightit2cbps.multi
weightit2npcbps <- function(formula, data, subset, verbose, ...) {
  A <- list(...)

  tt <- terms(formula)
  mf <- model.frame(tt, data[subset,], drop.unused.levels = TRUE)
  treat <- factor(model.response(mf))
  covs <- model.matrix(tt, data=mf)[,-1]
  new.data <- data.frame(treat = treat, covs)

  check.package("CBPS")
  fit <- CBPS::npCBPS(formula(new.data), data = new.data, print.level = verbose, ...)

  w <- cobalt::get.w(fit)

  obj <- list(w = w)

  return(obj)
}
weightit2npcbps.multi <- function(formula, data, subset, estimand, verbose, ...) {
  A <- list(...)

  tt <- terms(formula)
  mf <- model.frame(tt, data[subset,], drop.unused.levels = TRUE)
  treat <- model.response(mf)
  covs <- model.matrix(tt, data=mf)[,-1]
  new.data <- data.frame(.t = factor(treat), covs)

  check.package("CBPS")
  fit <- CBPS::npCBPS(formula(new.data), data = new.data, print.level = verbose, ...)

  w <- cobalt::get.w(fit)

  obj <- list(w = w)

  return(obj)
}
weightit2npcbps.cont <- function(formula, data, subset, estimand, verbose, ...) {
  A <- list(...)

  tt <- terms(formula)
  mf <- model.frame(tt, data[subset,], drop.unused.levels = TRUE)
  treat <- model.response(mf)
  covs <- model.matrix(tt, data=mf)[,-1]
  new.data <- data.frame(treat = treat, covs)

  check.package("CBPS")
  fit <- CBPS::npCBPS(formula(new.data), data = new.data, print.level = verbose, ...)

  w <- cobalt::get.w(fit)

  obj <- list(w = w)

  return(obj)
}

#Entropy balancing with ebal
weightit2ebal <- function(formula, data, s.weights, subset, estimand, stabilize, verbose,...) {
  A <- list(...)

  tt <- terms(formula)
  mf <- model.frame(tt, data[subset,], drop.unused.levels = TRUE)
  treat <- model.response(mf)
  covs <- apply(model.matrix(tt, data=mf)[,-1], 2, make.smaller)

  #covs <- covs * replicate(ncol(covs), s.weights)

  check.package("ebal")
  if (estimand %in% c("ATT", "ATC")) {
    if (estimand == "ATC" ) treat_ <- 1 - treat
    else treat_ <- treat

    ebal.out <- ebal::ebalance(Treatment = treat_, X = covs,
                               print.level = ifelse(verbose, 3, -1),
                               base.weight = s.weights[subset][treat_ == 0], ...)
    if (stabilize) ebal.out <- ebal::ebalance.trim(ebal.out,
                                                   print.level = ifelse(verbose, 3, -1),
                                                   ...)

    w <- cobalt::get.w(ebal.out, treat = treat_)
  }
  else if (estimand == "ATE") {
    w <- rep(1, length(treat))

    for (i in unique(treat)) {
      #Reweight controls to be like total (need treated to look like total)
      covs_i <- rbind(covs, covs[treat==i,])
      treat_i <- c(rep(1, nrow(covs)), rep(0, sum(treat==i)))

      ebal.out_i <- ebal::ebalance(Treatment = treat_i, X = covs_i,
                                   base.weight = s.weights[subset][treat==i],
                                   print.level = 3, ...)
      if (stabilize) ebal.out_i <- ebal::ebalance.trim(ebal.out_i, ...)



      w[treat == i] <- ebal.out_i$w
    }
  }

  #w <- w/s.weights
  obj <- list(w = w)
  return(obj)

}
weightit2ebal.multi <- function(formula, data, s.weights, subset, estimand, focal, stabilize, verbose, ...) {
  A <- list(...)

  tt <- terms(formula)
  mf <- model.frame(tt, data[subset,], drop.unused.levels = TRUE)
  treat <- model.response(mf)
  covs <- apply(model.matrix(tt, data=mf)[,-1], 2, make.smaller)

  check.package("ebal")
  if (estimand %in% c("ATT")) {
    w <- rep(1, length(treat))
    control.levels <- levels(treat)[levels(treat) != focal]

    for (i in control.levels) {
      treat_ <- ifelse(treat[treat %in% c(focal, i)] == i, 0, 1)
      covs_ <- covs[treat %in% c(focal, i),]
      ebal.out <- ebal::ebalance(Treatment = treat_, X = covs_,
                                 base.weight = s.weights[subset][treat == i],
                                 print.level = ifelse(verbose, 3, -1), ...)
      if (stabilize) ebal.out <- ebal::ebalance.trim(ebal.out,
                                                     print.level = ifelse(verbose, 3, -1), ...)
      w[treat == i] <- ebal.out$w
    }
  }
  else if (estimand == "ATE") {
    w <- rep(1, length(treat))

    for (i in levels(treat)) {
      covs_i <- rbind(covs, covs[treat==i,])
      treat_i <- c(rep(1, nrow(covs)), rep(0, sum(treat==i)))

      ebal.out_i <- ebal::ebalance(Treatment = treat_i, X = covs_i,
                                   base.weight = s.weights[subset][treat == i],
                                   print.level = ifelse(verbose, 3, -1), ...)
      if (stabilize) ebal.out_i <- ebal::ebalance.trim(ebal.out_i,
                                                       print.level = ifelse(verbose, 3, -1),
                                                       ...)

      w[treat == i] <- ebal.out_i$w
    }
  }

  obj <- list(w = w)
  return(obj)
}

#Stable balancing weights with sbw
weightit2sbw <- function() {

}

#Empirical Balancing Calibration weights with ATE
weightit2ebcw <- function(formula, data, subset, estimand, ...) {
  A <- list(...)

  tt <- terms(formula)
  mf <- model.frame(tt, data[subset,], drop.unused.levels = TRUE)
  treat <- model.response(mf)
  covs <- model.matrix(tt, data=mf)[,-1]
  Y <- rep(0, length(treat))

  check.package("ATE")
  if (estimand %in% c("ATT", "ATC")) {
    if (estimand == "ATC" ) treat_ <- 1 - treat
    else treat_ <- treat

    ate.out <- ATE::ATE(Y = Y, Ti = treat_, X = covs,
                   ATT = TRUE, ...)
    w <- ate.out$weights.q
    w[treat_ == 1] <- 1

  }
  else if (estimand == "ATE") {
    ate.out <- ATE::ATE(Y = Y, Ti = treat, X = covs,
                   ATT = FALSE, ...)
    w <- ate.out$weights.q + ate.out$weights.p
  }

  obj <- list(w = w)
  return(obj)
}
weightit2ebcw.multi <- function(formula, data, subset, estimand, ...) {
  A <- list(...)

  tt <- terms(formula)
  mf <- model.frame(tt, data[subset,], drop.unused.levels = TRUE)
  treat <- as.numeric(model.response(mf)) - 1
  covs <- model.matrix(tt, data=mf)[,-1]
  Y <- rep(0, length(treat))

  check.package("ATE")
  ate.out <- ATE::ATE(Y = Y, Ti = treat, X = covs,
                 ATT = FALSE, ...)

  w <- apply(ate.out$weights.mat, 2, sum)

  obj <- list(w = w)
  return(obj)
}
