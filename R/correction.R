at_boundary <- "phi-hat is at the boundary phi = 0, where the fit is equivalent to cloglog regression"
not_pd <- "the joint information of the parameters and phi is not positive definite"

glm_correction <- function(model, phi, boundary.tol) {
  k <- sum(!is.na(stats::coef(model)))
  V_na <- matrix(NA_real_, k, k)
  N <- model$prior.weights
  ll <- function(p) sum(stats::dbinom(round(N * model$y), round(N), p, log = TRUE))
  if (ll(stats::fitted(model)) - ll(stats::binomial("cloglog")$linkinv(model$linear.predictors)) < boundary.tol)
    return(list(status = at_boundary, V = V_na))
  info <- glm_joint_info(model, phi)
  if (!is.finite(info$den) || info$den <= 0) return(list(status = not_pd, V = V_na))
  list(status = "corrected", V = info$Sigma + outer(info$Sb_I, info$Sb_I) / info$den)
}

# Joint information of (beta, log phi) for a binomial gcloglog GLM
glm_joint_info <- function(object, phi) {
  cf <- stats::coef(object)
  X <- stats::model.matrix(object)[, !is.na(cf), drop = FALSE]
  Sigma <- stats::summary.glm(object)$cov.scaled
  N <- object$prior.weights
  y <- object$y * N
  eta <- object$linear.predictors

  e <- exp(eta)
  u <- 1 + phi * e
  log1mp <- -log1p(phi * e) / phi          # log(1 - p)
  p <- -expm1(log1mp)
  A <- e / u + log1mp                      # dp/dlog(phi) = (1 - p) * A
  r <- y / p - N
  dp_deta <- (1 - p) * e / u

  cross <- -y / p^2 * dp_deta * A + r * (e / u^2 - e / u)
  curv <- -y / p^2 * (1 - p) * A^2 + r * (-phi * (e / u)^2 - log1mp - e / u)

  I_bl <- -colSums(cross * X)              # I_{beta, log phi}
  I_ll <- -sum(curv) - min(-sum(r * A), 0) # I_{log phi, log phi}. min() applies a correction: below phi, log phi loses curvature much faster than phi, so the hessian returned is in terms of phi.
  Sb_I <- as.vector(Sigma %*% I_bl)
  list(Sigma = Sigma, Sb_I = Sb_I, den = I_ll - sum(I_bl * Sb_I))
}

glmm_correction <- function(model, phi, boundary.tol) {
  beta <- lme4::fixef(model)
  k <- length(beta)
  V_na <- matrix(NA_real_, k, k)
  nAGQ <- model@devcomp$dims[["nAGQ"]]
  theta <- lme4::getME(model, "theta")
  lower <- lme4::getME(model, "lower")
  free <- !is.finite(lower) | theta - lower > 1e-8
  nt <- sum(free)
  pars <- function(p) c(replace(theta, free, p[seq_len(nt)]), p[nt + seq_len(k)])
  devfuns <- list()
  nll <- function(p) {
    key <- sprintf("%.17g", p[nt + k + 1])
    if (is.null(devfuns[[key]]))
      devfuns[[key]] <<- glmm_devfun(model, stats::binomial(link = make.gcloglog(exp(p[nt + k + 1]))), nAGQ)
    devfuns[[key]](pars(p)) / 2
  }
  par0 <- c(theta[free], beta, log(phi))

  # phi-hat is at the boundary phi = 0 if the fit does not improve on the cloglog link,
  # the phi -> 0 limit, at the same parameters by more than boundary.tol; if PIRLS for the
  # cloglog link does not converge at these parameters, they are far from the cloglog fit
  nll0 <- tryCatch(glmm_devfun(model, stats::binomial("cloglog"), nAGQ)(pars(par0)) / 2, error = function(e) Inf)
  if (nll0 - nll(par0) < boundary.tol)
    return(list(status = at_boundary, V = V_na))

  # joint Hessian in (theta, beta, log phi). min() applies a correction: below phi, log phi
  # loses curvature much faster than phi, so the hessian returned is in terms of phi.
  H <- numDeriv::hessian(nll, par0)
  slope <- numDeriv::grad(function(lp) nll(replace(par0, nt + k + 1, lp)), log(phi))
  H[nt + k + 1, nt + k + 1] <- H[nt + k + 1, nt + k + 1] - min(slope, 0)
  incl <- nt + seq_len(k)
  incld <- c(seq_len(nt), nt + k + 1)
  A <- H[incl, incl, drop = FALSE]
  B <- H[incl, incld, drop = FALSE]
  D <- H[incld, incld, drop = FALSE]
  I <- A - B %*% solve(D, t(B))
  I <- 0.5 * (I + t(I))
  if (min(eigen(D, symmetric = TRUE, only.values = TRUE)$values) <= 0)
    return(list(status = not_pd, V = V_na))
  V <- solve(I)
  dimnames(V) <- list(names(beta), names(beta))
  list(status = "corrected", V = V)
}

glmm_devfun <- function(model, family, nAGQ) {
  reTrms <- lme4::getME(model, c("Zt", "theta", "Lambdat", "Lind", "Gp", "lower", "flist", "cnms", "Ztlist"))
  devfun <- lme4::mkGlmerDevfun(model@frame, lme4::getME(model, "X"), reTrms, family, nAGQ = 0L,
                                control = lme4::glmerControl(tolPwrss = 1e-12))
  lme4::updateGlmerDevfun(devfun, reTrms, nAGQ = nAGQ)
}
