as.gcloglog <- function(model, phi = NULL, boundary.tol = 1e-5) {
  if (is.null(phi)) phi <- environment(stats::family(model)$linkfun)$phi

  if (inherits(model, "glm")) {
    corr <- glm_correction(model, phi, boundary.tol)
    model$phi <- phi
    model$phi.status <- corr$status
    model$vcov.corrected <- corr$V
    class(model) <- c("gcloglog", class(model))
    model
  } else {
    corr <- glmm_correction(model, phi, boundary.tol)
    methods::new("gcloglogMerMod", model, phi = phi, phi.status = corr$status,
                 vcov.corrected = corr$V)
  }
}

#' @importClassesFrom lme4 glmerMod
#' @noRd
methods::setClass("gcloglog", representation("VIRTUAL"))
methods::setClass("gcloglogMerMod", contains = c("gcloglog", "glmerMod"),
                  slots = c(phi = "numeric", phi.status = "character", vcov.corrected = "matrix"))

#' @importFrom Matrix forceSymmetric
#' @export
#' @noRd
vcov.gcloglog <- function(object, correct = TRUE, ...) {
  if (isS4(object)) {
    phi <- object@phi; status <- object@phi.status; V <- object@vcov.corrected
  } else {
    phi <- object$phi; status <- object$phi.status; V <- object$vcov.corrected
  }
  if (!correct) return(NextMethod())
  if (status != "corrected") {
    warning(status, " (phi = ", signif(phi, 3), "); returning the covariance conditional on phi.",
            call. = FALSE)
    return(NextMethod())
  }
  if (isS4(object)) {
    # same form as lme4's vcov (dpoMatrix with cached correlation), which its
    # summary print method relies on
    V <- methods::as(Matrix::forceSymmetric(Matrix::Matrix(V)), "dpoMatrix")
    V@factors$correlation <- methods::as(V, "corMatrix")
    return(V)
  }
  # as vcov.glm(complete = TRUE): NA rows and columns for aliased coefficients
  cf <- stats::coef(object)
  Vc <- matrix(NA_real_, length(cf), length(cf), dimnames = list(names(cf), names(cf)))
  Vc[!is.na(cf), !is.na(cf)] <- V
  Vc
}

#' @export
#' @noRd
summary.gcloglog <- function(object, correct = TRUE, ...) {
  if (isS4(object)) {
    phi <- object@phi; status <- object@phi.status; V <- object@vcov.corrected
    # lme4's summary calls vcov() internally, which would dispatch to ours
    s <- summary(methods::as(object, "glmerMod"), ...)
  } else {
    phi <- object$phi; status <- object$phi.status; V <- object$vcov.corrected
    s <- NextMethod()
  }
  s$phi.link <- phi
  s$phi.status <- if (correct) status else "not requested"
  if (s$phi.status == "corrected") {
    se <- sqrt(diag(V))
    s$coefficients[, 2] <- se
    s$coefficients[, 3] <- s$coefficients[, 1] / se
    s$coefficients[, 4] <- 2 * stats::pnorm(-abs(s$coefficients[, 3]))
    if (isS4(object)) s$vcov <- stats::vcov(object) else s$cov.unscaled <- s$cov.scaled <- V
  }
  class(s) <- c("summary.gcloglog", class(s))
  s
}

#' @export
#' @noRd
print.summary.gcloglog <- function(x, ...) {
  NextMethod()
  note <- switch(x$phi.status,
    corrected = "standard errors account for its estimation.",
    `not requested` = "standard errors are conditional on it.",
    paste0(x$phi.status, "; standard errors are conditional on it."))
  cat("\nLink shape phi = ", signif(x$phi.link, 4), " (estimated); ", note, "\n\n", sep = "")
  invisible(x)
}

#' @export
#' @noRd
confint.gcloglog <- function(object, parm, level = 0.95, method = c("Wald", "profile", "boot"), ...) {
  method <- match.arg(method)
  if (method != "Wald") {
    if (isS4(object)) return(stats::confint(methods::as(object, "glmerMod"), parm, level, method = method, ...))
    if (method == "boot") stop("method = \"boot\" is only available for glmer fits.")
    class(object) <- class(object)[-1]
    return(stats::confint(object, parm, level, ...))
  }
  est <- if (isS4(object)) lme4::fixef(object) else stats::coef(object)
  se <- sqrt(diag(as.matrix(stats::vcov(object))))
  a <- (1 - level) / 2
  ci <- cbind(est + stats::qnorm(a) * se, est + stats::qnorm(1 - a) * se)
  dimnames(ci) <- list(names(est), paste(format(100 * c(a, 1 - a), trim = TRUE,
                                                 scientific = FALSE, digits = 3), "%"))
  if (!missing(parm)) ci <- ci[parm, , drop = FALSE]
  ci
}
