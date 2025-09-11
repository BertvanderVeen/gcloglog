#' @title MLE estimation of the shape of a binomial link function
#' @description This function takes a model object with a family argument, and estimates a parameter phi that controls the shape of the link function.
#'
#' @param model a model object which accepts a family argument (should correspond to "binomial") with a link function. The model should have corresponding \code{"\link{update}"} and \code{"\link{logLik}"} functionality implemented.
#' @param y response variable.
#' @param optimizer numerical optimisation routine used to estimate phi. Must accept the starting value, objective function, gradient, and control arugments (in that order). Defaults to \code{"\link{optim}"}.
#' @param optControl arguments passed to optimizer.
#' @param ... other arguments passed to \code{"\link{update}"}.
#'
#' @return The optimisation results of the profiling, with log(phi) as the estimated parameter.
#'
#' @author Bert van der Veen
#' @references
#' van der Veen and Hui (2025). In prep.
#'
#'@examples
#'# Define logit link via gcloglog
#'library(gcloglog)
#'gcloglog1 <- make.gcloglog(1)
#'
#'# Fit logistic regression with gcloglog
#'x <- rnorm(100)
#'b1 = 2.3
#'b0 = 0.1
#'eta = b0+b1*x
#'y = rbinom(length(x), 1, plogis(eta))
#'data <- data.frame(y = y,x = x)
#'model <- glm(y~x, family = binomial(link = gcloglog1), data = data)
#'
#'#' Find optimal link shape by profiling phi
#'res <- profile.phi(model)
#'exp(res$optr$par) # extract phi
#'
#' @importFrom nloptr bobyqa
#'
#' @rdname profile.phi
#' @export
#' @method profile.phi glm
profile.phi.glm <- function(model, y, optimizer = optim, optControl = list(method = "BFGS", maxit = 100), ...){
  N <- weights(model)

  nll0 <- -logLik(model)

  # objective in the GLM case
  gr <- function(logphi, model, y, N, nll0, ...){
    phi = exp(logphi)
    gcloglog <- make.gcloglog(phi)
    fit <- try(newmodelgr <- update(model, family = binomial(link = gcloglog), start = coef(model), ...), silent = TRUE)

    if(inherits(fit, "try-error")){
      # try without starting values
      fit <- try(newmodelgr <- update(model, family = binomial(link = gcloglog), ...), silent = TRUE)
    }

    if(inherits(fit, "try-error")){
      NA
    }else{
      eta = predict(newmodelgr)
      p = predict(newmodelgr, type="response")
      q = 1-p

      # gradient for nll
      a = eta-logphi
      -sum(phi * q * (-(pmax(0, a) + log1p(exp(-abs(a)))) + plogis(a)) *(y / p - (N - y) / (1 - p)))
    }
  }

  # could try to find a decent start
  # continue = TRUE
  # logphi.start = log(1.1)
  # while(continue){
  #   if(is.na(fn.glm(logphi.start, model = model, nll0 = nll0))){
  #   logphi.start = logphi.start+2
  #   }else{
  #     break
  #   }
  # }

  if(missing(optimizer)){
  if("method" %in% names(optControl)){
  method = optControl$method
  optControl <- optControl[names(optControl) != "method"]
  }else{
    method  = "BFGS"
  }

  optr <- try(optim(log(1.01), fn = fn.glm, gr = gr, method = method, control = optControl, model = model, y = y, N = N, nll0 = nll0), silent = TRUE)

  # Search a bit for a start
  if(inherits(optr,"try-error")){
    maxit = 20
    logphi.start = 1.1
    it <- 1
    while(inherits(optr,"try-error") && it  < maxit){
    optr <- try(optim(logphi.start, fn = fn.glm, gr = gr, method = method, control = optControl, model = model, y = y, N = N, nll0 = nll0), silent = TRUE)
    logphi.start <- logphi.start + 0.1
    it <- it + 1
    }
  }
  }else{
  optr <- try(optimizer(log(1.01), fn.glm, gr, control = optControl, model = model, y = y, N = N, nll0 = nll0), silent = TRUE)

  # Search a bit for a start
  if(inherits(optr,"try-error")){
    maxit = 20
    logphi.start = 1.1
    it <- 1
    while(inherits(optr,"try-error") && it  < maxit){
      optr <- try(optimizer(logphi.start, fn.glm, gr, control = optControl, model = model, y = y, N = N, nll0 = nll0), silent = TRUE)
      logphi.start <- logphi.start + 0.1
      it <- it + 1
    }
  }
  }
  return(list(optr = optr, start = coef(model)))
}

#' @rdname profile.phi
#' @export
#' @method profile.phi merMod
profile.phi.merMod <- function(model, y, optimizer = bobyqa, optControl = list(maxit = 100), ...){

  nll0 = -logLik(model)

  # Here we do gradient free optimisatin
  optr <- optimizer(log(1.01), fn.merMod, control = optControl, model = model, y = y, nll0 = nll0)

  return(list(optr = optr, start = list(fixef = fixef(model),
                                                     theta =  getME(model, "theta"))))
}

#' @rdname profile.phi
#' @export
#' @method profile.phi default
profile.phi.default <- function(model, y, optimizer = bobyqa, optControl = list(maxit = 100), ...){
  # still needs to be adjusted for N>1 case
  # Here we do gradient free optimisatin
  # The "general" class of models is harder to implement analytical derivatives for..
  optr <- optimizer(log(1.01), fn.generic, control = optControl, model = model, y = y)

  return(list(optr = optr))
}

#' @export
profile.phi <- function(x, ...) {
  UseMethod("profile.phi")
}

fn.glm <- function(logphi, model, nll0, ...){
  phi = exp(logphi)
  gcloglog <- make.gcloglog(phi)

  args <- list(...)
  args$y <- NULL
  args$N <- NULL

  args$object = model
  args$family = binomial(link = gcloglog)

  args$start = coef(model)

  fit <- try(newmodelfn <- do.call(update, args), silent = TRUE)

  if(inherits(fit, "try-error") && !is.null(args$start)){
    args$start <- NULL
    # try without starting values
    fit <- try(newmodelfn <- do.call(update, args), silent = TRUE)
  }

  if(inherits(fit, "try-error")){
    NA
  }else{
    nll <- -logLik(newmodelfn)
    if(nll < nll0){
      # improving starting values of glm
      assign("model", newmodelfn, envir = parent.frame())
      assign("nll0", nll, envir = parent.frame())
      }
    nll
  }
}

fn.merMod <- function(logphi, model, nll0, ...){
  phi = exp(logphi)
  gcloglog <- make.gcloglog(phi)

  args <- list(...)
  args$y <- NULL
  args$N <- NULL

  args$object = model
  args$family = binomial(link = gcloglog)

  args$start = list(fixef = fixef(model),
                    theta =  getME(model, "theta"))

  fit <- try(newmodelfn <- do.call(update, args), silent = TRUE)

  if(inherits(fit, "try-error") && !is.null(args$start)){
    args$start <- NULL
    # try without starting values
    fit <- try(newmodelfn <- do.call(update, args), silent = TRUE)
  }

  if(inherits(fit, "try-error")){
    NA
  }else{
    nll <- -logLik(newmodelfn)
    if(nll < nll0){
      # improving starting values of glm
      assign("model", newmodelfn, envir = parent.frame())
      assign("nll0", nll, envir = parent.frame())
      }
    nll
  }
}

fn.generic <- function(logphi, model, ...){
  phi = exp(logphi)
  gcloglog <- make.gcloglog(phi)

  args <- list(...)
  args$y <- NULL
  args$N <- NULL

  args$object = model
  args$family = binomial(link = gcloglog)

  fit <- try(newmodelfn <- do.call(update, args), silent = TRUE)

  if(inherits(fit, "try-error")){
    NA
  }else{
    -logLik(newmodelfn)
  }
}

#' @export profile.phi.CI
profile.phi.CI <- function(logphi.mle, model,
                           h = 0.02, ytol = 2, ystep = 0.1,
                           maxit = 100, adaptive = TRUE, trace = TRUE,...) {

  if(inherits(model, "glm")){
    fn <- function(logphi.mle, model, nll=0)fn.glm(logphi.mle, model, nll)
  }else  if(inherits(model, "merMod")){
    fn <- function(logphi.mle, model, nll=0)fn.merMod(logphi.mle, model, nll)
  }else{
    fn <- fn.default
  }

  nll.mle <- fn(logphi.mle, model = model,...)

  eval_along <- function(start, direction) {
    x <- start
    yvals <- nll.mle
    hcurrent <- h
    iter <- 0

    repeat {
      iter <- iter + 1
      if (iter > maxit) break

      xnext <- tail(x, 1) + direction * hcurrent
      nllnext <- tryCatch(fn(xnext, model=model,...), error=function(e) NA)

      # If NA, shrink step and try again
      if (is.na(nllnext)) {
        hcurrent <- hcurrent / 2
        if (hcurrent < 1e-6) break
        next
      }

      # Stop if drop exceeds tolerance
      if ((nllnext-nll.mle) > ytol) break

      # Append
      x <- c(x, xnext)
      yvals <- c(yvals, nllnext)

      # Adaptive step only if at least 2 differences
      if (adaptive && length(yvals) > 2) {
        delta <- abs(yvals[length(yvals)] - yvals[length(yvals)-1])
        if (delta > ystep) hcurrent <- hcurrent / 2
        if (delta < ystep/4) hcurrent <- hcurrent * 2
      }
    }

    data.frame(logphi = x, logLik = -yvals)
  }

  if (trace) cat("Profiling downwards\n")
  down <- eval_along(logphi.mle, -1)
  if (trace) cat("Profiling upwards\n")
  # Ensure that we step enough for a CI even if the likelihood is flat on the left
  ytol.down <- -nll.mle-min(down$logLik)
  if(ytol.down<ytol)ytol = ytol + 2-ytol.down

  up <- eval_along(logphi.mle, 1)

  out <- rbind(up, down)
  out <- out[order(out$logphi), ]
  rownames(out) <- NULL
  return(out)
}

#' @title MLE estimation of the shape of a binomial link function
#' @description This function takes a model object, estimates a parameter phi that controls the shape of the link function, and performs further profiling for retrieving a confidence interval. The implementation heavily borrows from \code{"\link{TMB::tmbprofile}"}.
#'
#' @param model a model object which accepts a family argument (should correspond to "binomial") with a link function. The model should have corresponding \code{"\link{update}"} and \code{"\link{logLik}"} functionality implemented.
#' @param CI logical, defaults to \code{TRUE}.
#' @param alpha the confidence level, defaults to 0.05.
#' @param plot logical, defaults to \code{TRUE}, so that a plot of the profiling is displayed. Only when CI is set to \code{TRUE}.
#' @param h initial adaptive stepsize.
#' @param ytol Adjusts the range of the likelihood values.
#' @param ystep Adjusts the reoslution of the likelihood profile.
#' @param maxit Maximum number of iterations for the adaptive algorithm.
#' @param adaptive logical, defaults to \code{TRUE}. Implements adaptive step size.
#' @param trace logical, defaults to \code{TRUE}. Prints progress.
#' @param ... other arguments passed to \code{"\link{profile.phi}"}.
#'
#' @return A list including the optimisation results, CI, and final model fit.
#'
#' @author Bert van der Veen
#' @references
#' van der Veen and Hui (2025). In prep.
#'
#'@examples
#'# Define logit link via gcloglog
#'library(gcloglog)
#'gcloglog1 <- make.gcloglog(1)
#'
#'# Profile the gcloglog link
#'x <- rnorm(100)
#'b1 = 2.3
#'b0 = 0.1
#'eta = b0+b1*x
#'y = rbinom(length(x), 1, plogis(eta))
#'data <- data.frame(y = y,x = x)
#'model <- glm(y~x, family = binomial(link = gcloglog1), data = data)
#'res <- profil.gcloglog(model)
#'
#'final.model <- res$final.model
#' @export profile.gcloglog
profile.gcloglog <- function(model, CI = TRUE, method = "profile", alpha = 0.05, plot = TRUE, h = 0.02, ytol = 2, ystep = 0.1,
                             maxit = 100, adaptive = TRUE, trace = TRUE, ...){

  if(!grepl("gcloglog", family(model)$link)){
    model <- update(model, family = binomial(link = make.gcloglog(1)))
  }

  res <- profile.phi(model = model, y = model.response(model.frame(model)), ...)
  logphi.mle  <- res$optr$par

  gcloglog = make.gcloglog(exp(logphi.mle))
  final.model <- try(update(model, family = binomial(link=gcloglog)), silent = TRUE)

  if(inherits(final.model, "try-error"))final.model <- try(update(model, family = binomial(link=gcloglog), start = res$start), silent = TRUE)

  logLik.mle <- logLik(final.model)
  prof <- NULL
  ans = NA

  if(CI){
    prof <- profile.phi.CI(logphi.mle, final.model, h = h, ytol = ytol, ystep = ystep, maxit = maxit, adaptive = adaptive, trace = trace)
    prof$phi <- exp(prof$logphi)
    prof <- prof[, -1]

    tmp <- try({
      li <- subset(prof, phi<exp(logphi.mle))
      ui <- subset(prof, phi>exp(logphi.mle))

      ans <- numeric(2)
      threshold <- logLik.mle-0.5 * qchisq(1-alpha, df = 1)

      sp <- spline(x = li[, "phi"], y = li[, "logLik"])
      ans[1] <- approx(sp$y, sp$x, xout = threshold)$y
      if(is.na(ans[1]))ans[1]<- 0
      sp <- spline(x = ui[, "phi"], y = ui[, "logLik"])
      ans[2] <- approx(sp$y, sp$x, xout = threshold)$y
    })
    if(inherits(tmp,"try-error")){ans <- NA;print(tmp)}

    if(plot){
      plot(y = prof$logLik, prof$phi, type = "l", ylab = "log-likelihood", xlab = expression(hat(phi)), xlim = range(prof$phi))

      abline(v = ans[1], lty = "dotted")
      abline(v = ans[2], lty = "dotted")
      abline(v = exp(logphi.mle), col = "blue", lty = "dashed");mtext(at = exp(logphi.mle), side = 1, padj = 1, expression(hat(phi[max])), col = "blue")
      abline(v = 1, col = "red", lty = "dashed");mtext(at = 1.1, padj = -1, "logit", col = "red")
      abline(v = 0, col = "red", lty = "dashed");mtext(at = 0, padj = -1, "cloglog", col = "red")

    }
  }


  return(list(phi.mle = exp(logphi.mle), final.model = final.model, CI = ans, prof = prof))
}


