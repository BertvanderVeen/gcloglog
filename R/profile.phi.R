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
#'model2 <- res$final.mod
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
      sum(phi * q * (log(phi / (phi + exp(eta))) + exp(eta) / (phi + exp(eta))) *(y / p - (N - y) / (1 - p)))
    }
  }

  if(missing(optimizer)){
  if("method" %in% names(optControl)){
  method = optControl$method
  optControl <- optControl[names(optControl) != "method"]
  }else{
    method  = "BFGS"
  }

  optr <- optim(log(10), fn = fn.glm, gr = gr, method = method, control = optControl, model = model, y = y, N = N, nll0 = nll0)
  }else{
  optr <- optimizer(log(10), fn.glm, gr, control = optControl, model = model, y = y, N = N, nll0 = nll0)
  }
  return(list(optr = optr))
}

#' @rdname profile.phi
#' @export
#' @method profile.phi merMod
profile.phi.merMod <- function(model, y, optimizer = bobyqa, optControl = list(maxit = 100), ...){

  nll0 = -logLik(model)

  # Here we do gradient free optimisatin
  optr <- optimizer(log(10), fn.merMod, control = optControl, model = model, y = y, nll0 = nll0)

  return(list(optr = optr))
}

#' @rdname profile.phi
#' @export
#' @method profile.phi default
profile.phi.default <- function(model, y, optimizer = bobyqa, optControl = list(maxit = 100), ...){
  # still needs to be adjusted for N>1 case
  # Here we do gradient free optimisatin
  # The "general" class of models is harder to implement analytical derivatives for..
  optr <- optimizer(log(10), fn.generic, control = optControl, model = model, y = y)

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
    9e9
  }else{
    nll <- -logLik(newmodelfn)
    if(nll < nll0)model <<- newmodelfn # improving starting values of glm
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
    9e9
  }else{
    nll <- -logLik(newmodelfn)
    if(nll < nll0)model <<- newmodelfn # improving starting values of glm
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
    9e9
  }else{
    -logLik(newmodelfn)
  }
}

