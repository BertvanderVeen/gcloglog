#' @title MLE estimation of the shape of a binomial link function
#' @description This function takes a model object with a family argument, and estimates a parameter phi that controls the shape of the link function.
#'
#' @param model a model object which accepts a family argument (should correspond to "binomial") with a link function. The model should have corresponding \code{"\link{update}"} and \code{"\link{logLik}"} functionality implemented.
#' @param opt.control arguments passed to \code{"\link{optim}"}.
#' @param ... other arguments passed to \code{"\link{update}"}.
#'
#' @return A list of length 2, including the optimisation results and the final model object.
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
#' @export
profile.phi <- function(model, opt.control = list(maxit = 100, method = "BFGS"), ...){
  newmodelfn <- model

  gr <- function(logphi, model){
    y = model$y
    phi = exp(logphi)
    gcloglog <- make.gcloglog(phi)
    fit <- try(newmodelgr <- update(model, family = binomial(link = gcloglog), ...), silent = TRUE)

    if(inherits(fit, "try-error")){
      NA
    }else{
      eta = predict(newmodelgr)
      p = predict(newmodelgr, type="response")
      q = 1-p

      # gradient for nll
      sum(phi * q * (log(phi / (phi + exp(eta))) + exp(eta) / (phi + exp(eta))) *(y / p - (1 - y) / (1 - p)))
    }
  }
  # needs to be adjusted for Ntrials >1
  fn <- function(logphi, model, ...){
    phi = exp(logphi)
    gcloglog <- make.gcloglog(phi)
    fit <- try(newmodelfn <<- update(model, family = binomial(link = gcloglog), ...), silent = TRUE)

    if(inherits(fit, "try-error")){
      NA
    }else{
      -logLik(newmodelfn)
    }
  }
  if("method" %in% names(opt.control)){
  method = opt.control$method
  opt.control <- opt.control[names(opt.control) != "method"]
  }else{
    method  = "BFGS"
  }

  optr <- optim(1e-12, fn, gr, method = method, control = opt.control, model = model)

  return(list(optr = optr, final.mod = newmodelfn))
}
