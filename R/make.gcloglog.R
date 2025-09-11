#' @title Create a Generalized complementary log-log link-function for GLM-type family
#' @description This function is used with \code{"\link{family}"} functions in \code{"\link{glm}"} type modeling functions.
#'
#' @param phi value for the dispersion parameter the link is evaluated for.
#'
#' @return An object of class "link-glm, a list with (amongst others) the components:
#'\describe{
#'  \item{\emph{linkfun}: }{ Link function \code{function(mu)} for fixed phi}
#'  \item{\emph{linkinv}: }{ Inverse link function \code{function(eta)} for fixed phi}
#'  \item{\emph{variance}: }{ variance function}
#'  \item{\emph{mu.eta}: }{ Derivative \code{function(eta)}}
#'  \item{\emph{valideta}: }{ \code{function(eta)}TRUE}
#'  \item{\emph{validmu}: }{ \code{function(mu)} that ensures 0<p<1}
#'  \item{\emph{name}: }{ the name for the link \code{gcloglog(phi)}}
#' }
#' @author Bert van der Veen
#' @references
#' van der Veen and Hui (2025). In prep.
#'
#' @examples
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
#' @export
make.gcloglog <- function(phi) {
  phi <- max(phi, .Machine$double.eps)
  link <- list(
    linkfun = function(mu) { -log(phi) - (phi) * log1p(-mu) + log1p(-exp(phi * log1p(-mu))) },
    linkinv = function(eta) { pmax(pmin(1 - (plogis(-log(phi)-eta))^(1/phi), 1 - .Machine$double.eps), .Machine$double.eps) },
    variance = binomial()$variance,
    dev.resids = binomial()$dev.resids,
    aic = binomial()$aic,
    mu.eta = function(eta) {
      exp(eta - (1/phi + 1) * (eta + log1p(exp(-eta-log(phi))))-(1/phi+1)*log(phi))
    },
    validmu = function(mu)all(is.finite(mu)) && all(mu > 0 & mu < 1),
    valideta = function(eta) TRUE,
    simulate = binomial()$simulate,
    initialize = binomial()$initialize,
    dispersion = 1,
    name = paste0("gcloglog(", phi, ")")
  )
  structure(link, class = "link-glm")
}
