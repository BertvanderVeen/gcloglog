#' @title Create a Generalized complementary log-log link-function for GLM-type family
#' @description This function is used with \code{"\link{family}"} functions in \code{"\link{glm}"} type modeling functions.
#'
#' @param phi value for the dispersion parameter the link is evaluated for.
#'
#' @return An object of class "link-glm, a list with components
#'\describe{
#'  \item{\emph{linkfun}: }{ Link function \code{function(mu)} for fixed phi}
#'  \item{\emph{linkinv}: }{ Inverse link function \code{function(eta)} for fixed phi}
#'  \item{\emph{mu.eta}: }{ Derivative \code{function(eta)}}
#'  \item{\emph{valideta}: }{ \code{function(eta)}TRUE}
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
  link <- list(
    linkfun = function(mu) { log(phi * ((1 - mu)^(-1/phi) - 1)) },
    linkinv = function(eta) { 1 - (phi / (phi + exp(eta)))^phi },
    mu.eta = function(eta) {
      phi^(phi+1)*(exp(eta)/(phi+exp(eta))^(phi+1))
    },
    valideta = function(eta) TRUE,
    name = paste0("gcloglog(", phi, ")")
  )
  structure(link, class = "link-glm")
}
