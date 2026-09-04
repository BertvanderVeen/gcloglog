# gcloglog

Package that facilitates fitting a generalized complementary log-log model. Software is appendix to van der Veen and Hui (2026) in prep.

The link is indexed by a dispersion parameter phi, which arises from a latent negative-binomial count. It defines a family of link functions that encompasses the logit (phi = 1) and cloglog (phi -> 0) links jointly as special cases, with other values of phi giving further shapes.

## Installation

```r
# install.packages("remotes")
remotes::install_github("BertvanderVeen/gcloglog")
```

## Usage

### A link function at fixed phi

`make.gcloglog(phi)` returns a `"link-glm"` object that can be passed to `binomial()`, and so to `glm()`, `glmer()`, or any other function accepting a family:

```r
library(gcloglog)

set.seed(1)
n <- 500
x <- rnorm(n)
y <- rbinom(n, 1, plogis(0.1 + 2.3 * x))
data <- data.frame(y = y, x = x)

# phi = 1 reproduces logistic regression
model <- glm(y ~ x, family = binomial(link = make.gcloglog(1)), data = data)
coef(model)
```

### Estimating the shape of the link

`profile.gcloglog()` estimates phi by profile likelihood, returns a confidence interval, and refits the model at the estimate:

```r
res <- profile.gcloglog(model)

res$phi.mle          # 0.458
res$CI               # 0.000 1.560, includes 1, so logit is not rejected
coef(res$final.model)
```

By default a plot of the profile likelihood is drawn, with the logit (phi = 1) and cloglog (phi = 0) special cases marked. Use `plot = FALSE` to suppress it, and `CI = FALSE` to skip the (more expensive) profiling step.

If only the point estimate is needed, `profile.phi()` optimises log(phi) directly for a fitted model:

```r
res <- profile.phi(model, y = data$y)
exp(res$optr$par)
```

Methods are available for `glm` (analytical gradient), `merMod` (gradient-free, via `nloptr::bobyqa`), and a default method for any model class implementing `update()` and `logLik()`.

## References

Aranda-Ordaz, F. J. (1981). On two families of transformations to additivity for binary response data. *Biometrika*, 68(2), 357-363. doi:10.1093/biomet/68.2.357

van der Veen, B. and Hui, F. K. C. (2026). In prep.
