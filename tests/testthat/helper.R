library(testthat)
library(RobustMediate)

# Bivariate standard normal draws with correlation r (avoids a MASS dependency)
MASS_free_bivnorm <- function(n, r) {
  u <- stats::rnorm(n); v <- stats::rnorm(n)
  cbind(u, r * u + sqrt(1 - r^2) * v)
}
