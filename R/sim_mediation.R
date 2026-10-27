#' Simulated Mediation Study Data
#'
#' @description
#' A synthetic dataset mimicking a clustered education study with a continuous
#' treatment (tutoring hours), a continuous mediator (mid-year test score),
#' and a continuous outcome (end-of-year test score). Designed to illustrate
#' **RobustMediate** with realistic effect sizes and non-trivial confounding.
#'
#' @format A data frame with 600 rows (30 schools x 20 students) and 7 columns:
#' \describe{
#'   \item{school}{Factor. School identifier (30 levels). Use as `cluster_var`.}
#'   \item{Y}{Numeric. End-of-year test score (outcome).}
#'   \item{X}{Numeric. Tutoring hours received (continuous treatment, >= 0).}
#'   \item{M}{Numeric. Mid-year test score (mediator).}
#'   \item{Z1}{Numeric. Prior achievement (continuous covariate).}
#'   \item{Z2}{Integer (0/1). Free-lunch status (binary covariate).}
#'   \item{Z3}{Numeric. Parental education index (continuous covariate).}
#' }
#'
#' @section True parameter targets:
#' The data-generating process sets:
#' \itemize{
#'   \item **NDE** (X → Y direct path) ~= **0.25**
#'   \item **NIE** (X → M → Y path)    ~= **0.35**
#'   \item **TE**                       ~= **0.60**
#'   \item **% mediated**               ~= **58%**
#' }
#' Use these as a ground truth to assess estimation accuracy.
#'
#' @source Generated via `data-raw/generate_sim_data.R`. See that script for
#'   the full data-generating process.
#'
#' @examples
#' data(sim_mediation)
#' str(sim_mediation)
#' summary(sim_mediation[, c("Y","X","M")])
#'
#' \donttest{
#' fit <- robustmediate(
#'   treatment_formula = X ~ Z1 + Z2 + Z3,
#'   mediator_formula  = M ~ X + Z1 + Z2 + Z3,
#'   outcome_formula   = Y ~ X + M + Z1 + Z2 + Z3,
#'   data        = sim_mediation,
#'   cluster_var = "school",
#'   R           = 10,
#'   verbose     = FALSE
#' )
#' diagnose(fit)
#' }
"sim_mediation"
