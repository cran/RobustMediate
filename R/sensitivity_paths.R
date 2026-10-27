#' Path-Wise and Indirect-Effect Robustness Values
#'
#' @description
#' Sensitivity analysis for the two paths of a linear mediation model and for
#' the indirect effect, in the partial \eqn{R^2} metric of Cinelli and Hazlett
#' (2020).
#'
#' One unmeasured confounder `U` is described by three partial \eqn{R^2}
#' values, all conditional on the observed covariates `C`:
#' \itemize{
#'   \item `r2_xu`: \eqn{R^2_{X \sim U \mid C}}, with the treatment;
#'   \item `r2_mu`: \eqn{R^2_{M \sim U \mid X, C}}, with the mediator;
#'   \item `r2_yu`: \eqn{R^2_{Y \sim U \mid M, X, C}}, with the outcome.
#' }
#' The a-path (treatment to mediator) depends on `r2_xu` and `r2_mu`. The
#' b-path (mediator to outcome) depends on `r2_mu` and `r2_yu`. Because
#' `r2_mu` enters both paths, one confounder related to the mediator weakens
#' both paths at once.
#'
#' A robustness value (RV) is the smallest \eqn{\rho} such that a confounder
#' with every partial \eqn{R^2} equal to \eqn{\rho} changes the conclusion.
#' The function reports:
#' \describe{
#'   \item{Path RVs}{`rv_q1` reduces the path estimate to zero; `rv_alpha`
#'     makes the path test non-significant (Cinelli and Hazlett, 2020).}
#'   \item{Indirect, joint test}{The indirect effect is declared significant
#'     when both paths are. Its RV equals the smaller of the two path RVs.}
#'   \item{Indirect, Sobel test}{The indirect effect is tested by
#'     \eqn{z = t_a t_b / \sqrt{t_a^2 + t_b^2}}. Its RV is smaller than the
#'     smaller path RV, because one confounder lowers both \eqn{t_a} and
#'     \eqn{t_b}.}
#' }
#' The ratio `overstatement` = RV(joint) / RV(Sobel) shows how much the
#' weaker-path value overstates robustness when the indirect effect is tested
#' as a product.
#'
#' @details
#' The adjusted t-statistic for a coefficient with t-statistic \eqn{t} and
#' residual degrees of freedom \eqn{df}, given treatment-side and outcome-side
#' partial \eqn{R^2} values \eqn{R^2_D} and \eqn{R^2_Y}, is
#' \deqn{t_{adj} = \frac{|t| - \sqrt{df}\sqrt{R^2_Y R^2_D / (1 - R^2_D)}}
#'   {\sqrt{(1 - R^2_Y)/(1 - R^2_D)}\sqrt{df/(df - 1)}}}
#' assuming the confounder biases the estimate toward zero. Path tests use the
#' t critical value with \eqn{df - 1} degrees of freedom; the Sobel test uses
#' the normal critical value.
#'
#' The result applies to linear (Gaussian, identity-link) mediator and outcome
#' models with the treatment and mediator entering as main effects.
#'
#' @param x A `robmedfit` object. Alternatively leave `x` as `NULL` and supply
#'   `mediator_model`, `outcome_model`, `treatment`, and `mediator`.
#' @param alpha Significance level. Default `0.05`.
#' @param mediator_model,outcome_model Fitted `lm` or Gaussian `glm` models
#'   for the mediator and the outcome. Used when `x` is `NULL`.
#' @param treatment,mediator Names of the treatment and mediator variables.
#'   Used when `x` is `NULL`.
#' @param rho_grid Values of \eqn{\rho} for the robustness curve stored in the
#'   result and used by [plot_paths()].
#'
#' @return An object of class `"robmed_paths"`: a list with
#' \describe{
#'   \item{`table`}{Data frame with one row per path and per indirect-effect
#'     test: `estimate`, `se`, `stat`, `df`, `rv_q1`, `rv_alpha`.}
#'   \item{`overstatement`}{RV(joint) / RV(Sobel).}
#'   \item{`curve`}{Data frame of adjusted `t_a`, `t_b`, and Sobel `z` over
#'     `rho_grid`.}
#'   \item{`inputs`}{t-statistics, degrees of freedom, and critical values.}
#' }
#'
#' @references
#' Cinelli, C., & Hazlett, C. (2020). Making sense of sensitivity: Extending
#' omitted variable bias. *Journal of the Royal Statistical Society: Series B*,
#' 82(1), 39--67. \doi{10.1111/rssb.12348}
#'
#' Sobel, M. E. (1982). Asymptotic confidence intervals for indirect effects
#' in structural equation models. *Sociological Methodology*, 13, 290--312.
#'
#' @seealso [adjust_paths()], [plot_paths()], [sensitivity_rho()]
#'
#' @examples
#' data(sim_mediation)
#' m_fit <- lm(M ~ X + Z1 + Z2 + Z3, data = sim_mediation)
#' y_fit <- lm(Y ~ X + M + Z1 + Z2 + Z3, data = sim_mediation)
#' sp <- sensitivity_paths(mediator_model = m_fit, outcome_model = y_fit,
#'                         treatment = "X", mediator = "M")
#' sp
#' plot_paths(sp)
#' @export
sensitivity_paths <- function(x = NULL, alpha = 0.05,
                              mediator_model = NULL, outcome_model = NULL,
                              treatment = NULL, mediator = NULL,
                              rho_grid = seq(0, 0.6, by = 0.005)) {
  if (!is.null(x)) {
    if (!inherits(x, "robmedfit")) rlang::abort("`x` must be a robmedfit object.")
    mediator_model <- x$models$mediator
    outcome_model  <- x$models$outcome
    treatment      <- x$meta$treat_var
    mediator       <- x$meta$med_var
  }
  if (is.null(mediator_model) || is.null(outcome_model) ||
      is.null(treatment) || is.null(mediator)) {
    rlang::abort(paste0(
      "Supply either a robmedfit object `x`, or all of `mediator_model`, ",
      "`outcome_model`, `treatment`, and `mediator`."
    ))
  }
  if (!is.numeric(alpha) || length(alpha) != 1L || alpha <= 0 || alpha >= 1)
    rlang::abort("`alpha` must be a single number in (0, 1).")
  if (any(rho_grid < 0 | rho_grid >= 1))
    rlang::abort("`rho_grid` must lie in [0, 1).")

  .check_linear_sem(mediator_model, outcome_model, treatment, mediator,
                    fn = "sensitivity_paths")

  pa <- .path_stats(mediator_model, treatment)
  pb <- .path_stats(outcome_model,  mediator)

  if (pa$df < 3 || pb$df < 3)
    rlang::abort("Too few residual degrees of freedom.")

  ta <- abs(pa$t); tb <- abs(pb$t)
  crit_a <- stats::qt(1 - alpha / 2, pa$df - 1)
  crit_b <- stats::qt(1 - alpha / 2, pb$df - 1)
  crit_z <- stats::qnorm(1 - alpha / 2)

  rv_a_q1 <- .rv_equal(ta, pa$df, 0)
  rv_b_q1 <- .rv_equal(tb, pb$df, 0)
  rv_a_al <- .rv_equal(ta, pa$df, crit_a)
  rv_b_al <- .rv_equal(tb, pb$df, crit_b)

  rv_joint <- min(rv_a_al, rv_b_al)

  z_obs <- .sobel(ta, tb)
  sobel_at <- function(rho) {
    .sobel(.tadj_equal(ta, pa$df, rho), .tadj_equal(tb, pb$df, rho))
  }
  upper <- min(rv_a_q1, rv_b_q1)
  rv_sobel <- if (sobel_at(0) <= crit_z) {
    0
  } else {
    stats::uniroot(function(r) sobel_at(r) - crit_z,
                   lower = 0, upper = upper, tol = 1e-12)$root
  }

  ab    <- pa$estimate * pb$estimate
  ab_se <- sqrt(pa$estimate^2 * pb$se^2 + pb$estimate^2 * pa$se^2)

  table <- data.frame(
    quantity = c("a-path", "b-path", "indirect (joint test)",
                 "indirect (Sobel test)"),
    estimate = c(pa$estimate, pb$estimate, ab, ab),
    se       = c(pa$se, pb$se, NA_real_, ab_se),
    stat     = c(pa$t, pb$t, NA_real_, sign(ab) * z_obs),
    df       = c(pa$df, pb$df, NA_real_, NA_real_),
    rv_q1    = c(rv_a_q1, rv_b_q1, min(rv_a_q1, rv_b_q1),
                 min(rv_a_q1, rv_b_q1)),
    rv_alpha = c(rv_a_al, rv_b_al, rv_joint, rv_sobel),
    stringsAsFactors = FALSE
  )

  curve <- data.frame(
    rho     = rho_grid,
    t_a     = .tadj_equal(ta, pa$df, rho_grid),
    t_b     = .tadj_equal(tb, pb$df, rho_grid),
    stringsAsFactors = FALSE
  )
  curve$z_sobel <- .sobel(curve$t_a, curve$t_b)

  structure(
    list(
      table         = table,
      overstatement = if (rv_sobel > 0) rv_joint / rv_sobel else NA_real_,
      curve         = curve,
      inputs        = list(t_a = pa$t, t_b = pb$t, df_a = pa$df, df_b = pb$df,
                           crit_a = crit_a, crit_b = crit_b, crit_z = crit_z,
                           alpha = alpha, treatment = treatment,
                           mediator = mediator)
    ),
    class = "robmed_paths"
  )
}


#' Adjust Both Paths for a Confounder of Given Strength
#'
#' @description
#' Returns the adjusted path t-statistics and the adjusted Sobel statistic for
#' a hypothetical confounder with the given partial \eqn{R^2} values. Use it
#' when the three values should differ, for example to explore a confounder
#' that is strongly related to the mediator but only weakly to the treatment.
#' See [sensitivity_paths()] for the definitions.
#'
#' @param x A `robmed_paths` object from [sensitivity_paths()].
#' @param r2_xu,r2_mu,r2_yu Partial \eqn{R^2} of the confounder with the
#'   treatment, the mediator, and the outcome. Vectors are recycled.
#'
#' @return A data frame with the inputs and columns `t_a`, `t_b`, `z_sobel`,
#'   `sig_joint`, and `sig_sobel`.
#'
#' @examples
#' data(sim_mediation)
#' m_fit <- lm(M ~ X + Z1 + Z2 + Z3, data = sim_mediation)
#' y_fit <- lm(Y ~ X + M + Z1 + Z2 + Z3, data = sim_mediation)
#' sp <- sensitivity_paths(mediator_model = m_fit, outcome_model = y_fit,
#'                         treatment = "X", mediator = "M")
#' adjust_paths(sp, r2_xu = 0.02, r2_mu = 0.10, r2_yu = 0.05)
#' @export
adjust_paths <- function(x, r2_xu, r2_mu, r2_yu) {
  if (!inherits(x, "robmed_paths"))
    rlang::abort("`x` must be a robmed_paths object from sensitivity_paths().")
  for (v in list(r2_xu, r2_mu, r2_yu)) {
    if (!is.numeric(v) || any(v < 0 | v >= 1))
      rlang::abort("Partial R2 values must lie in [0, 1).")
  }
  inp <- x$inputs
  out <- data.frame(r2_xu = r2_xu, r2_mu = r2_mu, r2_yu = r2_yu)
  out$t_a <- .tadj(abs(inp$t_a), inp$df_a, r2d = out$r2_xu, r2y = out$r2_mu)
  out$t_b <- .tadj(abs(inp$t_b), inp$df_b, r2d = out$r2_mu, r2y = out$r2_yu)
  out$z_sobel   <- .sobel(out$t_a, out$t_b)
  out$sig_joint <- out$t_a > inp$crit_a & out$t_b > inp$crit_b
  out$sig_sobel <- out$z_sobel > inp$crit_z
  out
}


#' Print a robmed_paths object
#'
#' @param x A `robmed_paths` object from [sensitivity_paths()].
#' @param digits Number of digits to print. Default `3`.
#' @param ... Ignored.
#' @return The input object, invisibly.
#' @export
print.robmed_paths <- function(x, digits = 3, ...) {
  cat("-- Path-wise and indirect-effect robustness values --------------\n")
  tab <- x$table
  num <- vapply(tab, is.numeric, logical(1L))
  tab[num] <- lapply(tab[num], function(v) round(v, digits))
  print(tab, row.names = FALSE)
  cat(sprintf(
    "\n  RV(joint) / RV(Sobel) = %.2f\n", x$overstatement
  ))
  cat("  rv_q1: confounder strength that reduces the estimate to zero.\n")
  cat(sprintf(
    "  rv_alpha: strength that makes the test non-significant (alpha = %.2f).\n",
    x$inputs$alpha
  ))
  cat("  If the indirect effect is tested as a product (Sobel, Monte Carlo,\n")
  cat("  bootstrap), report the Sobel value, not the smaller path value.\n")
  invisible(x)
}


#' Plot Adjusted Path and Sobel Statistics Against Confounder Strength
#'
#' @description
#' Plots the adjusted t-statistics of the two paths and the adjusted Sobel
#' statistic as the confounder strength \eqn{\rho} grows (all three partial
#' \eqn{R^2} values equal to \eqn{\rho}). Dashed vertical lines mark the
#' robustness values for the joint test and the Sobel test.
#'
#' @param x A `robmed_paths` object from [sensitivity_paths()].
#' @param ... Ignored.
#' @return A `ggplot2` object.
#'
#' @examples
#' data(sim_mediation)
#' m_fit <- lm(M ~ X + Z1 + Z2 + Z3, data = sim_mediation)
#' y_fit <- lm(Y ~ X + M + Z1 + Z2 + Z3, data = sim_mediation)
#' plot_paths(sensitivity_paths(mediator_model = m_fit, outcome_model = y_fit,
#'                              treatment = "X", mediator = "M"))
#' @importFrom ggplot2 ggplot aes geom_line geom_hline geom_vline labs
#'   theme_minimal theme scale_linetype_manual coord_cartesian
#' @export
plot_paths <- function(x, ...) {
  if (!inherits(x, "robmed_paths"))
    rlang::abort("`x` must be a robmed_paths object from sensitivity_paths().")
  cv <- x$curve
  long <- data.frame(
    rho  = rep(cv$rho, 3),
    stat = c(cv$t_a, cv$t_b, cv$z_sobel),
    line = factor(rep(c("a-path t", "b-path t", "Sobel z"), each = nrow(cv)),
                  levels = c("a-path t", "b-path t", "Sobel z")),
    stringsAsFactors = FALSE
  )
  rv <- x$table$rv_alpha
  ymax <- max(long$stat, na.rm = TRUE)

  ggplot2::ggplot(long, ggplot2::aes(x = rho, y = stat, linetype = line)) +
    ggplot2::geom_hline(yintercept = x$inputs$crit_z, colour = "grey50",
                        linewidth = 0.4) +
    ggplot2::geom_vline(xintercept = rv[3], linetype = "dashed",
                        colour = "grey40", linewidth = 0.4) +
    ggplot2::geom_vline(xintercept = rv[4], linetype = "dashed",
                        colour = "black", linewidth = 0.6) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::scale_linetype_manual(
      values = c("a-path t" = "dotted", "b-path t" = "dotdash",
                 "Sobel z" = "solid"),
      name = NULL
    ) +
    ggplot2::coord_cartesian(ylim = c(0, ymax)) +
    ggplot2::labs(
      x = expression("Confounder strength " * rho * " (each partial " * R^2 * ")"),
      y = "Adjusted test statistic",
      caption = paste0(
        "Horizontal line: critical value. Dashed lines: RV(joint) = ",
        round(rv[3], 3), " (grey), RV(Sobel) = ", round(rv[4], 3), " (black)."
      )
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(legend.position = "top")
}


# ── internal helpers ──────────────────────────────────────────────────────────

.path_stats <- function(model, var) {
  co <- stats::coef(model)
  se <- sqrt(diag(stats::vcov(model)))
  list(estimate = unname(co[var]), se = unname(se[var]),
       t = unname(co[var] / se[var]), df = stats::df.residual(model))
}

# General Cinelli-Hazlett adjusted t (bias toward zero), t >= 0
.tadj <- function(t, df, r2d, r2y) {
  bias <- sqrt(df) * sqrt(r2y * r2d / (1 - r2d))
  (t - bias) / (sqrt((1 - r2y) / (1 - r2d)) * sqrt(df / (df - 1)))
}

# Adjusted t when every partial R2 equals rho
.tadj_equal <- function(t, df, rho) {
  sqrt(df - 1) * (t / sqrt(df) - rho / sqrt(1 - rho))
}

# Smallest rho with .tadj_equal(t, df, rho) = crit; 0 if already below
.rv_equal <- function(t, df, crit) {
  c0 <- t / sqrt(df) - crit / sqrt(df - 1)
  if (c0 <= 0) return(0)
  (-c0^2 + c0 * sqrt(c0^2 + 4)) / 2
}

# Sobel statistic from two path t-statistics; 0 once a path reaches zero
.sobel <- function(ta, tb) {
  z <- ta * tb / sqrt(ta^2 + tb^2)
  z[ta <= 0 | tb <= 0] <- 0
  z
}
