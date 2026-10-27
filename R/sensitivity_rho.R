#' Sensitivity of the Indirect Effect to Mediator-Outcome Confounding
#'
#' @description
#' Implements the sensitivity analysis of Imai, Keele, and Yamamoto (2010) for
#' the linear structural equation model. The sensitivity parameter `rho` is the
#' correlation between the error terms of the mediator model and the outcome
#' model. Sequential ignorability implies `rho = 0`. A nonzero `rho` arises when
#' an unmeasured pre-treatment variable affects both the mediator and the
#' outcome.
#'
#' For each value of `rho`, the indirect effect per unit change in treatment is
#' \deqn{\delta(\rho) = \frac{\beta_2 \sigma_1}{\sigma_2}
#'   \left\{ \tilde\rho - \rho \sqrt{\frac{1 - \tilde\rho^2}{1 - \rho^2}} \right\}}
#' where \eqn{\beta_2} is the treatment coefficient in the mediator model,
#' \eqn{\sigma_2} is the residual standard deviation of the mediator model,
#' \eqn{\sigma_1} is the residual standard deviation of the outcome model with
#' the mediator removed, and \eqn{\tilde\rho} is the correlation between those
#' two sets of residuals. At `rho = 0` this equals the product of the treatment
#' coefficient in the mediator model and the mediator coefficient in the
#' outcome model. The indirect effect is zero exactly when `rho` equals
#' \eqn{\tilde\rho}.
#'
#' The same threshold can be stated in variance-explained terms. If an
#' unmeasured confounder explains a proportion \eqn{R^{2*}_M} of the residual
#' variance of the mediator and \eqn{R^{2*}_Y} of the residual variance of the
#' outcome, then \eqn{\rho^2 = R^{2*}_M R^{2*}_Y}.
#'
#' @details
#' The closed form holds only for the linear model. The function therefore
#' requires a Gaussian identity-link mediator and outcome model, a treatment
#' that enters both models as a single linear term, and no treatment-by-mediator
#' interaction. It treats the treatment assignment as ignorable given the
#' covariates; it addresses only confounding of the mediator-outcome
#' relationship. Results are point estimates; they do not include sampling
#' uncertainty.
#'
#' @param x A `robmedfit` object returned by [robustmediate()].
#' @param rho_seq Values of `rho` at which to evaluate the indirect effect.
#'   Default `seq(-0.9, 0.9, by = 0.01)`.
#' @param contrast Size of the treatment contrast. Default is the focal
#'   contrast reported by [robustmediate()] (focal dose minus reference dose).
#'
#' @return An object of class `"robmed_sens"`: a list with
#' \describe{
#'   \item{`curve`}{Data frame with columns `rho`, `nie_unit` (indirect effect
#'     per unit of treatment), and `nie` (indirect effect for `contrast`).}
#'   \item{`rho_zero`}{Value of `rho` at which the indirect effect is zero.}
#'   \item{`r2_product_zero`}{\eqn{R^{2*}_M R^{2*}_Y} at which the indirect
#'     effect is zero (equal to `rho_zero^2`).}
#'   \item{`nie_unit_0`, `contrast`}{Indirect effect per unit at `rho = 0`, and
#'     the contrast used.}
#' }
#'
#' @references
#' Imai, K., Keele, L., & Yamamoto, T. (2010). Identification, inference and
#' sensitivity analysis for causal mediation effects. *Statistical Science*,
#' 25(1), 51--71.
#'
#' Imai, K., Keele, L., & Tingley, D. (2010). A general approach to causal
#' mediation analysis. *Psychological Methods*, 15(4), 309--334.
#'
#' @examples
#' \donttest{
#' data(sim_mediation)
#' fit <- robustmediate(
#'   X ~ Z1 + Z2, M ~ X + Z1 + Z2, Y ~ X + M + Z1 + Z2,
#'   data = sim_mediation, R = 10, verbose = FALSE,
#'   dose_grid = seq(-2, 2, length.out = 25)
#' )
#' sens <- sensitivity_rho(fit)
#' sens
#' plot_sensitivity(fit)
#' }
#' @export
sensitivity_rho <- function(x, rho_seq = seq(-0.9, 0.9, by = 0.01),
                            contrast = NULL) {
  if (!inherits(x, "robmedfit")) rlang::abort("`x` must be a robmedfit object.")
  if (any(abs(rho_seq) >= 1)) rlang::abort("`rho_seq` must lie strictly inside (-1, 1).")

  tv <- x$meta$treat_var
  mv <- x$meta$med_var
  med_mod <- x$models$mediator
  out_mod <- x$models$outcome

  .check_linear_sem(med_mod, out_mod, tv, mv)

  # Outcome model with the mediator removed (same covariates, same data)
  red_formula <- stats::update(stats::formula(out_mod),
                               stats::as.formula(paste(". ~ . -", mv)))
  red_mod <- stats::lm(red_formula, data = stats::model.frame(out_mod))

  e_red <- stats::residuals(red_mod)
  e_med <- stats::residuals(med_mod)
  if (length(e_red) != length(e_med))
    rlang::abort("Mediator and outcome models were fit to different rows; remove missing values first.")

  beta2   <- unname(stats::coef(med_mod)[tv])
  sigma1  <- stats::sd(e_red)
  sigma2  <- stats::sd(e_med)
  rho_til <- stats::cor(e_red, e_med)

  nie_fun <- function(rho) {
    beta2 * sigma1 / sigma2 *
      (rho_til - rho * sqrt((1 - rho_til^2) / (1 - rho^2)))
  }

  if (is.null(contrast)) {
    contrast <- x$effects$summary$focal_dose - x$meta$ref_dose
  }

  curve <- data.frame(
    rho      = rho_seq,
    nie_unit = nie_fun(rho_seq),
    stringsAsFactors = FALSE
  )
  curve$nie <- curve$nie_unit * contrast

  structure(
    list(
      curve           = curve,
      rho_zero        = rho_til,
      r2_product_zero = rho_til^2,
      nie_unit_0      = nie_fun(0),
      contrast        = contrast
    ),
    class = "robmed_sens"
  )
}

# Check that the closed-form linear SEM result applies
.check_linear_sem <- function(med_mod, out_mod, tv, mv, fn = "sensitivity_rho") {
  is_lin <- function(m) {
    fam <- stats::family(m)
    fam$family == "gaussian" && fam$link == "identity"
  }
  if (!is_lin(med_mod) || !is_lin(out_mod)) {
    rlang::abort(paste0(
      fn, "() requires Gaussian identity-link mediator and outcome ",
      "models; the closed-form result does not hold for other families."
    ))
  }
  if (!tv %in% names(stats::coef(med_mod)) || !tv %in% names(stats::coef(out_mod)) ||
      !mv %in% names(stats::coef(out_mod))) {
    rlang::abort(paste0(
      fn, "() requires the treatment to enter the mediator and ",
      "outcome models, and the mediator to enter the outcome model, as single ",
      "linear terms."
    ))
  }
  terms_out <- attr(stats::terms(out_mod), "term.labels")
  terms_med <- attr(stats::terms(med_mod), "term.labels")
  uses <- function(lbls, v) any(vapply(lbls, function(l) v %in% all.vars(stats::as.formula(paste("~", l))), logical(1L)) & lbls != v)
  if (uses(terms_out, tv) || uses(terms_out, mv) || uses(terms_med, tv)) {
    rlang::abort(paste0(
      fn, "() requires the treatment and mediator to appear only ",
      "as main effects (no interactions, splines, or transformations)."
    ))
  }
  invisible(TRUE)
}

#' Print a robmed_sens object
#'
#' @param x A `robmed_sens` object from [sensitivity_rho()].
#' @param ... Ignored.
#' @return The input object, invisibly.
#' @export
print.robmed_sens <- function(x, ...) {
  cat("-- Sensitivity to mediator-outcome confounding (Imai et al., 2010) --\n")
  cat(sprintf("  Indirect effect per unit of treatment (rho = 0): %.4f\n", x$nie_unit_0))
  cat(sprintf("  Indirect effect for contrast of %.3f units:       %.4f\n",
              x$contrast, x$nie_unit_0 * x$contrast))
  cat(sprintf("  Indirect effect is zero at rho = %.3f\n", x$rho_zero))
  cat(sprintf("  Equivalent R2*_M x R2*_Y = %.4f\n", x$r2_product_zero))
  cat("  Point estimates only; assumes linear models and ignorable treatment.\n")
  invisible(x)
}

#' Plot the Indirect Effect as a Function of rho
#'
#' @description
#' Plots the indirect effect from [sensitivity_rho()] against the sensitivity
#' parameter `rho`, with the point where the indirect effect is zero marked.
#'
#' @param x A `robmedfit` object or a `robmed_sens` object.
#' @param ... Passed to [sensitivity_rho()] when `x` is a `robmedfit` object.
#'
#' @return A `ggplot2` object.
#'
#' @examples
#' \donttest{
#' data(sim_mediation)
#' fit <- robustmediate(
#'   X ~ Z1 + Z2, M ~ X + Z1 + Z2, Y ~ X + M + Z1 + Z2,
#'   data = sim_mediation, R = 10, verbose = FALSE,
#'   dose_grid = seq(-2, 2, length.out = 25)
#' )
#' plot_sensitivity(fit)
#' }
#' @importFrom ggplot2 ggplot aes geom_line geom_hline geom_vline labs
#'   theme_minimal
#' @export
plot_sensitivity <- function(x, ...) {
  if (inherits(x, "robmedfit")) {
    x <- if (!is.null(x$sensitivity) && length(list(...)) == 0L) x$sensitivity
         else sensitivity_rho(x, ...)
  }
  if (!inherits(x, "robmed_sens"))
    rlang::abort("`x` must be a robmedfit or robmed_sens object.")
  if (is.null(x$curve))
    rlang::abort("No sensitivity results available for this fit.")

  ggplot2::ggplot(x$curve, ggplot2::aes(x = rho, y = nie)) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.4) +
    ggplot2::geom_vline(xintercept = x$rho_zero, linetype = "dashed",
                        colour = "grey30", linewidth = 0.5) +
    ggplot2::geom_line(linewidth = 0.9) +
    ggplot2::labs(
      x = expression(rho ~ "(correlation of mediator and outcome errors)"),
      y = "Indirect effect",
      caption = paste0(
        "Contrast = ", round(x$contrast, 3), " units of treatment. ",
        "Dashed line: indirect effect is zero at rho = ", round(x$rho_zero, 3),
        ".\nImai, Keele & Yamamoto (2010); linear models; point estimates."
      )
    ) +
    ggplot2::theme_minimal(base_size = 11)
}
