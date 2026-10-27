# ==============================================================================
# effect_curvature
# ==============================================================================

#' Curvature of an Effect Curve, with Bootstrap Intervals
#'
#' @description
#' Describes the shape of an estimated effect curve (NIE, NDE, or TE) over the
#' treatment grid: the second derivative of the curve at each interior dose,
#' with percentile bootstrap intervals taken from the same bootstrap
#' replicates as the effect estimates.
#'
#' This is a description of the dose-response shape. It is not a sensitivity
#' analysis for unmeasured confounding; for that see [sensitivity_rho()] and
#' [sensitivity_paths()].
#'
#' @details
#' When the mediator and outcome models are linear (identity links) and the
#' treatment and mediator enter only as main effects, every effect curve is a
#' straight line by construction. The function then returns curvature exactly
#' zero, with no interval, and says so in a message, instead of differencing
#' rounding error. To let the curve bend, put a spline of the treatment in the
#' model formulas, for example `M ~ splines::ns(X, df = 3) + Z1`.
#'
#' Otherwise the second derivative is computed by second differences on the
#' (possibly unequally spaced) dose grid, for the point estimate and for each
#' bootstrap replicate. Values at the two end doses are `NA`.
#'
#' The curvature is that of the fitted models, so it depends on the spline
#' chosen. Natural splines are linear beyond their boundary knots, so with few
#' degrees of freedom the curvature is pulled toward zero near the ends of the
#' treatment range. The intervals reflect sampling variability only, not this
#' modelling choice; compare a few spline settings before interpreting shape. Fits created
#' before version 0.3.0 do not store replicate curves, so intervals are `NA`.
#'
#' @param x A `robmedfit` object.
#' @param estimand `"NIE"` (default), `"NDE"`, or `"TE"`.
#' @param alpha Level for the curvature intervals. Defaults to the level used
#'   in the fit.
#'
#' @return A data frame of class `"robmed_curvature"` with columns `dose`,
#'   `estimate`, `lower`, `upper` (the effect curve and its interval),
#'   `curvature`, `curv_lower`, `curv_upper`, and `ci_includes_zero` (whether
#'   the effect interval includes zero). Attributes `estimand` and
#'   `linear_by_construction` record the settings.
#'
#' @seealso [plot_curvature()]
#'
#' @examples
#' \donttest{
#' data(sim_mediation)
#' fit <- robustmediate(
#'   X ~ Z1 + Z2, M ~ splines::ns(X, df = 3) + Z1 + Z2, Y ~ X + M + Z1 + Z2,
#'   data = sim_mediation, R = 10, verbose = FALSE,
#'   dose_grid = seq(-2, 2, length.out = 25)
#' )
#' curv <- effect_curvature(fit, estimand = "NIE")
#' plot_curvature(curv, ref_dose = fit$meta$ref_dose)
#' }
#'
#' @export
effect_curvature <- function(x, estimand = c("NIE", "NDE", "TE"), alpha = NULL) {
  if (!inherits(x, "robmedfit")) rlang::abort("`x` must be a robmedfit object.")
  estimand <- match.arg(estimand)
  if (is.null(alpha)) alpha <- x$meta$alpha

  curves <- x$effects$curves
  d      <- curves[curves$estimand == estimand, , drop = FALSE]
  ord    <- order(d$dose)
  d      <- d[ord, c("dose", "estimate", "lower", "upper")]
  rownames(d) <- NULL
  n <- nrow(d)

  linear <- .linear_in_dose(x$models, x$meta$treat_var, x$meta$med_var)

  d$curvature <- NA_real_
  d$curv_lower <- NA_real_
  d$curv_upper <- NA_real_

  if (linear) {
    d$curvature <- 0
    message(
      "The ", estimand, " curve is a straight line by construction (linear ",
      "models, treatment and mediator as main effects), so its curvature is ",
      "zero. Add a spline of the treatment to the formulas to allow a curve."
    )
  } else if (n >= 3) {
    d$curvature <- .second_diff(d$dose, d$estimate)
    boot <- x$effects$boot
    if (!is.null(boot)) {
      reps <- boot[ord, estimand, , drop = TRUE]
      if (is.null(dim(reps))) reps <- matrix(reps, ncol = 1L)
      ok <- apply(reps, 2L, function(v) all(is.finite(v)))
      reps <- reps[, ok, drop = FALSE]
      if (ncol(reps) >= 2L) {
        k2 <- apply(reps, 2L, function(v) .second_diff(d$dose, v))
        d$curv_lower <- apply(k2, 1L, stats::quantile, probs = alpha / 2,
                              na.rm = TRUE, names = FALSE)
        d$curv_upper <- apply(k2, 1L, stats::quantile, probs = 1 - alpha / 2,
                              na.rm = TRUE, names = FALSE)
        d$curv_lower[c(1L, n)] <- NA_real_
        d$curv_upper[c(1L, n)] <- NA_real_
      }
    }
  }

  d$ci_includes_zero <- d$lower <= 0 & d$upper >= 0
  structure(d, class = c("robmed_curvature", "data.frame"),
            estimand = estimand, linear_by_construction = linear)
}


#' @rdname effect_curvature
#' @description
#' `sensitivity_curvature()` is deprecated; use `effect_curvature()`.
#' @export
sensitivity_curvature <- function(x, estimand = c("NIE", "NDE", "TE")) {
  .Deprecated("effect_curvature", package = "RobustMediate",
              msg = paste0("sensitivity_curvature() is deprecated; use ",
                           "effect_curvature(). It describes the shape of the ",
                           "effect curve and is not a sensitivity analysis."))
  effect_curvature(x, estimand = estimand)
}


# Second differences on a possibly unequal grid; NA at the two ends
.second_diff <- function(x, y) {
  n <- length(y)
  out <- rep(NA_real_, n)
  if (n < 3L) return(out)
  h1 <- x[2:(n - 1)] - x[1:(n - 2)]
  h2 <- x[3:n] - x[2:(n - 1)]
  out[2:(n - 1)] <- 2 * (y[3:n] / h2 - y[2:(n - 1)] * (1 / h1 + 1 / h2) +
                           y[1:(n - 2)] / h1) / (h1 + h2)
  out
}

# TRUE when every effect curve is linear in the dose by construction
.linear_in_dose <- function(models, treat_var, med_var) {
  idl <- function(m) stats::family(m)$link == "identity"
  only_main <- function(model, v) {
    lbls <- attr(stats::terms(model), "term.labels")
    uses <- vapply(lbls, function(l)
      v %in% all.vars(stats::as.formula(paste("~", l))), logical(1L))
    all(lbls[uses] == v)
  }
  idl(models$mediator) && idl(models$outcome) &&
    only_main(models$mediator, treat_var) &&
    only_main(models$outcome, treat_var) &&
    only_main(models$outcome, med_var)
}


# ==============================================================================
# plot_curvature
# ==============================================================================

#' Plot an Effect Curve and Its Curvature
#'
#' @description
#' Two panels: the effect curve with its bootstrap interval, and its second
#' derivative with a bootstrap interval. Doses where the effect interval
#' includes zero are shaded.
#'
#' @param x Output of [effect_curvature()].
#' @param estimand Label for the estimand. Defaults to the one stored in `x`.
#' @param ref_dose Optional reference dose, drawn as a vertical line.
#' @param ... Ignored.
#'
#' @return A `ggplot2` object.
#'
#' @examples
#' \donttest{
#' data(sim_mediation)
#' fit <- robustmediate(
#'   X ~ Z1 + Z2, M ~ splines::ns(X, df = 3) + Z1 + Z2, Y ~ X + M + Z1 + Z2,
#'   data = sim_mediation, R = 10, verbose = FALSE,
#'   dose_grid = seq(-2, 2, length.out = 25)
#' )
#' plot_curvature(effect_curvature(fit), ref_dose = fit$meta$ref_dose)
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_ribbon geom_line geom_rect geom_hline
#'   geom_vline labs theme_minimal theme facet_wrap element_text
#' @export
plot_curvature <- function(x, estimand = NULL, ref_dose = NULL, ...) {
  if (!is.data.frame(x) || !all(c("curvature", "curv_lower") %in% names(x))) {
    rlang::abort("`x` must be the output of effect_curvature().")
  }
  if (is.null(estimand)) estimand <- attr(x, "estimand") %||% "NIE"

  p1 <- paste0(estimand, " estimate")
  p2 <- paste0("Curvature of ", estimand, " (second derivative)")
  panel_df <- rbind(
    data.frame(dose = x$dose, y = x$estimate, lower = x$lower,
               upper = x$upper, panel = p1, stringsAsFactors = FALSE),
    data.frame(dose = x$dose, y = x$curvature, lower = x$curv_lower,
               upper = x$curv_upper, panel = p2, stringsAsFactors = FALSE)
  )
  panel_df$panel <- factor(panel_df$panel, levels = c(p1, p2))

  step <- if (nrow(x) > 1L) min(diff(sort(x$dose))) / 2 else 0.01
  zero_df <- x[x$ci_includes_zero %in% TRUE, , drop = FALSE]

  p <- ggplot2::ggplot(panel_df, ggplot2::aes(x = dose)) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed",
                        colour = "grey60", linewidth = 0.4)
  if (nrow(zero_df) > 0L) {
    p <- p + ggplot2::geom_rect(
      data = zero_df,
      ggplot2::aes(xmin = dose - step, xmax = dose + step,
                   ymin = -Inf, ymax = Inf),
      fill = "grey85", alpha = 0.5, inherit.aes = FALSE
    )
  }
  rib <- panel_df[is.finite(panel_df$lower) & is.finite(panel_df$upper), ]
  if (nrow(rib) > 0L) {
    p <- p + ggplot2::geom_ribbon(
      data = rib, ggplot2::aes(ymin = lower, ymax = upper),
      fill = "grey50", alpha = 0.25
    )
  }
  p <- p +
    ggplot2::geom_line(ggplot2::aes(y = y), linewidth = 0.9, na.rm = TRUE) +
    ggplot2::facet_wrap(~ panel, ncol = 1, scales = "free_y") +
    ggplot2::labs(
      x = "Treatment dose", y = NULL,
      caption = paste0(
        "Bands: bootstrap intervals. Grey shading: doses where the ",
        estimand, " interval includes zero.",
        if (isTRUE(attr(x, "linear_by_construction")))
          "\nThe curve is linear by construction, so its curvature is zero." else ""
      )
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      strip.text   = ggplot2::element_text(face = "bold", size = 9),
      plot.caption = ggplot2::element_text(size = 8, colour = "grey40")
    )

  if (!is.null(ref_dose)) {
    p <- p + ggplot2::geom_vline(xintercept = ref_dose, linetype = "dotted",
                                 colour = "grey40", linewidth = 0.5)
  }
  p
}

`%||%` <- function(a, b) if (is.null(a)) b else a
