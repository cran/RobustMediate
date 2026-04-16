# Suppress R CMD check NOTEs for ggplot2 aes() column-name variables.
# These are evaluated inside ggplot2's non-standard evaluation environment.
utils::globalVariables(c(
  # plot_mediation / compare_fits
  "dose", "estimate", "estimand", "lower", "upper", "model",
  # plot_balance
  "smd", "covariate", "timing",
  # plot_meditcv
  "corridor_lo", "corridor_hi", "r_crit", "r_obs",
  "x_pos", "exceeds_meditcv", "label_txt",
  # plot_meditcv_profile
  "delta", "r_attenuated", "pathway", "r_att",
  # plot_sensitivity
  "evalue", "rho", "effect",
  # plot_curvature
  "y",
  # internal package functions called across files
  "sensitivity_meditcv", "plot_meditcv", "sensitivity_meditcv_profile"
))


# ==============================================================================
# .meditcv_pathway — internal helper
# ==============================================================================

.meditcv_pathway <- function(model, predictor_var, alpha, label) {
  co <- stats::coef(model)
  vc <- stats::vcov(model)
  df <- stats::df.residual(model)

  if (!predictor_var %in% names(co)) {
    rlang::abort(
      paste0("Predictor '", predictor_var, "' not found in model coefficients.")
    )
  }

  b      <- unname(co[predictor_var])
  se     <- sqrt(diag(vc)[predictor_var])
  t_stat <- b / se

  r_obs  <- sign(t_stat) * sqrt(t_stat^2 / (t_stat^2 + df))
  t_crit <- stats::qt(1 - alpha / 2, df = df)
  r_crit <- t_crit / sqrt(t_crit^2 + df)

  meditcv     <- abs(r_obs) - r_crit
  meditcv_pct <- 100 * meditcv / (1 - r_crit)

  bm_r <- c(0.1, 0.2, 0.3, 0.5)
  benchmarks <- data.frame(
    r_confounder    = bm_r,
    impact          = bm_r^2,
    exceeds_meditcv = bm_r^2 > abs(meditcv),
    stringsAsFactors = FALSE
  )

  list(
    label       = label,
    predictor   = predictor_var,
    coef        = b,
    se          = se,
    t_stat      = t_stat,
    df          = df,
    t_crit      = t_crit,
    r_obs       = r_obs,
    r_crit      = r_crit,
    meditcv     = meditcv,
    meditcv_pct = meditcv_pct,
    significant = abs(t_stat) > t_crit,
    benchmarks  = benchmarks
  )
}


# ==============================================================================
# sensitivity_meditcv
# ==============================================================================

#' Mediation ITCV (medITCV) for pathway-specific robustness
#'
#' @description
#' Computes a mediation-specific extension of Kenneth Frank's (2000)
#' Impact Threshold for a Confounding Variable (ITCV) for both pathways of a
#' mediation model:
#'
#' - **a-path**: treatment -> mediator
#' - **b-path**: mediator -> outcome (controlling for treatment)
#'
#' The mediation ITCV (medITCV) quantifies how strong an unmeasured confounder
#' would need to be, in terms of the product \eqn{r_{XC} \cdot r_{YC}}, to
#' invalidate inference for each pathway.
#'
#' @param x     A `robmedfit` object returned by [robustmediate()].
#' @param alpha Significance level. Default is `0.05`.
#'
#' @return An object of class `"meditcv"`: a named list with elements
#'   `a_path`, `b_path`, `indirect`, and `alpha`. Each pathway element
#'   contains the observed partial correlation, critical partial correlation,
#'   medITCV value, and benchmark confounder impacts.
#'
#' @references
#' Frank, K. A. (2000). Impact of a confounding variable on a regression
#' coefficient. *Sociological Methods & Research*, 29(2), 147--194.
#'
#' @seealso `plot_meditcv()`, `print.meditcv()`
#'
#' @examples
#' \donttest{
#' data(sim_mediation)
#'   fit <- robustmediate(
#'     X ~ Z1 + Z2, M ~ X + Z1 + Z2, Y ~ X + M + Z1 + Z2,
#'     data = sim_mediation, R = 20, verbose = FALSE
#'   )
#' med <- sensitivity_meditcv(fit)
#' print(med)
#' plot_meditcv(med)
#' }
#'
#' @export
sensitivity_meditcv <- function(x, alpha = 0.05) {
  if (!inherits(x, "robmedfit")) {
    rlang::abort("`x` must be a `robmedfit` object.")
  }

  treat_var <- x$meta$treat_var
  med_var   <- x$meta$med_var

  a_res <- .meditcv_pathway(
    model         = x$models$mediator,
    predictor_var = treat_var,
    alpha         = alpha,
    label         = "a-path (X -> M)"
  )

  b_res <- .meditcv_pathway(
    model         = x$models$outcome,
    predictor_var = med_var,
    alpha         = alpha,
    label         = "b-path (M -> Y | X)"
  )

  med_indirect <- min(a_res$meditcv, b_res$meditcv)

  bottleneck <- if (a_res$meditcv <= b_res$meditcv) {
    "a-path (X -> M)"
  } else {
    "b-path (M -> Y | X)"
  }

  pct_a <- if (isTRUE(all.equal(a_res$r_obs, 0))) {
    NA_real_
  } else {
    100 * a_res$meditcv / abs(a_res$r_obs)
  }

  pct_b <- if (isTRUE(all.equal(b_res$r_obs, 0))) {
    NA_real_
  } else {
    100 * b_res$meditcv / abs(b_res$r_obs)
  }

  structure(
    list(
      a_path   = a_res,
      b_path   = b_res,
      indirect = list(
        meditcv_indirect = med_indirect,
        bottleneck       = bottleneck,
        pct_replace_a    = pct_a,
        pct_replace_b    = pct_b
      ),
      alpha = alpha
    ),
    class = "meditcv"
  )
}


# ==============================================================================
# print.meditcv
# ==============================================================================

#' Print a meditcv object
#'
#' @param x   A `meditcv` object from `sensitivity_meditcv()`.
#' @param ... Ignored.
#'
#' @return The input object, invisibly. Called for its side effect of printing
#'   a formatted medITCV report to the console.
#'
#' @export
print.meditcv <- function(x, ...) {
  cat("\u2554", strrep("\u2550", 62), "\u2557\n", sep = "")
  cat("\u2551   medITCV: Mediation Robustness Report                      \u2551\n")
  cat("\u255a", strrep("\u2550", 62), "\u255d\n\n", sep = "")

  for (pwy in list(x$a_path, x$b_path)) {
    sig_flag <- if (pwy$significant) "\u2713 significant" else "\u26a0 NOT significant"

    cat(sprintf("\u2500\u2500 %s (%s) \u2500\u2500\n", pwy$label, sig_flag))
    cat(sprintf(
      "  Coefficient : %+.4f  (SE = %.4f, t = %.3f, df = %d)\n",
      pwy$coef, pwy$se, pwy$t_stat, pwy$df
    ))
    cat(sprintf(
      "  Observed partial r : %.4f   |   Critical r (alpha = %.2f) : %.4f\n",
      pwy$r_obs, x$alpha, pwy$r_crit
    ))

    if (pwy$meditcv >= 0) {
      cat(sprintf(
        "  medITCV = %.4f  (%+.1f%% above threshold)\n",
        pwy$meditcv, pwy$meditcv_pct
      ))
      cat(sprintf(
        "  -> A confounder needs |r_XC x r_YC| > %.4f to invalidate this pathway.\n",
        pwy$meditcv
      ))
    } else {
      cat(sprintf(
        "  medITCV = %.4f  \u2014 inference is already below the threshold.\n",
        pwy$meditcv
      ))
    }

    bm <- pwy$benchmarks
    cat("\n  Benchmark impacts (r_confounder^2):\n")
    for (i in seq_len(nrow(bm))) {
      flag <- if (bm$exceeds_meditcv[i]) {
        "  \u26a0 WOULD INVALIDATE"
      } else {
        "  \u2713 would not invalidate"
      }
      cat(sprintf(
        "    r = %.1f  ->  impact = %.2f%s\n",
        bm$r_confounder[i], bm$impact[i], flag
      ))
    }
    cat("\n")
  }

  ind <- x$indirect
  cat("\u2500\u2500 Indirect effect (a x b) \u2500\u2500\n")
  cat(sprintf("  Bottleneck pathway    : %s\n", ind$bottleneck))
  cat(sprintf("  Minimum-path medITCV  : %.4f\n", ind$meditcv_indirect))
  cat("  Fragility of the indirect effect is governed by the weaker pathway.\n")
  cat(strrep("\u2500", 64), "\n", sep = "")

  invisible(x)
}


# ==============================================================================
# plot_meditcv
# ==============================================================================

#' Plot medITCV robustness corridors for both mediation pathways
#'
#' @description
#' Produces a two-panel pathway-specific robustness corridor plot showing the
#' observed partial correlation, critical partial correlation threshold, medITCV
#' corridor, and benchmark confounder impacts for each pathway.
#'
#' @param x   A `meditcv` object from `sensitivity_meditcv()`.
#' @param ... Ignored.
#'
#' @return A `ggplot2` object.
#'
#' @examples
#' \donttest{
#' data(sim_mediation)
#'   fit <- robustmediate(
#'     X ~ Z1 + Z2, M ~ X + Z1 + Z2, Y ~ X + M + Z1 + Z2,
#'     data = sim_mediation, R = 20, verbose = FALSE
#'   )
#' med <- sensitivity_meditcv(fit)
#' plot_meditcv(med)
#' }
#'
#' @importFrom ggplot2 aes facet_wrap geom_point geom_rect geom_segment
#'   geom_text geom_vline ggplot labs scale_colour_manual scale_x_continuous
#'   theme theme_minimal element_blank element_text
#' @export
plot_meditcv <- function(x, ...) {
  if (!inherits(x, "meditcv")) {
    rlang::abort("`x` must be a `meditcv` object.")
  }

  panels <- lapply(list(x$a_path, x$b_path), function(pwy) {
    r_obs         <- pwy$r_obs
    r_crit_signed <- pwy$r_crit * sign(r_obs)
    bm            <- pwy$benchmarks
    bm$x_pos      <- sqrt(bm$impact) * sign(r_obs)

    list(
      label       = pwy$label,
      r_obs       = r_obs,
      r_crit      = r_crit_signed,
      corridor_lo = min(r_crit_signed, r_obs),
      corridor_hi = max(r_crit_signed, r_obs),
      benchmarks  = bm
    )
  })

  path_df <- do.call(rbind, lapply(panels, function(p) {
    data.frame(
      pathway     = p$label,
      r_obs       = p$r_obs,
      r_crit      = p$r_crit,
      corridor_lo = p$corridor_lo,
      corridor_hi = p$corridor_hi,
      stringsAsFactors = FALSE
    )
  }))

  bm_df <- do.call(rbind, lapply(panels, function(p) {
    cbind(data.frame(pathway = p$label, stringsAsFactors = FALSE), p$benchmarks)
  }))
  bm_df$label_txt <- paste0("r=", bm_df$r_confounder)

  x_range <- range(c(path_df$r_obs, path_df$r_crit, bm_df$x_pos, 0), na.rm = TRUE)
  x_pad   <- diff(x_range) * 0.15
  if (!is.finite(x_pad) || x_pad == 0) x_pad <- 0.1
  x_lims  <- c(x_range[1] - x_pad, x_range[2] + x_pad)

  ggplot2::ggplot() +
    ggplot2::geom_rect(
      data = path_df,
      ggplot2::aes(
        xmin = corridor_lo, xmax = corridor_hi,
        ymin = -Inf,        ymax = Inf
      ),
      fill  = "#FFF3CD",
      alpha = 0.8
    ) +
    ggplot2::geom_vline(
      xintercept = 0,
      colour     = "grey70",
      linewidth  = 0.4
    ) +
    ggplot2::geom_vline(
      data      = path_df,
      ggplot2::aes(xintercept = r_crit),
      linetype  = "dashed",
      colour    = "#C0392B",
      linewidth = 0.7
    ) +
    ggplot2::geom_point(
      data   = path_df,
      ggplot2::aes(x = r_obs, y = 0.5),
      colour = "#2C3E7A",
      size   = 4,
      shape  = 18
    ) +
    ggplot2::geom_segment(
      data = path_df,
      ggplot2::aes(x = r_obs, xend = r_obs, y = 0.3, yend = 0.7),
      colour    = "#2C3E7A",
      linewidth = 1
    ) +
    ggplot2::geom_point(
      data = bm_df,
      ggplot2::aes(x = x_pos, y = -0.3, colour = exceeds_meditcv),
      size  = 3,
      shape = 19
    ) +
    ggplot2::geom_text(
      data = bm_df,
      ggplot2::aes(x = x_pos, y = -0.6, label = label_txt),
      size   = 2.8,
      colour = "grey30"
    ) +
    ggplot2::scale_colour_manual(
      values = c(`TRUE` = "#C0392B", `FALSE` = "#27AE60"),
      labels = c(`TRUE` = "Would invalidate", `FALSE` = "Would not invalidate"),
      name   = "Benchmark confounder"
    ) +
    ggplot2::scale_x_continuous(limits = x_lims) +
    ggplot2::facet_wrap(~ pathway, ncol = 1, scales = "free_x") +
    ggplot2::labs(
      title    = "Pathway-Specific medITCV Robustness Corridor",
      subtitle = paste0(
        "Diamond = observed r  |  Red dashed = critical r (alpha = ", x$alpha, ")",
        "  |  Yellow band = medITCV corridor\n",
        "Dots = benchmark confounders at r = 0.1, 0.2, 0.3, 0.5  (impact = r\u00b2)"
      ),
      caption = "Based on Frank's ITCV framework (Frank, 2000).",
      x = "Partial correlation",
      y = NULL
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      axis.text.y        = ggplot2::element_blank(),
      axis.ticks.y       = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor.y = ggplot2::element_blank(),
      legend.position    = "bottom",
      plot.subtitle      = ggplot2::element_text(size = 8, colour = "grey40"),
      strip.text         = ggplot2::element_text(face = "bold")
    )
}


# ==============================================================================
# sensitivity_meditcv_profile
# ==============================================================================

#' medITCV Robustness Profile: Pathway-Specific Fragility Framework
#'
#' @description
#' Implements the medITCV robustness profile, a formal extension of Frank's
#' ITCV to causal mediation. Computes pathway-specific fragility thresholds,
#' applies the minimum robustness principle, and identifies the bottleneck
#' pathway that governs indirect-effect fragility.
#'
#' @param x          A `robmedfit` object.
#' @param alpha      Significance level. Default `0.05`.
#' @param delta_grid Numeric vector of confounding impact values over which the
#'   robustness profile is evaluated. Default `seq(0, 0.5, by = 0.01)`.
#'
#' @return An object of class `"meditcv_profile"`: a named list with elements
#'   `a_path`, `b_path`, `meditcv_indirect`, `bottleneck`,
#'   `robustness_profile`, `fragility_ratio`, `meditcv_detail`, and `alpha`.
#'
#' @references
#' Frank, K. A. (2000). Impact of a confounding variable on a regression
#' coefficient. *Sociological Methods & Research*, 29(2), 147--194.
#'
#' @seealso `sensitivity_meditcv()`, `plot_meditcv_profile()`,
#'   `sensitivity_curvature()`
#'
#' @examples
#' \donttest{
#' data(sim_mediation)
#'   fit <- robustmediate(
#'     X ~ Z1 + Z2, M ~ X + Z1 + Z2, Y ~ X + M + Z1 + Z2,
#'     data = sim_mediation, R = 20, verbose = FALSE
#'   )
#' mp  <- sensitivity_meditcv_profile(fit)
#' print(mp)
#' plot_meditcv_profile(mp)
#' }
#'
#' @export
sensitivity_meditcv_profile <- function(x, alpha = 0.05,
                                        delta_grid = seq(0, 0.5, by = 0.01)) {
  if (!inherits(x, "robmedfit")) rlang::abort("`x` must be a robmedfit object.")

  meditcv_obj <- if (!is.null(x$meditcv)) {
    x$meditcv
  } else {
    sensitivity_meditcv(x, alpha = alpha)
  }

  a <- meditcv_obj$a_path
  b <- meditcv_obj$b_path

  meditcv_a <- a$meditcv
  meditcv_b <- b$meditcv

  meditcv_ind <- min(meditcv_a, meditcv_b)
  bottleneck  <- if (meditcv_a <= meditcv_b) "a-path (X \u2192 M)" else "b-path (M \u2192 Y | X)"
  frag_ratio  <- meditcv_a / max(meditcv_b, 1e-10)

  profile_rows <- lapply(delta_grid, function(d) {
    r_a_att <- a$r_obs - sign(a$r_obs) * d
    r_b_att <- b$r_obs - sign(b$r_obs) * d
    data.frame(
      delta        = d,
      pathway      = c("a-path (X\u2192M)", "b-path (M\u2192Y|X)"),
      r_attenuated = c(r_a_att, r_b_att),
      r_crit       = c(a$r_crit, b$r_crit),
      significant  = c(abs(r_a_att) > a$r_crit, abs(r_b_att) > b$r_crit),
      stringsAsFactors = FALSE
    )
  })
  profile <- do.call(rbind, profile_rows)
  rownames(profile) <- NULL

  structure(
    list(
      a_path             = meditcv_a,
      b_path             = meditcv_b,
      meditcv_indirect   = meditcv_ind,
      bottleneck         = bottleneck,
      robustness_profile = profile,
      fragility_ratio    = frag_ratio,
      meditcv_detail     = meditcv_obj,
      alpha              = alpha
    ),
    class = "meditcv_profile"
  )
}


# ==============================================================================
# print.meditcv_profile
# ==============================================================================

#' Print a meditcv_profile object
#'
#' @param x   A `meditcv_profile` object from `sensitivity_meditcv_profile()`.
#' @param ... Ignored.
#'
#' @return The input object, invisibly. Called for its side effect of printing
#'   a formatted medITCV robustness profile to the console.
#'
#' @export
print.meditcv_profile <- function(x, ...) {
  a_det <- x$meditcv_detail$a_path
  b_det <- x$meditcv_detail$b_path

  cat("\u2554", strrep("\u2550", 62), "\u2557\n", sep = "")
  cat("\u2551      medITCV: Mediation Pathway Fragility Framework        \u2551\n")
  cat("\u255a", strrep("\u2550", 62), "\u255d\n\n", sep = "")

  cat("\u2500\u2500 Minimum Robustness Principle \u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\n")
  cat(sprintf("  medITCV (a-path):        %.4f\n", x$a_path))
  cat(sprintf("  medITCV (b-path):        %.4f\n", x$b_path))
  cat(sprintf("  medITCV (indirect):      %.4f   [= min(a, b)]\n", x$meditcv_indirect))
  cat(sprintf("  Fragility ratio (a/b):   %.3f\n", x$fragility_ratio))

  frag_interp <- if (x$fragility_ratio < 0.5) {
    "  \u2192 Strong asymmetry: a-path is the clear bottleneck."
  } else if (x$fragility_ratio < 1.0) {
    "  \u2192 Moderate asymmetry: a-path is weaker but b-path is not far behind."
  } else if (x$fragility_ratio < 2.0) {
    "  \u2192 Moderate asymmetry: b-path is weaker \u2014 it governs overall fragility."
  } else {
    "  \u2192 Strong asymmetry: b-path is the clear bottleneck."
  }
  cat(frag_interp, "\n\n")

  cat("\u2500\u2500 Bottleneck Pathway \u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\n")
  cat(sprintf("  %s\n", x$bottleneck))
  cat("  Robustness of the indirect effect is bounded by this pathway.\n")
  cat("  Collect additional covariates / run balance checks here first.\n\n")

  cat("\u2500\u2500 Path Detail \u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\n")
  a_class <- if (x$a_path > 0.1) "robust" else if (x$a_path > 0) "fragile" else "not sig"
  b_class <- if (x$b_path > 0.1) "robust" else if (x$b_path > 0) "fragile" else "not sig"
  cat(sprintf("  a-path: r_obs = %+.4f  r_crit = %.4f  medITCV = %.4f  (%s)\n",
              a_det$r_obs, a_det$r_crit, x$a_path, a_class))
  cat(sprintf("  b-path: r_obs = %+.4f  r_crit = %.4f  medITCV = %.4f  (%s)\n",
              b_det$r_obs, b_det$r_crit, x$b_path, b_class))

  cat("\n\u2500\u2500 Tipping-Point Confounders \u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\n")
  for (nm in c("a-path", "b-path")) {
    val   <- if (nm == "a-path") x$a_path else x$b_path
    r_tip <- sqrt(max(val, 0))
    cat(sprintf("  %s: a confounder correlated r = %.3f with BOTH predictor and outcome\n",
                nm, r_tip))
    cat("       would be sufficient to invalidate this pathway.\n")
  }

  cat("\n\u2500\u2500 Publishable Language \u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\n")
  cat(sprintf(paste0(
    "  \"Extending ITCV to mediation reveals that robustness is not\n",
    "   uniform but pathway-specific (medITCV_a = %.3f, medITCV_b = %.3f),\n",
    "   with the %s determining overall indirect-effect\n",
    "   fragility (medITCV_indirect = %.3f). A confounding impact exceeding\n",
    "   this threshold would be sufficient to nullify the indirect effect,\n",
    "   regardless of how robust the complementary pathway is.\"\n"
  ), x$a_path, x$b_path, x$bottleneck, x$meditcv_indirect))

  cat(strrep("\u2500", 64), "\n", sep = "")
  invisible(x)
}


# ==============================================================================
# plot_meditcv_profile
# ==============================================================================

#' Plot the medITCV Robustness Profile
#'
#' @description
#' Visualises how each pathway's partial correlation is attenuated as
#' confounding impact delta increases, with tipping points and a fragility
#' zone marked.
#'
#' @param x   A `meditcv_profile` object from `sensitivity_meditcv_profile()`.
#' @param ... Ignored.
#'
#' @return A `ggplot2` object.
#'
#' @examples
#' \donttest{
#' data(sim_mediation)
#'   fit <- robustmediate(
#'     X ~ Z1 + Z2, M ~ X + Z1 + Z2, Y ~ X + M + Z1 + Z2,
#'     data = sim_mediation, R = 20, verbose = FALSE
#'   )
#' mp  <- sensitivity_meditcv_profile(fit)
#' plot_meditcv_profile(mp)
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_line geom_hline geom_vline geom_rect
#'   geom_point scale_color_manual labs theme_minimal theme element_text
#'   scale_x_continuous
#' @export
plot_meditcv_profile <- function(x, ...) {
  if (!inherits(x, "meditcv_profile")) {
    rlang::abort("`x` must be a `meditcv_profile` object.")
  }

  prof   <- x$robustness_profile
  r_crit <- x$meditcv_detail$a_path$r_crit

  pal <- c("a-path (X\u2192M)" = "#534AB7", "b-path (M\u2192Y|X)" = "#0F6E56")

  tip_a   <- x$a_path
  tip_b   <- x$b_path
  tip_min <- x$meditcv_indirect

  ggplot2::ggplot(prof,
                  ggplot2::aes(x = delta, y = abs(r_attenuated), colour = pathway)) +
    ggplot2::geom_rect(
      xmin = tip_min, xmax = max(prof$delta),
      ymin = -Inf,    ymax = Inf,
      fill = "#FDECEA", alpha = 0.4, inherit.aes = FALSE
    ) +
    ggplot2::geom_hline(yintercept = r_crit, linetype = "dashed",
                        colour = "#C0392B", linewidth = 0.6) +
    ggplot2::geom_vline(xintercept = tip_a, linetype = "dotted",
                        colour = pal["a-path (X\u2192M)"], linewidth = 0.8) +
    ggplot2::geom_vline(xintercept = tip_b, linetype = "dotted",
                        colour = pal["b-path (M\u2192Y|X)"], linewidth = 0.8) +
    ggplot2::geom_line(linewidth = 1.1) +
    ggplot2::geom_point(
      data = data.frame(
        delta   = c(tip_a, tip_b),
        r_att   = c(r_crit, r_crit),
        pathway = c("a-path (X\u2192M)", "b-path (M\u2192Y|X)")
      ),
      ggplot2::aes(x = delta, y = r_att, colour = pathway),
      size = 3.5, shape = 21, fill = "white", stroke = 1.5,
      inherit.aes = FALSE
    ) +
    ggplot2::scale_color_manual(values = pal, name = "Pathway") +
    ggplot2::scale_x_continuous(
      breaks = sort(unique(round(c(0, tip_a, tip_b, 0.1, 0.2, 0.3, 0.4, 0.5), 3)))
    ) +
    ggplot2::labs(
      x        = expression(paste(
        "Confounding impact  ", delta, " = ", r[XC], " \u00d7 ", r[YC]
      )),
      y        = "Attenuated |partial r|",
      title    = "medITCV Robustness Profile",
      subtitle = paste0(
        "medITCV_a = ", round(tip_a, 4),
        "   medITCV_b = ", round(tip_b, 4),
        "   medITCV_indirect = ", round(tip_min, 4),
        "   Bottleneck: ", x$bottleneck
      ),
      caption = paste0(
        "Red dashed = critical r (alpha = ", x$alpha, ").  ",
        "Open circles = tipping points.  ",
        "Red shaded zone = indirect effect invalidated."
      )
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      legend.position = "top",
      plot.subtitle   = ggplot2::element_text(size = 8.5, colour = "grey40")
    )
}


# ==============================================================================
# sensitivity_curvature
# ==============================================================================

#' Dose-Varying Fragility: Curvature-Based Sensitivity
#'
#' @description
#' Computes the fragility curvature of the mediation effect across the full
#' treatment dose grid. Returns local fragility index, numerical curvature of
#' the effect curve, and a fragility zone flag at each dose value.
#'
#' @param x        A `robmedfit` object.
#' @param estimand Which estimand to analyse: `"NIE"` (default), `"NDE"`,
#'   or `"TE"`.
#'
#' @return A data frame with columns `dose`, `estimate`, `lower`, `upper`,
#'   `se_approx`, `frag_local`, `curvature`, and `in_fragility_zone`.
#'
#' @seealso `plot_curvature()`, `sensitivity_meditcv_profile()`
#'
#' @examples
#' \donttest{
#' data(sim_mediation)
#'   fit <- robustmediate(
#'     X ~ Z1 + Z2, M ~ X + Z1 + Z2, Y ~ X + M + Z1 + Z2,
#'     data = sim_mediation, R = 20, verbose = FALSE
#'   )
#' curv <- sensitivity_curvature(fit, estimand = "NIE")
#' plot_curvature(curv)
#' }
#'
#' @export
sensitivity_curvature <- function(x, estimand = c("NIE", "NDE", "TE")) {
  if (!inherits(x, "robmedfit")) rlang::abort("`x` must be a robmedfit object.")
  estimand <- match.arg(estimand)

  curves <- x$effects$curves
  d      <- curves[curves$estimand == estimand, , drop = FALSE]
  d      <- d[order(d$dose), ]

  d$se_approx  <- (d$upper - d$lower) / (2 * 1.96)
  d$frag_local <- d$se_approx / pmax(abs(d$estimate), 1e-8)

  n  <- nrow(d)
  d2 <- rep(NA_real_, n)
  if (n >= 3) {
    h <- diff(d$dose)
    for (i in 2:(n - 1)) {
      h1 <- h[i - 1]; h2 <- h[i]
      d2[i] <- 2 * (d$estimate[i + 1] / h2 -
                      d$estimate[i] * (1/h1 + 1/h2) +
                      d$estimate[i - 1] / h1) / (h1 + h2)
    }
    d2[1] <- d2[2]; d2[n] <- d2[n - 1]
  }
  d$curvature         <- d2
  d$in_fragility_zone <- d$lower < 0 & d$upper > 0

  d[, c("dose", "estimate", "lower", "upper", "se_approx",
        "frag_local", "curvature", "in_fragility_zone")]
}


# ==============================================================================
# plot_curvature
# ==============================================================================

#' Plot Dose-Varying Fragility (Curvature-Based Sensitivity)
#'
#' @description
#' Three-panel visualisation of dose-varying fragility: (1) effect curve with
#' CI bands and fragility zones, (2) local fragility index, (3) normalised
#' curvature.
#'
#' @param x        Data frame from [sensitivity_curvature()].
#' @param estimand Label for the estimand. Default `"NIE"`.
#' @param ref_dose Optional reference dose vertical line.
#' @param ...      Ignored.
#'
#' @return A `ggplot2` object.
#'
#' @examples
#' \donttest{
#' data(sim_mediation)
#'   fit <- robustmediate(
#'     X ~ Z1 + Z2, M ~ X + Z1 + Z2, Y ~ X + M + Z1 + Z2,
#'     data = sim_mediation, R = 20, verbose = FALSE
#'   )
#' curv <- sensitivity_curvature(fit)
#' plot_curvature(curv, ref_dose = fit$meta$ref_dose)
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_ribbon geom_line geom_rect geom_hline
#'   geom_vline labs theme_minimal theme facet_wrap element_text
#' @export
plot_curvature <- function(x, estimand = "NIE", ref_dose = NULL, ...) {
  if (!is.data.frame(x) || !"frag_local" %in% names(x)) {
    rlang::abort("`x` must be a data frame from sensitivity_curvature().")
  }

  x$curv_norm <- abs(x$curvature) / max(abs(x$curvature), na.rm = TRUE)
  x$frag_norm <- pmin(x$frag_local, 3) / 3

  panel_df <- rbind(
    data.frame(dose = x$dose, y = x$estimate,
               lower = x$lower,   upper = x$upper,
               fz    = x$in_fragility_zone,
               panel = paste0(estimand, " estimate"),
               stringsAsFactors = FALSE),
    data.frame(dose = x$dose, y = x$frag_norm,
               lower = NA_real_, upper = NA_real_,
               fz    = x$in_fragility_zone,
               panel = "Local fragility index",
               stringsAsFactors = FALSE),
    data.frame(dose = x$dose, y = x$curv_norm,
               lower = NA_real_, upper = NA_real_,
               fz    = x$in_fragility_zone,
               panel = "Curvature (normalised |d\u00b2/dx\u00b2|)",
               stringsAsFactors = FALSE)
  )
  panel_df$panel <- factor(
    panel_df$panel,
    levels = c(paste0(estimand, " estimate"),
               "Local fragility index",
               "Curvature (normalised |d\u00b2/dx\u00b2|)")
  )

  p <- ggplot2::ggplot(panel_df, ggplot2::aes(x = dose)) +
    ggplot2::geom_rect(
      data = panel_df[panel_df$fz, ],
      ggplot2::aes(xmin = dose - 0.02, xmax = dose + 0.02,
                   ymin = -Inf, ymax = Inf),
      fill = "#FDECEA", alpha = 0.5, inherit.aes = FALSE
    ) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed",
                        colour = "grey60", linewidth = 0.4) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper),
                         fill = "#0F6E56", alpha = 0.15, na.rm = TRUE) +
    ggplot2::geom_line(ggplot2::aes(y = y), colour = "#0F6E56",
                       linewidth = 0.9, na.rm = TRUE) +
    ggplot2::facet_wrap(~ panel, ncol = 1, scales = "free_y") +
    ggplot2::labs(
      x       = "Treatment dose",
      y       = NULL,
      title   = paste0("Dose-Varying Fragility: ", estimand),
      caption = paste0(
        "Red shading = fragility zone (CI crosses zero).  ",
        "Fragility index = SE / |estimate|; capped at 3\u00d7.  ",
        "Curvature normalised to [0, 1]."
      )
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      strip.text   = ggplot2::element_text(face = "bold", size = 9),
      plot.caption = ggplot2::element_text(size = 7.5, colour = "grey40")
    )

  if (!is.null(ref_dose)) {
    p <- p + ggplot2::geom_vline(
      xintercept = ref_dose,
      linetype   = "dotted",
      colour     = "grey50",
      linewidth  = 0.5
    )
  }
  p
}


# ==============================================================================
# fragility_table
# ==============================================================================

#' Pathway Fragility Decomposition Table
#'
#' @description
#' Returns a publication-ready table decomposing indirect-effect robustness
#' into pathway-specific components. Columns follow the medITCV reporting
#' convention: pathway, coefficient, SE, t, df, observed partial r, critical r,
#' medITCV, medITCV%, fragility classification, and tipping-point confounder r.
#'
#' @param x     A `robmedfit` object.
#' @param alpha Significance level. Default `0.05`.
#'
#' @return A data frame with three rows (a-path, b-path, indirect effect) and
#'   columns `pathway`, `coefficient`, `SE`, `t_stat`, `df`, `r_obs`,
#'   `r_crit`, `medITCV`, `medITCV_pct`, `fragility`,
#'   `tipping_r_confounder`, and `bottleneck`.
#'
#' @examples
#' \donttest{
#' data(sim_mediation)
#'   fit <- robustmediate(
#'     X ~ Z1 + Z2, M ~ X + Z1 + Z2, Y ~ X + M + Z1 + Z2,
#'     data = sim_mediation, R = 20, verbose = FALSE
#'   )
#' fragility_table(fit)
#' }
#'
#' @export
fragility_table <- function(x, alpha = 0.05) {
  if (!inherits(x, "robmedfit")) rlang::abort("`x` must be a robmedfit object.")

  mp <- sensitivity_meditcv_profile(x, alpha = alpha)

  classify <- function(val) {
    if (val >= 0.20)      "Highly robust"
    else if (val >= 0.10) "Robust"
    else if (val >= 0.05) "Moderately fragile"
    else if (val >= 0)    "Fragile"
    else                  "Not significant"
  }

  a   <- mp$meditcv_detail$a_path
  b   <- mp$meditcv_detail$b_path
  ind <- mp$meditcv_indirect

  data.frame(
    pathway              = c("a-path (X\u2192M)", "b-path (M\u2192Y|X)", "Indirect (a\u00d7b)"),
    coefficient          = c(round(a$coef,   4), round(b$coef,   4), NA_real_),
    SE                   = c(round(a$se,     4), round(b$se,     4), NA_real_),
    t_stat               = c(round(a$t_stat, 3), round(b$t_stat, 3), NA_real_),
    df                   = c(a$df,               b$df,               NA_real_),
    r_obs                = c(round(a$r_obs,  4), round(b$r_obs,  4), NA_real_),
    r_crit               = c(round(a$r_crit, 4), round(b$r_crit, 4), NA_real_),
    medITCV              = c(round(mp$a_path, 4), round(mp$b_path, 4), round(ind, 4)),
    medITCV_pct          = c(round(a$meditcv_pct, 1), round(b$meditcv_pct, 1), NA_real_),
    fragility            = c(classify(mp$a_path), classify(mp$b_path), classify(ind)),
    tipping_r_confounder = c(
      round(sqrt(max(mp$a_path, 0)), 4),
      round(sqrt(max(mp$b_path, 0)), 4),
      round(sqrt(max(ind, 0)), 4)
    ),
    bottleneck = c(
      mp$bottleneck == "a-path (X \u2192 M)",
      mp$bottleneck == "b-path (M \u2192 Y | X)",
      TRUE
    ),
    stringsAsFactors = FALSE
  )
}
