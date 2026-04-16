#' Diagnose a robmedfit Object
#'
#' @description
#' Prints a formatted diagnostics report covering balance, mediation effects,
#' and sensitivity robustness. The output is structured so that it can be used
#' directly (or with minimal editing) in the Results section of an applied paper.
#' Returns the underlying results invisibly.
#'
#' @param x A `robmedfit` object.
#' @param ... Ignored.
#'
#' @return Invisibly returns a list with elements `balance`, `effects`,
#'   `sensitivity`, `meditcv`, and `meditcv_profile`.
#'
#' @export
diagnose <- function(x, ...) {

  if (!inherits(x, "robmedfit")) {
    rlang::abort("`x` must be a robmedfit object.")
  }

  bal  <- x$balance$summary_stats
  eff  <- x$effects$summary
  sens <- x$sensitivity$tipping
  m    <- x$meta

  # ── Header ─────────────────────────────────────────────
  cat("==============================================================\n")
  cat("           RobustMediate: Diagnostics Report\n")
  cat("==============================================================\n\n")

  # ── Balance ────────────────────────────────────────────
  cat("-- Balance diagnostics ---------------------------------------\n")
  cat(sprintf(
    "  Treatment pathway: max |SMD| = %.3f  (%d covariate(s) above 0.10)\n",
    bal$treatment$max_smd, bal$treatment$n_above
  ))
  cat(sprintf(
    "  Mediator pathway:  max |SMD| = %.3f  (%d covariate(s) above 0.10)\n",
    bal$mediator$max_smd, bal$mediator$n_above
  ))

  if (max(bal$treatment$max_smd, bal$mediator$max_smd) <= 0.10) {
    cat("  OK: Balance is satisfactory on both pathways (all |SMD| <= 0.10).\n")
  } else {
    cat("  WARNING: Balance exceeds 0.10 on at least one covariate.\n")
    cat("  Inspect plot_balance() carefully.\n")
  }

  # ── Mediation effects ───────────────────────────────────
  focal <- eff$focal_dose

  cat("\n-- Mediation effects (focal dose =", round(focal, 3),
      "vs ref =", round(m$ref_dose, 3), ") ----------------\n")

  for (e in c("NDE", "NIE", "TE")) {
    cat(sprintf(
      "  %-4s  %+.4f  [%+.4f, %+.4f]  (%d%% CI)\n",
      e,
      eff[[e]],
      eff[[paste0(e, "_lo")]],
      eff[[paste0(e, "_hi")]],
      round(100 * (1 - m$alpha))
    ))
  }

  cat(sprintf("  %% mediated: %.1f%%\n", eff$pct_mediated))

  # ── Sensitivity 1: E-value / rho ────────────────────────
  cat("\n-- Sensitivity 1: E-value x rho surface ----------------------\n")
  cat(sprintf("  E-value to nullify NIE:       %.2f\n", sens$evalue_NIE))
  cat(sprintf("  Minimum |rho| to nullify NIE: %.2f\n", abs(sens$rho_min)))

  cat("\n  Suggested text for Results section:\n")
  cat(sprintf(
    paste0(
      "  \"Sensitivity analyses indicated that an unmeasured confounder\n",
      "   would need to be associated %.1f times with both treatment and outcome\n",
      "   to fully explain away the estimated indirect effect\n",
      "   (NIE = %.3f, %d%% CI [%.3f, %.3f]).\n",
      "   The sequential-ignorability violation |rho| would need to exceed\n",
      "   %.2f to nullify this finding.\"\n"
    ),
    sens$evalue_NIE,
    eff$NIE,
    round(100 * (1 - m$alpha)),
    eff$NIE_lo,
    eff$NIE_hi,
    abs(sens$rho_min)
  ))

  # ── Sensitivity 2: medITCV ──────────────────────────────
  cat("\n-- Sensitivity 2: medITCV ------------------------------------\n")

  med_obj <- x$meditcv

  if (is.null(med_obj)) {

    cat("  medITCV not available (computation failed; check warnings).\n")

  } else {

    a   <- med_obj$a_path
    b   <- med_obj$b_path
    ind <- med_obj$indirect

    classify <- function(v) {
      if (v < 0) "not significant"
      else if (v < 0.10) "fragile"
      else "robust"
    }

    cat("  a-path (X -> M)\n")
    cat(sprintf(
      "    Observed r = %.4f  |  Critical r = %.4f  |  medITCV = %.4f  [%s]\n",
      a$r_obs, a$r_crit, a$meditcv, classify(a$meditcv)
    ))
    cat(sprintf(
      "    Need |r_XC * r_YC| > %.4f (%.1f%% above threshold)\n",
      a$meditcv, a$meditcv_pct
    ))

    cat("  b-path (M -> Y | X)\n")
    cat(sprintf(
      "    Observed r = %.4f  |  Critical r = %.4f  |  medITCV = %.4f  [%s]\n",
      b$r_obs, b$r_crit, b$meditcv, classify(b$meditcv)
    ))
    cat(sprintf(
      "    Need |r_XC * r_YC| > %.4f (%.1f%% above threshold)\n",
      b$meditcv, b$meditcv_pct
    ))

    cat(sprintf(
      "\n  Bottleneck: %s (minimum-path medITCV = %.4f)\n",
      ind$bottleneck,
      ind$meditcv_indirect
    ))

    cat("  Use plot(fit, type = \"meditcv\") for visualization.\n")
  }

  # ── Sensitivity 3: medITCV profile ──────────────────────
  cat("\n-- Sensitivity 3: medITCV robustness profile -----------------\n")

  if (!is.null(x$meditcv_profile)) {

    mp <- x$meditcv_profile

    cat(sprintf("  medITCV_a (X->M):        %.4f\n", mp$a_path))
    cat(sprintf("  medITCV_b (M->Y|X):      %.4f\n", mp$b_path))
    cat(sprintf("  medITCV_indirect:        %.4f\n", mp$meditcv_indirect))
    cat(sprintf("  Fragility ratio (a/b):   %.3f\n", mp$fragility_ratio))
    cat(sprintf("  Bottleneck:              %s\n", mp$bottleneck))

    cat(sprintf(
      "\n  Tipping-point confounder r = %.3f would nullify the indirect effect.\n",
      sqrt(max(mp$meditcv_indirect, 0))
    ))

    cat("  Use plot(fit, type = \"meditcv_profile\") for visualization.\n")

  } else {

    cat("  medITCV robustness profile not available.\n")
  }

  cat("\n--------------------------------------------------------------\n")

  invisible(list(
    balance         = bal,
    effects         = eff,
    sensitivity     = sens,
    meditcv         = x$meditcv,
    meditcv_profile = x$meditcv_profile
  ))
}
