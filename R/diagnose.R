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
#' @return Invisibly returns a list with elements `balance`, `effects`, and
#'   `sensitivity`.
#'
#' @export
diagnose <- function(x, ...) {

  if (!inherits(x, "robmedfit")) {
    rlang::abort("`x` must be a robmedfit object.")
  }

  bal  <- x$balance$summary_stats
  eff  <- x$effects$summary
  sens <- x$sensitivity
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

  # ── Sensitivity ─────────────────────────────────────────
  cat("\n-- Sensitivity to mediator-outcome confounding --------------\n")

  if (is.null(sens)) {
    cat("  Not available: requires Gaussian identity-link mediator and outcome\n")
    cat("  models with the treatment and mediator entering as main effects.\n")
  } else {
    cat(sprintf("  Indirect effect is zero at rho = %.3f\n", sens$rho_zero))
    cat(sprintf("  Equivalent R2*_M x R2*_Y = %.4f\n", sens$r2_product_zero))

    cat("\n  Suggested text for Results section:\n")
    cat(sprintf(
      paste0(
        "  \"Following Imai, Keele, and Yamamoto (2010), the estimated indirect\n",
        "   effect would be zero if the correlation between the errors of the\n",
        "   mediator and outcome models were %.2f. Equivalently, an unmeasured\n",
        "   confounder would need to explain proportions of the residual variance\n",
        "   of the mediator and outcome whose product is %.3f.\"\n"
      ),
      sens$rho_zero, sens$r2_product_zero
    ))
    cat("  Use plot(fit, type = \"sensitivity\") for visualization.\n")
  }

  cat("\n--------------------------------------------------------------\n")

  invisible(list(
    balance         = bal,
    effects         = eff,
    sensitivity     = sens
  ))
}
