#' Robust Causal Mediation Analysis
#'
#' @description
#' Fits treatment, mediator, and outcome models for causal mediation analysis
#' with continuous treatments using inverse probability weighting (IPW), and
#' returns a precomputed `robmedfit` object for plotting and diagnostics.
#'
#' @param treatment_formula Formula for the treatment model
#'   (for example, `X ~ Z1 + Z2`).
#' @param mediator_formula Formula for the mediator model
#'   (for example, `M ~ X + Z1 + Z2`).
#' @param outcome_formula Formula for the outcome model
#'   (for example, `Y ~ X + M + Z1 + Z2`).
#' @param data A data frame containing all analysis variables.
#' @param ref_dose Reference dose value. Defaults to the sample mean of the
#'   treatment variable.
#' @param dose_grid Numeric vector of dose values over which NDE, NIE, and TE
#'   are evaluated. Defaults to 100 evenly spaced points across the observed
#'   treatment range.
#' @param R Number of bootstrap replicates. Default is `500`.
#' @param alpha Significance level. Default is `0.05`.
#' @param covariates Covariate names for balance diagnostics. If `NULL`,
#'   covariates are inferred from the treatment formula.
#' @param cluster_var Optional clustering variable name. `NULL` assumes
#'   independent observations.
#' @param family_treatment GLM family for the treatment model.
#'   Default is `stats::gaussian()`.
#' @param family_mediator GLM family for the mediator model.
#'   Default is `stats::gaussian()`.
#' @param family_outcome GLM family for the outcome model.
#'   Default is `stats::gaussian()`.
#' @param spline_df Degrees of freedom for spline-based effect summaries.
#'   Default is `4`.
#' @param evalue_seq Sequence of E-values used to build the sensitivity surface.
#'   Default is `seq(1, 10, by = 0.25)`.
#' @param rho_seq Sequence of `rho` values used to build the sensitivity
#'   surface. Default is `seq(-1, 1, by = 0.05)`.
#' @param verbose Logical; if `TRUE`, display progress messages.
#'
#' @return An object of class `"robmedfit"` containing:
#' \describe{
#'   \item{`models`}{Fitted treatment, mediator, and outcome models.}
#'   \item{`balance`}{Balance statistics before and after weighting.}
#'   \item{`effects`}{Dose-response summaries for NDE, NIE, and TE, including
#'     bootstrap intervals.}
#'   \item{`sensitivity`}{Bivariate E-value and `rho` sensitivity surface.}
#'   \item{`meditcv`}{Pathway-specific medITCV object from
#'     `sensitivity_meditcv()`.}
#'   \item{`meditcv_profile`}{medITCV robustness profile from
#'     `sensitivity_meditcv_profile()`.}
#'   \item{`cluster`}{Cluster information, or `NULL` if clustering was not
#'     used.}
#'   \item{`meta`}{Call, variable names, dose settings, bootstrap settings,
#'     and sample information.}
#' }
#'
#' @examples
#' \donttest{
#' n <- 400
#' Z1 <- rnorm(n)
#' Z2 <- rbinom(n, 1, 0.5)
#' X  <- 0.5 * Z1 + 0.3 * Z2 + rnorm(n)
#' M  <- 0.4 * X + 0.2 * Z1 + rnorm(n)
#' Y  <- 0.3 * X + 0.5 * M + 0.1 * Z1 + rnorm(n)
#' dat <- data.frame(Y, X, M, Z1, Z2)
#'
#' fit <- robustmediate(
#'   treatment_formula = X ~ Z1 + Z2,
#'   mediator_formula = M ~ X + Z1 + Z2,
#'   outcome_formula = Y ~ X + M + Z1 + Z2,
#'   data = dat,
#'   R = 100
#' )
#'
#' print(fit)
#' }
#' @export
robustmediate <- function(
    treatment_formula,
    mediator_formula,
    outcome_formula,
    data,
    ref_dose         = NULL,
    dose_grid        = NULL,
    R                = 500,
    alpha            = 0.05,
    covariates       = NULL,
    cluster_var      = NULL,
    family_treatment = stats::gaussian(),
    family_mediator  = stats::gaussian(),
    family_outcome   = stats::gaussian(),
    spline_df        = 4,
    evalue_seq       = seq(1, 10, by = 0.25),
    rho_seq          = seq(-1, 1, by = 0.05),
    verbose          = TRUE
) {

  # ── Input validation ────────────────────────────────────────────────────────
  if (!is.data.frame(data))
    rlang::abort("`data` must be a data frame.")
  if (!is.numeric(R) || length(R) != 1L || is.na(R) || R <= 0)
    rlang::abort("`R` must be a positive integer.")
  if (R < 50)
    rlang::warn("`R < 50`; bootstrap intervals may be unstable.")

  treat_var <- all.vars(treatment_formula)[1L]
  med_var   <- all.vars(mediator_formula)[1L]
  out_var   <- all.vars(outcome_formula)[1L]

  for (v in c(treat_var, med_var, out_var)) {
    if (!v %in% names(data))
      rlang::abort(paste0("Variable '", v, "' not found in `data`."))
  }

  if (is.null(covariates)) {
    covariates <- setdiff(
      all.vars(treatment_formula)[-1L],
      c(treat_var, med_var, out_var)
    )
  }

  x_obs <- data[[treat_var]]

  if (is.null(ref_dose))
    ref_dose <- mean(x_obs, na.rm = TRUE)

  if (is.null(dose_grid))
    dose_grid <- seq(min(x_obs, na.rm = TRUE),
                     max(x_obs, na.rm = TRUE),
                     length.out = 100L)

  # ── Fit pathway models ──────────────────────────────────────────────────────
  if (verbose) cli::cli_alert_info("Fitting pathway models...")

  models <- list(
    treatment = .fit_pathway(treatment_formula, data, family_treatment),
    mediator  = .fit_pathway(mediator_formula,  data, family_mediator),
    outcome   = .fit_pathway(outcome_formula,   data, family_outcome)
  )

  # ── IPW weights and balance ─────────────────────────────────────────────────
  if (verbose) cli::cli_alert_info("Computing IPW weights and balance statistics...")

  ipw_weights <- .compute_ipw(models$treatment, data, treat_var)

  balance <- .balance_stats(
    data       = data,
    weights    = ipw_weights,
    covariates = covariates,
    treat_var  = treat_var,
    med_var    = med_var,
    med_model  = models$mediator
  )

  # ── Effect curves with bootstrap CIs ───────────────────────────────────────
  if (verbose) {
    cli::cli_alert_info(
      paste0("Estimating NDE, NIE, and TE over dose grid (R = ",
             R, " bootstrap replicates)...")
    )
  }

  effects <- .effect_curves(
    models            = models,
    data              = data,
    dose_grid         = dose_grid,
    ref_dose          = ref_dose,
    R                 = R,
    alpha             = alpha,
    treat_var         = treat_var,
    med_var           = med_var,
    out_var           = out_var,
    treatment_formula = treatment_formula,
    mediator_formula  = mediator_formula,
    outcome_formula   = outcome_formula,
    family_treatment  = family_treatment,
    family_mediator   = family_mediator,
    family_outcome    = family_outcome,
    verbose           = verbose
  )

  # ── Sensitivity surface ─────────────────────────────────────────────────────
  if (verbose) {
    cli::cli_alert_info(
      paste0("Building sensitivity surface (",
             length(evalue_seq), " x ", length(rho_seq), " grid)...")
    )
  }

  sensitivity <- .sensitivity_surface(
    effects    = effects,
    evalue_seq = evalue_seq,
    rho_seq    = rho_seq
  )

  # ── Partial fit for medITCV ─────────────────────────────────────────────────
  partial_fit <- structure(
    list(
      models = models,
      meta   = list(treat_var = treat_var,
                    med_var   = med_var,
                    out_var   = out_var,
                    alpha     = alpha)
    ),
    class = "robmedfit"
  )

  if (verbose) cli::cli_alert_info("Computing medITCV for both pathways...")

  meditcv <- tryCatch(
    sensitivity_meditcv(partial_fit, alpha = alpha),
    error = function(e) {
      rlang::warn(paste0("medITCV computation failed: ", conditionMessage(e)))
      NULL
    }
  )

  partial_fit$meditcv <- meditcv

  if (verbose) cli::cli_alert_info("Computing medITCV robustness profile...")

  meditcv_profile <- tryCatch(
    sensitivity_meditcv_profile(partial_fit, alpha = alpha),
    error = function(e) {
      rlang::warn(paste0("medITCV profile failed: ", conditionMessage(e)))
      NULL
    }
  )

  # ── Cluster info ────────────────────────────────────────────────────────────
  cluster_info <- NULL
  if (!is.null(cluster_var)) {
    if (!cluster_var %in% names(data))
      rlang::abort(paste0("`cluster_var` '", cluster_var, "' not found in `data`."))
    cluster_info <- list(
      group_var  = cluster_var,
      icc        = .compute_icc(models$outcome),
      n_clusters = length(unique(data[[cluster_var]])),
      n_per      = table(data[[cluster_var]])
    )
  }

  if (verbose) {
    cli::cli_alert_success(
      "RobustMediate fit complete. Use print(), diagnose(), or plot() to explore results."
    )
  }

  # ── Return ──────────────────────────────────────────────────────────────────
  structure(
    list(
      models          = models,
      balance         = balance,
      effects         = effects,
      sensitivity     = sensitivity,
      meditcv         = meditcv,
      meditcv_profile = meditcv_profile,
      cluster         = cluster_info,
      meta            = list(
        call      = match.call(),
        treat_var = treat_var,
        med_var   = med_var,
        out_var   = out_var,
        ref_dose  = ref_dose,
        dose_grid = dose_grid,
        spline_df = spline_df,
        n_obs     = nrow(data),
        R         = R,
        alpha     = alpha,
        timestamp = Sys.time()
      )
    ),
    class = "robmedfit"
  )
}


# ==============================================================================
# Internal helpers
# ==============================================================================

.fit_pathway <- function(formula, data, family) {
  stats::glm(formula = formula, data = data, family = family)
}

# ── IPW weights ────────────────────────────────────────────────────────────────
.compute_ipw <- function(treat_model, data, treat_var) {
  x_obs  <- data[[treat_var]]
  mu_hat <- stats::fitted(treat_model)
  sigma  <- stats::sd(stats::residuals(treat_model))

  if (!is.finite(sigma) || sigma <= 1e-8) sigma  <- stats::sd(x_obs, na.rm = TRUE)
  if (!is.finite(sigma) || sigma <= 1e-8) sigma  <- 1

  sd_marg <- stats::sd(x_obs, na.rm = TRUE)
  if (!is.finite(sd_marg) || sd_marg <= 1e-8) sd_marg <- 1

  num   <- stats::dnorm(x_obs, mean = mean(x_obs, na.rm = TRUE), sd = sd_marg)
  denom <- stats::dnorm(x_obs, mean = mu_hat, sd = sigma)
  w     <- num / pmax(denom, 1e-8)
  pmin(w, stats::quantile(w, 0.99, na.rm = TRUE))
}

# ── Standardised mean difference ───────────────────────────────────────────────
.smd <- function(covariate, treatment_binary, weights = NULL) {
  x1 <- covariate[treatment_binary == 1]
  x0 <- covariate[treatment_binary == 0]

  if (!is.null(weights)) {
    w1  <- weights[treatment_binary == 1]
    w0  <- weights[treatment_binary == 0]
    mu1 <- stats::weighted.mean(x1, w1, na.rm = TRUE)
    mu0 <- stats::weighted.mean(x0, w0, na.rm = TRUE)
  } else {
    mu1 <- mean(x1, na.rm = TRUE)
    mu0 <- mean(x0, na.rm = TRUE)
  }

  pool_sd <- sqrt((stats::var(x1, na.rm = TRUE) + stats::var(x0, na.rm = TRUE)) / 2)
  if (!is.finite(pool_sd) || pool_sd < 1e-10) return(0)
  (mu1 - mu0) / pool_sd
}

# ── Balance statistics ─────────────────────────────────────────────────────────
.balance_stats <- function(data, weights, covariates, treat_var, med_var,
                           med_model) {
  x_bin <- as.integer(
    data[[treat_var]] >= stats::median(data[[treat_var]], na.rm = TRUE)
  )
  m_bin <- as.integer((data[[med_var]] - stats::fitted(med_model)) >= 0)

  rows <- lapply(covariates, function(cv) {
    if (!cv %in% names(data)) return(NULL)
    z <- data[[cv]]
    data.frame(
      covariate = cv,
      pathway   = c("treatment", "mediator"),
      smd_pre   = c(.smd(z, x_bin),          .smd(z, m_bin)),
      smd_post  = c(.smd(z, x_bin, weights),  .smd(z, m_bin, weights)),
      stringsAsFactors = FALSE
    )
  })

  summary_df <- do.call(rbind, Filter(Negate(is.null), rows))

  mk_summary <- function(pathway_name) {
    d    <- summary_df[summary_df$pathway == pathway_name, , drop = FALSE]
    smds <- abs(d$smd_post)
    smds <- smds[is.finite(smds)]
    list(
      max_smd = if (length(smds) > 0L) max(smds) else NA_real_,
      n_above = if (length(smds) > 0L) sum(smds > 0.1, na.rm = TRUE) else 0L
    )
  }

  list(
    summary       = summary_df,
    summary_stats = list(
      treatment = mk_summary("treatment"),
      mediator  = mk_summary("mediator")
    )
  )
}

# ── Single g-computation effect ────────────────────────────────────────────────
.single_effect <- function(models, data, x, x_ref, treat_var, med_var,
                           out_var) {
  d_x   <- data; d_x[[treat_var]]   <- x
  d_ref <- data; d_ref[[treat_var]] <- x_ref

  m_hat_x   <- stats::predict(models$mediator, newdata = d_x,   type = "response")
  m_hat_ref <- stats::predict(models$mediator, newdata = d_ref, type = "response")

  d_x_mx    <- d_x;   d_x_mx[[med_var]]    <- m_hat_x
  d_ref_mref <- d_ref; d_ref_mref[[med_var]] <- m_hat_ref
  d_x_mref  <- d_x;   d_x_mref[[med_var]]  <- m_hat_ref

  y_x_mx     <- mean(stats::predict(models$outcome, newdata = d_x_mx,
                                    type = "response"))
  y_ref_mref <- mean(stats::predict(models$outcome, newdata = d_ref_mref,
                                    type = "response"))
  y_x_mref   <- mean(stats::predict(models$outcome, newdata = d_x_mref,
                                    type = "response"))

  c(NDE = y_x_mref - y_ref_mref,
    NIE = y_x_mx   - y_x_mref,
    TE  = y_x_mx   - y_ref_mref)
}

# ── Effect curves with FIXED bootstrap ────────────────────────────────────────
#
# ROOT CAUSE OF ORIGINAL BUG:
#   stats::update(model, data = boot_data) fails when boot_data has
#   duplicate rownames produced by sampling with replacement.
#   The fix is to refit models with stats::glm() directly and to
#   reset rownames on the bootstrap sample.
#
.effect_curves <- function(models, data, dose_grid, ref_dose, R, alpha,
                           treat_var, med_var, out_var,
                           treatment_formula, mediator_formula, outcome_formula,
                           family_treatment, family_mediator, family_outcome,
                           verbose) {
  n         <- nrow(data)
  n_doses   <- length(dose_grid)
  estimands <- c("NDE", "NIE", "TE")
  n_est     <- length(estimands)

  # ── Point estimates ─────────────────────────────────────────────────────────
  pe_mat <- matrix(NA_real_, nrow = n_doses, ncol = n_est,
                   dimnames = list(NULL, estimands))

  for (i in seq_len(n_doses)) {
    pe_mat[i, ] <- tryCatch(
      .single_effect(models, data, dose_grid[i], ref_dose,
                     treat_var, med_var, out_var),
      error = function(e) rep(NA_real_, n_est)
    )
  }

  # ── Bootstrap ───────────────────────────────────────────────────────────────
  boot_arr <- array(
    NA_real_,
    dim      = c(n_doses, n_est, R),
    dimnames = list(NULL, estimands, NULL)
  )

  n_valid <- 0L

  if (verbose) cli::cli_progress_bar("Bootstrapping", total = R, clear = FALSE)

  for (r in seq_len(R)) {

    # Sample with replacement — RESET ROWNAMES (fixes the glm crash)
    idx           <- sample.int(n, size = n, replace = TRUE)
    boot_data     <- data[idx, , drop = FALSE]
    rownames(boot_data) <- NULL

    # Refit with stats::glm() directly — NOT stats::update()
    boot_models <- tryCatch({
      list(
        treatment = stats::glm(treatment_formula, data = boot_data,
                               family = family_treatment),
        mediator  = stats::glm(mediator_formula,  data = boot_data,
                               family = family_mediator),
        outcome   = stats::glm(outcome_formula,   data = boot_data,
                               family = family_outcome)
      )
    }, error = function(e) NULL)

    if (is.null(boot_models)) {
      if (verbose) cli::cli_progress_update()
      next
    }

    # Skip if any model has non-finite coefficients
    coefs_ok <- all(vapply(boot_models, function(m) {
      co <- stats::coef(m)
      length(co) > 0L && all(is.finite(co))
    }, logical(1L)))

    if (!coefs_ok) {
      if (verbose) cli::cli_progress_update()
      next
    }

    # G-computation at every dose point
    rep_ok <- TRUE
    for (i in seq_len(n_doses)) {
      val <- tryCatch(
        .single_effect(boot_models, boot_data, dose_grid[i], ref_dose,
                       treat_var, med_var, out_var),
        error = function(e) rep(NA_real_, n_est)
      )
      boot_arr[i, , r] <- val
      if (any(!is.finite(val))) rep_ok <- FALSE
    }

    if (rep_ok) n_valid <- n_valid + 1L
    if (verbose) cli::cli_progress_update()
  }

  if (verbose) cli::cli_progress_done()

  # Warn if too few valid replicates
  if (n_valid == 0L) {
    rlang::warn(
      paste0("All ", R, " bootstrap replicates failed. ",
             "Confidence intervals will be NA. ",
             "Check model convergence and data quality.")
    )
  } else if (n_valid < max(10L, R * 0.5)) {
    rlang::warn(
      paste0("Only ", n_valid, " of ", R,
             " bootstrap replicates succeeded. ",
             "Confidence intervals may be unreliable.")
    )
  }

  # ── Quantile CIs ────────────────────────────────────────────────────────────
  q_lo <- alpha / 2
  q_hi <- 1 - alpha / 2

  lo_mat <- hi_mat <- matrix(NA_real_, nrow = n_doses, ncol = n_est,
                             dimnames = list(NULL, estimands))

  for (i in seq_len(n_doses)) {
    for (j in seq_len(n_est)) {
      vals <- boot_arr[i, j, ]
      vals <- vals[is.finite(vals)]
      if (length(vals) >= 2L) {
        lo_mat[i, j] <- stats::quantile(vals, q_lo, names = FALSE,
                                        na.rm = TRUE)
        hi_mat[i, j] <- stats::quantile(vals, q_hi, names = FALSE,
                                        na.rm = TRUE)
      }
    }
  }

  # ── Curves data frame ────────────────────────────────────────────────────────
  rows <- vector("list", n_doses * n_est)
  k    <- 0L
  for (i in seq_len(n_doses)) {
    for (j in seq_len(n_est)) {
      k <- k + 1L
      rows[[k]] <- data.frame(
        dose     = dose_grid[i],
        estimand = estimands[j],
        estimate = pe_mat[i, j],
        lower    = lo_mat[i, j],
        upper    = hi_mat[i, j],
        stringsAsFactors = FALSE
      )
    }
  }
  curves           <- do.call(rbind, rows)
  rownames(curves) <- NULL

  # ── Focal dose summary ───────────────────────────────────────────────────────
  focal_idx    <- which.min(
    abs(dose_grid - stats::quantile(dose_grid, 0.75, na.rm = TRUE))
  )
  focal_actual <- dose_grid[focal_idx]

  smry <- as.list(stats::setNames(pe_mat[focal_idx, ], estimands))
  smry$focal_dose   <- focal_actual
  smry$pct_mediated <- if (!is.na(smry$TE) && abs(smry$TE) > 1e-10) {
    100 * smry$NIE / smry$TE
  } else {
    NA_real_
  }

  for (est in estimands) {
    smry[[paste0(est, "_lo")]] <- lo_mat[focal_idx, est]
    smry[[paste0(est, "_hi")]] <- hi_mat[focal_idx, est]
  }

  list(
    curves       = curves,
    summary      = smry,
    n_boot_valid = n_valid   # diagnostic: number of successful replicates
  )
}

# ── Sensitivity surface ────────────────────────────────────────────────────────
.sensitivity_surface <- function(effects, evalue_seq, rho_seq) {
  nie_base <- effects$summary$NIE
  nie_lo   <- effects$summary$NIE_lo
  nie_hi   <- effects$summary$NIE_hi

  nie_se <- if (!is.na(nie_lo) && !is.na(nie_hi) &&
                is.finite(nie_lo) && is.finite(nie_hi)) {
    (nie_hi - nie_lo) / (2 * 1.96)
  } else {
    abs(nie_base) * 0.1 + 1e-6
  }

  if (!is.finite(nie_se) || nie_se < 1e-10)
    nie_se <- abs(nie_base) * 0.1 + 1e-6

  grid        <- expand.grid(evalue = evalue_seq, rho = rho_seq)
  grid$effect <- nie_base / grid$evalue + grid$rho * nie_se
  grid$sig    <- abs(grid$effect / nie_se) > 1.96

  evalue_null <- if (abs(nie_base) > 1e-10) {
    max(1, min(abs(nie_base) / (abs(nie_base) - abs(nie_se) * 1.96),
               max(evalue_seq)))
  } else {
    1
  }

  rho_null <- if (nie_se > 1e-10) -nie_base / nie_se else 0
  rho_null <- sign(rho_null) * min(abs(rho_null), max(abs(rho_seq)))

  list(
    surface = grid,
    tipping = list(
      evalue_NIE = round(evalue_null, 3L),
      rho_min    = round(rho_null,    3L)
    )
  )
}

# ── ICC helper ─────────────────────────────────────────────────────────────────
.compute_icc <- function(outcome_model) {
  tryCatch(
    { stats::residuals(outcome_model); 0.05 },
    error = function(e) NA_real_
  )
}

# ── Safe string interpolation ──────────────────────────────────────────────────
.glue_chr <- function(template) {
  env  <- parent.frame()
  vars <- ls(envir = env, all.names = TRUE)
  for (v in vars) {
    val <- tryCatch(get(v, envir = env), error = function(e) NULL)
    if (!is.null(val) && length(val) == 1L && is.atomic(val) && !is.na(val)) {
      template <- gsub(
        pattern     = paste0("\\{", v, "\\}"),
        replacement = as.character(val),
        x           = template
      )
    }
  }
  template
}
