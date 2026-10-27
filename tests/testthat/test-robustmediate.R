library(testthat)
library(RobustMediate)

# ── Shared fixture (runs once, reused across tests) ───────────────────────────
set.seed(42)
n  <- 200
Z1 <- rnorm(n); Z2 <- rbinom(n, 1, 0.5)
X  <- 0.5 * Z1 + 0.3 * Z2 + rnorm(n)
M  <- 0.4 * X  + 0.2 * Z1 + rnorm(n)
Y  <- 0.3 * X  + 0.5 * M  + 0.1 * Z1 + rnorm(n)
dat <- data.frame(Y, X, M, Z1, Z2)

fit <- suppressWarnings(robustmediate(
  treatment_formula = X ~ Z1 + Z2,
  mediator_formula  = M ~ X + Z1 + Z2,
  outcome_formula   = Y ~ X + M + Z1 + Z2,
  data    = dat,
  R       = 20,       # keep fast for testing
  verbose = FALSE
))


# ── 1. robustmediate() basic structure ────────────────────────────────────────
test_that("robustmediate() returns a robmedfit object with correct slots", {
  expect_s3_class(fit, "robmedfit")
  expect_true(all(c("models", "balance", "effects", "sensitivity",
                    "cluster", "meta") %in% names(fit)))
  expect_false(any(c("meditcv", "meditcv_profile") %in% names(fit)))
})

test_that("meta slot contains expected fields", {
  expect_named(fit$meta,
               c("call", "treat_var", "med_var", "out_var", "ref_dose",
                 "dose_grid", "n_obs", "R", "alpha", "timestamp"),
               ignore.order = TRUE)
  expect_equal(fit$meta$treat_var, "X")
  expect_equal(fit$meta$med_var,   "M")
  expect_equal(fit$meta$out_var,   "Y")
  expect_equal(fit$meta$n_obs,     n)
})

test_that("effects curves data frame has correct shape", {
  crv <- fit$effects$curves
  expect_s3_class(crv, "data.frame")
  expect_true(all(c("dose", "estimand", "estimate", "lower", "upper") %in% names(crv)))
  expect_setequal(unique(crv$estimand), c("NDE", "NIE", "TE"))
})

test_that("effects summary contains expected fields", {
  sm <- fit$effects$summary
  expect_true(all(c("NDE", "NIE", "TE", "NDE_lo", "NIE_lo", "TE_lo",
                    "NDE_hi", "NIE_hi", "TE_hi", "pct_mediated") %in% names(sm)))
})


# ── 2. Balance ────────────────────────────────────────────────────────────────
test_that("balance statistics are computed for both pathways", {
  bal <- fit$balance
  expect_true(all(c("summary", "summary_stats") %in% names(bal)))
  expect_true(all(c("treatment", "mediator") %in% names(bal$summary_stats)))
  # max_smd is a non-negative numeric
  expect_gte(fit$balance$summary_stats$treatment$max_smd, 0)
  expect_gte(fit$balance$summary_stats$mediator$max_smd,  0)
})


# ── 3. Sensitivity to mediator-outcome confounding ───────────────────────────
test_that("sensitivity_rho() returns a robmed_sens object", {
  sens <- fit$sensitivity
  expect_s3_class(sens, "robmed_sens")
  expect_true(all(c("rho", "nie_unit", "nie") %in% names(sens$curve)))
  expect_equal(sens$r2_product_zero, sens$rho_zero^2)
})

test_that("indirect effect at rho = 0 equals a * b from the fitted models", {
  ab <- unname(coef(fit$models$mediator)["X"] * coef(fit$models$outcome)["M"])
  expect_equal(fit$sensitivity$nie_unit_0, ab, tolerance = 1e-10)
  s0 <- sensitivity_rho(fit, rho_seq = 0)
  expect_equal(s0$curve$nie_unit, ab, tolerance = 1e-10)
})

test_that("indirect effect at rho = 0 matches g-computation for linear models", {
  contrast <- fit$effects$summary$focal_dose - fit$meta$ref_dose
  expect_equal(fit$sensitivity$nie_unit_0 * contrast,
               fit$effects$summary$NIE, tolerance = 1e-8)
})

test_that("indirect effect is zero at rho_zero", {
  s1 <- sensitivity_rho(fit, rho_seq = fit$sensitivity$rho_zero)
  expect_equal(s1$curve$nie_unit, 0, tolerance = 1e-10)
})

test_that("rho recovers the indirect effect under known error correlation", {
  set.seed(7)
  nn <- 20000
  z  <- rnorm(nn)
  x  <- 0.5 * z + rnorm(nn)
  e  <- MASS_free_bivnorm(nn, 0.4)
  m  <- 0.6 * x + 0.3 * z + e[, 1]
  y  <- 0.2 * x + 0.5 * m + 0.2 * z + e[, 2]
  d  <- data.frame(x, m, y, z)
  f  <- suppressWarnings(robustmediate(x ~ z, m ~ x + z, y ~ x + m + z,
                                       data = d, R = 2, verbose = FALSE))
  # naive estimate is biased; rho = 0.4 recovers the true a * b = 0.3
  expect_gt(abs(f$sensitivity$nie_unit_0 - 0.3), 0.1)
  s <- sensitivity_rho(f, rho_seq = 0.4)
  expect_equal(s$curve$nie_unit, 0.3, tolerance = 0.02)
})

test_that("sensitivity is skipped for non-linear models", {
  d2 <- dat
  d2$Yb <- as.integer(d2$Y > median(d2$Y))
  f2 <- suppressWarnings(robustmediate(
    X ~ Z1 + Z2, M ~ X + Z1 + Z2, Yb ~ X + M + Z1 + Z2,
    data = d2, family_outcome = binomial(), R = 5, verbose = FALSE
  ))
  expect_null(f2$sensitivity)
  expect_error(sensitivity_rho(f2), "Gaussian")
  f3 <- suppressWarnings(robustmediate(
    X ~ Z1 + Z2, M ~ X + Z1 + Z2, Y ~ X * M + Z1 + Z2,
    data = dat, R = 5, verbose = FALSE
  ))
  expect_null(f3$sensitivity)
})

test_that("removed ITCV functions are no longer exported", {
  ns <- getNamespaceExports("RobustMediate")
  expect_false(any(c("sensitivity_meditcv", "sensitivity_meditcv_profile",
                     "plot_meditcv", "plot_meditcv_profile",
                     "fragility_table", "tipping_table") %in% ns))
})


# ── 6. Bootstrap CIs ─────────────────────────────────────────────────────────
test_that("bootstrap CIs are non-degenerate (lower != upper) with R >= 20", {
  eff <- fit$effects$summary
  # With R=20 reps and real data, lo != hi
  expect_false(isTRUE(all.equal(eff$NDE_lo, eff$NDE_hi)))
  expect_false(isTRUE(all.equal(eff$NIE_lo, eff$NIE_hi)))
  expect_false(isTRUE(all.equal(eff$TE_lo,  eff$TE_hi)))
})


# ── 9. plot_curvature ────────────────────────────────────────────────────────
test_that("plot_curvature() returns a ggplot object", {
  curv <- suppressMessages(effect_curvature(fit, estimand = "NIE"))
  p2   <- plot_curvature(curv, ref_dose = fit$meta$ref_dose)
  expect_s3_class(p2, "gg")
})

# ── 10. diagnose ─────────────────────────────────────────────────────────────
test_that("diagnose() returns list invisibly and prints text", {
  result <- capture.output(out <- diagnose(fit))
  expect_true(length(result) > 0)
  expect_named(out, c("balance", "effects", "sensitivity"), ignore.order = TRUE)
  expect_false(any(grepl("ITCV", result)))
})


# ── 11. Print methods ─────────────────────────────────────────────────────────
test_that("print methods produce output", {
  expect_true(length(capture.output(print(fit))) > 0)
  expect_true(any(grepl("rho", capture.output(print(fit$sensitivity)))))
})


# ── 12. broom methods ─────────────────────────────────────────────────────────
test_that("tidy.robmedfit() returns a data frame with NDE/NIE/TE rows", {
  td <- tidy(fit)
  expect_s3_class(td, "data.frame")
  expect_setequal(td$term, c("NDE", "NIE", "TE"))
  expect_true(all(c("estimate", "conf.low", "conf.high") %in% names(td)))
})

test_that("glance.robmedfit() returns a one-row data frame", {
  gl <- glance(fit)
  expect_s3_class(gl, "data.frame")
  expect_equal(nrow(gl), 1)
  expect_true(all(c("n_obs", "R", "pct_mediated", "rho_zero_NIE") %in% names(gl)))
})


# ── 13. Mediation plots ───────────────────────────────────────────────────────
test_that("plot_mediation() returns a ggplot", {
  p <- plot_mediation(fit)
  expect_s3_class(p, "gg")
})

test_that("plot_balance() returns a ggplot", {
  p <- plot_balance(fit)
  expect_s3_class(p, "gg")
})

test_that("plot_sensitivity() returns a ggplot", {
  p <- plot_sensitivity(fit)
  expect_s3_class(p, "gg")
})

test_that("plot(fit) defaults to mediation plot", {
  p <- plot(fit)
  expect_s3_class(p, "gg")
})

test_that("plot(fit, type='balance') returns a ggplot", {
  p <- plot(fit, type = "balance")
  expect_s3_class(p, "gg")
})

test_that("plot(fit, type='curvature') returns a ggplot", {
  p <- suppressMessages(plot(fit, type = "curvature"))
  expect_s3_class(p, "gg")
})


# ── 14. sensitivity_paths() ──────────────────────────────────────────────────
test_that("adjusted path t-statistics match a regression that includes U", {
  set.seed(11)
  nn <- 3000
  c1 <- rnorm(nn); u <- rnorm(nn)
  x  <- 0.4 * c1 + 0.5 * u + rnorm(nn)
  m  <- 0.3 * x + 0.3 * c1 + 0.5 * u + rnorm(nn)
  y  <- 0.2 * x + 0.3 * m + 0.2 * c1 + 0.4 * u + rnorm(nn)
  d  <- data.frame(x, m, y, c1, u)

  sp <- sensitivity_paths(mediator_model = lm(m ~ x + c1, d),
                          outcome_model  = lm(y ~ x + m + c1, d),
                          treatment = "x", mediator = "m")

  pr2 <- function(fit, v) { t <- summary(fit)$coef[v, "t value"]; t^2 / (t^2 + fit$df.residual) }
  r2_xu <- pr2(lm(x ~ c1 + u, d), "u")
  r2_mu <- pr2(lm(m ~ x + c1 + u, d), "u")
  r2_yu <- pr2(lm(y ~ m + x + c1 + u, d), "u")
  adj <- adjust_paths(sp, r2_xu, r2_mu, r2_yu)

  t_a_true <- summary(lm(m ~ x + c1 + u, d))$coef["x", "t value"]
  t_b_true <- summary(lm(y ~ m + x + c1 + u, d))$coef["m", "t value"]
  expect_equal(adj$t_a, t_a_true, tolerance = 1e-8)
  expect_equal(adj$t_b, t_b_true, tolerance = 1e-8)
})

test_that("joint RV equals the smaller path RV; Sobel RV is smaller", {
  m_fit <- lm(M ~ X + Z1 + Z2, data = dat)
  y_fit <- lm(Y ~ X + M + Z1 + Z2, data = dat)
  sp <- sensitivity_paths(mediator_model = m_fit, outcome_model = y_fit,
                          treatment = "X", mediator = "M")
  rv <- sp$table$rv_alpha
  expect_equal(rv[3], min(rv[1], rv[2]))
  expect_lt(rv[4], rv[3])
  expect_gt(sp$overstatement, 1)
  # at RV(Sobel) the adjusted Sobel z equals the critical value
  adj <- adjust_paths(sp, rv[4], rv[4], rv[4])
  expect_equal(adj$z_sobel, sp$inputs$crit_z, tolerance = 1e-6)
  # at RV(joint) the weaker path sits at its critical value
  adj2 <- adjust_paths(sp, rv[3], rv[3], rv[3])
  expect_equal(min(adj2$t_a - sp$inputs$crit_a, adj2$t_b - sp$inputs$crit_b),
               0, tolerance = 1e-6)
})

test_that("balanced closed form holds", {
  # with equal t and df, g(RV_sobel) = f - sqrt(2) * z / sqrt(df - 1)
  set.seed(3)
  sp <- sensitivity_paths(fit)
  inp <- sp$inputs
  g <- function(r) r / sqrt(1 - r)
  rho <- sp$table$rv_alpha[4]
  ta <- sqrt(inp$df_a - 1) * (abs(inp$t_a) / sqrt(inp$df_a) - g(rho))
  tb <- sqrt(inp$df_b - 1) * (abs(inp$t_b) / sqrt(inp$df_b) - g(rho))
  expect_equal(ta * tb / sqrt(ta^2 + tb^2), inp$crit_z, tolerance = 1e-6)
})

test_that("sensitivity_paths() works from a robmedfit and plots", {
  sp <- sensitivity_paths(fit)
  expect_s3_class(sp, "robmed_paths")
  expect_equal(nrow(sp$table), 4)
  expect_true(any(grepl("Sobel", capture.output(print(sp)))))
  expect_s3_class(plot_paths(sp), "gg")
  expect_error(sensitivity_paths(), "Supply")
})


# ── 15. Estimation details added in 0.3.0 ────────────────────────────────────
test_that("cluster resampling keeps whole clusters", {
  set.seed(1)
  ids <- rep(c("a", "b", "c", "d"), times = c(3, 5, 2, 4))
  idx <- RobustMediate:::.cluster_resample(ids)
  drawn <- ids[idx]
  # every drawn cluster appears with all of its rows, a whole number of times
  tab <- table(drawn) / table(ids)[names(table(drawn))]
  expect_true(all(tab == round(tab)))
  expect_equal(length(unique(sapply(split(idx, ids[idx]), length) %% table(ids)[names(split(idx, ids[idx]))])), 1)
})

test_that("cluster bootstrap runs and is recorded", {
  data(sim_mediation)
  f <- suppressWarnings(robustmediate(
    X ~ Z1 + Z2, M ~ X + Z1 + Z2, Y ~ X + M + Z1 + Z2,
    data = sim_mediation, cluster_var = "school", R = 10, verbose = FALSE,
    dose_grid = seq(-1, 1, length.out = 5)
  ))
  expect_equal(f$cluster$n_clusters, 30)
  expect_null(f$cluster$icc)
  expect_true(all(is.finite(f$effects$curves$lower)))
})

test_that("spline_df is deprecated", {
  expect_warning(
    robustmediate(X ~ Z1 + Z2, M ~ X + Z1 + Z2, Y ~ X + M + Z1 + Z2,
                  data = dat, R = 50, spline_df = 3, verbose = FALSE,
                  dose_grid = c(-1, 0, 1)),
    "deprecated"
  )
})

test_that("splines in the formulas give a curved indirect effect", {
  set.seed(5)
  nn <- 800
  z <- rnorm(nn); x <- rnorm(nn)
  m <- sin(x) + 0.3 * z + rnorm(nn, sd = 0.5)
  y <- 0.5 * m + 0.2 * x + 0.2 * z + rnorm(nn)
  d <- data.frame(x, m, y, z)
  f <- suppressWarnings(robustmediate(
    x ~ z, m ~ splines::ns(x, df = 4) + z, y ~ x + m + z,
    data = d, R = 2, verbose = FALSE, ref_dose = 0,
    dose_grid = c(0.75, 1.5)
  ))
  nie <- f$effects$curves$estimate[f$effects$curves$estimand == "NIE"]
  # exact g-computation: b * (mean E[M | x] - mean E[M | 0])
  mm <- f$models$mediator
  em <- function(v) mean(predict(mm, newdata = transform(d, x = v)))
  b_hat <- coef(f$models$outcome)["m"]
  expect_equal(nie, unname(b_hat * (c(em(0.75), em(1.5)) - em(0))),
               tolerance = 1e-10)
  # curved: a straight line would give NIE(1.5) = 2 * NIE(0.75)
  expect_lt(nie[2], 1.8 * nie[1])
})

test_that("mediator is integrated out for a logit outcome", {
  set.seed(9)
  nn <- 1500
  z <- rnorm(nn); x <- rnorm(nn)
  m <- 0.6 * x + 0.3 * z + rnorm(nn)
  y <- rbinom(nn, 1, plogis(-0.3 + 0.4 * x + 0.8 * m + 0.2 * z))
  d <- data.frame(x, m, y, z)
  f <- suppressWarnings(robustmediate(
    x ~ z, m ~ x + z, y ~ x + m + z, data = d,
    family_outcome = binomial(), R = 2, verbose = FALSE,
    ref_dose = 0, dose_grid = 1
  ))
  # brute-force Monte Carlo over M from the same fitted models
  mm <- f$models$mediator; ym <- f$models$outcome
  sig <- sqrt(summary(mm)$dispersion)
  po <- function(xo, xm) {
    mu <- predict(mm, newdata = transform(d, x = xm))
    mean(sapply(1:400, function(k) {
      dd <- transform(d, x = xo, m = mu + sig * rnorm(nn))
      mean(predict(ym, newdata = dd, type = "response"))
    }))
  }
  set.seed(10)
  nie_mc <- po(1, 1) - po(1, 0)
  nie <- f$effects$curves$estimate[f$effects$curves$estimand == "NIE"]
  expect_equal(nie, nie_mc, tolerance = 0.01)
  # plug-in would differ
  plug <- mean(predict(ym, newdata = transform(d, x = 1, m = predict(mm, newdata = transform(d, x = 1))), type = "response")) -
          mean(predict(ym, newdata = transform(d, x = 1, m = predict(mm, newdata = transform(d, x = 0))), type = "response"))
  expect_gt(abs(plug - nie_mc), abs(nie - nie_mc))
})

test_that("binary mediator is averaged exactly", {
  set.seed(4)
  nn <- 1000
  z <- rnorm(nn); x <- rnorm(nn)
  m <- rbinom(nn, 1, plogis(0.8 * x + 0.3 * z))
  y <- rbinom(nn, 1, plogis(-0.2 + 0.3 * x + 1.0 * m + 0.2 * z))
  d <- data.frame(x, m, y, z)
  f <- suppressWarnings(robustmediate(
    x ~ z, m ~ x + z, y ~ x + m + z, data = d,
    family_mediator = binomial(), family_outcome = binomial(),
    R = 2, verbose = FALSE, ref_dose = 0, dose_grid = 1
  ))
  mm <- f$models$mediator; ym <- f$models$outcome
  e_y <- function(xo, xm) {
    p  <- predict(mm, newdata = transform(d, x = xm), type = "response")
    y1 <- predict(ym, newdata = transform(d, x = xo, m = 1), type = "response")
    y0 <- predict(ym, newdata = transform(d, x = xo, m = 0), type = "response")
    mean(p * y1 + (1 - p) * y0)
  }
  nie <- f$effects$curves$estimate[f$effects$curves$estimand == "NIE"]
  expect_equal(nie, e_y(1, 1) - e_y(1, 0), tolerance = 1e-10)
})


# ── 16. effect_curvature() ───────────────────────────────────────────────────
test_that("linear models give zero curvature, with a message", {
  expect_message(cv <- effect_curvature(fit), "straight line")
  expect_true(all(cv$curvature == 0))
  expect_true(attr(cv, "linear_by_construction"))
})

test_that("spline models give curvature with bootstrap intervals", {
  set.seed(12)
  nn <- 600
  z <- rnorm(nn); x <- rnorm(nn)
  m <- x^2 / 2 + 0.3 * z + rnorm(nn, sd = 0.5)
  y <- 0.5 * m + 0.2 * x + 0.2 * z + rnorm(nn)
  d <- data.frame(x, m, y, z)
  f <- suppressWarnings(robustmediate(
    x ~ z, m ~ splines::ns(x, df = 3) + z, y ~ x + m + z,
    data = d, R = 30, verbose = FALSE, ref_dose = 0,
    dose_grid = seq(-1.5, 1.5, length.out = 13)
  ))
  cv <- effect_curvature(f, "NIE")
  inner <- 2:12
  expect_false(attr(cv, "linear_by_construction"))
  expect_true(all(is.na(cv$curvature[c(1, 13)])))
  expect_true(all(is.finite(cv$curv_lower[inner])))
  expect_true(all(cv$curv_lower[inner] <= cv$curv_upper[inner]))
  # true NIE(x, 0) = 0.5 * x^2 / 2, so curvature is about 0.5 near the centre
  expect_equal(cv$curvature[7], 0.5, tolerance = 0.3)
  expect_s3_class(plot_curvature(cv, ref_dose = 0), "gg")
})

test_that("second differences are exact for a quadratic on an unequal grid", {
  xg <- c(0, 0.3, 1, 1.2, 2.5)
  expect_equal(RobustMediate:::.second_diff(xg, 3 * xg^2)[2:4], rep(6, 3))
})

test_that("sensitivity_curvature() is deprecated", {
  expect_warning(suppressMessages(sensitivity_curvature(fit)), "deprecated")
})
