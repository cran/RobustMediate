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

fit <- robustmediate(
  treatment_formula = X ~ Z1 + Z2,
  mediator_formula  = M ~ X + Z1 + Z2,
  outcome_formula   = Y ~ X + M + Z1 + Z2,
  data    = dat,
  R       = 20,       # keep fast for testing
  verbose = FALSE
)


# ── 1. robustmediate() basic structure ────────────────────────────────────────
test_that("robustmediate() returns a robmedfit object with correct slots", {
  expect_s3_class(fit, "robmedfit")
  expect_true(all(c("models", "balance", "effects", "sensitivity",
                    "meditcv", "meditcv_profile", "cluster", "meta")
                  %in% names(fit)))
})

test_that("meta slot contains expected fields", {
  expect_named(fit$meta,
               c("call", "treat_var", "med_var", "out_var", "ref_dose",
                 "dose_grid", "spline_df", "n_obs", "R", "alpha", "timestamp"),
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
  expect_true(bal$balance$summary_stats$treatment$max_smd >= 0 ||
                is.numeric(bal$balance$summary_stats$treatment$max_smd))
  # max_smd is a non-negative numeric
  expect_gte(fit$balance$summary_stats$treatment$max_smd, 0)
  expect_gte(fit$balance$summary_stats$mediator$max_smd,  0)
})


# ── 3. Sensitivity surface ────────────────────────────────────────────────────
test_that("sensitivity surface has evalue and rho axes", {
  surf <- fit$sensitivity$surface
  expect_true(all(c("evalue", "rho", "effect") %in% names(surf)))
  tip  <- fit$sensitivity$tipping
  expect_true(all(c("evalue_NIE", "rho_min") %in% names(tip)))
  expect_gte(tip$evalue_NIE, 1)
})


# ── 4. meditcv slot ──────────────────────────────────────────────────────────
test_that("sensitivity_meditcv() returns meditcv with correct structure", {
  med <- fit$meditcv
  expect_s3_class(med, "meditcv")
  expect_named(med, c("a_path", "b_path", "indirect", "alpha"), ignore.order = TRUE)
  expect_true(is.numeric(med$a_path$meditcv))
  expect_true(is.numeric(med$b_path$meditcv))
  # minimum robustness principle
  expect_equal(med$indirect$meditcv_indirect,
               min(med$a_path$meditcv, med$b_path$meditcv))
})

test_that("sensitivity_meditcv() can be called directly", {
  med2 <- sensitivity_meditcv(fit)
  expect_s3_class(med2, "meditcv")
})


# ── 5. plot_meditcv ───────────────────────────────────────────────────────────
test_that("plot_meditcv() returns a ggplot", {
  p <- plot_meditcv(fit$meditcv)
  expect_s3_class(p, "gg")
})

test_that("plot(fit, type='meditcv') works via S3 router", {
  p <- plot(fit, type = "meditcv")
  expect_s3_class(p, "gg")
})


# ── 6. Bootstrap CIs ─────────────────────────────────────────────────────────
test_that("bootstrap CIs are non-degenerate (lower != upper) with R >= 20", {
  eff <- fit$effects$summary
  # With R=20 reps and real data, lo != hi
  expect_false(isTRUE(all.equal(eff$NDE_lo, eff$NDE_hi)))
  expect_false(isTRUE(all.equal(eff$NIE_lo, eff$NIE_hi)))
  expect_false(isTRUE(all.equal(eff$TE_lo,  eff$TE_hi)))
})


# ── 7. meditcv_profile slot ──────────────────────────────────────────────────
test_that("sensitivity_meditcv_profile() returns correct class and minimum robustness principle holds", {
  mp <- fit$meditcv_profile
  expect_s3_class(mp, "meditcv_profile")
  expect_true(all(c("a_path", "b_path", "meditcv_indirect",
                    "bottleneck", "fragility_ratio", "robustness_profile",
                    "meditcv_detail") %in% names(mp)))
  # minimum robustness principle
  expect_equal(mp$meditcv_indirect, min(mp$a_path, mp$b_path))
  expect_true(mp$meditcv_indirect <= mp$a_path)
  expect_true(mp$meditcv_indirect <= mp$b_path)
})

test_that("sensitivity_meditcv_profile() can be called directly", {
  mp2 <- sensitivity_meditcv_profile(fit)
  expect_s3_class(mp2, "meditcv_profile")
})


# ── 8. fragility_table ───────────────────────────────────────────────────────
test_that("fragility_table() returns correct columns and 3 rows", {
  ft <- fragility_table(fit)
  expect_s3_class(ft, "data.frame")
  expect_equal(nrow(ft), 3)
  expect_true(all(c("medITCV", "medITCV_pct", "bottleneck", "fragility",
                    "tipping_r_confounder") %in% names(ft)))
  # exactly one bottleneck among a-path and b-path rows
  expect_equal(sum(ft$bottleneck[1:2]), 1)
  # indirect row is always bottleneck = TRUE
  expect_true(ft$bottleneck[3])
})


# ── 9. plot_meditcv_profile and plot_curvature ───────────────────────────────
test_that("plot_meditcv_profile() and plot_curvature() return ggplot objects", {
  p1 <- plot_meditcv_profile(fit$meditcv_profile)
  expect_s3_class(p1, "gg")

  curv <- sensitivity_curvature(fit, estimand = "NIE")
  p2   <- plot_curvature(curv, ref_dose = fit$meta$ref_dose)
  expect_s3_class(p2, "gg")
})

test_that("plot(fit, type='meditcv_profile') works via S3 router", {
  p <- plot(fit, type = "meditcv_profile")
  expect_s3_class(p, "gg")
})


# ── 10. diagnose ─────────────────────────────────────────────────────────────
test_that("diagnose() returns list invisibly and prints text", {
  result <- capture.output(out <- diagnose(fit))
  expect_true(length(result) > 0)
  expect_named(out, c("balance", "effects", "sensitivity",
                      "meditcv", "meditcv_profile"),
               ignore.order = TRUE)
})


# ── 11. Print methods ─────────────────────────────────────────────────────────
test_that("print.meditcv() produces output", {
  out <- capture.output(print(fit$meditcv))
  expect_true(length(out) > 0)
  expect_true(any(grepl("medITCV", out)))
})

test_that("print.meditcv_profile() produces output", {
  out <- capture.output(print(fit$meditcv_profile))
  expect_true(length(out) > 0)
  expect_true(any(grepl("medITCV", out)))
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
  expect_true(all(c("n_obs", "R", "pct_mediated", "evalue_NIE") %in% names(gl)))
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
  p <- plot(fit, type = "curvature")
  expect_s3_class(p, "gg")
})
