# End-to-end smoke tests with small sims (structure, not precision).

two_arm_design <- function(n = 60) {
  declare_recovery(
    target = target_estimand("ITT effect", "latent SMD", sesoi = 0.4),
    data_strategy = two_arm_trial(n),
    measurement = measured_outcome(0.70),
    missingness = attrition_model(0.15, mechanism = "differential"),
    answer_strategy = planned_analysis("linear_model",
                                       y_observed ~ treatment)
  )
}

test_that("two-arm recovery test runs the confirmatory grid end to end", {
  skip_if_not_installed("DeclareDesign")
  res <- recovery_test(two_arm_design(), sims = 80, seed = 42)
  expect_s3_class(res, "recovery_result")
  expect_named(res$runs, c("null_declared", "null_pessimistic",
                           "target_declared", "target_pessimistic"))
  diag <- res$runs$target_declared$diagnosands
  expect_setequal(diag$diagnosand,
                  c("rejection_rate", "target_bias", "estimator_bias",
                    "coverage", "analyzable_coverage", "type_s", "type_m",
                    "precision", "model_failure", "estimand_drift"))
  # lm on complete data should essentially never fail
  expect_lt(diag$value[diag$diagnosand == "model_failure"], 0.05)
  # Type S / Type M are n/a on null rows
  nd <- res$runs$null_declared$diagnosands
  expect_true(is.na(nd$value[nd$diagnosand == "type_m"]))
  # decomposition: target bias = estimator bias + drift (same sims, exact
  # up to floating point)
  tb <- diag$value[diag$diagnosand == "target_bias"]
  eb <- diag$value[diag$diagnosand == "estimator_bias"]
  dr <- diag$value[diag$diagnosand == "estimand_drift"]
  expect_equal(tb, eb + dr, tolerance = 1e-10)

  v <- verdict(res)
  expect_true(v$verdict %in% c("PASS", "RISK", "FAIL"))
  expect_true(v$verdict_strict %in% c("PASS", "RISK", "FAIL"))
  rep_lines <- capture.output(rl <- report(res))
  expect_true(any(grepl("PRE-DATA RECOVERY REPORT", rep_lines)))
  expect_true(any(grepl("SCENARIO GRID", rep_lines)))
  expect_true(any(grepl("SIGNED MARGINS", rep_lines)))
  expect_true(any(grepl("Not modeled", rep_lines)))
})

test_that("baseline-dependent differential attrition displaces the null", {
  skip_if_not_installed("DeclareDesign")
  # Strong selection, decent n: drift under the null must be positive
  # (retained treated over-represent high-baseline participants).
  d <- declare_recovery(
    target = target_estimand("ITT effect", "latent SMD", sesoi = 0.4),
    data_strategy = two_arm_trial(300, baseline_outcome_cor = 0.7),
    missingness = attrition_model(0.3, mechanism = "differential",
                                  baseline_slope_treated = -1.5),
    answer_strategy = planned_analysis("linear_model",
                                       y_observed ~ treatment)
  )
  res <- recovery_test(d, sims = 150, seed = 7)
  nd <- res$runs$null_declared$diagnosands
  drift <- nd[nd$diagnosand == "estimand_drift", ]
  expect_gt(drift$value, 2 * drift$mcse)
  # and the rejection-rate note flags it as target-null rejection
  expect_match(nd$note[nd$diagnosand == "rejection_rate"], "TARGET-NULL")
})

test_that("mcar attrition leaves the null undisplaced (test size)", {
  skip_if_not_installed("DeclareDesign")
  d <- declare_recovery(
    target = target_estimand("ITT effect", "latent SMD", sesoi = 0.4),
    data_strategy = two_arm_trial(100),
    missingness = attrition_model(0.15, mechanism = "mcar"),
    answer_strategy = planned_analysis("linear_model",
                                       y_observed ~ treatment)
  )
  res <- recovery_test(d, sims = 100, seed = 11)
  nd <- res$runs$null_declared$diagnosands
  expect_match(nd$note[nd$diagnosand == "rejection_rate"], "test size")
})

test_that("reproducibility: same seed gives identical diagnosands", {
  skip_if_not_installed("DeclareDesign")
  d <- two_arm_design(40)
  r1 <- recovery_test(d, sims = 30, seed = 7, scenarios = "target_grid")
  r2 <- recovery_test(d, sims = 30, seed = 7, scenarios = "target_grid")
  expect_equal(r1$runs$target_declared$diagnosands$value,
               r2$runs$target_declared$diagnosands$value)
})

test_that("MI baseline-adjusted estimator runs and beats complete-case drift", {
  skip_if_not_installed("DeclareDesign")
  d <- declare_recovery(
    target = target_estimand("ITT effect", "latent SMD", sesoi = 0.4),
    data_strategy = two_arm_trial(80),
    missingness = attrition_model(0.2, mechanism = "differential",
                                  baseline_slope_treated = -1),
    answer_strategy = planned_analysis(
      "mi_baseline_adjusted", y_observed ~ treatment + baseline,
      m_imputations = 10)
  )
  res <- recovery_test(d, sims = 40, seed = 5, scenarios = "target_grid")
  diag <- res$runs$target_declared$diagnosands
  expect_true(is.finite(diag$value[diag$diagnosand == "target_bias"]))
  expect_lt(diag$value[diag$diagnosand == "model_failure"], 0.1)
})

test_that("cluster answer strategies run and classify failures", {
  skip_if_not_installed("DeclareDesign")
  skip_if_not_installed("lme4")
  skip_if_not_installed("lmerTest")
  base <- function(est, inf = "satterthwaite") {
    suppressWarnings(declare_recovery(
      target = target_estimand("ITT pupil effect", "student SMD",
                               sesoi = 0.4),
      data_strategy = cluster_trial(8, 10, icc = 0.05,
                                    icc_pessimistic = 0.15),
      answer_strategy = planned_analysis(
        est,
        if (est == "cluster_mean_ttest") y_observed ~ treatment
        else y_observed ~ treatment + (1 | cluster),
        inference = inf)
    ))
  }
  res <- recovery_test(base("lmm_random_intercept"), sims = 40, seed = 3,
                       scenarios = "target_grid")
  diag <- res$runs$target_declared$diagnosands
  mf <- diag[diag$diagnosand == "model_failure", ]
  # degenerate fits detected and reported in the class breakdown ...
  expect_match(mf$note, "F/N/D/W")
  # ... but not counted by default (pre-specified in planned_analysis)
  expect_match(mf$note, "not counted")
  # singular/boundary fits are class (c) degenerate, NOT class (b)
  # nonconvergence: with this tiny low-ICC design singular fits are
  # common, and they must not leak into the always-counted classes.
  fc <- attr(diag, "failure_classes")
  expect_gt(fc[["degenerate"]], 0)
  expect_lt(fc[["nonconverged"]], fc[["degenerate"]])
  sim_data <- res$runs$target_declared$sim_data
  expect_equal(mf$value, mean(sim_data$counted_failure))
  expect_equal(mean(sim_data$fatal | sim_data$nonconverged), mf$value)

  res_t <- recovery_test(base("cluster_mean_ttest"), sims = 25, seed = 3,
                         scenarios = "target_grid")
  expect_true(is.finite(
    res_t$runs$target_declared$diagnosands$value[1]))

  res_z <- recovery_test(base("lmm_random_intercept", "wald_z"),
                         sims = 25, seed = 3, scenarios = "target_grid")
  expect_true(is.finite(
    res_z$runs$target_declared$diagnosands$value[1]))
})

test_that("fragility curves run and stay outside the verdict", {
  skip_if_not_installed("DeclareDesign")
  d <- two_arm_design(40)
  fr <- effect_fragility(d, effects = c(0.2, 0.4), sims = 25, seed = 1)
  expect_s3_class(fr, "recovery_fragility")
  expect_equal(nrow(fr), 2)
  nf <- nuisance_fragility(d, "attrition_rate", values = c(0.15, 0.3),
                           sims = 25, seed = 1)
  expect_s3_class(nf, "recovery_fragility")
  expect_equal(nf$attrition_rate, c(0.15, 0.3))
})
