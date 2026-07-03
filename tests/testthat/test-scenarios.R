design_31 <- function() {
  declare_recovery(
    target = target_estimand("ITT effect", "latent SMD", sesoi = 0.4),
    data_strategy = two_arm_trial(115),
    measurement = measured_outcome(0.70),
    missingness = attrition_model(0.15, mechanism = "differential"),
    answer_strategy = planned_analysis("linear_model",
                                       y_observed ~ treatment)
  )
}

test_that("the confirmatory grid crosses effect and nuisance", {
  scs <- recoverlite:::build_scenarios(design_31())
  expect_named(scs, c("null_declared", "null_pessimistic",
                      "target_declared", "target_pessimistic"))
  expect_equal(scs$null_declared$params$effect, 0)
  expect_equal(scs$null_pessimistic$params$effect, 0)
  expect_equal(scs$target_declared$params$effect, 0.4)
  expect_equal(scs$target_pessimistic$params$effect, 0.4)
  expect_equal(scs$null_declared$row_type, "null")
  expect_equal(scs$target_declared$counts_for, "declared")
  expect_equal(scs$target_pessimistic$counts_for, "pessimistic")
  # the null world keeps the declared missingness mechanism
  expect_equal(scs$null_declared$params$attrition$slope_treated, -0.5)
})

test_that("pessimistic rows perturb nuisance only, per the fallback table", {
  scs <- recoverlite:::build_scenarios(design_31())
  p <- scs$target_pessimistic$params
  expect_equal(p$attrition$rate_control, 0.15 * 1.5)
  expect_equal(p$attrition$rate_treated, 0.15 * 1.5)
  expect_equal(p$reliability, 0.70 - 0.10)
  expect_equal(p$effect, 0.4)   # the effect is NEVER shrunk
  d <- scs$target_declared$params
  expect_equal(d$attrition$rate_treated, 0.15)
  expect_equal(d$reliability, 0.70)
})

test_that("attrition perturbation respects the declared cap", {
  d <- declare_recovery(
    target = target_estimand("ITT effect", "SMD", sesoi = 0.4),
    data_strategy = two_arm_trial(100),
    missingness = attrition_model(0.5, max_rate = 0.6),
    answer_strategy = planned_analysis("linear_model",
                                       y_observed ~ treatment)
  )
  scs <- recoverlite:::build_scenarios(d)
  expect_equal(scs$target_pessimistic$params$attrition$rate_treated, 0.6)
})

test_that("declared effect above the SESOI: target rows at SESOI, expected-effect row informational", {
  d <- declare_recovery(
    target = target_estimand("ITT effect", "SMD", sesoi = 0.4),
    data_strategy = two_arm_trial(100),
    answer_strategy = planned_analysis("linear_model",
                                       y_observed ~ treatment),
    effect = 0.6
  )
  scs <- recoverlite:::build_scenarios(d)
  expect_equal(scs$target_declared$params$effect, 0.4)      # SESOI
  expect_equal(scs$target_pessimistic$params$effect, 0.4)
  expect_true("expected_effect" %in% names(scs))
  expect_equal(scs$expected_effect$params$effect, 0.6)
  expect_equal(scs$expected_effect$counts_for, "informational")
})

test_that("cluster ICC moves to the declared upper bound with its tier", {
  d <- suppressWarnings(declare_recovery(
    target = target_estimand("ITT effect", "SMD", sesoi = 0.4),
    data_strategy = cluster_trial(16, 30, icc = 0.05,
                                  icc_pessimistic = 0.15,
                                  evidence = "published ICC range"),
    answer_strategy = planned_analysis("lmm_random_intercept",
                                       y_observed ~ treatment + (1 | cluster))
  ))
  scs <- recoverlite:::build_scenarios(d)
  expect_equal(scs$target_pessimistic$params$icc, 0.15)
  expect_true(any(grepl("published ICC range", attr(scs, "tiers"))))
})

test_that("undeclared ICC upper bound falls back to icc + 0.10, labeled package default", {
  d <- declare_recovery(
    target = target_estimand("ITT effect", "SMD", sesoi = 0.4),
    data_strategy = cluster_trial(16, 30, icc = 0.05),
    answer_strategy = planned_analysis("lmm_random_intercept",
                                       y_observed ~ treatment + (1 | cluster))
  )
  scs <- recoverlite:::build_scenarios(d)
  expect_equal(scs$target_pessimistic$params$icc, 0.15)
  expect_true(any(grepl("package default", attr(scs, "tiers"))))
})

test_that("target_grid produces the two target rows only", {
  scs <- recoverlite:::build_scenarios(design_31(), "target_grid")
  expect_named(scs, c("target_declared", "target_pessimistic"))
})

test_that("dropout intercept calibration hits the marginal rate", {
  a <- recoverlite:::calibrate_dropout_intercept(0.15, -0.5)
  marg <- stats::integrate(function(b)
    stats::plogis(a - 0.5 * b) * stats::dnorm(b), -Inf, Inf)$value
  expect_equal(marg, 0.15, tolerance = 1e-4)
  expect_equal(recoverlite:::calibrate_dropout_intercept(0.15, 0),
               stats::qlogis(0.15))
})
