test_that("target_estimand validates and defaults the bias scale unit", {
  t <- target_estimand("ITT effect", "latent SMD", sesoi = 0.4)
  expect_s3_class(t, "recovery_target")
  expect_equal(t$bias_scale_unit, 0.4)
  expect_error(target_estimand("x", "SMD", sesoi = 0), "sesoi")
  expect_error(target_estimand(1, "SMD", sesoi = 0.4), "estimand")
})

test_that("data strategy constructors validate", {
  ta <- two_arm_trial(100)
  expect_s3_class(ta, "recovery_two_arm")
  expect_equal(ta$baseline_outcome_cor, 0.5)
  expect_error(two_arm_trial(1), "n_per_arm")
  expect_error(two_arm_trial(100, allocation = 1), "allocation")
  expect_error(two_arm_trial(100, baseline_outcome_cor = 1),
               "baseline_outcome_cor")
  expect_s3_class(cluster_trial(16, 30, icc = 0.05), "recovery_cluster")
  expect_error(cluster_trial(2, 30, icc = 0.05), "n_clusters")
  expect_error(cluster_trial(16, 30, icc = 0.05, icc_pessimistic = 0.01),
               "icc_pessimistic")
})

test_that("measurement and attrition constructors validate", {
  expect_s3_class(measured_outcome(0.7), "recovery_measurement")
  expect_error(measured_outcome(1.2), "reliability")
  a <- attrition_model(0.15)
  expect_equal(a$mechanism, "differential")
  expect_equal(a$rate_control, 0.15)
  expect_equal(a$rate_treated, 0.15)
  expect_equal(a$baseline_slope_treated, -0.5)
  expect_equal(a$baseline_slope_control, 0)
  expect_error(attrition_model(1.0), "rate")
  # mcar zeroes the slopes
  m <- attrition_model(0.15, mechanism = "mcar",
                       baseline_slope_treated = -1)
  expect_equal(m$baseline_slope_treated, 0)
})

test_that("planned_analysis records inference and failure-count choices", {
  pa <- planned_analysis("lmm_random_intercept",
                         y_observed ~ treatment + (1 | cluster),
                         inference = "kenward_roger")
  expect_equal(pa$inference, "kenward_roger")
  expect_false(pa$degenerate_counts)
  expect_error(planned_analysis("linear_model", "not a formula"), "formula")
})

test_that("declare_recovery records omissions and defaults the effect", {
  d <- declare_recovery(
    target = target_estimand("ITT effect", "latent SMD", sesoi = 0.4),
    data_strategy = two_arm_trial(100),
    answer_strategy = planned_analysis("linear_model",
                                       y_observed ~ treatment)
  )
  expect_equal(d$effect, 0.4)
  expect_true(any(grepl("perfectly reliable", d$omissions)))
  expect_true(any(grepl("Attrition was not modeled", d$omissions)))
})

test_that("cluster design with independent-observation analysis warns", {
  expect_warning(
    declare_recovery(
      target = target_estimand("ITT effect", "SMD", sesoi = 0.4),
      data_strategy = cluster_trial(16, 30, icc = 0.05,
                                    icc_pessimistic = 0.15),
      answer_strategy = planned_analysis("linear_model",
                                         y_observed ~ treatment)
    ),
    "independent"
  )
})

test_that("cluster estimators require a cluster data strategy", {
  expect_error(
    declare_recovery(
      target = target_estimand("ITT effect", "SMD", sesoi = 0.4),
      data_strategy = two_arm_trial(100),
      answer_strategy = planned_analysis("cluster_mean_ttest",
                                         y_observed ~ treatment)
    ),
    "cluster_trial"
  )
})
