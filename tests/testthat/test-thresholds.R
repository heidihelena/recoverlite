test_that("shipped profiles match the paper's section 2.5 exactly", {
  d <- recovery_thresholds("default")
  expect_equal(d$null_rejection_mult, 1.25)
  expect_equal(d$power, 0.80)
  expect_equal(d$target_bias, 0.05)
  expect_equal(d$coverage, 0.925)
  expect_equal(d$type_s, 0.01)
  expect_equal(d$type_m, 1.50)
  expect_equal(d$model_failure, 0.01)
  expect_equal(d$mcse_margin, 2)
  expect_equal(d$min_conditional_n, 200)
  expect_equal(d$overcoverage_flag, 0.975)
  expect_length(d$modified, 0)
  expect_match(d$version, "recoverlite-thresholds-")

  s <- recovery_thresholds("strict")
  expect_equal(s$null_rejection_mult, 1.10)
  expect_equal(s$power, 0.90)
  expect_equal(s$target_bias, 0.025)
  expect_equal(s$coverage, 0.940)
  expect_equal(s$type_s, 0.005)
  expect_equal(s$type_m, 1.25)
  expect_equal(s$model_failure, 0.005)

  l <- recovery_thresholds("lenient")
  expect_equal(l$null_rejection_mult, 1.50)
  expect_equal(l$power, 0.70)
  expect_equal(l$target_bias, 0.10)
  expect_equal(l$coverage, 0.900)
  expect_equal(l$type_s, 0.05)
  expect_equal(l$type_m, 2.00)
  expect_equal(l$model_failure, 0.05)
})

test_that("deviations from the shipped profile are recorded", {
  thr <- recovery_thresholds("default", power = 0.90, type_m = 2)
  expect_setequal(thr$modified, c("power", "type_m"))
})

test_that("the estimation profile carries a drift threshold", {
  thr <- recovery_thresholds("estimation")
  expect_equal(thr$drift, thr$target_bias)
  # regression: drift must track a *modified* bias threshold too
  thr2 <- recovery_thresholds("estimation", target_bias = 0.08)
  expect_equal(thr2$drift, 0.08)
  # unless set explicitly
  thr3 <- recovery_thresholds("estimation", target_bias = 0.08, drift = 0.03)
  expect_equal(thr3$drift, 0.03)
})
