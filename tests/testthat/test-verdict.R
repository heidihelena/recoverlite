# Unit tests of the verdict rule on synthetic diagnosand tables, so the
# PASS/RISK/FAIL logic is tested independently of the simulation engine.

make_diag <- function(row_type = "target",
                      rejection = if (row_type == "null") 0.05 else 0.90,
                      rej_mcse = 0.005,
                      bias = 0.00, bias_mcse = 0.002,
                      est_bias = 0.00, drift = 0.00, drift_mcse = 0.001,
                      coverage = 0.955, coverage_mcse = 0.005,
                      type_s = 0.000, type_m = 1.05, n_sig = 1800,
                      cond_mcse = 0.002, model_failure = 0.000,
                      mf_mcse = 0.0005, unstable_cond = FALSE, S = 2000) {
  is_target <- row_type == "target"
  data.frame(
    diagnosand = c("rejection_rate", "target_bias", "estimator_bias",
                   "coverage", "analyzable_coverage", "type_s", "type_m",
                   "precision", "model_failure", "estimand_drift"),
    value = c(rejection, bias, est_bias, coverage, 0.95,
              if (is_target) type_s else NA_real_,
              if (is_target) type_m else NA_real_,
              0.5, model_failure, drift),
    mcse = c(rej_mcse, bias_mcse, bias_mcse, coverage_mcse, 0.005,
             if (is_target) cond_mcse else NA_real_,
             if (is_target) cond_mcse else NA_real_,
             0.01, mf_mcse, drift_mcse),
    n_contributing = c(S, S, S, S, S,
                       if (is_target) n_sig else 0L,
                       if (is_target) n_sig else 0L, S, S, S),
    unstable = c(FALSE, FALSE, FALSE, FALSE, FALSE,
                 if (is_target) unstable_cond else FALSE,
                 if (is_target) unstable_cond else FALSE,
                 FALSE, FALSE, FALSE),
    # one-sided Wilson upper bound, populated only for a zero-count Type S
    upper = c(NA_real_, NA_real_, NA_real_, NA_real_, NA_real_,
              if (is_target && type_s == 0) wilson_upper(0, n_sig) else NA_real_,
              NA_real_, NA_real_, NA_real_, NA_real_),
    note = "", stringsAsFactors = FALSE
  )
}

make_run <- function(name, row_type, counts_for, diag) {
  list(scenario = list(name = name, label = name, rationale = "",
                       counts_for = counts_for, row_type = row_type),
       diagnosands = diag, theta = if (row_type == "null") 0 else 0.4)
}

make_result <- function(nd = make_diag("null"), np = make_diag("null"),
                        td = make_diag("target"),
                        tp = make_diag("target"),
                        thresholds = recovery_thresholds(),
                        drop = character(0)) {
  runs <- list(
    null_declared = make_run("null_declared", "null", "declared", nd),
    null_pessimistic = make_run("null_pessimistic", "null", "pessimistic", np),
    target_declared = make_run("target_declared", "target", "declared", td),
    target_pessimistic = make_run("target_pessimistic", "target",
                                  "pessimistic", tp)
  )
  runs <- runs[setdiff(names(runs), drop)]
  structure(list(design = NULL, runs = runs, sims = 2000L,
                 thresholds = thresholds, seed = 1, alpha = 0.05,
                 scenario_request = "confirmatory_grid",
                 evidence_tiers = character(0),
                 elapsed_secs = 0, session = "test"),
            class = "recovery_result")
}

test_that("all required rows passing with wide margins => PASS", {
  v <- verdict(make_result())
  expect_equal(v$verdict, "PASS")
  expect_true(is.na(v$binding))
  expect_equal(v$verdict_lenient, "PASS")
})

test_that("inflated target-null rejection under Null-declared => FAIL", {
  # 1.25 x 0.05 = 0.0625; 0.09 clearly exceeds it
  v <- verdict(make_result(nd = make_diag("null", rejection = 0.09)))
  expect_equal(v$verdict, "FAIL")
  expect_match(v$binding, "target_null_rejection")
})

test_that("failure under a declared-nuisance target row => FAIL", {
  v <- verdict(make_result(td = make_diag("target", rejection = 0.55)))
  expect_equal(v$verdict, "FAIL")
  expect_match(v$binding, "power")
  v2 <- verdict(make_result(td = make_diag("target", bias = -0.17)))
  expect_equal(v2$verdict, "FAIL")
  expect_match(v2$binding, "target_bias")
  v3 <- verdict(make_result(td = make_diag("target", coverage = 0.90)))
  expect_equal(v3$verdict, "FAIL")
  expect_match(v3$binding, "coverage")
})

test_that("pass declared rows but fail a pessimistic row => RISK", {
  v <- verdict(make_result(tp = make_diag("target", rejection = 0.70)))
  expect_equal(v$verdict, "RISK")
  expect_match(v$binding, "pessimistic")
})

test_that("margin within 2 MCSE of a threshold => RISK, not PASS", {
  # power 0.805 with MCSE 0.009: margin 0.005 < 2 * 0.009
  v <- verdict(make_result(td = make_diag("target", rejection = 0.805,
                                          rej_mcse = 0.009)))
  expect_equal(v$verdict, "RISK")
  expect_match(v$binding, "within 2 MCSE")
})

test_that("declared failure WITHIN 2 MCSE is capped at RISK, not FAIL", {
  # power 0.798 vs 0.80 with MCSE 0.009: fails the point estimate but the
  # margin (-0.002) is within 2 MCSE, so the stability guard caps at RISK.
  v <- verdict(make_result(td = make_diag("target", rejection = 0.798,
                                          rej_mcse = 0.009)))
  expect_equal(v$verdict, "RISK")
  expect_match(v$binding, "caps the verdict at RISK")
})

test_that("declared failure BEYOND 2 MCSE is a stable FAIL", {
  # power 0.78 vs 0.80 with MCSE 0.005: margin -0.02 is 4 MCSE => FAIL.
  v <- verdict(make_result(td = make_diag("target", rejection = 0.78,
                                          rej_mcse = 0.005)))
  expect_equal(v$verdict, "FAIL")
  expect_match(v$binding, "stable failure")
})

test_that("Type S zero-count is checked against the Wilson upper bound", {
  # 0 sign flips out of 210 significant sims: point estimate 0, but the
  # one-sided Wilson upper (~0.0127) exceeds the default 0.01 threshold, so
  # the target row cannot clear and the verdict is capped below PASS.
  v <- verdict(make_result(td = make_diag("target", type_s = 0,
                                          n_sig = 210)))
  expect_true(v$verdict %in% c("RISK", "FAIL"))
  # with many significant sims the Wilson upper is well under 0.01 => PASS
  v2 <- verdict(make_result(td = make_diag("target", type_s = 0,
                                           n_sig = 1800)))
  expect_equal(v2$verdict, "PASS")
})

test_that("overcoverage is flagged as inefficiency, not failure", {
  v <- verdict(make_result(td = make_diag("target", coverage = 0.985,
                                          rejection = 0.90)))
  expect_false(v$verdict == "FAIL")
  ev <- v$evaluations$target_declared
  expect_true(ev$pass[ev$criterion == "coverage"])
  expect_match(ev$note[ev$criterion == "coverage"], "inefficiency")
})

test_that("unstable conditional diagnosands block PASS but do not FAIL", {
  d <- make_diag("target", type_m = 2.4, n_sig = 20, unstable_cond = TRUE)
  v <- verdict(make_result(td = d))
  expect_equal(v$verdict, "RISK")
  expect_match(v$binding, "unstable")
})

test_that("counted model failure above threshold FAILs a declared row", {
  v <- verdict(make_result(td = make_diag("target", model_failure = 0.08)))
  expect_equal(v$verdict, "FAIL")
  expect_match(v$binding, "model_failure")
})

test_that("missing required rows cap the verdict at RISK", {
  v <- verdict(make_result(drop = c("null_declared", "null_pessimistic")))
  expect_equal(v$verdict, "RISK")
  expect_match(v$binding, "Required scenario rows")
})

test_that("strict and lenient recomputations can disagree", {
  # power 0.75: fails default (0.80), passes lenient (0.70)
  v <- verdict(make_result(td = make_diag("target", rejection = 0.75),
                           tp = make_diag("target", rejection = 0.75)))
  expect_equal(v$verdict, "FAIL")
  expect_equal(v$verdict_lenient, "PASS")
  expect_equal(v$verdict_strict, "FAIL")
})

test_that("the estimation profile evaluates target rows only, with drift", {
  thr <- recovery_thresholds("estimation")
  # drift 0.08 Delta exceeds the 0.05 drift threshold
  v <- verdict(make_result(
    td = make_diag("target", drift = 0.08),
    thresholds = thr,
    drop = c("null_declared", "null_pessimistic")))
  expect_equal(v$verdict, "FAIL")
  expect_match(v$binding, "estimand_drift")
  # and no rejection-rate criterion is evaluated
  expect_false("power" %in% v$evaluations$target_declared$criterion)
})
