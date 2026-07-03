#' Produce the standalone recovery report
#'
#' Renders the full recovery report (protocol Step 6): target estimand
#' and SESOI with scales, the declaration (including explicit statements
#' of what was *not* modeled), the scenario grid with the null-world
#' specification and the evidence tier of every pessimistic value, all
#' diagnosands under every scenario with MCSEs, contributing counts,
#' failure-class counts, and the bias decomposition, the threshold
#' profile with signed margins, the verdict under the selected and the
#' shipped strict/lenient profiles, the binding failure mode, and the
#' computation settings. The report is standalone: a collaborator,
#' reviewer, supervisor, or funder should be able to read it without
#' running the code. A FAIL report is as complete as a PASS report.
#'
#' @param result A `recovery_result` from [recovery_test()].
#' @param file Optional path; if supplied, the report is also written
#'   there as plain text.
#'
#' @return The report, invisibly, as a character vector of lines.
#' @export
report <- function(result, file = NULL) {
  stopifnot(
    "`result` must be a recovery_result from recovery_test()" =
      inherits(result, "recovery_result")
  )
  d <- result$design
  v <- verdict(result)
  thr <- result$thresholds
  L <- character(0)
  add <- function(...) L <<- c(L, paste0(...))

  add(strrep("=", 74))
  add("PRE-DATA RECOVERY REPORT  (recoverlite ",
      as.character(utils::packageVersion("recoverlite")), ")")
  add(strrep("=", 74))
  add("")
  add("1. TARGET")
  add("   Estimand: ", d$target$estimand)
  add("   Scale:    ", d$target$scale)
  add(sprintf("   SESOI:    %.3g   Declared expected effect: %.3g   Delta (bias/drift unit): %.3g",
              d$target$sesoi, d$effect, d$target$bias_scale_unit))
  if (!is.null(d$target$max_width)) {
    add(sprintf("   Declared maximum acceptable CI width: %.3g",
                d$target$max_width))
  }
  add("")
  add("2. DECLARATION")
  for (ln in describe_design(d)) add("   ", ln)
  if (length(d$omissions)) {
    add("   Not modeled (silence must not imply ideality):")
    for (om in d$omissions) add("     - ", om)
  }
  add("")
  add("3. SCENARIO GRID AND NULL-WORLD SPECIFICATION")
  for (run in result$runs) {
    add("   * ", run$scenario$label)
    add("     ", run$scenario$rationale)
    if (run$scenario$row_type == "null") {
      dr <- run$diagnosands
      drift <- dr[dr$diagnosand == "estimand_drift", ]
      lab <- if (is.finite(drift$mcse) &&
                 abs(drift$value) > 2 * drift$mcse) {
        sprintf(paste0("E[theta_obs] != 0 under this null world (drift ",
                       "%.4f [%.4f] Delta): rejections are FALSE CLAIMS ",
                       "ABOUT THE TARGET, partly induced by selection."),
                drift$value, drift$mcse)
      } else {
        "E[theta_obs] = 0 under this null world: the rejection rate is pure test size."
      }
      add("     Null world: theta = 0 with the declared missingness ",
          "mechanism persisting. ", lab)
    }
  }
  if (length(result$evidence_tiers)) {
    add("   Evidence tiers of the pessimistic values (hierarchy: empirical",
        " > prior-study > elicited > package default):")
    for (tier in result$evidence_tiers) add("     - ", tier)
  }
  add("")
  add("4. DIAGNOSANDS (value [MCSE]; n = contributing simulations)")
  add("   Inclusion rules: counted model failures are excluded from")
  add("   estimate-based diagnosands, which are conditional on successful")
  add("   analysis; failure classes are fatal/nonconvergence/degenerate/")
  add("   warnings (F/N/D/W); target bias = estimator bias + estimand drift.")
  for (nm in names(result$runs)) {
    run <- result$runs[[nm]]
    diag <- run$diagnosands
    add("")
    add("   -- ", run$scenario$label)
    add(sprintf("      analyzed n per sim: mean %.1f; realized attrition: %.3f",
                attr(diag, "mean_n_analyzed"),
                attr(diag, "mean_attrition_realized")))
    for (i in seq_len(nrow(diag))) {
      r <- diag[i, ]
      val <- if (is.na(r$value)) "  not estimable " else
        sprintf("%8.4f [%.4f]", r$value, r$mcse)
      add(sprintf("      %-20s %s  n=%-5d%s%s",
                  r$diagnosand, val, r$n_contributing,
                  if (r$unstable) "  UNSTABLE" else "",
                  if (nzchar(r$note)) paste0("  (", r$note, ")") else ""))
    }
  }
  add("")
  add("5. THRESHOLD PROFILE AND SIGNED MARGINS")
  add("   Profile: '", thr$profile, "' [", thr$version, "]",
      if (length(thr$modified))
        paste0("  DEVIATIONS from shipped profile: ",
               paste(thr$modified, collapse = ", "))
      else "  (shipped values)")
  add("   Signed margin to every threshold (positive = passing):")
  for (nm in names(v$evaluations)) {
    ev <- v$evaluations[[nm]]
    for (i in seq_len(nrow(ev))) {
      r <- ev[i, ]
      add(sprintf("      %-20s %-22s %s  margin %+.4f [MCSE %.4f]%s%s",
                  r$criterion, paste0("(", nm, ")"), r$requirement,
                  r$margin, r$mcse,
                  if (isTRUE(r$unstable)) "  UNSTABLE" else "",
                  if (nzchar(r$note)) paste0("  (", r$note, ")") else ""))
    }
  }
  add("")
  add("6. VERDICT: ", v$verdict, "  (profile '", thr$profile, "')")
  if (!is.na(v$verdict_strict)) {
    add(sprintf("   Under shipped profiles: strict %s | default-family %s | lenient %s",
                v$verdict_strict, v$verdict, v$verdict_lenient))
    if (length(unique(c(v$verdict_strict, v$verdict, v$verdict_lenient))) > 1) {
      add("   Profile disagreement is itself a finding; the RISK category exists to hold it.")
    }
  }
  if (!is.na(v$smallest_margin)) {
    add("   Smallest signed margin: ", v$smallest_margin)
  }
  add(sprintf(paste0("   Rule: PASS = all required rows pass with margins > ",
                     "%g MCSE; RISK = pessimistic-only failure, narrow ",
                     "margin, or unstable required diagnosand; FAIL = any ",
                     "failure under a declared-nuisance row."),
              thr$mcse_margin))
  add("   The verdict is a decision convention, not a validity classification.")
  add("")
  add("7. BINDING FAILURE MODE")
  add("   ", if (is.na(v$binding)) "None: all criteria passed with stable margins."
      else v$binding)
  if (v$verdict != "PASS") {
    add("   Fragility curves for the binding parameters: see",
        " nuisance_fragility() and effect_fragility().")
  }
  add("")
  add("8. DESIGN CHANGE")
  if (v$verdict == "PASS") {
    add("   No change required under the scenario rows this profile requires.")
  } else {
    add("   Rerun the recovery test on candidate repairs, each simulated as")
    add("   a full declaration (e.g. improved measurement reliability,")
    add("   larger recruitment, more clusters, a different inference method,")
    add("   or an estimator matched to the estimand such as")
    add("   mi_baseline_adjusted) until the binding criterion passes with")
    add("   stable margins. Note the repair must match the failure: a")
    add("   larger complete-case study estimates a displaced contrast more")
    add("   precisely - resources repair precision, not drift.")
  }
  add("")
  add(strrep("-", 74))
  add(sprintf("Computation: %d simulations per scenario row; alpha = %.3g; seed = %s;",
              result$sims, result$alpha,
              if (is.null(result$seed)) "not set (NOT reproducible)"
              else as.character(result$seed)))
  add(sprintf("elapsed %.1f s; %s.", result$elapsed_secs, result$session))
  add("A PASS is evidence about the instrument, not about the world.")
  add(strrep("=", 74))

  cat(L, sep = "\n")
  if (!is.null(file)) writeLines(L, file)
  invisible(L)
}

# Shared, human-readable summary of the declaration (used by report() and
# print.recovery_design()).
describe_design <- function(d) {
  ds <- d$data_strategy
  a <- d$answer_strategy
  lines <- character(0)
  if (inherits(ds, "recovery_two_arm")) {
    lines <- c(lines, sprintf(
      "Data strategy: two-arm randomized trial; %d recruited per arm (allocation %.2g); baseline-outcome correlation rho = %.3g (baseline observed at randomization for every participant).",
      ds$n_per_arm, ds$allocation, ds$baseline_outcome_cor))
    if (ds$noncompliance > 0) {
      lines <- c(lines, sprintf(
        "One-sided noncompliance: %.3g of treated receive no treatment.",
        ds$noncompliance))
    }
  } else {
    lines <- c(lines, sprintf(
      "Data strategy: cluster-randomized trial; %d clusters x %d individuals, sizes fixed by design (allocation %.2g); declared ICC %.3g%s.",
      ds$n_clusters, ds$n_per_cluster, ds$allocation, ds$icc,
      if (is.null(ds$icc_pessimistic))
        " (no upper plausible bound declared)"
      else sprintf(", upper plausible bound %.3g", ds$icc_pessimistic)))
  }
  if (!is.null(d$measurement)) {
    lines <- c(lines, sprintf(
      "Measurement: classical additive error, Var(e) = 1/r - 1 with declared reliability r = %.3g. The raw treatment contrast is not attenuated in expectation; the error inflates residual variance (charged to the variance account, not the bias account).",
      d$measurement$reliability))
  }
  if (!is.null(d$missingness) && d$missingness$rate > 0) {
    m <- d$missingness
    if (m$mechanism == "differential") {
      lines <- c(lines, sprintf(
        "Attrition: logit Pr(dropout | Z, B) = alpha_Z + gamma_Z B, with intercepts calibrated to marginal rates control %.3g / treated %.3g and baseline slopes gamma_0 = %.3g, gamma_1 = %.3g. Dropout depends only on observed quantities (Z, B): MAR given the baseline.",
        m$rate_control, m$rate_treated, m$baseline_slope_control,
        m$baseline_slope_treated))
    } else {
      lines <- c(lines, sprintf(
        "Attrition: MCAR at marginal rates control %.3g / treated %.3g.",
        m$rate_control, m$rate_treated))
    }
  }
  inference_txt <- if (a$estimator == "lmm_random_intercept") {
    paste0(" (inference: ", a$inference, ")")
  } else if (a$estimator == "mi_baseline_adjusted") {
    sprintf(" (m = %d imputations, Rubin's rules, Barnard-Rubin df)",
            a$m_imputations)
  } else ""
  lines <- c(lines, sprintf(
    "Answer strategy: %s%s, %s; two-sided alpha %.3g. Degenerate fits (singular/boundary) %s against the model-failure threshold (pre-specified); fatal errors and nonconvergence always count.",
    a$estimator, inference_txt, deparse(a$formula), a$alpha,
    if (a$degenerate_counts) "COUNT" else "do NOT count"))
  lines
}
