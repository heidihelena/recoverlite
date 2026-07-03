# Diagnosands (manuscript, section 2.3), each with Monte Carlo
# uncertainty and explicit inclusion rules (section 2.4).
#
# Inclusion rules, fixed in advance: simulations with COUNTED model
# failures (fatal, nonconvergence, and degenerate fits when the answer
# strategy pre-specifies that they count) contribute to the model-failure
# rate and are excluded from estimate-based diagnosands, which are
# therefore conditional on successful analysis. Degenerate fits that do
# not count are retained, and their marginal effect on coverage is
# reported separately. Conditional diagnosands (Type S, Type M) report
# contributing counts, use a nonparametric bootstrap MCSE that resamples
# over all included simulations and recomputes the significance
# conditioning within each resample, and are marked unstable below the
# pre-specified minimum count.

mcse_prop <- function(p, n) {
  if (is.na(p) || n < 1) return(NA_real_)
  sqrt(p * (1 - p) / n)
}

mcse_mean <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 2) return(NA_real_)
  stats::sd(x) / sqrt(length(x))
}

# Bootstrap MCSE for conditional diagnosands: resample the included
# simulations, recompute the conditioning, recompute the statistic.
mcse_boot <- function(stat_fun, n_included, B = 500) {
  vals <- vapply(seq_len(B), function(b) {
    idx <- sample.int(n_included, n_included, replace = TRUE)
    stat_fun(idx)
  }, numeric(1))
  stats::sd(vals, na.rm = TRUE)
}

# sim_df:   per-simulation results from run_scenario()
# theta:    true value of the target estimand in this scenario row
# unit:     declared substantive unit Delta (SESOI by default)
# row_type: "null" or "target" (Type S / Type M are target-only)
compute_diagnosands <- function(sim_df, theta, unit, thresholds, row_type,
                                has_attrition, alpha) {
  S <- nrow(sim_df)
  ok <- !sim_df$counted_failure & is.finite(sim_df$estimate)
  n_ok <- sum(ok)
  est <- sim_df$estimate[ok]
  sig <- sim_df$sig[ok]
  tobs <- sim_df$theta_obs[ok]
  n_sig <- sum(sig, na.rm = TRUE)
  min_n <- thresholds$min_conditional_n

  row <- function(name, value, mcse, n_contributing, unstable = FALSE,
                  note = "") {
    data.frame(diagnosand = name, value = value, mcse = mcse,
               n_contributing = n_contributing, unstable = unstable,
               note = note, stringsAsFactors = FALSE)
  }

  # Rejection rate: test size / target-null rejection under null rows,
  # power under target rows. The label depends on whether drift leaves
  # the analyzable contrast nonzero under the null (section 2.2, Step 3).
  rej <- mean(sig, na.rm = TRUE)
  drift_val <- mean(tobs - theta)
  drift_mcse <- mcse_mean(tobs - theta)
  rej_note <- if (row_type == "null") {
    if (is.finite(drift_mcse) && abs(drift_val) > 2 * drift_mcse) {
      "TARGET-NULL REJECTION RATE: E[theta_obs] != 0 under this null world (drift below); rejections are false claims about the target, partly induced by selection - not pure test size"
    } else {
      "test size: E[theta_obs] = 0 under this null world"
    }
  } else "power at this scenario's theta"
  d_rej <- row("rejection_rate", rej, mcse_prop(rej, n_ok), n_ok,
               note = rej_note)

  # Target bias and its exact decomposition -----------------------------
  d_tbias <- row("target_bias", mean(est - theta) / unit,
                 mcse_mean(est - theta) / unit, n_ok,
                 note = sprintf("(E[est] - theta) / %.3g; = estimator_bias + estimand_drift", unit))
  d_ebias <- row("estimator_bias", mean(est - tobs) / unit,
                 mcse_mean(est - tobs) / unit, n_ok,
                 note = "bias for the analyzable-data contrast (analysis problem)")

  # Coverage of theta, and of the analyzable contrast (not thresholded) --
  cov <- mean(sim_df$covered[ok], na.rm = TRUE)
  d_cov <- row("coverage", cov, mcse_prop(cov, n_ok), n_ok)
  cov_obs <- mean(sim_df$covered_obs[ok], na.rm = TRUE)
  d_cov_obs <- row("analyzable_coverage", cov_obs, mcse_prop(cov_obs, n_ok),
                   n_ok,
                   note = "Pr(theta_obs in CI): estimator calibration for the contrast it actually estimates (not thresholded)")

  # Type S / Type M: target rows only, conditional on significance -------
  if (row_type == "target") {
    if (n_sig > 0) {
      sig_idx <- which(sig)
      ts <- mean(sign(est[sig_idx]) != sign(theta))
      ts_mcse <- mcse_boot(function(idx) {
        s2 <- sig[idx]; e2 <- est[idx]
        if (!any(s2, na.rm = TRUE)) return(NA_real_)
        mean(sign(e2[which(s2)]) != sign(theta))
      }, n_ok)
      d_ts <- row("type_s", ts, ts_mcse, n_sig,
                  unstable = n_sig < min_n,
                  note = if (n_sig < min_n)
                    sprintf("only %d significant simulations (< %d); unstable",
                            n_sig, min_n) else "bootstrap MCSE")
      tm <- mean(abs(est[sig_idx])) / abs(theta)
      tm_mcse <- mcse_boot(function(idx) {
        s2 <- sig[idx]; e2 <- est[idx]
        if (!any(s2, na.rm = TRUE)) return(NA_real_)
        mean(abs(e2[which(s2)])) / abs(theta)
      }, n_ok)
      d_tm <- row("type_m", tm, tm_mcse, n_sig,
                  unstable = n_sig < min_n,
                  note = if (n_sig < min_n)
                    sprintf("only %d significant simulations (< %d); unstable",
                            n_sig, min_n) else "bootstrap MCSE")
    } else {
      d_ts <- row("type_s", NA_real_, NA_real_, 0L, unstable = TRUE,
                  note = "no significant simulations; not estimable")
      d_tm <- row("type_m", NA_real_, NA_real_, 0L, unstable = TRUE,
                  note = "no significant simulations; not estimable")
    }
  } else {
    d_ts <- row("type_s", NA_real_, NA_real_, 0L,
                note = "n/a: undefined under theta = 0")
    d_tm <- row("type_m", NA_real_, NA_real_, 0L,
                note = "n/a: undefined under theta = 0")
  }

  # Precision -------------------------------------------------------------
  width <- (sim_df$ci_hi - sim_df$ci_lo)[ok]
  prec_note <- "mean CI width"
  if (!is.null(thresholds$max_width)) {
    p_wide <- mean(width > thresholds$max_width, na.rm = TRUE)
    prec_note <- sprintf(
      "mean CI width; Pr(width > declared max %.3g) = %.3f",
      thresholds$max_width, p_wide)
  }
  d_prec <- row("precision", mean(width), mcse_mean(width), n_ok,
                note = prec_note)

  # Model failure: classified, not pooled (F/N/D/W) -----------------------
  mf <- mean(sim_df$counted_failure)
  cov_degen <- if (any(sim_df$degenerate & ok)) {
    sprintf("; coverage among degenerate fits %.3f vs %.3f among others",
            mean(sim_df$covered[sim_df$degenerate & ok], na.rm = TRUE),
            mean(sim_df$covered[!sim_df$degenerate & ok], na.rm = TRUE))
  } else ""
  d_mf <- row(
    "model_failure", mf, mcse_prop(mf, S), S,
    note = sprintf(
      "counted classes only; all classes F/N/D/W = %d/%d/%d/%d of %d (degenerate %s)%s",
      sum(sim_df$fatal), sum(sim_df$nonconverged), sum(sim_df$degenerate),
      sum(sim_df$warned), S,
      if (isTRUE(attr(sim_df, "degenerate_counts"))) "counted" else "not counted",
      cov_degen))

  # Estimand drift ---------------------------------------------------------
  d_drift <- row("estimand_drift", drift_val / unit, drift_mcse / unit, n_ok,
                 note = if (has_attrition)
                   "(E[theta_obs] - theta) / Delta: what the data strategy does to the question (design problem)"
                 else
                   "no attrition or exclusions declared; zero in expectation by construction")

  out <- rbind(d_rej, d_tbias, d_ebias, d_cov, d_cov_obs, d_ts, d_tm,
               d_prec, d_mf, d_drift)
  attr(out, "n_sims") <- S
  attr(out, "n_ok") <- n_ok
  attr(out, "n_sig") <- n_sig
  attr(out, "mean_n_analyzed") <- mean(sim_df$n_analyzed)
  attr(out, "mean_attrition_realized") <- mean(sim_df$attrition_realized)
  attr(out, "failure_classes") <- c(fatal = sum(sim_df$fatal),
                                    nonconverged = sum(sim_df$nonconverged),
                                    degenerate = sum(sim_df$degenerate),
                                    warned = sum(sim_df$warned))
  out
}
