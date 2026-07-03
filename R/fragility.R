#' Effect-size fragility curve
#'
#' Shows how power, Type M exaggeration, and precision behave as the true
#' effect shrinks below the declared target, with all nuisance
#' assumptions held at their declared values. This is deliberately
#' **outside the PASS/RISK/FAIL verdict**: shrinking the target effect
#' changes the inferential target, so pessimistic scenario rows perturb
#' nuisance assumptions only, and fragility to the effect size itself is
#' reported as a separate continuous curve.
#'
#' @param design A [declare_recovery()] object.
#' @param effects Numeric vector of true effects to evaluate. Defaults to
#'   an even grid from a quarter of the declared effect up to the declared
#'   effect, with the SESOI included.
#' @param sims Simulations per effect value (default 500; increase for
#'   smoother curves).
#' @param seed Optional integer seed.
#'
#' @return A data frame of class `recovery_fragility` with one row per
#'   effect: power, Type M, and precision, each with its MCSE.
#' @export
effect_fragility <- function(design, effects = NULL, sims = 500,
                             seed = NULL) {
  stopifnot(
    "`design` must be a declare_recovery() object" =
      inherits(design, "recovery_design")
  )
  s <- sign(design$effect)
  effects <- effects %||%
    sort(unique(c(seq(0.25, 1, length.out = 4) * design$effect,
                  s * design$target$sesoi)))
  stopifnot(is.numeric(effects), length(effects) >= 1, all(effects != 0))
  if (!is.null(seed)) set.seed(seed)

  out <- fragility_rows(design, lapply(effects, function(eff) {
    list(value = eff, overrides = list(effect = eff))
  }), sims, value_name = "effect")
  attr(out, "sims") <- sims
  attr(out, "sesoi") <- design$target$sesoi
  attr(out, "parameter") <- "effect"
  class(out) <- c("recovery_fragility", "data.frame")
  out
}

#' Nuisance-parameter fragility curve
#'
#' Fragility curves over a nuisance parameter the diagnosis identifies as
#' binding (manuscript, section 2.2, Step 3): the verdict uses the
#' pre-specified point scenarios; these curves show where the verdict
#' would change. The chosen diagnosands are power, target-null-relevant
#' quantities being outside scope here, target bias, coverage, and
#' estimand drift across the grid, at the target effect (the SESOI).
#'
#' @param design A [declare_recovery()] object.
#' @param parameter One of `"attrition_rate"`, `"reliability"`, `"icc"`.
#' @param values Numeric grid of parameter values to evaluate.
#' @param sims Simulations per grid point (default 500).
#' @param seed Optional integer seed.
#'
#' @return A data frame of class `recovery_fragility`.
#' @export
nuisance_fragility <- function(design,
                               parameter = c("attrition_rate",
                                             "reliability", "icc"),
                               values, sims = 500, seed = NULL) {
  stopifnot(
    "`design` must be a declare_recovery() object" =
      inherits(design, "recovery_design"),
    "`values` must be a numeric grid" =
      is.numeric(values) && length(values) >= 1
  )
  parameter <- match.arg(parameter)
  if (parameter == "icc" &&
      !inherits(design$data_strategy, "recovery_cluster")) {
    stop("`icc` fragility requires a cluster_trial() data strategy.",
         call. = FALSE)
  }
  if (parameter == "attrition_rate" && is.null(design$missingness)) {
    stop("`attrition_rate` fragility requires a declared attrition_model().",
         call. = FALSE)
  }
  if (!is.null(seed)) set.seed(seed)
  theta <- sign(design$effect) * design$target$sesoi

  specs <- lapply(values, function(vv) {
    ov <- switch(parameter,
      attrition_rate = list(rate_control = vv, rate_treated = vv),
      reliability = list(reliability = vv),
      icc = list(icc = vv))
    list(value = vv, overrides = c(list(effect = theta), ov))
  })
  out <- fragility_rows(design, specs, sims, value_name = parameter)
  attr(out, "sims") <- sims
  attr(out, "sesoi") <- design$target$sesoi
  attr(out, "parameter") <- parameter
  class(out) <- c("recovery_fragility", "data.frame")
  out
}

fragility_rows <- function(design, specs, sims, value_name) {
  unit <- design$target$bias_scale_unit
  rows <- lapply(specs, function(sp) {
    params <- scenario_params(design, overrides = sp$overrides)
    sim_df <- run_scenario(design, params, sims)
    theta <- params$effect
    ok <- !sim_df$counted_failure & is.finite(sim_df$estimate)
    est <- sim_df$estimate[ok]
    sig <- sim_df$sig[ok]
    n_sig <- sum(sig, na.rm = TRUE)
    pw <- mean(sig, na.rm = TRUE)
    exag <- abs(est[which(sig)]) / abs(theta)
    width <- (sim_df$ci_hi - sim_df$ci_lo)[ok]
    drift <- (sim_df$theta_obs[ok] - theta) / unit
    r <- data.frame(
      value = sp$value,
      power = pw, power_mcse = mcse_prop(pw, sum(ok)),
      target_bias = mean(est - theta) / unit,
      coverage = mean(sim_df$covered[ok], na.rm = TRUE),
      drift = mean(drift), drift_mcse = mcse_mean(drift),
      type_m = if (n_sig > 0) mean(exag) else NA_real_,
      type_m_n = n_sig,
      precision = mean(width)
    )
    names(r)[1] <- value_name
    r
  })
  do.call(rbind, rows)
}

#' @export
print.recovery_fragility <- function(x, ...) {
  cat("Fragility curve over `", attr(x, "parameter"), "` (",
      attr(x, "sims"), " simulations per point; SESOI = ",
      attr(x, "sesoi"), ")\n", sep = "")
  cat("Note: fragility curves are reported separately and do not enter",
      "the PASS/RISK/FAIL verdict; they show where the verdict would",
      "change.\n\n")
  y <- as.data.frame(x)
  y[] <- lapply(y, function(col) if (is.numeric(col)) signif(col, 3) else col)
  print(y, row.names = FALSE)
  invisible(x)
}
