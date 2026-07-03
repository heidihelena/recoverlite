#' Versioned verdict threshold profiles
#'
#' The recovery-test protocol ships threshold *profiles* matched to design
#' goals (manuscript, section 2.5). Thresholds are conventions: their
#' function is not to be correct but to prevent post hoc negotiation with
#' borderline results. A profile must be selected or justified before
#' simulation; every deviation from a shipped profile is recorded and
#' echoed in the report, and the report recomputes the verdict under the
#' shipped strict and lenient profiles with signed margins to every
#' threshold.
#'
#' The shipped confirmatory profiles are, exactly:
#'
#' | Criterion             | Lenient   | Default    | Strict     |
#' |-----------------------|-----------|------------|------------|
#' | Target-null rejection | <= 1.50 a | <= 1.25 a  | <= 1.10 a  |
#' | Power                 | >= 0.70   | >= 0.80    | >= 0.90    |
#' | Target bias           | <= 0.10 D | <= 0.05 D  | <= 0.025 D |
#' | Coverage              | >= 0.900  | >= 0.925   | >= 0.940   |
#' | Type S                | <= 0.05   | <= 0.01    | <= 0.005   |
#' | Type M                | <= 2.00   | <= 1.50    | <= 1.25    |
#' | Model failure         | <= 0.05   | <= 0.01    | <= 0.005   |
#'
#' The rejection-rate threshold is one-sided: excess false claims trigger
#' failure, while conservative behavior is flagged through power and
#' precision. The coverage threshold is a lower bound; overcoverage above
#' `overcoverage_flag` (default 0.975) is flagged as inefficiency and
#' evaluated through precision, not treated as failure.
#'
#' The `"estimation"` profile replaces the rejection-rate and Type S/M
#' criteria with target bias, coverage, precision against a declared
#' maximum acceptable width (see [target_estimand()]), and estimand
#' drift, evaluated on the target rows only.
#'
#' @param profile One of `"default"`, `"strict"`, `"lenient"`,
#'   `"estimation"`.
#' @param null_rejection_mult Maximum target-null rejection rate as a
#'   multiple of the nominal alpha.
#' @param power Minimum power at the SESOI.
#' @param target_bias Maximum absolute target bias in units of the
#'   declared substantive scale Delta.
#' @param coverage Minimum coverage of a nominal 95% interval (lower
#'   bound; one-sided).
#' @param type_s Maximum Type S rate among significant estimates, when
#'   stably estimable.
#' @param type_m Maximum Type M exaggeration among significant estimates,
#'   when stably estimable.
#' @param model_failure Maximum counted model-failure rate (see
#'   [planned_analysis()] for which classes count).
#' @param drift Maximum absolute estimand drift in Delta units;
#'   thresholded by the estimation profile only (reported, not
#'   thresholded, in confirmatory profiles, where its effect is captured
#'   by target bias and target-null rejection).
#' @param overcoverage_flag Coverage above this value is flagged as
#'   inefficiency (never failure).
#' @param mcse_margin All required threshold margins must exceed this many
#'   Monte Carlo standard errors for a PASS (default 2).
#' @param min_conditional_n Minimum contributing simulations for a
#'   conditional diagnosand (Type S, Type M) to be treated as stable
#'   (default 200 — a stability convention chosen so that the MCSE of a
#'   conditional proportion cannot exceed ~0.035).
#' @param max_width Maximum acceptable confidence-interval width for the
#'   estimation profile's precision criterion; usually inherited from
#'   [target_estimand()] by [recovery_test()].
#'
#' @return An object of class `recovery_thresholds`, carrying the
#'   threshold-set version, the profile name, and a record of any values
#'   changed from the shipped profile.
#' @export
recovery_thresholds <- function(profile = c("default", "strict", "lenient",
                                            "estimation"),
                                null_rejection_mult = NULL,
                                power = NULL,
                                target_bias = NULL,
                                coverage = NULL,
                                type_s = NULL,
                                type_m = NULL,
                                model_failure = NULL,
                                drift = NULL,
                                overcoverage_flag = 0.975,
                                mcse_margin = 2,
                                min_conditional_n = 200,
                                max_width = NULL) {
  profile <- match.arg(profile)
  base <- .shipped_profiles[[if (profile == "estimation") "default" else profile]]

  values <- list(
    null_rejection_mult = null_rejection_mult %||% base$null_rejection_mult,
    power = power %||% base$power,
    target_bias = target_bias %||% base$target_bias,
    coverage = coverage %||% base$coverage,
    type_s = type_s %||% base$type_s,
    type_m = type_m %||% base$type_m,
    model_failure = model_failure %||% base$model_failure
  )
  # The drift threshold tracks the resolved bias threshold unless set.
  values$drift <- drift %||% values$target_bias
  modified <- names(base)[!mapply(identical, values[names(base)], base)]
  if (!is.null(drift) && !identical(drift, values$target_bias)) {
    modified <- union(modified, "drift")
  }

  structure(
    c(values,
      list(profile = profile,
           overcoverage_flag = overcoverage_flag,
           mcse_margin = mcse_margin,
           min_conditional_n = min_conditional_n,
           max_width = max_width,
           version = .threshold_set_version,
           modified = modified)),
    class = "recovery_thresholds"
  )
}

# Bump when shipped profiles or verdict semantics change.
.threshold_set_version <- "recoverlite-thresholds-0.2"

.shipped_profiles <- list(
  lenient = list(null_rejection_mult = 1.50, power = 0.70,
                 target_bias = 0.10, coverage = 0.900, type_s = 0.05,
                 type_m = 2.00, model_failure = 0.05),
  default = list(null_rejection_mult = 1.25, power = 0.80,
                 target_bias = 0.05, coverage = 0.925, type_s = 0.01,
                 type_m = 1.50, model_failure = 0.01),
  strict = list(null_rejection_mult = 1.10, power = 0.90,
                target_bias = 0.025, coverage = 0.940, type_s = 0.005,
                type_m = 1.25, model_failure = 0.005)
)

#' @export
print.recovery_thresholds <- function(x, ...) {
  cat("Recovery-test threshold profile '", x$profile, "' [", x$version,
      "]\n", sep = "")
  cat(sprintf("  target-null rejection      <= %.3g x alpha\n",
              x$null_rejection_mult))
  cat(sprintf("  power at the SESOI         >= %.3g\n", x$power))
  cat(sprintf("  |target bias|              <= %.3g Delta\n", x$target_bias))
  cat(sprintf("  coverage (nominal 95%%)     >= %.3g (overcoverage > %.3g flagged, not failed)\n",
              x$coverage, x$overcoverage_flag))
  cat(sprintf("  Type S rate                <= %.3g (when stably estimable)\n",
              x$type_s))
  cat(sprintf("  Type M exaggeration        <= %.3g (when stably estimable)\n",
              x$type_m))
  cat(sprintf("  model failure (counted)    <= %.3g\n", x$model_failure))
  if (x$profile == "estimation") {
    cat(sprintf("  |estimand drift|           <= %.3g Delta (estimation profile)\n",
                x$drift))
    cat(sprintf("  max acceptable CI width    %s\n",
                if (is.null(x$max_width)) "not declared"
                else sprintf("%.3g", x$max_width)))
  }
  cat(sprintf("  PASS margin                >  %g MCSE\n", x$mcse_margin))
  cat(sprintf("  conditional stability      >= %g contributing simulations\n",
              x$min_conditional_n))
  if (length(x$modified)) {
    cat("  DEVIATIONS from the shipped profile:",
        paste(x$modified, collapse = ", "), "\n")
  } else {
    cat("  All values are the shipped profile's.\n")
  }
  invisible(x)
}
