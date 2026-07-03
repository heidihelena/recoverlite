# Scenario grid construction (protocol Step 3).
#
# Nuisance degradation and effect size are CROSSED, because an analysis
# that misbehaves under realistic nuisance conditions can show it in
# false-positive behavior as easily as in power. The minimum grid for a
# confirmatory hypothesis-testing design (manuscript, section 2.2):
#
#   null_declared       theta = 0      nuisance declared      target-null rejection / calibration
#   null_pessimistic    theta = 0      nuisance pessimistic   robustness of false-claim behavior
#   target_declared     theta = SESOI  nuisance declared      power, bias, coverage, precision, drift
#   target_pessimistic  theta = SESOI  nuisance pessimistic   robustness of recovery
#   expected_effect     declared effect > SESOI, nuisance declared -> informational, outside the verdict
#
# The pessimistic scenario perturbs NUISANCE assumptions only. It is not
# a worst case; it applies a pre-specified, evidence-anchored degradation.
# Fallback package defaults (lowest tier of the evidence hierarchy):
# attrition rates x1.5 capped at the declared maximum; reliability -0.10
# bounded within (0,1]; ICC at the declared upper bound (or icc + 0.10,
# labeled as a package default, when none was declared); noncompliance
# +50% capped. The target effect is NEVER shrunk automatically; fragility
# to the effect size is a separate curve (effect_fragility()).

build_scenarios <- function(design,
                            scenarios = c("confirmatory_grid",
                                          "target_grid")) {
  scenarios <- match.arg(scenarios)
  theta_target <- sign(design$effect) * design$target$sesoi

  pess <- pessimistic_overrides(design)

  sc <- function(name, label, theta, overrides, rationale, counts_for,
                 row_type) {
    list(name = name, label = label,
         params = scenario_params(design, c(list(effect = theta), overrides)),
         rationale = rationale, counts_for = counts_for, row_type = row_type)
  }

  declared_rationale <- "Nuisance assumptions exactly as declared."
  out <- list()
  if (scenarios == "confirmatory_grid") {
    out$null_declared <- sc(
      "null_declared", "Null-declared (theta = 0, nuisance declared)",
      0, list(), paste("Target-null rejection and calibration under the",
                       "planned design.", declared_rationale),
      "declared", "null")
    out$null_pessimistic <- sc(
      "null_pessimistic", "Null-pessimistic (theta = 0, nuisance pessimistic)",
      0, pess$overrides,
      paste("Robustness of false-claim behavior.", pess$rationale),
      "pessimistic", "null")
  }
  out$target_declared <- sc(
    "target_declared",
    sprintf("Target-declared (theta = %.3g, nuisance declared)", theta_target),
    theta_target, list(),
    paste("Power, bias, coverage, precision, drift at the SESOI.",
          declared_rationale),
    "declared", "target")
  out$target_pessimistic <- sc(
    "target_pessimistic",
    sprintf("Target-pessimistic (theta = %.3g, nuisance pessimistic)",
            theta_target),
    theta_target, pess$overrides,
    paste("Robustness of recovery at the SESOI.", pess$rationale),
    "pessimistic", "target")

  # Informational expected-effect row, outside the verdict.
  if (abs(design$effect) > design$target$sesoi) {
    out$expected_effect <- sc(
      "expected_effect",
      sprintf("Expected-effect (theta = %.3g, nuisance declared; informational)",
              design$effect),
      design$effect, list(),
      sprintf(paste("Secondary planning information: the declared expected",
                    "effect (%.3g) exceeds the SESOI (%.3g). Verdict rows",
                    "use the SESOI."), design$effect, design$target$sesoi),
      "informational", "target")
  }
  attr(out, "tiers") <- pess$tiers
  out
}

# Pessimistic nuisance perturbations with evidence-tier labels.
pessimistic_overrides <- function(design) {
  ds <- design$data_strategy
  overrides <- list()
  rationale <- character(0)
  tiers <- character(0)

  if (!is.null(design$missingness) && design$missingness$rate > 0) {
    a <- design$missingness
    overrides$rate_control <- min(a$rate_control * 1.5, a$max_rate)
    overrides$rate_treated <- min(a$rate_treated * 1.5, a$max_rate)
    rationale <- c(rationale, sprintf(
      "attrition rates (control %.3g, treated %.3g) -> (%.3g, %.3g) (x1.5, capped at %.3g)",
      a$rate_control, a$rate_treated, overrides$rate_control,
      overrides$rate_treated, a$max_rate))
    tiers <- c(tiers, sprintf("attrition perturbation: %s",
      a$evidence %||% "package default (x1.5)"))
  }
  if (!is.null(design$measurement)) {
    r <- design$measurement$reliability
    overrides$reliability <- max(r - 0.10, 0.01)
    rationale <- c(rationale, sprintf(
      "reliability %.3g -> %.3g (-0.10, bounded)", r, overrides$reliability))
    tiers <- c(tiers, sprintf("reliability perturbation: %s",
      design$measurement$evidence %||% "package default (-0.10)"))
  }
  if (inherits(ds, "recovery_cluster")) {
    if (!is.null(ds$icc_pessimistic)) {
      overrides$icc <- ds$icc_pessimistic
      rationale <- c(rationale, sprintf(
        "ICC %.3g -> %.3g (declared upper bound)", ds$icc, overrides$icc))
      tiers <- c(tiers, sprintf("ICC upper bound: %s",
        ds$evidence %||% "researcher-declared (no source stated)"))
    } else {
      overrides$icc <- min(ds$icc + 0.10, 0.99)
      rationale <- c(rationale, sprintf(
        "ICC %.3g -> %.3g (NO upper bound declared; icc + 0.10)",
        ds$icc, overrides$icc))
      tiers <- c(tiers,
        "ICC upper bound: package default (icc + 0.10) - declare `icc_pessimistic` from field evidence")
    }
  }
  if (inherits(ds, "recovery_two_arm") && ds$noncompliance > 0) {
    overrides$noncompliance <- min(ds$noncompliance * 1.5, 0.95)
    rationale <- c(rationale, sprintf(
      "noncompliance %.3g -> %.3g (+50%%, capped)",
      ds$noncompliance, overrides$noncompliance))
    tiers <- c(tiers, "noncompliance perturbation: package default (+50%)")
  }

  list(
    overrides = overrides,
    rationale = if (length(rationale)) {
      paste0("Nuisance perturbations only (target effect unchanged): ",
             paste(rationale, collapse = "; "), ".")
    } else {
      "No nuisance assumptions were declared that admit a perturbation; pessimistic rows equal declared rows."
    },
    tiers = tiers
  )
}

# Resolve the design + overrides into the concrete parameter set used by
# the simulation engine for one scenario row.
scenario_params <- function(design, overrides = list()) {
  ds <- design$data_strategy
  miss <- design$missingness
  p <- list(
    type = ds$type,
    effect = overrides$effect %||% design$effect,
    reliability = overrides$reliability %||%
      (if (is.null(design$measurement)) NULL else design$measurement$reliability),
    allocation = ds$allocation
  )
  if (!is.null(miss) && miss$rate > 0) {
    p$attrition <- list(
      mechanism = miss$mechanism,
      rate_control = overrides$rate_control %||% miss$rate_control,
      rate_treated = overrides$rate_treated %||% miss$rate_treated,
      slope_control = miss$baseline_slope_control,
      slope_treated = miss$baseline_slope_treated
    )
  } else {
    p$attrition <- NULL
  }
  if (inherits(ds, "recovery_two_arm")) {
    p$n_per_arm <- ds$n_per_arm
    p$rho <- ds$baseline_outcome_cor
    p$noncompliance <- overrides$noncompliance %||% ds$noncompliance
  } else {
    p$n_clusters <- ds$n_clusters
    p$n_per_cluster <- ds$n_per_cluster
    p$icc <- overrides$icc %||% ds$icc
  }
  p
}
