#' recoverlite: Pre-Data Recovery Tests for Planned Study Designs
#'
#' A prototype implementation of the recovery-test protocol (Andersen,
#' 2026, working draft): declare the target estimand and the planned
#' design-analysis pair, simulate it over a crossed scenario grid — null
#' and target effects, each under declared and pessimistically perturbed
#' nuisance assumptions — and convert the diagnosands into a
#' PASS/RISK/FAIL verdict under a pre-specified, versioned threshold
#' profile, with Monte Carlo uncertainty part of the verdict.
#'
#' The core workflow is:
#' 1. [target_estimand()] — state the estimand, its scale, and the SESOI.
#' 2. [declare_recovery()] — assemble the target, data strategy
#'    ([two_arm_trial()] or [cluster_trial()]), measurement
#'    ([measured_outcome()]), missingness ([attrition_model()]), and
#'    answer strategy ([planned_analysis()]).
#' 3. [recovery_test()] — simulate the scenario grid and compute
#'    diagnosands (rejection rate / power, target bias with its exact
#'    decomposition into estimator bias and estimand drift, coverage,
#'    Type S/M, precision, classified model failure) with MCSEs.
#' 4. [verdict()] — the PASS/RISK/FAIL rule under the selected profile,
#'    recomputed under the shipped strict and lenient profiles.
#' 5. [report()] — the standalone recovery report.
#'
#' Fragility is deliberately outside the verdict; see
#' [effect_fragility()] and [nuisance_fragility()].
#'
#' @keywords internal
#' @importFrom stats lm qt qnorm pnorm pt qlogis plogis dnorm rnorm rbinom
#'   runif sd var integrate uniroot rchisq coef vcov model.matrix lm.fit
#' @importFrom utils packageVersion modifyList
"_PACKAGE"

# Data-mask variables used inside DeclareDesign/fabricatr declarations
# (evaluated in the simulated data, not the package namespace).
utils::globalVariables(c(
  "N", "baseline", "eps", "complier", "treatment", "y_true", "y_observed",
  "p_drop", "retained", "u_c", "e", "cluster",
  "y_true_treatment_0", "y_true_treatment_1"
))

`%||%` <- function(x, y) if (is.null(x)) y else x
