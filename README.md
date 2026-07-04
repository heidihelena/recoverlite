# recoverlite

<!-- badges: start -->
[![R-CMD-check](https://github.com/heidihelena/recoverlite/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/heidihelena/recoverlite/actions/workflows/R-CMD-check.yaml)
[![recoverlite status badge](https://heidihelena.r-universe.dev/badges/recoverlite)](https://heidihelena.r-universe.dev/recoverlite)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21195920.svg)](https://doi.org/10.5281/zenodo.21195920)
[![License: Apache 2.0](https://img.shields.io/badge/License-Apache_2.0-blue.svg)](LICENSE.md)
<!-- badges: end -->

**Pre-data recovery tests for planned study designs.**

A planned study can be unable to support its intended inferential claim
even when the researcher's substantive assumptions are correct. Sampling,
measurement, missingness, assignment, and analysis may yield estimates
that are biased, poorly calibrated, unstable, or exaggerated conditional
on detection — and attrition and exclusions can change the estimand
itself: a study may retain precision while quietly answering a different
question. Conventional power analysis rarely diagnoses these failures.

`recoverlite` is a prototype implementation of the **recovery test**, a
standardized pre-data simulation protocol (Andersen, 2026, working
draft). Researchers declare the estimand, data-generating assumptions,
data strategy, missing-data process, and analysis strategy before
simulation. Data are generated under a **crossed scenario grid** — null
and target effects, each under declared and pessimistically perturbed
nuisance assumptions — and the planned analysis is applied to each
simulated dataset. The report converts the diagnosands into a
**PASS / RISK / FAIL** verdict under a pre-specified, versioned threshold
profile — a decision convention, not a validity classification.

Design declaration and data generation follow the
[DeclareDesign](https://declaredesign.org) grammar; mixed-model answer
strategies use `lme4` with `lmerTest` (Satterthwaite) and `pbkrtest`
(Kenward–Roger); `simr` is suggested for complementary GLMM power work.

## Installation

```r
# from the r-universe (recommended):
install.packages("recoverlite",
                 repos = c("https://heidihelena.r-universe.dev",
                           "https://cloud.r-project.org"))

# or the development version straight from GitHub:
# remotes::install_github("heidihelena/recoverlite")
```

Not yet on CRAN.

## The workflow in one block

```r
library(recoverlite)

design <- declare_recovery(
  target = target_estimand(
    estimand = "ITT mean difference at 12 weeks",
    scale    = "latent-outcome standardized mean difference",
    sesoi    = 0.40
  ),
  data_strategy   = two_arm_trial(n_per_arm = 115, allocation = 0.5),
  measurement     = measured_outcome(reliability = 0.70),
  missingness     = attrition_model(rate = 0.15, mechanism = "differential"),
  answer_strategy = planned_analysis(
    estimator = "linear_model",
    formula   = y_observed ~ treatment
  )
)

result <- recovery_test(design, sims = 2000,
                        scenarios = "confirmatory_grid", seed = 1)

verdict(result)   # PASS / RISK / FAIL under the selected threshold profile,
                  # recomputed under the shipped strict and lenient profiles
report(result)    # standalone recovery report; always travels with the verdict
```

Cluster-randomized designs use `cluster_trial()` with a mixed-model or
cluster-level answer strategy and an **explicit inference method** — with
few clusters, the inference method is not a detail, it is the design:

```r
design <- declare_recovery(
  target = target_estimand(
    estimand = "ITT mean difference in pupil outcome",
    scale    = "student-level standardized mean difference",
    sesoi    = 0.40
  ),
  data_strategy = cluster_trial(n_clusters = 16, n_per_cluster = 30,
                                icc = 0.05, icc_pessimistic = 0.15),
  answer_strategy = planned_analysis(
    estimator = "lmm_random_intercept",
    formula   = y_observed ~ treatment + (1 | cluster),
    inference = "kenward_roger"   # or "satterthwaite", "wald_z"
  )
)
```

## What the protocol fixes before results are known

* **A crossed scenario grid.** Null-declared, Null-pessimistic,
  Target-declared, Target-pessimistic. An analysis that misbehaves under
  realistic nuisance conditions can show it in false-positive behavior as
  easily as in power, so the null rows are required verdict rows. The
  null world is stated exactly: when the declared missingness mechanism
  displaces the analyzable contrast under the null, the rejection rate is
  reported as the **target-null rejection rate** (false claims about the
  target, partly induced by selection), not as pure test size.
* **Estimand drift as a diagnosand.** Target bias decomposes exactly into
  **estimator bias** (what the estimator does to the answer) plus
  **estimand drift** (what the data strategy does to the question). An
  unbiased estimator aimed at a displaced contrast is a design problem —
  and resources repair precision, not drift.
* **Pessimistic values by an evidence hierarchy.** Empirical ranges >
  prior-study ranges > elicited ranges > package defaults (attrition
  ×1.5 capped, reliability −0.10, ICC at its upper plausible bound,
  noncompliance +50%) — each labeled with its tier in the report. The
  target effect is never shrunk automatically; effect-size fragility is a
  separate curve (`effect_fragility()`), as are nuisance fragility curves
  (`nuisance_fragility()`), both outside the verdict.
* **A classified failure taxonomy.** Fatal errors / nonconvergence /
  degenerate (singular, boundary) fits / diagnostic warnings, reported
  separately. Fatal and nonconvergence always count against the failure
  threshold; whether degenerate fits count is pre-specified in
  `planned_analysis()`, and their marginal effect on coverage is
  reported.
* **Monte Carlo uncertainty in the verdict.** Every diagnosand carries an
  MCSE (bootstrap for conditional diagnosands); conditional diagnosands
  report contributing counts and are marked unstable below 200; a margin
  within 2 MCSE of its threshold caps the verdict at RISK.
* **Threshold profiles, not thresholds.** Shipped lenient / default /
  strict profiles (and an estimation profile); the report shows the
  signed margin to every threshold and recomputes the verdict under
  strict and lenient. A verdict that flips across profiles is itself a
  finding — the RISK category exists to hold it.

| Verdict | Rule |
|---|---|
| **PASS** | All required thresholds met under all scenario rows the profile requires, every margin > 2 MCSE. |
| **RISK** | Passes declared-nuisance rows but fails a pessimistic row, **or** any margin within 2 MCSE, **or** a required conditional diagnosand too unstable to confirm. |
| **FAIL** | Any required threshold fails under a declared-nuisance row — including an inflated target-null rejection rate. |

Default confirmatory profile: target-null rejection ≤ 1.25α, power ≥ .80
at the SESOI, |target bias| ≤ .05Δ, coverage ≥ .925 (overcoverage > .975
flagged as inefficiency, not failure), Type S ≤ .01, Type M ≤ 1.50,
model failure ≤ .01.

A PASS is evidence about the instrument, not about the world.

## Scope of the current prototype

**Supported:** two-arm parallel trials with an observed baseline,
classical additive measurement error, MCAR or baseline-dependent (MAR)
differential attrition, optional one-sided noncompliance; answer
strategies: complete-case linear model, baseline-adjusted multiple
imputation (`mi_baseline_adjusted`); cluster-randomized parallel trials
with random-intercept mixed models (Wald z / Satterthwaite /
Kenward–Roger inference) and the cluster-level t-test.

**Unsupported:** crossed and longitudinal random-effects structures,
Bayesian answer strategies, prediction models, and latent-variable
measurement models (manuscript §5.4).

The scripts that reproduce the paper's worked examples are in
[`inst/paper/`](inst/paper/). Agent-facing usage instructions are in
[`SKILL.md`](SKILL.md).

## Citation

```r
citation("recoverlite")
```

> Andersen, H. H. (2026). *Recovery before data: pre-data simulation
> diagnosis of planned study designs.* Working paper; preprint
> forthcoming. https://github.com/heidihelena/recoverlite

Versioned releases are archived on Zenodo:
[doi:10.5281/zenodo.21195920](https://doi.org/10.5281/zenodo.21195920)
(concept DOI, always resolves to the latest version).

The reusable methods sentence for preregistrations and grant
applications:

> "Design feasibility was evaluated using a pre-data recovery test, in
> which the planned design and analysis were simulated under declared
> assumptions and pessimistic perturbations to assess power, bias,
> coverage, precision, and model stability."

## License

[Apache License 2.0](LICENSE.md).
